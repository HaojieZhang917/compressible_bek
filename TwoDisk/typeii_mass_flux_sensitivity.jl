using LinearAlgebra
using Printf
using DelimitedFiles
using BSplineKit
using PyCall

include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")

const RE_H = 1000.0
const RO = -1.0
const MODE = 1
const N_CHEB = parse(Int, get(ENV, "R3C5_N_CHEB", "99"))
const EPSILON_AS = parse(Float64, get(ENV, "R3C5_EPSILON_AS", "0.01"))

struct CriticalMode
    name::String
    radius::Float64
    beta::Float64
    alpha::Float64
    omega_target::ComplexF64
end

const MODES = [
    CriticalMode("Type I", 283.21, 0.0824, 0.403, 0.0099 + 0.0im),
    CriticalMode("Type II", 91.08, -0.1135, 0.256, 0.0983 + 0.0im),
]

const AS_VALUES = [-0.4, -0.2, 0.0, 0.2, 0.4]
const RC_TYPE_I = [187.55, 228.41, 283.21, 353.16, 453.06]
const RC_TYPE_II = [83.10, 84.96, 91.08, 102.23, 121.75]

function reduced_indices(n_cheb)
    nfull = 4 * (n_cheb + 1)
    removed = (1, n_cheb + 1, n_cheb + 2, 2n_cheb + 2,
               2n_cheb + 3, 3n_cheb + 3)
    return setdiff(1:nfull, removed)
end

function base_profile(a_s, D, z)
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(RE_H, RO, a_s, MODE)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, N_CHEB, MODE)
    return vec(F), vec(G .- 1), vec(H)
end

function temporal_problem(F, G, H, mode::CriticalMode, D, D2)
    cof = CRC_STA.Spatial_mode_BEK1(
        reshape(F, :, 1), reshape(G, :, 1), reshape(H, :, 1),
        mode.radius, N_CHEB, D, D2, RE_H,
    )
    return CRC_STA.assemble_time_mat(
        cof, D, D2, mode.beta, mode.alpha, mode.radius, N_CHEB,
    )
end

function finite_eigenpairs(H0, H1)
    decomposition = eigen(H0, H1)
    keep = findall(v -> isfinite(real(v)) && isfinite(imag(v)), decomposition.values)
    return decomposition.values[keep], decomposition.vectors[:, keep]
end

function nearest_pair(H0, H1, target)
    values, vectors = finite_eigenpairs(H0, H1)
    index = argmin(abs.(values .- target))
    return values[index], vectors[:, index]
end

function paired_left_vector(H0, H1, omega)
    values, vectors = finite_eigenpairs(adjoint(H0), adjoint(H1))
    index = argmin(abs.(values .- conj(omega)))
    return values[index], vectors[:, index]
end

function block_range(component, n)
    first = (component - 1) * n + 1
    return first:(first + n - 1)
end

function profile_operator_components(F, G, H, mode::CriticalMode, D)
    n = length(F)
    total = 4n
    sqrt_re = sqrt(RE_H)
    R = mode.radius
    beta = mode.beta
    alpha = mode.alpha
    KD = kron(Matrix{Float64}(I, 4, 4), D)

    ru = block_range(1, n)
    rv = block_range(2, n)
    rw = block_range(3, n)

    f_value = zeros(ComplexF64, total, total)
    g_value = zeros(ComplexF64, total, total)
    h_value = zeros(ComplexF64, total, total)
    f_shear = zeros(ComplexF64, total, total)
    g_shear = zeros(ComplexF64, total, total)
    h_shear = zeros(ComplexF64, total, total)

    # F enters the streamwise-advection matrix A and the diagonal D1 blocks.
    for range in (ru, rv, rw)
        f_value[range, range] .+= im * alpha .* Diagonal(F)
    end
    f_value[ru, ru] .+= Diagonal(F ./ R)
    f_value[rv, rv] .+= Diagonal(F ./ R)

    # G enters the azimuthal-advection matrix B and the Coriolis/base-flow blocks.
    for range in (ru, rv, rw)
        g_value[range, range] .+= im * beta .* Diagonal(G)
    end
    g_value[ru, rv] .+= Diagonal(-2 .* G ./ R)
    g_value[rv, ru] .+= Diagonal(2 .* G ./ R)

    # H multiplies the wall-normal derivative in C*D.
    c_h = zeros(ComplexF64, total, total)
    for range in (ru, rv, rw)
        c_h[range, range] .+= Diagonal(H ./ (R * sqrt_re))
    end
    h_value .= c_h * KD

    # The remaining base-flow dependence is through F', G', and H'.
    f_shear[ru, rw] .+= Diagonal((D * F) ./ sqrt_re)
    g_shear[rv, rw] .+= Diagonal((D * G) ./ sqrt_re)
    h_shear[rw, rw] .+= Diagonal((D * H) ./ (R * sqrt_re))

    return Dict(
        "F value" => f_value,
        "G value" => g_value,
        "H value" => h_value,
        "F shear" => f_shear,
        "G shear" => g_shear,
        "H shear" => h_shear,
    )
end

function mode_sensitivity(mode::CriticalMode, epsilon, D, D2, z)
    F0, G0, H0profile = base_profile(0.0, D, z)
    Fp, Gp, Hp = base_profile(epsilon, D, z)
    Fm, Gm, Hm = base_profile(-epsilon, D, z)

    H00, H1 = temporal_problem(F0, G0, H0profile, mode, D, D2)
    H0p, H1p = temporal_problem(Fp, Gp, Hp, mode, D, D2)
    H0m, H1m = temporal_problem(Fm, Gm, Hm, mode, D, D2)

    omega, direct = nearest_pair(H00, H1, mode.omega_target)
    adjoint_value, left = paired_left_vector(H00, H1, omega)
    omega_p, _ = nearest_pair(H0p, H1p, omega)
    omega_m, _ = nearest_pair(H0m, H1m, omega)

    full_derivative = (H0p - H0m) / (2epsilon)
    plus_components = profile_operator_components(Fp, Gp, Hp, mode, D)
    minus_components = profile_operator_components(Fm, Gm, Hm, mode, D)
    indices = reduced_indices(N_CHEB)

    component_derivatives = Dict{String, Matrix{ComplexF64}}()
    for key in keys(plus_components)
        component_derivatives[key] =
            (plus_components[key][indices, indices] -
             minus_components[key][indices, indices]) / (2epsilon)
    end

    reconstructed = sum(values(component_derivatives))
    reconstruction_error = norm(full_derivative - reconstructed) / norm(full_derivative)
    denominator = dot(left, H1 * direct)
    sensitivities = Dict(
        key => dot(left, matrix * direct) / denominator
        for (key, matrix) in component_derivatives
    )
    adjoint_total = dot(left, full_derivative * direct) / denominator
    direct_fd = (omega_p - omega_m) / (2epsilon)

    value_total = sum(sensitivities[key] for key in ("F value", "G value", "H value"))
    shear_total = sum(sensitivities[key] for key in ("F shear", "G shear", "H shear"))

    return (
        mode = mode,
        omega = omega,
        adjoint_value = adjoint_value,
        omega_p = omega_p,
        omega_m = omega_m,
        direct_fd = direct_fd,
        adjoint_total = adjoint_total,
        value_total = value_total,
        shear_total = shear_total,
        components = sensitivities,
        reconstruction_error = reconstruction_error,
        eigen_residual = norm(H00 * direct - omega * H1 * direct) / norm(H00 * direct),
        adjoint_pair_error = abs(adjoint_value - conj(omega)),
    )
end

function critical_radius_summary()
    rows = NamedTuple[]
    for (name, values) in (("Type I", RC_TYPE_I), ("Type II", RC_TYPE_II))
        central_slope = (values[4] - values[2]) / (AS_VALUES[4] - AS_VALUES[2])
        endpoint_slope = (values[end] - values[1]) / (AS_VALUES[end] - AS_VALUES[1])
        relative_change = (values[end] - values[1]) / values[3]
        push!(rows, (
            name = name,
            central_slope = central_slope,
            endpoint_slope = endpoint_slope,
            relative_change = relative_change,
        ))
    end
    return rows
end

function main()
    output_dir = joinpath(@__DIR__, "r3_c5_validation_results",
                          "N$(N_CHEB)_eps$(replace(string(EPSILON_AS), "." => "p"))")
    mkpath(output_dir)
    D, D2, z_matrix = CRC_BF.Cheb(N_CHEB, MODE)
    z = vec(z_matrix)
    epsilon = EPSILON_AS

    results = [mode_sensitivity(mode, epsilon, D, D2, z) for mode in MODES]
    summary = critical_radius_summary()

    open(joinpath(output_dir, "mass_flux_sensitivity_summary.txt"), "w") do io
        println(io, "Reviewer 3 Comment 5: quantitative Type-II mass-flux sensitivity")
        println(io, "N_cheb = $N_CHEB; N_z = $(N_CHEB + 1); epsilon_as = $epsilon")
        println(io)
        println(io, "Critical-radius trends from the manuscript table")
        for row in summary
            @printf(io, "%s: central dRc/das = %.8f; endpoint dRc/das = %.8f; ",
                    row.name, row.central_slope, row.endpoint_slope)
            @printf(io, "endpoint DeltaRc/Rc(as=0) = %.8f\n", row.relative_change)
        end
        println(io)

        for result in results
            println(io, "[$(result.mode.name)]")
            @printf(io, "critical parameters: R=%.8f beta=%.8f alpha=%.8f\n",
                    result.mode.radius, result.mode.beta, result.mode.alpha)
            @printf(io, "omega(as=0) = %.12e %+.12ei\n", real(result.omega), imag(result.omega))
            @printf(io, "omega(+eps) = %.12e %+.12ei\n", real(result.omega_p), imag(result.omega_p))
            @printf(io, "omega(-eps) = %.12e %+.12ei\n", real(result.omega_m), imag(result.omega_m))
            @printf(io, "direct FD domega/das = %.12e %+.12ei\n",
                    real(result.direct_fd), imag(result.direct_fd))
            @printf(io, "adjoint domega/das = %.12e %+.12ei\n",
                    real(result.adjoint_total), imag(result.adjoint_total))
            @printf(io, "profile-value total = %.12e %+.12ei\n",
                    real(result.value_total), imag(result.value_total))
            @printf(io, "profile-shear total = %.12e %+.12ei\n",
                    real(result.shear_total), imag(result.shear_total))
            for key in ("F value", "G value", "H value", "F shear", "G shear", "H shear")
                value = result.components[key]
                @printf(io, "  %-8s = %.12e %+.12ei\n", key, real(value), imag(value))
            end
            @printf(io, "operator reconstruction error = %.12e\n", result.reconstruction_error)
            @printf(io, "direct eigen residual = %.12e\n", result.eigen_residual)
            @printf(io, "adjoint pairing error = %.12e\n\n", result.adjoint_pair_error)
        end
    end

    rows = Any[]
    for result in results
        push!(rows, [
            result.mode.name,
            real(result.omega), imag(result.omega),
            real(result.direct_fd), imag(result.direct_fd),
            real(result.adjoint_total), imag(result.adjoint_total),
            real(result.value_total), imag(result.value_total),
            real(result.shear_total), imag(result.shear_total),
        ])
    end
    open(joinpath(output_dir, "mass_flux_sensitivity.dat"), "w") do io
        println(io, "TITLE = \"Type-I and Type-II mass-flux sensitivity\"")
        println(io, "VARIABLES = \"Mode\" \"omega_r\" \"omega_i\" \"FD_domega_r\" \"FD_domega_i\" \"Adj_domega_r\" \"Adj_domega_i\" \"Value_domega_r\" \"Value_domega_i\" \"Shear_domega_r\" \"Shear_domega_i\"")
        for row in rows
            println(io, join(row, " "))
        end
    end

    println(read(joinpath(output_dir, "mass_flux_sensitivity_summary.txt"), String))
end

main()
