#!/usr/bin/env julia

"""Dynamical-systems analysis of the three Re_h=1000 similarity states.

Spatial eigenvalues describe the autonomous steady ODE in z. Temporal
eigenvalues describe axisymmetric, similarity-preserving perturbations; they
are not a replacement for the full three-dimensional stability problem.
"""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using CSV
using JSON3
using LinearAlgebra
using Plots
using Printf
using Statistics

const DEFAULT_INPUT = joinpath(@__DIR__, "..", "data", "boussinesq_singularity_results",
                               "three_solutions_Tw1.160.csv")
const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "data", "three_solution_dynamics")
const BRANCH_LABELS = Dict(1 => "upper / low-Pi", 2 => "middle",
                           3 => "principal / isothermal-connected")

struct BaseProfile
    branch::Int
    pressure::Float64
    z::Vector{Float64}
    h::Vector{Float64}
    f::Vector{Float64}
    g::Vector{Float64}
    temperature::Vector{Float64}
end

struct TemporalResult
    n::Int
    z::Vector{Float64}
    eigenvalues::Vector{ComplexF64}
    eigenvectors::Matrix{ComplexF64}
    f_map::Matrix{Float64}
    h_map::Matrix{Float64}
    g_map::Matrix{Float64}
    temperature_map::Matrix{Float64}
    n_f_reduced::Int
    n_g::Int
end

function load_profiles(path::AbstractString)
    rows = read_csv_rows(path)
    branches = sort(unique(Int(round(numeric_value(row["branch"]))) for row in rows))
    profiles = BaseProfile[]
    for branch in branches
        selected = filter(row -> Int(round(numeric_value(row["branch"]))) == branch, rows)
        sort!(selected; by=row -> numeric_value(row["z"]))
        column(name) = [numeric_value(row[name]) for row in selected]
        push!(profiles, BaseProfile(branch, numeric_value(first(selected)["pressure_gradient"]),
                                    column("z"), column("H"), column("F"), column("G"),
                                    column("T")))
    end
    return profiles
end

"""Chebyshev antiderivative matrix on [0,1], with every primitive zero at z=0."""
function chebyshev_integration_matrix(z::AbstractVector)
    n = length(z)
    x = 2 .* z .- 1
    theta = acos.(clamp.(x, -1.0, 1.0))
    vandermonde = hcat((cos.(degree .* theta) for degree in 0:n-1)...)
    integral_basis = zeros(n, n)
    integral_basis[:, 1] .= 0.5 .* (x .+ 1)
    if n >= 2
        integral_basis[:, 2] .= 0.25 .* (x.^2 .- 1)
    end
    for degree in 2:n-1
        tplus = cos.((degree + 1) .* theta)
        tminus = cos.((degree - 1) .* theta)
        value_at_minus_one = 0.5 * (((-1.0)^(degree + 1)) / (degree + 1) -
                                    ((-1.0)^(degree - 1)) / (degree - 1))
        integral_basis[:, degree + 1] .= 0.5 .* (
            0.5 .* (tplus ./ (degree + 1) .- tminus ./ (degree - 1)) .-
            value_at_minus_one)
    end
    return integral_basis / vandermonde
end

function interpolate_profile(profile::BaseProfile, z::AbstractVector,
                             D::AbstractMatrix, D2::AbstractMatrix)
    values = Dict{String,Vector{Float64}}(
        "h" => linear_interpolate(profile.z, profile.h, z),
        "f" => linear_interpolate(profile.z, profile.f, z),
        "g" => linear_interpolate(profile.z, profile.g, z),
        "temperature" => linear_interpolate(profile.z, profile.temperature, z),
    )
    for name in ("h", "f", "g", "temperature")
        values[name * "_z"] = D * values[name]
        values[name * "_zz"] = D2 * values[name]
    end
    return values
end

function temporal_operator(profile::BaseProfile, n::Int, re_h::Real, pr::Real)
    z, d1, d2, _ = chebyshev_operators(n - 1, 0.0, 1.0)
    integration = chebyshev_integration_matrix(z)
    weights = vec(integration[end, :])
    base = interpolate_profile(profile, z, d1, d2)
    interior = 2:n-1
    m = length(interior)
    sqrt_re = sqrt(re_h)

    extension = zeros(n, m)
    extension[interior, :] .= Matrix{Float64}(I, m, m)
    weights_i = weights[interior]
    q = nullspace(reshape(weights_i, 1, :))
    size(q) == (m, m - 1) || error("Unexpected zero-integral basis dimension")
    pressure_projection = Matrix{Float64}(I, m, m) .-
                          ones(m) * transpose(weights_i ./ sum(weights_i))
    f_map = extension * q
    g_map = extension
    temperature_map = extension
    h_map = -2sqrt_re .* integration * f_map

    f_i, g_i, h_i = (base[name][interior] for name in ("f", "g", "h"))
    f_z_i = base["f_z"][interior]
    g_z_i = base["g_z"][interior]
    temperature_z_i = base["temperature_z"][interior]
    d1_f, d2_f = (matrix * f_map for matrix in (d1, d2))
    d1_g, d2_g = (matrix * g_map for matrix in (d1, d2))
    d1_t, d2_t = (matrix * temperature_map for matrix in (d1, d2))
    f_i_map = f_map[interior, :]
    h_i_map = h_map[interior, :]

    radial_f = d2_f[interior, :] ./ re_h .-
               h_i .* d1_f[interior, :] ./ sqrt_re .-
               f_z_i .* h_i_map ./ sqrt_re .- 2f_i .* f_i_map
    radial_g = 2 .* Diagonal(g_i)
    radial_t = -Matrix{Float64}(I, m, m)
    azimuthal_f = -g_z_i .* h_i_map ./ sqrt_re .- 2g_i .* f_i_map
    azimuthal_g = d2_g[interior, :] ./ re_h .-
                   h_i .* d1_g[interior, :] ./ sqrt_re .- 2 .* Diagonal(f_i)
    thermal_f = -temperature_z_i .* h_i_map ./ sqrt_re
    thermal_t = d2_t[interior, :] ./ (re_h * pr) .-
                h_i .* d1_t[interior, :] ./ sqrt_re

    n_a = m - 1
    operator = [transpose(q) * pressure_projection * radial_f  transpose(q) * pressure_projection * radial_g  transpose(q) * pressure_projection * radial_t;
                azimuthal_f                                    azimuthal_g                                    zeros(m, m);
                thermal_f                                      zeros(m, m)                                    thermal_t]
    spectrum = eigen(operator)
    order = sortperm(eachindex(spectrum.values);
                     by=i -> (-real(spectrum.values[i]), abs(imag(spectrum.values[i]))))
    result = TemporalResult(n, z, ComplexF64.(spectrum.values[order]),
                            ComplexF64.(spectrum.vectors[:, order]), f_map, h_map,
                            g_map, temperature_map, n_a, m)
    return operator, result
end

function spatial_jacobian(state::AbstractVector, re_h::Real, pr::Real)
    h, f, fp, g, gp, _, tp = state
    sqrt_re = sqrt(re_h)
    matrix = zeros(7, 7)
    matrix[1, 2] = -2sqrt_re
    matrix[2, 3] = 1
    matrix[3, 1] = sqrt_re * fp
    matrix[3, 2] = 2re_h * f
    matrix[3, 3] = sqrt_re * h
    matrix[3, 4] = -2re_h * g
    matrix[3, 6] = re_h
    matrix[4, 5] = 1
    matrix[5, 1] = sqrt_re * gp
    matrix[5, 2] = 2re_h * g
    matrix[5, 4] = 2re_h * f
    matrix[5, 5] = sqrt_re * h
    matrix[6, 7] = 1
    matrix[7, 1] = pr * sqrt_re * tp
    matrix[7, 7] = pr * sqrt_re * h
    return matrix
end

trapz(values, z) = sum((values[1:end-1] .+ values[2:end]) .* diff(z)) / 2

function sign_changes(values; tolerance=2e-4)
    filtered = values[abs.(values) .> tolerance]
    length(filtered) < 2 && return 0
    return count(!=(0), diff(sign.(filtered)))
end

function branch_diagnostics(profile::BaseProfile, re_h::Real, pr::Real)
    z = collect(range(0.0, 1.0; length=1001))
    _, d1, d2, _ = chebyshev_operators(200, 0.0, 1.0)
    # Diagnostics use the dense stored profile directly to avoid differentiating
    # a 1001-by-1001 matrix. Central differences are sufficient for reporting.
    h = linear_interpolate(profile.z, profile.h, z)
    f = linear_interpolate(profile.z, profile.f, z)
    g = linear_interpolate(profile.z, profile.g, z)
    t = linear_interpolate(profile.z, profile.temperature, z)
    dz = z[2] - z[1]
    derivative(v) = vcat((v[2] - v[1]) / dz,
                         (v[3:end] .- v[1:end-2]) ./ (2dz),
                         (v[end] - v[end-1]) / dz)
    second(v) = vcat((v[3] - 2v[2] + v[1]) / dz^2,
                     (v[3:end] .- 2v[2:end-1] .+ v[1:end-2]) ./ dz^2,
                     (v[end] - 2v[end-1] + v[end-2]) / dz^2)
    fz, gz, tz = derivative(f), derivative(g), derivative(t)
    fzz, gzz, tzz = second(f), second(g), second(t)
    core = findall(@. 0.30 <= z <= 0.70)
    inviscid = g.^2 .- (t .- 1) .- profile.pressure
    radial = fzz ./ re_h .- h .* fz ./ sqrt(re_h) .+ g.^2 .- f.^2 .-
             (t .- 1) .- profile.pressure
    azimuthal = gzz ./ re_h .- h .* gz ./ sqrt(re_h) .- 2f .* g
    thermal = tzz ./ (re_h * pr) .- h .* tz ./ sqrt(re_h)
    mid = argmin(abs.(z .- 0.5))
    state = [h[mid], f[mid], fz[mid], g[mid], gz[mid], t[mid], tz[mid]]
    spatial = ComplexF64.(eigvals(spatial_jacobian(state, re_h, pr)))
    complex_spatial = filter(value -> abs(imag(value)) > 1e-8, spatial)
    min_wavelength = isempty(complex_spatial) ? NaN :
                     minimum(2pi ./ abs.(imag.(complex_spatial)))
    row = Dict{String,Any}(
        "branch" => profile.branch, "label" => BRANCH_LABELS[profile.branch],
        "pressure_gradient" => profile.pressure, "H_mid" => h[mid],
        "F_mid" => f[mid], "G_mid" => g[mid], "T_mid" => t[mid],
        "core_max_abs_F" => maximum(abs.(f[core])),
        "core_G_range" => maximum(g[core]) - minimum(g[core]),
        "core_rms_Gz" => sqrt(mean(gz[core].^2)),
        "core_rms_inviscid_residual" => sqrt(mean(inviscid[core].^2)),
        "mass_integral_F" => trapz(f, z),
        "steady_radial_residual_max" => maximum(abs.(radial[2:end-1])),
        "steady_azimuthal_residual_max" => maximum(abs.(azimuthal[2:end-1])),
        "steady_thermal_residual_max" => maximum(abs.(thermal[2:end-1])),
        "F_interior_sign_changes" => sign_changes(f[2:end-1]),
        "G_interior_extrema" => sign_changes(gz[2:end-1]; tolerance=2e-3),
        "mid_spatial_positive_real_count" => count(>(1e-8), real.(spatial)),
        "mid_spatial_negative_real_count" => count(<(-1e-8), real.(spatial)),
        "mid_spatial_neutral_count" => count(<=(1e-8), abs.(real.(spatial))),
        "shortest_mid_spatial_wavelength" => min_wavelength,
    )
    return row, spatial
end

function plot_temporal_spectra(results, output)
    panels = Any[]
    for branch in sort(collect(keys(results)))
        values = results[branch].eigenvalues
        mask = @. real(values) > -2.5 && abs(imag(values)) < 8
        lead = first(values)
        panel = scatter(real.(values[mask]), imag.(values[mask]); markersize=2.5,
                        label=false, xlabel="Re(lambda)/Omega", ylabel="Im(lambda)/Omega",
                        title=@sprintf("Branch %d\n%.4f%+.4fi", branch,
                                       real(lead), imag(lead)), gridalpha=0.25)
        vline!(panel, [0.0]; color=:black, linewidth=1, label=false)
        scatter!(panel, [real(lead)], [imag(lead)]; marker=:star5, markersize=8,
                 color=:red, label=false)
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(1, 3), size=(1500, 470)), output)
end

function plot_dominant_modes(results, output)
    panels = Any[]
    for branch in sort(collect(keys(results)))
        result = results[branch]
        vector = result.eigenvectors[:, 1]
        na, ng = result.n_f_reduced, result.n_g
        components = [result.f_map * vector[1:na],
                      result.g_map * vector[na+1:na+ng],
                      result.h_map * vector[1:na],
                      result.temperature_map * vector[na+ng+1:end]]
        scale = maximum(maximum(abs, item) for item in components)
        panel = plot(; xlabel="z/h", ylabel="normalised amplitude",
                     title="Branch $branch: $(BRANCH_LABELS[branch])", gridalpha=0.25)
        for (values, label) in zip(components, ("|f|", "|g|", "|h|", "|theta|"))
            plot!(panel, result.z, abs.(values) ./ scale; label)
        end
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(1, 3), size=(1500, 470)), output)
end

function plot_spatial_dynamics(profiles, output)
    p1 = plot(; xlabel="F", ylabel="F_z", title="Radial phase projection")
    p2 = plot(; xlabel="G", ylabel="G_z", title="Azimuthal phase projection")
    p3 = plot(; xlabel="z/h", ylabel="G^2-(T-1)-Pi",
              title="Departure from inviscid-core balance")
    for profile in profiles
        z = collect(range(0.0, 1.0; length=1001))
        f = linear_interpolate(profile.z, profile.f, z)
        g = linear_interpolate(profile.z, profile.g, z)
        t = linear_interpolate(profile.z, profile.temperature, z)
        dz = z[2] - z[1]
        fz = vcat(diff(f) ./ dz, (f[end] - f[end-1]) / dz)
        gz = vcat(diff(g) ./ dz, (g[end] - g[end-1]) / dz)
        plot!(p1, f, fz; label="B$(profile.branch): Pi=$(round(profile.pressure; digits=5))")
        plot!(p2, g, gz; label=false)
        plot!(p3, z, g.^2 .- (t .- 1) .- profile.pressure; label=false)
    end
    hline!(p3, [0.0]; color=:black, label=false)
    savefig(plot(p1, p2, p3; layout=(1, 3), size=(1550, 470)), output)
end

function scan_temporal_growth_along_branch(re_h, pr, order)
    config = RotorConfig(re_h=re_h, pr=pr, tolerance=2e-7,
                         collocation_degree=max(order - 1, 60))
    reference = joinpath(@__DIR__, "baseflow_Res1000.npz")
    isothermal = solve_rotor_isothermal(config, reference)
    folds_path = joinpath(@__DIR__, "boussinesq_singularity_results",
                          "folds_Re$(re_h)_traditional_centrifugal.json")
    folds = isfile(folds_path) ? JSON3.read(read(folds_path, String)) : Any[]
    fold_pressures = [Float64(item.pressure_gradient) for item in folds]
    pressures = sort(unique(vcat(collect(isothermal.pressure:-0.002:0.015),
                                 collect(0.014:-0.001:0.011),
                                 collect(0.0105:-0.0005:-0.001), fold_pressures)); rev=true)
    branch = continue_rotor_pressure(pressures, isothermal; profile_points=1001)
    rows = Dict{String,Any}[]
    for profile in branch.profiles
        z = collect(range(0.0, 1.0; length=1001))
        state = rotor_state(profile, z)
        base = BaseProfile(0, profile.pressure, z, vec(state[1, :]), vec(state[2, :]),
                           vec(state[4, :]), vec(state[6, :]))
        _, temporal = temporal_operator(base, order, re_h, pr)
        lead = first(temporal.eigenvalues)
        push!(rows, Dict("pressure_gradient" => profile.pressure, "Tw" => profile.tw,
                         "leading_real" => real(lead), "leading_imag" => imag(lead),
                         "positive_real_count" => count(>(1e-7), real.(temporal.eigenvalues)),
                         "fold_sample" => Int(any(abs.(fold_pressures .- profile.pressure) .< 1e-12))))
    end
    return rows
end

function plot_branch_growth(rows, output)
    pressure = [row["pressure_gradient"] for row in rows]
    tw = [row["Tw"] for row in rows]
    growth = [row["leading_real"] for row in rows]
    stable = growth .<= 0
    folds = [row["fold_sample"] == 1 for row in rows]
    p1 = plot(tw, pressure; color=:grey, label=false, xlabel="Tw", ylabel="Pi",
              title="Stability exchange on the S-curve")
    scatter!(p1, tw[stable], pressure[stable]; label="similarity-stable", color=:green)
    scatter!(p1, tw[.!stable], pressure[.!stable]; label="unstable", color=:red)
    scatter!(p1, tw[folds], pressure[folds]; marker=:star5, color=:black, label="folds")
    p2 = plot(pressure, growth; marker=:circle, markersize=2, label=false,
              xlabel="Pi", ylabel="max Re(lambda)/Omega", title="Leading growth rate")
    hline!(p2, [0.0]; color=:black, label=false)
    savefig(plot(p1, p2; layout=(1, 2), size=(1250, 480)), output)
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}(
        "input" => DEFAULT_INPUT, "output-dir" => DEFAULT_OUTPUT,
        "re-h" => 1000.0, "pr" => 0.72,
        "collocation-orders" => [61, 81, 101, 121],
        "branch-scan-order" => 81, "skip-branch-scan" => false)
    args = parse_cli(arguments, defaults)
    output = abspath(args["output-dir"])
    mkpath(output)
    profiles = load_profiles(abspath(args["input"]))

    diagnostic_rows = Dict{String,Any}[]
    spatial_rows = Dict{String,Any}[]
    for profile in profiles
        row, values = branch_diagnostics(profile, args["re-h"], args["pr"])
        push!(diagnostic_rows, row)
        ordered = sort(values; by=value -> (-real(value), abs(imag(value))))
        for (rank, value) in enumerate(ordered)
            push!(spatial_rows, Dict("branch" => profile.branch,
                 "label" => BRANCH_LABELS[profile.branch], "rank" => rank,
                 "real" => real(value), "imag" => imag(value),
                 "spatial_wavelength" => abs(imag(value)) > 1e-10 ? 2pi / abs(imag(value)) : NaN))
        end
    end

    temporal_rows = Dict{String,Any}[]
    final_results = Dict{Int,TemporalResult}()
    convergence = Dict{String,Any}()
    orders = args["collocation-orders"]
    for profile in profiles
        branch_convergence = Dict{String,Any}()
        for n in orders
            _, result = temporal_operator(profile, n, args["re-h"], args["pr"])
            lead = first(result.eigenvalues)
            branch_convergence[string(n)] = Dict("leading_real" => real(lead),
                "leading_imag" => imag(lead),
                "positive_real_count" => count(>(1e-7), real.(result.eigenvalues)))
            for (rank, value) in enumerate(result.eigenvalues[1:min(20, end)])
                push!(temporal_rows, Dict("branch" => profile.branch,
                    "label" => BRANCH_LABELS[profile.branch], "collocation_points" => n,
                    "rank" => rank, "real" => real(value), "imag" => imag(value),
                    "oscillation_frequency_over_Omega" => abs(imag(value)),
                    "period_Omega_t" => abs(imag(value)) > 1e-10 ? 2pi / abs(imag(value)) : NaN))
            end
            n == last(orders) && (final_results[profile.branch] = result)
        end
        convergence[string(profile.branch)] = branch_convergence
    end

    for row in diagnostic_rows
        result = final_results[Int(row["branch"])]
        spectrum_values = result.eigenvalues
        lead = first(spectrum_values)
        vector = result.eigenvectors[:, 1]
        na, ng = result.n_f_reduced, result.n_g
        components = Dict(
            "F" => result.f_map * vector[1:na],
            "H" => result.h_map * vector[1:na],
            "G" => result.g_map * vector[na+1:na+ng],
            "T" => result.temperature_map * vector[na+ng+1:end])
        maxima = Dict(name => maximum(abs, values) for (name, values) in components)
        scale = maximum(values(maxima))
        total = sum(abs2.(component ./ scale) for component in values(components))
        centroid = trapz(result.z .* total, result.z) / trapz(total, result.z)
        oscillatory = filter(value -> abs(imag(value)) > 1e-7, spectrum_values)
        least = oscillatory[argmax(real.(oscillatory))]
        row["leading_temporal_real"] = real(lead)
        row["leading_temporal_imag"] = imag(lead)
        row["temporal_positive_real_count"] = count(>(1e-7), real.(spectrum_values))
        for (name, maximum_value) in maxima
            row["dominant_mode_$(name)_relative_max"] = maximum_value / scale
        end
        row["dominant_mode_z_centroid"] = centroid
        row["least_damped_oscillatory_real"] = real(least)
        row["least_damped_oscillatory_imag_abs"] = abs(imag(least))
        row["least_damped_oscillatory_period"] = 2pi / abs(imag(least))
        row["least_damped_oscillatory_quality_factor"] =
            abs(imag(least)) / (2abs(real(least)))
        row["similarity_temporal_class"] = real(lead) > 1e-7 ?
            (abs(imag(lead)) > 1e-7 ? "unstable-oscillatory" : "unstable-monotonic") :
            (any(abs.(imag.(spectrum_values)) .> 1e-7) ? "stable-with-damped-oscillations" : "stable-monotonic")
    end

    write_csv_rows(joinpath(output, "branch_diagnostics.csv"), diagnostic_rows)
    write_csv_rows(joinpath(output, "midplane_spatial_eigenvalues.csv"), spatial_rows)
    write_csv_rows(joinpath(output, "temporal_eigenvalues_convergence.csv"), temporal_rows)
    write_json(joinpath(output, "temporal_convergence.json"), convergence)
    plot_temporal_spectra(final_results, joinpath(output, "temporal_spectra.png"))
    plot_dominant_modes(final_results, joinpath(output, "dominant_temporal_modes.png"))
    plot_spatial_dynamics(profiles, joinpath(output, "spatial_phase_portraits.png"))

    if !args["skip-branch-scan"]
        rows = scan_temporal_growth_along_branch(args["re-h"], args["pr"],
                                                 args["branch-scan-order"])
        write_csv_rows(joinpath(output, "branch_temporal_growth.csv"), rows)
        plot_branch_growth(rows, joinpath(output, "branch_stability_exchange.png"))
    end
    JSON3.pretty(stdout, diagnostic_rows; allow_inf=true)
    println()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
