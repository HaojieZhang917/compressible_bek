Base.@kwdef struct RotorConfig
    re_h::Float64 = 1000.0
    pr::Float64 = 0.72
    tolerance::Float64 = 2e-9
    collocation_degree::Int = 100
    model::String = "traditional_centrifugal"
end

struct RotorProfile
    z::Vector{Float64}
    D::Matrix{Float64}
    weights::Vector{Float64}
    fields::Matrix{Float64} # rows: H, F, G, T
    pressure::Float64
    tw::Float64
    residual::Float64
    config::RotorConfig
end

struct RotorBranch
    columns::Dict{String,Vector{Float64}}
    profiles::Vector{RotorProfile}
end

function _rotor_reference_fields(config::RotorConfig, reference_npz::AbstractString,
                                 z::AbstractVector)
    archive = NPZ.npzread(reference_npz)
    zref = vec(Float64.(archive["x"]))
    fields = zeros(4, length(z))
    fields[1, :] .= linear_interpolate(zref, vec(Float64.(archive["w0"])), z)
    fields[2, :] .= linear_interpolate(zref, vec(Float64.(archive["u0"])), z)
    fields[3, :] .= linear_interpolate(zref, vec(Float64.(archive["v0"])), z)
    fields[4, :] .= 1.0
    return fields
end

function _rotor_system(vector::AbstractVector, config::RotorConfig,
                       fixed::Symbol, fixed_value::Float64,
                       D::AbstractMatrix, D2::AbstractMatrix)
    config.model in ("traditional_centrifugal", "soong_rotating_forces") ||
        error("Unknown rotor--stator model: $(config.model)")
    fields_vector = view(vector, 1:length(vector)-1)
    parameter = Float64(vector[end])
    pressure = fixed == :temperature ? parameter : fixed_value
    tw = fixed == :pressure ? parameter : fixed_value
    nodes = size(D, 1)
    fields = unpack_fields(fields_vector, nodes)
    H, F, G, T = (view(fields, row, :) for row in 1:4)
    DF, DG, DT = D * F, D * G, D * T
    theta = T .- 1
    sqrt_re = sqrt(config.re_h)
    re_h = config.re_h
    ranges = ntuple(i -> ((i - 1) * nodes + 1):(i * nodes), 4)
    rH, rF, rG, rT = ranges

    if config.model == "traditional_centrifugal"
        radial = theta
        azimuthal = G
    else
        radial = theta .* (2G .- 1)
        azimuthal = G .- theta
    end
    residual = vcat(
        D * H + 2sqrt_re .* F,
        D2 * F - sqrt_re .* H .* DF - re_h .* (F.^2 .- G.^2 .+ radial .+ pressure),
        D2 * G - sqrt_re .* H .* DG - 2re_h .* F .* azimuthal,
        D2 * T - config.pr * sqrt_re .* H .* DT,
    )
    jacobian = zeros(4nodes + 1, 4nodes + 1)
    identity = Matrix{Float64}(I, nodes, nodes)
    jacobian[rH, rH] .= D
    jacobian[rH, rF] .= 2sqrt_re .* identity
    jacobian[rF, rH] .= -sqrt_re .* Diagonal(DF)
    jacobian[rF, rF] .= D2 - sqrt_re .* Diagonal(H) * D - 2re_h .* Diagonal(F)
    jacobian[rG, rH] .= -sqrt_re .* Diagonal(DG)
    jacobian[rT, rH] .= -config.pr * sqrt_re .* Diagonal(DT)
    jacobian[rT, rT] .= D2 - config.pr * sqrt_re .* Diagonal(H) * D

    if config.model == "traditional_centrifugal"
        jacobian[rF, rG] .= 2re_h .* Diagonal(G)
        jacobian[rF, rT] .= -re_h .* identity
        jacobian[rG, rF] .= -2re_h .* Diagonal(G)
        jacobian[rG, rG] .= D2 - sqrt_re .* Diagonal(H) * D - 2re_h .* Diagonal(F)
    else
        jacobian[rF, rG] .= 2re_h .* Diagonal(G .- theta)
        jacobian[rF, rT] .= -re_h .* Diagonal(2G .- 1)
        jacobian[rG, rF] .= -2re_h .* Diagonal(G .- theta)
        jacobian[rG, rG] .= D2 - sqrt_re .* Diagonal(H) * D - 2re_h .* Diagonal(F)
        jacobian[rG, rT] .= 2re_h .* Diagonal(F)
    end

    full_residual = vcat(residual, H[end])
    _set_identity_row!(full_residual, jacobian, first(rH), first(rH), H[1])
    _set_identity_row!(full_residual, jacobian, first(rF), first(rF), F[1])
    _set_identity_row!(full_residual, jacobian, last(rF), last(rF), F[end])
    _set_identity_row!(full_residual, jacobian, first(rG), first(rG), G[1] - 1)
    _set_identity_row!(full_residual, jacobian, last(rG), last(rG), G[end])
    _set_identity_row!(full_residual, jacobian, first(rT), first(rT), T[1] - tw)
    _set_identity_row!(full_residual, jacobian, last(rT), last(rT), T[end] - 1)
    jacobian[end, :] .= 0.0
    jacobian[end, nodes] = 1.0

    if fixed == :temperature
        jacobian[rF, end] .= -re_h
        jacobian[first(rF), end] = 0.0
        jacobian[last(rF), end] = 0.0
    elseif fixed == :pressure
        jacobian[first(rT), end] = -1.0
    else
        error("fixed must be :temperature or :pressure")
    end
    return full_residual, jacobian
end

function _solve_rotor(config::RotorConfig, seed_fields::AbstractMatrix,
                      seed_parameter::Float64, fixed::Symbol, fixed_value::Float64)
    z, D, D2, weights = chebyshev_operators(config.collocation_degree, 0.0, 1.0)
    size(seed_fields, 2) == length(z) || error("Rotor seed grid mismatch")
    initial = vcat(pack_fields(seed_fields), seed_parameter)
    system = state -> _rotor_system(state, config, fixed, fixed_value, D, D2)
    state, history, _ = _newton_linesearch(system, initial;
                                            tolerance=config.tolerance,
                                            max_iterations=24)
    fields = unpack_fields(view(state, 1:length(state)-1), length(z))
    pressure = fixed == :temperature ? state[end] : fixed_value
    tw = fixed == :pressure ? state[end] : fixed_value
    return RotorProfile(z, D, weights, fields, pressure, tw, last(history), config)
end

function solve_rotor_isothermal(config::RotorConfig, reference_npz::AbstractString)
    z, D, D2, _ = chebyshev_operators(config.collocation_degree, 0.0, 1.0)
    fields = _rotor_reference_fields(config, reference_npz, z)
    H, F, G = (view(fields, row, :) for row in 1:3)
    estimate = (D2 * F - sqrt(config.re_h) .* H .* (D * F)) ./ config.re_h .-
               F.^2 .+ G.^2
    pressure_guess = median(estimate[3:end-2])
    return _solve_rotor(config, fields, pressure_guess, :temperature, 1.0)
end

function rotor_state(profile::RotorProfile, targets::AbstractVector)
    derivative_nodes = permutedims(profile.D * permutedims(profile.fields))
    fields = barycentric_interpolate(profile.z, profile.weights, profile.fields, targets)
    derivatives = barycentric_interpolate(profile.z, profile.weights,
                                           derivative_nodes, targets)
    output = zeros(7, length(targets))
    output[1, :] .= fields[1, :]
    output[2, :] .= fields[2, :]
    output[3, :] .= derivatives[2, :]
    output[4, :] .= fields[3, :]
    output[5, :] .= derivatives[3, :]
    output[6, :] .= fields[4, :]
    output[7, :] .= derivatives[4, :]
    return output
end

rotor_state(profile::RotorProfile, target::Real) = rotor_state(profile, [Float64(target)])[:, 1]

function _rotor_branch_columns(profiles::Vector{RotorProfile}; profile_points::Int=401)
    zplot = collect(range(0.0, 1.0; length=profile_points))
    states = [rotor_state(profile, zplot) for profile in profiles]
    columns = Dict{String,Vector{Float64}}(
        "pressure_gradient" => [p.pressure for p in profiles],
        "Tw" => [p.tw for p in profiles],
        "thermal_rossby" => [p.tw - 1 for p in profiles],
        "rotor_heat_flux_Tz" => [state[7, 1] for state in states],
        "G_mid" => [rotor_state(p, 0.5)[4] for p in profiles],
        "H_min" => [minimum(state[1, :]) for state in states],
        "F_min" => [minimum(state[2, :]) for state in states],
        "F_max" => [maximum(state[2, :]) for state in states],
        "nodes" => fill(Float64(profiles[1].config.collocation_degree + 1), length(profiles)),
        "max_rms_residual" => [p.residual for p in profiles],
    )
    return columns
end

function continue_rotor_pressure(pressure_values, initial::RotorProfile;
                                 profile_points::Int=401)
    profiles = RotorProfile[]
    seed = initial
    for pressure in pressure_values
        seed = _solve_rotor(seed.config, seed.fields, seed.tw, :pressure,
                            Float64(pressure))
        push!(profiles, seed)
    end
    return RotorBranch(_rotor_branch_columns(profiles; profile_points=profile_points), profiles)
end

function continue_rotor_temperature(temperature_values, initial::RotorProfile;
                                    profile_points::Int=401)
    profiles = RotorProfile[]
    seed = initial
    for temperature in temperature_values
        seed = _solve_rotor(seed.config, seed.fields, seed.pressure, :temperature,
                            Float64(temperature))
        push!(profiles, seed)
    end
    return RotorBranch(_rotor_branch_columns(profiles; profile_points=profile_points), profiles)
end

function solve_rotor_fixed_temperature(tw::Real, seed::RotorProfile)
    return _solve_rotor(seed.config, seed.fields, seed.pressure, :temperature,
                        Float64(tw))
end

function _interpolate_branch_column(branch::RotorBranch, name::String, pressure::Float64)
    x = branch.columns["pressure_gradient"]
    y = branch.columns[name]
    return linear_interpolate(x, y, [pressure])[1]
end

function rotor_turning_points(branch::RotorBranch)
    extrema = locate_quadratic_extrema(branch.columns["pressure_gradient"],
                                       branch.columns["Tw"])
    output = Dict{String,Any}[]
    for point in extrema
        row = Dict{String,Any}(
            "kind" => point.kind,
            "pressure_gradient" => point.x,
            "Tw" => point.y,
            "thermal_rossby" => point.y - 1,
        )
        for name in ("rotor_heat_flux_Tz", "G_mid", "H_min", "F_min", "F_max")
            row[name] = _interpolate_branch_column(branch, name, point.x)
        end
        push!(output, row)
    end
    return output
end

rotor_roots_at_temperature(branch::RotorBranch, target_tw::Real) =
    bracketed_roots(branch.columns["pressure_gradient"], branch.columns["Tw"], target_tw)

function pressure_grid(initial_pressure::Real, fine_step::Real)
    coarse_high = collect(Float64(initial_pressure):-0.001:0.015)
    fine = collect((coarse_high[end] - fine_step):-fine_step:(-0.006 - fine_step / 2))
    coarse_low = collect((fine[end] - 0.0005):-0.0005:-0.0501)
    return vcat(coarse_high, fine, coarse_low)
end

function rotor_validation(profile::RotorProfile, reference_npz::AbstractString)
    archive = NPZ.npzread(reference_npz)
    applicable = isapprox(profile.config.re_h, 1000.0)
    output = Dict{String,Any}(
        "reference_re_h" => 1000.0,
        "reference_comparison_applicable" => applicable,
        "pressure_gradient" => profile.pressure,
        "max_rms_residual" => profile.residual,
        "nodes" => length(profile.z),
    )
    if applicable
        zref = vec(Float64.(archive["x"]))
        computed = rotor_state(profile, zref)
        for (name, row, key) in (("H", 1, "w0"), ("F", 2, "u0"), ("G", 4, "v0"))
            delta = computed[row, :] .- vec(Float64.(archive[key]))
            output["$(name)_max_abs"] = maximum(abs, delta)
            output["$(name)_rms"] = sqrt(mean(abs2, delta))
        end
    end
    return output
end
