const VK_PR = 0.72
const VK_LAMBDA_C = 1.0

struct VKProfile
    z::Vector{Float64}
    D::Matrix{Float64}
    weights::Vector{Float64}
    fields::Matrix{Float64} # rows: H, F, G, T
    tw::Float64
    hinf::Float64
    residual::Float64
end

struct VKBranch
    hinf::Vector{Float64}
    tw::Vector{Float64}
    profiles::Vector{VKProfile}
end

pack_fields(fields::AbstractMatrix) = vcat((collect(view(fields, row, :)) for row in axes(fields, 1))...)

function unpack_fields(vector::AbstractVector, nodes::Int)
    return permutedims(reshape(vector, nodes, 4))
end

function _vk_initial_fields(z::AbstractVector)
    fields = zeros(4, length(z))
    fields[1, :] .= -0.8845 .* (1 .- exp.(-0.8 .* z))
    fields[2, :] .= 0.5102 .* z .* exp.(-0.8 .* z)
    fields[3, :] .= 1 .- exp.(-0.8 .* z)
    fields[4, :] .= 1.0
    return fields
end

function _set_identity_row!(residual, jacobian, row::Int, column::Int, value)
    residual[row] = value
    jacobian[row, :] .= 0.0
    jacobian[row, column] = 1.0
end

function _vk_fixed_tw_system(vector::AbstractVector, tw::Float64,
                             D::AbstractMatrix, D2::AbstractMatrix)
    nodes = size(D, 1)
    fields = unpack_fields(vector, nodes)
    H, F, G, T = (view(fields, row, :) for row in 1:4)
    DF, DG, DT = D * F, D * G, D * T
    ranges = ntuple(i -> ((i - 1) * nodes + 1):(i * nodes), 4)
    rH, rF, rG, rT = ranges
    residual = vcat(
        D * H + 2F,
        D2 * F - F.^2 - H .* DF + (G .- 1).^2 - VK_LAMBDA_C .* (T .- 1),
        D2 * G - 2F .* G - H .* DG + 2F,
        D2 * T - VK_PR .* H .* DT,
    )
    jacobian = zeros(4nodes, 4nodes)
    identity = Matrix{Float64}(I, nodes, nodes)
    jacobian[rH, rH] .= D
    jacobian[rH, rF] .= 2identity
    jacobian[rF, rH] .= -Diagonal(DF)
    jacobian[rF, rF] .= D2 - Diagonal(2F) - Diagonal(H) * D
    jacobian[rF, rG] .= Diagonal(2 .* (G .- 1))
    jacobian[rF, rT] .= -VK_LAMBDA_C .* identity
    jacobian[rG, rH] .= -Diagonal(DG)
    jacobian[rG, rF] .= Diagonal(2 .- 2G)
    jacobian[rG, rG] .= D2 - Diagonal(2F) - Diagonal(H) * D
    jacobian[rT, rH] .= -VK_PR .* Diagonal(DT)
    jacobian[rT, rT] .= D2 - VK_PR .* Diagonal(H) * D

    _set_identity_row!(residual, jacobian, first(rH), first(rH), H[1])
    _set_identity_row!(residual, jacobian, first(rF), first(rF), F[1])
    _set_identity_row!(residual, jacobian, last(rF), last(rF), F[end])
    _set_identity_row!(residual, jacobian, first(rG), first(rG), G[1])
    _set_identity_row!(residual, jacobian, last(rG), last(rG), G[end] - 1)
    _set_identity_row!(residual, jacobian, first(rT), first(rT), T[1] - tw)
    _set_identity_row!(residual, jacobian, last(rT), last(rT), T[end] - 1)
    return residual, jacobian
end

function _vk_fixed_h_system(vector::AbstractVector, hinf::Float64,
                            D::AbstractMatrix, D2::AbstractMatrix)
    fields_vector = view(vector, 1:length(vector)-1)
    tw = Float64(vector[end])
    residual, base_jacobian = _vk_fixed_tw_system(fields_vector, tw, D, D2)
    nodes = size(D, 1)
    jacobian = zeros(length(vector), length(vector))
    jacobian[1:end-1, 1:end-1] .= base_jacobian
    temperature_row = 3nodes + 1
    jacobian[temperature_row, end] = -1.0
    full_residual = vcat(residual, fields_vector[nodes] - hinf)
    jacobian[end, nodes] = 1.0
    return full_residual, jacobian
end

function _newton_linesearch(system, initial::AbstractVector;
                            tolerance::Float64=2e-10, max_iterations::Int=20)
    state = collect(Float64, initial)
    history = Float64[]
    for _ in 1:max_iterations
        residual, jacobian = system(state)
        norm0 = norm(residual, Inf)
        push!(history, norm0)
        norm0 < tolerance && return state, history, jacobian
        step = -(jacobian \ residual)
        damping = 1.0
        accepted = false
        while damping >= 2.0^-18
            trial = state .+ damping .* step
            trial_norm = norm(first(system(trial)), Inf)
            if isfinite(trial_norm) && trial_norm < norm0
                state = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        accepted || error("Newton line search failed at residual $norm0")
    end
    residual, jacobian = system(state)
    error("Newton iteration failed: residual=$(norm(residual, Inf))")
end

function solve_vk_isothermal(; zmax::Float64=20.0, degree::Int=100,
                             tolerance::Float64=2e-10)
    z, D, D2, weights = chebyshev_operators(degree, 0.0, zmax)
    initial = pack_fields(_vk_initial_fields(z))
    system = state -> _vk_fixed_tw_system(state, 1.0, D, D2)
    state, history, _ = _newton_linesearch(system, initial; tolerance=tolerance)
    fields = unpack_fields(state, length(z))
    return VKProfile(z, D, weights, fields, 1.0, fields[1, end], last(history))
end

function solve_vk_fixed_h(hinf::Real, seed::VKProfile;
                          degree::Int=length(seed.z)-1,
                          zmax::Float64=seed.z[end],
                          tolerance::Float64=2e-10)
    z, D, D2, weights = chebyshev_operators(degree, 0.0, zmax)
    initial_fields = barycentric_interpolate(seed.z, seed.weights, seed.fields, z)
    initial = vcat(pack_fields(initial_fields), seed.tw)
    system = state -> _vk_fixed_h_system(state, Float64(hinf), D, D2)
    state, history, _ = _newton_linesearch(system, initial; tolerance=tolerance)
    fields = unpack_fields(view(state, 1:length(state)-1), length(z))
    return VKProfile(z, D, weights, fields, state[end], Float64(hinf), last(history))
end

function vk_state(profile::VKProfile, targets::AbstractVector)
    derivatives = profile.D * permutedims(profile.fields)
    derivatives = permutedims(derivatives)
    fields = barycentric_interpolate(profile.z, profile.weights, profile.fields, targets)
    slopes = barycentric_interpolate(profile.z, profile.weights, derivatives, targets)
    output = zeros(7, length(targets))
    output[1, :] .= fields[1, :]
    output[2, :] .= slopes[2, :]
    output[3, :] .= fields[2, :]
    output[4, :] .= slopes[3, :]
    output[5, :] .= fields[3, :]
    output[6, :] .= slopes[4, :]
    output[7, :] .= fields[4, :]
    return output
end

vk_state(profile::VKProfile, target::Real) = vk_state(profile, [Float64(target)])[:, 1]

function trace_vk_h_branch(; zmax::Float64=20.0, degree::Int=100,
                           h_start::Float64=-0.75, h_stop::Float64=0.10,
                           dh::Float64=0.005, tolerance::Float64=2e-9)
    seed = solve_vk_isothermal(; zmax=zmax, degree=degree, tolerance=tolerance)
    profiles = VKProfile[]
    hs = collect(h_start:dh:(h_stop + dh / 4))
    for hinf in hs
        seed = solve_vk_fixed_h(hinf, seed; degree=degree, zmax=zmax,
                                tolerance=tolerance)
        push!(profiles, seed)
    end
    return VKBranch(hs, [profile.tw for profile in profiles], profiles)
end

function vk_turning_points(branch::VKBranch)
    return locate_quadratic_extrema(branch.hinf, branch.tw)
end

vk_roots_at_tw(branch::VKBranch, target_tw::Real) =
    bracketed_roots(branch.hinf, branch.tw, target_tw)

function vk_newton_diagnostics(profile::VKProfile)
    D2 = profile.D * profile.D
    residual, jacobian = _vk_fixed_tw_system(pack_fields(profile.fields), profile.tw,
                                             profile.D, D2)
    row_norm = [norm(view(jacobian, row, :)) for row in axes(jacobian, 1)]
    scaled = jacobian ./ max.(row_norm, eps(Float64))
    singular = svdvals(scaled)
    return (residual=norm(residual, Inf),
            sigma_ratio=singular[end] / singular[1])
end

function vk_similarity_eigenvalues(profile::VKProfile; degree::Int=60,
                                   zmax::Float64=profile.z[end])
    z, D, D2, weights = chebyshev_operators(degree, 0.0, zmax)
    fields = barycentric_interpolate(profile.z, profile.weights, profile.fields, z)
    H, F, G, T = (view(fields, row, :) for row in 1:4)
    Fp, Gp, Tp = D * F, D * G, D * T
    nodes = length(z)
    ranges = ntuple(i -> ((i - 1) * nodes + 1):(i * nodes), 4)
    rf, rg, rh, rt = ranges
    identity = Matrix{Float64}(I, nodes, nodes)
    A = zeros(4nodes, 4nodes)
    B = zeros(4nodes, 4nodes)
    A[rf, rf] .= D2 - Diagonal(H) * D - Diagonal(2F)
    A[rf, rg] .= Diagonal(2 .* (G .- 1))
    A[rf, rh] .= -Diagonal(Fp)
    A[rf, rt] .= -identity
    A[rg, rf] .= Diagonal(2 .* (1 .- G))
    A[rg, rg] .= D2 - Diagonal(H) * D - Diagonal(2F)
    A[rg, rh] .= -Diagonal(Gp)
    A[rh, rf] .= 2identity
    A[rh, rh] .= D
    A[rt, rh] .= -Diagonal(Tp)
    A[rt, rt] .= D2 / VK_PR - Diagonal(H) * D
    B[rf, rf] .= identity
    B[rg, rg] .= identity
    B[rt, rt] .= identity

    function set_bc(row, column)
        A[row, :] .= 0.0
        B[row, :] .= 0.0
        A[row, column] = 1.0
    end
    set_bc(first(rf), first(rf))
    set_bc(last(rf), last(rf))
    set_bc(first(rg), first(rg))
    set_bc(last(rg), last(rg))
    set_bc(first(rh), first(rh))
    set_bc(first(rt), first(rt))
    set_bc(last(rt), last(rt))
    values = eigen(A, B).values
    values = values[isfinite.(real.(values)) .& isfinite.(imag.(values)) .& (abs.(values) .< 1e5)]
    return values[sortperm(real.(values); rev=true)]
end

struct VKIVP
    z::Vector{Float64}
    state::Matrix{Float64}
end

function (solution::VKIVP)(targets::AbstractVector)
    output = zeros(size(solution.state, 1), length(targets))
    for row in axes(output, 1)
        output[row, :] .= linear_interpolate(solution.z, view(solution.state, row, :), targets)
    end
    return output
end

(solution::VKIVP)(target::Real) = solution([Float64(target)])[:, 1]

function vk_shoot(tw::Float64, slopes::AbstractVector; zmax::Float64=20.0,
                  step::Float64=0.002)
    function rhs(u)
        du = similar(u)
        H, Fp, F, Gp, G, Tp, T = u
        du[1] = -2F
        du[2] = F^2 + H * Fp - (G - 1)^2 + VK_LAMBDA_C * (T - 1)
        du[3] = Fp
        du[4] = 2F * G + H * Gp - 2F
        du[5] = Gp
        du[6] = VK_PR * H * Tp
        du[7] = Tp
        return du
    end
    function integrate(current)
        initial = [0.0, current[1], 0.0, current[2], 0.0, current[3], tw]
        intervals = ceil(Int, zmax / step)
        z = collect(range(0.0, zmax; length=intervals + 1))
        h = z[2] - z[1]
        state = zeros(7, length(z))
        state[:, 1] .= initial
        for index in 1:intervals
            u = view(state, :, index)
            k1 = rhs(u)
            k2 = rhs(u .+ 0.5h .* k1)
            k3 = rhs(u .+ 0.5h .* k2)
            k4 = rhs(u .+ h .* k3)
            state[:, index + 1] .= u .+ h .* (k1 .+ 2k2 .+ 2k3 .+ k4) ./ 6
        end
        return VKIVP(z, state)
    end
    function residual(current)
        solution = integrate(current)
        edge = solution.state[:, end]
        return [edge[3], edge[5] - 1, edge[7] - 1]
    end
    current = Float64.(slopes)
    converged = false
    for _ in 1:16
        value = residual(current)
        if norm(value, Inf) < 2e-9
            converged = true
            break
        end
        jacobian = zeros(3, 3)
        for column in 1:3
            delta = 2e-6 * max(abs(current[column]), 1.0)
            trial = copy(current)
            trial[column] += delta
            jacobian[:, column] .= (residual(trial) .- value) ./ delta
        end
        update = -(jacobian \ value)
        damping = 1.0
        while damping >= 2.0^-12
            trial = current .+ damping .* update
            if norm(residual(trial), Inf) < norm(value, Inf)
                current = trial
                break
            end
            damping *= 0.5
        end
    end
    solution = integrate(current)
    value = residual(current)
    answer = (zero=current, converged=converged || norm(value, Inf) < 2e-9)
    return answer, value, solution
end
