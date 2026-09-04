module BEKIsothermal

using LinearAlgebra
using DelimitedFiles

export BEKOperators, BEKSolution, bek_coefficients, solve_isothermal,
       residual_norm, write_solution

struct BEKOperators
    degree::Int
    a::Float64
    b::Float64
    c::Float64
    x::Vector{Float64}
    eta::Vector{Float64}
    dx::Matrix{Float64}
    deta::Matrix{Float64}
    deta2::Matrix{Float64}
end

struct BEKSolution
    Ro::Float64
    Co::Float64
    operators::BEKOperators
    fields::Matrix{Float64} # rows: H, F, G
    residual::Float64
    iterations::Int
end

function bek_coefficients(Ro::Real)
    r = Float64(Ro)
    return (Ro=r, Co=2 - r - r^2)
end

function _operators(degree::Int, a::Float64, b::Float64, c::Float64)
    degree >= 8 || error("degree must be at least 8")
    j = collect(0:degree)
    x = cos.(pi .* j ./ degree)
    weights = vcat(2.0, ones(degree - 1), 2.0) .* ((-1.0) .^ j)
    dx = (weights * transpose(1.0 ./ weights)) ./ (x .- transpose(x) + Matrix{Float64}(I, degree + 1, degree + 1))
    dx .-= Diagonal(vec(sum(dx; dims=2)))
    xi = -x
    q = b .* xi .+ (1 - b) .* (xi.^3 .+ c .* (1 .- xi.^2))
    eta = [k == degree + 1 ? Inf : a * (1 + q[k]) / (1 - q[k]) for k in eachindex(q)]
    dq = b .+ (1 - b) .* (3 .* xi.^2 .- 2c .* xi)
    dxi = -dx
    deta = Diagonal((1 .- q).^2 ./ (2a .* dq)) * dxi
    return BEKOperators(degree, a, b, c, x, eta, dx, deta, deta * deta)
end

function _initial_state(op::BEKOperators)
    z = op.eta
    H = zeros(length(z)); F = zeros(length(z)); G = zeros(length(z))
    finite = isfinite.(z)
    zz = z[finite]
    H[finite] .= -0.8845 .* (1 .- exp.(-0.8 .* zz))
    F[finite] .= 0.5102 .* zz .* exp.(-0.8 .* zz)
    G[finite] .= 1 .- exp.(-0.8 .* zz)
    H[end] = -0.5
    G[end] = 1.0
    return vcat(H, F, G)
end

function _residual_jacobian(state, Ro, op::BEKOperators)
    n = length(op.x); H = view(state, 1:n); F = view(state, n+1:2n); G = view(state, 2n+1:3n)
    D, D2 = op.deta, op.deta2
    Fp, Gp = D * F, D * G
    Co = 2 - Ro - Ro^2
    r = vcat(D * H + 2F,
             D2 * F + Ro .* (F.^2 + H .* Fp .- (G.^2 .- 1)) .- Co .* (G .- 1),
             D2 * G + Ro .* (2 .* F .* G + H .* Gp) .+ Co .* F)
    J = zeros(3n, 3n); I_n = Matrix{Float64}(I, n, n)
    rH, rF, rG = 1:n, n+1:2n, 2n+1:3n
    J[rH,rH] .= D; J[rH,rF] .= 2I_n
    J[rF,rH] .= Ro .* Diagonal(Fp); J[rF,rF] .= D2 + Ro .* Diagonal(2 .* F) ; J[rF,rF] .+= Ro .* Diagonal(H) * D
    J[rF,rG] .= -Ro .* Diagonal(2 .* G) - Co .* I_n
    J[rG,rF] .= Ro .* Diagonal(2 .* G) + Co .* I_n; J[rG,rH] .= Ro .* Diagonal(Gp); J[rG,rG] .= D2 + Ro .* Diagonal(H) * D
    # Boundary rows: H(0)=F(0)=G(0)=0, F(inf)=0, G(inf)=1.
    for (row, col, value) in ((first(rH), first(rH), H[1]), (first(rF), first(rF), F[1]),
                              (first(rG), first(rG), G[1]), (last(rF), last(rF), F[end]),
                              (last(rG), last(rG), G[end]-1))
        r[row] = value; J[row,:] .= 0; J[row,col] = 1
    end
    # Replace one redundant far-field continuity equation by H'(infinity)=0.
    row = last(rH); r[row] = dot(op.dx[end,:], H); J[row,:] .= 0; J[row,1:n] .= op.dx[end,:]
    return r, J
end

function _newton(state, Ro, op; tolerance=1e-10, max_iterations=30)
    u = Float64.(state)
    for it in 1:max_iterations
        r, J = _residual_jacobian(u, Ro, op); nr = norm(r, Inf)
        nr < tolerance && return u, nr, it
        scale = 1 ./ max.([norm(view(J,i,:)) for i in axes(J,1)], 1e-14)
        du = -((scale .* J) \ (scale .* r)); damping = 1.0
        accepted = false
        while damping >= 2.0^-20
            trial = u .+ damping .* du
            nt = norm(first(_residual_jacobian(trial, Ro, op)), Inf)
            if isfinite(nt) && nt < nr; u = trial; accepted = true; break end
            damping *= 0.5
        end
        accepted || error("BEK Newton line search failed: Ro=$Ro residual=$nr")
    end
    error("BEK Newton did not converge: Ro=$Ro")
end

function solve_isothermal(Ro::Real; degree=120, a=2.0, b=0.6, c=0.5, initial=nothing, tolerance=1e-10)
    op = _operators(degree, Float64(a), Float64(b), Float64(c))
    u0 = initial === nothing ? _initial_state(op) : initial
    u, r, it = _newton(u0, Float64(Ro), op; tolerance)
    n = degree + 1
    return BEKSolution(Float64(Ro), 2 - Float64(Ro) - Float64(Ro)^2, op,
                       permutedims(reshape(u, n, 3)), r, it)
end

residual_norm(sol::BEKSolution) = sol.residual

function write_solution(path::AbstractString, sol::BEKSolution)
    mkpath(dirname(path)); n = length(sol.operators.x)
    data = hcat(sol.operators.eta, sol.fields')
    open(path, "w") do io
        println(io, "# eta H F G; Ro=$(sol.Ro) Co=$(sol.Co) residual=$(sol.residual)")
        writedlm(io, data)
    end
    path
end

end
