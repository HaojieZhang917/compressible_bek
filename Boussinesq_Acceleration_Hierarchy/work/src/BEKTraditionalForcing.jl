module BEKTraditionalForcing

using LinearAlgebra
using ..BEKIsothermal

export ForcingSolution, lambda_cf, wall_temperature, physical_temperature,
       solve_forcing_isothermal, solve_forcing_fixed_h, solve_forcing_fixed_b,
       solve_forcing_pseudoarclength,
       fixed_b_condition, thermal_tail_length, farfield_decay_rates

const PR = 0.72

"""Traditional BEK solution using the regular forcing `B=Lambda_cf*(Tw-1)`.

The fourth field is the normalized response `Theta=(T-1)/(Tw-1)`, with
`Theta(0)=1` and `Theta(infinity)=0`.  This representation remains regular
when the wall-temperature conversion becomes singular at `Ro=0` or degenerate
at `Ro=1`.
"""
struct ForcingSolution
    Ro::Float64
    Co::Float64
    operators::BEKOperators
    fields::Matrix{Float64} # rows: H, F, G, Theta
    B::Float64
    Hinf::Float64
    residual::Float64
    iterations::Int
end

pack(fields) = vec(permutedims(fields))
unpack(state, n) = permutedims(reshape(state, n, 4))

function lambda_cf(Ro::Real)
    r = Float64(Ro)
    abs(r) > 1e-14 || throw(DomainError(r,
        "Lambda_cf=Co^2/(4Ro) is non-uniform at the Ekman limit Ro=0"))
    Co = 2 - r - r^2
    Co^2 / (4r)
end

function wall_temperature(B::Real, Ro::Real)
    lam = lambda_cf(Ro)
    if abs(lam) <= 1e-14
        abs(B) <= 1e-12 && return 1.0
        return copysign(Inf, Float64(B))
    end
    1 + Float64(B) / lam
end

function physical_temperature(solution::ForcingSolution)
    tw = wall_temperature(solution.B, solution.Ro)
    isfinite(tw) || error("physical temperature is undefined for B != 0 when Lambda_cf=0")
    1 .+ (tw - 1) .* view(solution.fields, 4, :)
end

function _residual_jacobian(state, Ro, B, op::BEKOperators)
    n = length(op.x)
    H = view(state, 1:n)
    F = view(state, n+1:2n)
    G = view(state, 2n+1:3n)
    Theta = view(state, 3n+1:4n)
    D, D2 = op.deta, op.deta2
    Fp, Gp, Thetap = D * F, D * G, D * Theta
    Co = 2 - Ro - Ro^2
    I_n = Matrix{Float64}(I, n, n)

    r = vcat(D * H + 2F,
             D2 * F + Ro .* (F.^2 + H .* Fp - (G.^2 .- 1)) .-
                 Co .* (G .- 1) + B .* Theta,
             D2 * G + Ro .* (2 .* F .* G + H .* Gp) + Co .* F,
             D2 * Theta - PR .* H .* Thetap)

    J = zeros(4n, 4n)
    rH, rF, rG, rTheta = 1:n, n+1:2n, 2n+1:3n, 3n+1:4n
    J[rH,rH] .= D
    J[rH,rF] .= 2I_n
    J[rF,rH] .= Ro .* Diagonal(Fp)
    J[rF,rF] .= D2 + Ro .* Diagonal(2 .* F) + Ro .* Diagonal(H) * D
    J[rF,rG] .= -Ro .* Diagonal(2 .* G) - Co .* I_n
    J[rF,rTheta] .= B .* I_n
    J[rG,rF] .= Ro .* Diagonal(2 .* G) + Co .* I_n
    J[rG,rH] .= Ro .* Diagonal(Gp)
    J[rG,rG] .= D2 + Ro .* Diagonal(H) * D
    J[rTheta,rH] .= -PR .* Diagonal(Thetap)
    J[rTheta,rTheta] .= D2 - PR .* Diagonal(H) * D

    for (row, col, value) in
        ((first(rH), first(rH), H[1]),
         (first(rF), first(rF), F[1]),
         (first(rG), first(rG), G[1]),
         (last(rF), last(rF), F[end]),
         (last(rG), last(rG), G[end] - 1),
         (first(rTheta), first(rTheta), Theta[1] - 1),
         (last(rTheta), last(rTheta), Theta[end]))
        r[row] = value
        J[row,:] .= 0
        J[row,col] = 1
    end
    row = last(rH)
    r[row] = dot(op.dx[end,:], H)
    J[row,:] .= 0
    J[row,1:n] .= op.dx[end,:]
    r, J
end

function _residual_b_derivative(state, op::BEKOperators)
    n = length(op.x)
    Theta = view(state, 3n+1:4n)
    d = vcat(zeros(n), collect(Theta), zeros(2n))
    d[first(n+1:2n)] = 0
    d[last(n+1:2n)] = 0
    d
end

function _newton(system, state; tolerance=1e-10, max_iterations=40)
    u = Float64.(state)
    for it in 1:max_iterations
        r, J = system(u)
        nr = norm(r, Inf)
        nr < tolerance && return u, nr, it
        scales = 1 ./ max.([norm(view(J,i,:)) for i in axes(J,1)], 1e-14)
        du = -((scales .* J) \ (scales .* r))
        damping = 1.0
        accepted = false
        while damping >= 2.0^-20
            trial = u .+ damping .* du
            nt = norm(first(system(trial)), Inf)
            if isfinite(nt) && nt < nr
                u = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        accepted || error("traditional-forcing Newton line search failed: residual=$nr")
    end
    error("traditional-forcing Newton did not converge")
end

function _initial_state(Ro, op::BEKOperators; tolerance=1e-10)
    iso = BEKIsothermal.solve_isothermal(Ro; degree=op.degree,
        a=op.a, b=op.b, c=op.c, tolerance=tolerance)
    z = op.eta
    Theta = zeros(length(z))
    finite = isfinite.(z)
    rate = max(0.2, -PR * iso.fields[1,end])
    Theta[finite] .= exp.(-rate .* z[finite])
    Theta[1] = 1
    Theta[end] = 0
    vcat(vec(iso.fields[1,:]), vec(iso.fields[2,:]),
         vec(iso.fields[3,:]), Theta)
end

function solve_forcing_isothermal(Ro::Real=-1.0; degree=120,
                                  a=2.0, b=0.6, c=0.5,
                                  tolerance=1e-10)
    op = BEKIsothermal._operators(degree, Float64(a), Float64(b), Float64(c))
    r = Float64(Ro)
    u0 = _initial_state(r, op; tolerance=max(tolerance, 1e-10))
    u, nr, it = _newton(s -> _residual_jacobian(s, r, 0.0, op), u0;
                         tolerance=tolerance)
    fields = unpack(u, degree + 1)
    ForcingSolution(r, 2-r-r^2, op, fields, 0.0, fields[1,end], nr, it)
end

function solve_forcing_fixed_h(hinf::Real, seed::ForcingSolution;
                               tolerance=1e-10)
    op = seed.operators
    n = length(op.x)
    r = seed.Ro
    hfix = Float64(hinf)
    u0 = vcat(pack(seed.fields), seed.B)
    function system(u)
        fields = view(u, 1:4n)
        B = u[end]
        residual, J = _residual_jacobian(fields, r, B, op)
        dB = _residual_b_derivative(fields, op)
        out = vcat(residual, fields[n] - hfix)
        JJ = zeros(4n+1, 4n+1)
        JJ[1:4n,1:4n] .= J
        JJ[1:4n,end] .= dB
        JJ[end,n] = 1
        out, JJ
    end
    u, nr, it = _newton(system, u0; tolerance=tolerance)
    fields = unpack(view(u,1:4n), n)
    ForcingSolution(r, 2-r-r^2, op, fields, u[end], hfix, nr, it)
end

function solve_forcing_fixed_b(B::Real, seed::ForcingSolution;
                               tolerance=1e-10)
    op = seed.operators
    n = length(op.x)
    bfix = Float64(B)
    u, nr, it = _newton(
        s -> _residual_jacobian(s, seed.Ro, bfix, op), pack(seed.fields);
        tolerance=tolerance)
    fields = unpack(u, n)
    ForcingSolution(seed.Ro, seed.Co, op, fields, bfix, fields[1,end], nr, it)
end

"""Pseudo-arclength step in the `(Hinf,B)` projection at fixed `Ro`."""
function solve_forcing_pseudoarclength(previous::ForcingSolution,
                                       current::ForcingSolution;
                                       step::Real=0.002,
                                       tolerance=1e-10)
    previous.Ro == current.Ro || error("pseudo-arclength seeds must have equal Ro")
    previous.operators.degree == current.operators.degree || error("grid mismatch")
    op = current.operators
    n = length(op.x)
    dH = current.Hinf-previous.Hinf
    dB = current.B-previous.B
    ds = hypot(dH,dB)
    ds > 1e-14 || error("pseudo-arclength seeds are indistinguishable")
    tH, tB = dH/ds, dB/ds
    scale = Float64(step)/ds
    qprevious = vcat(pack(previous.fields), previous.B)
    qcurrent = vcat(pack(current.fields), current.B)
    prediction = qcurrent + scale .* (qcurrent-qprevious)
    hpred = current.Hinf + scale*dH
    bpred = current.B + scale*dB
    function system(u)
        fields = view(u,1:4n)
        B = u[end]
        residual, J = _residual_jacobian(fields,current.Ro,B,op)
        dBd = _residual_b_derivative(fields,op)
        out = vcat(residual, tH*(fields[n]-hpred)+tB*(B-bpred))
        JJ = zeros(4n+1,4n+1)
        JJ[1:4n,1:4n] .= J
        JJ[1:4n,end] .= dBd
        JJ[end,n] = tH
        JJ[end,end] = tB
        out, JJ
    end
    u, nr, it = _newton(system,prediction;tolerance=tolerance)
    fields = unpack(view(u,1:4n),n)
    ForcingSolution(current.Ro,current.Co,op,fields,u[end],fields[1,end],nr,it)
end

"""Smallest singular value of the row-scaled fixed-`B` steady Jacobian."""
function fixed_b_condition(solution::ForcingSolution)
    _, J = _residual_jacobian(pack(solution.fields), solution.Ro,
                              solution.B, solution.operators)
    rownorm = max.([norm(view(J,i,:)) for i in axes(J,1)], 1e-14)
    singular = svdvals(J ./ rownorm)
    (sigma_min=singular[end], sigma_max=singular[1],
     ratio=singular[end]/singular[1])
end

thermal_tail_length(solution::ForcingSolution) =
    solution.Hinf < 0 ? -1 / (PR * solution.Hinf) : Inf

function farfield_decay_rates(Ro::Real, Hinf::Real)
    r = Float64(Ro)
    h = Float64(Hinf)
    Co = 2-r-r^2
    coupling = 2r + Co
    q = r*h
    roots = ComplexF64[]
    for sign in (-1.0, 1.0)
        root = sqrt(complex(q^2, 4sign*coupling))
        push!(roots, (q+root)/2, (q-root)/2)
    end
    decaying = sort([real(z) for z in roots if real(z) > 1e-12])
    (thermal=-PR*h, velocity=decaying, coupling=coupling, roots=roots)
end

farfield_decay_rates(solution::ForcingSolution) =
    farfield_decay_rates(solution.Ro, solution.Hinf)

end
