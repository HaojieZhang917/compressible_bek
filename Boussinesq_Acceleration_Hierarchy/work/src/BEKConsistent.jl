module BEKConsistent

using LinearAlgebra
using DelimitedFiles
using ..BEKIsothermal

export ConsistentSolution, solve_consistent_isothermal,
       ConsistentFold, solve_consistent_fixed_tw, solve_consistent_fixed_h,
       solve_consistent_fixed_h_at_ro,
       solve_consistent_pseudoarclength,
       solve_consistent_pseudoarclength_ro_full, fixed_tw_condition,
       solve_consistent_fold, continue_consistent_fold,
       write_solution

const DEFAULT_PR = 0.72

struct ConsistentSolution
    Ro::Float64
    Co::Float64
    gamma::Float64
    Pr::Float64
    operators::BEKOperators
    fields::Matrix{Float64} # rows: H, F, G, T
    Tw::Float64
    Hinf::Float64
    residual::Float64
    iterations::Int
end

struct ConsistentFold
    solution::ConsistentSolution
    nullvector::Vector{Float64}
    residual::Float64
    iterations::Int
end

pack(fields) = vec(permutedims(fields))
unpack(state, n) = permutedims(reshape(state, n, 4))

function _residual_jacobian(state, Ro, tw, gamma, Pr,
                            op::BEKOperators)
    abs(Ro) > 1e-12 || error(
        "the present differential-rotation form is non-uniform at Ro=0")
    n = length(op.x)
    H = view(state, 1:n)
    F = view(state, n+1:2n)
    G = view(state, 2n+1:3n)
    T = view(state, 3n+1:4n)
    D, D2 = op.deta, op.deta2
    Fp, Gp, Tp = D * F, D * G, D * T
    Co = 2 - Ro - Ro^2
    s = Co / 2
    f = s + Ro
    chi = 1 .- gamma .* (T .- 1)
    q = s .+ Ro .* G
    Ar = Ro .* (F.^2 + H .* Fp) .- q.^2 ./ Ro
    Atheta = Ro .* (2 .* F .* G + H .* Gp) .+ Co .* F
    identity = Matrix{Float64}(I, n, n)

    residual = vcat(
        D * H + 2F,
        D2 * F + chi .* Ar .+ f^2 / Ro,
        D2 * G + chi .* Atheta,
        D2 * T - Pr .* H .* Tp,
    )

    jacobian = zeros(4n, 4n)
    rH, rF, rG, rT = 1:n, n+1:2n, 2n+1:3n, 3n+1:4n
    jacobian[rH, rH] .= D
    jacobian[rH, rF] .= 2identity

    jacobian[rF, rH] .= Diagonal(chi .* Ro .* Fp)
    jacobian[rF, rF] .= D2 + Diagonal(chi .* 2Ro .* F) +
                         Diagonal(chi .* Ro .* H) * D
    jacobian[rF, rG] .= Diagonal(chi .* (-2 .* q))
    jacobian[rF, rT] .= Diagonal(-gamma .* Ar)

    jacobian[rG, rH] .= Diagonal(chi .* Ro .* Gp)
    jacobian[rG, rF] .= Diagonal(chi .* (2Ro .* G .+ Co))
    jacobian[rG, rG] .= D2 + Diagonal(chi .* 2Ro .* F) +
                         Diagonal(chi .* Ro .* H) * D
    jacobian[rG, rT] .= Diagonal(-gamma .* Atheta)

    jacobian[rT, rH] .= -Pr .* Diagonal(Tp)
    jacobian[rT, rT] .= D2 - Pr .* Diagonal(H) * D

    for (row, column, value) in
        ((first(rH), first(rH), H[1]),
         (first(rF), first(rF), F[1]),
         (first(rG), first(rG), G[1]),
         (last(rF), last(rF), F[end]),
         (last(rG), last(rG), G[end] - 1),
         (first(rT), first(rT), T[1] - tw),
         (last(rT), last(rT), T[end] - 1))
        residual[row] = value
        jacobian[row, :] .= 0
        jacobian[row, column] = 1
    end
    row = last(rH)
    residual[row] = dot(op.dx[end, :], H)
    jacobian[row, :] .= 0
    jacobian[row, 1:n] .= op.dx[end, :]
    residual, jacobian
end

function _tw_derivative(state, op::BEKOperators)
    n = length(op.x)
    derivative = zeros(4n)
    derivative[3n+1] = -1
    derivative
end

function _residual_ro_derivative(state, Ro, tw, gamma, Pr,
                                 op::BEKOperators)
    n = length(op.x)
    H = view(state,1:n)
    F = view(state,n+1:2n)
    G = view(state,2n+1:3n)
    T = view(state,3n+1:4n)
    D = op.deta
    Fp, Gp = D*F, D*G
    Co = 2-Ro-Ro^2
    dCo = -1-2Ro
    s, ds = Co/2, dCo/2
    f, df = s+Ro, ds+1
    chi = 1 .- gamma.*(T.-1)
    q = s .+ Ro.*G
    dq = ds .+ G
    core_r = F.^2 + H.*Fp
    core_t = 2 .* F .* G + H .* Gp
    dAr = core_r .- 2 .* q .* dq ./ Ro .+ q.^2 ./ Ro^2
    dAtheta = core_t .+ dCo.*F
    derivative = vcat(
        zeros(n),
        chi.*dAr .+ 2f*df/Ro .- f^2/Ro^2,
        chi.*dAtheta,
        zeros(n),
    )
    rH, rF, rG, rT = 1:n, n+1:2n, 2n+1:3n, 3n+1:4n
    for row in (first(rH),last(rH),first(rF),last(rF),
                first(rG),last(rG),first(rT),last(rT))
        derivative[row] = 0
    end
    derivative
end

"""Directional Hessian `d(J*v)/du` for the consistent residual."""
function _jacobian_directional_derivative(state, Ro, tw, gamma, Pr,
                                          op::BEKOperators, vector)
    n = length(op.x)
    H = view(state,1:n); F = view(state,n+1:2n)
    G = view(state,2n+1:3n); T = view(state,3n+1:4n)
    vH = view(vector,1:n); vF = view(vector,n+1:2n)
    vG = view(vector,2n+1:3n); vT = view(vector,3n+1:4n)
    D = op.deta
    Fp, Gp = D*F, D*G
    DvF, DvG, DvT = D*vF, D*vG, D*vT
    Co = 2-Ro-Ro^2
    s = Co/2
    chi = 1 .- gamma.*(T.-1)
    q = s .+ Ro.*G
    dArv = Ro.*(2 .* F .* vF + vH.*Fp + H.*DvF) .-
            2 .* q .* vG
    dAthetav = Ro.*(2 .* vF .* G + 2 .* F .* vG +
                    vH.*Gp + H.*DvG) .+ Co.*vF

    K = zeros(4n,4n)
    rH, rF, rG, rT = 1:n, n+1:2n, 2n+1:3n, 3n+1:4n

    K[rF,rH] .= Diagonal(chi.*Ro.*DvF) .-
                 gamma.*Diagonal(vT.*Ro.*Fp)
    K[rF,rF] .= Diagonal(chi.*2Ro.*vF) .+
                 Diagonal(chi.*Ro.*vH)*D .-
                 gamma.*Diagonal(vT) *
                 (Diagonal(2Ro.*F) + Diagonal(Ro.*H)*D)
    K[rF,rG] .= Diagonal(chi.*(-2Ro.*vG)) .+
                 gamma.*Diagonal(2 .* vT .* q)
    K[rF,rT] .= Diagonal(-gamma.*dArv)

    K[rG,rH] .= Diagonal(chi.*Ro.*DvG) .-
                 gamma.*Diagonal(vT.*Ro.*Gp)
    K[rG,rF] .= Diagonal(chi.*2Ro.*vG) .-
                 gamma.*Diagonal(vT.*(2Ro.*G .+ Co))
    K[rG,rG] .= Diagonal(chi.*2Ro.*vF) .+
                 Diagonal(chi.*Ro.*vH)*D .-
                 gamma.*Diagonal(vT) *
                 (Diagonal(2Ro.*F) + Diagonal(Ro.*H)*D)
    K[rG,rT] .= Diagonal(-gamma.*dAthetav)

    K[rT,rH] .= Diagonal(-Pr.*DvT)
    K[rT,rT] .= -Pr.*Diagonal(vH)*D

    for row in (first(rH),last(rH),first(rF),last(rF),
                first(rG),last(rG),first(rT),last(rT))
        K[row,:] .= 0
    end
    K
end

function _newton(system, initial; tolerance=1e-10, max_iterations=40)
    state = Float64.(initial)
    for iteration in 1:max_iterations
        residual, jacobian = system(state)
        norm0 = norm(residual, Inf)
        norm0 < tolerance && return state, norm0, iteration
        scales = 1 ./ max.([norm(view(jacobian, row, :))
                            for row in axes(jacobian, 1)], 1e-14)
        step = -((scales .* jacobian) \ (scales .* residual))
        damping = 1.0
        accepted = false
        while damping >= 2.0^-22
            trial = state .+ damping .* step
            trial_norm = norm(first(system(trial)), Inf)
            if isfinite(trial_norm) && trial_norm < norm0
                state = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        if !accepted
            norm0 <= 5tolerance && return state, norm0, iteration
            error("consistent-model Newton line search failed at residual=$norm0")
        end
    end
    error("consistent-model Newton iteration did not converge")
end

function solve_consistent_isothermal(Ro::Real=-1.0; gamma=1.0,
                                     Pr=DEFAULT_PR, degree=120,
                                     a=2.0, b=0.6, c=0.5,
                                     tolerance=1e-10)
    isothermal = BEKIsothermal.solve_isothermal(Ro; degree=degree,
        a=a, b=b, c=c, tolerance=tolerance)
    n = degree + 1
    fields = vcat(isothermal.fields, ones(1, n))
    state, residual, iterations = _newton(
        u -> _residual_jacobian(u, Float64(Ro), 1.0, Float64(gamma),
                                Float64(Pr), isothermal.operators),
        pack(fields); tolerance=tolerance)
    output = unpack(state, n)
    ConsistentSolution(Float64(Ro), 2-Float64(Ro)-Float64(Ro)^2,
        Float64(gamma), Float64(Pr), isothermal.operators, output, 1.0,
        output[1,end], residual, iterations)
end

function solve_consistent_fixed_tw(tw::Real, seed::ConsistentSolution;
                                   tolerance=1e-10)
    value = Float64(tw)
    state, residual, iterations = _newton(
        u -> _residual_jacobian(u, seed.Ro, value, seed.gamma, seed.Pr,
                                seed.operators),
        pack(seed.fields); tolerance=tolerance)
    fields = unpack(state, length(seed.operators.x))
    ConsistentSolution(seed.Ro, seed.Co, seed.gamma, seed.Pr, seed.operators,
        fields, value, fields[1,end], residual, iterations)
end

function solve_consistent_fixed_h(hinf::Real, seed::ConsistentSolution;
                                  tolerance=1e-10)
    op = seed.operators
    n = length(op.x)
    target = Float64(hinf)
    initial = vcat(pack(seed.fields), seed.Tw)
    function system(state)
        fields = view(state, 1:4n)
        tw = state[end]
        residual, jacobian = _residual_jacobian(
            fields, seed.Ro, tw, seed.gamma, seed.Pr, op)
        output = vcat(residual, fields[n] - target)
        full_jacobian = zeros(4n+1, 4n+1)
        full_jacobian[1:4n, 1:4n] .= jacobian
        full_jacobian[1:4n, end] .= _tw_derivative(fields, op)
        full_jacobian[end, n] = 1
        output, full_jacobian
    end
    state, residual, iterations = _newton(system, initial;
                                           tolerance=tolerance)
    fields = unpack(view(state, 1:4n), n)
    ConsistentSolution(seed.Ro, seed.Co, seed.gamma, seed.Pr, op, fields,
        state[end], target, residual, iterations)
end

"""Solve at prescribed `Ro` and `Hinf` from a neighbouring consistent state."""
function solve_consistent_fixed_h_at_ro(hinf::Real, Ro::Real,
                                        seed::ConsistentSolution;
                                        tolerance=1e-10)
    op = seed.operators
    n = length(op.x)
    target = Float64(hinf)
    rnew = Float64(Ro)
    initial = vcat(pack(seed.fields),seed.Tw)
    function system(state)
        fields = view(state,1:4n)
        tw = state[end]
        residual,jacobian = _residual_jacobian(
            fields,rnew,tw,seed.gamma,seed.Pr,op)
        output = vcat(residual,fields[n]-target)
        full = zeros(4n+1,4n+1)
        full[1:4n,1:4n] .= jacobian
        full[1:4n,end] .= _tw_derivative(fields,op)
        full[end,n] = 1
        output,full
    end
    state,residual,iterations = _newton(system,initial;
                                         tolerance=tolerance)
    fields = unpack(view(state,1:4n),n)
    ConsistentSolution(rnew,2-rnew-rnew^2,seed.gamma,seed.Pr,op,
        fields,state[end],target,residual,iterations)
end

"""Pseudo-arclength step using the complete state and wall temperature."""
function solve_consistent_pseudoarclength(previous::ConsistentSolution,
                                          current::ConsistentSolution;
                                          step=0.01, tolerance=1e-10)
    previous.Ro == current.Ro || error("Ro mismatch")
    previous.gamma == current.gamma || error("gamma mismatch")
    op = current.operators
    n = length(op.x)
    q0 = vcat(pack(previous.fields), previous.Tw)
    q1 = vcat(pack(current.fields), current.Tw)
    tangent = q1 - q0
    tangent ./= norm(tangent)
    prediction = q1 .+ Float64(step) .* tangent
    function system(state)
        fields = view(state, 1:4n)
        tw = state[end]
        residual, jacobian = _residual_jacobian(
            fields, current.Ro, tw, current.gamma, current.Pr, op)
        output = vcat(residual, dot(tangent, state - prediction))
        full_jacobian = zeros(4n+1, 4n+1)
        full_jacobian[1:4n, 1:4n] .= jacobian
        full_jacobian[1:4n, end] .= _tw_derivative(fields, op)
        full_jacobian[end, :] .= tangent
        output, full_jacobian
    end
    state, residual, iterations = _newton(system, prediction;
                                           tolerance=tolerance)
    fields = unpack(view(state, 1:4n), n)
    ConsistentSolution(current.Ro, current.Co, current.gamma, current.Pr, op,
        fields, state[end], fields[1,end], residual, iterations)
end

"""Full-state pseudo-arclength continuation at fixed `Hinf`, varying `Ro`."""
function solve_consistent_pseudoarclength_ro_full(
        previous::ConsistentSolution,current::ConsistentSolution;
        step=0.005,tolerance=1e-9)
    previous.operators.degree == current.operators.degree ||
        error("grid mismatch")
    abs(previous.Hinf-current.Hinf) < 1e-10 ||
        error("Hinf must remain fixed")
    previous.gamma == current.gamma || error("gamma mismatch")
    op = current.operators
    n = length(op.x)
    m = 4n
    qprevious = vcat(pack(previous.fields),previous.Tw,previous.Ro)
    qcurrent = vcat(pack(current.fields),current.Tw,current.Ro)
    difference = qcurrent-qprevious
    metric = vcat(fill(1/m,m),1.0,1.0)
    distance = sqrt(sum(metric.*difference.^2))
    distance > 1e-13 || error("continuation seeds coincide")
    direction = difference./distance
    tangent = metric.*direction
    prediction = qcurrent .+ Float64(step).*direction
    hfix = current.Hinf

    function system(state)
        fields = view(state,1:m)
        tw = state[m+1]
        ro = state[m+2]
        residual,jacobian = _residual_jacobian(fields,ro,tw,
            current.gamma,current.Pr,op)
        dro = _residual_ro_derivative(fields,ro,tw,current.gamma,
                                      current.Pr,op)
        output = vcat(residual,fields[n]-hfix,
                      dot(tangent,state-prediction))
        full = zeros(m+2,m+2)
        full[1:m,1:m] .= jacobian
        full[1:m,m+1] .= _tw_derivative(fields,op)
        full[1:m,m+2] .= dro
        full[m+1,n] = 1
        full[m+2,:] .= tangent
        output,full
    end
    state,residual,iterations = _newton(system,prediction;
        tolerance=tolerance,max_iterations=25)
    ro = state[m+2]
    fields = unpack(view(state,1:m),n)
    ConsistentSolution(ro,2-ro-ro^2,current.gamma,current.Pr,op,
        fields,state[m+1],hfix,residual,iterations)
end

function fixed_tw_condition(solution::ConsistentSolution)
    _, jacobian = _residual_jacobian(pack(solution.fields), solution.Ro,
        solution.Tw, solution.gamma, solution.Pr, solution.operators)
    row_norm = max.([norm(view(jacobian, row, :))
                     for row in axes(jacobian, 1)], 1e-14)
    singular = svdvals(jacobian ./ row_norm)
    (sigma_min=singular[end], sigma_max=singular[1],
     ratio=singular[end]/singular[1])
end

function _fold_system(state, Ro, gamma, Pr, op::BEKOperators)
    n = length(op.x)
    m = 4n
    fields = view(state,1:m)
    tw = state[m+1]
    vector = view(state,m+2:2m+1)
    residual,jacobian = _residual_jacobian(fields,Ro,tw,gamma,Pr,op)
    K = _jacobian_directional_derivative(fields,Ro,tw,gamma,Pr,op,
                                         vector)
    output = vcat(residual,jacobian*vector,dot(vector,vector)-1)
    full = zeros(2m+1,2m+1)
    full[1:m,1:m] .= jacobian
    full[1:m,m+1] .= _tw_derivative(fields,op)
    full[m+1:2m,1:m] .= K
    full[m+1:2m,m+2:2m+1] .= jacobian
    full[end,m+2:2m+1] .= 2 .* vector
    output,full
end

function _initial_nullvector(solution::ConsistentSolution)
    state = pack(solution.fields)
    _,jacobian = _residual_jacobian(state,solution.Ro,solution.Tw,
        solution.gamma,solution.Pr,solution.operators)
    rownorm = max.([norm(view(jacobian,row,:))
                    for row in axes(jacobian,1)],1e-14)
    decomposition = svd(jacobian ./ rownorm)
    vector = decomposition.V[:,end]
    vector ./ norm(vector)
end

"""Correct a nearby state to a simple fixed-`Tw` fold."""
function solve_consistent_fold(seed::ConsistentSolution;
                               null_seed=nothing,tolerance=1e-9)
    op = seed.operators
    n = length(op.x)
    vector = null_seed === nothing ? _initial_nullvector(seed) :
             Float64.(null_seed)
    vector ./= norm(vector)
    initial = vcat(pack(seed.fields),seed.Tw,vector)
    state,residual,iterations = _newton(
        z -> _fold_system(z,seed.Ro,seed.gamma,seed.Pr,op),initial;
        tolerance=tolerance,max_iterations=20)
    fields = unpack(view(state,1:4n),n)
    tw = state[4n+1]
    output = ConsistentSolution(seed.Ro,seed.Co,seed.gamma,seed.Pr,op,
        fields,tw,fields[1,end],residual,iterations)
    ConsistentFold(output,collect(view(state,4n+2:8n+1)),
                   residual,iterations)
end

function _fold_ro_column(fields,tw,vector,Ro,gamma,Pr,op)
    residual_ro = _residual_ro_derivative(fields,Ro,tw,gamma,Pr,op)
    # J_Ro*v is evaluated as the state directional derivative of R_Ro.
    # This avoids differencing the large collocation Jacobian itself.
    epsilon = sqrt(eps(Float64))*(1+norm(fields))/max(norm(vector),1e-14)
    plus = _residual_ro_derivative(fields.+epsilon.*vector,Ro,tw,
                                   gamma,Pr,op)
    minus = _residual_ro_derivative(fields.-epsilon.*vector,Ro,tw,
                                    gamma,Pr,op)
    vpart = (plus-minus)./(2epsilon)
    vcat(residual_ro,vpart,0.0)
end

"""One full-state pseudo-arclength step along the two-parameter fold locus."""
function continue_consistent_fold(previous::ConsistentFold,
                                  current::ConsistentFold;
                                  step=0.01,tolerance=1e-9)
    previous.solution.operators.degree == current.solution.operators.degree ||
        error("fold grid mismatch")
    op = current.solution.operators
    n = length(op.x)
    m = 4n
    prevv = previous.nullvector
    curv = copy(current.nullvector)
    dot(prevv,curv) < 0 && (curv .*= -1)
    qprev = vcat(pack(previous.solution.fields),previous.solution.Tw,
                 prevv,previous.solution.Ro)
    qcur = vcat(pack(current.solution.fields),current.solution.Tw,
                curv,current.solution.Ro)
    difference = qcur-qprev
    metric = vcat(fill(1/m,m),1.0,fill(1/m,m),1.0)
    distance = sqrt(sum(metric.*difference.^2))
    distance > 1e-13 || error("fold continuation seeds coincide")
    direction = difference./distance
    tangent = metric.*direction
    prediction = qcur .+ Float64(step).*direction
    gamma,Pr = current.solution.gamma,current.solution.Pr

    function system(state)
        fields = view(state,1:m)
        tw = state[m+1]
        vector = view(state,m+2:2m+1)
        ro = state[end]
        base,jacobian = _fold_system(view(state,1:2m+1),ro,
                                     gamma,Pr,op)
        rocolumn = _fold_ro_column(fields,tw,vector,ro,gamma,Pr,op)
        output = vcat(base,dot(tangent,state-prediction))
        full = zeros(2m+2,2m+2)
        full[1:2m+1,1:2m+1] .= jacobian
        full[1:2m+1,end] .= rocolumn
        full[end,:] .= tangent
        output,full
    end

    state,residual,iterations = _newton(system,prediction;
        tolerance=tolerance,max_iterations=18)
    ro = state[end]
    fields = unpack(view(state,1:m),n)
    solution = ConsistentSolution(ro,2-ro-ro^2,gamma,Pr,op,fields,
        state[m+1],fields[1,end],residual,iterations)
    ConsistentFold(solution,collect(view(state,m+2:2m+1)),
                   residual,iterations)
end

function write_solution(path::AbstractString, solution::ConsistentSolution)
    mkpath(dirname(path))
    data = hcat(solution.operators.eta, solution.fields')
    open(path, "w") do stream
        println(stream, "eta,H,F,G,T")
        writedlm(stream, data, ',')
    end
    path
end

end
