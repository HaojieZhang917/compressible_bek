module BEKThermal

using LinearAlgebra
using ..BEKIsothermal

export ThermalSolution, solve_thermal_isothermal, solve_thermal_fixed_h,
       solve_thermal_fixed_h_at_ro, solve_thermal_fixed_tw,
       solve_thermal_pseudoarclength, solve_thermal_pseudoarclength_ro,
       solve_thermal_pseudoarclength_ro_full,
       thermal_fixed_tw_condition, thermal_tail_length,
       farfield_decay_rates

const PR = 0.72

struct ThermalSolution
    Ro::Float64
    Co::Float64
    operators::BEKOperators
    fields::Matrix{Float64} # rows: H, F, G, T
    Tw::Float64
    Hinf::Float64
    residual::Float64
    iterations::Int
end

pack(fields) = vec(permutedims(fields))
unpack(state,n) = permutedims(reshape(state,n,4))

function _initial_state(op::BEKOperators; tw=1.0)
    n=length(op.x); z=op.eta; finite=isfinite.(z); zz=z[finite]
    H=zeros(n); F=zeros(n); G=zeros(n); T=ones(n)
    H[finite].=-0.8845 .* (1 .- exp.(-0.8 .* zz)); F[finite].=0.5102 .* zz .* exp.(-0.8 .* zz); G[finite].=1 .- exp.(-0.8 .* zz)
    H[end]=-0.8845; G[end]=1.0; T[1]=tw
    pack(vcat(H',F',G',T'))
end

function _residual_jacobian(state, Ro, tw, op::BEKOperators)
    n=length(op.x); H=view(state,1:n); F=view(state,n+1:2n); G=view(state,2n+1:3n); T=view(state,3n+1:4n)
    D,D2=op.deta,op.deta2; Fp,Gp,Tp=D*F,D*G,D*T; Co=2-Ro-Ro^2
    abs(Ro) > 1e-12 || error("traditional centrifugal thermal closure is singular at Ro=0; use a separately scaled Ekman-limit model")
    thermal_coeff = Co^2/(4Ro)
    I_n=Matrix{Float64}(I,n,n)
    r=vcat(D*H+2F,
           D2*F + Ro.*(F.^2 + H.*Fp - (G.^2 .- 1)) - Co.*(G .- 1) + thermal_coeff .* (T .- 1),
           D2*G + Ro.*(2 .* F .* G + H .* Gp) + Co .* F,
           D2*T - PR .* H .* Tp)
    J=zeros(4n,4n); rH=1:n; rF=n+1:2n; rG=2n+1:3n; rT=3n+1:4n
    J[rH,rH].=D; J[rH,rF].=2I_n
    J[rF,rH].=Ro.*Diagonal(Fp); J[rF,rF].=D2+Ro.*Diagonal(2 .* F)+Ro.*Diagonal(H)*D; J[rF,rG].=-Ro.*Diagonal(2 .* G)-Co.*I_n; J[rF,rT].=thermal_coeff .* I_n
    J[rG,rF].=Ro.*Diagonal(2 .* G)+Co.*I_n; J[rG,rH].=Ro.*Diagonal(Gp); J[rG,rG].=D2+Ro.*Diagonal(H)*D
    J[rT,rH].=-PR.*Diagonal(Tp); J[rT,rT].=D2-PR.*Diagonal(H)*D
    for (row,col,val) in ((first(rH),first(rH),H[1]),(first(rF),first(rF),F[1]),(first(rG),first(rG),G[1]),(last(rF),last(rF),F[end]),(last(rG),last(rG),G[end]-1),(first(rT),first(rT),T[1]-tw),(last(rT),last(rT),T[end]-1))
        r[row]=val; J[row,:].=0; J[row,col]=1
    end
    row=last(rH); r[row]=dot(op.dx[end,:],H); J[row,:].=0; J[row,1:n].=op.dx[end,:]
    r,J
end

function _newton(system, state; tolerance=2e-10, max_iterations=35)
    u=Float64.(state)
    for it in 1:max_iterations
        r,J=system(u); nr=norm(r,Inf); nr<tolerance && return u,nr,it
        scales=1 ./ max.([norm(view(J,i,:)) for i in axes(J,1)],1e-14); du=-((scales .* J)\(scales .* r)); damp=1.0; accepted=false
        while damp >= 2.0^-20
            trial=u.+damp.*du; nt=norm(first(system(trial)),Inf)
            if isfinite(nt) && nt<nr; u=trial; accepted=true; break end
            damp*=0.5
        end
        accepted || error("thermal BEK Newton line search failed: residual=$nr")
    end
    error("thermal BEK Newton did not converge")
end

function solve_thermal_isothermal(Ro::Real=-1.0; degree=120, a=2.0,b=0.6,c=0.5,tolerance=1e-9)
    op=BEKIsothermal._operators(degree,Float64(a),Float64(b),Float64(c)); u0=_initial_state(op)
    u,r,it=_newton(s->_residual_jacobian(s,Float64(Ro),1.0,op),u0;tolerance)
    f=unpack(u,degree+1); ThermalSolution(Float64(Ro),2-Float64(Ro)-Float64(Ro)^2,op,f,1.0,f[1,end],r,it)
end

function solve_thermal_fixed_h(hinf::Real, seed::ThermalSolution; tolerance=1e-9)
    op=seed.operators; n=length(op.x); u0=vcat(pack(seed.fields),seed.Tw); Ro=seed.Ro
    function system(u)
        fields= view(u,1:4n); tw=u[end]; r,J=_residual_jacobian(fields,Ro,tw,op); out=vcat(r,fields[n]-Float64(hinf)); JJ=zeros(4n+1,4n+1); JJ[1:4n,1:4n].=J; JJ[3n+1,end]=-1; JJ[end,n]=1; out,JJ
    end
    u,r,it=_newton(system,u0;tolerance); f=unpack(view(u,1:4n),n)
    ThermalSolution(Ro,2-Ro-Ro^2,op,f,u[end],Float64(hinf),r,it)
end

"""Solve at prescribed `Ro` and `Hinf`, using a neighbouring solution as seed."""
function solve_thermal_fixed_h_at_ro(hinf::Real, Ro::Real,
                                     seed::ThermalSolution; tolerance=1e-9)
    op=seed.operators; n=length(op.x); rnew=Float64(Ro)
    u0=vcat(pack(seed.fields),seed.Tw)
    function system(u)
        fields=view(u,1:4n); tw=u[end]
        r,J=_residual_jacobian(fields,rnew,tw,op)
        out=vcat(r,fields[n]-Float64(hinf))
        JJ=zeros(4n+1,4n+1)
        JJ[1:4n,1:4n].=J; JJ[3n+1,end]=-1; JJ[end,n]=1
        out,JJ
    end
    u,r,it=_newton(system,u0;tolerance)
    f=unpack(view(u,1:4n),n)
    ThermalSolution(rnew,2-rnew-rnew^2,op,f,u[end],Float64(hinf),r,it)
end

function _residual_ro_derivative(state, Ro, tw, op::BEKOperators)
    n=length(op.x); H=view(state,1:n); F=view(state,n+1:2n); G=view(state,2n+1:3n); T=view(state,3n+1:4n)
    D=op.deta; Co=2-Ro-Ro^2; dCo=-1-2Ro
    abs(Ro)>1e-12 || error("traditional centrifugal thermal closure is singular at Ro=0")
    lam=Co^2/(4Ro); dlam=Co*dCo/(2Ro)-Co^2/(4Ro^2)
    d=vcat(zeros(n), F.^2 + H.*(D*F) .- (G.^2 .- 1) .- dCo.*(G .- 1) .+ dlam.*(T .- 1),
            2 .* F .* G + H .* (D*G) .+ dCo .* F, zeros(n))
    d[first(1:n)]=0; d[first(n+1:2n)]=0; d[last(n+1:2n)]=0
    d[first(2n+1:3n)]=0; d[last(2n+1:3n)]=0
    d[last(1:n)]=0
    d
end

"""Pseudo-arclength continuation with fixed `Hinf` and varying `Ro`.

This is used to cross a turning point of the long-tail branch in the
`(Ro,Tw)` projection; it is not a replacement for a bordered cusp solve.
"""
function solve_thermal_pseudoarclength_ro(previous::ThermalSolution,
                                           current::ThermalSolution;
                                           step::Real=0.002, tolerance=1e-9)
    previous.operators.degree == current.operators.degree || error("grid mismatch")
    abs(previous.Hinf-current.Hinf)<1e-8 || error("Hinf must be fixed")
    op=current.operators; n=length(op.x); hfix=current.Hinf
    dRo=current.Ro-previous.Ro; dT=current.Tw-previous.Tw
    ds=hypot(dRo,dT); ds>1e-13 || error("continuation seeds are indistinguishable")
    tRo,tT=dRo/ds,dT/ds; fac=Float64(step)/ds
    ropred=current.Ro+fac*dRo; tpred=current.Tw+fac*dT
    prediction=vcat(pack(current.fields),current.Tw,current.Ro) .+
               fac .* (vcat(pack(current.fields),current.Tw,current.Ro)-
                       vcat(pack(previous.fields),previous.Tw,previous.Ro))
    function system(u)
        fields=view(u,1:4n); tw=u[4n+1]; ro=u[4n+2]
        r,J=_residual_jacobian(fields,ro,tw,op); dro=_residual_ro_derivative(fields,ro,tw,op)
        out=vcat(r,fields[n]-hfix,tRo*(ro-ropred)+tT*(tw-tpred))
        JJ=zeros(4n+2,4n+2); JJ[1:4n,1:4n].=J; JJ[1:4n,4n+2].=dro
        JJ[3n+1,4n+1]=-1; JJ[4n+1,n]=1
        JJ[4n+2,4n+1]=tT; JJ[4n+2,4n+2]=tRo
        out,JJ
    end
    u,r,it=_newton(system,prediction;tolerance)
    f=unpack(view(u,1:4n),n); ro=u[4n+2]
    ThermalSolution(ro,2-ro-ro^2,op,f,u[4n+1],hfix,r,it)
end

"""Full-state weighted pseudo-arclength step at fixed `Hinf`.

The continuation tangent contains every collocation field as well as `Tw`
and `Ro`.  Field components are averaged in the tangent norm so that the
large collocation block does not overwhelm the two scalar parameters.  This
is useful when the `(Ro,Tw)` projection develops a turning point; it remains
a continuation diagnostic rather than a bordered fold/cusp solve.
"""
function solve_thermal_pseudoarclength_ro_full(previous::ThermalSolution,
                                                current::ThermalSolution;
                                                step::Real=0.002,
                                                tolerance::Real=1e-9)
    previous.operators.degree == current.operators.degree || error("grid mismatch")
    abs(previous.Hinf-current.Hinf)<1e-8 || error("Hinf must be fixed")
    op=current.operators; n=length(op.x); hfix=current.Hinf
    qprev=vcat(pack(previous.fields),previous.Tw,previous.Ro)
    qcur=vcat(pack(current.fields),current.Tw,current.Ro)
    dq=qcur-qprev
    # Average the field block in the metric, while retaining unit weights for
    # the scalar continuation coordinates.
    metric=vcat(fill(1.0/(4n),4n),1.0,1.0)
    ds=sqrt(sum(metric .* dq.^2)); ds>1e-13 || error("continuation seeds are indistinguishable")
    tangent=metric .* dq ./ ds
    prediction=qcur + Float64(step) .* dq ./ ds

    function system(u)
        fields=view(u,1:4n); tw=u[4n+1]; ro=u[4n+2]
        r,J=_residual_jacobian(fields,ro,tw,op)
        dro=_residual_ro_derivative(fields,ro,tw,op)
        out=vcat(r,fields[n]-hfix,sum(tangent .* (u-prediction)))
        JJ=zeros(4n+2,4n+2)
        JJ[1:4n,1:4n].=J; JJ[1:4n,4n+2].=dro
        JJ[3n+1,4n+1]=-1; JJ[4n+1,n]=1; JJ[4n+2,:].=tangent
        out,JJ
    end
    u,r,it=_newton(system,prediction;tolerance)
    f=unpack(view(u,1:4n),n); ro=u[4n+2]
    ThermalSolution(ro,2-ro-ro^2,op,f,u[4n+1],hfix,r,it)
end

"""Solve the steady BEK equations at a prescribed wall temperature."""
function solve_thermal_fixed_tw(tw::Real, seed::ThermalSolution; tolerance=1e-9)
    op=seed.operators; n=length(op.x); Ro=seed.Ro
    u,r,it=_newton(s->_residual_jacobian(s,Ro,Float64(tw),op),pack(seed.fields); tolerance)
    f=unpack(u,n)
    ThermalSolution(Ro,2-Ro-Ro^2,op,f,Float64(tw),f[1,end],r,it)
end

"""
Take one pseudo-arclength step from `previous` through `current`.

The arclength constraint is imposed in the physically transparent
`(Hinf,Tw)` projection. This crosses ordinary `Tw` folds while retaining the
full collocation residual as the primary system.
"""
function solve_thermal_pseudoarclength(previous::ThermalSolution,
                                       current::ThermalSolution;
                                       step::Real=0.002, tolerance=1e-9)
    previous.Ro == current.Ro || error("pseudo-arclength seeds must have equal Ro")
    previous.operators.degree == current.operators.degree || error("pseudo-arclength grid mismatch")
    op=current.operators; n=length(op.x); Ro=current.Ro
    dH=current.Hinf-previous.Hinf; dT=current.Tw-previous.Tw
    ds=hypot(dH,dT); ds > 1e-13 || error("pseudo-arclength seeds are indistinguishable")
    tH,tT=dH/ds,dT/ds; scale=Float64(step)/ds
    prediction=vcat(pack(current.fields),current.Tw) .+
               scale .* (vcat(pack(current.fields),current.Tw) .-
                         vcat(pack(previous.fields),previous.Tw))
    hpred=current.Hinf+scale*dH; tpred=current.Tw+scale*dT
    function system(u)
        fields=view(u,1:4n); tw=u[end]
        r,J=_residual_jacobian(fields,Ro,tw,op)
        out=vcat(r,tH*(fields[n]-hpred)+tT*(tw-tpred))
        JJ=zeros(4n+1,4n+1)
        JJ[1:4n,1:4n].=J
        JJ[3n+1,end]=-1
        JJ[end,n]=tH; JJ[end,end]=tT
        out,JJ
    end
    u,r,it=_newton(system,prediction;tolerance)
    f=unpack(view(u,1:4n),n)
    ThermalSolution(Ro,2-Ro-Ro^2,op,f,u[end],f[1,end],r,it)
end

"""Smallest singular value of the row-scaled fixed-`Tw` Jacobian."""
function thermal_fixed_tw_condition(solution::ThermalSolution)
    r,J=_residual_jacobian(pack(solution.fields),solution.Ro,solution.Tw,solution.operators)
    rownorm=max.([norm(view(J,i,:)) for i in axes(J,1)],1e-14)
    singular=svdvals(J ./ rownorm)
    (sigma_min=singular[end], sigma_max=singular[1], ratio=singular[end]/singular[1])
end

"""Far-field thermal e-folding length from `T''-Pr*Hinf*T'=0`."""
thermal_tail_length(solution::ThermalSolution) =
    solution.Hinf < 0 ? -1/(PR*solution.Hinf) : Inf

"""Linear far-field decay rates for the thermal and coupled velocity modes."""
function farfield_decay_rates(Ro::Real, Hinf::Real)
    r=Float64(Ro); h=Float64(Hinf); Co=2-r-r^2
    coupling=2r+Co # = (2-r)(1+r)
    q=r*h
    roots=ComplexF64[]
    for sign in (-1.0,1.0)
        disc=complex(q^2,4sign*coupling)
        root=sqrt(disc)
        push!(roots,(q+root)/2,(q-root)/2)
    end
    decaying=sort([real(z) for z in roots if real(z)>1e-12])
    (thermal=-PR*h, velocity=decaying, coupling=coupling, roots=roots)
end

farfield_decay_rates(solution::ThermalSolution) =
    farfield_decay_rates(solution.Ro,solution.Hinf)

end
