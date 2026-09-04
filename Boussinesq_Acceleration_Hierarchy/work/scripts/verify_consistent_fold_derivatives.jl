#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using Printf

function main()
    solution = solve_consistent_isothermal(-0.5; degree=40,
        a=2.0,b=0.6,c=0.5,gamma=1.0,Pr=0.72,tolerance=1e-10)
    solution = solve_consistent_fixed_tw(1.2,solution;tolerance=1e-10)
    state = BEKConsistent.pack(solution.fields)
    m = length(state)
    vector = normalize(sin.(collect(1.0:m)))
    direction = normalize(cos.(collect(1.0:m)))
    epsilon = 2e-6
    K = BEKConsistent._jacobian_directional_derivative(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators,vector)
    _,Jplus = BEKConsistent._residual_jacobian(
        state.+epsilon.*direction,solution.Ro,solution.Tw,
        solution.gamma,solution.Pr,solution.operators)
    _,Jminus = BEKConsistent._residual_jacobian(
        state.-epsilon.*direction,solution.Ro,solution.Tw,
        solution.gamma,solution.Pr,solution.operators)
    finite = ((Jplus-Jminus)*vector)./(2epsilon)
    analytic = K*direction
    hessian_error = norm(analytic-finite)/max(norm(finite),1e-14)

    roepsilon = 2e-6
    rplus,_ = BEKConsistent._residual_jacobian(state,
        solution.Ro+roepsilon,solution.Tw,solution.gamma,solution.Pr,
        solution.operators)
    rminus,_ = BEKConsistent._residual_jacobian(state,
        solution.Ro-roepsilon,solution.Tw,solution.gamma,solution.Pr,
        solution.operators)
    finite_ro = (rplus-rminus)./(2roepsilon)
    analytic_ro = BEKConsistent._residual_ro_derivative(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators)
    ro_error = norm(analytic_ro-finite_ro)/max(norm(finite_ro),1e-14)
    @printf("directional Hessian relative error = %.6e\n",hessian_error)
    @printf("Ro derivative relative error        = %.6e\n",ro_error)
    hessian_error < 2e-7 || error("directional Hessian verification failed")
    ro_error < 2e-7 || error("Ro derivative verification failed")
end

main()
