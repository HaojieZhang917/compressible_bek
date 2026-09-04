#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N=120; const HFIX=-0.02; const TOL=1e-8

function seed_at_h(Ro)
    s=solve_thermal_isothermal(Ro; degree=N, a=2.0, b=0.6, c=0.5, tolerance=TOL)
    h=s.Hinf; step=0.005
    while h < HFIX-1e-12
        hn=min(h+step,HFIX)
        try
            s=solve_thermal_fixed_h(hn,s;tolerance=TOL); h=hn
            step=min(0.01,1.2step)
        catch
            step*=0.5; step>=1e-5 || error("cannot construct initial tail seed")
        end
    end
    s
end

function main()
    outdir=joinpath(@__DIR__,"..","results","traditional_bek_solvability"); mkpath(outdir)
    state=seed_at_h(-0.5)
    open(joinpath(outdir,"long_tail_ro_continuation.csv"),"w") do io
        println(io,"Ro,Hinf,Tw,K,alpha_T,alpha_V_min,residual")
        ros=vcat(collect(-0.50:-0.01:-0.99),
                  collect(-0.991:-0.001:-0.999),
                  collect(-0.9991:-0.0001:-1.0))
        for Ro in ros
            state=solve_thermal_fixed_h_at_ro(HFIX,Ro,state;tolerance=TOL)
            rates=farfield_decay_rates(state); av=minimum(rates.velocity)
            @printf(io,"%.8f,%.12f,%.12f,%.12e,%.12e,%.12e,%.5e\n",
                    Ro,state.Hinf,state.Tw,rates.coupling,rates.thermal,av,state.residual)
            @printf("Ro=%+.2f Tw=%.9f alphaV=%.6f residual=%.2e\n",Ro,state.Tw,av,state.residual)
        end
    end
end

main()
