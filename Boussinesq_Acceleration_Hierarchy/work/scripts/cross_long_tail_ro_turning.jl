#!/usr/bin/env julia
include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKThermal.jl"))
using .BEKIsothermal,.BEKThermal
using Printf

const N=120; const HFIX=-0.02; const TOL=1e-8
function seed(ro)
    s=solve_thermal_isothermal(ro;degree=N,tolerance=TOL); h=s.Hinf; step=.005
    while h<HFIX; hn=min(h+step,HFIX); s=solve_thermal_fixed_h(hn,s;tolerance=TOL); h=hn; step=min(.01,1.2step); end
    s
end
function main()
    out=joinpath(@__DIR__,"..","results","traditional_bek_solvability"); mkpath(out)
    p=seed(-.49); c=solve_thermal_fixed_h_at_ro(HFIX,-.50,p;tolerance=TOL)
    open(joinpath(out,"long_tail_ro_turning.csv"),"w") do io
        println(io,"Ro,Hinf,Tw,alpha_T,alpha_V_min,residual")
        for k in 1:180
            rates=farfield_decay_rates(c); av=minimum(rates.velocity)
            @printf(io,"%.10f,%.12f,%.12f,%.12e,%.12e,%.5e\n",c.Ro,c.Hinf,c.Tw,rates.thermal,av,c.residual)
            @printf("k=%d Ro=%+.8f Tw=%.9f H=%+.5f\n",k,c.Ro,c.Tw,c.Hinf)
            try
                n=solve_thermal_pseudoarclength_ro(p,c;step=.005,tolerance=TOL)
                p,c=c,n
            catch e
                println("stop: ",sprint(showerror,e)); break
            end
        end
    end
end
main()
