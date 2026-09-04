#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKThermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKTraditionalForcing.jl"))
using .BEKIsothermal
using .BEKThermal
using .BEKTraditionalForcing
using Printf

function option(prefix,default)
    for arg in ARGS
        startswith(arg,prefix) && return split(arg,"=";limit=2)[2]
    end
    default
end

const N=parse(Int,option("--degree=","120"))
const A=parse(Float64,option("--a=","2.0"))
const MAP_B=parse(Float64,option("--b=","0.6"))
const C=parse(Float64,option("--c=","0.5"))
const TARGET_RO=parse(Float64,option("--ro=","-0.99"))
const TOL=parse(Float64,option("--tol=","1e-8"))
const HSEED=-0.02

function high_tail_seed()
    if TARGET_RO >= -0.98
        s=BEKThermal.solve_thermal_isothermal(TARGET_RO;degree=N,a=A,b=MAP_B,c=C,tolerance=TOL)
        h=s.Hinf
        step=0.005
        while h < HSEED-1e-13
            hn=min(h+step,HSEED)
            try
                s=BEKThermal.solve_thermal_fixed_h(hn,s;tolerance=TOL)
                h=hn
                step=min(0.01,1.2step)
            catch
                step*=0.5
                step>=1e-5 || error("failed to reach direct high-tail seed")
            end
        end
        theta=(s.fields[4,:] .- 1) ./ (s.Tw-1)
        fields=vcat(s.fields[1:3,:],permutedims(theta))
        bforce=BEKTraditionalForcing.lambda_cf(s.Ro)*(s.Tw-1)
        return BEKTraditionalForcing.ForcingSolution(s.Ro,s.Co,s.operators,fields,
                   bforce,s.Hinf,s.residual,s.iterations)
    end
    s=BEKThermal.solve_thermal_isothermal(-0.5;degree=N,a=A,b=MAP_B,c=C,tolerance=TOL)
    h=s.Hinf
    step=0.005
    while h < HSEED-1e-13
        hn=min(h+step,HSEED)
        try
            s=BEKThermal.solve_thermal_fixed_h(hn,s;tolerance=TOL)
            h=hn
            step=min(0.01,1.2step)
        catch
            step*=0.5
            step>=1e-5 || error("failed to reach high-tail seed at Ro=-0.5")
        end
    end
    ro=-0.5
    while ro > TARGET_RO+1e-13
        rn=max(TARGET_RO,ro-0.01)
        s=BEKThermal.solve_thermal_fixed_h_at_ro(HSEED,rn,s;tolerance=TOL)
        ro=rn
    end
    theta=(s.fields[4,:] .- 1) ./ (s.Tw-1)
    fields=vcat(s.fields[1:3,:],permutedims(theta))
    bforce=BEKTraditionalForcing.lambda_cf(s.Ro)*(s.Tw-1)
    BEKTraditionalForcing.ForcingSolution(s.Ro,s.Co,s.operators,fields,bforce,
                                           s.Hinf,s.residual,s.iterations)
end

function main()
    current=high_tail_seed()
    # On the high-tail sheet Hinf is a poor local coordinate.  A slightly
    # weaker (less negative) forcing gives the preceding state reliably.
    previous=BEKTraditionalForcing.solve_forcing_fixed_b(current.B+0.001,current;
                                                          tolerance=TOL)
    states=BEKTraditionalForcing.ForcingSolution[previous,current]
    step=0.001
    turned=false
    for _ in 1:500
        trial=nothing
        local accepted=false
        for _ in 1:12
            try
                trial=BEKTraditionalForcing.solve_forcing_pseudoarclength(
                    previous,current;step=step,tolerance=TOL)
                accepted=true
                break
            catch
                step*=0.5
            end
        end
        accepted || break
        push!(states,trial)
        if trial.Hinf < current.Hinf
            turned=true
        end
        previous,current=current,trial
        step=min(0.003,1.15step)
        turned && current.Hinf < -0.03 && break
    end
    idx=argmax([s.Hinf for s in states])
    closest=states[idx]
    outdir=joinpath(@__DIR__,"..","results","full_bek_traditional_availability",
                    "tail_turning_validation")
    mkpath(outdir)
    token=@sprintf("Ro_%+.5f_N%d_a%.2f_b%.2f_c%.2f",TARGET_RO,N,A,MAP_B,C)
    token=replace(token,"+"=>"p","-"=>"m","."=>"p")
    path=joinpath(outdir,token*".csv")
    open(path,"w") do io
        println(io,"point,Hinf,B,Tw,residual,ell_T")
        for (i,s) in enumerate(states)
            @printf(io,"%d,%.12e,%.12e,%.12e,%.5e,%.12e\n",i,s.Hinf,s.B,
                    BEKTraditionalForcing.wall_temperature(s.B,s.Ro),s.residual,
                    BEKTraditionalForcing.thermal_tail_length(s))
        end
    end
    summary=joinpath(outdir,"summary.csv")
    newfile=!isfile(summary)
    open(summary,"a") do io
        newfile && println(io,"Ro,N,a,b,c,npoints,turned,Hinf_closest,Tw_closest,B_closest,ell_T,residual")
        @printf(io,"%.8f,%d,%.8f,%.8f,%.8f,%d,%s,%.12e,%.12e,%.12e,%.12e,%.5e\n",
                TARGET_RO,N,A,MAP_B,C,length(states),turned,closest.Hinf,
                BEKTraditionalForcing.wall_temperature(closest.B,closest.Ro),
                closest.B,BEKTraditionalForcing.thermal_tail_length(closest),closest.residual)
    end
    @printf("Ro=%+.5f N=%d map=(%.2f,%.2f,%.2f) closest H=%+.8e Tw=%.10f ellT=%.3e turned=%s points=%d\n",
            TARGET_RO,N,A,MAP_B,C,closest.Hinf,
            BEKTraditionalForcing.wall_temperature(closest.B,closest.Ro),
            BEKTraditionalForcing.thermal_tail_length(closest),turned,length(states))
    println("Output: $path")
end

main()
