#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N=120; const A=2.0; const B=0.6; const C=0.5; const TOL=1e-8

function extrema(H,T)
    out=NamedTuple[]
    for i in 2:length(H)-1
        (T[i]-T[i-1])*(T[i+1]-T[i]) <= 0 || continue
        q=hcat(ones(3),H[i-1:i+1],H[i-1:i+1].^2) \ T[i-1:i+1]
        abs(q[3])>1e-14 || continue
        hc=-q[2]/(2q[3]); H[i-1] <= hc <= H[i+1] || continue
        push!(out,(kind=q[3]<0 ? "maximum" : "minimum",Hinf=hc,Tw=q[1]+q[2]*hc+q[3]*hc^2))
    end
    out
end

function trace(Ro)
    s=solve_thermal_isothermal(Ro;degree=N,a=A,b=B,c=C,tolerance=TOL)
    H=Float64[]; T=Float64[]; h=s.Hinf; step=0.002
    while h < -0.05-1e-12
        hn=min(h+step,-0.05)
        try
            s=solve_thermal_fixed_h(hn,s;tolerance=TOL); h=hn
            push!(H,h); push!(T,s.Tw); step=min(0.004,1.2step)
        catch
            step*=0.5; step >= 1e-4 || break
        end
    end
    H,T,extrema(H,T)
end

outdir=joinpath(@__DIR__,"..","results","negative_ro_fold_motion"); mkpath(outdir)
open(joinpath(outdir,"fold_locus.csv"),"w") do io
    println(io,"Ro,kind,Hinf,Tw")
    for Ro in -1.0:0.01:-0.90
        H,T,fs=trace(Ro)
        @printf("Ro=%+.2f points=%d folds=%d\n",Ro,length(H),length(fs))
        for f in fs
            @printf("  %-7s Hinf=%+.10f Tw=%.10f\n",f.kind,f.Hinf,f.Tw)
            @printf(io,"%.4f,%s,%.12f,%.12f\n",Ro,f.kind,f.Hinf,f.Tw)
        end
    end
end
println("Outputs: ",outdir)
