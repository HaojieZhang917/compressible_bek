#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N=120; const A=2.0; const B=0.6; const C=0.5; const TOL=1e-9

function extrema(H,T)
    out=NamedTuple[]
    for i in 2:length(H)-1
        (T[i]-T[i-1])*(T[i+1]-T[i]) <= 0 || continue
        q=hcat(ones(3),H[i-1:i+1],H[i-1:i+1].^2) \ T[i-1:i+1]
        abs(q[3])>1e-13 || continue
        hc=-q[2]/(2q[3]); H[i-1] <= hc <= H[i+1] || continue
        push!(out,(kind=q[3]<0 ? "max" : "min",Hinf=hc,Tw=q[1]+q[2]*hc+q[3]*hc^2))
    end
    out
end

function trace(Ro)
    s=solve_thermal_isothermal(Ro;degree=N,a=A,b=B,c=C,tolerance=TOL)
    H=Float64[]; T=Float64[]; h=s.Hinf; step=0.0005
    while h < -0.05-1e-12
        hn=min(h+step,-0.05)
        try
            s=solve_thermal_fixed_h(hn,s;tolerance=TOL); h=hn
            push!(H,h); push!(T,s.Tw); step=min(0.001,1.15step)
        catch
            step*=0.5; step >= 2e-5 || break
        end
    end
    H,T,extrema(H,T)
end

function scan(path, ros)
    open(path,"w") do io
        println(io,"Ro,kind,Hinf,Tw,deltaH")
        for Ro in ros
            H,T,fs=trace(Ro)
            @printf("Ro=%+.5f points=%d folds=%d\n",Ro,length(H),length(fs))
            if length(fs)==2
                sep=abs(fs[1].Hinf-fs[2].Hinf)
                @printf("  max H=%+.10f Tw=%.10f; min H=%+.10f Tw=%.10f; dH=%.6e\n",fs[1].Hinf,fs[1].Tw,fs[2].Hinf,fs[2].Tw,sep)
            end
            for f in fs
                @printf(io,"%.6f,%s,%.12f,%.12f,%.12e\n",Ro,f.kind,f.Hinf,f.Tw, length(fs)==2 ? abs(fs[1].Hinf-fs[2].Hinf) : NaN)
            end
        end
    end
end

out=joinpath(@__DIR__,"..","results","refined_fold_merger"); mkpath(out)
scan(joinpath(out,"local.csv"), [-0.999,-0.998,-0.997,-0.996,-0.995,-0.994])
println("Outputs: ",out)
