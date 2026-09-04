#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N=60; const A=2.0; const B=0.6; const C=0.5; const TOL=1e-8

function extrema(H,T)
    out=NamedTuple[]
    for i in 2:length(H)-1
        (T[i]-T[i-1])*(T[i+1]-T[i]) <= 0 || continue
        q=hcat(ones(3),H[i-1:i+1],H[i-1:i+1].^2) \ T[i-1:i+1]
        abs(q[3]) > 1e-13 || continue
        hc=-q[2]/(2q[3]); H[i-1] <= hc <= H[i+1] || continue
        push!(out,(kind=q[3]<0 ? "max" : "min",Hinf=hc,Tw=q[1]+q[2]*hc+q[3]*hc^2))
    end
    out
end

function trace(Ro)
    s=solve_thermal_isothermal(Ro;degree=N,a=A,b=B,c=C,tolerance=TOL)
    h=s.Hinf; h0=h; H=Float64[]; T=Float64[]; step=0.01; status="reached_H0"
    while h < -0.001
        hn=min(h+step,-0.001)
        try
            s=solve_thermal_fixed_h(hn,s;tolerance=TOL); h=hn
            push!(H,h); push!(T,s.Tw); step=min(0.02,1.2step)
        catch
            if step > 2e-4
                step*=0.5
            else
                status="newton_failure"; break
            end
        end
    end
    fs=extrema(H,T)
    return h0,H,T,fs,status
end

out=joinpath(@__DIR__,"..","results","systematic_thermal_availability"); mkpath(out)
ros=[-1.0,-0.999,-0.998,-0.995,-0.994,-0.99,-0.95,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1]
open(joinpath(out,"summary.csv"),"w") do io
    println(io,"Ro,Co,Lambda_cf,Hiso,status,Hlast,Twlast,nfolds,fold1_H,fold1_Tw,fold2_H,fold2_Tw")
    for Ro in ros
        h0,H,T,fs,status=trace(Ro); Co=2-Ro-Ro^2; lam=Co^2/(4Ro)
        hl=isempty(H) ? h0 : H[end]; tl=isempty(T) ? 1.0 : T[end]
        f1=length(fs)>=1 ? fs[1] : (Hinf=NaN,Tw=NaN); f2=length(fs)>=2 ? fs[2] : (Hinf=NaN,Tw=NaN)
        @printf("Ro=%+.3f lambda=%+.6f status=%s Hlast=%+.5f Twlast=%.6f folds=%d\n",Ro,lam,status,hl,tl,length(fs))
        @printf(io,"%.4f,%.8f,%.8f,%.12f,%s,%.12f,%.12f,%d,%.12f,%.12f,%.12f,%.12f\n",Ro,Co,lam,h0,status,hl,tl,length(fs),f1.Hinf,f1.Tw,f2.Hinf,f2.Tw)
    end
end
println("Outputs: ",out)
