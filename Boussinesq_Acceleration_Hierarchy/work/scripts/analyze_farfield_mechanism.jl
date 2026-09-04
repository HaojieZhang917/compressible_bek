#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

states = [
    (label="vk_fold_1", Ro=-1.0, Hinf=-0.53276166, Tw=1.048021731),
    (label="vk_fold_2", Ro=-1.0, Hinf=-0.11330573, Tw=1.030144664),
    (label="near_cusp_upper", Ro=-0.9938, Hinf=-0.44873720, Tw=1.055343990),
    (label="near_cusp_lower", Ro=-0.9938, Hinf=-0.44194023, Tw=1.055343766),
    (label="Ro_m095_tail", Ro=-0.95, Hinf=-0.02, Tw=1.869486417),
    (label="Ro_m090_tail", Ro=-0.90, Hinf=-0.02, Tw=1.899375237),
    (label="Ro_m050_tail", Ro=-0.50, Hinf=-0.02, Tw=1.674120716),
]

outdir=joinpath(@__DIR__,"..","results","traditional_bek_solvability")
mkpath(outdir)
open(joinpath(outdir,"farfield_rates.csv"),"w") do io
    println(io,"label,Ro,Hinf,Tw,K,alpha_T,alpha_V_min,ell_T,ell_V_max,alpha_T_over_alpha_V")
    for s in states
        rates=farfield_decay_rates(s.Ro,s.Hinf)
        av=isempty(rates.velocity) ? NaN : minimum(rates.velocity)
        @printf(io,"%s,%.8f,%.12f,%.12f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                s.label,s.Ro,s.Hinf,s.Tw,rates.coupling,rates.thermal,av,
                1/rates.thermal,1/av,rates.thermal/av)
        @printf("%-18s K=%+.5e alpha_T=%.5e alpha_V=%.5e ell_T=%.2f ell_V=%.2f\n",
                s.label,rates.coupling,rates.thermal,av,1/rates.thermal,1/av)
    end
end
