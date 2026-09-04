#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N = 120
const A = 2.0
const B = 0.6
const C = 0.5
const TOL = 1e-8
const H_TARGET = -0.02
const ROS = [-1.0, -0.999, -0.998, -0.994, -0.990, -0.950, -0.900, -0.750, -0.500]

function quadratic_extrema(H, T)
    extrema = NamedTuple[]
    for i in 2:length(H)-1
        (T[i]-T[i-1]) * (T[i+1]-T[i]) <= 0 || continue
        q = hcat(ones(3), H[i-1:i+1], H[i-1:i+1].^2) \ T[i-1:i+1]
        abs(q[3]) > 1e-13 || continue
        h = -q[2]/(2q[3])
        H[i-1] <= h <= H[i+1] || continue
        push!(extrema, (kind=q[3] < 0 ? "maximum" : "minimum",
                        Hinf=h, Tw=q[1]+q[2]*h+q[3]*h^2))
    end
    extrema
end

function trace_branch(Ro)
    state = solve_thermal_isothermal(Ro; degree=N, a=A, b=B, c=C, tolerance=TOL)
    H = [state.Hinf]
    T = [state.Tw]
    residual = [state.residual]
    condition = [thermal_fixed_tw_condition(state).ratio]
    h = state.Hinf
    step = 0.004
    status = "reached_H_target"
    while h < H_TARGET - 1e-12
        hnext = min(h + step, H_TARGET)
        try
            state = solve_thermal_fixed_h(hnext, state; tolerance=TOL)
            h = hnext
            push!(H, state.Hinf)
            push!(T, state.Tw)
            push!(residual, state.residual)
            push!(condition, thermal_fixed_tw_condition(state).ratio)
            step = min(0.010, 1.20step)
        catch err
            step *= 0.5
            if step < 1e-5
                status = "fixed_H_continuation_limit"
                @printf("Ro=%+.3f stopped at Hinf=%+.8f, Tw=%.9f: %s\n",
                        Ro, h, state.Tw, sprint(showerror, err))
                break
            end
        end
    end
    return H, T, residual, condition, quadratic_extrema(H, T), status
end

function main()
    outdir = joinpath(@__DIR__, "..", "results", "traditional_bek_solvability")
    mkpath(outdir)
    open(joinpath(outdir, "branch_points.csv"), "w") do points
        open(joinpath(outdir, "summary.csv"), "w") do summary
            println(points, "Ro,Hinf,Tw,tail_length,residual,fixed_Tw_sigma_ratio")
            println(summary, "Ro,status,points,H_last,Tw_last,tail_last,nfolds,fold1_kind,fold1_Hinf,fold1_Tw,fold2_kind,fold2_Hinf,fold2_Tw")
            for Ro in ROS
                H,T,residual,condition,folds,status = trace_branch(Ro)
                for i in eachindex(H)
                    @printf(points, "%.6f,%.12f,%.12f,%.12f,%.5e,%.5e\n", Ro, H[i], T[i], -1/(0.72H[i]), residual[i], condition[i])
                end
                f1 = length(folds) >= 1 ? folds[1] : (kind="", Hinf=NaN, Tw=NaN)
                f2 = length(folds) >= 2 ? folds[2] : (kind="", Hinf=NaN, Tw=NaN)
                @printf(summary, "%.6f,%s,%d,%.12f,%.12f,%.12f,%d,%s,%.12f,%.12f,%s,%.12f,%.12f\n",
                        Ro,status,length(H),H[end],T[end],-1/(0.72H[end]),length(folds),
                        f1.kind,f1.Hinf,f1.Tw,f2.kind,f2.Hinf,f2.Tw)
                @printf("Ro=%+.3f status=%s, Hlast=%+.6f, Twlast=%.7f, tail=%.1f, folds=%d\n",
                        Ro,status,H[end],T[end],-1/(0.72H[end]),length(folds))
            end
        end
    end
    println("Outputs: ", outdir)
end

main()
