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

function extrema(H, Tw)
    output = NamedTuple[]
    for i in 2:length(H)-1
        (Tw[i] - Tw[i-1]) * (Tw[i+1] - Tw[i]) <= 0 || continue
        q = hcat(ones(3), H[i-1:i+1], H[i-1:i+1].^2) \ Tw[i-1:i+1]
        hc = -q[2] / (2q[3])
        H[i-1] <= hc <= H[i+1] || continue
        push!(output, (kind=q[3] < 0 ? "maximum" : "minimum",
                       Hinf=hc, Tw=q[1] + q[2] * hc + q[3] * hc^2))
    end
    output
end

function branch_at_ro(Ro)
    solution = solve_thermal_isothermal(Ro; degree=N, a=A, b=B, c=C, tolerance=TOL)
    h = solution.Hinf
    H = Float64[]
    Tw = Float64[]
    step = 0.002
    while h < -0.05 - 1e-12
        hnext = min(h + step, -0.05)
        try
            solution = solve_thermal_fixed_h(hnext, solution; tolerance=TOL)
            h = hnext
            push!(H, h)
            push!(Tw, solution.Tw)
            step = min(0.004, 1.2step)
        catch
            step *= 0.5
            step >= 1e-4 || break
        end
    end
    H, Tw, extrema(H, Tw)
end

function main()
    outdir = joinpath(@__DIR__, "..", "results", "traditional_bek_cusp")
    mkpath(outdir)
    # First bracket the cusp from the von Karman endpoint toward Ro=-0.9.
    ros = collect(-0.9940:0.0002:-0.9930)
    rows = NamedTuple[]
    for Ro in ros
        H, Tw, folds = branch_at_ro(Ro)
        @printf("Ro=%+.3f points=%d folds=%d\n", Ro, length(H), length(folds))
        for fold in folds
            @printf("  %s Hinf=%+.10f Tw=%.10f\n", fold.kind, fold.Hinf, fold.Tw)
            push!(rows, (Ro=Ro, kind=fold.kind, Hinf=fold.Hinf, Tw=fold.Tw))
        end
    end
    open(joinpath(outdir, "traditional_cusp_bracket.csv"), "w") do io
        println(io, "Ro,kind,Hinf,Tw")
        for row in rows
            @printf(io, "%.8f,%s,%.12f,%.12f\n", row.Ro, row.kind, row.Hinf, row.Tw)
        end
    end
end

main()
