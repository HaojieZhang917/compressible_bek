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

function quadratic_extrema(x, y)
    out = NamedTuple[]
    for i in 2:length(x)-1
        d1, d2 = y[i] - y[i-1], y[i+1] - y[i]
        d1 * d2 <= 0 || continue
        M = hcat(ones(3), x[i-1:i+1], x[i-1:i+1].^2)
        q = M \ y[i-1:i+1]
        abs(q[3]) > 1e-14 || continue
        xe = -q[2] / (2q[3])
        x[i-1] <= xe <= x[i+1] || continue
        ye = q[1] + q[2] * xe + q[3] * xe^2
        push!(out, (kind=q[3] < 0 ? "maximum" : "minimum", Hinf=xe, Tw=ye))
    end
    out
end

function trace_ro(Ro)
    seed = solve_thermal_isothermal(Ro; degree=N, a=A, b=B, c=C, tolerance=TOL)
    # Start exactly from the converged isothermal state; avoid an artificial
    # H-infinity jump that can exceed the Newton basin for weak rotation.
    h0 = seed.Hinf
    H = Float64[]; Tw = Float64[]; residual = Float64[]
    h = h0
    step = 0.01
    # The first fold is reached well before the long-tail region. Keep this
    # production scan focused on the fold locus; the second fold is retained
    # from the dedicated Ro=-1 continuation.
    while h < -0.25 - 1e-12
        hnext = min(h + step, -0.25)
        try
            candidate = solve_thermal_fixed_h(hnext, seed; tolerance=TOL)
            seed = candidate
            h = hnext
            push!(H, h); push!(Tw, seed.Tw); push!(residual, seed.residual)
            step = min(0.01, 1.25step)
        catch err
            step *= 0.5
            if step < 2e-4
                @printf("Ro=%+.3f stopped at Hinf=%+.5f: %s\n", Ro, h, sprint(showerror, err))
                break
            end
        end
    end
    return H, Tw, residual, quadratic_extrema(H, Tw)
end

function main()
    outdir = joinpath(@__DIR__, "..", "results", "traditional_bek_folds")
    mkpath(outdir)
    ros = collect(-1.0:0.1:1.0)
    fold_rows = NamedTuple[]
    branch_io = open(joinpath(outdir, "traditional_branches.csv"), "w")
    println(branch_io, "Ro,Hinf,Tw,residual")
    for Ro in ros
        H, Tw, residual, folds = trace_ro(Ro)
        for i in eachindex(H)
            @printf(branch_io, "%.8f,%.12f,%.12f,%.5e\n", Ro, H[i], Tw[i], residual[i])
        end
        for f in folds
            push!(fold_rows, (Ro=Ro, kind=f.kind, Hinf=f.Hinf, Tw=f.Tw))
            @printf("Ro=%+.3f %-7s Hinf=%+.8f Tw=%.10f\n", Ro, f.kind, f.Hinf, f.Tw)
        end
    end
    close(branch_io)
    open(joinpath(outdir, "traditional_fold_locus.csv"), "w") do io
        println(io, "Ro,kind,Hinf,Tw")
        for r in fold_rows
            @printf(io, "%.8f,%s,%.12f,%.12f\n", r.Ro, r.kind, r.Hinf, r.Tw)
        end
    end
    println("Outputs: $outdir")
end

main()
