#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
using .BEKIsothermal, .BEKThermal
using Printf

const N = 120
const HFIX = -0.02
const TOL = 1e-8

function seed_at_h(ro)
    s = solve_thermal_isothermal(ro; degree=N, a=2.0, b=0.6, c=0.5,
                                 tolerance=TOL)
    h = s.Hinf
    step = 0.005
    while h < HFIX - 1e-12
        hn = min(h + step, HFIX)
        try
            s = solve_thermal_fixed_h(hn, s; tolerance=TOL)
            h = hn
            step = min(0.01, 1.2 * step)
        catch
            step *= 0.5
            step >= 1e-5 || error("cannot construct long-tail seed")
        end
    end
    s
end

function high_tail_seed()
    s = seed_at_h(-0.5)
    for ro in collect(-0.51:-0.01:-0.99)
        s = solve_thermal_fixed_h_at_ro(HFIX, ro, s; tolerance=TOL)
    end
    for ro in collect(-0.991:-0.001:-0.995)
        s = solve_thermal_fixed_h_at_ro(HFIX, ro, s; tolerance=TOL)
    end
    s
end

function main()
    outdir = joinpath(@__DIR__, "..", "results", "traditional_bek_solvability")
    mkpath(outdir)
    seed = high_tail_seed()
    previous = solve_thermal_fixed_h_at_ro(HFIX, -0.9948, seed; tolerance=TOL)
    current = solve_thermal_fixed_h_at_ro(HFIX, -0.9949, previous; tolerance=TOL)
    path = joinpath(outdir, "long_tail_ro_fullstate_toward_vk.csv")
    open(path, "w") do io
        println(io, "k,Ro,Hinf,Tw,K,alpha_T,alpha_V_min,ell_T,residual")
        arc_step = 0.001
        for k in 1:707
            rates = farfield_decay_rates(current)
            av = minimum(rates.velocity)
            @printf(io, "%d,%.10f,%.12f,%.12f,%.12e,%.12e,%.12e,%.12e,%.5e\n",
                    k, current.Ro, current.Hinf, current.Tw, rates.coupling,
                    rates.thermal, av, 1 / rates.thermal, current.residual)
            @printf("k=%d Ro=%+.10f Tw=%.10f H=%+.5f residual=%.2e\n",
                    k, current.Ro, current.Tw, current.Hinf, current.residual)
            k == 707 && break
            try
                next = solve_thermal_pseudoarclength_ro_full(previous, current;
                                                               step=arc_step,
                                                               tolerance=TOL)
                previous, current = current, next
                arc_step = min(0.001, 1.15 * arc_step)
            catch err
                arc_step *= 0.5
                if arc_step < 2.5e-5
                    println("stop: ", sprint(showerror, err))
                    break
                end
                @printf("retry at step=%.3e after: %s\n", arc_step,
                        sprint(showerror, err))
            end
        end
    end
    println("Outputs: ", path)
end

main()
