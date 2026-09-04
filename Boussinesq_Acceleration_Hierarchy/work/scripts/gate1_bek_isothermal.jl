#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
using .BEKIsothermal
using Printf
include(joinpath(@__DIR__, "..", "..", "baselines", "saddle_node", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode

function main()
    outdir = joinpath(@__DIR__, "..", "results", "gate1_isothermal")
    mkpath(outdir)
    vk_ref = solve_vk_isothermal(; zmax=80.0, degree=160, tolerance=2e-11)
    ros = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
    solutions = BEKSolution[]
    seed = nothing
    # Seed the first mapped solve from the existing finite-domain von Karman
    # solution; the mapped solver then continues this state in Ro.
    seed_op = BEKIsothermal._operators(120, 2.0, 0.6, 0.5)
    seed_targets = min.(seed_op.eta[1:end-1], vk_ref.z[end])
    mapped = barycentric_interpolate(vk_ref.z, vk_ref.weights, vk_ref.fields,
                                     seed_targets)[1:3, :]
    mapped = hcat(mapped, [vk_ref.fields[1,end], 0.0, 1.0])
    seed = vec(mapped')
    for Ro in ros
        sol = solve_isothermal(Ro; degree=120, a=2.0, b=0.6, c=0.5, initial=seed, tolerance=2e-10)
        push!(solutions, sol)
        write_solution(joinpath(outdir, @sprintf("bek_Ro_%+.2f.dat", Ro)), sol)
        @printf("Ro=%+.2f Co=%+.8f residual=%.3e iterations=%d Hinf=%+.8f\n",
                Ro, sol.Co, sol.residual, sol.iterations, sol.fields[1, end])
        seed = vec(sol.fields')
    end

    open(joinpath(outdir, "summary.csv"), "w") do io
        println(io, "Ro,Co,residual,iterations,Hinf")
        for sol in solutions
            @printf(io, "%.8f,%.8f,%.5e,%d,%.12f\n", sol.Ro, sol.Co, sol.residual,
                    sol.iterations, sol.fields[1, end])
        end
    end

    # Regression target: the BEK equations at Ro=-1, Co=2 are exactly the
    # isothermal von Karman equations after U=-F and W=-H.
    bek_vk = solve_isothermal(-1.0; degree=120, a=2.0, b=0.6, c=0.5, tolerance=2e-11)
    targets = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
    max_error = 0.0
    for zq in targets
        k = argmin(abs.(bek_vk.operators.eta .- zq))
        ref = vk_state(vk_ref, zq)
        # BEK fields are [H,F,G]; vk_state returns [H,F',F,G',G,T',T].
        max_error = max(max_error, maximum(abs.([bek_vk.fields[1,k]-ref[1],
                                                   bek_vk.fields[2,k]-ref[3],
                                                   bek_vk.fields[3,k]-ref[5]])))
    end
    println("von Karman regression: residual=$(residual_norm(bek_vk)) max_field_error=$(max_error)")
    println("Outputs: $outdir")
end

main()
