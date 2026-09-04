#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using Printf
using DelimitedFiles
using Plots

const DEGREE = 120
const MAP_A = 2.0
const MAP_B = 0.6
const MAP_C = 0.5
const TOLERANCE = 1e-10
const PROFILE_TEMPERATURES = [1.0, 1.05, 1.1, 1.2, 1.5, 1.8, 1.95, 1.99]

function temperature_targets()
    unique(sort(vcat(
        collect(1.0:0.01:1.2),
        collect(1.22:0.02:1.8),
        collect(1.81:0.01:1.95),
        collect(1.96:0.01:1.99),
        PROFILE_TEMPERATURES,
    )))
end

function profile_metrics(solution)
    eta = solution.operators.eta
    finite = findall(isfinite, eta)
    H, F, G, T = (view(solution.fields, row, :) for row in 1:4)
    D = solution.operators.deta
    Fp, Gp, Tp = D*F, D*G, D*T
    peak_local = argmax(F[finite])
    peak = finite[peak_local]
    condition = fixed_tw_condition(solution)
    (Tw=solution.Tw, Hinf=solution.Hinf,
     Fmax=F[peak], eta_Fmax=eta[peak], Fp0=Fp[1], Gp0=Gp[1],
     minus_Tp0=-Tp[1], chi_wall=2-solution.Tw,
     residual=solution.residual, sigma_min=condition.sigma_min,
     sigma_ratio=condition.ratio)
end

function trace_fixed_tw()
    seed = solve_consistent_isothermal(-1.0; gamma=1.0, Pr=0.72,
        degree=DEGREE, a=MAP_A, b=MAP_B, c=MAP_C,
        tolerance=TOLERANCE)
    solutions = Dict{Float64,ConsistentSolution}(1.0 => seed)
    for tw in temperature_targets()[2:end]
        seed = solve_consistent_fixed_tw(tw, seed; tolerance=TOLERANCE)
        solutions[tw] = seed
        @printf("Tw=%.3f Hinf=%+.10f residual=%.2e\n",
                tw, seed.Hinf, seed.residual)
    end
    solutions
end

function pseudoarclength_crosscheck(solutions)
    targets = temperature_targets()
    previous = solutions[targets[1]]
    current = solutions[targets[2]]
    rows = NamedTuple[]
    push!(rows, (step=0, Tw=previous.Tw, Hinf=previous.Hinf,
                 dTw=NaN, residual=previous.residual))
    push!(rows, (step=1, Tw=current.Tw, Hinf=current.Hinf,
                 dTw=current.Tw-previous.Tw, residual=current.residual))
    index = 1
    while current.Tw < 1.989 && index < 1000
        candidate = solve_consistent_pseudoarclength(previous, current;
            step=0.05, tolerance=5TOLERANCE)
        index += 1
        push!(rows, (step=index, Tw=candidate.Tw, Hinf=candidate.Hinf,
                     dTw=candidate.Tw-current.Tw,
                     residual=candidate.residual))
        previous, current = current, candidate
    end
    rows
end

function write_outputs(solutions, arc_rows)
    output = joinpath(@__DIR__, "..", "results",
                      "vonkarman_consistent_baseflow")
    profiles = joinpath(output, "profiles")
    mkpath(profiles)
    targets = temperature_targets()
    metrics = [profile_metrics(solutions[tw]) for tw in targets]

    open(joinpath(output, "branch.csv"), "w") do stream
        println(stream, "Tw,Hinf,Fmax,eta_Fmax,Fp0,Gp0,minus_Tp0,chi_wall,residual,sigma_min,sigma_ratio")
        for row in metrics
            @printf(stream, "%.8f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.5e,%.5e,%.5e\n",
                    row.Tw,row.Hinf,row.Fmax,row.eta_Fmax,row.Fp0,row.Gp0,
                    row.minus_Tp0,row.chi_wall,row.residual,row.sigma_min,
                    row.sigma_ratio)
        end
    end
    open(joinpath(output, "pseudoarclength.csv"), "w") do stream
        println(stream, "step,Tw,Hinf,dTw,residual")
        for row in arc_rows
            @printf(stream, "%d,%.12f,%.12f,%.12e,%.5e\n",
                    row.step,row.Tw,row.Hinf,row.dTw,row.residual)
        end
    end
    for tw in PROFILE_TEMPERATURES
        tag = replace(@sprintf("%.2f", tw), "." => "p")
        BEKConsistent.write_solution(
            joinpath(profiles, "profile_Tw_$(tag).csv"), solutions[tw])
    end

    default(; linewidth=2, framestyle=:box, gridalpha=0.25,
            legendfontsize=8, guidefontsize=10, tickfontsize=8)
    branch_plot = plot([row.Tw for row in metrics],
                       [row.Hinf for row in metrics],
                       xlabel="T_w", ylabel="H_infinity",
                       label="consistent branch", marker=:circle,
                       markersize=2.2)
    savefig(branch_plot, joinpath(output, "Hinf_Tw.png"))

    condition_plot = plot([row.Tw for row in metrics],
                          [row.sigma_ratio for row in metrics],
                          yscale=:log10, xlabel="T_w",
                          ylabel="scaled sigma_min/sigma_max",
                          label="fixed-T_w Jacobian", marker=:circle,
                          markersize=2.2)
    savefig(condition_plot, joinpath(output, "jacobian_condition.png"))

    selected = [solutions[tw] for tw in PROFILE_TEMPERATURES]
    palette = cgrad(:viridis, length(selected), categorical=true)
    panels = [plot() for _ in 1:4]
    labels = ("H", "F", "G", "T")
    for (index, solution) in enumerate(selected)
        finite = isfinite.(solution.operators.eta)
        eta = solution.operators.eta[finite]
        keep = eta .<= 20
        for row in 1:4
            plot!(panels[row], eta[keep], solution.fields[row,finite][keep],
                  color=palette[index],
                  label=row == 1 ? @sprintf("T_w=%.2f",solution.Tw) : "")
        end
    end
    for row in 1:4
        xlabel!(panels[row], "eta")
        ylabel!(panels[row], labels[row])
    end
    profile_plot = plot(panels...; layout=(2,2), size=(1000,720),
                        legend=:best)
    savefig(profile_plot, joinpath(output, "profiles.png"))
    output, metrics
end

function main()
    solutions = trace_fixed_tw()
    arc_rows = pseudoarclength_crosscheck(solutions)
    output, metrics = write_outputs(solutions, arc_rows)

    ratios = [row.sigma_ratio for row in metrics]
    minimum_index = argmin(ratios)
    d_tw = [row.dTw for row in arc_rows if isfinite(row.dTw)]
    h = [row.Hinf for row in metrics]
    dh = diff(h)
    h_extrema = count(index -> dh[index]*dh[index+1] <= 0,
                      1:length(dh)-1)
    @printf("minimum fixed-Tw sigma ratio %.6e at Tw=%.6f\n",
            ratios[minimum_index], metrics[minimum_index].Tw)
    @printf("pseudoarclength minimum forward dTw %.6e\n", minimum(d_tw))
    @printf("Hinf extrema detected: %d\n", h_extrema)
    println("Outputs: $output")
end

main()
