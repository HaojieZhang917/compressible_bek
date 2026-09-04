#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKThermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKThermal, .BEKConsistent
using LinearAlgebra
using Printf
using Plots

const ROSSBY_NUMBERS = [-1.0, -0.75, -0.5, -0.25, -0.1, 0.25, 0.5, 1.0]
const DEGREE = 120
const MAP_A = 2.0
const MAP_B = 0.6
const MAP_C = 0.5
const PRANDTL = 0.72
const GAMMA = 1.0
const TOLERANCE = 1e-9
const MAX_TW = 1.99
const MAX_STEP = 0.01
const MIN_STEP = 1e-5

pack_solution(solution::ThermalSolution) =
    vcat(BEKThermal.pack(solution.fields), solution.Tw)
pack_solution(solution::ConsistentSolution) =
    vcat(BEKConsistent.pack(solution.fields), solution.Tw)

function isothermal(model::Symbol, ro)
    model == :traditional && return solve_thermal_isothermal(ro;
        degree=DEGREE, a=MAP_A, b=MAP_B, c=MAP_C,
        tolerance=TOLERANCE)
    model == :consistent && return solve_consistent_isothermal(ro;
        gamma=GAMMA, Pr=PRANDTL, degree=DEGREE,
        a=MAP_A, b=MAP_B, c=MAP_C, tolerance=TOLERANCE)
    error("unknown model $model")
end

fixed_tw(model::Symbol, tw, seed) = model == :traditional ?
    solve_thermal_fixed_tw(tw, seed; tolerance=TOLERANCE) :
    solve_consistent_fixed_tw(tw, seed; tolerance=TOLERANCE)

fixed_h(model::Symbol, h, seed) = model == :traditional ?
    solve_thermal_fixed_h(h, seed; tolerance=TOLERANCE) :
    solve_consistent_fixed_h(h, seed; tolerance=TOLERANCE)

condition(model::Symbol, solution) = model == :traditional ?
    thermal_fixed_tw_condition(solution) : fixed_tw_condition(solution)

pseudoarc(model::Symbol, previous, current, step) = model == :traditional ?
    solve_thermal_pseudoarclength(previous, current; step=step,
                                  tolerance=5TOLERANCE) :
    solve_consistent_pseudoarclength(previous, current; step=step,
                                     tolerance=5TOLERANCE)

function trace_fixed_tw(model::Symbol, ro)
    seed = isothermal(model, ro)
    solutions = [seed]
    tw = 1.0
    step = MAX_STEP
    failures = 0
    while tw < MAX_TW - 1e-12 && step >= MIN_STEP
        target = min(MAX_TW, tw + step)
        try
            candidate = fixed_tw(model, target, seed)
            candidate.residual <= 5TOLERANCE ||
                error("residual $(candidate.residual) exceeds tolerance")
            push!(solutions, candidate)
            seed = candidate
            tw = target
            step = min(MAX_STEP, 1.25step)
        catch
            failures += 1
            step *= 0.5
        end
    end
    (solutions=solutions, reached=tw >= MAX_TW-1e-10,
     failures=failures, final_step=step)
end

function trace_pseudoarc(model::Symbol, solutions; max_steps=240)
    length(solutions) >= 2 || return typeof(first(solutions))[]
    previous, current = solutions[end-1], solutions[end]
    states = [previous, current]
    step = model == :traditional ? 0.002 : 0.02
    minimum_step = model == :traditional ? 2e-6 : 2e-4
    negative = 0
    for _ in 1:max_steps
        accepted = false
        candidate = current
        while step >= minimum_step
            try
                candidate = pseudoarc(model, previous, current, step)
                accepted = true
                break
            catch
                step *= 0.5
            end
        end
        accepted || break
        push!(states, candidate)
        negative = candidate.Tw < current.Tw ? negative+1 : 0
        previous, current = current, candidate
        step = min(model == :traditional ? 0.01 : 0.03, 1.15step)
        negative >= 10 && break
    end
    states
end

function trace_through_fold_in_h(model::Symbol, solutions; max_steps=80)
    current = solutions[end]
    states = [current]
    step = model == :traditional ? 5e-4 : 1e-3
    negative = 0
    for _ in 1:max_steps
        accepted = false
        candidate = current
        local_step = step
        while local_step >= 2e-7
            try
                candidate = fixed_h(model, current.Hinf+local_step, current)
                accepted = true
                break
            catch
                local_step *= 0.5
            end
        end
        accepted || break
        push!(states,candidate)
        negative = candidate.Tw < current.Tw ? negative+1 : 0
        current = candidate
        step = min(model == :traditional ? 0.003 : 0.006,
                   1.2local_step)
        negative >= 10 && break
    end
    states
end

function quadratic_vertex(x, y)
    coefficients = hcat(ones(length(x)), x, x.^2) \ y
    abs(coefficients[3]) > 1e-14 || error("degenerate quadratic")
    xv = -coefficients[2]/(2coefficients[3])
    yv = coefficients[1]+coefficients[2]*xv+coefficients[3]*xv^2
    xv, yv, coefficients[3]
end

function refine_fold(model::Symbol, arc)
    length(arc) >= 3 || error("too few pseudo-arclength states")
    index = argmax(getfield.(arc, :Tw))
    1 < index < length(arc) || error("pseudo-arclength did not bracket fold")
    triplet = arc[index-1:index+1]
    h = getfield.(triplet, :Hinf)
    tw = getfield.(triplet, :Tw)
    order = sortperm(h)
    h, tw = h[order], tw[order]
    hstar, _, curvature = quadratic_vertex(h, tw)
    curvature < 0 || error("bracketed extremum is not a maximum")
    minimum(h) <= hstar <= maximum(h) || error("fold vertex outside bracket")
    halfwidth = max((maximum(h)-minimum(h))/4, 1e-5)
    seed = arc[index]
    center = fixed_h(model, hstar, seed)
    for _ in 1:6
        left = fixed_h(model, hstar-halfwidth, center)
        right = fixed_h(model, hstar+halfwidth, center)
        hs = [left.Hinf, center.Hinf, right.Hinf]
        tws = [left.Tw, center.Tw, right.Tw]
        trial, _, q2 = quadratic_vertex(hs, tws)
        # At widths comparable to the nonlinear tolerance, the fitted
        # curvature becomes roundoff dominated.  Retain the previous
        # converged vertex rather than corrupting it with that final fit.
        q2 < 0 || break
        minimum(hs) <= trial <= maximum(hs) || break
        hstar = trial
        center = fixed_h(model, hstar, center)
        halfwidth /= 3
    end
    delta = max(halfwidth, 2e-7)
    left = fixed_h(model, hstar-delta, center)
    right = fixed_h(model, hstar+delta, center)
    ds = norm(pack_solution(right)-pack_solution(left))
    dTw_ds = (right.Tw-left.Tw)/ds
    cond = condition(model, center)
    (solution=center, condition=cond, dTw_ds=dTw_ds,
     bracket_positive=(center.Tw-left.Tw)/
        norm(pack_solution(center)-pack_solution(left)),
     bracket_negative=(right.Tw-center.Tw)/
        norm(pack_solution(right)-pack_solution(center)))
end

function tail_limit(solutions)
    count = min(8, length(solutions))
    selected = solutions[end-count+1:end]
    h = getfield.(selected, :Hinf)
    tw = getfield.(selected, :Tw)
    coefficients = hcat(ones(count), h, h.^2) \ tw
    coefficients[1]
end

function classify(model::Symbol, ro, trace)
    last = trace.solutions[end]
    if trace.reached
        topology = model == :traditional && ro == 1.0 ? "passive" : "smooth"
        cond = condition(model, last)
        return (topology=topology, Tw_c=last.Tw, H_c=last.Hinf,
                sigma_min=cond.sigma_min, sigma_ratio=cond.ratio,
                dTw_ds=NaN, sigma_location="Tw=1.99",
                last=last, arc=typeof(last)[], fold=nothing)
    end
    if -0.05 < last.Hinf <= 0
        cond = condition(model, last)
        return (topology="tail", Tw_c=tail_limit(trace.solutions), H_c=0.0,
                sigma_min=cond.sigma_min, sigma_ratio=cond.ratio,
                dTw_ds=NaN, sigma_location="last_resolved",
                last=last, arc=typeof(last)[], fold=nothing)
    end
    arc = trace_pseudoarc(model, trace.solutions)
    crossed = any(i -> arc[i].Tw < arc[i-1].Tw, 2:length(arc))
    if !crossed
        # Fixed Hinf is regular at an ordinary Tw fold and supplies separated
        # full-state seeds when the final fixed-Tw points are nearly identical.
        # The resulting curve is still parameterised and differentiated by
        # full-state arclength below.
        arc = trace_through_fold_in_h(model, trace.solutions)
        crossed = any(i -> arc[i].Tw < arc[i-1].Tw, 2:length(arc))
    end
    crossed || error("$model Ro=$ro stopped away from tail without a fold crossing")
    fold = refine_fold(model, arc)
    return (topology="fold", Tw_c=fold.solution.Tw,
            H_c=fold.solution.Hinf,
            sigma_min=fold.condition.sigma_min,
            sigma_ratio=fold.condition.ratio,
            dTw_ds=fold.dTw_ds, sigma_location="refined_fold",
            last=last, arc=arc, fold=fold)
end

function token(ro)
    replace(@sprintf("%+.2f",ro), "+"=>"p", "-"=>"m", "."=>"p")
end

function write_branch(path, solutions)
    open(path, "w") do io
        println(io, "point,Ro,Tw,Hinf,residual,thermal_length")
        for (i, solution) in enumerate(solutions)
            ell = solution.Hinf < 0 ? -1/(PRANDTL*solution.Hinf) : Inf
            @printf(io, "%d,%.8f,%.12f,%.12f,%.5e,%.12e\n",
                    i,solution.Ro,solution.Tw,solution.Hinf,
                    solution.residual,ell)
        end
    end
end

function write_arc(path, arc)
    open(path, "w") do io
        println(io, "point,Ro,Tw,Hinf,residual")
        for (i, solution) in enumerate(arc)
            @printf(io, "%d,%.8f,%.12f,%.12f,%.5e\n",
                    i,solution.Ro,solution.Tw,solution.Hinf,
                    solution.residual)
        end
    end
end

function make_plot(output, rows, branches)
    default(; linewidth=2, framestyle=:box, gridalpha=0.25,
            legendfontsize=7, guidefontsize=10, tickfontsize=8)
    p = plot(; xlabel="T_w", ylabel="H_infinity")
    colors = Dict(:traditional=>:darkorange, :consistent=>:navy)
    for model in (:traditional, :consistent)
        first_curve = true
        for ro in ROSSBY_NUMBERS
            states = branches[(model,ro)]
            label = first_curve ? String(model) : ""
            plot!(p, getfield.(states,:Tw), getfield.(states,:Hinf);
                  color=colors[model], alpha=0.72, label=label)
            first_curve = false
        end
    end
    for row in rows
        row.topology == "fold" || continue
        scatter!(p,[row.Tw_c],[row.H_c];marker=:diamond,markersize=5,
                 color=colors[row.model],label="")
    end
    savefig(p,joinpath(output,"unified_Hinf_Tw.png"))
end

function main()
    output = joinpath(@__DIR__,"..","results",
                      "two_model_topology_comparison")
    mkpath(output)
    rows = NamedTuple[]
    branches = Dict{Tuple{Symbol,Float64},Any}()
    for model in (:traditional, :consistent)
        for ro in ROSSBY_NUMBERS
            @printf("%s Ro=%+.2f: tracing ...\n",model,ro)
            trace = trace_fixed_tw(model,ro)
            result = classify(model,ro,trace)
            branches[(model,ro)] = trace.solutions
            write_branch(joinpath(output,"branch_$(model)_Ro_$(token(ro)).csv"),
                         trace.solutions)
            !isempty(result.arc) && write_arc(
                joinpath(output,"arc_$(model)_Ro_$(token(ro)).csv"),
                result.arc)
            push!(rows,(model=model,Ro=ro,topology=result.topology,
                Tw_c=result.Tw_c,H_c=result.H_c,
                sigma_min=result.sigma_min,sigma_ratio=result.sigma_ratio,
                dTw_ds=result.dTw_ds,
                sigma_location=result.sigma_location,
                Tw_last=result.last.Tw,H_last=result.last.Hinf,
                failures=trace.failures,residual=result.last.residual))
            @printf("  %-7s Tw_c=%.10f H_c=%+.9f sigma=%.3e dTw/ds=%+.3e\n",
                    result.topology,result.Tw_c,result.H_c,
                    result.sigma_min,result.dTw_ds)
            if result.fold !== nothing
                @printf("  fold crossing slopes: %+.3e, %+.3e\n",
                        result.fold.bracket_positive,
                        result.fold.bracket_negative)
            end
        end
    end
    open(joinpath(output,"topology_table.csv"),"w") do io
        println(io,"model,Ro,topology,Tw_c,Hinf_c,scaled_sigma_min,scaled_sigma_ratio,dTw_ds,sigma_location,Tw_last,Hinf_last,failures,last_residual")
        for row in rows
            @printf(io,"%s,%.8f,%s,%.12f,%.12f,%.12e,%.12e,%.12e,%s,%.12f,%.12f,%d,%.5e\n",
                    row.model,row.Ro,row.topology,row.Tw_c,row.H_c,
                    row.sigma_min,row.sigma_ratio,row.dTw_ds,
                    row.sigma_location,row.Tw_last,row.H_last,
                    row.failures,row.residual)
        end
    end
    make_plot(output,rows,branches)
    println("Outputs: $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
