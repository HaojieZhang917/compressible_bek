#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
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

solution_vector(solution) = vcat(BEKConsistent.pack(solution.fields),
                                  solution.Tw)

function trace_fixed_tw(ro; degree=DEGREE)
    seed = solve_consistent_isothermal(ro; gamma=GAMMA, Pr=PRANDTL,
        degree=degree, a=MAP_A, b=MAP_B, c=MAP_C,
        tolerance=TOLERANCE)
    solutions = ConsistentSolution[seed]
    tw = 1.0
    step = MAX_STEP
    failures = 0
    while tw < MAX_TW - 1e-12 && step >= MIN_STEP
        target = min(MAX_TW, tw + step)
        try
            candidate = solve_consistent_fixed_tw(target, seed;
                                                   tolerance=TOLERANCE)
            candidate.residual <= TOLERANCE || error(
                "stored residual $(candidate.residual) exceeds tolerance")
            push!(solutions, candidate)
            seed = candidate
            tw = target
            step = min(MAX_STEP, 1.25step)
        catch
            failures += 1
            step *= 0.5
        end
    end
    (solutions=solutions, final_step=step, failures=failures,
     reached=tw >= MAX_TW - 1e-10)
end

function branch_metrics(solutions)
    initial_ratio = fixed_tw_condition(first(solutions)).ratio
    rows = NamedTuple[]
    for solution in solutions
        condition = fixed_tw_condition(solution)
        length_scale = solution.Hinf < 0 ?
            -1 / (solution.Pr * solution.Hinf) : Inf
        push!(rows, (Ro=solution.Ro, Tw=solution.Tw,
                     Hinf=solution.Hinf, thermal_length=length_scale,
                     residual=solution.residual,
                     sigma_ratio=condition.ratio,
                     sigma_relative=condition.ratio/initial_ratio))
    end
    rows
end

function trace_pseudoarclength(solutions; max_steps=180)
    length(solutions) >= 2 || return NamedTuple[]
    previous, current = solutions[end-1], solutions[end]
    rows = NamedTuple[]
    cumulative = 0.0
    push!(rows, (step=0, arclength=cumulative, Tw=previous.Tw,
                 Hinf=previous.Hinf, dTw=NaN,
                 residual=previous.residual))
    increment = norm(solution_vector(current)-solution_vector(previous))
    cumulative += increment
    push!(rows, (step=1, arclength=cumulative, Tw=current.Tw,
                 Hinf=current.Hinf, dTw=current.Tw-previous.Tw,
                 residual=current.residual))
    arc_step = 0.02
    negative_steps = 0
    for index in 2:max_steps
        candidate = nothing
        accepted = false
        while arc_step >= 2e-4
            try
                candidate = solve_consistent_pseudoarclength(
                    previous, current; step=arc_step,
                    tolerance=5TOLERANCE)
                accepted = true
                break
            catch
                arc_step *= 0.5
            end
        end
        accepted || break
        increment = norm(solution_vector(candidate)-solution_vector(current))
        cumulative += increment
        delta_tw = candidate.Tw-current.Tw
        push!(rows, (step=index, arclength=cumulative, Tw=candidate.Tw,
                     Hinf=candidate.Hinf, dTw=delta_tw,
                     residual=candidate.residual))
        negative_steps = delta_tw < 0 ? negative_steps+1 : 0
        previous, current = current, candidate
        arc_step = min(0.03, 1.1arc_step)
        negative_steps >= 8 && break
    end
    rows
end

function quadratic_extremum(x, y, index; maximum=false)
    1 < index < length(x) || return (x=x[index], y=y[index], resolved=false)
    xx = x[index-1:index+1]
    yy = y[index-1:index+1]
    coefficients = hcat(ones(3), xx, xx.^2) \ yy
    abs(coefficients[3]) > 1e-14 ||
        return (x=x[index], y=y[index], resolved=false)
    xe = -coefficients[2]/(2coefficients[3])
    ye = coefficients[1]+coefficients[2]*xe+coefficients[3]*xe^2
    correct_curvature = maximum ? coefficients[3] < 0 : coefficients[3] > 0
    (x=xe, y=ye, resolved=correct_curvature && xx[1] <= xe <= xx[end])
end

function fold_from_arc(rows)
    isempty(rows) && return nothing
    tw = [row.Tw for row in rows]
    s = [row.arclength for row in rows]
    h = [row.Hinf for row in rows]
    index = argmax(tw)
    extremum = quadratic_extremum(s, tw, index; maximum=true)
    extremum.resolved || return nothing
    hcoeff = hcat(ones(3), s[index-1:index+1],
                  s[index-1:index+1].^2) \ h[index-1:index+1]
    hfold = hcoeff[1]+hcoeff[2]*extremum.x+hcoeff[3]*extremum.x^2
    (Tw=extremum.y, Hinf=hfold, index=index)
end

function tail_limit(rows)
    usable = [row for row in rows if isfinite(row.Hinf)]
    length(usable) >= 5 || return NaN
    selected = usable[max(1,end-7):end]
    h = [row.Hinf for row in selected]
    tw = [row.Tw for row in selected]
    coefficients = hcat(ones(length(h)), h, h.^2) \ tw
    coefficients[1]
end

function classify(trace, metrics)
    last = metrics[end]
    minimum_relative = minimum(row.sigma_relative for row in metrics)
    if trace.reached
        if abs(last.Hinf) > 0.05 && minimum_relative > 1e-3
            return (mechanism="smooth_to_density_boundary",
                    endpoint_tw=last.Tw, endpoint_h=last.Hinf,
                    arc=NamedTuple[], fold=nothing)
        end
        return (mechanism="unresolved", endpoint_tw=last.Tw,
                endpoint_h=last.Hinf, arc=NamedTuple[], fold=nothing)
    end
    if -0.05 < last.Hinf <= 0
        return (mechanism="thermal_tail", endpoint_tw=tail_limit(metrics),
                endpoint_h=0.0, arc=NamedTuple[], fold=nothing)
    end
    arc = trace_pseudoarclength(trace.solutions)
    fold = fold_from_arc(arc)
    if fold !== nothing && any(row -> isfinite(row.dTw) && row.dTw < 0, arc)
        return (mechanism="fold", endpoint_tw=fold.Tw,
                endpoint_h=fold.Hinf, arc=arc, fold=fold)
    end
    (mechanism="unresolved", endpoint_tw=last.Tw,
     endpoint_h=last.Hinf, arc=arc, fold=nothing)
end

function ro_tag(ro)
    replace(@sprintf("%+.2f", ro), "+"=>"p", "-"=>"m", "."=>"p")
end

function write_branch(path, rows)
    open(path, "w") do stream
        println(stream,
            "Ro,Tw,Hinf,thermal_length,residual,sigma_ratio,sigma_relative")
        for row in rows
            @printf(stream,
                "%.8f,%.12f,%.12f,%.12e,%.5e,%.5e,%.5e\n",
                row.Ro,row.Tw,row.Hinf,row.thermal_length,row.residual,
                row.sigma_ratio,row.sigma_relative)
        end
    end
end

function write_arc(path, rows)
    open(path, "w") do stream
        println(stream,"step,arclength,Tw,Hinf,dTw,residual")
        for row in rows
            @printf(stream,"%d,%.12e,%.12f,%.12f,%.12e,%.5e\n",
                    row.step,row.arclength,row.Tw,row.Hinf,row.dTw,
                    row.residual)
        end
    end
end

function make_plots(output, branches, summaries)
    default(; linewidth=2, framestyle=:box, gridalpha=0.25,
            legendfontsize=8, guidefontsize=10, tickfontsize=8)
    p1 = plot(; xlabel="T_w", ylabel="H_infinity")
    p2 = plot(; xlabel="T_w",
              ylabel="sigma ratio / isothermal ratio", yscale=:log10)
    p3 = plot(; xlabel="T_w", ylabel="thermal decay length",
              yscale=:log10)
    colors = cgrad(:viridis, length(ROSSBY_NUMBERS), categorical=true)
    for (index, ro) in enumerate(ROSSBY_NUMBERS)
        rows = branches[ro]
        label = @sprintf("Ro=%+.2f",ro)
        plot!(p1,[row.Tw for row in rows],[row.Hinf for row in rows];
              label=label,color=colors[index])
        plot!(p2,[row.Tw for row in rows],
              [row.sigma_relative for row in rows];
              label=label,color=colors[index])
        plot!(p3,[row.Tw for row in rows],
              [row.thermal_length for row in rows];
              label=label,color=colors[index])
        summary = summaries[index]
        summary.mechanism == "fold" && scatter!(p1,[summary.Tw_endpoint],
            [summary.H_endpoint];marker=:diamond,color=colors[index],label="")
    end
    savefig(p1,joinpath(output,"Hinf_Tw.png"))
    savefig(p2,joinpath(output,"jacobian_relative.png"))
    savefig(p3,joinpath(output,"thermal_length.png"))
end

function main()
    output = joinpath(@__DIR__,"..","results",
                      "consistent_representative_ro_sweep")
    mkpath(output)
    branches = Dict{Float64,Vector{NamedTuple}}()
    summaries = NamedTuple[]
    for ro in ROSSBY_NUMBERS
        @printf("Tracing Ro=%+.2f ...\n",ro)
        trace = trace_fixed_tw(ro)
        metrics = branch_metrics(trace.solutions)
        branches[ro] = metrics
        result = classify(trace,metrics)
        if !isempty(result.arc)
            write_arc(joinpath(output,"arc_Ro_$(ro_tag(ro)).csv"),result.arc)
        end
        write_branch(joinpath(output,"branch_Ro_$(ro_tag(ro)).csv"),metrics)
        minimum_relative = minimum(row.sigma_relative for row in metrics)
        maximum_length = maximum(row.thermal_length for row in metrics)
        push!(summaries,(Ro=ro,mechanism=result.mechanism,
            Tw_endpoint=result.endpoint_tw,H_endpoint=result.endpoint_h,
            Tw_last=metrics[end].Tw,H_last=metrics[end].Hinf,
            minimum_sigma_relative=minimum_relative,
            maximum_thermal_length=maximum_length,
            failures=trace.failures,residual=metrics[end].residual))
        @printf("  %-28s endpoint Tw=%.9f H=%+.8f\n",
                result.mechanism,result.endpoint_tw,result.endpoint_h)
    end
    open(joinpath(output,"summary.csv"),"w") do stream
        println(stream,"Ro,mechanism,Tw_endpoint,H_endpoint,Tw_last,H_last,minimum_sigma_relative,maximum_thermal_length,failures,last_residual")
        for row in summaries
            @printf(stream,
                "%.8f,%s,%.12f,%.12f,%.12f,%.12f,%.5e,%.12e,%d,%.5e\n",
                row.Ro,row.mechanism,row.Tw_endpoint,row.H_endpoint,
                row.Tw_last,row.H_last,row.minimum_sigma_relative,
                row.maximum_thermal_length,row.failures,row.residual)
        end
    end
    make_plots(output,branches,summaries)
    println("Outputs: $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
