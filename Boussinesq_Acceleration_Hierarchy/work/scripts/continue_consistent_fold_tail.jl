#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using Printf
using Plots

function option(prefix,default)
    for argument in ARGS
        startswith(argument,prefix) &&
            return split(argument,"=";limit=2)[2]
    end
    default
end

const DEGREE = parse(Int,option("--degree=","120"))
const TAG = option("--tag=","production")
const TOLERANCE = parse(Float64,option("--tol=","1e-9"))
const MAX_FOLD_POINTS = parse(Int,option("--fold-points=","180"))
const MAP_A = 2.0
const MAP_B = 0.6
const MAP_C = 0.5
const PRANDTL = 0.72
const GAMMA = 1.0
const TAIL_LEVELS = [-0.02,-0.01,-0.005,-0.0025]

function approach_tw(ro,targets)
    solution = solve_consistent_isothermal(ro;degree=DEGREE,a=MAP_A,
        b=MAP_B,c=MAP_C,Pr=PRANDTL,gamma=GAMMA,
        tolerance=TOLERANCE)
    for tw in targets
        solution = solve_consistent_fixed_tw(tw,solution;
                                               tolerance=TOLERANCE)
    end
    solution
end

function initial_fold()
    targets = [1.10,1.20,1.30,1.40,1.50,1.60,1.65,1.68,1.70,
               1.71,1.715,1.72,1.721,1.7215,1.72165]
    seed = approach_tw(-0.5,targets)
    seed = solve_consistent_fixed_h(-0.297584,seed;
                                     tolerance=TOLERANCE)
    solve_consistent_fold(seed;tolerance=TOLERANCE)
end

function fold_diagnostics(fold)
    solution = fold.solution
    state = BEKConsistent.pack(solution.fields)
    _,jacobian = BEKConsistent._residual_jacobian(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators)
    rownorm = max.([norm(view(jacobian,row,:))
                    for row in axes(jacobian,1)],1e-14)
    decomposition = svd(jacobian ./ rownorm)
    left = (1 ./ rownorm).*decomposition.U[:,end]
    left ./= norm(left)
    vector = fold.nullvector
    K = BEKConsistent._jacobian_directional_derivative(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators,vector)
    quadratic = dot(left,K*vector)
    transversality = dot(left,BEKConsistent._tw_derivative(
        state,solution.operators))
    (scaled_sigma_min=decomposition.S[end],
     scaled_sigma_ratio=decomposition.S[end]/decomposition.S[1],
     quadratic=quadratic,transversality=transversality,
     cusp_indicator=quadratic/transversality)
end

function neighbouring_fold(center,target_ro)
    seed = solve_consistent_fixed_h_at_ro(center.solution.Hinf,target_ro,
        center.solution;tolerance=TOLERANCE)
    solve_consistent_fold(seed;null_seed=center.nullvector,
                          tolerance=TOLERANCE)
end

function trace_fold_direction(center,direction;max_points=MAX_FOLD_POINTS)
    neighbour = neighbouring_fold(center,center.solution.Ro+direction*0.005)
    folds = ConsistentFold[center,neighbour]
    previous,current = center,neighbour
    step = 0.012
    failures = 0
    for _ in 3:max_points
        accepted = false
        candidate = current
        local_step = step
        while local_step >= 2e-5
            try
                candidate = continue_consistent_fold(previous,current;
                    step=local_step,tolerance=TOLERANCE)
                accepted = true
                break
            catch
                local_step *= 0.5
                failures += 1
            end
        end
        accepted || break
        if dot(candidate.nullvector,current.nullvector) < 0
            candidate = ConsistentFold(candidate.solution,
                -candidate.nullvector,candidate.residual,candidate.iterations)
        end
        push!(folds,candidate)
        previous,current = current,candidate
        step = min(0.022,1.15local_step)
        ro,h,tw = current.solution.Ro,current.solution.Hinf,current.solution.Tw
        (ro < -0.95 || ro > -0.015 || h > -0.0015 ||
         tw < 0.98 || tw > 2.02) && break
    end
    folds,failures
end

function seed_tail_levels()
    solution = approach_tw(-0.25,[1.10,1.20,1.30,1.35,1.38,1.40])
    h = solution.Hinf
    while h < TAIL_LEVELS[1]-1e-12
        target = min(TAIL_LEVELS[1],h+0.025)
        solution = solve_consistent_fixed_h(target,solution;
                                             tolerance=TOLERANCE)
        h = target
    end
    seeds = Dict{Float64,ConsistentSolution}(TAIL_LEVELS[1]=>solution)
    for target in TAIL_LEVELS[2:end]
        solution = solve_consistent_fixed_h(target,solution;
                                             tolerance=TOLERANCE)
        seeds[target] = solution
    end
    seeds
end

function solve_ro_target(target,current;depth=0)
    try
        return solve_consistent_fixed_h_at_ro(current.Hinf,target,current;
                                               tolerance=TOLERANCE)
    catch
        depth >= 7 && rethrow()
        midpoint = (current.Ro+target)/2
        intermediate = solve_ro_target(midpoint,current;depth=depth+1)
        return solve_ro_target(target,intermediate;depth=depth+1)
    end
end

function trace_tail_level(seed)
    # Follow the regular low-temperature arm toward the Ekman limit.
    right = ConsistentSolution[seed]
    current = seed
    for target in (-0.24:0.01:-0.02)
        try
            current = solve_ro_target(target,current)
            push!(right,current)
            (current.Tw > 2.05 || current.Tw < 0.98) && break
        catch
            break
        end
    end

    # The opposite direction develops a turning point in the Ro projection.
    # Seed the full-state arclength tangent from the regular right neighbour,
    # then continue through that projection turning point onto the high-tail
    # sheet.
    neighbour = solve_ro_target(-0.24,seed)
    backward = ConsistentSolution[seed]
    previous,current = neighbour,seed
    step = 0.008
    for _ in 1:180
        accepted = false
        candidate = current
        local_step = step
        while local_step >= 2e-5
            try
                candidate = solve_consistent_pseudoarclength_ro_full(
                    previous,current;step=local_step,tolerance=TOLERANCE)
                accepted = true
                break
            catch
                local_step *= 0.5
            end
        end
        accepted || break
        push!(backward,candidate)
        previous,current = current,candidate
        step = min(0.015,1.15local_step)
        (current.Ro < -0.9 || current.Ro > -0.015 ||
         current.Tw > 2.08 || current.Tw < 0.98) && break
    end
    vcat(reverse(backward[2:end]),right)
end

function tail_from_fold(fold)
    solution = fold.solution
    hs = solution.Hinf .* collect(1.0:-0.1:0.1)
    states = ConsistentSolution[solution]
    current = solution
    for h in hs[2:end]
        try
            current = solve_consistent_fixed_h(h,current;
                                                tolerance=TOLERANCE)
            push!(states,current)
        catch
            return nothing
        end
    end
    selected = states[max(1,end-5):end]
    hh = getfield.(selected,:Hinf)
    tw = getfield.(selected,:Tw)
    coefficients = hcat(ones(length(hh)),hh,hh.^2) \ tw
    linear = hcat(ones(2),hh[end-1:end]) \ tw[end-1:end]
    (Ro=solution.Ro,fold_Tw=solution.Tw,fold_H=solution.Hinf,
     tail_Tw=coefficients[1],linear_Tw=linear[1],
     estimate_difference=abs(coefficients[1]-linear[1]),
     gap=coefficients[1]-solution.Tw)
end

function tail_extrapolation(curves)
    dictionaries = Dict{Float64,Dict{Int,ConsistentSolution}}()
    for h in TAIL_LEVELS
        regular = [s for s in curves[h]
                   if abs(100s.Ro-round(100s.Ro)) < 1e-8 &&
                      -0.2500001 <= s.Ro <= -0.0199999]
        dictionaries[h] = Dict(round(Int,100s.Ro)=>s for s in regular)
    end
    common = reduce(intersect,[Set(keys(dictionaries[h]))
                               for h in TAIL_LEVELS])
    rows = NamedTuple[]
    for key in sort(collect(common))
        hs = TAIL_LEVELS
        tws = [dictionaries[h][key].Tw for h in hs]
        coefficients = hcat(ones(length(hs)),hs,hs.^2) \ tws
        linear = hcat(ones(2),hs[end-1:end]) \ tws[end-1:end]
        push!(rows,(Ro=key/100,Tw=coefficients[1],
                    linear_Tw=linear[1],
                    estimate_difference=abs(coefficients[1]-linear[1])))
    end
    rows
end

function write_folds(path,label,folds)
    open(path,"w") do io
        println(io,"branch,point,Ro,Tw,Hinf,residual,scaled_sigma_min,scaled_sigma_ratio,fold_quadratic,tw_transversality,cusp_indicator")
        for (index,fold) in enumerate(folds)
            d = fold_diagnostics(fold)
            s = fold.solution
            @printf(io,"%s,%d,%.12f,%.12f,%.12f,%.5e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    label,index,s.Ro,s.Tw,s.Hinf,fold.residual,
                    d.scaled_sigma_min,d.scaled_sigma_ratio,
                    d.quadratic,d.transversality,d.cusp_indicator)
        end
    end
end

function write_tail_curve(path,curve)
    open(path,"w") do io
        println(io,"Ro,Tw,Hinf,residual,thermal_length")
        for solution in curve
            ell = -1/(solution.Pr*solution.Hinf)
            @printf(io,"%.12f,%.12f,%.12f,%.5e,%.12e\n",
                    solution.Ro,solution.Tw,solution.Hinf,
                    solution.residual,ell)
        end
    end
end

function make_plot(output,left,right,tail)
    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8)
    p = plot(;xlabel="Ro",ylabel="T_w",legend=:best)
    plot!(p,getfield.(getfield.(left,:solution),:Ro),
          getfield.(getfield.(left,:solution),:Tw);
          color=:crimson,label="fold locus: left")
    plot!(p,getfield.(getfield.(right,:solution),:Ro),
          getfield.(getfield.(right,:solution),:Tw);
          color=:darkorange,label="fold locus: right")
    plot!(p,getfield.(tail,:Ro),getfield.(tail,:Tw);
          color=:navy,linestyle=:dash,label="Hinf -> 0 tail boundary")
    savefig(p,joinpath(output,"fold_tail_parameter_plane.png"))
end

function main()
    output = joinpath(@__DIR__,"..","results",
        "consistent_fold_tail_continuation","$(TAG)_N$(DEGREE)")
    mkpath(output)
    println("Correcting Ro=-0.5 seed fold ...")
    center = initial_fold()
    diagnostics = fold_diagnostics(center)
    @printf("seed: Ro=%.12f Tw=%.12f Hinf=%.12f sigma=%.3e a=%.3e\n",
            center.solution.Ro,center.solution.Tw,center.solution.Hinf,
            diagnostics.scaled_sigma_min,diagnostics.quadratic)

    println("Continuing fold locus toward decreasing Ro ...")
    left,left_failures = trace_fold_direction(center,-1.0)
    println("Continuing fold locus toward increasing Ro ...")
    right,right_failures = trace_fold_direction(center,+1.0)
    write_folds(joinpath(output,"fold_left.csv"),"left",left)
    write_folds(joinpath(output,"fold_right.csv"),"right",right)

    println("Continuing finite-Hinf tail proxies ...")
    seeds = seed_tail_levels()
    curves = Dict{Float64,Vector{ConsistentSolution}}()
    for h in TAIL_LEVELS
        @printf("  Hinf=%+.5f ...\n",h)
        curve = trace_tail_level(seeds[h])
        curves[h] = curve
        token = replace(@sprintf("%.4f",abs(h)),"."=>"p")
        write_tail_curve(joinpath(output,"tail_H_m$(token).csv"),curve)
    end
    tail = tail_extrapolation(curves)
    open(joinpath(output,"tail_boundary.csv"),"w") do io
        println(io,"Ro,Tw_tail,linear_Tw_tail,estimate_difference")
        for row in tail
            @printf(io,"%.12f,%.12f,%.12f,%.12e\n",row.Ro,row.Tw,
                    row.linear_Tw,row.estimate_difference)
        end
    end
    approach = NamedTuple[]
    for (index,fold) in enumerate(right)
        fold.solution.Hinf > -0.025 || continue
        index % 3 == 0 || continue
        estimate = tail_from_fold(fold)
        estimate === nothing || push!(approach,estimate)
    end
    open(joinpath(output,"fold_tail_approach.csv"),"w") do io
        println(io,"Ro,fold_Tw,fold_Hinf,Tw_tail,linear_Tw_tail,estimate_difference,Tw_gap")
        for row in approach
            @printf(io,"%.12f,%.12f,%.12f,%.12f,%.12f,%.12e,%.12e\n",
                    row.Ro,row.fold_Tw,row.fold_H,row.tail_Tw,
                    row.linear_Tw,row.estimate_difference,row.gap)
        end
    end
    open(joinpath(output,"run_summary.txt"),"w") do io
        @printf(io,"degree=%d\na=%.6f\nb=%.6f\nc=%.6f\nPr=%.6f\ngamma=%.6f\n",
                DEGREE,MAP_A,MAP_B,MAP_C,PRANDTL,GAMMA)
        @printf(io,"seed_Ro=%.12f\nseed_Tw=%.12f\nseed_Hinf=%.12f\n",
                center.solution.Ro,center.solution.Tw,center.solution.Hinf)
        @printf(io,"left_points=%d\nright_points=%d\nleft_failures=%d\nright_failures=%d\n",
                length(left),length(right),left_failures,right_failures)
        @printf(io,"tail_boundary_points=%d\n",length(tail))
        @printf(io,"fold_tail_approach_points=%d\n",length(approach))
    end
    make_plot(output,left,right,tail)
    println("Outputs: $output")
end

main()
