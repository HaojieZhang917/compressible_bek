#!/usr/bin/env julia

using Printf
using Plots

const ROOT = joinpath(@__DIR__,"..","results",
    "consistent_fold_tail_continuation")
const RUNS = [(80,"convergence_N80"),(120,"production_N120"),
              (160,"convergence_N160")]
const LEVELS = [-0.10,-0.05,-0.02,-0.01,-0.0075,-0.005]

function numeric_csv(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines),','))
    rows = NamedTuple[]
    for line in lines[2:end]
        values = split(line,',')
        parsed = Tuple(something(tryparse(Float64,value),NaN)
                       for value in values)
        push!(rows,NamedTuple{Tuple(header)}(parsed))
    end
    rows
end

function primary_arm(rows)
    # Starting at Ro=-0.5, the converged arm initially moves monotonically
    # toward larger Ro and Hinf.  The first reversal is the resolution-limited
    # near-tail structure whose migration is being tested here.
    endpoint = length(rows)
    for index in 3:length(rows)-1
        if rows[index+1].Ro < rows[index].Ro &&
           rows[index].Ro > rows[index-1].Ro
            endpoint = index
            break
        end
    end
    rows[1:endpoint]
end

function at_h_level(rows,target)
    for index in 1:length(rows)-1
        left,right = rows[index],rows[index+1]
        left.Hinf <= target <= right.Hinf || continue
        weight = (target-left.Hinf)/(right.Hinf-left.Hinf)
        return (Ro=left.Ro+weight*(right.Ro-left.Ro),
                Tw=left.Tw+weight*(right.Tw-left.Tw))
    end
    nothing
end

function main()
    intersections = NamedTuple[]
    migrations = NamedTuple[]
    arms = Dict{Int,Vector{NamedTuple}}()
    for (degree,directory) in RUNS
        all = numeric_csv(joinpath(ROOT,directory,"fold_right.csv"))
        arm = primary_arm(all)
        arms[degree] = arm
        endpoint = arm[end]
        push!(migrations,(degree=degree,Ro=endpoint.Ro,Tw=endpoint.Tw,
                          Hinf=endpoint.Hinf,point=Int(endpoint.point)))
        for level in LEVELS
            value = at_h_level(arm,level)
            value === nothing && continue
            push!(intersections,(degree=degree,Hinf=level,
                Ro=value.Ro,Tw=value.Tw))
        end
    end

    open(joinpath(ROOT,"primary_fold_level_convergence.csv"),"w") do io
        println(io,"N,Hinf,Ro,Tw")
        for row in intersections
            @printf(io,"%d,%.12f,%.12f,%.12f\n",
                    row.degree,row.Hinf,row.Ro,row.Tw)
        end
    end
    open(joinpath(ROOT,"apparent_cusp_migration.csv"),"w") do io
        println(io,"N,last_primary_point,Ro,Tw,Hinf")
        for row in migrations
            @printf(io,"%d,%d,%.12f,%.12f,%.12f\n",
                    row.degree,row.point,row.Ro,row.Tw,row.Hinf)
        end
    end

    # Extrapolate only the highest-resolution primary arm, using its last
    # points that remain before the projection reversal.  This is a numerical
    # estimate of a fold--tail limiting coordinate, not an exact Hinf=0 solve.
    arm = arms[160]
    selected = arm[max(1,end-7):end]
    h = getfield.(selected,:Hinf)
    ro = getfield.(selected,:Ro)
    tw = getfield.(selected,:Tw)
    design = hcat(ones(length(h)),h,h.^2)
    rofit = design\ro
    twfit = design\tw
    open(joinpath(ROOT,"fold_tail_limit_estimate.txt"),"w") do io
        @printf(io,"source_N=160\npoints=%d\n",length(selected))
        @printf(io,"Hinf_range=[%.12e,%.12e]\n",minimum(h),maximum(h))
        @printf(io,"Ro_at_Hinf_0_quadratic=%.12f\n",rofit[1])
        @printf(io,"Tw_at_Hinf_0_quadratic=%.12f\n",twfit[1])
        @printf(io,"interpretation=resolution_limited_extrapolation_not_exact_tail_solution\n")
    end


    tail = numeric_csv(joinpath(ROOT,"production_N120","tail_boundary.csv"))
    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    p = plot(;xlabel="Ro",ylabel="T_w",legend=:best)
    colors = Dict(80=>:gray55,120=>:darkorange,160=>:crimson)
    for (degree,_) in RUNS
        branch = arms[degree]
        plot!(p,getfield.(branch,:Ro),getfield.(branch,:Tw);
              color=colors[degree],label="primary fold, N=$degree")
        scatter!(p,[branch[end].Ro],[branch[end].Tw];
                 color=colors[degree],marker=:circle,markersize=4,label="")
    end
    plot!(p,getfield.(tail,:Ro),getfield.(tail,:Tw_tail);
          color=:navy,linestyle=:dash,label="Hinf -> 0 boundary, N=120")
    savefig(p,joinpath(ROOT,"fold_tail_resolution_convergence.png"))

    println("Apparent primary-arm reversal migration:")
    for row in migrations
        @printf("  N=%3d Ro=%+.8f Tw=%.8f Hinf=%+.4e\n",
                row.degree,row.Ro,row.Tw,row.Hinf)
    end
    @printf("N=160 quadratic Hinf->0 estimate: Ro=%+.8f Tw=%.8f\n",
            rofit[1],twfit[1])
end

main()
