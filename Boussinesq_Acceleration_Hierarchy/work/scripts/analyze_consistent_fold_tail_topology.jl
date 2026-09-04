#!/usr/bin/env julia

using Printf
using Plots

function option(prefix,default)
    for argument in ARGS
        startswith(argument,prefix) &&
            return split(argument,"=";limit=2)[2]
    end
    default
end

const DIRECTORY = option("--directory=",joinpath(@__DIR__,"..","results",
    "consistent_fold_tail_continuation","production_N120"))

function numeric_csv(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines),','))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = split(line,',')
        parsed = Tuple(something(tryparse(Float64,value),NaN)
                       for value in values)
        push!(rows,NamedTuple{Tuple(header)}(parsed))
    end
    rows
end

function cusp_brackets(folds)
    output = NamedTuple[]
    for index in 1:length(folds)-1
        left,right = folds[index],folds[index+1]
        y1,y2 = left.cusp_indicator,right.cusp_indicator
        y1*y2 < 0 || continue
        max(abs(y1),abs(y2)) < 5e-4 || continue
        weight = -y1/(y2-y1)
        interpolate(field) = getfield(left,field) + weight *
                             (getfield(right,field)-getfield(left,field))
        push!(output,(point_left=Int(left.point),
            point_right=Int(right.point),Ro=interpolate(:Ro),
            Tw=interpolate(:Tw),Hinf=interpolate(:Hinf),
            bracket_indicator=max(abs(y1),abs(y2))))
    end
    output
end

function interpolate_tail(tail,ro)
    ro < first(tail).Ro && return nothing
    ro > last(tail).Ro && return nothing
    for index in 1:length(tail)-1
        left,right = tail[index],tail[index+1]
        left.Ro <= ro <= right.Ro || continue
        weight = (ro-left.Ro)/(right.Ro-left.Ro)
        return left.Tw_tail + weight*(right.Tw_tail-left.Tw_tail)
    end
    nothing
end

function fold_tail_gaps(folds,tail)
    output = NamedTuple[]
    for fold in folds
        fold.Hinf > -0.015 || continue
        value = interpolate_tail(tail,fold.Ro)
        value === nothing && continue
        push!(output,(Ro=fold.Ro,fold_Tw=fold.Tw,Hinf=fold.Hinf,
                      tail_Tw=value,gap=fold.Tw-value))
    end
    output
end

function main()
    left = numeric_csv(joinpath(DIRECTORY,"fold_left.csv"))
    right = numeric_csv(joinpath(DIRECTORY,"fold_right.csv"))
    tail = numeric_csv(joinpath(DIRECTORY,"tail_boundary.csv"))
    cusps = cusp_brackets(right)
    gaps = fold_tail_gaps(right,tail)

    open(joinpath(DIRECTORY,"cusp_brackets.csv"),"w") do io
        println(io,"point_left,point_right,Ro,Tw,Hinf,max_bracket_indicator")
        for row in cusps
            @printf(io,"%d,%d,%.12f,%.12f,%.12f,%.12e\n",
                    row.point_left,row.point_right,row.Ro,row.Tw,row.Hinf,
                    row.bracket_indicator)
        end
    end
    open(joinpath(DIRECTORY,"fold_tail_gap.csv"),"w") do io
        println(io,"Ro,fold_Tw,fold_Hinf,Tw_tail,fold_minus_tail")
        for row in gaps
            @printf(io,"%.12f,%.12f,%.12f,%.12f,%.12e\n",
                    row.Ro,row.fold_Tw,row.Hinf,row.tail_Tw,row.gap)
        end
    end

    closest = isempty(gaps) ? nothing : gaps[argmin(abs.(getfield.(gaps,:Hinf)))]
    open(joinpath(DIRECTORY,"topology_analysis.txt"),"w") do io
        @printf(io,"fold_left_points=%d\nfold_right_points=%d\n",length(left),length(right))
        @printf(io,"cusp_indicator_zero_brackets=%d\n",length(cusps))
        for (index,row) in enumerate(cusps)
            @printf(io,"cusp_%d_Ro=%.12f\ncusp_%d_Tw=%.12f\ncusp_%d_Hinf=%.12f\n",
                    index,row.Ro,index,row.Tw,index,row.Hinf)
        end
        if closest !== nothing
            @printf(io,"closest_fold_tail_Ro=%.12f\n",closest.Ro)
            @printf(io,"closest_fold_tail_Hinf=%.12f\n",closest.Hinf)
            @printf(io,"closest_fold_tail_Tw_gap=%.12e\n",closest.gap)
        end
    end

    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    p = plot(getfield.(right,:Ro),getfield.(right,:Tw);
             color=:crimson,label="fold locus")
    plot!(p,getfield.(left,:Ro),getfield.(left,:Tw);
          color=:crimson,label="")
    plot!(p,getfield.(tail,:Ro),getfield.(tail,:Tw_tail);
          color=:navy,linestyle=:dash,label="Hinf -> 0 boundary")
    scatter!(p,getfield.(cusps,:Ro),getfield.(cusps,:Tw);
             color=:black,marker=:diamond,markersize=6,label="cusp brackets")
    xlabel!(p,"Ro")
    ylabel!(p,"T_w")
    savefig(p,joinpath(DIRECTORY,"fold_tail_topology.png"))

    @printf("cusp brackets: %d\n",length(cusps))
    for (index,row) in enumerate(cusps)
        @printf("  %d: Ro=%.9f Tw=%.9f Hinf=%.9f\n",
                index,row.Ro,row.Tw,row.Hinf)
    end
    closest === nothing || @printf(
        "closest fold-tail comparison: Hinf=%.3e Tw gap=%.3e\n",
        closest.Hinf,closest.gap)
end

main()
