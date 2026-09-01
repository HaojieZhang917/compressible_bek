#!/usr/bin/env julia
"""Compute finite-domain von Karman branches for several outer boundaries."""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using Plots
using Printf

const VK_ROOT = normpath(joinpath(@__DIR__, ".."))

function trace_branch(zmax, hmin, hmax, dh; degree=100)
    branch = trace_vk_h_branch(; zmax=Float64(zmax), degree=degree,
                               h_start=Float64(hmin), h_stop=Float64(hmax),
                               dh=Float64(dh), tolerance=2e-8)
    rows = Dict{String,Any}[]
    for profile in branch.profiles
        wall = vk_state(profile, 0.0)
        push!(rows, Dict(
            "zmax" => zmax, "Hinf" => profile.hinf, "Tw" => profile.tw,
            "Fp0" => wall[2], "Gp0" => wall[4], "Tp0" => wall[6],
            "nodes" => length(profile.z),
        ))
    end
    return branch, rows
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}(
        "zmax" => [15.0, 20.0, 25.0, 30.0, 40.0],
        "hmin" => -0.75, "hmax" => -0.05, "dh" => 0.0025,
        "degree" => 100,
        "out" => joinpath(VK_ROOT, "data", "finite_domain_branches"),
    )
    # Multiple zmax values are accepted as a comma-separated list.
    args = parse_cli(arguments, defaults)
    zmax_values = args["zmax"] isa Vector ? args["zmax"] : parse.(Float64, split(string(args["zmax"]), ','))
    output = abspath(args["out"])
    mkpath(output)
    turns = Dict{String,Any}[]
    all_rows = Dict{String,Any}[]
    branches = Tuple{Float64,VKBranch}[]
    for zmax in zmax_values
        branch, rows = trace_branch(zmax, args["hmin"], args["hmax"], args["dh"];
                                    degree=args["degree"])
        append!(all_rows, rows)
        push!(branches, (zmax, branch))
        tag = replace(@sprintf("%g", zmax), "." => "p")
        write_csv_rows(joinpath(output, "Hinf_Tw_zmax_$(tag).csv"), rows;
                       columns=["zmax", "Hinf", "Tw", "Fp0", "Gp0", "Tp0", "nodes"])
        for (index, point) in enumerate(vk_turning_points(branch)[1:min(2, end)])
            push!(turns, Dict("zmax" => zmax, "turning_point" => index,
                              "Hinf" => point.x, "Tw" => point.y))
        end
    end
    columns = ["zmax", "Hinf", "Tw", "Fp0", "Gp0", "Tp0", "nodes"]
    write_csv_rows(joinpath(output, "Hinf_Tw_all_zmax.csv"), all_rows; columns)
    write_csv_rows(joinpath(output, "turning_points_by_zmax.csv"), turns;
                   columns=["zmax", "turning_point", "Hinf", "Tw"])
    main_plot = plot(xlabel="Hinf", ylabel="Tw", gridalpha=0.25,
                     title="Finite-domain Boussinesq branches: Tw(Hinf)")
    zoom_plot = plot(xlabel="Hinf", ylabel="Tw", gridalpha=0.25,
                     title="Turning-point region", xlim=(-0.58, -0.10),
                     ylim=(1.025, 1.057))
    for (zmax, branch) in branches
        plot!(main_plot, branch.hinf, branch.tw; label="zmax=$zmax", linewidth=2)
        plot!(zoom_plot, branch.hinf, branch.tw; label="zmax=$zmax", linewidth=2)
        selected = filter(row -> row["zmax"] == zmax, turns)
        for panel in (main_plot, zoom_plot)
            scatter!(panel, [row["Hinf"] for row in selected],
                     [row["Tw"] for row in selected]; color=:black,
                     markersize=4, label=false)
        end
    end
    savefig(main_plot, joinpath(output, "Hinf_Tw_curves.png"))
    savefig(main_plot, joinpath(output, "Hinf_Tw_curves.pdf"))
    savefig(zoom_plot, joinpath(output, "Hinf_Tw_turning_region.png"))
    savefig(zoom_plot, joinpath(output, "Hinf_Tw_turning_region.pdf"))
    open(joinpath(output, "summary.md"), "w") do io
        println(io, "# Domain-dependent Boussinesq turning points\n")
        println(io, "Branches were parameterized by Hinf on [$(args["hmin"]), $(args["hmax"])] with dHinf=$(args["dh"]).\n")
        println(io, "| zmax | point | Hinf | Tw |\n|---:|---:|---:|---:|")
        for row in turns
            @printf(io, "| %g | %d | %.9f | %.9f |\n", row["zmax"],
                    row["turning_point"], row["Hinf"], row["Tw"])
        end
    end
    println("Wrote finite-domain branch data to ", output)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
