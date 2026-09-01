#!/usr/bin/env julia
"""Compute principal-branch rotor--stator base flows and export Tecplot DAT."""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using Printf

const ROTOR_ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE_NPZ = joinpath(ROTOR_ROOT, "reference", "baseflow_Res1000.npz")
const DEFAULT_OUTPUT = joinpath(ROTOR_ROOT, "tecplot", "data", "baseflow_profiles",
                                "principal_branch_Tw1.00_1.16_step0.04.dat")

function parameter_grid(start, stop, step)
    count = round(Int, (stop - start) / step)
    grid = start .+ step .* collect(0:count)
    isapprox(grid[end], stop; atol=1e-12, rtol=0) || error("Temperature interval is not divisible by step")
    grid[end] = stop
    return grid
end

function write_dat(path, temperatures, profiles, points, re_h, pr)
    mkpath(dirname(path))
    variables = ("z", "H", "F", "F_z", "G", "G_z", "T", "T_z", "pressure_gradient", "Tw")
    open(path, "w") do io
        println(io, "TITLE = \"Rotor-stator principal-branch base flows\"")
        println(io, "VARIABLES = ", join(("\"$name\"" for name in variables), ' '))
        println(io, "DATASETAUXDATA Model = \"traditional_centrifugal\"")
        @printf(io, "DATASETAUXDATA Re_h = \"%.16g\"\n", re_h)
        @printf(io, "DATASETAUXDATA Pr = \"%.16g\"\n", pr)
        for (temperature, profile) in zip(temperatures, profiles)
            values = rotor_state(profile, points)
            @printf(io, "ZONE T=\"Tw=%.2f\", I=%d, DATAPACKING=POINT\n", temperature, length(points))
            @printf(io, "AUXDATA Tw = \"%.16g\"\n", temperature)
            @printf(io, "AUXDATA Pressure_Gradient = \"%.16g\"\n", profile.pressure)
            for index in eachindex(points)
                row = (points[index], values[1,index], values[2,index], values[3,index],
                       values[4,index], values[5,index], values[6,index], values[7,index],
                       profile.pressure, temperature)
                println(io, join((@sprintf("%.16g", value) for value in row), ' '))
            end
        end
    end
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}(
        "tw-start" => 1.0, "tw-stop" => 1.16, "tw-step" => 0.04,
        "internal-step" => 0.005, "profile-points" => 2001,
        "re-h" => 1000.0, "pr" => 0.72, "tol" => 2e-9, "degree" => 100,
        "output" => DEFAULT_OUTPUT)
    args = parse_cli(arguments, defaults)
    isodd(args["profile-points"]) && args["profile-points"] >= 3 ||
        error("--profile-points must be an odd integer >= 3")
    requested = parameter_grid(args["tw-start"], args["tw-stop"], args["tw-step"])
    internal = parameter_grid(args["tw-start"], args["tw-stop"], args["internal-step"])
    temperatures = sort(unique(vcat(requested, internal)))
    config = RotorConfig(re_h=args["re-h"], pr=args["pr"], tolerance=args["tol"],
                         collocation_degree=args["degree"], model="traditional_centrifugal")
    isothermal = solve_rotor_isothermal(config, REFERENCE_NPZ)
    branch = continue_rotor_temperature(temperatures, isothermal;
                                        profile_points=args["profile-points"])
    indices = [argmin(abs.(branch.columns["Tw"] .- temperature)) for temperature in requested]
    selected = branch.profiles[indices]
    points = collect(range(0.0, 1.0; length=args["profile-points"]))
    write_dat(abspath(args["output"]), requested, selected, points, args["re-h"], args["pr"])
    validation = Dict{String,Any}[]
    for (temperature, profile) in zip(requested, selected)
        values = rotor_state(profile, points)
        mass_simpson = simpson_uniform(view(values, 2, :), points)
        mass_end = -(values[1,end] - values[1,1]) / (2sqrt(args["re-h"]))
        boundary = maximum(abs.((values[1,1], values[2,1], values[4,1]-1,
                                values[6,1]-temperature, values[1,end], values[2,end],
                                values[4,end], values[6,end]-1)))
        push!(validation, Dict("Tw" => temperature, "pressure_gradient" => profile.pressure,
              "profile_points" => length(points), "bvp_nodes" => length(profile.z),
              "max_rms_residual" => profile.residual, "boundary_residual_max" => boundary,
              "mass_integral_F_simpson" => mass_simpson,
              "mass_integral_F_from_H_endpoints" => mass_end))
    end
    summary = splitext(abspath(args["output"]))[1] * "_validation.json"
    write_json(summary, validation)
    println("Wrote $(args["output"]): $(length(requested)) zones, $(length(points)) points/zone")
    println("Wrote ", summary)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
