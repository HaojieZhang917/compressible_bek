#!/usr/bin/env julia
"""Continuation analysis of thermally stratified rotor--stator similarity flow."""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using JSON3
using Plots
using Printf

const ROTOR_ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE_NPZ = joinpath(ROTOR_ROOT, "reference", "baseflow_Res1000.npz")
const DEFAULT_OUTPUT = joinpath(ROTOR_ROOT, "data", "boussinesq_singularity_results")

function branch_rows(branch::RotorBranch)
    names = ["pressure_gradient", "Tw", "thermal_rossby", "rotor_heat_flux_Tz",
             "G_mid", "H_min", "F_min", "F_max", "nodes", "max_rms_residual"]
    return [Dict(name => branch.columns[name][i] for name in names)
            for i in eachindex(branch.columns["Tw"])], names
end

function plot_branch(branch, folds, re_h, path)
    pressure = branch.columns["pressure_gradient"]
    tw = branch.columns["Tw"]
    figure = plot(pressure, tw; linewidth=2, label="branch", xlabel="radial-pressure parameter Pi",
                  ylabel="rotor temperature ratio T_w", title="Traditional Boussinesq two-disk branch (Re_h=$(re_h))",
                  gridalpha=0.25)
    for fold in folds
        scatter!(figure, [fold["pressure_gradient"]], [fold["Tw"]]; label=fold["kind"])
    end
    savefig(figure, path)
end

function plot_three_solutions(solutions, pressures, tw, path)
    z = collect(range(0.0, 1.0; length=1001))
    figure = plot(layout=(2, 2), size=(920, 700), plot_title="Three steady similarity solutions at T_w=$(tw)")
    for (solution, pressure) in zip(solutions, pressures)
        values = rotor_state(solution, z)
        for (panel, row, name) in ((1, 2, "F"), (2, 4, "G"), (3, 1, "H"), (4, 6, "T"))
            plot!(figure[panel], z, values[row, :]; linewidth=1.8,
                  label=@sprintf("Pi=%.6f", pressure), xlabel="z/h", ylabel=name, gridalpha=0.22)
        end
    end
    savefig(figure, path)
end

function run_analysis(args)
    output = abspath(args["output-dir"])
    mkpath(output)
    config = RotorConfig(re_h=args["re-h"], pr=args["pr"],
                         tolerance=args["tol"], collocation_degree=args["degree"],
                         model=args["model"])
    isothermal = solve_rotor_isothermal(config, REFERENCE_NPZ)
    validation = rotor_validation(isothermal, REFERENCE_NPZ)
    if args["model"] == "traditional_centrifugal"
        grid = pressure_grid(isothermal.pressure, args["fine-pressure-step"])
        branch = continue_rotor_pressure(grid, isothermal)
        folds = rotor_turning_points(branch)
    else
        temperatures = collect(1.0:0.01:(args["alternative-max-temperature"] + 0.005))
        branch = continue_rotor_temperature(temperatures, isothermal)
        folds = Dict{String,Any}[]
    end
    rows, names = branch_rows(branch)
    tag = @sprintf("%g", args["re-h"])
    write_csv_rows(joinpath(output, "branch_Re$(tag)_$(args["model"]).csv"), rows; columns=names)
    write_json(joinpath(output, "folds_Re$(tag)_$(args["model"]).json"), folds)
    write_json(joinpath(output, "validation_Re$(tag)_$(args["model"]).json"), validation)

    if args["model"] == "traditional_centrifugal"
        plot_branch(branch, folds, args["re-h"], joinpath(output, "traditional_branch_Re$(tag).png"))
        target = args["profile-temperature"]
        pressures = rotor_roots_at_temperature(branch, target)
        solutions = RotorProfile[]
        for pressure in pressures
            nearest = argmin(abs.(branch.columns["pressure_gradient"] .- pressure))
            seed = branch.profiles[nearest]
            fixed_pressure_seed = RotorProfile(seed.z, seed.D, seed.weights, seed.fields,
                                               pressure, seed.tw, seed.residual, seed.config)
            push!(solutions, solve_rotor_fixed_temperature(target, fixed_pressure_seed))
        end
        if !isempty(solutions)
            plot_three_solutions(solutions, pressures, target,
                                 joinpath(output, @sprintf("three_solutions_Tw%.3f.png", target)))
            z = collect(range(0.0, 1.0; length=1001))
            profile_rows = Dict{String,Any}[]
            for (branch_index, (solution, pressure)) in enumerate(zip(solutions, pressures))
                values = rotor_state(solution, z)
                for j in eachindex(z)
                    push!(profile_rows, Dict("branch" => branch_index,
                         "pressure_gradient" => pressure, "z" => z[j], "H" => values[1,j],
                         "F" => values[2,j], "G" => values[4,j], "T" => values[6,j]))
                end
            end
            write_csv_rows(joinpath(output, @sprintf("three_solutions_Tw%.3f.csv", target)),
                           profile_rows; columns=["branch", "pressure_gradient", "z", "H", "F", "G", "T"])
        end
    else
        figure = plot(branch.columns["Tw"], branch.columns["pressure_gradient"];
                      linewidth=2, xlabel="T_w", ylabel="Pi", label="Pi",
                      title="Soong rotating-force continuation (Re_h=$(args["re-h"]))")
        savefig(figure, joinpath(output, "soong_temperature_continuation_Re$(tag).png"))
    end
    return Dict("config" => args, "validation" => validation, "folds" => folds,
                "branch_points" => length(branch.profiles))
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}(
        "re-h" => 1000.0, "pr" => 0.72, "tol" => 2e-9, "degree" => 100,
        "fine-pressure-step" => 5e-5, "profile-temperature" => 1.16,
        "alternative-max-temperature" => 2.0,
        "model" => "traditional_centrifugal", "output-dir" => DEFAULT_OUTPUT,
    )
    args = parse_cli(arguments, defaults;
                     choices=Dict("model" => ["traditional_centrifugal", "soong_rotating_forces"]))
    summary = run_analysis(args)
    JSON3.pretty(stdout, summary; allow_inf=true)
    println()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
