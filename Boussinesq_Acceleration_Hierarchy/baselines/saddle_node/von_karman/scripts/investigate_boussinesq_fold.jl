#!/usr/bin/env julia

"""Cross-validate the Boussinesq rotating-disk similarity folds near Tw=1.05."""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using LinearAlgebra
using Plots
using Printf

const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "data", "finite_fold_analysis")
const PR = 0.72

function nearest_profile(branch::VKBranch, hinf::Real)
    return branch.profiles[argmin(abs.(branch.hinf .- hinf))]
end

function exact_profile(branch::VKBranch, hinf::Real; degree::Int, zmax::Float64,
                       tolerance::Float64=2e-10)
    return solve_vk_fixed_h(hinf, nearest_profile(branch, hinf);
                            degree, zmax, tolerance)
end

function domain_fold_estimate(zmax::Real; degree::Int=90)
    branch = trace_vk_h_branch(; zmax=Float64(zmax), degree,
        h_start=-0.70, h_stop=-0.04, dh=0.01, tolerance=5e-8)
    return vk_turning_points(branch)
end

function upper_branch_domain_stability(; base_degree::Int=70)
    targets = (-0.18, -0.14, -0.10, -0.08)
    rows = Dict{String,Any}[]
    for zmax in (20.0, 40.0, 60.0, 80.0)
        degree = max(base_degree, round(Int, base_degree + zmax / 2))
        seed = solve_vk_isothermal(; zmax, degree, tolerance=2e-8)
        for hinf in collect(-0.70:0.02:-0.08)
            seed = solve_vk_fixed_h(hinf, seed; zmax, degree, tolerance=8e-8)
            if any(abs(hinf - target) < 1e-10 for target in targets)
                values = vk_similarity_eigenvalues(seed; degree, zmax)
                leading = first(values)
                push!(rows, Dict("zmax" => zmax, "Hinf" => hinf,
                    "Tw" => seed.tw, "leading_real" => real(leading),
                    "leading_imag" => imag(leading),
                    "unstable_eigenvalues" => count(>(1e-8), real.(values)),
                    "thermal_decay_length" => 1 / (PR * abs(hinf))))
            end
        end
    end
    return rows
end

function save_branch_plot(branch, turns, target_tw, output)
    physical = branch.hinf .< 0
    p1 = plot(branch.tw[physical], branch.hinf[physical]; color=:black,
              linewidth=2, label=false, xlabel="Tw", ylabel="Hinf",
              title="Negative-H branch", gridalpha=0.25)
    p2 = plot(branch.tw[physical], branch.hinf[physical]; color=:black,
              linewidth=2, label=false, xlabel="Tw", ylabel="Hinf",
              title="Fold region", xlim=(1.0375, 1.0510),
              ylim=(-0.72, -0.14), gridalpha=0.25)
    for panel in (p1, p2)
        scatter!(panel, [item.y for item in turns], [item.x for item in turns];
                 color=:crimson, marker=:star5, markersize=7,
                 label="turning points")
        vline!(panel, [target_tw]; color=:grey, linestyle=:dash,
               label="Tw=$target_tw")
    end
    figure = plot(p1, p2; layout=(1, 2), size=(1080, 440))
    savefig(figure, joinpath(output, "branch_diagram.png"))
    savefig(figure, joinpath(output, "branch_diagram.pdf"))
end

function save_profile_plot(profiles, output)
    panels = [plot(; xlabel=label, ylabel="eta", gridalpha=0.25)
              for label in ("F", "G", "H")]
    z = collect(range(0.0, 8.0; length=600))
    for (branch_id, profile) in profiles
        state = vk_state(profile, z)
        for (panel, row) in zip(panels, (3, 5, 1))
            plot!(panel, state[row, :], z; linewidth=2,
                  label=panel === panels[1] ? "branch $branch_id" : false)
        end
    end
    figure = plot(panels...; layout=(1, 3), size=(1200, 380))
    savefig(figure, joinpath(output, "profiles_Tw1p045.png"))
    savefig(figure, joinpath(output, "profiles_Tw1p045.pdf"))
end

function write_report(path, turning_rows, domain_rows, validation,
                      fixed_stability, upper_rows)
    open(path, "w") do io
        println(io, "# Boussinesq similarity-fold investigation\n")
        println(io, "## Turning points\n")
        println(io, "| point | Hinf | Tw | scaled Jacobian sigma_min/sigma_max | thermal decay length |")
        println(io, "|---:|---:|---:|---:|---:|")
        for row in turning_rows
            @printf(io, "| %d | %.9f | %.9f | %.3e | %.4f |\n",
                    row["point"], row["Hinf"], row["Tw"],
                    row["newton_sigma_min_over_max"], row["thermal_decay_length"])
        end
        println(io, "\nThe benchmark-connected branch reaches a saddle-node; continuation in Tw is singular there, whereas Hinf continuation remains regular and exposes the S-shaped branch.\n")
        println(io, "## Domain sensitivity\n")
        println(io, "| zmax | point | Hinf | Tw |\n|---:|---:|---:|---:|")
        for row in domain_rows
            @printf(io, "| %.0f | %d | %.9f | %.9f |\n", row["zmax"],
                    row["turning_point"], row["Hinf"], row["Tw"])
        end
        println(io, "\nThe first fold converges rapidly; the second is more sensitive because small |Hinf| creates a long thermal tail.\n")
        println(io, "## Cross-validation at Tw=1.045\n")
        println(io, "| branch | Hinf | shooting profile error | spectral-grid error | Newton residual |")
        println(io, "|---:|---:|---:|---:|---:|")
        for row in validation
            @printf(io, "| %d | %.9f | %.3e | %.3e | %.3e |\n", row["branch"],
                    row["Hinf"], row["shoot_profile_error"],
                    row["newton_profile_error"], row["newton_residual"])
        end
        println(io, "\nIndependent IVP shooting and two Chebyshev resolutions converge to the same profiles.\n")
        println(io, "## Similarity-preserving temporal stability at Tw=1.045\n")
        println(io, "| branch | Hinf | leading eigenvalue | unstable eigenvalues |")
        println(io, "|---:|---:|---:|---:|")
        for row in fixed_stability
            @printf(io, "| %d | %.9f | %+.6e%+.6ei | %d |\n", row["branch"],
                    row["Hinf"], row["leading_real"], row["leading_imag"],
                    row["unstable_eigenvalues"])
        end
        println(io, "\nThis spectrum is restricted to axisymmetric similarity-preserving disturbances; it is not the full physical stability spectrum.\n")
        println(io, "## Infinite-domain condition\n")
        println(io, "For a non-isothermal solution, T' behaves asymptotically as exp(Pr Hinf eta). Thus Hinf<0 is required for exponential thermal decay. As Hinf approaches zero from below, the upper branch requires progressively larger domains.\n")
        println(io, "## Upper-branch domain stability\n")
        println(io, "| zmax | Hinf | Tw | leading eigenvalue | unstable eigenvalues |")
        println(io, "|---:|---:|---:|---:|---:|")
        for row in upper_rows
            @printf(io, "| %.0f | %.2f | %.9f | %+.6e%+.6ei | %d |\n",
                    row["zmax"], row["Hinf"], row["Tw"], row["leading_real"],
                    row["leading_imag"], row["unstable_eigenvalues"])
        end
        println(io, "\n## Interpretation\n")
        println(io, "The failure near Tw=1.05 is a fold of the similarity solution, not disappearance of every mathematical solution. Heating reduces radial outflow, weakens axial entrainment, thickens the thermal layer and increases integrated inward buoyancy. This positive feedback produces the saddle-node; physical selection remains a separate stability question.")
    end
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}("output-dir" => DEFAULT_OUTPUT, "zmax" => 20.0,
        "degree" => 100, "dh" => 0.005, "target-tw" => 1.045,
        "skip-domain-study" => false)
    args = parse_cli(arguments, defaults)
    output = abspath(args["output-dir"])
    mkpath(output)
    zmax, degree = args["zmax"], args["degree"]
    branch = trace_vk_h_branch(; zmax, degree, h_start=-0.75, h_stop=0.10,
                               dh=args["dh"], tolerance=2e-9)
    turns = vk_turning_points(branch)
    branch_rows = Dict{String,Any}[]
    for profile in branch.profiles
        wall = vk_state(profile, 0.0)
        push!(branch_rows, Dict("Hinf" => profile.hinf, "Tw" => profile.tw,
            "Fp0" => wall[2], "Gp0" => wall[4], "Tp0" => wall[6],
            "nodes" => length(profile.z)))
    end
    write_csv_rows(joinpath(output, "branch_by_Hinf.csv"), branch_rows)

    turning_rows = Dict{String,Any}[]
    for (point_id, item) in enumerate(turns[1:min(2, end)])
        profile = exact_profile(branch, item.x; degree=max(degree, 70), zmax)
        diagnostics = vk_newton_diagnostics(profile)
        push!(turning_rows, Dict("point" => point_id, "Hinf" => item.x,
            "Tw" => profile.tw, "newton_residual" => diagnostics.residual,
            "newton_sigma_min_over_max" => diagnostics.sigma_ratio,
            "thermal_decay_length" => 1 / (PR * abs(item.x))))
    end
    write_csv_rows(joinpath(output, "turning_points.csv"), turning_rows)

    target_tw = args["target-tw"]
    roots = vk_roots_at_tw(branch, target_tw)
    validation = Dict{String,Any}[]
    profile_rows = Dict{String,Any}[]
    selected_profiles = Tuple{Int,VKProfile}[]
    fixed_stability = Dict{String,Any}[]
    z_compare = collect(range(0.0, zmax; length=1001))
    for (branch_id, hinf) in enumerate(roots)
        profile = exact_profile(branch, hinf; degree, zmax)
        slopes = vk_state(profile, 0.0)[[2, 4, 6]]
        answer, shoot_residual, ivp = vk_shoot(target_tw, slopes; zmax, step=0.002)
        reference = vk_state(profile, z_compare)[[1, 3, 5, 7], :]
        shooting = ivp(z_compare)[[1, 3, 5, 7], :]
        low_degree = solve_vk_fixed_h(hinf, profile;
            degree=max(50, degree - 20), zmax, tolerance=2e-9)
        low_values = barycentric_interpolate(low_degree.z, low_degree.weights,
                                              low_degree.fields, profile.z)
        diagnostic = vk_newton_diagnostics(profile)
        push!(validation, Dict("branch" => branch_id, "Hinf" => hinf,
            "Tw" => profile.tw, "Fp0" => slopes[1], "Gp0" => slopes[2],
            "Tp0" => slopes[3], "shoot_bc_residual" => norm(shoot_residual, Inf),
            "shoot_profile_error" => maximum(abs.(shooting .- reference)),
            "newton_residual" => diagnostic.residual,
            "newton_profile_error" => maximum(abs.(low_values .- profile.fields)),
            "newton_sigma_min_over_max" => diagnostic.sigma_ratio,
            "shoot_converged" => answer.converged))
        push!(selected_profiles, (branch_id, profile))
        for index in eachindex(z_compare)
            push!(profile_rows, Dict("branch" => branch_id, "z" => z_compare[index],
                "H" => reference[1, index], "F" => reference[2, index],
                "G" => reference[3, index], "T" => reference[4, index]))
        end
        values = vk_similarity_eigenvalues(profile; degree=70, zmax)
        leading = first(values)
        push!(fixed_stability, Dict("branch" => branch_id, "Hinf" => hinf,
            "Tw" => profile.tw, "leading_real" => real(leading),
            "leading_imag" => imag(leading),
            "unstable_eigenvalues" => count(>(1e-8), real.(values))))
    end
    write_csv_rows(joinpath(output, "cross_validation.csv"), validation)
    write_csv_rows(joinpath(output, "profiles_Tw1p045.csv"), profile_rows)
    write_csv_rows(joinpath(output, "similarity_stability_Tw1p045.csv"), fixed_stability)

    stability_branch = Dict{String,Any}[]
    for hinf in collect(-0.70:0.02:-0.06)
        profile = nearest_profile(branch, hinf)
        values = vk_similarity_eigenvalues(profile; degree=60, zmax)
        leading = first(values)
        push!(stability_branch, Dict("Hinf" => profile.hinf, "Tw" => profile.tw,
            "leading_real" => real(leading), "leading_imag" => imag(leading),
            "unstable_eigenvalues" => count(>(1e-8), real.(values))))
    end
    write_csv_rows(joinpath(output, "similarity_stability_branch.csv"), stability_branch)

    domain_rows = Dict{String,Any}[]
    upper_rows = Dict{String,Any}[]
    if !args["skip-domain-study"]
        for local_zmax in (15.0, 20.0, 25.0, 30.0, 40.0)
            local_turns = domain_fold_estimate(local_zmax;
                                               degree=max(70, round(Int, 2local_zmax + 40)))
            for (index, item) in enumerate(local_turns[1:min(2, end)])
                push!(domain_rows, Dict("zmax" => local_zmax,
                    "turning_point" => index, "Hinf" => item.x, "Tw" => item.y))
            end
        end
        upper_rows = upper_branch_domain_stability()
    end
    write_csv_rows(joinpath(output, "domain_sensitivity.csv"), domain_rows;
                   columns=["zmax", "turning_point", "Hinf", "Tw"])
    write_csv_rows(joinpath(output, "upper_branch_domain_stability.csv"), upper_rows;
                   columns=["zmax", "Hinf", "Tw", "leading_real", "leading_imag",
                            "unstable_eigenvalues", "thermal_decay_length"])

    save_branch_plot(branch, turns, target_tw, output)
    save_profile_plot(selected_profiles, output)
    p = plot([row["Hinf"] for row in stability_branch],
             [row["leading_real"] for row in stability_branch]; marker=:circle,
             xlabel="Hinf", ylabel="leading Re(lambda)", label=false,
             title="Similarity-preserving temporal stability", gridalpha=0.25)
    hline!(p, [0.0]; color=:black, label=false)
    for item in turns[1:min(2, end)]
        vline!(p, [item.x]; color=:crimson, linestyle=:dash, label=false)
    end
    savefig(p, joinpath(output, "similarity_stability_branch.png"))
    savefig(p, joinpath(output, "similarity_stability_branch.pdf"))
    write_report(joinpath(output, "report.md"), turning_rows, domain_rows,
                 validation, fixed_stability, upper_rows)
    println("Output written to $output")
    for item in turns[1:min(2, end)]
        @printf("  Hinf=%.9f, Tw=%.9f\n", item.x, item.y)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
