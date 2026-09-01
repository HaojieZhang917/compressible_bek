#!/usr/bin/env julia

"""Validate and compare von Karman and rotor--stator saddle-node dynamics."""

include(joinpath(@__DIR__, "..", "..", "rotor_stator", "scripts",
                 "analyze_three_solution_dynamics.jl"))
using JSON3
using LinearAlgebra
using Plots
using Printf

const CROSS_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TWODISK = joinpath(CROSS_PROJECT_ROOT, "rotor_stator")
const TWODISK_DATA = joinpath(TWODISK, "data")
const VONKARMAN = joinpath(CROSS_PROJECT_ROOT, "von_karman")
const CROSS_OUTPUT = joinpath(CROSS_PROJECT_ROOT, "cross_model", "data",
                              "dynamical_singularity_comparison")

function first_fold_by_domain()
    rows = read_csv_rows(joinpath(VONKARMAN, "data", "finite_domain_branches",
                                  "turning_points_by_zmax.csv"))
    return Dict(numeric_value(row["zmax"]) =>
                (numeric_value(row["Hinf"]), numeric_value(row["Tw"]))
                for row in rows
                if Int(round(numeric_value(row["turning_point"]))) == 1)
end

function continue_vk(zmax::Real, target_h_values; tolerance::Float64=2e-8,
                     degree::Int=max(100, round(Int, 2zmax + 60)))
    targets = sort(unique(Float64.(target_h_values)))
    hmax = maximum(targets)
    stepping = collect(-0.75:0.01:(hmax + 0.0101))
    values = filter(h -> -0.7500001 <= h <= hmax + 1e-12,
                    sort(unique(vcat(stepping, targets))))
    seed = solve_vk_isothermal(; zmax=Float64(zmax), degree, tolerance=2e-9)
    output = Dict{Float64,VKProfile}()
    for hinf in values
        seed = solve_vk_fixed_h(hinf, seed; zmax=Float64(zmax), degree, tolerance)
        if any(abs.(targets .- hinf) .< 2e-12)
            output[round(hinf; digits=12)] = seed
        end
    end
    return output
end

function vk_similarity_spectrum(profile::VKProfile; zmax::Float64,
                                degree::Int)
    z, D, D2, _ = chebyshev_operators(degree, 0.0, zmax)
    fields = barycentric_interpolate(profile.z, profile.weights, profile.fields, z)
    H, F, G, T = (view(fields, row, :) for row in 1:4)
    Fp, Gp, Tp = D * F, D * G, D * T
    nodes = length(z)
    rf, rg, rh, rt = ntuple(i -> ((i - 1) * nodes + 1):(i * nodes), 4)
    identity = Matrix{Float64}(I, nodes, nodes)
    A, B = zeros(4nodes, 4nodes), zeros(4nodes, 4nodes)
    A[rf, rf] .= D2 - Diagonal(H) * D - Diagonal(2F)
    A[rf, rg] .= Diagonal(2 .* (G .- 1)); A[rf, rh] .= -Diagonal(Fp)
    A[rf, rt] .= -identity
    A[rg, rf] .= Diagonal(2 .* (1 .- G))
    A[rg, rg] .= D2 - Diagonal(H) * D - Diagonal(2F)
    A[rg, rh] .= -Diagonal(Gp)
    A[rh, rf] .= 2identity; A[rh, rh] .= D
    A[rt, rh] .= -Diagonal(Tp); A[rt, rt] .= D2 / 0.72 - Diagonal(H) * D
    B[rf, rf] .= identity; B[rg, rg] .= identity; B[rt, rt] .= identity
    function set_bc(row, column)
        A[row, :] .= 0; B[row, :] .= 0; A[row, column] = 1
    end
    set_bc(first(rf), first(rf)); set_bc(last(rf), last(rf))
    set_bc(first(rg), first(rg)); set_bc(last(rg), last(rg))
    set_bc(first(rh), first(rh)); set_bc(first(rt), first(rt)); set_bc(last(rt), last(rt))
    spectrum = eigen(A, B)
    keep = isfinite.(real.(spectrum.values)) .& isfinite.(imag.(spectrum.values)) .&
           (abs.(spectrum.values) .< 1e5)
    values = ComplexF64.(spectrum.values[keep])
    vectors = ComplexF64.(spectrum.vectors[:, keep])
    order = sortperm(eachindex(values); by=i -> (-real(values[i]), abs(imag(values[i]))))
    return z, values[order], vectors[:, order]
end

leading_vk(profile; zmax, degree) = first(vk_similarity_spectrum(profile; zmax, degree)[2])

function through_origin_fit(x, y)
    slope = dot(x, y) / dot(x, x)
    residual = y .- slope .* x
    return slope, 1 - dot(residual, residual) / dot(y, y)
end

function vk_fold_scan(zmax, h_fold, tw_fold; degree::Int=130, quick::Bool=false)
    offsets = quick ? [-0.04, -0.01, -0.0025, 0.0, 0.0025, 0.01, 0.04] :
        [-0.12, -0.10, -0.08, -0.06, -0.04, -0.03, -0.02, -0.015,
         -0.010, -0.0075, -0.005, -0.0025, 0.0, 0.0025, 0.005,
         0.0075, 0.010, 0.015, 0.020, 0.030, 0.040, 0.060, 0.080]
    targets = h_fold .+ offsets
    solutions = continue_vk(zmax, targets; tolerance=8e-9,
                            degree=quick ? 80 : max(120, degree))
    rows = Dict{String,Any}[]
    for hinf in targets
        solution = solutions[round(hinf; digits=12)]
        leading = leading_vk(solution; zmax=Float64(zmax), degree)
        push!(rows, Dict("Hinf" => hinf, "Tw" => solution.tw,
            "leading_real" => real(leading), "leading_imag" => imag(leading),
            "unstable_count_class" => Int(real(leading) >= 0),
            "fold_sample" => Int(abs(hinf - h_fold) < 1e-12)))
    end
    fits = Dict{String,Any}[]
    for (name, side) in (("benchmark_connected_stable", (<)(h_fold)),
                         ("post_fold_unstable", (>)(h_fold)))
        selected = filter(row -> side(row["Hinf"]), rows)
        mu = [tw_fold - row["Tw"] for row in selected]
        lam2 = [row["leading_real"]^2 for row in selected]
        use = (mu .> 1e-9) .& (mu .< 2e-5)
        if count(use) >= 1
            slope, r2 = through_origin_fit(mu[use], lam2[use])
            push!(fits, Dict("model" => "von_Karman", "fold" => "principal",
                "side" => name, "points" => count(use),
                "temperature_distance_max" => 2e-5,
                "lambda_squared_slope" => slope, "r2_through_origin" => r2))
        end
    end
    return rows, fits, solutions[round(h_fold; digits=12)]
end

function steady_jacobian_checks(profile::VKProfile, tw_fold::Real)
    state = BoussinesqSaddleNode.pack_fields(profile.fields)
    residual, jacobian = BoussinesqSaddleNode._vk_fixed_tw_system(
        state, Float64(tw_fold), profile.D, profile.D * profile.D)
    row_norm = [norm(view(jacobian, row, :)) for row in axes(jacobian, 1)]
    scale = 1.0 ./ max.(row_norm, 1e-300)
    scaled = scale .* jacobian
    factor = svd(scaled; full=true)
    psi, phi = factor.U[:, end], factor.V[:, end]
    delta_tw = 1e-6
    plus = first(BoussinesqSaddleNode._vk_fixed_tw_system(
        state, Float64(tw_fold + delta_tw), profile.D, profile.D * profile.D))
    minus = first(BoussinesqSaddleNode._vk_fixed_tw_system(
        state, Float64(tw_fold - delta_tw), profile.D, profile.D * profile.D))
    parameter_direction = scale .* (plus .- minus) ./ (2delta_tw)
    transversality = dot(psi, parameter_direction)
    coefficients = Dict{String,Any}[]
    for epsilon in (1e-3, 1e-4, 1e-5)
        _, jplus = BoussinesqSaddleNode._vk_fixed_tw_system(
            state .+ epsilon .* phi, Float64(tw_fold), profile.D, profile.D * profile.D)
        _, jminus = BoussinesqSaddleNode._vk_fixed_tw_system(
            state .- epsilon .* phi, Float64(tw_fold), profile.D, profile.D * profile.D)
        second = scale .* (((jplus .- jminus) ./ (2epsilon)) * phi)
        push!(coefficients, Dict("epsilon" => epsilon,
            "transversality_coefficient" => transversality,
            "quadratic_coefficient" => 0.5dot(psi, second)))
    end
    checks = Dict("degree" => length(profile.z) - 1,
        "newton_residual" => norm(residual, Inf), "newton_iterations" => missing,
        "sigma_min_over_max" => factor.S[end] / factor.S[1],
        "sigma_second_min_over_max" => factor.S[end-1] / factor.S[1],
        "nullity_gap_sigma_second_over_sigma_min" => factor.S[end-1] / factor.S[end],
        "transversality_coefficient" => transversality)
    nodes = length(profile.z)
    null_fields = permutedims(reshape(phi, nodes, 4))
    return null_fields, checks, coefficients
end

trapz_cross(values, z) = sum((values[1:end-1] .+ values[2:end]) .* diff(z)) / 2

function temporal_mode_metrics(z, vector)
    modes = ndims(vector) == 1 ? permutedims(reshape(vector, length(z), 4)) : vector
    maxima = [maximum(abs, view(modes, row, :)) for row in 1:4]
    scale = maximum(maxima)
    normalised = modes ./ scale
    amplitude_sq = vec(sum(abs2.(normalised); dims=1))
    active_sq = vec(sum(abs2.(normalised[[1, 2, 4], :]); dims=1))
    peaks = [z[argmax(abs.(view(normalised, row, :)))] for row in 1:4]
    metrics = Dict{String,Any}(
        "F_relative_max" => maxima[1] / scale, "G_relative_max" => maxima[2] / scale,
        "H_relative_max" => maxima[3] / scale, "T_relative_max" => maxima[4] / scale,
        "coordinate_centroid" => trapz_cross(z .* amplitude_sq, z) / trapz_cross(amplitude_sq, z),
        "active_fields_centroid" => trapz_cross(z .* active_sq, z) / trapz_cross(active_sq, z),
        "F_peak_coordinate" => peaks[1], "G_peak_coordinate" => peaks[2],
        "H_peak_coordinate" => peaks[3], "T_peak_coordinate" => peaks[4])
    return metrics, normalised
end

function read_rotor_folds()
    path = joinpath(TWODISK_DATA, "boussinesq_singularity_results",
                    "folds_Re1000_traditional_centrifugal.json")
    return JSON3.read(read(path, String))
end

function rotor_branch_at_pressures(target_pressures; degree::Int=100)
    config = RotorConfig(re_h=1000.0, pr=0.72, tolerance=2e-7,
                         collocation_degree=degree,
                         model="traditional_centrifugal")
    reference = joinpath(TWODISK, "reference", "baseflow_Res1000.npz")
    isothermal = solve_rotor_isothermal(config, reference)
    coarse = collect(isothermal.pressure:-0.002:0.015)
    moderate = collect((last(coarse) - 0.001):-0.001:0.011)
    fine = collect((last(moderate) - 0.0005):-0.0005:-0.001)
    path = sort(unique(vcat(coarse, moderate, fine, target_pressures)); rev=true)
    return continue_rotor_pressure(path, isothermal; profile_points=1001)
end

function base_profile(profile::RotorProfile)
    z = collect(range(0.0, 1.0; length=1001))
    state = rotor_state(profile, z)
    return BaseProfile(0, profile.pressure, z, vec(state[1, :]), vec(state[2, :]),
                       vec(state[4, :]), vec(state[6, :]))
end

function scan_rotor_stator_near_folds(; order::Int=101, quick::Bool=false)
    folds = read_rotor_folds()
    offsets = Dict(
        "maximum" => (quick ? [-0.00035, -0.00008, 0.0, 0.00008, 0.00035] :
            [-0.0012, -0.0009, -0.0007, -0.0005, -0.00035, -0.00024,
             -0.00016, -0.00008, 0.0, 0.00008, 0.00016, 0.00024,
             0.00035, 0.0005, 0.0007, 0.0009, 0.0012]),
        "minimum" => (quick ? [-0.00015, -0.00003, 0.0, 0.00003, 0.00015] :
            [-0.00070, -0.00055, -0.00040, -0.00030, -0.00022, -0.00015,
             -0.00010, -0.00005, -0.00003, -0.00002, -0.00001, 0.0,
             0.00001, 0.00002, 0.00003, 0.00005, 0.00010, 0.00015,
             0.00022, 0.00030, 0.00040, 0.00055, 0.00070]))
    metadata = Tuple{Float64,String,Float64}[]
    for fold in folds
        pressure = Float64(fold.pressure_gradient)
        for offset in offsets[String(fold.kind)]
            push!(metadata, (pressure + offset, String(fold.kind), offset))
        end
    end
    branch = rotor_branch_at_pressures(first.(metadata); degree=quick ? 70 : 100)
    branch_pressures = branch.columns["pressure_gradient"]
    rows = Dict{String,Any}[]
    for (pressure, kind, offset) in metadata
        index = argmin(abs.(branch_pressures .- pressure))
        profile = base_profile(branch.profiles[index])
        _, temporal = temporal_operator(profile, order, 1000.0, 0.72)
        leading = first(temporal.eigenvalues)
        push!(rows, Dict("fold" => kind, "pressure_gradient" => pressure,
            "pressure_offset" => offset, "Tw" => branch.profiles[index].tw,
            "leading_real" => real(leading), "leading_imag" => imag(leading),
            "positive_real_count" => count(>(1e-7), real.(temporal.eigenvalues)),
            "fold_sample" => Int(abs(offset) < 1e-15), "collocation_order" => order))
    end
    sort!(rows; by=row -> (row["fold"], -row["pressure_gradient"]))
    return rows
end

function two_disk_fold_fits(growth)
    rows = Dict{String,Any}[]
    for fold in read_rotor_folds()
        kind = String(fold.kind); twc = Float64(fold.Tw); pic = Float64(fold.pressure_gradient)
        is_maximum = kind == "maximum"; cap = is_maximum ? 1e-4 : 5e-5
        for (name, test) in (("higher_Pi", pi -> pi > pic), ("lower_Pi", pi -> pi < pic))
            selected = Tuple{Float64,Float64}[]
            for item in growth
                item["fold"] == kind || continue
                mu = is_maximum ? twc - item["Tw"] : item["Tw"] - twc
                test(item["pressure_gradient"]) && 1e-9 < mu < cap &&
                    push!(selected, (mu, item["leading_real"]^2))
            end
            isempty(selected) && continue
            slope, r2 = through_origin_fit(first.(selected), last.(selected))
            push!(rows, Dict("model" => "rotor_stator", "fold" => kind,
                "side" => name, "points" => length(selected),
                "temperature_distance_max" => cap,
                "lambda_squared_slope" => slope, "r2_through_origin" => r2))
        end
    end
    return rows
end

function rotor_stator_principal_fold_mode(; order::Int=121, quick::Bool=false)
    fold = only(filter(item -> String(item.kind) == "maximum", read_rotor_folds()))
    pressure = Float64(fold.pressure_gradient)
    branch = rotor_branch_at_pressures([pressure]; degree=quick ? 70 : 100)
    index = argmin(abs.(branch.columns["pressure_gradient"] .- pressure))
    profile = base_profile(branch.profiles[index])
    _, result = temporal_operator(profile, order, 1000.0, 0.72)
    vector = result.eigenvectors[:, 1]
    na, ng = result.n_f_reduced, result.n_g
    modes = [transpose(result.f_map * vector[1:na]);
             transpose(result.g_map * vector[na+1:na+ng]);
             transpose(result.h_map * vector[1:na]);
             transpose(result.temperature_map * vector[na+ng+1:end])]
    metrics, normalised = temporal_mode_metrics(result.z, modes)
    metrics["leading_real"] = real(first(result.eigenvalues))
    metrics["leading_imag"] = imag(first(result.eigenvalues))
    metrics["Tw"] = branch.profiles[index].tw
    metrics["pressure_gradient"] = pressure
    return result.z, normalised, metrics
end

function plot_vk_fold(rows, tw_fold, h_fold, output)
    h = [row["Hinf"] for row in rows]; tw = [row["Tw"] for row in rows]
    growth = [row["leading_real"] for row in rows]; stable = growth .< 0
    p1 = plot(tw, h; color=:grey, label=false, xlabel="Tw", ylabel="Hinf",
              title="Fold and stability exchange")
    scatter!(p1, tw[stable], h[stable]; color=:green, label="stable")
    scatter!(p1, tw[.!stable], h[.!stable]; color=:red, label="unstable")
    scatter!(p1, [tw_fold], [h_fold]; marker=:star5, color=:black, label="fold")
    p2 = plot(h, growth; marker=:circle, label=false, xlabel="Hinf",
              ylabel="max Re(lambda)", title="Real eigenvalue crossing")
    hline!(p2, [0.0]; color=:black, label=false)
    mu = tw_fold .- tw
    p3 = scatter(mu[stable], growth[stable].^2; color=:green, label="stable",
                 xlabel="Tw,c-Tw", ylabel="[Re(lambda)]^2",
                 title="Square-root scaling")
    scatter!(p3, mu[.!stable], growth[.!stable].^2; color=:red, label="unstable")
    savefig(plot(p1, p2, p3; layout=(1, 3), size=(1600, 470)), output)
end

function plot_zero_modes(z_vk, vk_mode, z_td, td_mode, output)
    p1 = plot(; xlabel="eta", ylabel="normalised amplitude",
              title="von Karman principal-fold zero mode", xlim=(0, min(20, last(z_vk))))
    p2 = plot(; xlabel="z/h", ylabel="normalised amplitude",
              title="rotor-stator principal-fold zero mode")
    for (row, label) in enumerate(("F", "G", "H", "temperature"))
        plot!(p1, z_vk, abs.(vk_mode[row, :]); label)
        plot!(p2, z_td, abs.(td_mode[row, :]); label)
    end
    savefig(plot(p1, p2; layout=(1, 2), size=(1250, 480)), output)
end

function plot_normal_form_scaling(vk_rows, td_rows, output)
    fold_map = Dict(String(item.kind) => item for item in read_rotor_folds())
    panels = Any[]
    for (model, rows, kind) in (("von Karman", vk_rows, "principal"),
                                ("rotor-stator", td_rows, "minimum"),
                                ("rotor-stator", td_rows, "maximum"))
        if model == "von Karman"
            selected = rows; twc = maximum(row["Tw"] for row in rows); cap = 2e-5
            mu = [twc - row["Tw"] for row in selected]
        else
            selected = filter(row -> row["fold"] == kind, rows)
            twc = Float64(fold_map[kind].Tw); cap = kind == "minimum" ? 5e-5 : 1e-4
            mu = kind == "minimum" ? [row["Tw"] - twc for row in selected] :
                                      [twc - row["Tw"] for row in selected]
        end
        lam2 = [row["leading_real"]^2 for row in selected]
        stable = [row["leading_real"] <= 0 for row in selected]
        panel = scatter(mu[stable], lam2[stable]; color=:green, label="stable",
                        xlabel="temperature distance from fold",
                        ylabel="[Re(lambda)]^2", title="$model $kind fold")
        scatter!(panel, mu[.!stable], lam2[.!stable]; color=:red, label="unstable")
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(1, 3), size=(1580, 460)), output)
end

function main(arguments=ARGS)
    args = parse_cli(arguments, Dict{String,Any}(
        "output-dir" => CROSS_OUTPUT, "quick" => false))
    output = abspath(args["output-dir"]); mkpath(output)
    quick = args["quick"]
    fold_domains = first_fold_by_domain()
    zmax = 60.0
    h_fold, tw_fold = fold_domains[zmax]
    vk_degree = quick ? 70 : 130
    near_rows, vk_fits, fold_profile = vk_fold_scan(zmax, h_fold, tw_fold;
                                                    degree=vk_degree, quick)
    write_csv_rows(joinpath(output, "vonkarman_near_fold_temporal.csv"), near_rows)
    plot_vk_fold(near_rows, tw_fold, h_fold,
                 joinpath(output, "vonkarman_fold_dynamics.png"))

    convergence = Dict{String,Any}[]
    domains = quick ? (60.0,) : (20.0, 40.0, 60.0, 80.0)
    for domain in domains
        hc, twc = fold_domains[domain]
        solve_degree = quick ? 80 : max(100, round(Int, 2domain + 60))
        solution = continue_vk(domain, [hc]; tolerance=1e-8,
                               degree=solve_degree)[round(hc; digits=12)]
        degrees = quick ? (60, 70) : Dict(20.0 => (60, 80, 100),
            40.0 => (80, 110, 140), 60.0 => (100, 130, 160),
            80.0 => (120, 150, 180))[domain]
        for degree in degrees
            leading = leading_vk(solution; zmax=domain, degree)
            push!(convergence, Dict("zmax" => domain, "degree" => degree,
                "Hinf_fold" => hc, "Tw_fold" => twc,
                "leading_real" => real(leading), "leading_imag" => imag(leading)))
        end
    end
    write_csv_rows(joinpath(output, "vonkarman_fold_spectrum_convergence.csv"), convergence)

    steady_null, jacobian_checks, coefficients = steady_jacobian_checks(fold_profile, tw_fold)
    write_csv_rows(joinpath(output, "vonkarman_nondegeneracy_coefficients.csv"), coefficients)
    write_json(joinpath(output, "vonkarman_jacobian_checks.json"), jacobian_checks)
    z_vk, vk_values, vk_vectors = vk_similarity_spectrum(fold_profile;
        zmax, degree=quick ? 70 : 120)
    vk_metrics, vk_mode = temporal_mode_metrics(z_vk, vk_vectors[:, 1])
    vk_metrics["leading_real"] = real(first(vk_values)); vk_metrics["leading_imag"] = imag(first(vk_values))
    vk_metrics["Tw"] = tw_fold; vk_metrics["Hinf"] = h_fold
    steady_on_temporal_grid = barycentric_interpolate(
        fold_profile.z, fold_profile.weights, steady_null, z_vk)
    steady_reordered = steady_on_temporal_grid[[2, 3, 1, 4], :]
    steady_reordered ./= maximum(abs, steady_reordered)
    correlation = abs(dot(vec(steady_reordered), vec(real.(vk_mode)))) /
                  (norm(steady_reordered) * norm(real.(vk_mode)))
    vk_metrics["steady_temporal_zero_mode_correlation"] = correlation
    write_json(joinpath(output, "vonkarman_zero_mode_metrics.json"), vk_metrics)
    write_csv_rows(joinpath(output, "vonkarman_principal_fold_zero_mode.csv"),
        [Dict("eta" => z_vk[i], "F_mode_abs" => abs(vk_mode[1, i]),
              "G_mode_abs" => abs(vk_mode[2, i]), "H_mode_abs" => abs(vk_mode[3, i]),
              "T_mode_abs" => abs(vk_mode[4, i])) for i in eachindex(z_vk)])

    td_z, td_mode, td_metrics = rotor_stator_principal_fold_mode(;
        order=quick ? 61 : 121, quick)
    write_json(joinpath(output, "rotor_stator_zero_mode_metrics.json"), td_metrics)
    write_csv_rows(joinpath(output, "rotor_stator_principal_fold_zero_mode.csv"),
        [Dict("z_over_h" => td_z[i], "F_mode_abs" => abs(td_mode[1, i]),
              "G_mode_abs" => abs(td_mode[2, i]), "H_mode_abs" => abs(td_mode[3, i]),
              "T_mode_abs" => abs(td_mode[4, i])) for i in eachindex(td_z)])
    plot_zero_modes(z_vk, vk_mode, td_z, td_mode,
                    joinpath(output, "principal_fold_zero_modes.png"))

    td_rows = scan_rotor_stator_near_folds(; order=quick ? 61 : 101, quick)
    write_csv_rows(joinpath(output, "rotor_stator_near_fold_temporal.csv"), td_rows)
    all_fits = vcat(vk_fits, two_disk_fold_fits(td_rows))
    write_csv_rows(joinpath(output, "normal_form_square_root_fits.csv"), all_fits)
    plot_normal_form_scaling(near_rows, td_rows,
                             joinpath(output, "both_models_square_root_scaling.png"))
    comparison = [
        Dict("model" => "von_Karman", "domain" => "eta_max_60",
             "principal_fold_control" => "Hinf", "principal_fold_coordinate" => h_fold,
             "principal_fold_temperature" => tw_fold,
             "fold_temporal_real" => real(first(vk_values)),
             "fold_temporal_imag" => imag(first(vk_values)),
             "zero_mode_centroid" => vk_metrics["coordinate_centroid"],
             "active_fields_centroid" => vk_metrics["active_fields_centroid"],
             "global_pressure_parameter" => 0),
        Dict("model" => "rotor_stator", "domain" => "finite_gap_Re_h_1000",
             "principal_fold_control" => "Pi",
             "principal_fold_coordinate" => td_metrics["pressure_gradient"],
             "principal_fold_temperature" => td_metrics["Tw"],
             "fold_temporal_real" => td_metrics["leading_real"],
             "fold_temporal_imag" => td_metrics["leading_imag"],
             "zero_mode_centroid" => td_metrics["coordinate_centroid"],
             "active_fields_centroid" => td_metrics["active_fields_centroid"],
             "global_pressure_parameter" => 1)]
    write_csv_rows(joinpath(output, "model_dynamical_comparison.csv"), comparison)
    JSON3.pretty(stdout, Dict("vonkarman_fold" => Dict("Hinf" => h_fold,
        "Tw" => tw_fold, "lambda" => [real(first(vk_values)), imag(first(vk_values))],
        "mode_metrics" => vk_metrics, "jacobian" => jacobian_checks,
        "nondegeneracy" => coefficients), "rotor_stator_fold" => td_metrics,
        "normal_form_fits" => all_fits); allow_inf=true)
    println()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
