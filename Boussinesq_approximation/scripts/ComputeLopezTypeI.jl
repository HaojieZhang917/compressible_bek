const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(WORKSPACE_ROOT, "NeutralCurveRunner.jl"))

using .NeutralCurveRunner
using DelimitedFiles
using Printf

const DEFAULT_TYPEI_SEEDS = Dict(
    1.08 => (beta=0.1202, alpha=0.6630),
)

function typeI_seed(Tw::Real)
    key = round(Float64(Tw); digits=2)
    haskey(DEFAULT_TYPEI_SEEDS, key) || throw(ArgumentError(
        "No default Type-I seed is available for Tw=$Tw",
    ))
    return DEFAULT_TYPEI_SEEDS[key]
end

function typeI_config(
    Tw::Real, output_dir::AbstractString;
    N_cheb::Integer=69, beta_step::Real=-8.0e-4,
)
    seed = typeI_seed(Tw)
    return CurveConfig(
        Tw=Float64(Tw), omega=0.0, R_initial=500.0,
        beta_initial=seed.beta, alpha_target=seed.alpha,
        num_modes=2, model=:lopez, Mr=0.3, Ro=-1.0,
        N_cheb=Int(N_cheb), beta_step=Float64(beta_step),
        beta_scan_step=2.0e-4, beta_bounds=(0.04, 0.16),
        min_beta_step=2.5e-5, step_recovery_successes=4,
        neutral_tol=1.0e-7, R_tol=1.0e-4,
        corrector_R_step=0.5, max_scan_steps=160,
        max_prediction_step=20.0, min_mode_overlap=0.50,
        max_curve_points=200, min_valid_points=8,
        minimum_complete_points=50, minimum_complete_beta=0.08,
        output_dir=String(output_dir), keep_logs=true,
    )
end

function neutral_columns(data::AbstractMatrix)
    return [abs(data[index, 5]) <= abs(data[index, 7]) ? 1 : 2
            for index in axes(data, 1)]
end

function typeI_turning_issues(data::AbstractMatrix)
    issues = String[]
    minimum_index = argmin(data[:, 2])
    minimum_index > 1 || push!(issues, "Type-I minimum R is at the first point")
    minimum_index < size(data, 1) || push!(
        issues, "Type-I minimum R is at the last point",
    )
    if 1 < minimum_index < size(data, 1)
        all(diff(data[begin:minimum_index, 2]) .< 0) || push!(
            issues, "R is not decreasing before the Type-I critical point",
        )
        all(diff(data[minimum_index:end, 2]) .> 0) || push!(
            issues, "R is not increasing after the Type-I critical point",
        )
    end
    return issues
end

function validate_typeI(result, config::CurveConfig)
    data = result.data
    issues = copy(result.validation.issues)
    size(data, 1) >= 50 || push!(issues, "fewer than 50 Type-I points")
    maximum(data[:, 3]) >= 0.11 || push!(issues, "branch did not reach beta=0.11")
    minimum(data[:, 3]) <= 0.055 || push!(issues, "branch did not reach beta=0.055")
    abs(data[1, 2] - config.R_initial) <= 1.0e-3 || push!(
        issues, "first point is not the R=500 neutral crossing",
    )
    residuals = min.(abs.(data[:, 5]), abs.(data[:, 7]))
    maximum(residuals) <= config.neutral_tol || push!(
        issues, "neutral residual exceeds $(config.neutral_tol)",
    )
    columns = neutral_columns(data)
    switch_count = count(!=(0), diff(columns))
    switch_count == 0 || push!(
        issues, "active neutral eigenvalue column switched $switch_count times",
    )
    append!(issues, typeI_turning_issues(data))
    return (ok=isempty(issues), issues=issues, columns=columns)
end

function compute_lopez_typeI(
    Tw::Real;
    output_dir::AbstractString=joinpath(WORKSPACE_ROOT, "neutral_curve_batch"),
    N_cheb::Integer=69,
)
    mkpath(output_dir)
    return mktempdir() do work_dir
        config = typeI_config(Tw, work_dir; N_cheb=N_cheb)
        result = compute_neutral_curve(config)
        check = validate_typeI(result, config)
        check.ok || error(
            "Type-I validation failed for Tw=$Tw: " * join(check.issues, "; "),
        )

        destination = joinpath(
            output_dir,
            NeutralCurveRunner.case_tag(config) * ".dat",
        )
        cp(result.path, destination; force=true)
        @printf(
            "Tw=%.2f Type-I: points=%d beta=[%.8f, %.8f] R=[%.8f, %.8f] max_residual=%.3e stop=%s\n",
            Tw, size(result.data, 1), minimum(result.data[:, 3]),
            maximum(result.data[:, 3]), minimum(result.data[:, 2]),
            maximum(result.data[:, 2]),
            maximum(min.(abs.(result.data[:, 5]), abs.(result.data[:, 7]))),
            String(result.stop_reason),
        )
        println("saved $destination")
        return merge(result, (path=destination, typeI_validation=check))
    end
end

function resume_lopez_typeI(
    Tw::Real, prefix_path::AbstractString;
    output_dir::AbstractString=joinpath(WORKSPACE_ROOT, "neutral_curve_batch"),
    N_cheb::Integer=69,
)
    prefix = NeutralCurveRunner.load_curve(prefix_path)
    size(prefix, 1) >= 2 || error("Type-I prefix contains fewer than two points")
    last_row = prefix[end, :]
    active_column = abs(last_row[5]) <= abs(last_row[7]) ? 1 : 2
    alpha_target = active_column == 1 ? last_row[4] : last_row[6]

    return mktempdir() do work_dir
        config = CurveConfig(
            Tw=Float64(Tw), omega=0.0, R_initial=last_row[2],
            beta_initial=last_row[3], alpha_target=alpha_target,
            num_modes=2, model=:lopez, Mr=0.3, Ro=-1.0,
            N_cheb=Int(N_cheb), beta_step=-4.0e-4,
            beta_scan_step=1.0e-4, beta_bounds=(0.04, 0.08),
            min_beta_step=2.5e-5, step_recovery_successes=4,
            neutral_tol=1.0e-7, R_tol=1.0e-4,
            corrector_R_step=0.5, max_scan_steps=160,
            max_prediction_step=15.0, min_mode_overlap=0.50,
            max_curve_points=100, min_valid_points=1,
            minimum_complete_points=1, minimum_complete_beta=0.055,
            output_dir=work_dir, keep_logs=true,
        )
        tail = compute_neutral_curve(config)
        combined = vcat(prefix, tail.data[2:end, :])
        final_config = typeI_config(Tw, output_dir; N_cheb=N_cheb)
        combined_result = (
            data=combined, validation=NeutralCurveRunner.validate_curve_data(
                combined, final_config,
            ),
        )
        check = validate_typeI(combined_result, final_config)
        check.ok || error(
            "Resumed Type-I validation failed for Tw=$Tw: " *
            join(check.issues, "; "),
        )
        all(diff(combined[:, 3]) .< 0) || error(
            "Resumed Type-I beta is not strictly decreasing",
        )

        destination = joinpath(
            output_dir,
            NeutralCurveRunner.case_tag(final_config) * ".dat",
        )
        NeutralCurveRunner.write_curve(destination, final_config, combined)
        residual = maximum(min.(abs.(combined[:, 5]), abs.(combined[:, 7])))
        @printf(
            "Tw=%.2f resumed Type-I: prefix=%d tail=%d total=%d beta=[%.8f, %.8f] R=[%.8f, %.8f] max_residual=%.3e stop=%s\n",
            Tw, size(prefix, 1), size(tail.data, 1), size(combined, 1),
            minimum(combined[:, 3]), maximum(combined[:, 3]),
            minimum(combined[:, 2]), maximum(combined[:, 2]), residual,
            String(tail.stop_reason),
        )
        println("saved $destination")
        return (
            data=combined, path=destination, tail=tail,
            typeI_validation=check,
        )
    end
end

function truncate_after_mode_jump(data::AbstractMatrix)
    minimum_index = argmin(data[:, 2])
    for index in minimum_index + 1:size(data, 1)
        if data[index, 2] <= data[index - 1, 2]
            return Matrix(data[begin:index - 1, :]), index
        end
    end
    return Matrix(data), nothing
end

function refine_lopez_typeI_endpoint(
    Tw::Real, curve_path::AbstractString;
    output_dir::AbstractString=joinpath(WORKSPACE_ROOT, "neutral_curve_batch"),
    N_cheb::Integer=69,
)
    raw = NeutralCurveRunner.load_curve(curve_path)
    prefix, jump_index = truncate_after_mode_jump(raw)
    jump_index === nothing && error("No Type-I mode jump was detected")
    size(prefix, 1) >= 3 || error("Type-I prefix is too short to refine")
    prefix = Matrix(prefix[begin:end - 1, :])
    seed = prefix[end, :]
    active_column = abs(seed[5]) <= abs(seed[7]) ? 1 : 2
    alpha_target = active_column == 1 ? seed[4] : seed[6]

    return mktempdir() do work_dir
        config = CurveConfig(
            Tw=Float64(Tw), omega=0.0, R_initial=seed[2],
            beta_initial=seed[3], alpha_target=alpha_target,
            num_modes=2, model=:lopez, Mr=0.3, Ro=-1.0,
            N_cheb=Int(N_cheb), beta_step=-1.0e-4,
            beta_scan_step=2.5e-5, beta_bounds=(0.045, 0.06),
            min_beta_step=1.25e-5, step_recovery_successes=4,
            neutral_tol=1.0e-7, R_tol=1.0e-4,
            corrector_R_step=0.25, max_scan_steps=200,
            max_prediction_step=8.0, min_mode_overlap=0.80,
            max_curve_points=100, min_valid_points=1,
            minimum_complete_points=1, minimum_complete_beta=0.05,
            output_dir=work_dir, keep_logs=true,
        )
        tail = compute_neutral_curve(config)
        candidate = vcat(prefix, tail.data[2:end, :])
        combined, second_jump = truncate_after_mode_jump(candidate)
        second_jump === nothing || @printf(
            "discarded %d post-jump endpoint rows\n",
            size(candidate, 1) - size(combined, 1),
        )
        final_config = typeI_config(Tw, output_dir; N_cheb=N_cheb)
        combined_result = (
            data=combined, validation=NeutralCurveRunner.validate_curve_data(
                combined, final_config,
            ),
        )
        check = validate_typeI(combined_result, final_config)
        check.ok || error(
            "Refined Type-I validation failed for Tw=$Tw: " *
            join(check.issues, "; "),
        )

        destination = joinpath(
            output_dir,
            NeutralCurveRunner.case_tag(final_config) * ".dat",
        )
        NeutralCurveRunner.write_curve(destination, final_config, combined)
        residual = maximum(min.(abs.(combined[:, 5]), abs.(combined[:, 7])))
        @printf(
            "Tw=%.2f refined Type-I: coarse_prefix=%d refined_tail=%d total=%d beta=[%.8f, %.8f] R=[%.8f, %.8f] max_residual=%.3e stop=%s\n",
            Tw, size(prefix, 1), size(tail.data, 1), size(combined, 1),
            minimum(combined[:, 3]), maximum(combined[:, 3]),
            minimum(combined[:, 2]), maximum(combined[:, 2]), residual,
            String(tail.stop_reason),
        )
        println("saved $destination")
        return (
            data=combined, path=destination, tail=tail,
            typeI_validation=check,
        )
    end
end

function main(args=ARGS)
    if !isempty(args) && first(args) == "--resume"
        length(args) == 3 || throw(ArgumentError(
            "Usage: ComputeLopezTypeI.jl --resume Tw prefix_path",
        ))
        return resume_lopez_typeI(parse(Float64, args[2]), args[3])
    end
    if !isempty(args) && first(args) == "--refine-endpoint"
        length(args) == 3 || throw(ArgumentError(
            "Usage: ComputeLopezTypeI.jl --refine-endpoint Tw curve_path",
        ))
        return refine_lopez_typeI_endpoint(
            parse(Float64, args[2]), args[3],
        )
    end
    temperatures = isempty(args) ? (1.08,) : parse.(Float64, args)
    return [compute_lopez_typeI(Tw) for Tw in temperatures]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
