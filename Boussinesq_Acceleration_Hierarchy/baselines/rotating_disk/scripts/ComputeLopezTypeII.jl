const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(WORKSPACE_ROOT, "NeutralCurveRunner.jl"))

using .NeutralCurveRunner
using Printf

const DEFAULT_TYPEII_SEEDS = Dict(
    1.08 => (beta=0.0278, alpha=0.0715),
    1.18 => (beta=0.0191, alpha=0.0530),
    1.20 => (beta=0.0178, alpha=0.0500),
)

function typeII_seed(Tw::Real)
    key = round(Float64(Tw); digits=2)
    haskey(DEFAULT_TYPEII_SEEDS, key) || throw(ArgumentError(
        "No default Type-II seed is available for Tw=$Tw",
    ))
    return DEFAULT_TYPEII_SEEDS[key]
end

function typeII_config(
    Tw::Real, output_dir::AbstractString;
    N_cheb::Integer=69, beta_step::Real=4.0e-4,
)
    seed = typeII_seed(Tw)
    return CurveConfig(
        Tw=Float64(Tw), omega=0.0, R_initial=500.0,
        beta_initial=seed.beta, alpha_target=seed.alpha,
        num_modes=2, model=:lopez, Mr=0.3, Ro=-1.0,
        N_cheb=Int(N_cheb), beta_step=Float64(beta_step),
        beta_scan_step=1.0e-4, beta_bounds=(0.005, 0.08),
        min_beta_step=2.5e-5, step_recovery_successes=4,
        neutral_tol=1.0e-7, R_tol=1.0e-4,
        corrector_R_step=0.5, max_scan_steps=160,
        max_prediction_step=20.0, min_mode_overlap=0.50,
        max_curve_points=200, min_valid_points=8,
        minimum_complete_points=8, minimum_complete_beta=0.04,
        output_dir=String(output_dir), keep_logs=true,
    )
end

function neutral_columns(data::AbstractMatrix)
    return [abs(data[index, 5]) <= abs(data[index, 7]) ? 1 : 2
            for index in axes(data, 1)]
end

function validate_typeII(result, config::CurveConfig)
    data = result.data
    issues = copy(result.validation.issues)
    size(data, 1) >= 20 || push!(issues, "fewer than 20 Type-II points")
    maximum(data[:, 3]) >= 0.04 || push!(issues, "branch did not reach beta=0.04")
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
    return (ok=isempty(issues), issues=issues, columns=columns)
end

function compute_lopez_typeII(
    Tw::Real;
    output_dir::AbstractString=joinpath(WORKSPACE_ROOT, "neutral_curve_batch"),
    N_cheb::Integer=69,
)
    mkpath(output_dir)
    return mktempdir() do work_dir
        config = typeII_config(Tw, work_dir; N_cheb=N_cheb)
        result = compute_neutral_curve(config)
        check = validate_typeII(result, config)
        check.ok || error(
            "Type-II validation failed for Tw=$Tw: " * join(check.issues, "; "),
        )

        destination = joinpath(
            output_dir,
            NeutralCurveRunner.case_tag(config) * "_branch=typeII.dat",
        )
        cp(result.path, destination; force=true)
        @printf(
            "Tw=%.2f Type-II: points=%d beta=[%.8f, %.8f] R=[%.8f, %.8f] max_residual=%.3e stop=%s\n",
            Tw, size(result.data, 1), minimum(result.data[:, 3]),
            maximum(result.data[:, 3]), minimum(result.data[:, 2]),
            maximum(result.data[:, 2]),
            maximum(min.(abs.(result.data[:, 5]), abs.(result.data[:, 7]))),
            String(result.stop_reason),
        )
        println("saved $destination")
        return merge(result, (path=destination, typeII_validation=check))
    end
end

function main(args=ARGS)
    temperatures = isempty(args) ? (1.18, 1.20) : parse.(Float64, args)
    return [compute_lopez_typeII(Tw) for Tw in temperatures]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
