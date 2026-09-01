const WORKSPACE_ROOT = @__DIR__
include(joinpath(WORKSPACE_ROOT, "BlackburnNeutralCurveRunner.jl"))

using .BlackburnNeutralCurveRunner
using Printf

const DEFAULT_TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20)

function production_config(
    Tw::Real, output_dir::AbstractString;
    N_cheb::Integer=69,
)
    return CurveConfig(
        Tw=Float64(Tw),
        omega=0.0,
        R_initial=500.0,
        beta_initial=0.04,
        alpha_target=0.10,
        num_modes=2,
        model=:blackburn,
        Mr=0.3,
        Ro=-1.0,
        N_cheb=Int(N_cheb),
        beta_step=8.0e-4,
        beta_scan_step=5.0e-4,
        beta_bounds=(1.0e-3, 0.20),
        min_beta_step=5.0e-5,
        neutral_tol=1.0e-7,
        R_tol=1.0e-4,
        corrector_R_step=1.0,
        max_scan_steps=80,
        max_prediction_step=30.0,
        min_mode_overlap=0.60,
        max_curve_points=500,
        min_valid_points=8,
        minimum_complete_points=50,
        minimum_complete_beta=0.08,
        output_dir=String(output_dir),
        keep_logs=true,
    )
end

function compute_temperature(
    Tw::Real;
    output_dir::AbstractString=joinpath(
        WORKSPACE_ROOT, "blackburn_neutral_curve_batch",
    ),
    N_cheb::Integer=69,
)
    mkpath(output_dir)
    config = production_config(Tw, output_dir; N_cheb=N_cheb)
    @printf(
        "Computing Blackburn neutral curve: Tw=%.2f N_cheb=%d output=%s\n",
        Tw, N_cheb, output_dir,
    )
    result = run_case_with_retries(config)
    residual = maximum(min.(abs.(result.data[:, 5]), abs.(result.data[:, 7])))
    @printf(
        "Completed Tw=%.2f points=%d beta=[%.9f, %.9f] R=[%.9f, %.9f] max_residual=%.3e stop=%s\n",
        Tw, size(result.data, 1),
        minimum(result.data[:, 3]), maximum(result.data[:, 3]),
        minimum(result.data[:, 2]), maximum(result.data[:, 2]),
        residual, String(result.stop_reason),
    )
    println("Saved ", result.path)
    return result
end

function main(args=ARGS)
    temperatures = isempty(args) ? collect(DEFAULT_TEMPERATURES) :
        parse.(Float64, args)
    output_dir = get(
        ENV,
        "BLACKBURN_OUTPUT_DIR",
        joinpath(WORKSPACE_ROOT, "blackburn_neutral_curve_batch"),
    )
    N_cheb = parse(Int, get(ENV, "BLACKBURN_N_CHEB", "69"))
    return [
        compute_temperature(
            Tw;
            output_dir=output_dir,
            N_cheb=N_cheb,
        )
        for Tw in temperatures
    ]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
