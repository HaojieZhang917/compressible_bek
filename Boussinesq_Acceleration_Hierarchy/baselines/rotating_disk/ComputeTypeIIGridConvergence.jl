include(joinpath(@__DIR__, "NeutralCurveRunner.jl"))

using .NeutralCurveRunner

function parse_env(name, default, parser)
    return parser(get(ENV, name, string(default)))
end

function main()
    Tw = parse_env("TYPEII_GRID_TW", 1.12, x -> parse(Float64, x))
    N_cheb = parse_env("TYPEII_GRID_N_CHEB", 69, x -> parse(Int, x))
    beta_step = parse_env(
        "TYPEII_GRID_BETA_STEP", 2.0e-4, x -> parse(Float64, x),
    )
    output_root = get(
        ENV,
        "TYPEII_GRID_OUTPUT_ROOT",
        joinpath(@__DIR__, "typeII_grid_convergence_results", "raw"),
    )
    output_dir = joinpath(output_root, "N=$(N_cheb)")
    # Keep the initial shift centered on the converged Type-II eigenvalue.
    # A fixed target of 0.11 can select a nearby secondary Ritz value at high
    # resolution for Tw=1.12 even though the desired branch remains present.
    alpha_target = 0.1033 + 0.13 * (Tw - 1.12)

    config = NeutralCurveRunner.CurveConfig(
        Tw=Tw,
        omega=0.0,
        R_initial=440.0,
        beta_initial=0.038,
        alpha_target=alpha_target,
        num_modes=2,
        model=:compressible,
        Mr=0.3,
        Ro=-1.0,
        N_cheb=N_cheb,
        property_perturbations=true,
        base_property_variation=true,
        beta_step=beta_step,
        beta_scan_step=1.0e-4,
        beta_bounds=(0.036, 0.0475),
        min_beta_step=2.5e-5,
        step_recovery_successes=8,
        step_growth_factor=2.0,
        neutral_tol=1.0e-8,
        beta_tol=1.0e-10,
        R_tol=1.0e-5,
        corrector_R_step=0.25,
        max_scan_steps=160,
        max_prediction_step=8.0,
        min_mode_overlap=0.70,
        max_curve_points=100,
        min_valid_points=20,
        minimum_complete_points=40,
        minimum_complete_beta=0.0465,
        output_dir=output_dir,
        keep_logs=true,
    )

    println(
        "Computing Type-II grid check: Tw=$(Tw), N=$(N_cheb), " *
        "beta_step=$(beta_step)",
    )
    flush(stdout)
    result = NeutralCurveRunner.compute_neutral_curve(config)
    println(
        "Completed Tw=$(Tw), N=$(N_cheb): points=$(size(result.data, 1)), " *
        "stop=$(result.stop_reason), valid=$(result.validation.ok), " *
        "path=$(result.path)",
    )
end

main()
