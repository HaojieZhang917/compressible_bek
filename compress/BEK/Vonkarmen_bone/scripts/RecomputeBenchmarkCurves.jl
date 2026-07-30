const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(WORKSPACE_ROOT,"NeutralCurveRunner.jl"))

using .NeutralCurveRunner

function benchmark_config(model::Symbol)
    return CurveConfig(
        Tw=1.0,
        omega=0.0,
        R_initial=500.0,
        beta_initial=0.04,
        alpha_target=0.1,
        num_modes=1,
        model=model,
        Mr=0.3,
        Ro=-1.0,
        N_cheb=69,
        property_perturbations=true,
        base_property_variation=true,
        beta_step=8.0e-4,
        neutral_tol=1.0e-7,
        output_dir=joinpath(WORKSPACE_ROOT,"neutral_curve_batch"),
        keep_logs=false,
    )
end

function main()
    for model in (:lopez,:compressible)
        config = benchmark_config(model)
        println("Computing Tw=1.0 model=$model")
        result = NeutralCurveRunner.run_case_with_retries(config)
        validation = NeutralCurveRunner.validate_curve_file(result.path,config)
        validation.ok || error(
            "Validation failed for model=$model: $(join(validation.issues, "; "))",
        )
        println("Validated $(result.path) with $(size(result.data,1)) points")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
