const WORKSPACE_ROOT = @__DIR__
include(joinpath(WORKSPACE_ROOT, "ComputeBlackburnNeutralCurves.jl"))

using .BlackburnNeutralCurveRunner
using Printf

const DEFAULT_EXTENDED_TEMPERATURES = (1.40, 1.60, 1.80)

function extended_config(model::Symbol, Tw::Real; N_cheb::Integer=69)
    model in (:blackburn, :compressible) || error(
        "ZERO_FREQUENCY_MODEL must be blackburn or compressible",
    )
    output_dir = model === :blackburn ?
        joinpath(WORKSPACE_ROOT, "blackburn_neutral_curve_batch") :
        joinpath(WORKSPACE_ROOT, "neutral_curve_batch")
    config = production_config(Tw, output_dir; N_cheb=N_cheb)
    if model === :compressible
        config = BlackburnNeutralCurveRunner.with_config(
            config;
            model=:compressible,
            property_perturbations=true,
            base_property_variation=true,
        )
    end
    # At Tw >= 1.4 the low-beta R=500 seed no longer brackets a neutral
    # crossing.  Start from the high-beta Type-I endpoint and continue toward
    # lower beta; this direction still traverses every connected fold.
    if Tw >= 1.4
        config = BlackburnNeutralCurveRunner.typeI_fallback_config(config)
        return BlackburnNeutralCurveRunner.with_config(
            config;
            beta_scan_step=1.0e-3,
            beta_bounds=(1.0e-3, 0.14),
            corrector_R_step=1.0,
            min_mode_overlap=0.60,
            minimum_complete_beta=0.055,
        )
    end
    return config
end

function completed_curve(config::CurveConfig)
    path = curve_path(config)
    isfile(path) || return false
    validation = validate_curve_file(path, config)
    validation.ok || return false
    issues = BlackburnNeutralCurveRunner.file_completion_issues(
        validation.data, config,
    )
    isempty(issues) || return false
    if size(validation.data, 1) >= 2 && validation.data[1, 3] > validation.data[end, 3]
        BlackburnNeutralCurveRunner.write_curve(
            path, config, reverse(validation.data; dims=1),
        )
    end
    return true
end

function compute_case(model::Symbol, Tw::Real; N_cheb::Integer=69)
    config = extended_config(model, Tw; N_cheb=N_cheb)
    if completed_curve(config)
        @printf("Skipping validated case model=%s Tw=%.2f path=%s\n",
            String(model), Tw, curve_path(config))
        return (status=:skipped, path=curve_path(config))
    end
    @printf("Computing zero-frequency curve model=%s Tw=%.2f N_cheb=%d\n",
        String(model), Tw, N_cheb)
    flush(stdout)
    result = run_case_with_retries(config)
    residual = maximum(min.(abs.(result.data[:, 5]), abs.(result.data[:, 7])))
    @printf(
        "Completed model=%s Tw=%.2f points=%d beta=[%.9f, %.9f] R=[%.9f, %.9f] max_residual=%.3e stop=%s\n",
        String(model), Tw, size(result.data, 1),
        minimum(result.data[:, 3]), maximum(result.data[:, 3]),
        minimum(result.data[:, 2]), maximum(result.data[:, 2]),
        residual, String(result.stop_reason),
    )
    println("Saved ", result.path)
    return result
end

function main(args=ARGS)
    model = Symbol(lowercase(get(ENV, "ZERO_FREQUENCY_MODEL", "blackburn")))
    N_cheb = parse(Int, get(ENV, "ZERO_FREQUENCY_N_CHEB", "69"))
    temperatures = isempty(args) ? collect(DEFAULT_EXTENDED_TEMPERATURES) :
        parse.(Float64, args)
    return [compute_case(model, Tw; N_cheb=N_cheb) for Tw in temperatures]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
