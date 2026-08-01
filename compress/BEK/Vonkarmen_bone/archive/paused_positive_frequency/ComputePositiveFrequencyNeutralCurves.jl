const WORKSPACE_ROOT = @__DIR__
include(joinpath(WORKSPACE_ROOT, "BlackburnNeutralCurveRunner.jl"))

using .BlackburnNeutralCurveRunner
using Printf

const DEFAULT_TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20)

function positive_frequency_config(
    Tw::Real,
    model::Symbol,
    omega::Real,
    output_dir::AbstractString;
    N_cheb::Integer=69,
    R_initial::Real=400.0,
)
    model in (:blackburn, :compressible) ||
        throw(ArgumentError("model must be :blackburn or :compressible"))
    return CurveConfig(
        Tw=Float64(Tw),
        omega=Float64(omega),
        R_initial=Float64(R_initial),
        beta_initial=0.035,
        # This shift follows the same Type-II spatial branch continuously:
        # alpha_r≈0.152, 0.159, 0.167 at Tw=1.00, 1.04, 1.08.
        alpha_target=0.16,
        num_modes=2,
        model=model,
        Mr=0.3,
        Ro=-1.0,
        N_cheb=Int(N_cheb),
        property_perturbations=true,
        base_property_variation=true,
        beta_step=1.2e-3,
        beta_scan_step=6.0e-4,
        # Lingwood's omega=0.008 Type-II lobe closes at beta≈0.040.
        # Stop immediately before the Type-I/Type-II neutral branch point so
        # that mode switching there cannot contaminate the Type-II comparison.
        beta_bounds=(-0.02, 0.0398),
        min_beta_step=2.5e-5,
        step_recovery_successes=8,
        step_growth_factor=2.0,
        neutral_tol=1.0e-8,
        beta_tol=1.0e-10,
        R_tol=1.0e-5,
        corrector_R_step=0.5,
        max_scan_steps=160,
        max_prediction_step=15.0,
        min_mode_overlap=0.65,
        allow_mode_switch=false,
        max_curve_points=180,
        min_valid_points=20,
        minimum_complete_points=20,
        minimum_complete_beta=0.0385,
        output_dir=String(output_dir),
        keep_logs=true,
    )
end

function compute_case(
    Tw::Real,
    model::Symbol,
    omega::Real;
    output_dir::AbstractString,
    N_cheb::Integer=69,
)
    mkpath(output_dir)
    # Heating moves the whole positive-frequency Type-II lobe downward in R.
    # Start from a temperature-dependent section and try nearby sections only
    # when that section does not intersect the lobe.
    nominal_R = 400.0 - 100.0 * (Float64(Tw) - 1.0)
    R_candidates = unique([
        nominal_R,
        nominal_R - 20.0,
        nominal_R + 20.0,
        nominal_R - 40.0,
        nominal_R + 40.0,
    ])
    last_exception = nothing
    for R_initial in R_candidates
        config = positive_frequency_config(
            Tw, model, omega, output_dir;
            N_cheb=N_cheb, R_initial=R_initial,
        )
        @printf(
            "Computing positive-frequency curve: model=%s omega=%.6f Tw=%.2f N_cheb=%d R_initial=%.1f\n",
            String(model), omega, Tw, N_cheb, R_initial,
        )
        flush(stdout)
        result = try
            # This calculation deliberately isolates Type II and stops before
            # its branch point.  Do not invoke the generic Type-I fallback.
            compute_neutral_curve(config)
        catch exception
            exception isa InterruptException && rethrow()
            message = sprint(showerror, exception)
            if occursin("No initial neutral crossing", message)
                last_exception = exception
                @printf(
                    "No Type-II crossing at R_initial=%.1f; trying the next section\n",
                    R_initial,
                )
                flush(stdout)
                continue
            end
            rethrow()
        end
        residual = maximum(min.(
            abs.(result.data[:, 5]), abs.(result.data[:, 7]),
        ))
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
    last_exception === nothing && error("No R_initial candidates were generated")
    throw(last_exception)
end

function main(args=ARGS)
    temperatures = isempty(args) ? collect(DEFAULT_TEMPERATURES) :
        parse.(Float64, args)
    omega = parse(Float64, get(ENV, "POSITIVE_FREQUENCY_OMEGA", "0.008"))
    model = Symbol(lowercase(get(ENV, "POSITIVE_FREQUENCY_MODEL", "blackburn")))
    N_cheb = parse(Int, get(ENV, "POSITIVE_FREQUENCY_N_CHEB", "69"))
    output_dir = get(
        ENV,
        "POSITIVE_FREQUENCY_OUTPUT_DIR",
        joinpath(
            WORKSPACE_ROOT,
            "positive_frequency_neutral_curve_batch",
            @sprintf("omega=%.3f", omega),
        ),
    )
    return [
        compute_case(
            Tw, model, omega;
            output_dir=output_dir,
            N_cheb=N_cheb,
        )
        for Tw in temperatures
    ]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
