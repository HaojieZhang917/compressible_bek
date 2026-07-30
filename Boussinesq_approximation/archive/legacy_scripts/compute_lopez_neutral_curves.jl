const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(WORKSPACE_ROOT, "LopezStability.jl"))
include(joinpath(WORKSPACE_ROOT, "LopezBaseflow.jl"))
include(joinpath(WORKSPACE_ROOT, "CRD_STA.jl"))

using .LopezStability
using .LopezBaseflow
using .CRD_BF
using BSplineKit
using Dates
using Printf

const OMEGA = 0.0
const N_CHEB = 99
const BETA_START = 0.1185
const BETA_STEP = -0.0008
const NEUTRAL_TOL = 1.0e-7
const RESIDUAL_LIMIT = 1.0e-8
const OVERLAP_LIMIT = 0.60
const ALPHA_JUMP_LIMIT = 0.05

const OUTPUT_DIR = joinpath(WORKSPACE_ROOT, "lopez_neutral_curves")
const STATUS_PATH = joinpath(OUTPUT_DIR, "status.log")

function wall_temperature_tag(Tw)
    return replace(@sprintf("%.2f", Tw), "." => "p")
end

function append_status(message)
    open(STATUS_PATH, "a") do stream
        println(stream, "$(Dates.now())  $message")
        flush(stream)
    end
end

function sample_on_cheb(z, values, points)
    spline = BSplineKit.interpolate(z, values, BSplineOrder(4))
    return vec([spline(point) for point in points])
end

function calculate_curve(Tw, beta_end)
    tag = wall_temperature_tag(Tw)
    output_path = joinpath(OUTPUT_DIR, "neutral_Tw$(tag).dat")
    append_status("START Tw=$Tw beta_end=$beta_end")

    profile = solve_lopez_baseflow(Tw)
    D, D2, z_cheb = CRD_BF.Cheb(N_CHEB)
    F = sample_on_cheb(profile.z, profile.F, z_cheb)
    G = sample_on_cheb(profile.z, profile.G, z_cheb)
    H = sample_on_cheb(profile.z, profile.H, z_cheb)
    T = sample_on_cheb(profile.z, profile.T, z_cheb)

    previous_point = Ref{Any}(nothing)
    result = open(output_path, "w") do stream
        println(
            stream,
            "Variables=\"omega\" \"R\" \"beta\" \"alpha_r\" " *
            "\"alpha_i\" \"residual\" \"overlap\"",
        )
        println(stream, "Zone T=\"Lopez Tw=$Tw\"")
        flush(stream)

        on_point = function (point, index)
            alpha_jump = previous_point[] === nothing ? 0.0 :
                         abs(point.alpha - previous_point[].alpha)
            overlap_ok = !isfinite(point.overlap) || point.overlap >= OVERLAP_LIMIT
            checks_ok = (
                abs(imag(point.alpha)) <= 2NEUTRAL_TOL &&
                point.residual <= RESIDUAL_LIMIT &&
                overlap_ok &&
                alpha_jump <= ALPHA_JUMP_LIMIT
            )
            checks_ok || error(
                "live validation failed at Tw=$Tw, beta=$(point.beta): " *
                "alpha_i=$(imag(point.alpha)), residual=$(point.residual), " *
                "overlap=$(point.overlap), alpha_jump=$alpha_jump",
            )
            @printf(
                stream,
                "%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                OMEGA, point.R, point.beta, real(point.alpha), imag(point.alpha),
                point.residual, point.overlap,
            )
            flush(stream)
            append_status(
                @sprintf(
                    "POINT Tw=%.2f i=%d beta=%.7f R=%.8f alpha=(%.9f,%+.3e) res=%.3e overlap=%.8f jump=%.3e OK",
                    Tw, index, point.beta, point.R, real(point.alpha),
                    imag(point.alpha), point.residual, point.overlap, alpha_jump,
                ),
            )
            previous_point[] = point
        end

        track_neutral_curve(
            F, G, H, T, OMEGA, N_CHEB, D, D2;
            beta_start=BETA_START,
            beta_end=beta_end,
            beta_step=BETA_STEP,
            R_start=495.0,
            alpha_start=0.6716502816766712 + 0im,
            initial_R_delta=-10.8,
            initial_alpha_delta=-0.00350 + 0im,
            R_step=1.0,
            neutral_tol=NEUTRAL_TOL,
            num_candidates=2,
            min_overlap=OVERLAP_LIMIT,
            max_alpha_jump=ALPHA_JUMP_LIMIT,
            residual_tol=RESIDUAL_LIMIT,
            on_point=on_point,
            verbose=true,
        )
    end

    if result.failure === nothing
        append_status("COMPLETE Tw=$Tw points=$(length(result.points)) file=$output_path")
    else
        append_status(
            "FAILED Tw=$Tw points=$(length(result.points)) beta=$(result.failure.beta) " *
            "message=$(result.failure.message)",
        )
    end
    return result
end

mkpath(OUTPUT_DIR)
open(STATUS_PATH, "w") do stream
    println(stream, "Lopez neutral-curve run started at $(Dates.now())")
end

beta_end = isempty(ARGS) ? 0.0401 : parse(Float64, ARGS[1])
temperatures = length(ARGS) <= 1 ? [1.05, 1.10, 1.20] : parse.(Float64, ARGS[2:end])

results = Dict{Float64, Any}()
for Tw in temperatures
    try
        results[Tw] = calculate_curve(Tw, beta_end)
    catch exception
        append_status("ABORTED Tw=$Tw message=$(sprint(showerror, exception))")
        showerror(stderr, exception, catch_backtrace())
        println(stderr)
    end
end

all_complete = all(
    haskey(results, Tw) && results[Tw].failure === nothing for Tw in temperatures
)
all_complete || exit(1)
