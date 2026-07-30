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
const NEUTRAL_TOL = 1.0e-7
const OUTPUT_DIR = joinpath(WORKSPACE_ROOT, "lopez_neutral_curves")

Tw = isempty(ARGS) ? 1.10 : parse(Float64, ARGS[1])
tag = replace(@sprintf("%.2f", Tw), "." => "p")
output_path = joinpath(OUTPUT_DIR, "neutral_Tw$(tag).dat")
status_path = joinpath(OUTPUT_DIR, "status_Tw$(tag).log")
mkpath(OUTPUT_DIR)

function log_status(message)
    open(status_path, "a") do stream
        println(stream, "$(Dates.now())  $message")
        flush(stream)
    end
end

profile = solve_lopez_baseflow(Tw)
D, D2, z_cheb = CRD_BF.Cheb(N_CHEB)
sample(values) = begin
    spline = BSplineKit.interpolate(profile.z, values, BSplineOrder(4))
    vec([spline(point) for point in z_cheb])
end
F, G, H, T = sample(profile.F), sample(profile.G), sample(profile.H), sample(profile.T)

open(status_path, "w") do stream
    println(stream, "Adaptive Lopez neutral curve started at $(Dates.now())")
end

all_points = NamedTuple[]
last_written = Ref{Any}(nothing)

result = open(output_path, "w") do stream
    println(
        stream,
        "Variables=\"omega\" \"R\" \"beta\" \"alpha_r\" " *
        "\"alpha_i\" \"residual\" \"overlap\" \"beta_step\"",
    )
    println(stream, "Zone T=\"Lopez Tw=$Tw adaptive\"")
    flush(stream)

    function write_point(point, index, beta_step)
        alpha_jump = last_written[] === nothing ? 0.0 :
                     abs(point.alpha - last_written[].alpha)
        checks_ok = (
            abs(imag(point.alpha)) <= 2NEUTRAL_TOL &&
            point.residual <= 1.0e-8 &&
            (!isfinite(point.overlap) || point.overlap >= 0.60) &&
            alpha_jump <= 0.05
        )
        checks_ok || error(
            "live validation failed at beta=$(point.beta): " *
            "alpha_i=$(imag(point.alpha)), residual=$(point.residual), " *
            "overlap=$(point.overlap), alpha_jump=$alpha_jump",
        )
        @printf(
            stream,
            "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            OMEGA, point.R, point.beta, real(point.alpha), imag(point.alpha),
            point.residual, point.overlap, beta_step,
        )
        flush(stream)
        log_status(
            @sprintf(
                "POINT i=%d beta=%.7f R=%.8f alpha=(%.9f,%+.3e) res=%.3e overlap=%.8f jump=%.3e step=%.1e OK",
                index, point.beta, point.R, real(point.alpha), imag(point.alpha),
                point.residual, point.overlap, alpha_jump, beta_step,
            ),
        )
        last_written[] = point
    end

    segments = [
        (beta_start=0.1185, beta_end=0.0553, beta_step=-0.0008),
        (beta_start=0.0552, beta_end=0.0441, beta_step=-0.0001),
        (beta_start=0.0437, beta_end=0.0301, beta_step=-0.0004),
    ]

    last_result = nothing
    for (segment_index, segment) in enumerate(segments)
        if isempty(all_points)
            R_start = 495.0
            alpha_start = 0.6716502816766712 + 0im
            vector_start = nothing
            initial_R_delta = -10.8
            initial_alpha_delta = -0.00350 + 0im
        else
            previous = all_points[end - 1]
            current = all_points[end]
            ratio = segment.beta_step / (current.beta - previous.beta)
            initial_R_delta = ratio * (current.R - previous.R)
            initial_alpha_delta = ratio * (current.alpha - previous.alpha)
            R_start = current.R + initial_R_delta
            alpha_start = current.alpha + initial_alpha_delta
            vector_start = current.vector
        end
        offset = length(all_points)
        log_status(
            "SEGMENT $segment_index start=$(segment.beta_start) " *
            "end=$(segment.beta_end) step=$(segment.beta_step)",
        )
        last_result = track_neutral_curve(
            F, G, H, T, OMEGA, N_CHEB, D, D2;
            beta_start=segment.beta_start,
            beta_end=segment.beta_end,
            beta_step=segment.beta_step,
            R_start=R_start,
            alpha_start=alpha_start,
            vector_start=vector_start,
            initial_R_delta=initial_R_delta,
            initial_alpha_delta=initial_alpha_delta,
            R_step=segment_index == 2 ? 0.5 : 1.0,
            neutral_tol=NEUTRAL_TOL,
            num_candidates=2,
            on_point=(point, index) -> write_point(
                point, offset + index, segment.beta_step,
            ),
            verbose=true,
        )
        append!(all_points, last_result.points)
        last_result.failure === nothing || break
    end
    last_result
end

if result.failure === nothing
    log_status("COMPLETE points=$(length(all_points)) file=$output_path")
else
    log_status(
        "FAILED points=$(length(all_points)) beta=$(result.failure.beta) " *
        "message=$(result.failure.message)",
    )
    error("adaptive neutral curve failed: $(result.failure)")
end
