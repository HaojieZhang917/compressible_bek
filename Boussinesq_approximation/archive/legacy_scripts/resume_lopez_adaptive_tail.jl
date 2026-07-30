const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(WORKSPACE_ROOT, "LopezStability.jl"))
include(joinpath(WORKSPACE_ROOT, "LopezBaseflow.jl"))
include(joinpath(WORKSPACE_ROOT, "CRD_STA.jl"))

using .LopezStability
using .LopezBaseflow
using .CRD_BF
using BSplineKit
using DelimitedFiles
using Dates
using Printf

const N_CHEB = 99
const NEUTRAL_TOL = 1.0e-7
const OUTPUT_DIR = joinpath(WORKSPACE_ROOT, "lopez_neutral_curves")

Tw = isempty(ARGS) ? 1.10 : parse(Float64, ARGS[1])
resume_beta = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0483
tag = replace(@sprintf("%.2f", Tw), "." => "p")
source_path = joinpath(OUTPUT_DIR, "neutral_Tw$(tag).dat")
output_path = joinpath(OUTPUT_DIR, "tail_Tw$(tag).dat")
status_path = joinpath(OUTPUT_DIR, "status_tail_Tw$(tag).log")

raw = readdlm(source_path, Float64; skipstart=2)
row_index = argmin(abs.(raw[:, 3] .- resume_beta))
abs(raw[row_index, 3] - resume_beta) <= 1.0e-8 || error(
    "resume beta $resume_beta is not present in $source_path",
)
row_index >= 2 || error("the resume point needs a preceding point")
previous_row = raw[row_index - 1, :]
resume_row = raw[row_index, :]

profile = solve_lopez_baseflow(Tw)
D, D2, z_cheb = CRD_BF.Cheb(N_CHEB)
sample(values) = begin
    spline = BSplineKit.interpolate(profile.z, values, BSplineOrder(4))
    vec([spline(point) for point in z_cheb])
end
F, G, H, T = sample(profile.F), sample(profile.G), sample(profile.H), sample(profile.T)

open(status_path, "w") do stream
    println(stream, "Tail continuation started at $(Dates.now()) from beta=$resume_beta")
end
log_status(message) = open(status_path, "a") do stream
    println(stream, "$(Dates.now())  $message")
    flush(stream)
end

all_points = NamedTuple[]
last_written = Ref{Any}(nothing)
open(output_path, "w") do stream
    println(
        stream,
        "Variables=\"omega\" \"R\" \"beta\" \"alpha_r\" " *
        "\"alpha_i\" \"residual\" \"overlap\" \"beta_step\"",
    )
    println(stream, "Zone T=\"Lopez Tw=$Tw adaptive tail\"")
    flush(stream)

    function write_point(point, index, beta_step)
        jump = last_written[] === nothing ? 0.0 : abs(point.alpha - last_written[].alpha)
        (
            abs(imag(point.alpha)) <= 2NEUTRAL_TOL &&
            point.residual <= 1.0e-8 &&
            (!isfinite(point.overlap) || point.overlap >= 0.60) &&
            jump <= 0.05
        ) || error("tail live check failed at beta=$(point.beta)")
        @printf(
            stream,
            "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            0.0, point.R, point.beta, real(point.alpha), imag(point.alpha),
            point.residual, point.overlap, beta_step,
        )
        flush(stream)
        log_status(
            @sprintf(
                "POINT i=%d beta=%.7f R=%.8f alpha=(%.9f,%+.3e) res=%.3e overlap=%.8f jump=%.3e step=%.1e OK",
                index, point.beta, point.R, real(point.alpha), imag(point.alpha),
                point.residual, point.overlap, jump, beta_step,
            ),
        )
        last_written[] = point
    end

    fine_step = -0.0001
    ratio = fine_step / (resume_row[3] - previous_row[3])
    initial_R_delta = ratio * (resume_row[2] - previous_row[2])
    initial_alpha_delta = ratio * (
        complex(resume_row[4], resume_row[5]) -
        complex(previous_row[4], previous_row[5])
    )
    fine_start = resume_beta + fine_step
    fine = track_neutral_curve(
        F, G, H, T, 0.0, N_CHEB, D, D2;
        beta_start=fine_start,
        beta_end=0.0441,
        beta_step=fine_step,
        R_start=resume_row[2] + initial_R_delta,
        alpha_start=complex(resume_row[4], resume_row[5]) + initial_alpha_delta,
        initial_R_delta=initial_R_delta,
        initial_alpha_delta=initial_alpha_delta,
        R_step=0.5,
        neutral_tol=NEUTRAL_TOL,
        num_candidates=2,
        on_point=(point, index) -> write_point(point, index, fine_step),
        verbose=true,
    )
    append!(all_points, fine.points)
    fine.failure === nothing || error("fine tail failed: $(fine.failure)")

    previous = all_points[end - 1]
    current = all_points[end]
    low_step = -0.0004
    low_ratio = low_step / (current.beta - previous.beta)
    low_R_delta = low_ratio * (current.R - previous.R)
    low_alpha_delta = low_ratio * (current.alpha - previous.alpha)
    low = track_neutral_curve(
        F, G, H, T, 0.0, N_CHEB, D, D2;
        beta_start=0.0437,
        beta_end=0.0301,
        beta_step=low_step,
        R_start=current.R + low_R_delta,
        alpha_start=current.alpha + low_alpha_delta,
        vector_start=current.vector,
        initial_R_delta=low_R_delta,
        initial_alpha_delta=low_alpha_delta,
        R_step=1.0,
        neutral_tol=NEUTRAL_TOL,
        num_candidates=2,
        on_point=(point, index) -> write_point(
            point, length(all_points) + index, low_step,
        ),
        verbose=true,
    )
    append!(all_points, low.points)
    low.failure === nothing || error("low tail failed: $(low.failure)")
end

log_status("COMPLETE points=$(length(all_points)) file=$output_path")
