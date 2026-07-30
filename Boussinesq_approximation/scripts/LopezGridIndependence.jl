const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(WORKSPACE_ROOT, "NeutralCurveRunner.jl"))

using .NeutralCurveRunner
using BSplineKit
using DelimitedFiles
using LinearAlgebra
using Printf

const TW = 1.08
const GRID_SIZES = (49, 59, 69, 79, 89)
const COMMON_Z = collect(range(0.0, 20.0; length=801))

const MODE_REFERENCES = (
    (
        name="Type-I", R=263.707514331873,
        beta=0.069884194636872, alpha=0.390278130496983,
        beta_half_width=4.0e-4,
    ),
    (
        name="Type-II", R=303.521565914879,
        beta=0.0450380158825237, alpha=0.145743738926843,
        beta_half_width=4.0e-4,
    ),
)

function grid_config(N::Integer)
    return CurveConfig(
        Tw=TW, omega=0.0, R_initial=300.0,
        beta_initial=0.06, alpha_target=0.3,
        num_modes=2, model=:lopez, N_cheb=Int(N),
        neutral_tol=1.0e-8, R_tol=1.0e-6,
    )
end

function initial_mode_state(solve_at, reference)
    values, vectors = solve_at(
        reference.R, reference.beta, reference.alpha, 2,
    )
    active_index = argmin(abs.(values .- reference.alpha))
    return (
        R=reference.R, beta=reference.beta,
        values=values, vectors=vectors, active_index=active_index,
    )
end

function neutral_root(solve_at, beta, R_guess, seed)
    return NeutralCurveRunner.NeutralContinuation.find_neutral_R(
        solve_at, beta, R_guess, seed.values, seed.vectors,
        seed.active_index;
        R_step=0.2, preferred_direction=0,
        max_scan_steps=100, max_refine=50,
        neutral_tol=1.0e-8, R_tol=1.0e-6,
        R_bounds=(100.0, 600.0), max_R_deviation=15.0,
        min_overlap=0.75,
    )
end

function parabola_vertex(betas, reynolds)
    center = betas[2]
    scale = maximum(abs.(betas .- center))
    x = (betas .- center) ./ scale
    coefficients = hcat(ones(3), x, x .^ 2) \ reynolds
    coefficients[3] > 0 || error(
        "Local neutral curve is not convex: coefficients=$coefficients",
    )
    x_vertex = clamp(-coefficients[2] / (2coefficients[3]), -1.25, 1.25)
    return center + scale * x_vertex
end

function local_critical_point(solve_at, reference)
    initial = initial_mode_state(solve_at, reference)
    center = neutral_root(
        solve_at, reference.beta, reference.R, initial,
    )
    half_width = reference.beta_half_width

    for _ in 1:2
        left = neutral_root(
            solve_at, center.beta - half_width, center.R, center,
        )
        right = neutral_root(
            solve_at, center.beta + half_width, center.R, center,
        )
        beta_vertex = parabola_vertex(
            [left.beta, center.beta, right.beta],
            [left.R, center.R, right.R],
        )
        nearest = abs(beta_vertex - left.beta) < abs(beta_vertex - right.beta) ?
            left : right
        abs(beta_vertex - center.beta) < abs(beta_vertex - nearest.beta) &&
            (nearest = center)
        center = neutral_root(
            solve_at, beta_vertex, minimum((left.R, center.R, right.R)), nearest,
        )
        half_width /= 2
    end
    return center
end

function full_mode(config::CurveConfig, prepared, critical)
    z, H0, F0, G0, T0, _, _, _, _ =
        NeutralCurveRunner.LopezBaseflow.get_baseflow(config.Tw)
    F = NeutralCurveRunner.sample_lopez_profile(z, F0, prepared.x)
    G = NeutralCurveRunner.sample_lopez_profile(z, G0, prepared.x)
    H = NeutralCurveRunner.sample_lopez_profile(z, H0, prepared.x)
    T = NeutralCurveRunner.sample_lopez_profile(z, T0, prepared.x)
    active = critical.active_index
    values, vectors, info = NeutralCurveRunner.LopezStability.eigsol_lopez(
        F, G, H, T, critical.R, config.omega, critical.beta,
        config.N_cheb, prepared.D, prepared.D2,
        real(critical.values[active]), 1;
        initial_vector=critical.vectors[:, active],
        full_eigenvectors=true, return_info=true,
    )
    return (
        alpha=values[1], vector=vectors[:, 1],
        polynomial_residual=info.residuals[1], x=vec(prepared.x),
    )
end

function unique_coordinates(x, values)
    keep = Int[]
    for index in eachindex(x)
        (isempty(keep) || x[index] > x[last(keep)] + 1.0e-12) && push!(keep, index)
    end
    return x[keep], values[keep]
end

function spectral_sample(x, values, points)
    x_unique, value_unique = unique_coordinates(x, values)
    real_interpolation = BSplineKit.interpolate(
        x_unique, real.(value_unique), BSplineOrder(4),
    )
    imaginary_interpolation = BSplineKit.interpolate(
        x_unique, imag.(value_unique), BSplineOrder(4),
    )
    return ComplexF64.(
        real_interpolation.(points) .+ im .* imaginary_interpolation.(points),
    )
end

function sampled_physical_mode(result)
    n = result.N + 1
    fields = Matrix{ComplexF64}(undef, length(COMMON_Z), 4)
    for component in 1:4
        indices = (component - 1) * n + 1:component * n
        fields[:, component] .= spectral_sample(
            result.x, result.vector[indices], COMMON_Z,
        )
    end
    scale = norm(fields[:, 1:3])
    scale > eps(Float64) || error("Velocity eigenfunction norm is zero")
    return fields ./ scale
end

function aligned_shape_error(candidate, reference)
    phase = dot(vec(reference), vec(candidate))
    aligned = abs(phase) > eps(Float64) ?
        candidate .* (conj(phase) / abs(phase)) : candidate
    velocity_error = norm(aligned[:, 1:3] - reference[:, 1:3]) /
        norm(reference[:, 1:3])
    combined_error = norm(aligned - reference) / norm(reference)
    return velocity_error, combined_error
end

function compute_case(N::Integer, reference)
    config = grid_config(N)
    prepared = NeutralCurveRunner.prepare_solver(config)
    critical = local_critical_point(prepared.solve_at, reference)
    mode = full_mode(config, prepared, critical)
    active = critical.active_index
    @printf(
        "%s N=%d: R=%.10f beta=%.10f alpha=(%.10f,%+.3e) pep=%.3e\n",
        reference.name, N, critical.R, critical.beta,
        real(mode.alpha), imag(mode.alpha), mode.polynomial_residual,
    )
    flush(stdout)
    return (
        mode=reference.name, N=Int(N), R=critical.R,
        beta=critical.beta, alpha=mode.alpha,
        neutral_residual=abs(imag(mode.alpha)),
        polynomial_residual=mode.polynomial_residual,
        x=mode.x, vector=mode.vector,
    )
end

function write_results(results, output_dir)
    mkpath(output_dir)
    rows = NamedTuple[]
    for mode_name in ("Type-I", "Type-II")
        mode_results = filter(result -> result.mode == mode_name, results)
        sort!(mode_results; by=result -> result.N)
        reference = sampled_physical_mode(mode_results[end])
        R_reference = mode_results[end].R
        beta_reference = mode_results[end].beta
        alpha_reference = real(mode_results[end].alpha)
        for result in mode_results
            sampled = sampled_physical_mode(result)
            velocity_error, combined_error = aligned_shape_error(sampled, reference)
            push!(rows, (
                mode=result.mode, N=result.N, R=result.R, beta=result.beta,
                alpha_r=real(result.alpha), alpha_i=imag(result.alpha),
                neutral_residual=result.neutral_residual,
                polynomial_residual=result.polynomial_residual,
                R_relative_error=abs(result.R - R_reference) / abs(R_reference),
                beta_relative_error=abs(result.beta - beta_reference) /
                    abs(beta_reference),
                alpha_relative_error=abs(real(result.alpha) - alpha_reference) /
                    abs(alpha_reference),
                velocity_shape_error=velocity_error,
                combined_shape_error=combined_error,
            ))
        end
    end

    path = joinpath(output_dir, "lopez_grid_independence_Tw1.08.tsv")
    open(path, "w") do io
        names = propertynames(first(rows))
        println(io, join(names, '\t'))
        for row in rows
            println(io, join((getproperty(row, name) for name in names), '\t'))
        end
    end
    println("saved $path")
    return (rows=rows, path=path)
end

function main()
    results = NamedTuple[]
    for reference in MODE_REFERENCES
        for N in GRID_SIZES
            push!(results, compute_case(N, reference))
        end
    end
    return write_results(results, joinpath(WORKSPACE_ROOT, "grid_independence"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
