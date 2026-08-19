using LinearAlgebra
using BSplineKit
using DelimitedFiles
using Printf
using PyCall
using NonlinearEigenproblems

include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")

Base.@kwdef mutable struct Config
    target_radius::Float64 = 570.0
    Re_h::Int = 1000
    mass_flux::Float64 = 0.0
    mode::Int = 1
    Ro::Float64 = -1.0
    azimuthal_mode::Int = 30
    omega_bar::Float64 = 0.0
    N_cheb::Int = 99
    scan_start::Float64 = 200.0
    scan_step::Float64 = 2.0
    integration_step::Float64 = 2.0
    alpha_guess_real::Float64 = 1.05
    alpha_guess_imag::Float64 = 0.0
    candidate_count::Int = 1
    eig_tolerance::Float64 = 1e-13
    neutral_alpha_tolerance::Float64 = 1e-9
    neutral_radius_tolerance::Float64 = 1e-5
    neutral_max_iterations::Int = 30
    roughness_height::Union{Nothing, Float64} = nothing
    # Fourier parameter corresponding to the manuscript Gaussian variance
    # c^2=1: l_s=1/(2c^2)=1/2. This reproduces the published Cr curves.
    roughness_localization::Float64 = 0.5
    output_prefix::String = "total_amplitude"
end

function usage()
    println("""
Compute the total downstream amplitude of one spatial instability branch:

    N(R)       = -integral(alpha_i dR)
    gain(R)    = exp(N(R))
    |A(R)|     = |C_r(R_l)| exp(N(R))

The lower neutral radius R_l is found automatically unless --lower-radius is
supplied. The initial receptivity coefficient is evaluated at R_l using the
same roughness projection and normalization as receptivity.ipynb.

Usage:
  julia total_amplitude.jl [target_radius] [options]

Examples:
  julia total_amplitude.jl 570
  julia total_amplitude.jl --target-radius 570 --n 30 --omega-bar 0
  julia total_amplitude.jl 570 --lower-radius 283.0 --output-prefix baseline

Options:
  --target-radius R       Radius at which the total amplitude is requested.
  --lower-radius R        Use a known lower neutral radius; otherwise find it.
  --Re N                  Cavity Reynolds number Re_h (default: 1000).
  --mass-flux A           Injection/suction coefficient a_s (default: 0).
  --n N                   Integer azimuthal mode number (default: 30).
  --omega-bar W           Scaled frequency R*omega (default: 0).
  --N-cheb N              Chebyshev polynomial order (default: 99).
  --scan-start R          First radius used to find R_l (default: 200).
  --scan-step DR          Neutral-point scan spacing (default: 2).
  --integration-step DR   Radial spacing for the N-factor integral (default: 2).
  --alpha-real X          Real part of initial eigenvalue shift (default: 1.05).
  --alpha-imag X          Imaginary part of initial shift (default: 0).
  --candidates N          Candidate roots used for branch tracking (default: 1).
                          Use 2 or 3 near a modal-switching region; the script
                          then also uses eigenfunction overlap for selection.
  --roughness-height H    Wall-displacement height (default: 1/Re_h).
  --roughness-ls L        Fourier localization parameter (default: 0.5,
                          corresponding to manuscript c^2=1).
  --output-prefix NAME    Output filename prefix (default: total_amplitude).
  --help                  Show this message.

Outputs:
  <prefix>_curve.dat      R, alpha, growth rate, N, gain, and total amplitude.
  <prefix>_summary.txt    Parameters and the requested-radius result.
""")
end

function parse_cli(args)
    cfg = Config()
    lower_radius = nothing
    positional_target_seen = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif !startswith(arg, "--")
            positional_target_seen && error("Only one positional target radius is allowed")
            cfg.target_radius = parse(Float64, arg)
            positional_target_seen = true
            i += 1
            continue
        end

        key, value = if occursin('=', arg)
            split(arg[3:end], '='; limit = 2)
        else
            i == length(args) && error("Missing value after $arg")
            (arg[3:end], args[i + 1])
        end
        if !occursin('=', arg)
            i += 1
        end

        if key == "target-radius"
            cfg.target_radius = parse(Float64, value)
        elseif key == "lower-radius"
            lower_radius = parse(Float64, value)
        elseif key == "Re"
            cfg.Re_h = parse(Int, value)
        elseif key == "mass-flux"
            cfg.mass_flux = parse(Float64, value)
        elseif key == "n"
            cfg.azimuthal_mode = parse(Int, value)
        elseif key == "omega-bar"
            cfg.omega_bar = parse(Float64, value)
        elseif key == "N-cheb"
            cfg.N_cheb = parse(Int, value)
        elseif key == "scan-start"
            cfg.scan_start = parse(Float64, value)
        elseif key == "scan-step"
            cfg.scan_step = parse(Float64, value)
        elseif key == "integration-step"
            cfg.integration_step = parse(Float64, value)
        elseif key == "alpha-real"
            cfg.alpha_guess_real = parse(Float64, value)
        elseif key == "alpha-imag"
            cfg.alpha_guess_imag = parse(Float64, value)
        elseif key == "candidates"
            cfg.candidate_count = parse(Int, value)
        elseif key == "roughness-height"
            cfg.roughness_height = parse(Float64, value)
        elseif key == "roughness-ls"
            cfg.roughness_localization = parse(Float64, value)
        elseif key == "output-prefix"
            cfg.output_prefix = value
        else
            error("Unknown option --$key. Run with --help for valid options.")
        end
        i += 1
    end

    # Keep the notebook convention h_r=1/Re_h unless the user overrides it.
    if cfg.roughness_height === nothing
        cfg.roughness_height = 1 / cfg.Re_h
    end
    return cfg, lower_radius
end

function validate(cfg, lower_radius)
    cfg.mode in (1, 3) || error("Only mode=1 (cavity) or mode=3 (single disk) is supported")
    cfg.N_cheb >= 20 || error("N_cheb must be at least 20")
    cfg.scan_step > 0 || error("scan_step must be positive")
    cfg.integration_step > 0 || error("integration_step must be positive")
    cfg.target_radius > cfg.scan_start || error("target_radius must exceed scan_start")
    cfg.candidate_count >= 1 || error("candidate_count must be positive")
    cfg.roughness_height > 0 || error("roughness_height must be positive")
    cfg.roughness_localization > 0 || error("roughness_localization must be positive")
    if lower_radius !== nothing
        lower_radius < cfg.target_radius || error("lower_radius must be smaller than target_radius")
    end
end

function cheb_mass_matrix(N_cheb::Int, mode::Int)
    weights = zeros(N_cheb + 1)
    for j in 0:N_cheb
        series_sum = 0.0
        for k in 1:floor(Int, N_cheb / 2)
            term = 2 / (1 - (2k)^2) * cos(2k * j * pi / N_cheb)
            series_sum += 2k == N_cheb ? 0.5 * term : term
        end
        endpoint_factor = (j == 0 || j == N_cheb) ? 1.0 : 2.0
        weights[j + 1] = endpoint_factor / N_cheb * (1 + series_sum)
    end

    full_mass = kron(I(4), Diagonal(0.5 .* weights))
    constrained = if mode == 1
        (1, N_cheb + 1, N_cheb + 2, 2N_cheb + 2,
         2N_cheb + 3, 3N_cheb + 3)
    else
        (1, N_cheb + 1, N_cheb + 2, 2N_cheb + 2,
         2N_cheb + 3, 3N_cheb + 3, 4N_cheb + 4)
    end
    free = setdiff(1:size(full_mass, 1), constrained)
    return full_mass[free, free]
end

function full_mode(vector, N_cheb)
    matrix = reshape(copy(vector), :, 1)
    return CRC_STA.eig_full(matrix, N_cheb, 1)
end

function normalized_overlap(a, b)
    denominator = norm(a) * norm(b)
    denominator == 0 && return 0.0
    return abs(dot(a, b)) / denominator
end

function choose_candidate(values, vectors, shift, previous_vector)
    if previous_vector === nothing
        return argmin(abs.(values .- shift))
    end
    spectral_scale = max(abs(shift), 0.1)
    scores = similar(real.(values))
    for j in eachindex(values)
        spectral_distance = abs(values[j] - shift) / spectral_scale
        overlap_penalty = 1 - normalized_overlap(previous_vector, vectors[:, j])
        scores[j] = spectral_distance + 0.25 * overlap_penalty
    end
    return argmin(scores)
end

function direct_solution(
    radius, shift, previous_vector, cfg, F, G, H, D, D2
)
    beta = cfg.azimuthal_mode / radius
    omega = cfg.omega_bar / radius
    cof = CRC_STA.Spatial_mode_BEK1(
        F, G .- 1, H, radius, cfg.N_cheb, D, D2, cfg.Re_h
    )
    L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(
        cof, D, D2, beta, omega, radius
    )
    L0, L1, L2 = CRC_STA.boudary_condition(
        L0_raw, L1_raw, L2_raw, cfg.N_cheb, cfg.mode
    )
    problem = PEP([L0, L1, L2])
    values, vectors = try
        iar(
            problem;
            σ = shift,
            neigs = cfg.candidate_count,
            maxit = 500,
            tol = cfg.eig_tolerance,
        )
    catch exception
        cfg.candidate_count == 1 && rethrow(exception)
        @warn "Multiple-root solve failed; retrying the target root only" radius exception
        iar(
            problem;
            σ = shift,
            neigs = 1,
            maxit = 750,
            tol = cfg.eig_tolerance,
        )
    end
    index = choose_candidate(values, vectors, shift, previous_vector)
    return (
        R = radius,
        alpha = values[index],
        vector = copy(vectors[:, index]),
        beta = beta,
        omega = omega,
        cof = cof,
        L1 = L1,
        L2 = L2,
    )
end

function radial_points(start_radius, end_radius, step)
    points = collect(start_radius:step:end_radius)
    if isempty(points) || points[end] < end_radius - 10eps(end_radius)
        push!(points, end_radius)
    elseif points[end] > end_radius + 10eps(end_radius)
        points[end] = end_radius
    end
    return points
end

function find_lower_neutral(cfg, F, G, H, D, D2)
    shift = complex(cfg.alpha_guess_real, cfg.alpha_guess_imag)
    previous = direct_solution(cfg.scan_start, shift, nothing, cfg, F, G, H, D, D2)
    @printf("neutral scan: R=%9.4f  alpha=% .9f %+.9fi\n",
            previous.R, real(previous.alpha), imag(previous.alpha))
    imag(previous.alpha) >= 0 || error(
        "The selected branch is already unstable at scan_start=$(cfg.scan_start). " *
        "Use a smaller --scan-start or a different initial alpha shift."
    )

    scan_radii = radial_points(
        cfg.scan_start + cfg.scan_step, cfg.target_radius, cfg.scan_step
    )
    for radius in scan_radii
        current = direct_solution(
            radius, previous.alpha, previous.vector, cfg, F, G, H, D, D2
        )
        @printf("neutral scan: R=%9.4f  alpha=% .9f %+.9fi\n",
                current.R, real(current.alpha), imag(current.alpha))
        if imag(previous.alpha) >= 0 && imag(current.alpha) <= 0
            return refine_neutral(previous, current, cfg, F, G, H, D, D2)
        end
        previous = current
    end
    error("No lower-branch neutral crossing was found before target_radius=$(cfg.target_radius)")
end

function refine_neutral(stable, unstable, cfg, F, G, H, D, D2)
    lo, hi = stable, unstable
    best = abs(imag(lo.alpha)) < abs(imag(hi.alpha)) ? lo : hi
    for iteration in 1:cfg.neutral_max_iterations
        ai_lo, ai_hi = imag(lo.alpha), imag(hi.alpha)
        fraction = ai_lo / (ai_lo - ai_hi)
        fraction = clamp(fraction, 0.1, 0.9)
        radius = lo.R + fraction * (hi.R - lo.R)
        shift = lo.alpha + fraction * (hi.alpha - lo.alpha)
        reference = fraction <= 0.5 ? lo : hi
        current = direct_solution(
            radius, shift, reference.vector, cfg, F, G, H, D, D2
        )
        if abs(imag(current.alpha)) < abs(imag(best.alpha))
            best = current
        end
        @printf("neutral refine %2d: R=%12.7f  alpha_i=%+.4e\n",
                iteration, radius, imag(current.alpha))
        if abs(imag(current.alpha)) <= cfg.neutral_alpha_tolerance ||
           hi.R - lo.R <= cfg.neutral_radius_tolerance
            return current
        elseif imag(current.alpha) > 0
            lo = current
        else
            hi = current
        end
    end
    @warn "Neutral refinement reached the iteration limit" best.R imag(best.alpha)
    return best
end

function solution_at_known_lower(lower_radius, cfg, F, G, H, D, D2)
    shift = complex(cfg.alpha_guess_real, cfg.alpha_guess_imag)
    solution = direct_solution(
        lower_radius, shift, nothing, cfg, F, G, H, D, D2
    )
    if abs(imag(solution.alpha)) > 5e-4
        @warn "The supplied lower radius is not close to neutral" lower_radius imag(solution.alpha)
    end
    return solution
end

function normalize_pair(direct_vector, adjoint_vector, N_cheb)
    direct_fields = full_mode(direct_vector, N_cheb)
    adjoint_fields = full_mode(adjoint_vector, N_cheb)
    direct_scale = maximum(abs.(direct_fields[2]))
    adjoint_scale = maximum(abs.(adjoint_fields[1]))
    direct_scale > 0 || error("Direct-mode normalization scale is zero")
    adjoint_scale > 0 || error("Adjoint-mode normalization scale is zero")
    return direct_vector / direct_scale, adjoint_vector / adjoint_scale
end

function initial_receptivity(neutral, cfg, F, G, D, D2, mass_matrix)
    A0_raw, A1_raw, A2_raw = CRC_STA.assemble_adjmat(
        neutral.cof, D, D2, neutral.beta, neutral.omega, neutral.R
    )
    A0, A1, A2 = CRC_STA.boudary_condition(
        A0_raw, A1_raw, A2_raw, cfg.N_cheb, cfg.mode
    )
    adjoint_values, adjoint_vectors = iar(
        PEP([A0, A1, A2]);
        σ = neutral.alpha,
        neigs = 1,
        maxit = 500,
        tol = cfg.eig_tolerance,
    )

    # The transpose problem preserves the analytic polynomial in alpha. Its
    # complex conjugate is the Hermitian-adjoint mode used in the inner product.
    direct_vector, adjoint_transpose = normalize_pair(
        neutral.vector, adjoint_vectors[:, 1], cfg.N_cheb
    )
    alpha_adjoint = conj(adjoint_values[1])
    adjoint_hermitian = conj.(adjoint_transpose)

    Q = adjoint_hermitian' * mass_matrix *
        (neutral.L1 + (neutral.alpha + conj(alpha_adjoint)) * neutral.L2) *
        direct_vector

    adjoint_fields = full_mode(adjoint_hermitian, cfg.N_cheb)
    roughness_transform = cfg.roughness_height *
        exp(-neutral.alpha^2 / (4cfg.roughness_localization)) *
        sqrt(pi / cfg.roughness_localization)
    u_wall = -(D * F)[1] * roughness_transform
    v_wall = -(D * G)[1] * roughness_transform
    boundary_projection = (
        conj((D * adjoint_fields[1])[1]) * u_wall +
        conj((D * adjoint_fields[2])[1]) * v_wall
    ) / (neutral.R * sqrt(cfg.Re_h))

    C_r = -im * boundary_projection / Q
    return (
        C_r = C_r,
        Q = Q,
        boundary_projection = boundary_projection,
        alpha_adjoint = alpha_adjoint,
    )
end

function integrate_branch(neutral, cfg, F, G, H, D, D2, C_r)
    radii = radial_points(neutral.R, cfg.target_radius, cfg.integration_step)
    solutions = Vector{typeof(neutral)}(undef, length(radii))
    solutions[1] = neutral
    N_factor = zeros(length(radii))
    phase_integral = zeros(length(radii))
    complex_integral = zeros(ComplexF64, length(radii))

    for j in 2:length(radii)
        previous = solutions[j - 1]
        current = direct_solution(
            radii[j], previous.alpha, previous.vector, cfg, F, G, H, D, D2
        )
        solutions[j] = current
        delta_R = current.R - previous.R
        alpha_average = 0.5 * (previous.alpha + current.alpha)
        complex_integral[j] = complex_integral[j - 1] + alpha_average * delta_R
        N_factor[j] = -imag(complex_integral[j])
        phase_integral[j] = real(complex_integral[j])
        @printf(
            "branch: R=%9.4f  alpha=% .9f %+.9fi  -alpha_i=%+.4e  N=%+.6f\n",
            current.R, real(current.alpha), imag(current.alpha),
            -imag(current.alpha), N_factor[j]
        )
    end

    gains = exp.(N_factor)
    total_complex = C_r .* exp.(im .* complex_integral)
    total_magnitude = abs.(C_r) .* gains
    return solutions, N_factor, phase_integral, gains, total_complex, total_magnitude
end

function write_outputs(
    cfg, neutral, receptivity, solutions, N_factor, phase_integral,
    gains, total_complex, total_magnitude
)
    prefix = abspath(cfg.output_prefix)
    mkpath(dirname(prefix))
    curve_path = prefix * "_curve.dat"
    summary_path = prefix * "_summary.txt"

    open(curve_path, "w") do io
        println(io, "# R beta omega alpha_r alpha_i growth_rate N_factor gain phase_integral C_r_initial_abs A_abs A_real A_imag")
        for j in eachindex(solutions)
            sol = solutions[j]
            @printf(
                io,
                "%.12g %.12g %.12g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                sol.R, sol.beta, sol.omega, real(sol.alpha), imag(sol.alpha),
                -imag(sol.alpha), N_factor[j], gains[j], phase_integral[j],
                abs(receptivity.C_r), total_magnitude[j],
                real(total_complex[j]), imag(total_complex[j])
            )
        end
    end

    final = solutions[end]
    open(summary_path, "w") do io
        println(io, "Total-amplitude calculation")
        println(io, "Re_h = ", cfg.Re_h)
        println(io, "a_s = ", cfg.mass_flux)
        println(io, "n = ", cfg.azimuthal_mode)
        println(io, "omega_bar = ", cfg.omega_bar)
        println(io, "N_cheb = ", cfg.N_cheb, " (", cfg.N_cheb + 1, " collocation points)")
        println(io, "R_lower = ", neutral.R)
        println(io, "alpha_lower = ", repr(neutral.alpha))
        println(io, "alpha_adjoint_lower = ", repr(receptivity.alpha_adjoint))
        println(io, "C_r_lower = ", repr(receptivity.C_r))
        println(io, "abs_C_r_lower = ", abs(receptivity.C_r))
        println(io, "Q_lower = ", repr(receptivity.Q))
        println(io, "boundary_projection_lower = ", repr(receptivity.boundary_projection))
        println(io, "R_target = ", final.R)
        println(io, "alpha_target = ", repr(final.alpha))
        println(io, "growth_rate_target = ", -imag(final.alpha))
        println(io, "N_target = ", N_factor[end])
        println(io, "gain_target = ", gains[end])
        println(io, "log_abs_A_target = ", log(abs(receptivity.C_r)) + N_factor[end])
        println(io, "abs_A_target = ", total_magnitude[end])
        println(io, "A_target_raw_phase = ", repr(total_complex[end]))
        println(io, "curve_file = ", curve_path)
        println(io, "note = A_target_raw_phase depends on the arbitrary eigenvector phase; abs_A_target is the phase-invariant result.")
    end
    return curve_path, summary_path
end

function main(args)
    cfg, supplied_lower_radius = parse_cli(args)
    validate(cfg, supplied_lower_radius)

    println("Solving base flow...")
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(
        u0, v0, w0, z, cfg.N_cheb, cfg.mode
    )
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    neutral = supplied_lower_radius === nothing ?
        find_lower_neutral(cfg, F, G, H, D, D2) :
        solution_at_known_lower(supplied_lower_radius, cfg, F, G, H, D, D2)
    @printf("\nLower neutral point: R_l=%.10f, alpha=% .12f %+.12fi\n",
            neutral.R, real(neutral.alpha), imag(neutral.alpha))

    receptivity = initial_receptivity(neutral, cfg, F, G, D, D2, mass_matrix)
    @printf("Initial receptivity: |C_r(R_l)|=%.12e\n\n", abs(receptivity.C_r))

    solutions, N_factor, phase_integral, gains, total_complex,
    total_magnitude = integrate_branch(
        neutral, cfg, F, G, H, D, D2, receptivity.C_r
    )
    curve_path, summary_path = write_outputs(
        cfg, neutral, receptivity, solutions, N_factor, phase_integral,
        gains, total_complex, total_magnitude
    )

    final = solutions[end]
    println("\nRequested-radius result")
    @printf("  R_target       = %.10f\n", final.R)
    @printf("  growth rate    = %.12e\n", -imag(final.alpha))
    @printf("  N factor       = %.12e\n", N_factor[end])
    @printf("  gain exp(N)    = %.12e\n", gains[end])
    @printf("  |C_r(R_l)|     = %.12e\n", abs(receptivity.C_r))
    @printf("  |A(R_target)|  = %.12e\n", total_magnitude[end])
    println("  curve output   = ", curve_path)
    println("  summary output = ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
