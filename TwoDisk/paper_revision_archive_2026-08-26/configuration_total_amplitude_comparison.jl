using LinearAlgebra
using Printf
using NonlinearEigenproblems

include(joinpath(@__DIR__, "configuration_effect_validation.jl"))

const AMP_OUTPUT_DIR = joinpath(@__DIR__, "configuration_total_amplitude_results")
const R_TARGET = 500.0

function single_direct_point(radius, shift; N=99, n=30.0)
    params = FlowParameters(1000, -1.0, 0.0, radius, n, 0.0, N, 2)
    baseflow, grid = solve_baseflow(params)
    F, G, H = baseflow.F, baseflow.G, baseflow.H
    D, D2 = grid.D, grid.D2
    beta = n / radius
    Ro = -1.0
    Co = 2 - Ro - Ro^2
    cof = assemble_coeff_matrix_BEK2(F, G, H, radius, N, D, D2, Ro, Co)
    L0, L1, L2 = assemble_direct_matrices(cof, D, D2, beta, 0.0, radius)
    L0, L1, L2 = apply_boundary_conditions!(L0, L1, L2, N)
    alpha, vector = solve_eigenvalue_problem(L0, L1, L2, ComplexF64(shift), 1)
    return (; R=radius, alpha=alpha[1], vector=copy(vector[:, 1]))
end

function find_single_neutral(; N=99, n=30.0)
    lo = single_direct_point(300.0, 0.51 + 0.01im; N=N, n=n)
    hi = single_direct_point(330.0, lo.alpha; N=N, n=n)
    imag(lo.alpha) > 0 || error("Single-disk neutral bracket starts unstable")
    imag(hi.alpha) < 0 || error("Single-disk neutral bracket ends stable")
    best = abs(imag(lo.alpha)) < abs(imag(hi.alpha)) ? lo : hi
    for _ in 1:30
        fraction = imag(lo.alpha) / (imag(lo.alpha) - imag(hi.alpha))
        fraction = clamp(fraction, 0.1, 0.9)
        radius = lo.R + fraction * (hi.R - lo.R)
        shift = lo.alpha + fraction * (hi.alpha - lo.alpha)
        current = single_direct_point(radius, shift; N=N, n=n)
        abs(imag(current.alpha)) < abs(imag(best.alpha)) && (best = current)
        if abs(imag(current.alpha)) < 1e-10 || hi.R - lo.R < 1e-6
            return current
        elseif imag(current.alpha) > 0
            lo = current
        else
            hi = current
        end
    end
    return best
end

function single_amplitude_curve(neutral, C_r; N=99, n=30.0, step=1.0)
    radii = radial_points(neutral.R, R_TARGET, step)
    solutions = Vector{typeof(neutral)}(undef, length(radii))
    solutions[1] = neutral
    N_factor = zeros(length(radii))
    for j in 2:length(radii)
        previous = solutions[j - 1]
        current = single_direct_point(radii[j], previous.alpha; N=N, n=n)
        solutions[j] = current
        dR = current.R - previous.R
        N_factor[j] = N_factor[j - 1] -
                      0.5 * (imag(previous.alpha) + imag(current.alpha)) * dR
        if mod(j, 25) == 0 || j == length(radii)
            @printf("single branch: R=%8.3f alpha=% .8f%+.8fi N=% .7f\n",
                    current.R, real(current.alpha), imag(current.alpha), N_factor[j])
        end
    end
    gain = exp.(N_factor)
    amplitude = abs(C_r) .* gain
    return solutions, N_factor, gain, amplitude
end

function integrate_subset(radii, alpha_i, stride)
    ids = collect(1:stride:length(radii))
    ids[end] == length(radii) || push!(ids, length(radii))
    value = 0.0
    for k in 2:length(ids)
        i0, i1 = ids[k - 1], ids[k]
        value -= 0.5 * (alpha_i[i0] + alpha_i[i1]) *
                 (radii[i1] - radii[i0])
    end
    return value
end

function interpolate_value(radii, values, target)
    target <= radii[1] && return values[1]
    target >= radii[end] && return values[end]
    j = searchsortedfirst(radii, target)
    fraction = (target - radii[j - 1]) / (radii[j] - radii[j - 1])
    return values[j - 1] + fraction * (values[j] - values[j - 1])
end

function order_one_radius(radii, amplitude)
    j = findfirst(x -> x >= 1, amplitude)
    j === nothing && return NaN
    j == 1 && return radii[1]
    log0, log1 = log(amplitude[j - 1]), log(amplitude[j])
    fraction = -log0 / (log1 - log0)
    return radii[j - 1] + fraction * (radii[j] - radii[j - 1])
end

function write_zone(io, name, solutions, N_factor, gain, amplitude, C_r)
    println(io, "ZONE T=\"", name, "\", I=", length(solutions), ", F=POINT")
    for j in eachindex(solutions)
        solution = solutions[j]
        @printf(io, "%.12f %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
                solution.R, real(solution.alpha), imag(solution.alpha),
                -imag(solution.alpha), N_factor[j], gain[j], abs(C_r), amplitude[j])
    end
end

function main()
    mkpath(AMP_OUTPUT_DIR)
    N = 99
    n = 30.0

    cavity_ctx = cavity_context(N)
    cavity_ctx.cfg.target_radius = R_TARGET
    cavity_ctx.cfg.scan_start = 200.0
    cavity_ctx.cfg.scan_step = 2.0
    cavity_ctx.cfg.integration_step = 1.0

    println("Finding the cavity lower neutral point...")
    cavity_neutral = find_lower_neutral(
        cavity_ctx.cfg, cavity_ctx.F, cavity_ctx.G, cavity_ctx.H,
        cavity_ctx.D, cavity_ctx.D2
    )
    cavity_receptivity_native = cavity_point(
        cavity_neutral.R, cavity_neutral.alpha, cavity_neutral.vector, cavity_ctx
    )
    cavity_Cr = cavity_receptivity_native.Cr / FOURIER_GAUSSIAN_FACTOR

    println("Integrating the cavity branch...")
    cavity_solutions, cavity_N, _, cavity_gain, _, cavity_A = integrate_branch(
        cavity_neutral, cavity_ctx.cfg, cavity_ctx.F, cavity_ctx.G,
        cavity_ctx.H, cavity_ctx.D, cavity_ctx.D2, cavity_Cr
    )

    println("Finding the isolated-disk lower neutral point...")
    single_neutral = find_single_neutral(N=N, n=n)
    single_receptivity = solve_point(
        R=single_neutral.R, n=n, N=N, shift=single_neutral.alpha
    )
    single_Cr = single_receptivity.Cr_thomas

    println("Integrating the isolated-disk branch...")
    single_solutions, single_N, single_gain, single_A = single_amplitude_curve(
        single_neutral, single_Cr; N=N, n=n, step=1.0
    )

    curve_path = joinpath(AMP_OUTPUT_DIR, "total_amplitude_comparison_N99.dat")
    open(curve_path, "w") do io
        println(io, "TITLE=\"Total amplitude from each lower neutral point\"")
        println(io, "VARIABLES=\"R\",\"alpha_r\",\"alpha_i\",\"growth_rate\",\"N_factor\",\"gain\",\"Cr_initial\",\"A_abs\"")
        write_zone(io, "rotor-stator cavity", cavity_solutions,
                   cavity_N, cavity_gain, cavity_A, cavity_Cr)
        write_zone(io, "isolated rotating disk", single_solutions,
                   single_N, single_gain, single_A, single_Cr)
    end

    cavity_radii = [s.R for s in cavity_solutions]
    cavity_ai = [imag(s.alpha) for s in cavity_solutions]
    single_radii = [s.R for s in single_solutions]
    single_ai = [imag(s.alpha) for s in single_solutions]

    summary_path = joinpath(AMP_OUTPUT_DIR, "total_amplitude_comparison_summary.txt")
    open(summary_path, "w") do io
        println(io, "Parameters: a_s=0, n=30, omega_bar=0, c^2=1, N=99")
        println(io, "Both coefficients use the same unitary Gaussian Fourier normalization.")
        println(io)
        println(io, "Cavity lower neutral R = ", cavity_neutral.R)
        println(io, "Cavity alpha at lower neutral = ", repr(cavity_neutral.alpha))
        println(io, "Cavity corrected |Cr(R_l)| = ", cavity_Cr)
        println(io, "Cavity N(R_target) = ", cavity_N[end])
        println(io, "Cavity |A(R_target)| = ", cavity_A[end])
        println(io, "Cavity order-one radius = ", order_one_radius(cavity_radii, cavity_A))
        println(io)
        println(io, "Single-disk lower neutral R = ", single_neutral.R)
        println(io, "Single-disk alpha at lower neutral = ", repr(single_neutral.alpha))
        println(io, "Single-disk corrected |Cr(R_l)| = ", single_Cr)
        println(io, "Single-disk N(R_target) = ", single_N[end])
        println(io, "Single-disk |A(R_target)| = ", single_A[end])
        println(io, "Single-disk order-one radius = ", order_one_radius(single_radii, single_A))
        println(io)
        println(io, "Common target radius = ", R_TARGET)
        println(io, "Target amplitude ratio disk/cavity = ", single_A[end] / cavity_A[end])
        println(io)
        println(io, "Integration-step check derived from the Delta R=1 curves")
        for stride in (1, 2, 4)
            Nc = integrate_subset(cavity_radii, cavity_ai, stride)
            Ns = integrate_subset(single_radii, single_ai, stride)
            println(io, "Delta R approximately ", stride,
                    ": cavity Ntarget=", Nc, " Atarget=", cavity_Cr*exp(Nc),
                    "; single Ntarget=", Ns, " Atarget=", single_Cr*exp(Ns))
        end
        println(io)
        println(io, "curve_file = ", curve_path)
    end

    println("Wrote ", curve_path)
    println("Wrote ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
