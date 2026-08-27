using LinearAlgebra
using DelimitedFiles
using Printf

include("total_amplitude.jl")

const C6_RADIUS_TARGET = 570.0
const C6_LOCAL_REFERENCE = 360.0

function refine_extremum(xs, ys, index)
    if index == 1 || index == length(xs)
        return xs[index], ys[index]
    end
    x1, x2, x3 = xs[index - 1:index + 1]
    y1, y2, y3 = ys[index - 1:index + 1]
    spacing = x2 - x1
    denominator = y1 - 2y2 + y3
    abs(denominator) < eps(Float64) && return x2, y2
    offset = 0.5 * (y1 - y3) / denominator
    x_peak = x2 + offset * spacing
    y_peak = y2 - 0.25 * (y1 - y3) * offset
    return x_peak, y_peak
end

function interpolate_complex(xs, ys, x)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    upper = searchsortedfirst(xs, x)
    lower = upper - 1
    weight = (x - xs[lower]) / (xs[upper] - xs[lower])
    return (1 - weight) * ys[lower] + weight * ys[upper]
end

function cumulative_integral(xs, ys)
    values = zeros(ComplexF64, length(xs))
    for j in 2:length(xs)
        values[j] = values[j - 1] +
            0.5 * (ys[j - 1] + ys[j]) * (xs[j] - xs[j - 1])
    end
    return values
end

function solution_interpolated(radius, branch, cfg, F, G, H, D, D2)
    radii = getfield.(branch, :R)
    index = argmin(abs.(radii .- radius))
    reference = branch[index]
    return direct_solution(
        radius, reference.alpha, reference.vector, cfg, F, G, H, D, D2
    )
end

function calculate_receptivity_curve(branch, cfg, F, G, H, D, D2, mass_matrix)
    curve = zeros(length(branch), 5)
    for j in eachindex(branch)
        result = initial_receptivity(branch[j], cfg, F, G, D, D2, mass_matrix)
        curve[j, :] .= (
            branch[j].R,
            abs(result.C_r),
            abs(result.Q),
            abs(result.boundary_projection),
            abs(result.alpha_adjoint - conj(branch[j].alpha)),
        )
        @printf(
            "receptivity: R=%9.4f  |Cr|=%.9e  |Q|=%.9e  |BC|=%.9e  adj_err=%.3e\n",
            curve[j, 1], curve[j, 2], curve[j, 3], curve[j, 4], curve[j, 5]
        )
    end
    return curve
end

function case_result(
    label, radius, branch, integral, cfg, F, G, H, D, D2, mass_matrix
)
    solution = solution_interpolated(radius, branch, cfg, F, G, H, D, D2)
    receptivity = initial_receptivity(solution, cfg, F, G, D, D2, mass_matrix)
    radii = getfield.(branch, :R)
    integral_start = interpolate_complex(radii, integral, radius)
    downstream_integral = integral[end] - integral_start
    n_factor = -imag(downstream_integral)
    gain = exp(n_factor)
    total_amplitude = abs(receptivity.C_r) * gain
    return (
        label = label,
        R = radius,
        alpha = solution.alpha,
        growth_rate = -imag(solution.alpha),
        C_r = receptivity.C_r,
        Q = receptivity.Q,
        BC = receptivity.boundary_projection,
        N = n_factor,
        gain = gain,
        A = total_amplitude,
        adjoint_error = abs(receptivity.alpha_adjoint - conj(solution.alpha)),
    )
end

function main()
    cfg = Config(
        target_radius = C6_RADIUS_TARGET,
        N_cheb = 99,
        scan_start = 240.0,
        scan_step = 2.0,
        integration_step = 2.0,
        candidate_count = 1,
        roughness_localization = 0.5,
        output_prefix = joinpath("c6_validation_results", "unused"),
    )
    cfg.roughness_height = 1 / cfg.Re_h

    output_directory = abspath("c6_validation_results")
    mkpath(output_directory)

    println("Solving base flow for C6 three-case validation...")
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    neutral = find_lower_neutral(cfg, F, G, H, D, D2)
    unit_receptivity = 1.0 + 0.0im
    branch, _, _, _, _, _ = integrate_branch(
        neutral, cfg, F, G, H, D, D2, unit_receptivity
    )
    radii = getfield.(branch, :R)
    alphas = getfield.(branch, :alpha)
    growth_rates = -imag.(alphas)
    integral = cumulative_integral(radii, alphas)

    growth_index = argmax(growth_rates)
    growth_peak_radius, growth_peak_value = refine_extremum(
        radii, growth_rates, growth_index
    )

    # A two-unit receptivity scan is sufficient to locate the broad extrema.
    scan_indices = collect(1:length(branch))
    receptivity_curve = calculate_receptivity_curve(
        branch[scan_indices], cfg, F, G, H, D, D2, mass_matrix
    )
    global_index = argmax(receptivity_curve[:, 2])
    global_peak_radius, global_peak_value = refine_extremum(
        receptivity_curve[:, 1], receptivity_curve[:, 2], global_index
    )

    cases = [
        case_result(
            "lower_neutral", neutral.R, branch, integral,
            cfg, F, G, H, D, D2, mass_matrix
        ),
        case_result(
            "maximum_growth_local_case", growth_peak_radius, branch, integral,
            cfg, F, G, H, D, D2, mass_matrix
        ),
        case_result(
            "global_receptivity_peak", global_peak_radius, branch, integral,
            cfg, F, G, H, D, D2, mass_matrix
        ),
    ]
    check_360 = case_result(
        "fixed_R_360_check", C6_LOCAL_REFERENCE, branch, integral,
        cfg, F, G, H, D, D2, mass_matrix
    )

    curve_path = joinpath(output_directory, "c6_branch_and_receptivity.dat")
    open(curve_path, "w") do io
        println(io, "# R alpha_r alpha_i growth_rate cumulative_N_from_Rl Cr_abs Q_abs BC_abs adjoint_eigenvalue_error")
        for j in eachindex(branch)
            @printf(
                io, "%.12g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                radii[j], real(alphas[j]), imag(alphas[j]), growth_rates[j],
                -imag(integral[j]), receptivity_curve[j, 2],
                receptivity_curve[j, 3], receptivity_curve[j, 4],
                receptivity_curve[j, 5]
            )
        end
    end

    case_path = joinpath(output_directory, "c6_three_cases.dat")
    open(case_path, "w") do io
        println(io, "# label R alpha_r alpha_i growth_rate Cr_abs Q_abs BC_abs N_to_570 gain_to_570 A_abs_at_570 adjoint_eigenvalue_error")
        for result in [cases; check_360]
            @printf(
                io, "%s %.12g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                result.label, result.R, real(result.alpha), imag(result.alpha),
                result.growth_rate, abs(result.C_r), abs(result.Q), abs(result.BC),
                result.N, result.gain, result.A, result.adjoint_error
            )
        end
    end

    summary_path = joinpath(output_directory, "c6_validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "C6 three-case program-validation summary")
        println(io, "Re_h = ", cfg.Re_h)
        println(io, "a_s = ", cfg.mass_flux)
        println(io, "n = ", cfg.azimuthal_mode)
        println(io, "omega_bar = ", cfg.omega_bar)
        println(io, "N_cheb = ", cfg.N_cheb)
        println(io, "roughness_height = ", cfg.roughness_height)
        println(io, "roughness_localization = ", cfg.roughness_localization)
        println(io, "R_target = ", cfg.target_radius)
        println(io, "R_lower = ", neutral.R)
        println(io, "growth_peak_radius = ", growth_peak_radius)
        println(io, "growth_peak_value = ", growth_peak_value)
        println(io, "local_case_definition = maximum growth-rate location near the reviewer-discussed local feature")
        println(io, "global_Cr_peak_radius = ", global_peak_radius)
        println(io, "global_Cr_peak_value = ", global_peak_value)
        println(io, "case_file = ", case_path)
        println(io, "curve_file = ", curve_path)
    end

    println("\nC6 landmarks")
    @printf("  lower neutral R       = %.8f\n", neutral.R)
    @printf("  maximum growth R      = %.8f, -alpha_i=%.10e\n", growth_peak_radius, growth_peak_value)
    println("  strict local Cr peak near R=360 was not detected; Cr remains increasing there")
    @printf("  global Cr peak R      = %.8f, |Cr|=%.10e\n", global_peak_radius, global_peak_value)
    println("\nC6 cases propagated to R=", cfg.target_radius)
    for result in [cases; check_360]
        @printf(
            "  %-26s R=%10.5f |Cr|=%.9e N=%+.9f gain=%.9e |A|=%.9e\n",
            result.label, result.R, abs(result.C_r), result.N,
            result.gain, result.A
        )
    end
    println("Results written to ", output_directory)
end

main()
