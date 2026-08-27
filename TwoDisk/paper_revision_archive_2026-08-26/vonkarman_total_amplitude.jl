using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "configuration_total_amplitude_comparison.jl"))

const VONKARMAN_OUTPUT_DIR = joinpath(@__DIR__, "vonkarman_total_amplitude_results")
const VONKARMAN_TARGET_RADIUS = 500.0
const VONKARMAN_RADIAL_STEP = 2.0

function validate_malik_point()
    reference_R = 285.36
    reference_beta = 0.07759
    reference_alpha = 0.38482
    azimuthal_mode = reference_beta * reference_R
    results = NamedTuple[]
    for N in (79, 99, 129)
        result = solve_point(
            R=reference_R,
            n=azimuthal_mode,
            N=N,
            shift=reference_alpha + 0im,
        )
        push!(results, (
            N=N,
            alpha=result.alpha,
            pairing_error=result.pairing_error,
            direct_residual=result.direct_residual,
            adjoint_residual=result.adjoint_residual,
        ))
    end
    primary = only(filter(row -> row.N == 99, results))
    relative_error = abs(real(primary.alpha) - reference_alpha) / reference_alpha
    relative_error < 5e-3 || error("Malik alpha_r benchmark failed")
    abs(imag(primary.alpha)) < 1e-3 || error("Malik neutral offset benchmark failed")
    primary.direct_residual < 1e-10 || error("Malik direct residual benchmark failed")
    return (; reference_R, reference_beta, reference_alpha, results, relative_error)
end

function write_vonkarman_curve(path, solutions, N_factor, gain, amplitude, C_r)
    open(path, "w") do io
        println(io, "TITLE=\"von Karman total amplitude from the lower neutral point\"")
        println(io,
            "VARIABLES=\"R\",\"beta\",\"alpha_r\",\"alpha_i\"," *
            "\"growth_rate\",\"N_factor\",\"gain\",\"Cr_initial\",\"A_abs\""
        )
        println(io, "DATASETAUXDATA n=\"30\"")
        println(io, "DATASETAUXDATA omega=\"0\"")
        println(io, "DATASETAUXDATA c_squared=\"1\"")
        println(io, "DATASETAUXDATA N_cheb=\"99\"")
        println(io, "DATASETAUXDATA radial_step=\"$(VONKARMAN_RADIAL_STEP)\"")
        println(io, "ZONE T=\"isolated von Karman disk\", I=$(length(solutions)), F=POINT")
        for index in eachindex(solutions)
            solution = solutions[index]
            @printf(io,
                "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                solution.R, 30.0 / solution.R, real(solution.alpha),
                imag(solution.alpha), -imag(solution.alpha), N_factor[index],
                gain[index], abs(C_r), amplitude[index],
            )
        end
    end
end

function main_vonkarman_amplitude()
    LinearAlgebra.BLAS.set_num_threads(1)
    mkpath(VONKARMAN_OUTPUT_DIR)

    benchmark = validate_malik_point()
    println("Malik Type-I point accepted by the project benchmark criteria.")

    neutral = find_single_neutral(N=99, n=30.0)
    receptivity = solve_point(R=neutral.R, n=30.0, N=99, shift=neutral.alpha)
    C_r = receptivity.Cr
    solutions, N_factor, gain, amplitude = single_amplitude_curve(
        neutral, C_r; N=99, n=30.0, step=VONKARMAN_RADIAL_STEP
    )

    curve_path = joinpath(VONKARMAN_OUTPUT_DIR, "vonkarman_A_N99.dat")
    write_vonkarman_curve(curve_path, solutions, N_factor, gain, amplitude, C_r)

    radii = [solution.R for solution in solutions]
    alpha_i = [imag(solution.alpha) for solution in solutions]
    summary_path = joinpath(VONKARMAN_OUTPUT_DIR, "vonkarman_A_N99_summary.txt")
    open(summary_path, "w") do io
        println(io, "Malik Type-I stationary benchmark")
        println(io, "reference R = ", benchmark.reference_R)
        println(io, "reference beta = ", benchmark.reference_beta)
        println(io, "reference alpha_r = ", benchmark.reference_alpha)
        for row in benchmark.results
            println(io, "N=", row.N,
                    " alpha=", repr(row.alpha),
                    " pairing_error=", row.pairing_error,
                    " direct_residual=", row.direct_residual,
                    " adjoint_residual=", row.adjoint_residual)
        end
        println(io, "N=99 relative alpha_r error = ", benchmark.relative_error)
        println(io)
        println(io, "von Karman amplitude evolution")
        println(io, "parameters: n=30, omega=0, c^2=1, N=99, Delta R=2")
        println(io, "Fourier convention: hhat(alpha)=integral h(s) exp(-i alpha s) ds")
        println(io, "lower neutral R = ", neutral.R)
        println(io, "lower neutral beta = ", 30.0 / neutral.R)
        println(io, "lower neutral alpha = ", repr(neutral.alpha))
        println(io, "initial |Cr| = ", abs(C_r))
        println(io, "target R = ", VONKARMAN_TARGET_RADIUS)
        println(io, "N(target) = ", N_factor[end])
        println(io, "gain(target) = ", gain[end])
        println(io, "|A(target)| = ", amplitude[end])
        println(io, "order-one radius = ", order_one_radius(radii, amplitude))
        println(io, "direct residual at lower neutral = ", receptivity.direct_residual)
        println(io, "adjoint residual at lower neutral = ", receptivity.adjoint_residual)
        println(io, "pairing error at lower neutral = ", receptivity.pairing_error)
        println(io)
        println(io, "integration-step check")
        for stride in (1, 2, 4)
            value = integrate_subset(radii, alpha_i, stride)
            println(io, "Delta R approximately ", VONKARMAN_RADIAL_STEP * stride,
                    ": Ntarget=", value,
                    " Atarget=", abs(C_r) * exp(value))
        end
        println(io, "curve_file = ", curve_path)
    end

    println("curve: ", curve_path)
    println("summary: ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_vonkarman_amplitude()
end
