using DelimitedFiles
using LinearAlgebra
using Printf

include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")

const RES = parse(Int, get(ENV, "GROUP_VELOCITY_RES", "1000"))
const N_CHEB = parse(Int, get(ENV, "GROUP_VELOCITY_N", "199"))
const MODE = 1
const RO = -1.0
const NEWTON_STEP = 1.0e-4
const MATRIX_DERIVATIVE_STEP = 1.0e-3
const DERIVATIVE_STEPS = (2.0e-4, 1.0e-4, 5.0e-5, 2.5e-5)

const TABLE_V_CASES = [
    (a_s=-0.4, radius=362.0, beta=0.1595, alpha=0.225-0.130im, omega=-0.0318+0.0im),
    (a_s=-0.2, radius=452.0, beta=0.1445, alpha=0.233-0.137im, omega=-0.0251+0.0im),
    (a_s= 0.0, radius=570.0, beta=0.1520, alpha=0.262-0.143im, omega=-0.0247+0.0im),
    (a_s= 0.2, radius=732.0, beta=0.1550, alpha=0.270-0.148im, omega=-0.0238+0.0im),
    (a_s= 0.4, radius=962.0, beta=0.1525, alpha=0.290-0.159im, omega=-0.0220+0.0im),
]

function inverse_iteration(matrix_a, matrix_b, shift, initial_vector; tolerance=1.0e-13, maximum_iterations=40)
    eigenvalue = ComplexF64(shift)
    eigenvector = initial_vector / norm(initial_vector)
    for _ in 1:maximum_iterations
        previous = eigenvalue
        updated = (matrix_a - eigenvalue * matrix_b) \ (matrix_b * eigenvector)
        eigenvector = updated / norm(updated)
        eigenvalue = dot(eigenvector, matrix_a * eigenvector) / dot(eigenvector, matrix_b * eigenvector)
        abs(eigenvalue - previous) < tolerance && return eigenvalue, eigenvector
    end
    error("Inverse iteration failed to converge")
end

function initial_eigenpair(coefficient, differentiation, differentiation2, beta, alpha, radius, target)
    matrix_a, matrix_b = CRC_STA.assemble_time_mat(
        coefficient, differentiation, differentiation2, beta, alpha, radius, N_CHEB
    )
    decomposition = eigen(matrix_a, matrix_b)
    finite_indices = findall(index -> isfinite(decomposition.values[index]), eachindex(decomposition.values))
    selected = finite_indices[argmin(abs.(decomposition.values[finite_indices] .- target))]
    return decomposition.values[selected], decomposition.vectors[:, selected]
end

function continued_eigenpair(coefficient, differentiation, differentiation2, beta, alpha, radius, target, vector)
    matrix_a, matrix_b = CRC_STA.assemble_time_mat(
        coefficient, differentiation, differentiation2, beta, alpha, radius, N_CHEB
    )
    eigenvalue, eigenvector = inverse_iteration(matrix_a, matrix_b, target, vector)
    return eigenvalue, eigenvector, matrix_a, matrix_b
end

function finite_difference_derivatives(coefficient, differentiation, differentiation2, beta, alpha, radius, omega, vector, step)
    omega_plus, _, _, _ = continued_eigenpair(
        coefficient, differentiation, differentiation2, beta, alpha + step, radius, omega, vector
    )
    omega_minus, _, _, _ = continued_eigenpair(
        coefficient, differentiation, differentiation2, beta, alpha - step, radius, omega, vector
    )
    first = (omega_plus - omega_minus) / (2step)
    second = (omega_plus - 2omega + omega_minus) / step^2
    return first, second
end

function refine_saddle(coefficient, differentiation, differentiation2, beta, radius, alpha_initial, omega_initial)
    omega, vector = initial_eigenpair(
        coefficient, differentiation, differentiation2, beta, alpha_initial, radius, omega_initial
    )
    alpha = ComplexF64(alpha_initial)
    for iteration in 1:20
        omega, vector, _, _ = continued_eigenpair(
            coefficient, differentiation, differentiation2, beta, alpha, radius, omega, vector
        )
        first, second = finite_difference_derivatives(
            coefficient, differentiation, differentiation2, beta, alpha, radius, omega, vector, NEWTON_STEP
        )
        @printf("  iteration %2d: alpha=%+.12e%+.12ei omega=%+.12e%+.12ei |cg|=%.3e\n",
            iteration, real(alpha), imag(alpha), real(omega), imag(omega), abs(first))
        abs(first) < 1.0e-10 && return alpha, omega, vector
        abs(second) > 1.0e-12 || error("Vanishing second derivative during saddle refinement")
        correction = first / second
        abs(correction) < 0.1 || error("Newton correction left the target branch")
        alpha -= correction
    end
    error("Saddle refinement failed to converge")
end

function adjoint_group_velocity(coefficient, differentiation, differentiation2, beta, alpha, radius, omega, right_vector, step)
    _, _, matrix_a, matrix_b = continued_eigenpair(
        coefficient, differentiation, differentiation2, beta, alpha, radius, omega, right_vector
    )
    matrix_a_plus, matrix_b_plus = CRC_STA.assemble_time_mat(
        coefficient, differentiation, differentiation2, beta, alpha + step, radius, N_CHEB
    )
    matrix_a_minus, matrix_b_minus = CRC_STA.assemble_time_mat(
        coefficient, differentiation, differentiation2, beta, alpha - step, radius, N_CHEB
    )
    left_decomposition = eigen(adjoint(matrix_a), adjoint(matrix_b))
    finite_indices = findall(index -> isfinite(left_decomposition.values[index]), eachindex(left_decomposition.values))
    selected = finite_indices[argmin(abs.(left_decomposition.values[finite_indices] .- conj(omega)))]
    left_vector = left_decomposition.vectors[:, selected]
    derivative_a = (matrix_a_plus - matrix_a_minus) / (2step)
    derivative_b = (matrix_b_plus - matrix_b_minus) / (2step)
    numerator = dot(left_vector, (derivative_a - omega * derivative_b) * right_vector)
    denominator = dot(left_vector, matrix_b * right_vector)
    return numerator / denominator
end

function compute_case(case)
    velocity_r, velocity_theta, velocity_z, _, _, _ = CRC_BF.BaseFlow(RES, RO, case.a_s, MODE)
    differentiation, differentiation2, collocation = CRC_BF.Cheb(N_CHEB, MODE)
    base_r, base_theta, base_z = CRC_BF.interp(
        velocity_r, velocity_theta, velocity_z, collocation, N_CHEB, MODE
    )
    coefficient = CRC_STA.Spatial_mode_BEK1(
        base_r, base_theta .- 1, base_z, case.radius, N_CHEB,
        differentiation, differentiation2, RES
    )
    @printf("a_s=%+.1f R=%.1f beta=%.4f\n", case.a_s, case.radius, case.beta)
    table_omega, table_vector = initial_eigenpair(
        coefficient, differentiation, differentiation2, case.beta, case.alpha, case.radius, case.omega
    )
    table_finite_value, _ = finite_difference_derivatives(
        coefficient, differentiation, differentiation2, case.beta, case.alpha, case.radius,
        table_omega, table_vector, NEWTON_STEP
    )
    table_adjoint_value = adjoint_group_velocity(
        coefficient, differentiation, differentiation2, case.beta, case.alpha, case.radius,
        table_omega, table_vector, MATRIX_DERIVATIVE_STEP
    )
    alpha, omega, vector = refine_saddle(
        coefficient, differentiation, differentiation2, case.beta, case.radius, case.alpha, case.omega
    )
    finite_values = ComplexF64[]
    for step in DERIVATIVE_STEPS
        first, _ = finite_difference_derivatives(
            coefficient, differentiation, differentiation2, case.beta, alpha, case.radius, omega, vector, step
        )
        push!(finite_values, first)
    end
    adjoint_value = adjoint_group_velocity(
        coefficient, differentiation, differentiation2, case.beta, alpha, case.radius, omega, vector,
        MATRIX_DERIVATIVE_STEP
    )
    return (
        a_s=case.a_s, radius=case.radius, beta=case.beta, alpha=alpha, omega=omega,
        table_alpha=case.alpha, table_omega=table_omega,
        table_finite_value=table_finite_value, table_adjoint_value=table_adjoint_value,
        finite_values=finite_values, adjoint_value=adjoint_value
    )
end

results = map(compute_case, TABLE_V_CASES)
output_directory = joinpath(@__DIR__, "group_velocity_results")
mkpath(output_directory)
output_file = joinpath(output_directory, "table_v_critical_group_velocity_Re$(RES)_N$(N_CHEB).dat")

open(output_file, "w") do io
    println(io, "VARIABLES = \"a_s\" \"R\" \"beta\" \"table_alpha_r\" \"table_alpha_i\" \"table_omega_r\" \"table_omega_i\" \"table_cg_abs_fd\" \"table_cg_abs_adjoint\" \"saddle_alpha_r\" \"saddle_alpha_i\" \"saddle_omega_r\" \"saddle_omega_i\" \"cg_r_fd\" \"cg_i_fd\" \"cg_abs_fd\" \"cg_r_adjoint\" \"cg_i_adjoint\" \"cg_abs_adjoint\"")
    for result in results
        finite_value = result.finite_values[end]
        @printf(io, "%.1f %.8f %.8f %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
            result.a_s, result.radius, result.beta, real(result.table_alpha), imag(result.table_alpha),
            real(result.table_omega), imag(result.table_omega), abs(result.table_finite_value),
            abs(result.table_adjoint_value), real(result.alpha), imag(result.alpha),
            real(result.omega), imag(result.omega), real(finite_value), imag(finite_value),
            abs(finite_value), real(result.adjoint_value), imag(result.adjoint_value), abs(result.adjoint_value))
    end
end

summary_file = joinpath(output_directory, "table_v_critical_group_velocity_Re$(RES)_N$(N_CHEB)_summary.txt")
open(summary_file, "w") do io
    for result in results
        @printf(io, "a_s=%+.1f R=%.1f beta=%.4f\n", result.a_s, result.radius, result.beta)
        @printf(io, "  table alpha = %.14e %+.14ei\n", real(result.table_alpha), imag(result.table_alpha))
        @printf(io, "  omega at table alpha = %.14e %+.14ei\n", real(result.table_omega), imag(result.table_omega))
        @printf(io, "  table-alpha FD cg = %.14e %+.14ei, |cg|=%.14e\n",
            real(result.table_finite_value), imag(result.table_finite_value), abs(result.table_finite_value))
        @printf(io, "  table-alpha adjoint cg = %.14e %+.14ei, |cg|=%.14e\n",
            real(result.table_adjoint_value), imag(result.table_adjoint_value), abs(result.table_adjoint_value))
        @printf(io, "  alpha = %.14e %+.14ei\n", real(result.alpha), imag(result.alpha))
        @printf(io, "  omega = %.14e %+.14ei\n", real(result.omega), imag(result.omega))
        for (step, value) in zip(DERIVATIVE_STEPS, result.finite_values)
            @printf(io, "  FD h=%.1e: cg = %.14e %+.14ei, |cg|=%.14e\n",
                step, real(value), imag(value), abs(value))
        end
        @printf(io, "  adjoint: cg = %.14e %+.14ei, |cg|=%.14e\n\n",
            real(result.adjoint_value), imag(result.adjoint_value), abs(result.adjoint_value))
    end
end

println("Wrote ", output_file)
println("Wrote ", summary_file)
