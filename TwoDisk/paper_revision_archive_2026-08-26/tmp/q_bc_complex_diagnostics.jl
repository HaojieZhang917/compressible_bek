using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "total_amplitude.jl"))

function mass_norm(vector, mass_matrix)
    value = real(vector' * mass_matrix * vector)
    return sqrt(max(value, 0.0))
end

function detailed_receptivity(solution, cfg, F, G, D, D2, mass_matrix)
    A0_raw, A1_raw, A2_raw = CRC_STA.assemble_adjmat(
        solution.cof, D, D2, solution.beta, solution.omega, solution.R
    )
    A0, A1, A2 = CRC_STA.boudary_condition(
        A0_raw, A1_raw, A2_raw, cfg.N_cheb, cfg.mode
    )
    adjoint_values, adjoint_vectors = iar(
        PEP([A0, A1, A2]);
        σ = solution.alpha,
        neigs = 1,
        maxit = 500,
        tol = cfg.eig_tolerance,
    )

    direct, adjoint_transpose = normalize_pair(
        solution.vector, adjoint_vectors[:, 1], cfg.N_cheb
    )
    alpha_adjoint = conj(adjoint_values[1])
    adjoint = conj.(adjoint_transpose)
    polynomial_derivative = (
        solution.L1 + (solution.alpha + conj(alpha_adjoint)) * solution.L2
    )
    derivative_vector = polynomial_derivative * direct
    q_linear = adjoint' * mass_matrix * solution.L1 * direct
    q_quadratic = adjoint' * mass_matrix *
        ((solution.alpha + conj(alpha_adjoint)) * solution.L2) * direct
    q_total = q_linear + q_quadratic

    adjoint_fields = full_mode(adjoint, cfg.N_cheb)
    roughness_transform = cfg.roughness_height *
        exp(-solution.alpha^2 / (4cfg.roughness_localization)) *
        sqrt(pi / cfg.roughness_localization)
    u_wall = -(D * F)[1] * roughness_transform
    v_wall = -(D * G)[1] * roughness_transform
    bc_u = conj((D * adjoint_fields[1])[1]) * u_wall /
        (solution.R * sqrt(cfg.Re_h))
    bc_v = conj((D * adjoint_fields[2])[1]) * v_wall /
        (solution.R * sqrt(cfg.Re_h))
    bc_total = bc_u + bc_v
    cr = -im * bc_total / q_total

    norm_adjoint = mass_norm(adjoint, mass_matrix)
    norm_direct = mass_norm(direct, mass_matrix)
    norm_derivative = mass_norm(derivative_vector, mass_matrix)
    q_angle = abs(q_total) / max(norm_adjoint * norm_derivative, eps())
    q_condition = norm_adjoint * norm_direct / max(abs(q_total), eps())

    return (
        alpha_adjoint = alpha_adjoint,
        q_linear = q_linear,
        q_quadratic = q_quadratic,
        q_total = q_total,
        bc_u = bc_u,
        bc_v = bc_v,
        bc_total = bc_total,
        cr = cr,
        norm_direct = norm_direct,
        norm_adjoint = norm_adjoint,
        norm_derivative = norm_derivative,
        q_angle = q_angle,
        q_condition = q_condition,
    )
end

function write_complex(io, value)
    @printf(io, " %.16e %.16e", real(value), imag(value))
end

function main()
    cfg = Config(
        N_cheb = 99,
        candidate_count = 3,
        roughness_localization = 0.5,
    )
    cfg.roughness_height = 1 / cfg.Re_h
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    output_path = joinpath(@__DIR__, "q_bc_complex_diagnostics.dat")
    open(output_path, "w") do io
        println(
            io,
            "# R branch alpha_r alpha_i " *
            "Q1_r Q1_i Q2_r Q2_i Q_r Q_i " *
            "BCu_r BCu_i BCv_r BCv_i BC_r BC_i Cr_r Cr_i " *
            "norm_x norm_y norm_Palpha_x q_angle q_condition adjoint_error"
        )
        previous_solutions = nothing
        for radius in 420.0:1.0:500.0
            solutions = if previous_solutions === nothing
                (
                    direct_solution(
                        radius, 0.23 - 0.023im, nothing,
                        cfg, F, G, H, D, D2
                    ),
                    direct_solution(
                        radius, 0.189 + 0.034im, nothing,
                        cfg, F, G, H, D, D2
                    ),
                )
            else
                (
                    direct_solution(
                        radius, previous_solutions[1].alpha,
                        previous_solutions[1].vector,
                        cfg, F, G, H, D, D2
                    ),
                    direct_solution(
                        radius, previous_solutions[2].alpha,
                        previous_solutions[2].vector,
                        cfg, F, G, H, D, D2
                    ),
                )
            end
            for (branch_index, solution) in enumerate(solutions)
                result = detailed_receptivity(
                    solution, cfg, F, G, D, D2, mass_matrix
                )
                @printf(
                    io, "%.8f %d %.16e %.16e",
                    radius, branch_index, real(solution.alpha), imag(solution.alpha)
                )
                write_complex(io, result.q_linear)
                write_complex(io, result.q_quadratic)
                write_complex(io, result.q_total)
                write_complex(io, result.bc_u)
                write_complex(io, result.bc_v)
                write_complex(io, result.bc_total)
                write_complex(io, result.cr)
                @printf(
                    io, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
                    result.norm_direct,
                    result.norm_adjoint,
                    result.norm_derivative,
                    result.q_angle,
                    result.q_condition,
                    abs(result.alpha_adjoint - conj(solution.alpha)),
                )
            end
            previous_solutions = solutions
            @printf("R = %.1f complete\n", radius)
        end
    end
    println("Wrote ", output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
