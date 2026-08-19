using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "q_bc_complex_diagnostics.jl"))

function paired_modes(solution, cfg, F, G, D, D2)
    A0_raw, A1_raw, A2_raw = CRC_STA.assemble_adjmat(
        solution.cof, D, D2, solution.beta, solution.omega, solution.R
    )
    A0, A1, A2 = CRC_STA.boudary_condition(
        A0_raw, A1_raw, A2_raw, cfg.N_cheb, cfg.mode
    )
    values, vectors = iar(
        PEP([A0, A1, A2]);
        σ = solution.alpha,
        neigs = 1,
        maxit = 500,
        tol = cfg.eig_tolerance,
    )
    direct, transpose_adjoint = normalize_pair(
        solution.vector, vectors[:, 1], cfg.N_cheb
    )
    hermitian_adjoint = conj.(transpose_adjoint)
    tilde_alpha = conj(values[1])
    return direct, hermitian_adjoint, tilde_alpha
end

function normalized_pairing(adjoint, operator_vector, mass_matrix)
    numerator = abs(adjoint' * mass_matrix * operator_vector)
    denominator = mass_norm(adjoint, mass_matrix) *
        mass_norm(operator_vector, mass_matrix)
    return numerator / denominator
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

    previous = nothing
    targets = Set([420.0, 460.0, 470.0, 480.0, 500.0])
    println("# R rho_11 rho_12 rho_21 rho_22 abs_B12 abs_B21")
    for radius in 420.0:1.0:500.0
        solutions = if previous === nothing
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
                    radius, previous[1].alpha, previous[1].vector,
                    cfg, F, G, H, D, D2
                ),
                direct_solution(
                    radius, previous[2].alpha, previous[2].vector,
                    cfg, F, G, H, D, D2
                ),
            )
        end
        if radius in targets
            x1, y1, t1 = paired_modes(solutions[1], cfg, F, G, D, D2)
            x2, y2, t2 = paired_modes(solutions[2], cfg, F, G, D, D2)
            xs = (x1, x2)
            ys = (y1, y2)
            ts = (t1, t2)
            rho = zeros(2, 2)
            raw = zeros(ComplexF64, 2, 2)
            for m in 1:2, n in 1:2
                operator = solutions[n].L1 +
                    (solutions[n].alpha + conj(ts[m])) * solutions[n].L2
                operator_vector = operator * xs[n]
                raw[m, n] = ys[m]' * mass_matrix * operator_vector
                rho[m, n] = normalized_pairing(
                    ys[m], operator_vector, mass_matrix
                )
            end
            @printf(
                "%.1f %.8e %.8e %.8e %.8e %.8e %.8e\n",
                radius, rho[1, 1], rho[1, 2], rho[2, 1], rho[2, 2],
                abs(raw[1, 2]), abs(raw[2, 1])
            )
        end
        previous = solutions
    end
end

main()
