using LinearAlgebra
using NonlinearEigenproblems
using Printf

include(joinpath(@__DIR__, "..", "total_amplitude.jl"))

function shape_diagnostic(solution, cfg, F, G, D, D2, mass_matrix, z)
    A0r, A1r, A2r = CRC_STA.assemble_adjmat(
        solution.cof, D, D2, solution.beta, solution.omega, solution.R
    )
    A0, A1, A2 = CRC_STA.boudary_condition(
        A0r, A1r, A2r, cfg.N_cheb, cfg.mode
    )
    values, vectors = iar(
        PEP([A0, A1, A2]); σ = solution.alpha, neigs = 1,
        maxit = 500, tol = cfg.eig_tolerance
    )
    direct, adjoint_t = normalize_pair(
        solution.vector, vectors[:, 1], cfg.N_cheb
    )
    adjoint_h = conj.(adjoint_t)
    adjoint_fields = full_mode(adjoint_h, cfg.N_cheb)
    u_adj = abs.(adjoint_fields[1])
    z_peak = z[argmax(u_adj)]

    # Stability_Cavity.jl receives G-1 as the rotating-frame azimuthal base
    # velocity, so use the same convention for physical layer diagnostics.
    G_rotating = G .- 1
    critical_measure = abs.(real(solution.alpha) .* vec(F) .+
        solution.beta .* vec(G_rotating) .- solution.omega)
    shear_measure = abs.(real(solution.alpha) .* vec(D * F) .+
        solution.beta .* vec(D * G_rotating))
    interior = 2:length(z)-1
    z_critical = z[interior[argmin(critical_measure[interior])]]
    z_shear = z[interior[argmin(shear_measure[interior])]]
    return (
        alpha = solution.alpha,
        z_peak = only(z_peak),
        z_critical = only(z_critical),
        z_shear = only(z_shear),
        distance_critical = abs(only(z_peak) - only(z_critical)),
        distance_shear = abs(only(z_peak) - only(z_shear)),
    )
end

cfg = Config(N_cheb = 99, candidate_count = 1)
cfg.roughness_height = 1 / cfg.Re_h
u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode)
D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

open(joinpath(@__DIR__, "..", "c7_validation_results", "c7_shape_diagnostics.dat"), "w") do io
    println(io, "# R branch alpha_r alpha_i z_adj_u_peak z_critical z_shear distance_to_critical distance_to_shear")
    for radius in (420.0, 440.0, 460.0, 470.0, 480.0, 500.0)
        for (label, shift) in (("I", 0.20 - 0.015im), ("II", 0.18 + 0.030im))
            solution = direct_solution(radius, shift, nothing, cfg, F, G, H, D, D2)
            diag = shape_diagnostic(solution, cfg, F, G, D, D2, mass_matrix, z)
            @printf(
                io, "%.6f %s %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                radius, label, real(diag.alpha), imag(diag.alpha),
                diag.z_peak, diag.z_critical, diag.z_shear,
                diag.distance_critical, diag.distance_shear
            )
            @printf(
                "R=%5.0f %s zpeak=%.5f zcrit=%.5f zshear=%.5f dcrit=%.5f dshear=%.5f\n",
                radius, label, diag.z_peak, diag.z_critical, diag.z_shear,
                diag.distance_critical, diag.distance_shear
            )
        end
    end
end
