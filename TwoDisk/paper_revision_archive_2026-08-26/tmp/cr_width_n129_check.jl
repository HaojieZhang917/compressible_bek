using Printf

include(joinpath(@__DIR__, "..", "total_amplitude.jl"))

for order in (99, 129)
    cfg = Config(N_cheb = order, candidate_count = 1)
    cfg.roughness_height = 1 / cfg.Re_h
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)
    solution = direct_solution(
        470.0, 0.178 - 0.0074im, nothing, cfg, F, G, H, D, D2
    )
    for localization in (0.5, 0.125)
        cfg.roughness_localization = localization
        result = initial_receptivity(
            solution, cfg, F, G, D, D2, mass_matrix
        )
        @printf(
            "N=%d ls=%.3f alpha=% .12f%+.12fi Cr=%.12f\n",
            order, localization, real(solution.alpha), imag(solution.alpha),
            abs(result.C_r)
        )
    end
end
