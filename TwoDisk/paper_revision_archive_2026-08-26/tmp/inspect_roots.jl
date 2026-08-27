using LinearAlgebra
using NonlinearEigenproblems

include(joinpath(@__DIR__, "..", "total_amplitude.jl"))

cfg = Config(N_cheb = 99, candidate_count = 1)
u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode)
D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)

for (radius, shift) in [
    (360.0, 0.318 - 0.031im),
    (400.0, 0.255 - 0.029im),
    (420.0, 0.230 - 0.024im),
    (440.0, 0.208 - 0.017im),
    (460.0, 0.188 - 0.010im),
    (480.0, 0.169 - 0.006im),
    (500.0, 0.153 - 0.005im),
]
    beta = cfg.azimuthal_mode / radius
    omega = cfg.omega_bar / radius
    cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, radius, cfg.N_cheb, D, D2, cfg.Re_h)
    L0r, L1r, L2r = CRC_STA.assemble_mat(cof, D, D2, beta, omega, radius)
    L0, L1, L2 = CRC_STA.boudary_condition(L0r, L1r, L2r, cfg.N_cheb, cfg.mode)
    values, vectors = iar(PEP([L0, L1, L2]); σ = shift, neigs = 8, maxit = 750, tol = 1e-13)
    order = sortperm(abs.(values .- shift))
    println("R = ", radius, ", shift = ", shift)
    for j in order
        println("  ", values[j], " distance=", abs(values[j] - shift))
    end
end
