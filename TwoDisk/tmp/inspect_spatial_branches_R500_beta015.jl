include(joinpath(@__DIR__, "..", "zarf_spatial_growth_scan.jl"))

LinearAlgebra.BLAS.set_num_threads(1)
cfg = Config(N_cheb=99, candidate_count=12, eig_tolerance=1e-12, max_iterations=700)
radius = 500.0
beta = 0.15
omega = 0.0
mass_flux = 0.0

u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, -1, mass_flux, 1)
D, D2, z = CRC_BF.Cheb(cfg.N_cheb, 1)
F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, 1)
cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, radius, cfg.N_cheb, D, D2, cfg.Re_h)
L0, L1, L2, _ = spatial_matrices(cof, D, D2, beta, omega, radius, cfg)
problem = PEP([L0, L1, L2])

roots = NamedTuple[]
for alpha_shift in 0.05:0.05:1.0
    shift = complex(alpha_shift, -0.08)
    values, vectors = iar(
        problem;
        σ=shift,
        neigs=cfg.candidate_count,
        maxit=cfg.max_iterations,
        tol=cfg.eig_tolerance,
    )
    for index in eachindex(values)
        alpha = values[index]
        isfinite(real(alpha)) && isfinite(imag(alpha)) || continue
        real(alpha) > 0 || continue
        abs(imag(alpha)) < 0.5 || continue
        vector = vectors[:, index]
        residual = norm((L0 + alpha * L1 + alpha^2 * L2) * vector) / norm(vector)
        any(abs(alpha - row.alpha) < 1e-7 for row in roots) && continue
        push!(roots, (alpha=alpha, growth=-imag(alpha), residual=residual, shift=shift))
    end
end

sort!(roots; by=row -> -row.growth)
println("# alpha_r alpha_i growth residual shift_r shift_i")
for row in roots
    @printf(
        "%.12e %.12e %.12e %.12e %.6f %.6f\n",
        real(row.alpha), imag(row.alpha), row.growth, row.residual,
        real(row.shift), imag(row.shift),
    )
end
