include(joinpath(@__DIR__, "..", "zarf_spatial_growth_scan.jl"))

using Printf

LinearAlgebra.BLAS.set_num_threads(1)

cfg = Config(
    N_cheb=99,
    candidate_count=2,
    eig_tolerance=1e-12,
    max_iterations=700,
    newton_iterations=8,
    newton_tolerance=1e-9,
)

mass_flux = -0.4
radius = 500.0
beta = 0.2
temporal_path = joinpath(temporal_directory(cfg.zarf_root, mass_flux), "R=500.dat")
temporal_rows = parse_temporal_file(temporal_path)
temporal_rows = filter(row -> abs(row.beta - beta) < 1e-10, temporal_rows)
isempty(temporal_rows) && error("No temporal rows found at R=500, beta=0.2")
temporal_seed = temporal_rows[argmin(abs(row.omega_r) for row in temporal_rows)]
base_sample = merge(temporal_seed, (omega_r=0.0, branch_id=1))

u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, -1, mass_flux, 1)
D, D2, z = CRC_BF.Cheb(cfg.N_cheb, 1)
F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, 1)
cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, radius, cfg.N_cheb, D, D2, cfg.Re_h)

zero_roots = solve_spatial(base_sample, cof, D, D2, cfg)
sort!(zero_roots; by=root -> -imag(root.alpha))
length(zero_roots) >= 2 || error("Fewer than two admissible candidate roots at omega=0")
zero_roots = zero_roots[1:2]

function record_frequency(sample, evaluation, omega)
    merge(sample, (
        omega_r=omega,
        alpha_r=real(evaluation.alpha),
        alpha_i=imag(evaluation.alpha),
        growth_rate=evaluation.growth,
        residual=evaluation.residual,
        derivative=evaluation.derivative,
    ))
end

frequency_values = collect(range(-0.02, 0.02; length=9))
branch_results = Vector{Vector{NamedTuple}}()
newton_results = NamedTuple[]
for (branch_index, root) in enumerate(zero_roots)
    local_results = NamedTuple[]
    for omega in frequency_values
        evaluation = try
            evaluate_frequency(base_sample, omega, root.alpha, cof, D, D2, cfg)
        catch exception
            @warn "Frequency evaluation failed" branch_index omega exception
            nothing
        end
        evaluation === nothing || push!(local_results, record_frequency(base_sample, evaluation, omega))
    end
    sort!(local_results; by=row -> row.omega_r)
    for local_index in 2:(length(local_results) - 1)
        left = local_results[local_index - 1]
        seed = local_results[local_index]
        right = local_results[local_index + 1]
        seed.growth_rate >= left.growth_rate || continue
        seed.growth_rate >= right.growth_rate || continue
        peak = try
            newton_peak(base_sample, left, seed, right, cof, D, D2, cfg)
        catch exception
            @warn "Newton peak failed" branch_index omega=seed.omega_r exception
            nothing
        end
        peak === nothing || push!(newton_results, merge(base_sample, (
            branch_id=branch_index,
            omega_r=peak.omega,
            alpha_r=real(peak.alpha),
            alpha_i=imag(peak.alpha),
            growth_rate=peak.growth,
            residual=peak.residual,
        )))
    end
    push!(branch_results, local_results)
end

output_path = joinpath(
    @__DIR__, "..", "zarf_spatial_growth_results", "Ts_m0p4_N99_newton",
    "R500_beta020_two_root_newton.dat",
)
open(output_path, "w") do io
    println(io, "# a_s R beta temporal_seed_omega candidate_count")
    @printf(io, "# %.6f %.6f %.6f %.12e %d\n", mass_flux, radius, beta, temporal_seed.omega_r, cfg.candidate_count)
    println(io, "# zero-frequency roots: branch alpha_r alpha_i growth residual")
    for (branch_index, root) in enumerate(zero_roots)
        @printf(io, "# %d %.12e %.12e %.12e %.12e\n",
                branch_index, real(root.alpha), imag(root.alpha), -imag(root.alpha), root.residual)
    end
    println(io, "# scan: branch omega alpha_r alpha_i growth residual dgrowth_domega_imag")
    for (branch_index, rows) in enumerate(branch_results)
        for row in rows
            @printf(io, "scan %d %.12e %.12e %.12e %.12e %.12e %.12e\n",
                    branch_index, row.omega_r, row.alpha_r, row.alpha_i,
                    row.growth_rate, row.residual, imag(row.derivative))
        end
    end
    println(io, "# Newton peaks: branch omega alpha_r alpha_i growth residual")
    for row in newton_results
        @printf(io, "newton %d %.12e %.12e %.12e %.12e %.12e\n",
                row.branch_id, row.omega_r, row.alpha_r, row.alpha_i,
                row.growth_rate, row.residual)
    end
end

println("output: $output_path")
println(@sprintf("temporal seed omega=%.12e alpha_t=%.12e omega_i=%.12e",
                 temporal_seed.omega_r, temporal_seed.alpha_t, temporal_seed.omega_i))
for (branch_index, root) in enumerate(zero_roots)
    println(@sprintf("zero root %d: alpha=%.12e%+.12ei growth=%.12e residual=%.6e",
                     branch_index, real(root.alpha), imag(root.alpha),
                     -imag(root.alpha), root.residual))
end
for row in newton_results
    println(@sprintf("Newton branch %d: omega=%.12e growth=%.12e alpha=%.12e%+.12ei residual=%.6e",
                     row.branch_id, row.omega_r, row.growth_rate,
                     row.alpha_r, row.alpha_i, row.residual))
end
