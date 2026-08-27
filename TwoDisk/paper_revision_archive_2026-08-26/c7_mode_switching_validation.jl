using LinearAlgebra
using DelimitedFiles
using Printf

include("total_amplitude.jl")

function mode_diagnostics(solution, cfg, F, G, D, D2, mass_matrix)
    result = initial_receptivity(solution, cfg, F, G, D, D2, mass_matrix)
    return (
        C_r = abs(result.C_r),
        Q = abs(result.Q),
        BC = abs(result.boundary_projection),
        adjoint_error = abs(result.alpha_adjoint - conj(solution.alpha)),
    )
end

function scan_pair(radius, cfg, F, G, H, D, D2, mass_matrix)
    # The two shifts deliberately target the two roots visible in the switching
    # interval: the unstable Type-I continuation and the stable Type-II root.
    type_i_shift = 0.20 - 0.015im
    type_ii_shift = 0.18 + 0.030im
    type_i = direct_solution(
        radius, type_i_shift, nothing, cfg, F, G, H, D, D2
    )
    type_ii = direct_solution(
        radius, type_ii_shift, nothing, cfg, F, G, H, D, D2
    )
    diag_i = mode_diagnostics(type_i, cfg, F, G, D, D2, mass_matrix)
    diag_ii = mode_diagnostics(type_ii, cfg, F, G, D, D2, mass_matrix)
    return type_i, diag_i, type_ii, diag_ii
end

function main()
    cfg = Config(
        N_cheb = 99,
        candidate_count = 1,
        roughness_localization = 0.5,
    )
    cfg.roughness_height = 1 / cfg.Re_h
    output_directory = abspath("c7_validation_results")
    mkpath(output_directory)

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    radii = collect(420.0:2.0:500.0)
    rows = []
    for radius in radii
        type_i, diag_i, type_ii, diag_ii = scan_pair(
            radius, cfg, F, G, H, D, D2, mass_matrix
        )
        push!(rows, (
            R = radius,
            alpha_i_mode = type_i.alpha,
            Cr_i = diag_i.C_r,
            Q_i = diag_i.Q,
            BC_i = diag_i.BC,
            adjerr_i = diag_i.adjoint_error,
            alpha_ii_mode = type_ii.alpha,
            Cr_ii = diag_ii.C_r,
            Q_ii = diag_ii.Q,
            BC_ii = diag_ii.BC,
            adjerr_ii = diag_ii.adjoint_error,
        ))
        @printf(
            "R=%6.1f I=% .8f%+.8fi CrI=%.6f II=% .8f%+.8fi CrII=%.6f dRe=%.3e\n",
            radius, real(type_i.alpha), imag(type_i.alpha), diag_i.C_r,
            real(type_ii.alpha), imag(type_ii.alpha), diag_ii.C_r,
            abs(real(type_i.alpha) - real(type_ii.alpha))
        )
    end

    differences_real = [
        abs(real(row.alpha_i_mode) - real(row.alpha_ii_mode)) for row in rows
    ]
    differences_full = [
        abs(row.alpha_i_mode - row.alpha_ii_mode) for row in rows
    ]
    closest_index = argmin(differences_real)
    closest = rows[closest_index]
    cr_cross_index = argmin([abs(row.Cr_i - row.Cr_ii) for row in rows])
    cr_cross = rows[cr_cross_index]

    output_path = joinpath(output_directory, "c7_two_branches.dat")
    open(output_path, "w") do io
        println(io, "# R alphaI_r alphaI_i growthI CrI QI BCI alphaII_r alphaII_i growthII CrII QII BCII delta_alpha_r delta_alpha_i delta_alpha_abs adjerrI adjerrII")
        for row in rows
            @printf(
                io, "%.12g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                row.R, real(row.alpha_i_mode), imag(row.alpha_i_mode),
                -imag(row.alpha_i_mode), row.Cr_i, row.Q_i, row.BC_i,
                real(row.alpha_ii_mode), imag(row.alpha_ii_mode),
                -imag(row.alpha_ii_mode), row.Cr_ii, row.Q_ii, row.BC_ii,
                abs(real(row.alpha_i_mode) - real(row.alpha_ii_mode)),
                abs(imag(row.alpha_i_mode) - imag(row.alpha_ii_mode)),
                abs(row.alpha_i_mode - row.alpha_ii_mode),
                row.adjerr_i, row.adjerr_ii
            )
        end
    end

    summary_path = joinpath(output_directory, "c7_validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "C7 two-branch diagnostic summary")
        println(io, "N_cheb = ", cfg.N_cheb)
        println(io, "R_min_delta_alpha_r = ", closest.R)
        println(io, "alpha_I = ", repr(closest.alpha_i_mode))
        println(io, "alpha_II = ", repr(closest.alpha_ii_mode))
        println(io, "delta_alpha_r = ", differences_real[closest_index])
        println(io, "delta_alpha_i = ", abs(imag(closest.alpha_i_mode) - imag(closest.alpha_ii_mode)))
        println(io, "delta_alpha_abs = ", differences_full[closest_index])
        println(io, "Cr_I_at_closest = ", closest.Cr_i)
        println(io, "Cr_II_at_closest = ", closest.Cr_ii)
        println(io, "Q_I_at_closest = ", closest.Q_i)
        println(io, "Q_II_at_closest = ", closest.Q_ii)
        println(io, "BC_I_at_closest = ", closest.BC_i)
        println(io, "BC_II_at_closest = ", closest.BC_ii)
        println(io, "R_min_Cr_difference = ", cr_cross.R)
        println(io, "Cr_I_at_cross = ", cr_cross.Cr_i)
        println(io, "Cr_II_at_cross = ", cr_cross.Cr_ii)
        println(io, "data_file = ", output_path)
    end

    println("\nClosest real parts at R=", closest.R)
    println("  alpha_I  = ", closest.alpha_i_mode)
    println("  alpha_II = ", closest.alpha_ii_mode)
    println("  delta alpha_r = ", differences_real[closest_index])
    println("  delta alpha_i = ", abs(imag(closest.alpha_i_mode) - imag(closest.alpha_ii_mode)))
    println("  full complex separation = ", differences_full[closest_index])
    println("Closest Cr values at R=", cr_cross.R)
    println("Results written to ", output_directory)
end

main()
