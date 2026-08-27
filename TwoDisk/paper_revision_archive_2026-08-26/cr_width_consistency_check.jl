using Printf

include("total_amplitude.jl")

function receptivity_at(radius, shift, localization, F, G, H, D, D2, mass_matrix)
    cfg = Config(
        N_cheb = 99,
        candidate_count = 1,
        roughness_localization = localization,
    )
    cfg.roughness_height = 1 / cfg.Re_h
    solution = direct_solution(radius, shift, nothing, cfg, F, G, H, D, D2)
    result = initial_receptivity(solution, cfg, F, G, D, D2, mass_matrix)
    return solution, result
end

cfg_base = Config(N_cheb = 99)
u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
    cfg_base.Re_h, cfg_base.Ro, cfg_base.mass_flux, cfg_base.mode
)
D, D2, z = CRC_BF.Cheb(cfg_base.N_cheb, cfg_base.mode)
F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg_base.N_cheb, cfg_base.mode)
mass_matrix = cheb_mass_matrix(cfg_base.N_cheb, cfg_base.mode)

points = [
    (287.821608895, 0.449677971 - 0.0im, "lower_neutral"),
    (360.0, 0.317205897 - 0.031420707im, "R360"),
    (372.988887508, 0.295803 - 0.032262im, "maximum_growth"),
    (468.909531590, 0.179054 - 0.007661im, "published_peak_region"),
    (470.0, 0.177999 - 0.007407im, "R470"),
]

output_dir = abspath("cr_width_consistency_results")
mkpath(output_dir)
output_path = joinpath(output_dir, "cr_width_comparison.dat")

open(output_path, "w") do io
    println(io, "# label R alpha_r alpha_i Cr_ls_0p5 Cr_ls_0p125 ratio new_expected_from_fourier_ratio")
    for (radius, shift, label) in points
        solution_old, old_result = receptivity_at(
            radius, shift, 0.5, F, G, H, D, D2, mass_matrix
        )
        solution_new, new_result = receptivity_at(
            radius, solution_old.alpha, 0.125, F, G, H, D, D2, mass_matrix
        )
        cr_old = abs(old_result.C_r)
        cr_new = abs(new_result.C_r)
        fourier_ratio = 2 * exp(-1.5 * real(solution_old.alpha^2))
        @printf(
            io, "%s %.12g %.16g %.16g %.16g %.16g %.16g %.16g\n",
            label, radius, real(solution_old.alpha), imag(solution_old.alpha),
            cr_old, cr_new, cr_new / cr_old, cr_old * fourier_ratio
        )
        @printf(
            "%-24s R=%9.4f alpha=% .9f%+.9fi Cr(ls=.5)=%.9f Cr(ls=.125)=%.9f ratio=%.6f\n",
            label, radius, real(solution_old.alpha), imag(solution_old.alpha),
            cr_old, cr_new, cr_new / cr_old
        )
    end
end

println("comparison file = ", output_path)
