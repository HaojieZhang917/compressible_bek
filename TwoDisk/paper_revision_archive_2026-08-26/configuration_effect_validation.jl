using LinearAlgebra
using DelimitedFiles
using Printf
using NonlinearEigenproblems

include(joinpath(@__DIR__, "tmp", "test_bek2_configuration.jl"))
include(joinpath(@__DIR__, "total_amplitude.jl"))

const OUTPUT_DIR = joinpath(@__DIR__, "configuration_effect_validation_results")
const FOURIER_GAUSSIAN_FACTOR = sqrt(2pi)

function cavity_context(N::Int)
    cfg = Config(
        target_radius=600.0,
        Re_h=1000,
        mass_flux=0.0,
        mode=1,
        Ro=-1.0,
        azimuthal_mode=30,
        omega_bar=0.0,
        N_cheb=N,
        candidate_count=1,
        eig_tolerance=1e-13,
        roughness_height=1 / 1000,
        roughness_localization=0.5,
    )
    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    W = cheb_mass_matrix(cfg.N_cheb, cfg.mode)
    return (; cfg, F, G, H, D, D2, z, W)
end

function cavity_point(radius, shift, previous_vector, ctx)
    solution = direct_solution(
        radius, shift, previous_vector, ctx.cfg,
        ctx.F, ctx.G, ctx.H, ctx.D, ctx.D2
    )
    A0_raw, A1_raw, A2_raw = CRC_STA.assemble_adjmat(
        solution.cof, ctx.D, ctx.D2, solution.beta, solution.omega, solution.R
    )
    A0, A1, A2 = CRC_STA.boudary_condition(
        A0_raw, A1_raw, A2_raw, ctx.cfg.N_cheb, ctx.cfg.mode
    )
    alpha_t, y_t = iar(
        PEP([A0, A1, A2]); σ=solution.alpha, neigs=1,
        maxit=500, tol=ctx.cfg.eig_tolerance
    )
    x, yt = normalize_pair(solution.vector, y_t[:, 1], ctx.cfg.N_cheb)
    y = conj.(yt)
    alpha_a = conj(alpha_t[1])
    Q = y' * ctx.W *
        (solution.L1 + (solution.alpha + conj(alpha_a)) * solution.L2) * x
    direct_fields = full_mode(x, ctx.cfg.N_cheb)
    adjoint_fields = full_mode(y, ctx.cfg.N_cheb)
    hhat = ctx.cfg.roughness_height *
           exp(-solution.alpha^2 / (4ctx.cfg.roughness_localization)) *
           sqrt(pi / ctx.cfg.roughness_localization)
    u_wall = -(ctx.D * ctx.F)[1] * hhat
    v_wall = -(ctx.D * ctx.G)[1] * hhat
    BC_u = conj((ctx.D * adjoint_fields[1])[1]) * u_wall /
           (solution.R * sqrt(ctx.cfg.Re_h))
    BC_v = conj((ctx.D * adjoint_fields[2])[1]) * v_wall /
           (solution.R * sqrt(ctx.cfg.Re_h))
    BC = BC_u + BC_v
    Cr = abs(-im * BC / Q)
    P = solution.cof
    L0_raw, _, _ = CRC_STA.assemble_mat(
        P, ctx.D, ctx.D2, solution.beta, solution.omega, solution.R
    )
    L0, _, _ = CRC_STA.boudary_condition(
        L0_raw, solution.L1, solution.L2, ctx.cfg.N_cheb, ctx.cfg.mode
    )
    direct_residual = norm(
        (L0 + solution.alpha * solution.L1 + solution.alpha^2 * solution.L2) * x
    ) / ((norm(L0) + abs(solution.alpha) * norm(solution.L1) +
           abs(solution.alpha)^2 * norm(solution.L2)) * norm(x))
    adjoint_residual = norm(
        (A0 + alpha_t[1] * A1 + alpha_t[1]^2 * A2) * yt
    ) / ((norm(A0) + abs(alpha_t[1]) * norm(A1) +
           abs(alpha_t[1])^2 * norm(A2)) * norm(yt))
    return (; R=radius, alpha=solution.alpha, alpha_a, Cr, Q, BC, BC_u, BC_v,
            vector=solution.vector, direct_fields, adjoint_fields,
            pairing_error=abs(alpha_t[1] - solution.alpha),
            direct_residual, adjoint_residual)
end

function scan_configuration(; N=99, radii=collect(300.0:10.0:600.0))
    ctx = cavity_context(N)
    cavity_results = []
    single_results = []

    cavity_shift = 0.43 - 0.01im
    cavity_vector = nothing
    single_shift = 0.38 - 0.01im
    for R in radii
        cavity = cavity_point(R, cavity_shift, cavity_vector, ctx)
        push!(cavity_results, cavity)
        cavity_shift = cavity.alpha
        cavity_vector = cavity.vector

        single = solve_point(R=R, n=30.0, N=N, shift=single_shift)
        push!(single_results, single)
        single_shift = single.alpha

        @printf("R=%6.1f  cavity alpha=% .7f%+.7fi Cr=% .7f  single alpha=% .7f%+.7fi Cr=% .7f\n",
                R, real(cavity.alpha), imag(cavity.alpha), cavity.Cr,
                real(single.alpha), imag(single.alpha), single.Cr_thomas)
    end
    return ctx, cavity_results, single_results
end

function write_curve(radii, cavity, single)
    path = joinpath(OUTPUT_DIR, "configuration_comparison_N99.dat")
    open(path, "w") do io
        println(io, "TITLE=\"Cavity and isolated-disk Type-I comparison\"")
        println(io, "VARIABLES=\"R\",\"alpha_r_cavity\",\"alpha_i_cavity\",\"growth_cavity\",\"Cr_cavity_native\",\"Cr_cavity_unitary\",\"alpha_r_single\",\"alpha_i_single\",\"growth_single\",\"Cr_single_native_Thomas\",\"Cr_single_unnormalized\",\"pair_error_cavity\",\"pair_error_single\"")
        println(io, "ZONE T=\"n=30, omega_bar=0, N=99\", I=$(length(radii)), F=POINT")
        for j in eachindex(radii)
            c, s = cavity[j], single[j]
            @printf(io, "%.10f %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
                    radii[j], real(c.alpha), imag(c.alpha), -imag(c.alpha),
                    c.Cr, c.Cr / FOURIER_GAUSSIAN_FACTOR,
                    real(s.alpha), imag(s.alpha), -imag(s.alpha),
                    s.Cr_thomas, s.Cr, c.pairing_error, s.pairing_error)
        end
    end
    return path
end

function write_profiles(ctx, cavity, single)
    path = joinpath(OUTPUT_DIR, "baseflow_modes_R470_N99.dat")
    open(path, "w") do io
        println(io, "TITLE=\"Base flow and Type-I mode at R=470\"")
        println(io, "VARIABLES=\"eta\",\"F\",\"G\",\"H\",\"abs_u\",\"abs_v\",\"abs_w\",\"abs_ua\",\"abs_va\",\"abs_wa\"")
        println(io, "ZONE T=\"rotor-stator cavity\", I=$(length(ctx.z)), F=POINT")
        for j in eachindex(ctx.z)
            @printf(io, "%.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
                    ctx.z[j] * sqrt(ctx.cfg.Re_h), ctx.F[j], ctx.G[j], ctx.H[j],
                    abs(cavity.direct_fields[1][j]), abs(cavity.direct_fields[2][j]),
                    abs(cavity.direct_fields[3][j]), abs(cavity.adjoint_fields[1][j]),
                    abs(cavity.adjoint_fields[2][j]), abs(cavity.adjoint_fields[3][j]))
        end
        eta = single.grid.z
        finite_ids = findall(isfinite, eta)
        println(io, "ZONE T=\"isolated disk\", I=$(length(finite_ids)), F=POINT")
        for j in finite_ids
            @printf(io, "%.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e\n",
                    eta[j], single.baseflow.F[j], single.baseflow.G[j],
                    single.baseflow.H[j], abs(single.direct_fields[1][j]),
                    abs(single.direct_fields[2][j]), abs(single.direct_fields[3][j]),
                    abs(single.adjoint_fields[1][j]), abs(single.adjoint_fields[2][j]),
                    abs(single.adjoint_fields[3][j]))
        end
    end
    return path
end

function nearest_index(radii, target)
    return argmin(abs.(radii .- target))
end

function main()
    mkpath(OUTPUT_DIR)
    radii = collect(300.0:10.0:600.0)
    ctx, cavity, single = scan_configuration(N=99, radii=radii)
    curve_path = write_curve(radii, cavity, single)
    j470 = nearest_index(radii, 470.0)
    profile_path = write_profiles(ctx, cavity[j470], single[j470])

    grid_results = Dict{Int, Any}()
    for N in (79, 99, 129)
        grid_results[N] = solve_point(
            R=470.0, n=30.0, N=N, shift=0.28 - 0.06im
        )
    end

    c = cavity[j470]
    s = single[j470]
    summary_path = joinpath(OUTPUT_DIR, "configuration_validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "Configuration-effect diagnostic validation")
        println(io, "Parameters: R=470, n=30, omega_bar=0, a_s=0, N=99")
        println(io)
        println(io, "Cavity alpha = ", repr(c.alpha))
        println(io, "Cavity |Cr| (current manuscript convention) = ", c.Cr)
        println(io, "Cavity pairing error = ", c.pairing_error)
        println(io, "Cavity direct residual = ", c.direct_residual)
        println(io, "Cavity adjoint residual = ", c.adjoint_residual)
        println(io)
        println(io, "Single-disk alpha = ", repr(s.alpha))
        println(io, "Single-disk |Cr| (Thomas/Fig. Cr_ref convention) = ", s.Cr_thomas)
        println(io, "Single-disk |Cr| (same unnormalized Gaussian transform as cavity code) = ", s.Cr)
        println(io, "Single-disk |Q| = ", abs(s.Q))
        println(io, "Single-disk |K| = ", abs(s.K))
        println(io, "Relative |Q|-|K| difference = ", abs(abs(s.Q)-abs(s.K))/abs(s.K))
        println(io, "Single-disk pairing error = ", s.pairing_error)
        println(io, "Single-disk direct residual = ", s.direct_residual)
        println(io, "Single-disk adjoint residual = ", s.adjoint_residual)
        println(io)
        println(io, "Growth-rate ratio (-alpha_i single)/(-alpha_i cavity) = ",
                (-imag(s.alpha))/(-imag(c.alpha)))
        println(io, "Native plotted Cr ratio cavity/single = ", c.Cr/s.Cr_thomas)
        println(io, "Same-Fourier-convention Cr ratio cavity/single = ", c.Cr/s.Cr)
        println(io)
        println(io, "Single-disk grid convergence at R=470")
        for N in (79, 99, 129)
            r = grid_results[N]
            println(io, "N=", N, " alpha=", repr(r.alpha),
                    " Cr_Thomas=", r.Cr_thomas,
                    " pairing_error=", r.pairing_error)
        end
        println(io)
        println(io, "Malik/Thomas stationary neutral benchmark")
        neutral = solve_point(R=285.36, n=0.07759*285.36, N=99,
                              shift=0.38482 + 0im)
        println(io, "reference alpha approximately 0.38482+0i")
        println(io, "computed alpha = ", repr(neutral.alpha))
        println(io, "curve_file = ", curve_path)
        println(io, "profile_file = ", profile_path)
        println(io)
        println(io, "Interpretation warning: the current cavity and single-disk plotting paths use Gaussian Fourier amplitudes that differ by sqrt(2*pi). The raw BC and Q are also formulation-dependent. A physical configuration mechanism must therefore be based on a common forcing convention and scale-invariant modal diagnostics, not on the native Cr ratio alone.")
    end
    println("Wrote ", curve_path)
    println("Wrote ", profile_path)
    println("Wrote ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
