using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "total_amplitude.jl"))

const SOURCE_CURVE = joinpath(
    @__DIR__, "c_6_validation_results", "c6_lower_neutral_A_N.dat",
)
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "fig21_forcing_scan_step25_results")
const TARGET_RADIUS = 570.0
const FORCING_END_RADIUS = 500.0
const FORCING_STEP = 25.0
const MAXIMUM_GROWTH_RADIUS = 372.988887508
const RECEPTIVITY_PEAK_RADIUS = 468.909531590
const MAXIMUM_GROWTH_CURVE = joinpath(
    @__DIR__, "c_6_validation_results", "c6_maximum_growth_Rg_A_N.dat",
)
const RECEPTIVITY_PEAK_CURVE = joinpath(
    @__DIR__, "c_6_validation_results", "c6_global_Cr_peak_A_N.dat",
)

function read_existing_branch(path)
    rows = NamedTuple[]
    for line in eachline(path)
        values = tryparse.(Float64, split(strip(line)))
        if length(values) == 8 && all(!isnothing, values)
            data = Float64[value for value in values]
            push!(rows, (
                R=data[1], N_from_lower=data[2],
                alpha=complex(data[5], data[6]),
            ))
        end
    end
    isempty(rows) && error("No numeric branch data found in $path")
    issorted(getfield.(rows, :R)) || error("Source radii are not sorted")
    return rows
end

function interpolate_branch(rows, radius)
    radii = getfield.(rows, :R)
    radius < radii[1] - 1.0e-10 && error("Radius lies below source curve")
    radius > radii[end] + 1.0e-10 && error("Radius lies above source curve")
    index = searchsortedfirst(radii, radius)
    if index <= length(rows) && abs(rows[index].R - radius) < 1.0e-10
        return rows[index]
    end
    lower = rows[index - 1]
    upper = rows[index]
    fraction = (radius - lower.R) / (upper.R - lower.R)
    return (
        R=radius,
        N_from_lower=lower.N_from_lower +
            fraction * (upper.N_from_lower - lower.N_from_lower),
        alpha=lower.alpha + fraction * (upper.alpha - lower.alpha),
    )
end

function forcing_radii(lower_radius)
    radii = collect(lower_radius:FORCING_STEP:FORCING_END_RADIUS)
    if abs(radii[end] - FORCING_END_RADIUS) > 100eps(FORCING_END_RADIUS)
        push!(radii, FORCING_END_RADIUS)
    end
    append!(radii, (MAXIMUM_GROWTH_RADIUS, RECEPTIVITY_PEAK_RADIUS))
    sort!(radii)
    return unique(radii)
end

function location_label(radius)
    if abs(radius - MAXIMUM_GROWTH_RADIUS) < 1.0e-8
        return "R_g maximum spatial growth"
    elseif abs(radius - RECEPTIVITY_PEAK_RADIUS) < 1.0e-8
        return "R_p maximum receptivity"
    end
    return @sprintf("R_f=%.6f", radius)
end

function special_curve_path(radius)
    if abs(radius - MAXIMUM_GROWTH_RADIUS) < 1.0e-8
        return MAXIMUM_GROWTH_CURVE
    elseif abs(radius - RECEPTIVITY_PEAK_RADIUS) < 1.0e-8
        return RECEPTIVITY_PEAK_CURVE
    end
    return nothing
end

function solve_forcing_points(radii, rows, cfg, F, G, H, D, D2)
    solutions = []
    previous_vector = nothing
    for (index, radius) in enumerate(radii)
        reference = interpolate_branch(rows, radius)
        @printf(
            "Direct mode %d/%d at R_f=%.12f, existing alpha=% .9f%+.9fi\n",
            index, length(radii), radius,
            real(reference.alpha), imag(reference.alpha),
        )
        solution = direct_solution(
            radius, reference.alpha, previous_vector,
            cfg, F, G, H, D, D2,
        )
        abs(solution.alpha - reference.alpha) < 1.0e-6 || error(
            "Solved mode at R_f=$radius does not match the existing branch",
        )
        push!(solutions, solution)
        previous_vector = solution.vector
    end
    return solutions
end

function propagate_with_existing_integral(
    forcing_solution, rows, cfg, F, G, D, D2, mass_matrix,
)
    forcing_radius = forcing_solution.R
    receptivity = initial_receptivity(
        forcing_solution, cfg, F, G, D, D2, mass_matrix,
    )

    special_path = special_curve_path(forcing_radius)
    if special_path === nothing
        forcing_reference = interpolate_branch(rows, forcing_radius)
        downstream = [forcing_reference]
        for row in rows
            if row.R > forcing_radius + 100eps(forcing_radius)
                push!(downstream, row)
            end
        end
    else
        downstream = read_existing_branch(special_path)
        forcing_reference = downstream[1]
        abs(forcing_reference.R - forcing_radius) < 1.0e-10 ||
            error("Special curve does not start at its specified forcing radius")
    end
    abs(downstream[end].R - TARGET_RADIUS) < 1.0e-10 ||
        error("Existing alpha_i curve does not end at R=570")

    radius = getfield.(downstream, :R)
    alpha = getfield.(downstream, :alpha)
    n_factor = getfield.(downstream, :N_from_lower) .-
        forcing_reference.N_from_lower
    n_factor[1] = 0.0
    gain = exp.(n_factor)
    amplitude = abs(receptivity.C_r) .* gain

    return (
        forcing_radius=forcing_radius, receptivity=receptivity,
        radius=radius, alpha=alpha, n_factor=n_factor,
        gain=gain, amplitude=amplitude,
        adjoint_error=abs(
            receptivity.alpha_adjoint - conj(forcing_solution.alpha),
        ),
    )
end

function write_curve_zone(io, result)
    @printf(
        io, "ZONE T=\"%s\", I=%d, F=POINT\n",
        location_label(result.forcing_radius), length(result.radius),
    )
    @printf(io, "AUXDATA forcing_radius=\"%.12e\"\n", result.forcing_radius)
    @printf(io, "AUXDATA Cr_initial_abs=\"%.12e\"\n", abs(result.receptivity.C_r))
    @printf(io, "AUXDATA adjoint_pairing_error=\"%.12e\"\n", result.adjoint_error)
    for index in eachindex(result.radius)
        alpha = result.alpha[index]
        @printf(
            io,
            "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
            result.radius[index], result.forcing_radius,
            real(alpha), imag(alpha), -imag(alpha),
            result.n_factor[index], result.gain[index],
            abs(result.receptivity.C_r), result.amplitude[index],
        )
    end
end

function write_endpoint_zone(io, results)
    println(io, "ZONE T=\"R=570 endpoint map\", I=$(length(results)), F=POINT")
    println(io, "AUXDATA zone_role=\"A(570;R_f) and N(570;R_f) versus forcing radius\"")
    for result in results
        alpha = result.alpha[end]
        @printf(
            io,
            "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
            result.radius[end], result.forcing_radius,
            real(alpha), imag(alpha), -imag(alpha),
            result.n_factor[end], result.gain[end],
            abs(result.receptivity.C_r), result.amplitude[end],
        )
    end
end

function write_outputs(results, cfg, lower_radius)
    mkpath(OUTPUT_DIRECTORY)
    curve_path = joinpath(OUTPUT_DIRECTORY, "fig21_N_A_forcing_scan_step25.dat")
    open(curve_path, "w") do io
        println(io, "TITLE=\"Fig. 21 forcing-location scan: N-factor and total amplitude\"")
        println(
            io,
            "VARIABLES=\"R\" \"R_f\" \"alpha_r\" \"alpha_i\" " *
            "\"growth_rate\" \"N_factor\" \"gain\" \"Cr_initial_abs\" \"A_abs\"",
        )
        println(io, "DATASETAUXDATA a_s=\"$(cfg.mass_flux)\"")
        println(io, "DATASETAUXDATA Re_h=\"$(cfg.Re_h)\"")
        println(io, "DATASETAUXDATA omega_bar=\"$(cfg.omega_bar)\"")
        println(io, "DATASETAUXDATA n=\"$(cfg.azimuthal_mode)\"")
        println(io, "DATASETAUXDATA c_squared=\"1\"")
        println(io, "DATASETAUXDATA N_cheb=\"$(cfg.N_cheb)\"")
        @printf(io, "DATASETAUXDATA lower_neutral_radius=\"%.12e\"\n", lower_radius)
        println(io, "DATASETAUXDATA forcing_radius_step=\"$(FORCING_STEP)\"")
        println(io, "DATASETAUXDATA forcing_radius_end=\"$(FORCING_END_RADIUS)\"")
        println(io, "DATASETAUXDATA R_g=\"$(MAXIMUM_GROWTH_RADIUS)\"")
        println(io, "DATASETAUXDATA R_p=\"$(RECEPTIVITY_PEAK_RADIUS)\"")
        println(io, "DATASETAUXDATA downstream_target_radius=\"$(TARGET_RADIUS)\"")
        println(io, "DATASETAUXDATA alpha_i_source=\"c6_lower_neutral_A_N.dat\"")
        println(io, "DATASETAUXDATA R_g_curve_source=\"c6_maximum_growth_Rg_A_N.dat\"")
        println(io, "DATASETAUXDATA R_p_curve_source=\"c6_global_Cr_peak_A_N.dat\"")
        for result in results
            write_curve_zone(io, result)
        end
        write_endpoint_zone(io, results)
    end

    summary_path = joinpath(OUTPUT_DIRECTORY, "fig21_N_A_forcing_scan_step25_summary.txt")
    open(summary_path, "w") do io
        println(io, "Fig. 21 forcing-location scan")
        println(io, "Parameters: a_s=0, Re_h=1000, omega_bar=0, n=30, c^2=1, N=99")
        println(io, "Existing alpha_i source = ", SOURCE_CURVE)
        println(io, "Lower neutral radius = ", lower_radius)
        println(io, "Forcing-radius step = ", FORCING_STEP)
        println(io, "Forcing-radius final point = ", FORCING_END_RADIUS)
        println(io, "Maximum-growth radius R_g = ", MAXIMUM_GROWTH_RADIUS)
        println(io, "Maximum-receptivity radius R_p = ", RECEPTIVITY_PEAK_RADIUS)
        println(io, "Downstream target radius = ", TARGET_RADIUS)
        println(io)
        println(io, "# label R_f Cr_abs N_at_570 gain_at_570 A_abs_at_570 points adjoint_error")
        for result in results
            @printf(
                io,
                "%s %.12e %.12e %.12e %.12e %.12e %d %.12e\n",
                replace(location_label(result.forcing_radius), " " => "_"),
                result.forcing_radius, abs(result.receptivity.C_r),
                result.n_factor[end], result.gain[end], result.amplitude[end],
                length(result.radius), result.adjoint_error,
            )
        end
        println(io)
        println(io, "Tecplot file = ", curve_path)
    end
    return curve_path, summary_path
end

function main()
    LinearAlgebra.BLAS.set_num_threads(1)
    cfg = Config(
        target_radius=TARGET_RADIUS, N_cheb=99,
        candidate_count=1, roughness_localization=0.5,
    )
    cfg.roughness_height = 1 / cfg.Re_h

    rows = read_existing_branch(SOURCE_CURVE)
    lower_radius = rows[1].R
    abs(rows[1].N_from_lower) < 1.0e-12 || error("Source N must start at zero")
    abs(rows[end].R - TARGET_RADIUS) < 1.0e-10 || error("Source must end at R=570")

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode,
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    radii = forcing_radii(lower_radius)
    println("Using existing alpha_i and N curve: ", SOURCE_CURVE)
    println("Only the direct/adjoint receptivity problem is solved at each R_f.")
    forcing_solutions = solve_forcing_points(radii, rows, cfg, F, G, H, D, D2)

    results = []
    for (index, solution) in enumerate(forcing_solutions)
        @printf(
            "Receptivity %d/%d at R_f=%.12f\n",
            index, length(forcing_solutions), solution.R,
        )
        push!(
            results,
            propagate_with_existing_integral(
                solution, rows, cfg, F, G, D, D2, mass_matrix,
            ),
        )
    end

    maximum(result.adjoint_error for result in results) < 1.0e-8 ||
        error("Adjoint eigenvalue pairing validation failed")
    all(abs(result.n_factor[1]) < 1.0e-14 for result in results) ||
        error("N(R_f;R_f) must be zero")
    all(
        abs(result.amplitude[1] - abs(result.receptivity.C_r)) < 1.0e-12
        for result in results
    ) || error("A(R_f;R_f) must equal |C_r(R_f)|")
    all(abs(result.radius[end] - TARGET_RADIUS) < 1.0e-10 for result in results) ||
        error("Every downstream curve must end at R=570")

    curve_path, summary_path = write_outputs(results, cfg, lower_radius)
    println("Tecplot data: ", curve_path)
    println("Summary: ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
