using Printf

include("total_amplitude.jl")

const C6_OUTPUT_DIRECTORY = "c_6_validation_results"
const C6_TARGET_RADIUS = 570.0
const C6_MAXIMUM_GROWTH_RADIUS = 372.988887508
const C6_GLOBAL_PEAK_RADIUS = 468.90953159
const C6_RADIAL_STEP = 1.0

function tracked_curve(seed, cfg, F, G, H, D, D2)
    radii = radial_points(seed.R, cfg.target_radius, cfg.integration_step)
    solutions = [seed]
    for radius in radii[2:end]
        previous = solutions[end]
        push!(
            solutions,
            direct_solution(
                radius, previous.alpha, previous.vector,
                cfg, F, G, H, D, D2
            )
        )
    end
    return solutions
end

function solution_on_reference(radius, reference_branch, cfg, F, G, H, D, D2)
    reference_radii = getfield.(reference_branch, :R)
    reference = reference_branch[argmin(abs.(reference_radii .- radius))]
    return direct_solution(
        radius, reference.alpha, reference.vector,
        cfg, F, G, H, D, D2,
    )
end

function downstream_segment(seed, reference_branch, target_radius)
    solutions = [seed]
    for solution in reference_branch
        if solution.R > seed.R + 100eps(seed.R) &&
           solution.R <= target_radius + 100eps(target_radius)
            push!(solutions, solution)
        end
    end
    if abs(solutions[end].R - target_radius) > 100eps(target_radius)
        error("The reference branch does not end at the target radius.")
    end
    return solutions
end

function propagate_case(
    label, seed, reference_branch, cfg, F, G, D, D2, mass_matrix
)
    receptivity = initial_receptivity(
        seed, cfg, F, G, D, D2, mass_matrix
    )
    solutions = downstream_segment(seed, reference_branch, cfg.target_radius)
    count = length(solutions)
    radius = getfield.(solutions, :R)
    alpha = getfield.(solutions, :alpha)
    n_factor = zeros(count)
    for j in 2:count
        delta_radius = radius[j] - radius[j - 1]
        n_factor[j] = n_factor[j - 1] -
            0.5 * (imag(alpha[j - 1]) + imag(alpha[j])) * delta_radius
    end
    gain = exp.(n_factor)
    amplitude = abs(receptivity.C_r) .* gain
    return (
        label = label,
        start_radius = seed.R,
        initial_Cr = receptivity.C_r,
        solutions = solutions,
        n_factor = n_factor,
        gain = gain,
        amplitude = amplitude,
        adjoint_error = abs(
            receptivity.alpha_adjoint - conj(seed.alpha)
        ),
    )
end

function write_header(io, title, cfg)
    println(io, "TITLE = \"", title, "\"")
    println(
        io,
        "VARIABLES = \"R\" \"N_factor\" \"A_abs\" \"gain\" " *
        "\"alpha_r\" \"alpha_i\" \"growth_rate\" \"Cr_initial_abs\""
    )
    println(io, "DATASETAUXDATA N_cheb = \"", cfg.N_cheb, "\"")
    println(io, "DATASETAUXDATA R_target = \"", cfg.target_radius, "\"")
    println(io, "DATASETAUXDATA radial_step = \"", cfg.integration_step, "\"")
    println(io, "DATASETAUXDATA n = \"", cfg.azimuthal_mode, "\"")
    println(io, "DATASETAUXDATA omega_bar = \"", cfg.omega_bar, "\"")
    println(io, "DATASETAUXDATA roughness_ls = \"", cfg.roughness_localization, "\"")
end

function write_zone(io, result)
    println(
        io,
        "ZONE T=\"", result.label, "\", I=",
        length(result.solutions), ", F=POINT"
    )
    initial_amplitude = abs(result.initial_Cr)
    for j in eachindex(result.solutions)
        solution = result.solutions[j]
        @printf(
            io,
            "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
            solution.R,
            result.n_factor[j],
            result.amplitude[j],
            result.gain[j],
            real(solution.alpha),
            imag(solution.alpha),
            -imag(solution.alpha),
            initial_amplitude,
        )
    end
end

function write_outputs(results, cfg)
    output_directory = abspath(C6_OUTPUT_DIRECTORY)
    mkpath(output_directory)

    combined_path = joinpath(
        output_directory, "c6_three_start_A_N_evolution.dat"
    )
    open(combined_path, "w") do io
        write_header(
            io,
            "C6 three-start downstream amplitude and N-factor evolution",
            cfg,
        )
        for result in results
            write_zone(io, result)
        end
    end

    individual_paths = String[]
    for result in results
        path = joinpath(
            output_directory, "c6_" * result.label * "_A_N.dat"
        )
        open(path, "w") do io
            write_header(
                io,
                "C6 downstream amplitude and N-factor: " * result.label,
                cfg,
            )
            write_zone(io, result)
        end
        push!(individual_paths, path)
    end

    summary_path = joinpath(
        output_directory, "c6_three_start_A_N_summary.txt"
    )
    open(summary_path, "w") do io
        println(io, "C6 three-start downstream evolution")
        println(io, "N_cheb = ", cfg.N_cheb)
        println(io, "R_target = ", cfg.target_radius)
        println(io, "radial_step = ", cfg.integration_step)
        println(io, "a_s = ", cfg.mass_flux)
        println(io, "Re_h = ", cfg.Re_h)
        println(io, "omega_bar = ", cfg.omega_bar)
        println(io, "n = ", cfg.azimuthal_mode)
        println(io, "roughness_height = ", cfg.roughness_height)
        println(io, "roughness_localization = ", cfg.roughness_localization)
        println(
            io,
            "middle_case_note = R=372.988887508 is the maximum spatial " *
            "growth-rate location on the continuously tracked branch."
        )
        println(io)
        println(
            io,
            "# label R_start alpha_r_start alpha_i_start Cr_initial_abs " *
            "N_at_570 gain_at_570 A_abs_at_570 points adjoint_error"
        )
        for result in results
            seed = result.solutions[1]
            @printf(
                io,
                "%s %.12e %.12e %.12e %.12e %.12e %.12e %.12e %d %.12e\n",
                result.label,
                result.start_radius,
                real(seed.alpha),
                imag(seed.alpha),
                abs(result.initial_Cr),
                result.n_factor[end],
                result.gain[end],
                result.amplitude[end],
                length(result.solutions),
                result.adjoint_error,
            )
        end
        println(io)
        println(io, "combined_tecplot_file = ", combined_path)
        for path in individual_paths
            println(io, "individual_tecplot_file = ", path)
        end
    end

    return combined_path, individual_paths, summary_path
end

function main()
    cfg = Config(
        target_radius = C6_TARGET_RADIUS,
        N_cheb = 99,
        scan_start = 240.0,
        scan_step = 2.0,
        integration_step = C6_RADIAL_STEP,
        candidate_count = 1,
        roughness_localization = 0.5,
    )
    cfg.roughness_height = 1 / cfg.Re_h

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(
        u0, v0, w0, z, cfg.N_cheb, cfg.mode
    )
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    println("Locating the lower neutral point...")
    neutral = find_lower_neutral(cfg, F, G, H, D, D2)

    cfg.candidate_count = 3
    println("Building the continuously tracked N=99 reference branch...")
    reference_branch = tracked_curve(neutral, cfg, F, G, H, D, D2)
    growth_seed = solution_on_reference(
        C6_MAXIMUM_GROWTH_RADIUS,
        reference_branch,
        cfg, F, G, H, D, D2,
    )
    global_seed = solution_on_reference(
        C6_GLOBAL_PEAK_RADIUS,
        reference_branch,
        cfg, F, G, H, D, D2,
    )

    seeds = (
        ("lower_neutral", neutral),
        ("maximum_growth_Rg", growth_seed),
        ("global_Cr_peak", global_seed),
    )
    results = []
    for (label, seed) in seeds
        println("Integrating ", label, " from R=", seed.R, "...")
        push!(
            results,
            propagate_case(
                label, seed, reference_branch,
                cfg, F, G, D, D2, mass_matrix
            )
        )
    end

    combined_path, individual_paths, summary_path = write_outputs(
        results, cfg
    )
    println("Combined Tecplot file: ", combined_path)
    for path in individual_paths
        println("Individual Tecplot file: ", path)
    end
    println("Summary: ", summary_path)
end

main()
