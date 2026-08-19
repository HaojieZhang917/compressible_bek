using Printf

include("total_amplitude.jl")

const OUTPUT_DIRECTORY = "c_6_epsilon_results"
const TARGET_RADIUS = 570.0
const RADIAL_STEP = 2.0
const N_CHEB = 99

function normalized_velocity_fields(solution, n_cheb)
    fields = full_mode(copy(solution.vector), n_cheb)
    v_scale = maximum(abs.(fields[2]))
    v_scale > 0 || error("The direct-mode azimuthal velocity scale is zero.")
    return fields[1] ./ v_scale, fields[2] ./ v_scale, fields[3] ./ v_scale
end

function velocity_magnitude(u, v, w)
    return sqrt.(abs2.(u) .+ abs2.(v) .+ abs2.(w))
end

function threshold_radius(radius, epsilon, threshold)
    index = findfirst(epsilon .>= threshold)
    index === nothing && return nothing
    index == 1 && return radius[1]

    fraction = (threshold - epsilon[index - 1]) /
        (epsilon[index] - epsilon[index - 1])
    return radius[index - 1] + fraction * (radius[index] - radius[index - 1])
end

function write_results(path, cfg, neutral, receptivity, radius, n_factor, epsilon,
                       mode_peak, u_ref, component_peaks)
    open(path, "w") do io
        println(io, "TITLE = \"C6 lower-neutral epsilon evolution\"")
        println(io, "VARIABLES = \"R\" \"N_factor\" \"A_abs\" \"epsilon\" " *
                    "\"mode_velocity_peak\" \"U_ref\" \"u_peak\" \"v_peak\" \"w_peak\"")
        println(io, "DATASETAUXDATA N_cheb = \"", cfg.N_cheb, "\"")
        println(io, "DATASETAUXDATA radial_step = \"", cfg.integration_step, "\"")
        println(io, "DATASETAUXDATA R_lower = \"", neutral.R, "\"")
        println(io, "DATASETAUXDATA abs_Cr_lower = \"", abs(receptivity.C_r), "\"")
        println(io, "DATASETAUXDATA U_ref_definition = \"max_z sqrt(F^2+G^2+H^2)\"")
        println(io, "DATASETAUXDATA epsilon_definition = \"abs(Cr) exp(N) max_z sqrt(abs(u_hat)^2+abs(v_hat)^2+abs(w_hat)^2) / U_ref\"")
        println(io, "ZONE T=\"lower_neutral\", I=", length(radius), ", F=POINT")

        for index in eachindex(radius)
            amplitude = abs(receptivity.C_r) * exp(n_factor[index])
            @printf(io, "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                    radius[index], n_factor[index], amplitude, epsilon[index],
                    mode_peak[index], u_ref[index], component_peaks[index, 1],
                    component_peaks[index, 2], component_peaks[index, 3])
        end
    end
end

function write_summary(path, cfg, neutral, receptivity, radius, n_factor, epsilon)
    open(path, "w") do io
        println(io, "C6 lower-neutral epsilon evolution")
        println(io, "Re_h = ", cfg.Re_h)
        println(io, "a_s = ", cfg.mass_flux)
        println(io, "omega_bar = ", cfg.omega_bar)
        println(io, "n = ", cfg.azimuthal_mode)
        println(io, "N_cheb = ", cfg.N_cheb)
        println(io, "radial_step = ", cfg.integration_step)
        println(io, "R_lower = ", neutral.R)
        println(io, "abs_Cr_lower = ", abs(receptivity.C_r))
        println(io, "N_at_target = ", n_factor[end])
        println(io, "epsilon_at_lower = ", epsilon[1])
        println(io, "epsilon_at_target = ", epsilon[end])
        println(io, "U_ref_definition = max_z sqrt(F^2 + G^2 + H^2)")
        println(io, "epsilon_definition = abs(Cr) exp(N) max_z sqrt(abs(u_hat)^2 + abs(v_hat)^2 + abs(w_hat)^2) / U_ref")
        for threshold in (0.01, 0.05, 0.1)
            crossing = threshold_radius(radius, epsilon, threshold)
            println(io, "R_at_epsilon_", threshold, " = ",
                    crossing === nothing ? "not reached" : crossing)
        end
    end
end

function main()
    cfg = Config(
        target_radius = TARGET_RADIUS,
        N_cheb = N_CHEB,
        scan_start = 240.0,
        scan_step = 2.0,
        integration_step = RADIAL_STEP,
        candidate_count = 1,
        roughness_localization = 0.5,
    )
    cfg.roughness_height = 1 / cfg.Re_h

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(
        cfg.Re_h, cfg.Ro, cfg.mass_flux, cfg.mode
    )
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, cfg.mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, cfg.mode)
    mass_matrix = cheb_mass_matrix(cfg.N_cheb, cfg.mode)

    println("Locating the lower neutral point...")
    neutral = find_lower_neutral(cfg, F, G, H, D, D2)
    receptivity = initial_receptivity(neutral, cfg, F, G, D, D2, mass_matrix)
    base_velocity = velocity_magnitude(F, G, H)
    reference_velocity = maximum(base_velocity)
    radius = radial_points(neutral.R, cfg.target_radius, cfg.integration_step)
    count = length(radius)
    epsilon = zeros(count)
    mode_peak = zeros(count)
    u_ref = fill(reference_velocity, count)
    component_peaks = zeros(count, 3)
    n_factor = zeros(count)
    previous = neutral

    println("Following the lower-neutral branch and reconstructing epsilon...")
    for index in eachindex(radius)
        current = if index == 1
            neutral
        else
            direct_solution(
                radius[index], previous.alpha, previous.vector,
                cfg, F, G, H, D, D2,
            )
        end
        if index > 1
            delta_radius = radius[index] - radius[index - 1]
            n_factor[index] = n_factor[index - 1] -
                0.5 * (imag(previous.alpha) + imag(current.alpha)) * delta_radius
        end
        u_hat, v_hat, w_hat = normalized_velocity_fields(current, cfg.N_cheb)
        component_peaks[index, :] .= (maximum(abs.(u_hat)), maximum(abs.(v_hat)), maximum(abs.(w_hat)))
        mode_peak[index] = maximum(velocity_magnitude(u_hat, v_hat, w_hat))
        amplitude = abs(receptivity.C_r) * exp(n_factor[index])
        epsilon[index] = amplitude * mode_peak[index] / reference_velocity
        previous = current
        index % 20 == 0 && GC.gc()
    end

    output_directory = abspath(OUTPUT_DIRECTORY)
    mkpath(output_directory)
    data_path = joinpath(output_directory, "c6_lower_neutral_epsilon_N99_dr2.dat")
    summary_path = joinpath(output_directory, "c6_lower_neutral_epsilon_N99_dr2_summary.txt")
    write_results(data_path, cfg, neutral, receptivity, radius, n_factor, epsilon,
                  mode_peak, u_ref, component_peaks)
    write_summary(summary_path, cfg, neutral, receptivity, radius, n_factor, epsilon)

    println("epsilon at R_lower = ", epsilon[1])
    println("epsilon at R_target = ", epsilon[end])
    for threshold in (0.01, 0.05, 0.1)
        crossing = threshold_radius(radius, epsilon, threshold)
        println("R at epsilon=", threshold, ": ",
                crossing === nothing ? "not reached" : crossing)
    end
    println("data output: ", data_path)
    println("summary output: ", summary_path)
end

main()
