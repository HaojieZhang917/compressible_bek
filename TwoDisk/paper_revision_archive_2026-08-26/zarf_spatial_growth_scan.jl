using LinearAlgebra
using NonlinearEigenproblems
using Printf
using Statistics

const DEFAULT_ZARF_ROOT = raw"/home/zhj/main/code/compress/compressible_bek/TwoDisk/data/stability/Zarf_neutral"
include(joinpath(DEFAULT_ZARF_ROOT, "BaseFlow_cavity.jl"))
include(joinpath(DEFAULT_ZARF_ROOT, "Stability_Cavity.jl"))

Base.@kwdef struct Config
    zarf_root::String = DEFAULT_ZARF_ROOT
    output_root::String = joinpath(@__DIR__, "zarf_spatial_growth_results")
    mass_fluxes::Vector{Float64} = [0.0]
    N_cheb::Int = 99
    Re_h::Int = 1000
    mode::Int = 1
    candidate_count::Int = 1
    eig_tolerance::Float64 = 1e-10
    max_iterations::Int = 400
    imaginary_shift::Float64 = -0.005
    frequency_samples::Int = 9
    full_frequency::Bool = false
    max_points::Int = 0
    use_newton::Bool = true
    newton_iterations::Int = 6
    newton_tolerance::Float64 = 1e-7
    stationarity_step::Float64 = 1e-5
    radius_min::Float64 = -Inf
    radius_max::Float64 = Inf
    beta_min::Float64 = -Inf
    beta_max::Float64 = Inf
    points_file::String = ""
end

function parse_cli(args)
    mass_fluxes = [0.0]
    zarf_root = DEFAULT_ZARF_ROOT
    output_root = joinpath(@__DIR__, "zarf_spatial_growth_results")
    N_cheb = 99
    frequency_samples = 9
    candidate_count = 1
    full_frequency = false
    max_points = 0
    use_newton = true
    radius_min = -Inf
    radius_max = Inf
    beta_min = -Inf
    beta_max = Inf
    points_file = ""
    i = 1
    while i <= length(args)
        key = args[i]
        if key == "--all"
            mass_fluxes = [-0.4, -0.2, 0.0, 0.2, 0.4]
            i += 1
            continue
        elseif key == "--remaining"
            mass_fluxes = [-0.4, -0.2, 0.2, 0.4]
            i += 1
            continue
        elseif key == "--help"
            println("""
Usage:
  julia -t auto zarf_spatial_growth_scan.jl [options]

Options:
  --ts VALUE          One mass-flux value (default: 0.0).
  --all               Calculate -0.4, -0.2, 0.0, 0.2, and 0.4.
  --remaining         Calculate -0.4, -0.2, 0.2, and 0.4.
  --zarf-root PATH    Directory containing Ts=*/R=*.dat.
  --output-root PATH  Output directory.
  --N VALUE           Chebyshev order (default: 99).
  --frequency-samples VALUE
                     Samples per continuous unstable frequency segment (default: 9).
  --candidate-count VALUE
                     Polynomial roots requested around each shift (default: 1).
  --full-frequency   Use every omega_i>0 sample (brute-force reference).
  --max-points VALUE  Limit rows for a runtime test; 0 means all rows.
  --no-newton        Disable local frequency optimization.
  --R VALUE          Restrict the scan to one radius.
  --R-min VALUE      Minimum radius included in the scan.
  --R-max VALUE      Maximum radius included in the scan.
  --beta VALUE       Restrict the scan to one beta value.
  --points-file PATH Restrict the scan to R beta pairs listed in a text file.
""")
            exit()
        end
        if key == "--full-frequency"
            full_frequency = true
            i += 1
            continue
        elseif key == "--no-newton"
            use_newton = false
            i += 1
            continue
        end
        i == length(args) && error("Missing value after $key")
        value = args[i + 1]
        if key == "--ts"
            mass_fluxes = [parse(Float64, value)]
        elseif key == "--zarf-root"
            zarf_root = value
        elseif key == "--output-root"
            output_root = value
        elseif key == "--N"
            N_cheb = parse(Int, value)
        elseif key == "--frequency-samples"
            frequency_samples = parse(Int, value)
            frequency_samples >= 2 || error("--frequency-samples must be at least 2")
        elseif key == "--candidate-count"
            candidate_count = parse(Int, value)
            candidate_count >= 1 || error("--candidate-count must be positive")
        elseif key == "--max-points"
            max_points = parse(Int, value)
        elseif key == "--R"
            radius_min = parse(Float64, value)
            radius_max = radius_min
        elseif key == "--R-min"
            radius_min = parse(Float64, value)
        elseif key == "--R-max"
            radius_max = parse(Float64, value)
        elseif key == "--beta"
            beta_min = parse(Float64, value)
            beta_max = beta_min
        elseif key == "--points-file"
            points_file = value
        else
            error("Unknown option $key")
        end
        i += 2
    end
    return Config(
        zarf_root=zarf_root,
        output_root=output_root,
        mass_fluxes=mass_fluxes,
        N_cheb=N_cheb,
        frequency_samples=frequency_samples,
        candidate_count=candidate_count,
        full_frequency=full_frequency,
        max_points=max_points,
        use_newton=use_newton,
        radius_min=radius_min,
        radius_max=radius_max,
        beta_min=beta_min,
        beta_max=beta_max,
        points_file=points_file,
    )
end

function temporal_directory(root, mass_flux)
    label = @sprintf("Ts=%.1f", mass_flux)
    return joinpath(root, label)
end

function parse_temporal_file(path)
    samples = NamedTuple[]
    for line in eachline(path)
        values = split(strip(line))
        length(values) < 5 && continue
        parsed = try
            parse.(Float64, values[1:5])
        catch
            continue
        end
        parsed[5] > 0 || continue
        push!(samples, (
            R=parsed[1],
            alpha_t=parsed[2],
            beta=parsed[3],
            omega_r=parsed[4],
            omega_i=parsed[5],
        ))
    end
    return samples
end

function temporal_paths(root, mass_flux)
    directory = temporal_directory(root, mass_flux)
    isdir(directory) || error("Temporal data directory does not exist: $directory")
    paths = filter(path -> startswith(basename(path), "R=") && endswith(path, ".dat"),
                   readdir(directory; join=true))
    isempty(paths) && error("No R=*.dat files found in $directory")
    sort!(paths; by=path -> parse(Float64,
        replace(basename(path), r"^R=|\.dat$" => "")))
    return paths
end

function frequency_segments(rows)
    sort!(rows; by=x -> x.omega_r)
    length(rows) <= 1 && return [rows]
    spacings = diff([row.omega_r for row in rows])
    positive_spacings = filter(>(0.0), spacings)
    threshold = isempty(positive_spacings) ? Inf : 2.5 * median(positive_spacings)
    segments = Vector{typeof(rows)}()
    first_index = 1
    for index in 2:length(rows)
        if rows[index].omega_r - rows[index - 1].omega_r > threshold
            push!(segments, rows[first_index:index - 1])
            first_index = index
        end
    end
    push!(segments, rows[first_index:end])
    return segments
end

function evenly_select(rows, count)
    length(rows) <= count && return rows
    indices = unique(round.(Int, range(1, length(rows); length=count)))
    return rows[indices]
end

function select_frequency_samples(samples, cfg)
    grouped = Dict{Tuple{Float64, Float64}, Vector{NamedTuple}}()
    for sample in samples
        push!(get!(grouped, (sample.R, sample.beta), NamedTuple[]), sample)
    end
    selected = NamedTuple[]
    branch_count = 0
    for key in sort(collect(keys(grouped)))
        for segment in frequency_segments(grouped[key])
            branch_count += 1
            chosen = cfg.full_frequency ? segment : evenly_select(segment, cfg.frequency_samples)
            append!(selected, [merge(row, (branch_id=branch_count,)) for row in chosen])
        end
    end
    sort!(selected; by=x -> (x.R, x.beta, x.omega_r, x.alpha_t))
    if cfg.max_points > 0 && length(selected) > cfg.max_points
        resize!(selected, cfg.max_points)
    end
    return selected, branch_count
end

function load_temporal_samples(root, mass_flux, cfg)
    paths = temporal_paths(root, mass_flux)
    samples = NamedTuple[]
    for path in paths
        append!(samples, parse_temporal_file(path))
    end
    filter!(sample ->
        cfg.radius_min <= sample.R <= cfg.radius_max &&
        cfg.beta_min <= sample.beta <= cfg.beta_max,
        samples,
    )
    if !isempty(cfg.points_file)
        point_keys = Set{Tuple{Float64, Float64}}()
        for line in eachline(cfg.points_file)
            values = split(strip(line))
            length(values) < 2 && continue
            point = try
                (round(parse(Float64, values[1]); digits=10),
                 round(parse(Float64, values[2]); digits=10))
            catch
                continue
            end
            push!(point_keys, point)
        end
        isempty(point_keys) && error("No R beta pairs found in $(cfg.points_file)")
        filter!(sample ->
            (round(sample.R; digits=10), round(sample.beta; digits=10)) in point_keys,
            samples,
        )
    end
    selected, branch_count = select_frequency_samples(samples, cfg)
    isempty(selected) && error("No unstable temporal samples found in $(temporal_directory(root, mass_flux))")
    return selected, branch_count, paths
end

function group_by_radius(samples)
    groups = Dict{Float64, Vector{Int}}()
    for (index, sample) in enumerate(samples)
        push!(get!(groups, sample.R, Int[]), index)
    end
    return sort(collect(groups); by=first)
end

function root_indices(values, shift, alpha_t)
    candidates = findall(value ->
        isfinite(real(value)) && isfinite(imag(value)) &&
        real(value) > 0 && abs(imag(value)) < 0.5,
        values,
    )
    isempty(candidates) && (candidates = collect(eachindex(values)))
    sort!(candidates; by=index -> begin
        value = values[index]
        abs(real(value) - alpha_t) + 0.35abs(imag(value) - imag(shift))
    end)
    return candidates
end

function solve_spatial(sample, cof, D, D2, cfg)
    L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(
        cof, D, D2, sample.beta, sample.omega_r, sample.R
    )
    L0, L1, L2 = CRC_STA.boudary_condition(
        L0_raw, L1_raw, L2_raw, cfg.N_cheb, cfg.mode
    )
    problem = PEP([L0, L1, L2])
    shift = complex(sample.alpha_t, cfg.imaginary_shift)
    values, vectors = try
        iar(
            problem;
            σ=shift,
            neigs=cfg.candidate_count,
            maxit=cfg.max_iterations,
            tol=cfg.eig_tolerance,
        )
    catch
        iar(
            problem;
            σ=complex(sample.alpha_t, 0.0),
            neigs=1,
            maxit=2cfg.max_iterations,
            tol=1e-8,
        )
    end
    roots = NamedTuple[]
    for index in root_indices(values, shift, sample.alpha_t)
        alpha = values[index]
        vector = vectors[:, index]
        residual = norm((L0 + alpha * L1 + alpha^2 * L2) * vector) /
                   max(norm(vector), eps())
        push!(roots, (alpha=alpha, residual=residual))
    end
    isempty(roots) && error("No finite spatial roots found")
    return roots
end

function spatial_matrices(cof, D, D2, beta, omega, radius, cfg)
    L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(
        cof, D, D2, beta, omega, radius
    )
    L0, L1, L2 = CRC_STA.boudary_condition(
        L0_raw, L1_raw, L2_raw, cfg.N_cheb, cfg.mode
    )
    Tomega_raw = -im * cof.Ta
    Tzero_raw = zeros(eltype(Tomega_raw), size(Tomega_raw))
    Tomega, _, _ = CRC_STA.boudary_condition(
        Tomega_raw, Tzero_raw, Tzero_raw, cfg.N_cheb, cfg.mode
    )
    return L0, L1, L2, Tomega
end

function polynomial_pair(L0, L1, L2, shift, cfg)
    problem = PEP([L0, L1, L2])
    values, vectors = iar(
        problem;
        σ=shift,
        neigs=cfg.candidate_count,
        maxit=cfg.max_iterations,
        tol=cfg.eig_tolerance,
    )
    index = first(root_indices(values, shift, real(shift)))
    alpha = values[index]
    vector = vectors[:, index]
    residual = norm((L0 + alpha * L1 + alpha^2 * L2) * vector) /
               max(norm(vector), eps())
    return alpha, vector, residual
end

function spatial_pair(sample, omega, shift, cof, D, D2, cfg)
    L0, L1, L2, Tomega = spatial_matrices(
        cof, D, D2, sample.beta, omega, sample.R, cfg
    )
    alpha, q, residual = try
        polynomial_pair(L0, L1, L2, shift, cfg)
    catch
        polynomial_pair(L0, L1, L2, complex(real(shift), 0.0), cfg)
    end
    return alpha, q, residual, L0, L1, L2, Tomega
end

function closest_root(values, shift)
    candidates = findall(value ->
        isfinite(real(value)) && isfinite(imag(value)), values
    )
    isempty(candidates) && error("No finite polynomial eigenvalue found")
    return candidates[argmin(abs.(values[candidates] .- shift))]
end

function spatial_derivative(alpha, q, L0, L1, L2, Tomega, cfg)
    adjoint_problem = PEP([adjoint(L0), adjoint(L1), adjoint(L2)])
    adjoint_values, adjoint_vectors = iar(
        adjoint_problem;
        σ=conj(alpha),
        neigs=cfg.candidate_count,
        maxit=cfg.max_iterations,
        tol=cfg.eig_tolerance,
    )
    adjoint_index = closest_root(adjoint_values, conj(alpha))
    p = adjoint_vectors[:, adjoint_index]
    Palpha = L1 + 2alpha * L2
    numerator = dot(p, Tomega * q)
    denominator = dot(p, Palpha * q)
    abs(denominator) > eps() || error("Singular polynomial derivative denominator")
    return -numerator / denominator
end

function evaluate_frequency(sample, omega, shift, cof, D, D2, cfg)
    alpha, q, residual, L0, L1, L2, Tomega = spatial_pair(
        sample, omega, shift, cof, D, D2, cfg
    )
    derivative = spatial_derivative(
        alpha, q, L0, L1, L2, Tomega, cfg
    )
    return (
        omega=omega,
        alpha=alpha,
        growth=-imag(alpha),
        derivative=derivative,
        residual=residual,
    )
end

function newton_peak(sample, left, seed, right, cof, D, D2, cfg)
    left.omega_r < seed.omega_r < right.omega_r || return nothing
    design = [
        1.0 left.omega_r left.omega_r^2
        1.0 seed.omega_r seed.omega_r^2
        1.0 right.omega_r right.omega_r^2
    ]
    coefficients = design \ [
        left.growth_rate, seed.growth_rate, right.growth_rate
    ]
    curvature = 2coefficients[3]
    isfinite(curvature) && curvature < 0 || return nothing
    vertex = -coefficients[2] / curvature
    vertex = clamp(vertex, left.omega_r, right.omega_r)
    candidate = try
        evaluate_frequency(
            sample, vertex, complex(seed.alpha_r, seed.alpha_i),
            cof, D, D2, cfg,
        )
    catch
        return nothing
    end
    best = candidate
    abs(imag(candidate.derivative)) <= cfg.newton_tolerance && return best

    corrected_omega = vertex + imag(candidate.derivative) / curvature
    if !isfinite(corrected_omega) ||
       corrected_omega <= left.omega_r || corrected_omega >= right.omega_r
        return best
    end
    corrected = try
        alpha, _, residual, _, _, _, _ = spatial_pair(
            sample, corrected_omega, candidate.alpha, cof, D, D2, cfg
        )
        (
            omega=corrected_omega,
            alpha=alpha,
            growth=-imag(alpha),
            derivative=candidate.derivative,
            residual=residual,
        )
    catch
        return best
    end
    corrected.growth > best.growth && (best = corrected)
    return best
end

function calculate_case(mass_flux, cfg)
    samples, branch_count, paths = load_temporal_samples(cfg.zarf_root, mass_flux, cfg)
    groups = group_by_radius(samples)

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, -1, mass_flux, 1)
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, 1)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, 1)

    group_results = Vector{Vector{NamedTuple}}(undef, length(groups))
    completed = Threads.Atomic{Int}(0)
    total_groups = length(groups)
    Threads.@threads for group_index in eachindex(groups)
        radius, indices = groups[group_index]
        cof = CRC_STA.Spatial_mode_BEK1(
            F, G .- 1, H, radius, cfg.N_cheb, D, D2, cfg.Re_h
        )
        branch_indices = Dict{Int, Vector{Int}}()
        for index in indices
            push!(get!(branch_indices, samples[index].branch_id, Int[]), index)
        end
        local_results = NamedTuple[]
        for branch_index in sort(collect(keys(branch_indices)))
            branch_rows = sort(branch_indices[branch_index];
                               by=index -> samples[index].omega_r)
            coarse = NamedTuple[]
            for index in branch_rows
                sample = samples[index]
                try
                    roots = solve_spatial(sample, cof, D, D2, cfg)
                    for root in roots
                        alpha = root.alpha
                        push!(coarse, merge(sample, (
                            alpha_r=real(alpha),
                            alpha_i=imag(alpha),
                            growth_rate=-imag(alpha),
                            omega_bar=sample.R * sample.omega_r,
                            residual=root.residual,
                            status=1,
                        )))
                    end
                catch exception
                    @warn "Spatial solve failed" mass_flux radius index exception
                end
            end
            append!(local_results, coarse)
            if cfg.use_newton && length(coarse) >= 3
                for local_index in 2:(length(coarse) - 1)
                    left = coarse[local_index - 1]
                    seed = coarse[local_index]
                    right = coarse[local_index + 1]
                    seed.growth_rate >= left.growth_rate || continue
                    seed.growth_rate >= right.growth_rate || continue
                    peak = try
                        newton_peak(
                            seed, left, seed, right, cof, D, D2, cfg,
                        )
                    catch exception
                        @warn "Newton frequency optimization failed" mass_flux radius branch_index local_index exception
                        nothing
                    end
                    peak === nothing && continue
                    push!(local_results, merge(seed, (
                        omega_r=peak.omega,
                        omega_bar=seed.R * peak.omega,
                        alpha_r=real(peak.alpha),
                        alpha_i=imag(peak.alpha),
                        growth_rate=peak.growth,
                        residual=peak.residual,
                        status=1,
                    )))
                end
            end
        end
        group_results[group_index] = local_results
        done = Threads.atomic_add!(completed, 1) + 1
        if done == 1 || done % 5 == 0 || done == total_groups
            @printf("a_s=%+.1f: completed %d/%d radius groups\n",
                    mass_flux, done, total_groups)
            flush(stdout)
        end
    end
    results = reduce(vcat, group_results; init=NamedTuple[])
    sort!(results; by=x -> (x.R, x.beta, x.omega_r, x.alpha_t))
    return join(paths, ";"), branch_count, results
end

function envelope(results)
    best = Dict{Tuple{Float64, Float64}, NamedTuple}()
    for row in results
        row.status == 1 || continue
        key = (row.R, row.beta)
        if !haskey(best, key) || row.growth_rate > best[key].growth_rate
            best[key] = row
        end
    end
    return sort(collect(values(best)); by=x -> (x.R, x.beta))
end

function write_samples(path, source, mass_flux, cfg, results)
    open(path, "w") do io
        println(io, "TITLE=\"Frequency-resolved spatial growth inside temporal Zarf\"")
        println(io,
            "VARIABLES=\"R\",\"beta\",\"omega\",\"omega_bar\",\"branch_id\"," *
            "\"alpha_temporal\",\"omega_i_temporal\",\"alpha_r\",\"alpha_i\"," *
            "\"growth_rate\",\"residual\",\"status\""
        )
        println(io, "DATASETAUXDATA source=\"$(replace(source, '"' => '\''))\"")
        println(io, "DATASETAUXDATA mass_flux=\"$mass_flux\"")
        println(io, "DATASETAUXDATA N_cheb=\"$(cfg.N_cheb)\"")
        println(io, "DATASETAUXDATA frequency_samples=\"$(cfg.frequency_samples)\"")
        println(io, "DATASETAUXDATA full_frequency=\"$(cfg.full_frequency)\"")
        println(io, "ZONE T=\"a_s=$(mass_flux)\", I=$(length(results)), F=POINT")
        for row in results
            @printf(io,
                "%.12e %.12e %.12e %.12e %d %.12e %.12e %.12e %.12e %.12e %.12e %d\n",
                row.R, row.beta, row.omega_r, row.omega_bar, row.branch_id,
                row.alpha_t, row.omega_i, row.alpha_r, row.alpha_i,
                row.growth_rate, row.residual, row.status,
            )
        end
    end
end

function write_envelope(path, source, mass_flux, cfg, rows)
    open(path, "w") do io
        println(io, "TITLE=\"Maximum local spatial growth inside temporal Zarf\"")
        println(io,
            "VARIABLES=\"R\",\"beta\",\"growth_rate_max\",\"omega_opt\",\"branch_id\"," *
            "\"omega_bar_opt\",\"alpha_r\",\"alpha_i\",\"alpha_temporal\"," *
            "\"omega_i_temporal\",\"residual\""
        )
        println(io, "DATASETAUXDATA source=\"$(replace(source, '"' => '\''))\"")
        println(io, "DATASETAUXDATA mass_flux=\"$mass_flux\"")
        println(io, "DATASETAUXDATA N_cheb=\"$(cfg.N_cheb)\"")
        println(io, "ZONE T=\"a_s=$(mass_flux)\", I=$(length(rows)), F=POINT")
        for row in rows
            @printf(io,
                "%.12e %.12e %.12e %.12e %d %.12e %.12e %.12e %.12e %.12e %.12e\n",
                row.R, row.beta, row.growth_rate, row.omega_r, row.branch_id,
                row.omega_bar, row.alpha_r, row.alpha_i, row.alpha_t,
                row.omega_i, row.residual,
            )
        end
    end
end

function main()
    cfg = parse_cli(ARGS)
    isdir(cfg.zarf_root) || error("Zarf root does not exist: $(cfg.zarf_root)")
    mkpath(cfg.output_root)
    LinearAlgebra.BLAS.set_num_threads(1)

    for mass_flux in cfg.mass_fluxes
        label = replace(@sprintf("%+.1f", mass_flux), "+" => "p", "-" => "m", "." => "p")
        output_directory = joinpath(cfg.output_root, "Ts_$(label)_N$(cfg.N_cheb)")
        mkpath(output_directory)
        source, branch_count, results = calculate_case(mass_flux, cfg)
        envelope_rows = envelope(results)
        sample_path = joinpath(output_directory, "spatial_growth_samples.dat")
        envelope_path = joinpath(output_directory, "spatial_growth_envelope.dat")
        write_samples(sample_path, source, mass_flux, cfg, results)
        write_envelope(envelope_path, source, mass_flux, cfg, envelope_rows)
        successful = count(row -> row.status == 1, results)
        @printf("a_s=%+.1f complete: %d/%d samples, %d envelope points\n",
                mass_flux, successful, length(results), length(envelope_rows))
        println("frequency branches: ", branch_count)
        println("samples:  ", sample_path)
        println("envelope: ", envelope_path)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
