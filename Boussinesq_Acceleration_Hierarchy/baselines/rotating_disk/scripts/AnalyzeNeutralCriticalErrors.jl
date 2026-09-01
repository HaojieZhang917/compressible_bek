using LinearAlgebra
using Printf

const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT_DIR = get(
    ENV, "NEUTRAL_CURVE_DIR", joinpath(WORKSPACE_ROOT, "neutral_curve_integrated"),
)
const OUTPUT_DIR = get(ENV, "CRITICAL_ERROR_DIR", INPUT_DIR)
const TYPE_BETA_SPLIT = 0.055
const CONFIG_ORDER = Dict(
    "lopez" => 1,
    "compressible_full" => 2,
    "compressible_propPert_off_baseProp_variable" => 3,
    "compressible_propPert_off_baseProp_frozen" => 4,
)

function load_curve(path::AbstractString)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        text = strip(line)
        (isempty(text) || startswith(lowercase(text), "variables") ||
         startswith(lowercase(text), "zone") || startswith(text, "#")) && continue
        fields = split(text)
        length(fields) == 7 || error("Expected seven columns in $path: $text")
        push!(rows, parse.(Float64, fields))
    end
    isempty(rows) && return zeros(Float64, 0, 7)
    return reduce(vcat, permutedims.(rows))
end

function parse_case(path::AbstractString)
    name = basename(path)
    tw_match = match(r"Tw=([0-9]+(?:\.[0-9]+)?)", name)
    tw_match === nothing && error("Cannot parse Tw from $name")
    Tw = parse(Float64, tw_match.captures[1])
    configuration = if occursin("model=lopez", name)
        "lopez"
    elseif occursin("model=compressible", name) &&
           occursin("propPert=on", name) && occursin("baseProp=variable", name)
        "compressible_full"
    elseif occursin("model=compressible", name) &&
           occursin("propPert=off", name) && occursin("baseProp=variable", name)
        "compressible_propPert_off_baseProp_variable"
    elseif occursin("model=compressible", name) &&
           occursin("propPert=off", name) && occursin("baseProp=frozen", name)
        "compressible_propPert_off_baseProp_frozen"
    else
        error("Unknown model/switch configuration in $name")
    end
    return (Tw=Tw, configuration=configuration, file=name)
end

function split_merged_segments(data::AbstractMatrix; tolerance=1.0e-10)
    size(data, 1) >= 3 || return [Matrix(data)]
    # Integrated files store each source branch with decreasing beta. A positive
    # beta jump marks the boundary between an appended Type-I and Type-II source.
    breaks = findall(diff(data[:, 3]) .> tolerance)
    starts = vcat(1, breaks .+ 1)
    stops = vcat(breaks, size(data, 1))
    return [Matrix(data[first:last, :]) for (first, last) in zip(starts, stops)]
end

function refine_minimum(segment::AbstractMatrix, index::Int)
    beta0 = segment[index, 3]
    R0 = segment[index, 2]
    indices = index-1:index+1
    x = segment[indices, 3] .- beta0
    y = segment[indices, 2]
    coefficients = hcat(x.^2, x, ones(3)) \ y
    a, b, c = coefficients
    if !(isfinite(a) && isfinite(b) && isfinite(c) && a > 0)
        return (R=R0, beta=beta0, fitted=false)
    end
    vertex = -b / (2a)
    lower, upper = extrema(x)
    if !(lower <= vertex <= upper)
        return (R=R0, beta=beta0, fitted=false)
    end
    R_vertex = a * vertex^2 + b * vertex + c
    isfinite(R_vertex) || return (R=R0, beta=beta0, fitted=false)
    return (R=R_vertex, beta=beta0 + vertex, fitted=true)
end

function critical_candidates(data::AbstractMatrix)
    candidates = NamedTuple[]
    for (segment_index, segment) in enumerate(split_merged_segments(data))
        size(segment, 1) >= 3 || continue
        for index in 2:size(segment, 1)-1
            R_left = segment[index-1, 2]
            R_center = segment[index, 2]
            R_right = segment[index+1, 2]
            (R_center < R_left && R_center < R_right) || continue
            refined = refine_minimum(segment, index)
            mode = refined.beta >= TYPE_BETA_SPLIT ? "Type-I" : "Type-II"
            push!(candidates, (
                mode=mode, R=refined.R, beta=refined.beta,
                R_discrete=R_center, beta_discrete=segment[index, 3],
                fitted=refined.fitted, segment=segment_index, index=index,
            ))
        end
    end
    return candidates
end

function select_critical_points(data::AbstractMatrix)
    candidates = critical_candidates(data)
    selected = Dict{String, NamedTuple}()
    for mode in ("Type-I", "Type-II")
        mode_candidates = filter(candidate -> candidate.mode == mode, candidates)
        isempty(mode_candidates) && continue
        selected[mode] = mode_candidates[argmin(getfield.(mode_candidates, :R))]
    end
    return selected
end

signed_percent(value, reference) = 100 * (value - reference) / reference

function format_value(value)
    value === missing && return ""
    value isa Bool && return string(value)
    value isa AbstractString && return value
    return @sprintf("%.12g", value)
end

function write_tsv(path, rows)
    columns = (
        :Tw, :mode, :configuration, :status, :R_c, :beta_c,
        :R_c_discrete, :beta_c_discrete, :quadratic_fit,
        :R_c_reference, :beta_c_reference, :delta_R_c,
        :R_c_error_percent, :abs_R_c_error_percent, :delta_beta_c,
        :beta_c_error_percent, :abs_beta_c_error_percent, :source_file,
    )
    open(path, "w") do io
        println(io, join(string.(columns), '\t'))
        for row in rows
            println(io, join((format_value(getproperty(row, column)) for column in columns), '\t'))
        end
    end
end

function write_dat(path, rows)
    available = filter(row -> row.status == "available", rows)
    variables = (
        "Tw", "R_c", "beta_c", "R_c_reference", "beta_c_reference",
        "delta_R_c", "R_c_error_percent", "abs_R_c_error_percent",
        "delta_beta_c", "beta_c_error_percent", "abs_beta_c_error_percent",
    )
    open(path, "w") do io
        println(io, "TITLE = \"Neutral-curve critical-point errors\"")
        println(io, "VARIABLES = ", join(("\"$name\"" for name in variables), ' '))
        groups = Dict{Tuple{String, String}, Vector{NamedTuple}}()
        for row in available
            push!(get!(groups, (row.configuration, row.mode), NamedTuple[]), row)
        end
        keys_ordered = sort(collect(keys(groups)); by=key -> (
            CONFIG_ORDER[key[1]], key[2] == "Type-I" ? 1 : 2,
        ))
        for key in keys_ordered
            group = sort(groups[key]; by=row -> row.Tw)
            println(
                io,
                "ZONE T=\"$(key[1]) $(key[2])\", I=$(length(group)), DATAPACKING=POINT",
            )
            for row in group
                values = (
                    row.Tw, row.R_c, row.beta_c, row.R_c_reference,
                    row.beta_c_reference, row.delta_R_c, row.R_c_error_percent,
                    row.abs_R_c_error_percent, row.delta_beta_c,
                    row.beta_c_error_percent, row.abs_beta_c_error_percent,
                )
                println(io, join((@sprintf("%.12e", value) for value in values), ' '))
            end
        end
    end
end

function write_report(path, rows, skipped_typeII_temperatures)
    available = filter(row -> row.status == "available", rows)
    open(path, "w") do io
        println(io, "# Neutral-curve critical-point error summary")
        println(io)
        println(io, "The Lopez model is the reference at each wall temperature.")
        println(io, "A critical point is an interior local minimum of R(beta). Each")
        println(io, "reported minimum is refined with a three-point quadratic fit.")
        println(io)
        println(io, "Signed errors are defined as")
        println(io)
        println(io, "`100 * (model - reference) / reference`.")
        println(io)
        if isempty(skipped_typeII_temperatures)
            println(io, "A full-compressible Type-II minimum exists at every sampled Tw.")
        else
            values = join((@sprintf("%.2f", Tw) for Tw in skipped_typeII_temperatures), ", ")
            println(io, "Type-II comparisons were omitted at Tw = $values because the")
            println(io, "full-compressible curve has no interior Type-II minimum.")
        end
        println(io)
        println(io, "## Error ranges")
        println(io)
        println(io, "| Configuration | Mode | Tw range | max abs Rc error | max abs beta error |")
        println(io, "|---|---|---:|---:|---:|")
        groups = Dict{Tuple{String, String}, Vector{NamedTuple}}()
        for row in available
            row.configuration == "lopez" && continue
            push!(get!(groups, (row.configuration, row.mode), NamedTuple[]), row)
        end
        keys_ordered = sort(collect(keys(groups)); by=key -> (
            CONFIG_ORDER[key[1]], key[2] == "Type-I" ? 1 : 2,
        ))
        for key in keys_ordered
            group = groups[key]
            tw_min = minimum(getfield.(group, :Tw))
            tw_max = maximum(getfield.(group, :Tw))
            max_R = maximum(getfield.(group, :abs_R_c_error_percent))
            max_beta = maximum(getfield.(group, :abs_beta_c_error_percent))
            @printf(
                io, "| `%s` | %s | %.2f-%.2f | %.3f%% | %.3f%% |\n",
                key[1], key[2], tw_min, tw_max, max_R, max_beta,
            )
        end
        missing_models = filter(row -> row.status == "model_critical_missing", rows)
        if !isempty(missing_models)
            println(io)
            println(io, "## Missing model critical points")
            println(io)
            for row in missing_models
                @printf(io, "- Tw=%.2f, %s, `%s`\n", row.Tw, row.mode, row.configuration)
            end
        end
        println(io)
        println(io, "The Lopez source curves predate the latest removal of the")
        println(io, "O(R^-2) axial thermal-feedback matrix entry. Recompute those")
        println(io, "curves before using final numerical values in a publication.")
    end
end

function main()
    isdir(INPUT_DIR) || error("Neutral-curve directory does not exist: $INPUT_DIR")
    paths = sort(filter(readdir(INPUT_DIR; join=true)) do path
        name = basename(path)
        isfile(path) && startswith(name, "ome=") && endswith(name, ".dat")
    end)
    isempty(paths) && error("No integrated per-case neutral curves found in $INPUT_DIR")

    cases = NamedTuple[]
    for path in paths
        metadata = parse_case(path)
        data = load_curve(path)
        critical = select_critical_points(data)
        push!(cases, merge(metadata, (path=path, critical=critical)))
    end

    references = Dict{Tuple{Int, String}, NamedTuple}()
    for case in cases
        case.configuration == "lopez" || continue
        tw_key = round(Int, 100 * case.Tw)
        for (mode, point) in case.critical
            references[(tw_key, mode)] = point
        end
    end

    rows = NamedTuple[]
    full_typeII_keys = Set(
        round(Int, 100 * case.Tw) for case in cases
        if case.configuration == "compressible_full" &&
           haskey(case.critical, "Type-II")
    )
    skipped_typeII_temperatures = Float64[]
    temperatures = sort(unique(getfield.(cases, :Tw)))
    for Tw in temperatures
        tw_key = round(Int, 100 * Tw)
        tw_key in full_typeII_keys || push!(skipped_typeII_temperatures, Tw)
    end

    ordered_cases = sort(cases; by=case -> (
        case.Tw, CONFIG_ORDER[case.configuration],
    ))
    for case in ordered_cases
        tw_key = round(Int, 100 * case.Tw)
        for mode in ("Type-I", "Type-II")
            mode == "Type-II" && !(tw_key in full_typeII_keys) && continue
            reference = get(references, (tw_key, mode), nothing)
            reference === nothing && continue
            point = get(case.critical, mode, nothing)
            if point === nothing
                push!(rows, (
                    Tw=case.Tw, mode=mode, configuration=case.configuration,
                    status="model_critical_missing", R_c=missing, beta_c=missing,
                    R_c_discrete=missing, beta_c_discrete=missing,
                    quadratic_fit=missing, R_c_reference=reference.R,
                    beta_c_reference=reference.beta, delta_R_c=missing,
                    R_c_error_percent=missing, abs_R_c_error_percent=missing,
                    delta_beta_c=missing, beta_c_error_percent=missing,
                    abs_beta_c_error_percent=missing, source_file=case.file,
                ))
                continue
            end
            delta_R = point.R - reference.R
            delta_beta = point.beta - reference.beta
            R_error = signed_percent(point.R, reference.R)
            beta_error = signed_percent(point.beta, reference.beta)
            push!(rows, (
                Tw=case.Tw, mode=mode, configuration=case.configuration,
                status="available", R_c=point.R, beta_c=point.beta,
                R_c_discrete=point.R_discrete,
                beta_c_discrete=point.beta_discrete,
                quadratic_fit=point.fitted, R_c_reference=reference.R,
                beta_c_reference=reference.beta, delta_R_c=delta_R,
                R_c_error_percent=R_error,
                abs_R_c_error_percent=abs(R_error), delta_beta_c=delta_beta,
                beta_c_error_percent=beta_error,
                abs_beta_c_error_percent=abs(beta_error), source_file=case.file,
            ))
        end
    end

    mkpath(OUTPUT_DIR)
    tsv_path = joinpath(OUTPUT_DIR, "neutral_critical_point_errors.tsv")
    dat_path = joinpath(OUTPUT_DIR, "neutral_critical_point_errors.dat")
    report_path = joinpath(OUTPUT_DIR, "neutral_critical_point_errors.md")
    write_tsv(tsv_path, rows)
    write_dat(dat_path, rows)
    write_report(report_path, rows, skipped_typeII_temperatures)

    available = count(row -> row.status == "available", rows)
    missing_count = count(row -> row.status == "model_critical_missing", rows)
    println("Processed $(length(cases)) integrated neutral curves")
    println("Available comparisons: $available; missing model minima: $missing_count")
    println("Type-II omitted at Tw: $(skipped_typeII_temperatures)")
    println("TSV: $tsv_path")
    println("Tecplot DAT: $dat_path")
    println("Report: $report_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
