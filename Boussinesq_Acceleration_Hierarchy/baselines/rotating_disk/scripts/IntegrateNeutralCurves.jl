using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BATCH_DIR = joinpath(ROOT,"neutral_curve_batch")
const OUTPUT_DIR = joinpath(ROOT,"neutral_curve_integrated")
const VARIABLES_LINE =
    "Variables=\"omega\" \"R\" \"beta\" \"alpha_r_1\" " *
    "\"alpha_i_1\" \"alpha_r_2\" \"alpha_i_2\""

function load_curve(path::AbstractString)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        text = strip(line)
        (isempty(text) || startswith(lowercase(text),"variables") ||
         startswith(lowercase(text),"zone") || startswith(text,"#")) && continue
        fields = split(text)
        length(fields) == 7 || error("Expected seven columns in $path: $text")
        push!(rows,parse.(Float64,fields))
    end
    isempty(rows) && return zeros(Float64,0,7)
    return reduce(vcat,permutedims.(rows))
end

function source_quality(path::AbstractString,data::AbstractMatrix)
    points = size(data,1)
    points > 0 || return (
        quality="invalid",direction="none",max_residual=NaN,
    )
    beta_difference = diff(data[:,3])
    direction = all(beta_difference .> 0) ? "increasing" :
        all(beta_difference .< 0) ? "decreasing" : "mixed"
    max_residual = maximum(min.(abs.(data[:,5]),abs.(data[:,7])))
    is_branch = occursin("_branch=",basename(path))
    quality = if is_branch && points >= 10 && max_residual <= 1.0e-6 &&
                 direction != "mixed"
        "validated_branch"
    elseif points >= 50 && max_residual <= 1.0e-6 &&
           direction != "mixed" && maximum(data[:,3]) >= 0.08
        "validated_strict"
    else
        "invalid"
    end
    return (
        quality=quality,direction=direction,max_residual=max_residual,
    )
end

case_stem(name::AbstractString) = replace(name,r"_branch=[^.]+" => "")

function extract_temperature(name::AbstractString)
    matched = match(r"Tw=([0-9.]+)",name)
    matched === nothing && return NaN
    return parse(Float64,matched.captures[1])
end

function endpoint_R(data::AbstractMatrix,which::Symbol)
    index = which === :minimum ? argmin(data[:,3]) : argmax(data[:,3])
    return data[index,2]
end

function case_quality(names,records)
    data = reduce(vcat,[records[name].data for name in names])
    beta_min = minimum(data[:,3])
    beta_max = maximum(data[:,3])
    low_R = endpoint_R(data,:minimum)
    high_R = endpoint_R(data,:maximum)
    max_residual = maximum(record.quality.max_residual for record in
                           (records[name] for name in names))
    strict_sources = all(
        startswith(records[name].quality.quality,"validated_") for name in names
    )
    has_separate_typeII = any(occursin("_branch=typeII",name) for name in names)

    # Complete cases reach the outer R approximately 500 boundary on both sides.
    # A disconnected Type-II mode is accepted only when its explicit source file
    # is present; otherwise low-beta endpoint coverage must be in the main curve.
    has_typeII = has_separate_typeII || (beta_min <= 0.045 && low_R >= 430.0)
    has_typeI = beta_max >= 0.10 && high_R >= 430.0
    status = if has_typeI && has_typeII
        strict_sources ? "complete_strict" : "complete_invalid_source"
    elseif has_typeI
        "missing_typeII"
    elseif has_typeII
        "missing_typeI"
    else
        "missing_both"
    end
    return (
        status=status,has_typeI=has_typeI,has_typeII=has_typeII,
        separate_typeII=has_separate_typeII,zones=length(names),
        points=size(data,1),beta_min=beta_min,beta_max=beta_max,
        low_R=low_R,high_R=high_R,max_residual=max_residual,
    )
end

function ordered_zone_names(names)
    return sort(names; by=name -> occursin("_branch=typeII",name) ? 0 : 1)
end

function decreasing_beta(data::AbstractMatrix)
    size(data,1) <= 1 && return Matrix(data)
    return data[1,3] >= data[end,3] ? Matrix(data) :
        Matrix(reverse(data; dims=1))
end

function merged_case_data(names,records,quality)
    if !quality.separate_typeII
        return decreasing_beta(records[only(names)].data)
    end
    typeII_name = only(filter(name -> occursin("_branch=typeII",name),names))
    typeI_name = only(filter(name -> !occursin("_branch=",name),names))
    typeI = decreasing_beta(records[typeI_name].data)
    typeII = decreasing_beta(records[typeII_name].data)
    return vcat(typeI,typeII)
end

function write_zone(io::IO,title::AbstractString,data::AbstractMatrix)
    println(io,"Zone T=\"$title\", I=$(size(data,1)), F=POINT")
    for row in eachrow(data)
        println(io,join(row,'\t'))
    end
end

function write_complete_case(path,stem,names,records,quality)
    title_stem = chop(stem; tail=4)
    data = merged_case_data(names,records,quality)
    open(path,"w") do io
        println(io,VARIABLES_LINE)
        write_zone(io,"$title_stem complete-Type-I-Type-II",data)
    end
end

function raw_source_paths()
    return sort(filter(readdir(BATCH_DIR; join=true)) do path
        name = basename(path)
        isfile(path) && startswith(name,"ome=") && endswith(name,".dat") &&
            !endswith(name,"_allbranches.dat")
    end)
end

function main()
    source_paths = raw_source_paths()
    isempty(source_paths) && error("No neutral-curve source files found in $BATCH_DIR")

    records = Dict{String,NamedTuple}()
    for path in source_paths
        name = basename(path)
        data = load_curve(path)
        quality = source_quality(path,data)
        quality.quality == "invalid" && error(
            "Invalid source $name: direction=$(quality.direction), " *
            "residual=$(quality.max_residual)",
        )
        records[name] = (path=path,data=data,quality=quality)
    end

    case_groups = Dict{String,Vector{String}}()
    for name in keys(records)
        push!(get!(case_groups,case_stem(name),String[]),name)
    end
    ordered_cases = sort(
        collect(keys(case_groups)); by=name -> (extract_temperature(name),name),
    )
    qualities = Dict(
        stem => case_quality(case_groups[stem],records) for stem in ordered_cases
    )
    incomplete = filter(
        stem -> qualities[stem].status != "complete_strict",ordered_cases,
    )
    isempty(incomplete) || error(
        "Refusing to publish incomplete cases: $(join(incomplete, ", "))",
    )

    mkpath(OUTPUT_DIR)
    for path in readdir(OUTPUT_DIR; join=true)
        name = basename(path)
        isfile(path) && startswith(name,"ome=") && endswith(name,".dat") &&
            rm(path; force=true)
    end

    for stem in ordered_cases
        write_complete_case(
            joinpath(OUTPUT_DIR,stem),stem,case_groups[stem],records,qualities[stem],
        )
    end

    combined_path = joinpath(OUTPUT_DIR,"neutral_curves_all.dat")
    open(combined_path,"w") do io
        println(io,VARIABLES_LINE)
        for stem in ordered_cases
            quality = qualities[stem]
            title_stem = chop(stem; tail=4)
            data = merged_case_data(case_groups[stem],records,quality)
            write_zone(io,"$title_stem complete-Type-I-Type-II",data)
        end
    end
    complete_path = joinpath(OUTPUT_DIR,"neutral_curves_complete.dat")
    cp(combined_path,complete_path; force=true)

    manifest_path = joinpath(OUTPUT_DIR,"manifest.tsv")
    open(manifest_path,"w") do io
        println(
            io,
            "file\tTw\tstatus\toutput_zones\tsource_segments\tpoints\t" *
            "beta_min\tbeta_max\t" *
            "R_at_beta_min\tR_at_beta_max\tmax_neutral_residual\tsources",
        )
        for stem in ordered_cases
            quality = qualities[stem]
            fields = (
                stem,extract_temperature(stem),quality.status,1,quality.zones,
                quality.points,quality.beta_min,quality.beta_max,
                quality.low_R,quality.high_R,quality.max_residual,
                join(ordered_zone_names(case_groups[stem]),','),
            )
            println(io,join(fields,'\t'))
        end
    end

    completeness_path = joinpath(OUTPUT_DIR,"completeness_manifest.tsv")
    open(completeness_path,"w") do io
        println(
            io,
            "case\tTw\tstatus\thas_typeI\thas_typeII\tseparate_typeII\t" *
            "output_zones\tpoints\tbeta_min\tbeta_max\tR_at_beta_min\t" *
            "R_at_beta_max\tmax_neutral_residual",
        )
        for stem in ordered_cases
            quality = qualities[stem]
            fields = (
                chop(stem; tail=4),extract_temperature(stem),quality.status,
                quality.has_typeI,quality.has_typeII,quality.separate_typeII,1,
                quality.points,quality.beta_min,quality.beta_max,
                quality.low_R,quality.high_R,quality.max_residual,
            )
            println(io,join(fields,'\t'))
        end
    end

    source_manifest_path = joinpath(OUTPUT_DIR,"source_manifest.tsv")
    open(source_manifest_path,"w") do io
        println(io,"source\tcase\tquality\tpoints\tmax_neutral_residual\tdirection")
        for name in sort(collect(keys(records)))
            record = records[name]
            fields = (
                relpath(record.path,ROOT),chop(case_stem(name); tail=4),
                record.quality.quality,size(record.data,1),
                record.quality.max_residual,record.quality.direction,
            )
            println(io,join(fields,'\t'))
        end
    end

    println("Published $(length(ordered_cases)) complete per-case files in $OUTPUT_DIR")
    println("No standalone Type-II files were copied to the integrated directory")
    println("Manifest: $manifest_path")
    println("Completeness manifest: $completeness_path")
    println("Tecplot multi-case file: $combined_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
