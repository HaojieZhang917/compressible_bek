using DelimitedFiles
using Printf

const VARIABLES_LINE =
    "Variables=\"omega\" \"R\" \"beta\" \"alpha_r_1\" \"alpha_i_1\" \"alpha_r_2\" \"alpha_i_2\""

function number_tag(value::Real; digits::Int=4)
    result = @sprintf("%.*f", digits, Float64(value))
    result = replace(result, r"0+$" => "")
    result = endswith(result, ".") ? result[begin:end-1] : result
    return occursin('.', result) ? result : result * ".0"
end

function lopez_stem(Tw::Real)
    return "ome=0.0_Tw=$(number_tag(Tw))_model=lopez"
end

function read_single_zone(path::AbstractString)
    isfile(path) || throw(ArgumentError("Curve file does not exist: $path"))
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        text = strip(line)
        (isempty(text) || startswith(lowercase(text), "variables") ||
         startswith(lowercase(text), "zone")) && continue
        fields = split(text)
        length(fields) == 7 || throw(ArgumentError(
            "Expected seven columns in $path, got: $text",
        ))
        push!(rows, parse.(Float64, fields))
    end
    isempty(rows) && return zeros(Float64, 0, 7)
    return reduce(vcat, permutedims.(rows))
end

function active_neutral_columns(data::AbstractMatrix)
    return [abs(data[index, 5]) <= abs(data[index, 7]) ? 1 : 2
            for index in axes(data, 1)]
end

function split_main_curve(main_data::AbstractMatrix)
    columns = active_neutral_columns(main_data)
    switches = findall(!=(0), diff(columns))
    isempty(switches) && return (typeI=Matrix(main_data), typeII=zeros(0, 7))
    length(switches) == 1 || error(
        "Expected at most one active-mode switch, found $(length(switches))",
    )
    split_index = only(switches)
    first_part = Matrix(main_data[begin:split_index, :])
    second_part = Matrix(main_data[split_index + 1:end, :])
    if columns[1] == 2
        return (typeI=second_part, typeII=first_part)
    end
    return (typeI=first_part, typeII=second_part)
end

function increasing_beta(data::AbstractMatrix)
    size(data, 1) <= 1 && return Matrix(data)
    return data[1, 3] <= data[end, 3] ? Matrix(data) : Matrix(reverse(data; dims=1))
end

function validate_zone(data::AbstractMatrix, name::AbstractString)
    size(data, 1) >= 2 || error("$name contains fewer than two points")
    all(isfinite, data) || error("$name contains non-finite values")
    all(diff(data[:, 3]) .> 0) || error("$name beta is not strictly increasing")
    residual = maximum(min.(abs.(data[:, 5]), abs.(data[:, 7])))
    residual <= 1.0e-6 || error("$name neutral residual is $residual")
    return residual
end

function write_zone(io::IO, title::AbstractString, data::AbstractMatrix)
    println(io, "Zone T=\"$title\", I=$(size(data, 1)), F=POINT")
    writedlm(io, data)
end

function branch_data(input_dir::AbstractString, Tw::Real)
    stem = lopez_stem(Tw)
    main_path = joinpath(input_dir, stem * ".dat")
    branch_path = joinpath(input_dir, stem * "_branch=typeII.dat")
    split = split_main_curve(read_single_zone(main_path))
    typeI = increasing_beta(split.typeI)
    typeII = if isfile(branch_path)
        increasing_beta(read_single_zone(branch_path))
    else
        increasing_beta(split.typeII)
    end
    size(typeII, 1) > 0 || error("No Type-II data found for Tw=$Tw")
    residual_I = validate_zone(typeI, "Tw=$Tw Type-I")
    residual_II = validate_zone(typeII, "Tw=$Tw Type-II")
    return (
        Tw=Float64(Tw), stem=stem, typeI=typeI, typeII=typeII,
        residual_I=residual_I, residual_II=residual_II,
    )
end

function merge_temperature(
    input_dir::AbstractString, Tw::Real;
    output_dir::AbstractString=input_dir,
)
    curves = branch_data(input_dir, Tw)
    mkpath(output_dir)
    destination = joinpath(output_dir, curves.stem * "_allbranches.dat")
    open(destination, "w") do io
        println(io, VARIABLES_LINE)
        write_zone(io, "Tw=$(number_tag(Tw)) Type-II", curves.typeII)
        write_zone(io, "Tw=$(number_tag(Tw)) Type-I", curves.typeI)
    end
    @printf(
        "Tw=%.2f merged: Type-II=%d Type-I=%d residuals=[%.3e, %.3e]\n",
        Tw, size(curves.typeII, 1), size(curves.typeI, 1),
        curves.residual_II, curves.residual_I,
    )
    return merge(curves, (path=destination,))
end

function merge_all(
    temperatures=(1.08, 1.10, 1.12, 1.14, 1.16, 1.18, 1.20);
    input_dir::AbstractString=joinpath(@__DIR__, "..", "neutral_curve_batch"),
    output_dir::AbstractString=input_dir,
)
    results = [merge_temperature(input_dir, Tw; output_dir=output_dir)
               for Tw in temperatures]
    first_tag = number_tag(minimum(temperatures))
    last_tag = number_tag(maximum(temperatures))
    master_path = joinpath(
        output_dir, "lopez_allbranches_Tw$(first_tag)-$(last_tag).dat",
    )
    open(master_path, "w") do io
        println(io, VARIABLES_LINE)
        for result in results
            tag = number_tag(result.Tw)
            write_zone(io, "Tw=$tag Type-II", result.typeII)
            write_zone(io, "Tw=$tag Type-I", result.typeI)
        end
    end
    println("saved combined multi-temperature file $master_path")
    return (results=results, master_path=master_path)
end

function main(args=ARGS)
    temperatures = isempty(args) ?
        (1.08, 1.10, 1.12, 1.14, 1.16, 1.18, 1.20) :
        Tuple(parse.(Float64, args))
    return merge_all(temperatures)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
