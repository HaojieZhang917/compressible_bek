using Printf

const EXTENDED_TEMPERATURES = (1.4, 1.6, 1.8)

function number_tag(value::Real)
    text = replace(@sprintf("%.4f", value), r"0+$" => "")
    return endswith(text, ".") ? text * "0" : text
end

function paths_for_temperature(Tw::Real)
    tag = number_tag(Tw)
    return (
        joinpath(
            @__DIR__, "blackburn_neutral_curve_batch",
            "ome=0.0_Tw=$(tag)_model=blackburn.dat",
        ),
        joinpath(
            @__DIR__, "neutral_curve_batch",
            "ome=0.0_Tw=$(tag)_model=compressible_Mr=0.3_" *
            "propPert=on_baseProp=variable.dat",
        ),
    )
end

function normalize_file(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 3 || error("curve has no data rows: $path")
    rows = lines[3:end]
    fields = split.(strip.(rows))
    all(length(row) == 7 for row in fields) || error(
        "curve does not have seven columns: $path",
    )
    beta = parse.(Float64, getindex.(fields, 3))
    differences = diff(beta)
    if all(differences .< 0)
        rows = reverse(rows)
        beta = reverse(beta)
        open(path, "w") do io
            println(io, lines[1])
            println(io, lines[2])
            foreach(line -> println(io, line), rows)
        end
        action = "reversed"
    elseif all(differences .> 0)
        action = "unchanged"
    else
        error("beta is not strictly monotone: $path")
    end
    all(diff(beta) .> 0) || error("beta normalization failed: $path")
    @printf("%s points=%d beta=[%.9f, %.9f] %s\n",
        basename(path), length(beta), first(beta), last(beta), action)
end

for Tw in EXTENDED_TEMPERATURES, path in paths_for_temperature(Tw)
    normalize_file(path)
end
