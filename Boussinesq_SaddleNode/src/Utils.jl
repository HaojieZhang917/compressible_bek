const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

project_root() = PROJECT_ROOT

"""Parse GNU-style `--name value` arguments using typed defaults."""
function parse_cli(arguments::Vector{String}, defaults::Dict{String,Any};
                   choices::Dict{String,Vector{String}}=Dict{String,Vector{String}}())
    values = copy(defaults)
    index = 1
    while index <= length(arguments)
        token = arguments[index]
        startswith(token, "--") || error("Unexpected positional argument: $token")
        raw = token[3:end]
        if occursin('=', raw)
            name, text = split(raw, '='; limit=2)
        else
            name = raw
            haskey(defaults, name) || error("Unknown option --$name")
            if defaults[name] isa Bool
                text = "true"
            else
                index == length(arguments) && error("Missing value after --$name")
                index += 1
                text = arguments[index]
            end
        end
        haskey(defaults, name) || error("Unknown option --$name")
        default = defaults[name]
        value = if default isa Bool
            lowercase(text) in ("1", "true", "yes", "on")
        elseif default isa Int
            parse(Int, text)
        elseif default isa AbstractFloat
            parse(Float64, text)
        elseif default isa Vector{Int}
            parse.(Int, split(text, ','))
        elseif default isa Vector{Float64}
            parse.(Float64, split(text, ','))
        else
            text
        end
        if haskey(choices, name) && !(string(value) in choices[name])
            error("Invalid --$name=$(value); choose one of $(join(choices[name], ", "))")
        end
        values[name] = value
        index += 1
    end
    return values
end

function read_csv_rows(path::AbstractString)
    rows = Vector{Dict{String,Any}}()
    for row in CSV.File(path; normalizenames=false)
        item = Dict{String,Any}()
        for name in propertynames(row)
            item[string(name)] = getproperty(row, name)
        end
        push!(rows, item)
    end
    return rows
end

function write_csv_rows(path::AbstractString, rows::Vector{<:AbstractDict};
                        columns::Union{Nothing,Vector{String}}=nothing)
    mkpath(dirname(path))
    if isempty(rows)
        open(path, "w") do io
            columns === nothing || println(io, join(columns, ','))
        end
        return path
    end
    names = columns === nothing ? string.(collect(keys(first(rows)))) : columns
    table = (; (Symbol(name) => [get(row, name, missing) for row in rows] for name in names)...)
    CSV.write(path, table)
    return path
end

function write_csv_rows(path::AbstractString, columns::AbstractDict)
    names = string.(collect(keys(columns)))
    rows = [Dict(name => columns[name][i] for name in names)
            for i in eachindex(columns[first(names)])]
    return write_csv_rows(path, rows; columns=names)
end

function write_json(path::AbstractString, value)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.pretty(io, value; allow_inf=true)
        println(io)
    end
    return path
end

function chebyshev_operators(degree::Int, left::Float64=0.0, right::Float64=1.0)
    degree >= 1 || error("Chebyshev degree must be positive")
    j = collect(0:degree)
    x = cos.(pi .* j ./ degree)
    c = vcat(2.0, ones(degree - 1), 2.0) .* ((-1.0) .^ j)
    X = repeat(x, 1, degree + 1)
    dX = X .- X'
    D = (c * (1.0 ./ c)') ./ (dX + Matrix{Float64}(I, degree + 1, degree + 1))
    D .-= Diagonal(vec(sum(D; dims=2)))
    z = left .+ 0.5 * (right - left) .* (1.0 .- x)
    Dz = -2.0 .* D ./ (right - left)
    weights = (-1.0) .^ j
    weights[1] *= 0.5
    weights[end] *= 0.5
    return z, Dz, Dz * Dz, weights
end

function barycentric_interpolate(nodes::AbstractVector, weights::AbstractVector,
                                 values::AbstractVector, targets::AbstractVector)
    result = Vector{promote_type(eltype(values), Float64)}(undef, length(targets))
    scale = max(maximum(abs, nodes), 1.0)
    for (i, target) in pairs(targets)
        exact = findfirst(x -> abs(x - target) <= 8eps(Float64) * scale, nodes)
        if exact !== nothing
            result[i] = values[exact]
        else
            terms = weights ./ (target .- nodes)
            result[i] = dot(terms, values) / sum(terms)
        end
    end
    return result
end

function barycentric_interpolate(nodes::AbstractVector, weights::AbstractVector,
                                 values::AbstractMatrix, targets::AbstractVector)
    output = Matrix{promote_type(eltype(values), Float64)}(undef, size(values, 1), length(targets))
    for row in axes(values, 1)
        output[row, :] .= barycentric_interpolate(nodes, weights, view(values, row, :), targets)
    end
    return output
end

function linear_interpolate(x::AbstractVector, y::AbstractVector, xq::AbstractVector)
    order = sortperm(x)
    xs = x[order]
    ys = y[order]
    output = similar(collect(xq), promote_type(eltype(y), Float64))
    for (k, q) in pairs(xq)
        if q <= xs[1]
            output[k] = ys[1]
        elseif q >= xs[end]
            output[k] = ys[end]
        else
            hi = searchsortedfirst(xs, q)
            lo = hi - 1
            alpha = (q - xs[lo]) / (xs[hi] - xs[lo])
            output[k] = (1 - alpha) * ys[lo] + alpha * ys[hi]
        end
    end
    return output
end

"""Locate smooth local extrema with a three-point nonuniform quadratic fit."""
function locate_quadratic_extrema(x::AbstractVector, y::AbstractVector)
    order = sortperm(x)
    xs, ys = collect(x[order]), collect(y[order])
    output = NamedTuple[]
    slopes = diff(ys) ./ diff(xs)
    for i in 2:length(xs)-1
        slopes[i-1] * slopes[i] <= 0 || continue
        A = hcat(ones(3), xs[i-1:i+1], xs[i-1:i+1].^2)
        a, b, c = A \ ys[i-1:i+1]
        abs(c) > eps(Float64) || continue
        root = -b / (2c)
        xs[i-1] < root < xs[i+1] || continue
        value = a + b * root + c * root^2
        kind = c > 0 ? "minimum" : "maximum"
        push!(output, (kind=kind, x=root, y=value, curvature=2c))
    end
    sort!(output; by=item -> item.x)
    unique_output = NamedTuple[]
    for item in output
        if isempty(unique_output) || abs(item.x - last(unique_output).x) > 1e-10
            push!(unique_output, item)
        end
    end
    return unique_output
end

function bracketed_roots(x::AbstractVector, y::AbstractVector, target::Real)
    order = sortperm(x)
    xs, values = collect(x[order]), collect(y[order] .- target)
    roots = Float64[]
    for i in 1:length(xs)-1
        if values[i] == 0
            push!(roots, xs[i])
        elseif values[i] * values[i+1] < 0
            alpha = -values[i] / (values[i+1] - values[i])
            push!(roots, xs[i] + alpha * (xs[i+1] - xs[i]))
        end
    end
    return unique(round.(roots; digits=12))
end

function simpson_uniform(values::AbstractVector, coordinates::AbstractVector)
    intervals = length(values) - 1
    iseven(intervals) || error("Simpson integration requires an even interval count")
    h = coordinates[2] - coordinates[1]
    return h * (values[1] + values[end] +
                4sum(values[2:2:end-1]) + 2sum(values[3:2:end-2])) / 3
end

clean_name(value) = begin
    text = replace(strip(string(value)), r"\s+" => "_")
    text = replace(text, r"[^0-9A-Za-z_]+" => "_")
    text = strip(text, '_')
    isempty(text) && return "value"
    isdigit(first(text)) ? "v_" * text : text
end

quote_tecplot(value) = replace(replace(string(value), '\\' => '/'), '"' => "\\\"")

function numeric_value(value)
    value === missing && return NaN
    value === nothing && return NaN
    value isa Bool && return value ? 1.0 : 0.0
    value isa Number && return Float64(value)
    text = lowercase(strip(string(value)))
    isempty(text) && return NaN
    text == "true" && return 1.0
    text == "false" && return 0.0
    return parse(Float64, text)
end

format_number(value) = begin
    number = numeric_value(value)
    isnan(number) && return "NaN"
    isinf(number) && return number > 0 ? "Inf" : "-Inf"
    @sprintf("%.16g", number)
end
