using Printf
using Dates
using PyCall

const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(WORKSPACE_ROOT, "CRD_STA.jl"))

const DEFAULT_TW_VALUES = "1.0:0.1:2.0"

function parse_float_env(name::String, default::Float64)
    return parse(Float64, get(ENV, name, string(default)))
end

function parse_int_env(name::String, default::Int)
    return parse(Int, get(ENV, name, string(default)))
end

function float_range(a::Float64, step::Float64, b::Float64)
    if step == 0.0
        error("Tw step cannot be zero")
    end
    vals = Float64[]
    x = a
    if step > 0
        while x <= b + 0.5 * abs(step)
            push!(vals, round(x; digits=12))
            x += step
        end
    else
        while x >= b - 0.5 * abs(step)
            push!(vals, round(x; digits=12))
            x += step
        end
    end
    return vals
end

function parse_tw_values(spec::String)
    s = strip(spec)
    if occursin(":", s)
        parts = parse.(Float64, split(s, ":"))
        if length(parts) == 2
            return float_range(parts[1], 0.1, parts[2])
        elseif length(parts) == 3
            return float_range(parts[1], parts[2], parts[3])
        end
        error("TW_VALUES with ':' must be start:stop or start:step:stop")
    end
    return parse.(Float64, split(s, ","))
end

function fmt(x)
    if x isa Integer
        return string(x)
    elseif x isa Real
        return @sprintf("%.16e", Float64(x))
    else
        return string(x)
    end
end

function write_header(path::String, header)
    open(path, "w") do io
        println(io, join(header, ","))
    end
end

function append_csv_row(io, row)
    println(io, join(fmt.(row), ","))
end

function vec64(x)
    return vec(Float64.(x))
end

function unique_sorted_xy(xsrc, ysrc)
    x = vec64(xsrc)
    y = vec64(ysrc)
    p = sortperm(x)
    x = x[p]
    y = y[p]
    xu = Float64[]
    yu = Float64[]
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[j + 1] == x[i]
            j += 1
        end
        push!(xu, x[j])
        push!(yu, y[j])
        i = j + 1
    end
    return xu, yu
end

function linear_interp(xsrc, ysrc, xdst)
    x, y = unique_sorted_xy(xsrc, ysrc)
    out = similar(vec64(xdst))
    n = length(x)
    for (i, xi) in pairs(vec64(xdst))
        if xi <= x[1]
            out[i] = y[1]
        elseif xi >= x[end]
            out[i] = y[end]
        else
            k = searchsortedlast(x, xi)
            t = (xi - x[k]) / (x[k + 1] - x[k])
            out[i] = (1 - t) * y[k] + t * y[k + 1]
        end
    end
    return out
end

function py_to_dict(info)
    d = Dict{String, Float64}()
    for key in ["Tw", "Pr", "lambda_c", "Hinf", "Fp0", "Gp0", "Tp0", "N", "zmax"]
        try
            d[key] = Float64(info[key])
        catch
        end
    end
    return d
end

function import_bone()
    sys = pyimport("sys")
    path_list = PyVector(sys."path")
    thisdir = String(WORKSPACE_ROOT)
    if !(thisdir in String.(path_list))
        pushfirst!(path_list, thisdir)
    end
    return pyimport("Bone")
end

function compressible_base_raw(Tw, Mr, gamma, u0, v0, w0, f, q)
    T_raw = vec64(T_ca(Mr, f, q, w0, gamma, Tw)[2])
    H_raw = vec64(w0) .* T_raw
    rho_raw = 1.0 ./ T_raw
    z_raw = CRD_BF.Physical_Interpretation(T_raw, 0.004, 10001)
    return vec64(u0), vec64(v0), H_raw, T_raw, rho_raw, vec64(z_raw)
end

function compressible_base_cheb(Tw, Mr, gamma, N_cheb, x_cheb, u0, v0, w0, f, q)
    H_raw, T_raw = T_ca(Mr, f, q, w0, gamma, Tw)
    F, G, H, T, rho, z_phys = interp(u0, v0, H_raw, T_raw, x_cheb, N_cheb, "phy")
    return vec64(F), vec64(G), vec64(H), vec64(T), vec64(rho), vec64(z_phys)
end

function maxabs(x)
    return maximum(abs.(vec64(x)))
end

function main()
    tw_values = parse_tw_values(get(ENV, "TW_VALUES", DEFAULT_TW_VALUES))
    n_cheb = parse_int_env("N_CHEB", 199)
    n_common = parse_int_env("N_COMMON", 2001)
    z_common_max = parse_float_env("Z_COMMON_MAX", 20.0)
    ro = parse_float_env("RO", -1.0)
    mr = parse_float_env("MR", 0.3)
    gamma = parse_float_env("GAMMA", 1.4)
    sigma = parse_float_env("SIGMA", 0.72)
    out_dir = get(ENV, "OUT_DIR", joinpath(WORKSPACE_ROOT, "baseflow_comparison_data"))

    mkpath(out_dir)

    raw_inc_path = joinpath(out_dir, "incompressible_boussinesq_physical_raw.csv")
    raw_comp_path = joinpath(out_dir, "compressible_physical_raw.csv")
    common_path = joinpath(out_dir, "physical_common_grid_interpolated.csv")
    cheb_path = joinpath(out_dir, "cheb_grid_interpolated.csv")
    summary_path = joinpath(out_dir, "summary.csv")
    metadata_path = joinpath(out_dir, "metadata.txt")

    write_header(raw_inc_path, [
        "Tw", "z", "F_i", "G_i", "minus_G_i", "H_i", "T_i", "dF_i", "dG_i", "dT_i",
        "Hinf_i", "Fp0_i", "Gp0_i", "Tp0_i"
    ])
    write_header(raw_comp_path, ["Tw", "z_physical", "F_c", "G_c", "H_c", "T_c", "rho_c"])
    write_header(common_path, [
        "Tw", "z",
        "F_i", "G_i", "minus_G_i", "H_i", "T_i",
        "F_c", "G_c", "H_c", "T_c", "rho_c",
        "dF_c_minus_i", "dG_c_minus_i_raw", "dG_c_minus_minus_Gi", "dH_c_minus_i", "dT_c_minus_i"
    ])
    write_header(cheb_path, [
        "Tw", "x_cheb",
        "F_i", "G_i", "minus_G_i", "H_i", "T_i",
        "F_c", "G_c", "H_c", "T_c", "rho_c",
        "dF_c_minus_i", "dG_c_minus_i_raw", "dG_c_minus_minus_Gi", "dH_c_minus_i", "dT_c_minus_i"
    ])
    write_header(summary_path, [
        "Tw", "Ro", "Mr", "N_cheb", "N_common", "z_common_max",
        "Hinf_i", "Fp0_i", "Gp0_i", "Tp0_i",
        "Hinf_c_at_zmax", "F_c_at_zmax", "G_c_at_zmax", "T_c_at_zmax",
        "maxabs_dF_common", "maxabs_dG_raw_common", "maxabs_dG_signmatched_common",
        "maxabs_dH_common", "maxabs_dT_common",
        "maxabs_dF_cheb", "maxabs_dG_raw_cheb", "maxabs_dG_signmatched_cheb",
        "maxabs_dH_cheb", "maxabs_dT_cheb"
    ])

    open(metadata_path, "w") do io
        println(io, "Generated at: $(Dates.now())")
        println(io, "Script: $(joinpath(@__DIR__, "generate_baseflow_comparison.jl"))")
        println(io, "Tw values: $(join(tw_values, ", "))")
        println(io, "Ro: $ro")
        println(io, "Mr: $mr")
        println(io, "gamma: $gamma")
        println(io, "sigma: $sigma")
        println(io, "N_cheb: $n_cheb")
        println(io, "N_common: $n_common")
        println(io, "Z_COMMON_MAX: $z_common_max")
        println(io, "Compressible path follows Bone.ipynb: baseflow_var -> T_ca -> interp(..., \"phy\").")
        println(io, "Incompressible path uses Bone.py get_baseflow(Tw), with ideal-gas Boussinesq lambda_c = 1.")
        println(io, "The column minus_G_i is included because the notebook passes -G into the incompressible stability matrix.")
    end

    println("Preparing compressible base-flow operators: N_cheb=$n_cheb, Ro=$ro")
    co = 2 - ro - ro^2
    u0, v0, w0, f, q, D, D2, x_cheb_mat = baseflow_var(n_cheb, ro, co)
    x_cheb = vec64(x_cheb_mat)
    z_common = collect(range(0.0, z_common_max; length=n_common))
    bone = import_bone()

    open(raw_inc_path, "a") do raw_inc_io
        open(raw_comp_path, "a") do raw_comp_io
            open(common_path, "a") do common_io
                open(cheb_path, "a") do cheb_io
                    open(summary_path, "a") do summary_io
                        for Tw in tw_values
                            println("Computing Tw=$(Tw)")

                            z_i, H_i, F_i, G_i, T_i, dF_i, dG_i, dT_i, info_py = bone.get_baseflow(Tw)
                            z_i = vec64(z_i)
                            H_i = vec64(H_i)
                            F_i = vec64(F_i)
                            G_i = vec64(G_i)
                            T_i = vec64(T_i)
                            dF_i = vec64(dF_i)
                            dG_i = vec64(dG_i)
                            dT_i = vec64(dT_i)
                            info = py_to_dict(info_py)

                            for j in eachindex(z_i)
                                append_csv_row(raw_inc_io, [
                                    Tw, z_i[j], F_i[j], G_i[j], -G_i[j], H_i[j], T_i[j],
                                    dF_i[j], dG_i[j], dT_i[j],
                                    get(info, "Hinf", NaN), get(info, "Fp0", NaN),
                                    get(info, "Gp0", NaN), get(info, "Tp0", NaN)
                                ])
                            end

                            F_c_raw, G_c_raw, H_c_raw, T_c_raw, rho_c_raw, z_c_raw =
                                compressible_base_raw(Tw, mr, gamma, u0, v0, w0, f, q)
                            for j in eachindex(z_c_raw)
                                append_csv_row(raw_comp_io, [
                                    Tw, z_c_raw[j], F_c_raw[j], G_c_raw[j],
                                    H_c_raw[j], T_c_raw[j], rho_c_raw[j]
                                ])
                            end

                            F_i_common = linear_interp(z_i, F_i, z_common)
                            G_i_common = linear_interp(z_i, G_i, z_common)
                            H_i_common = linear_interp(z_i, H_i, z_common)
                            T_i_common = linear_interp(z_i, T_i, z_common)
                            F_c_common = linear_interp(z_c_raw, F_c_raw, z_common)
                            G_c_common = linear_interp(z_c_raw, G_c_raw, z_common)
                            H_c_common = linear_interp(z_c_raw, H_c_raw, z_common)
                            T_c_common = linear_interp(z_c_raw, T_c_raw, z_common)
                            rho_c_common = linear_interp(z_c_raw, rho_c_raw, z_common)

                            for j in eachindex(z_common)
                                append_csv_row(common_io, [
                                    Tw, z_common[j],
                                    F_i_common[j], G_i_common[j], -G_i_common[j], H_i_common[j], T_i_common[j],
                                    F_c_common[j], G_c_common[j], H_c_common[j], T_c_common[j], rho_c_common[j],
                                    F_c_common[j] - F_i_common[j],
                                    G_c_common[j] - G_i_common[j],
                                    G_c_common[j] + G_i_common[j],
                                    H_c_common[j] - H_i_common[j],
                                    T_c_common[j] - T_i_common[j]
                                ])
                            end

                            F_c_cheb, G_c_cheb, H_c_cheb, T_c_cheb, rho_c_cheb, z_c_phys =
                                compressible_base_cheb(Tw, mr, gamma, n_cheb, x_cheb_mat, u0, v0, w0, f, q)
                            F_i_cheb = linear_interp(z_i, F_i, x_cheb)
                            G_i_cheb = linear_interp(z_i, G_i, x_cheb)
                            H_i_cheb = linear_interp(z_i, H_i, x_cheb)
                            T_i_cheb = linear_interp(z_i, T_i, x_cheb)

                            for j in eachindex(x_cheb)
                                append_csv_row(cheb_io, [
                                    Tw, x_cheb[j],
                                    F_i_cheb[j], G_i_cheb[j], -G_i_cheb[j], H_i_cheb[j], T_i_cheb[j],
                                    F_c_cheb[j], G_c_cheb[j], H_c_cheb[j], T_c_cheb[j], rho_c_cheb[j],
                                    F_c_cheb[j] - F_i_cheb[j],
                                    G_c_cheb[j] - G_i_cheb[j],
                                    G_c_cheb[j] + G_i_cheb[j],
                                    H_c_cheb[j] - H_i_cheb[j],
                                    T_c_cheb[j] - T_i_cheb[j]
                                ])
                            end

                            append_csv_row(summary_io, [
                                Tw, ro, mr, n_cheb, n_common, z_common_max,
                                get(info, "Hinf", NaN), get(info, "Fp0", NaN),
                                get(info, "Gp0", NaN), get(info, "Tp0", NaN),
                                H_c_common[end], F_c_common[end], G_c_common[end], T_c_common[end],
                                maxabs(F_c_common - F_i_common),
                                maxabs(G_c_common - G_i_common),
                                maxabs(G_c_common + G_i_common),
                                maxabs(H_c_common - H_i_common),
                                maxabs(T_c_common - T_i_common),
                                maxabs(F_c_cheb - F_i_cheb),
                                maxabs(G_c_cheb - G_i_cheb),
                                maxabs(G_c_cheb + G_i_cheb),
                                maxabs(H_c_cheb - H_i_cheb),
                                maxabs(T_c_cheb - T_i_cheb)
                            ])
                        end
                    end
                end
            end
        end
    end

    println("Done. Data written to: $out_dir")
end

main()
