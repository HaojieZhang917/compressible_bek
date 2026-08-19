using Printf

include(joinpath(@__DIR__, "zarf_spatial_growth_scan.jl"))

const RESULT_PATH = joinpath(
    get(ENV, "ZARF_RESULT_DIRECTORY", joinpath(
        @__DIR__,
        "zarf_spatial_growth_results",
        "Ts_p0p0_N99",
    )),
    "spatial_growth_samples.dat",
)
const OUTPUT_PATH = joinpath(
    get(ENV, "ZARF_RESULT_DIRECTORY", joinpath(
        @__DIR__,
        "zarf_spatial_growth_results",
        "Ts_p0p0_N99",
    )),
    "representative_high_accuracy_check.dat",
)

function load_result_rows(path)
    rows = NamedTuple[]
    for line in eachline(path)
        values = split(strip(line))
        (length(values) == 11 || length(values) == 12) || continue
        parsed = try
            parse.(Float64, values)
        catch
            continue
        end
        push!(rows, (
            R=parsed[1],
            beta=parsed[2],
            omega_r=parsed[3],
            omega_bar=parsed[4],
            branch_id=length(parsed) == 12 ? parsed[5] : 0.0,
            alpha_t=length(parsed) == 12 ? parsed[6] : parsed[5],
            omega_i=length(parsed) == 12 ? parsed[7] : parsed[6],
            alpha_r_old=length(parsed) == 12 ? parsed[8] : parsed[7],
            alpha_i_old=length(parsed) == 12 ? parsed[9] : parsed[8],
            growth_old=length(parsed) == 12 ? parsed[10] : parsed[9],
            residual_old=length(parsed) == 12 ? parsed[11] : parsed[10],
        ))
    end
    return rows
end

function representative_indices(rows)
    selected = Set{Int}()
    selectors = (
        sortperm(rows; by=row -> row.omega_i),
        sortperm(rows; by=row -> -row.growth_old),
        sortperm(rows; by=row -> -abs(row.alpha_r_old - row.alpha_t)),
        sortperm(rows; by=row -> -row.residual_old),
    )
    for order in selectors
        foreach(index -> push!(selected, index), order[1:min(5, length(order))])
    end
    radial_order = sortperm(rows; by=row -> (row.R, row.beta, row.omega_r))
    for fraction in range(0.0, 1.0; length=12)
        position = 1 + round(Int, fraction * (length(radial_order) - 1))
        push!(selected, radial_order[position])
    end
    return sort(collect(selected); by=index -> (
        rows[index].R, rows[index].beta, rows[index].omega_r
    ))
end

function main_validation()
    rows = load_result_rows(RESULT_PATH)
    indices = representative_indices(rows)
    cfg = Config(
        mass_fluxes=[0.0],
        N_cheb=99,
        candidate_count=4,
        eig_tolerance=1e-13,
        max_iterations=700,
        imaginary_shift=-0.005,
    )

    u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(cfg.Re_h, -1, 0.0, 1)
    D, D2, z = CRC_BF.Cheb(cfg.N_cheb, 1)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, cfg.N_cheb, 1)
    cof_by_radius = Dict{Float64, Any}()
    checks = NamedTuple[]

    for (counter, index) in enumerate(indices)
        row = rows[index]
        cof = get!(cof_by_radius, row.R) do
            CRC_STA.Spatial_mode_BEK1(
                F, G .- 1, H, row.R, cfg.N_cheb, D, D2, cfg.Re_h
            )
        end
        sample = (
            R=row.R,
            alpha_t=row.alpha_t,
            beta=row.beta,
            omega_r=row.omega_r,
            omega_i=row.omega_i,
        )
        alpha, residual = solve_spatial(sample, cof, D, D2, cfg)
        push!(checks, merge(row, (
            alpha_r_new=real(alpha),
            alpha_i_new=imag(alpha),
            growth_new=-imag(alpha),
            residual_new=residual,
            delta_alpha=abs(alpha - complex(row.alpha_r_old, row.alpha_i_old)),
        )))
        @printf("validated %d/%d: R=%.1f beta=%+.3f delta_alpha=%.3e\n",
                counter, length(indices), row.R, row.beta, checks[end].delta_alpha)
        flush(stdout)
    end

    open(OUTPUT_PATH, "w") do io
        println(io,
            "# R beta omega alpha_t omega_i alpha_r_old alpha_i_old " *
            "alpha_r_new alpha_i_new growth_new residual_old residual_new delta_alpha"
        )
        for row in checks
            @printf(io,
                "%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                row.R, row.beta, row.omega_r, row.alpha_t, row.omega_i,
                row.alpha_r_old, row.alpha_i_old,
                row.alpha_r_new, row.alpha_i_new, row.growth_new,
                row.residual_old, row.residual_new, row.delta_alpha,
            )
        end
    end

    @printf("checks=%d max_delta_alpha=%.6e max_residual_new=%.6e\n",
            length(checks), maximum(row.delta_alpha for row in checks),
            maximum(row.residual_new for row in checks))
    println("output: ", OUTPUT_PATH)
end

main_validation()
