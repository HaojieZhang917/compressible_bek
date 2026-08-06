#!/usr/bin/env julia

"""Rational-Chebyshev continuation of the heated von Karman branch.

The compactification eta=L(1-x)/(1+x) places the disk at x=1 and infinity
at x=-1. H_infinity is the continuation coordinate, so temperature folds do
not make the nonlinear solve singular.
"""

include(joinpath(@__DIR__, "..", "..", "src", "BoussinesqSaddleNode.jl"))
using .BoussinesqSaddleNode
using JSON3
using LinearAlgebra
using Plots
using Printf

const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "data", "infinite_mapping")
const PR = 0.72

struct MappedOperators
    degree::Int
    scale::Float64
    x::Vector{Float64}
    eta::Vector{Float64}
    dx::Matrix{Float64}
    deta::Matrix{Float64}
    deta2::Matrix{Float64}
    weights::Vector{Float64}
end

struct MappedSolution
    hinf::Float64
    tw::Float64
    fields::Matrix{Float64}
    residual::Float64
    iterations::Int
end

function rational_chebyshev(degree::Int, scale::Real)
    degree >= 2 || error("Mapped Chebyshev degree must be at least 2")
    j = collect(0:degree)
    x = cos.(pi .* j ./ degree)
    c = vcat(2.0, ones(degree - 1), 2.0) .* ((-1.0) .^ j)
    X = repeat(x, 1, degree + 1)
    delta = X .- transpose(X)
    dx = (c * transpose(1.0 ./ c)) ./
         (delta + Matrix{Float64}(I, degree + 1, degree + 1))
    dx .-= Diagonal(vec(sum(dx; dims=2)))
    eta = [index == degree + 1 ? Inf :
           Float64(scale) * (1 - x[index]) / (1 + x[index])
           for index in eachindex(x)]
    dx_deta = -(1 .+ x).^2 ./ (2Float64(scale))
    deta = Diagonal(dx_deta) * dx
    weights = (-1.0) .^ j
    weights[1] *= 0.5
    weights[end] *= 0.5
    return MappedOperators(degree, Float64(scale), x, eta, dx, deta,
                           deta * deta, weights)
end

function mapped_residual_jacobian(state::AbstractVector, h_target::Real,
                                  operators::MappedOperators)
    nodes = length(operators.x)
    field_state = view(state, 1:length(state)-1)
    tw = Float64(state[end])
    residual, base_jacobian = BoussinesqSaddleNode._vk_fixed_tw_system(
        field_state, tw, operators.deta, operators.deta2)
    h = view(field_state, 1:nodes)
    far_continuity_row = nodes
    residual[far_continuity_row] = dot(view(operators.dx, nodes, :), h)
    base_jacobian[far_continuity_row, :] .= 0
    base_jacobian[far_continuity_row, 1:nodes] .= operators.dx[nodes, :]

    augmented_residual = vcat(residual, h[end] - h_target)
    jacobian = zeros(4nodes + 1, 4nodes + 1)
    jacobian[1:end-1, 1:end-1] .= base_jacobian
    jacobian[3nodes + 1, end] = -1
    jacobian[end, nodes] = 1
    return augmented_residual, jacobian
end

function mapped_newton(h_target::Real, initial_state::AbstractVector,
                       operators::MappedOperators; tolerance::Float64=2e-10,
                       max_iterations::Int=18)
    state = Float64.(initial_state)
    history = Float64[]
    for _ in 1:max_iterations
        residual, jacobian = mapped_residual_jacobian(state, h_target, operators)
        norm0 = norm(residual, Inf)
        push!(history, norm0)
        if norm0 < tolerance
            nodes = length(operators.x)
            return MappedSolution(Float64(h_target), state[end],
                permutedims(reshape(view(state, 1:length(state)-1), nodes, 4)),
                norm0, length(history))
        end
        row_norm = [norm(view(jacobian, row, :)) for row in axes(jacobian, 1)]
        scaling = 1.0 ./ max.(row_norm, 1e-300)
        step = -((scaling .* jacobian) \ (scaling .* residual))
        damping = 1.0
        accepted = false
        while damping >= 2e-6
            trial = state .+ damping .* step
            trial_norm = norm(first(mapped_residual_jacobian(trial, h_target, operators)), Inf)
            if isfinite(trial_norm) && trial_norm < norm0
                state = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        accepted || error("Mapped Newton line search failed at Hinf=$h_target, residual=$norm0")
    end
    error("Mapped Newton did not converge at Hinf=$h_target; residual=$(last(history))")
end

function finite_domain_initial_state(h_start::Real, operators::MappedOperators)
    zmax = 60.0
    degree = max(140, operators.degree)
    seed = solve_vk_isothermal(; zmax, degree, tolerance=2e-9)
    source = solve_vk_fixed_h(h_start, seed; zmax, degree, tolerance=2e-9)
    sample_eta = min.(operators.eta, zmax)
    fields = barycentric_interpolate(source.z, source.weights, source.fields, sample_eta)
    far = operators.eta .>= zmax
    fields[1, far] .= h_start
    fields[2, far] .= 0
    fields[3, far] .= 1
    fields[4, far] .= 1
    return vcat(BoussinesqSaddleNode.pack_fields(fields), source.tw)
end

function default_h_values(h_stop::Real)
    pieces = (collect(-0.75:0.01:-0.60), collect(-0.60:0.002:-0.47),
              collect(-0.47:0.005:-0.16), collect(-0.16:0.001:-0.06),
              collect(-0.06:0.002:(Float64(h_stop) + 1e-10)))
    values = sort(unique(round.(vcat(pieces...); digits=12)))
    return filter(value -> -0.7500001 <= value <= h_stop + 1e-10, values)
end

function trace_mapped_branch(degree::Int, scale::Real, h_stop::Real,
                             tolerance::Real)
    operators = rational_chebyshev(degree, scale)
    h_values = default_h_values(h_stop)
    state = finite_domain_initial_state(first(h_values), operators)
    solutions = MappedSolution[]
    for (index, hinf) in enumerate(h_values)
        solution = mapped_newton(hinf, state, operators;
                                 tolerance=Float64(tolerance))
        push!(solutions, solution)
        state = vcat(BoussinesqSaddleNode.pack_fields(solution.fields), solution.tw)
        if index == 1 || index % 40 == 0
            @printf("N=%d L=%g: %d/%d, Hinf=%.4f, Tw=%.8f, res=%.2e\n",
                    degree, scale, index, length(h_values), hinf,
                    solution.tw, solution.residual)
            flush(stdout)
        end
    end
    return operators, solutions
end

function locate_folds(solutions::Vector{MappedSolution})
    h = [item.hinf for item in solutions]
    tw = [item.tw for item in solutions]
    extrema = locate_quadratic_extrema(h, tw)
    rows = Dict{String,Any}[]
    for (index, item) in enumerate(extrema)
        push!(rows, Dict("fold_index" => index, "kind" => item.kind,
                        "Hinf" => item.x, "Tw" => item.y,
                        "d2Tw_dHinf2" => item.curvature))
    end
    return rows
end

function interpolate_fields(operators::MappedOperators, solution::MappedSolution,
                            eta_query::AbstractVector)
    x_query = (operators.scale .- eta_query) ./ (operators.scale .+ eta_query)
    return barycentric_interpolate(operators.x, operators.weights,
                                   solution.fields, x_query)
end

function branch_rows(operators, solutions, folds)
    fold_h = [row["Hinf"] for row in folds]
    rows = Dict{String,Any}[]
    for solution in solutions
        h, f, g, temperature = (view(solution.fields, row, :) for row in 1:4)
        push!(rows, Dict(
            "degree" => operators.degree, "mapping_scale" => operators.scale,
            "Hinf" => solution.hinf, "Tw" => solution.tw,
            "Fp0" => dot(operators.deta[1, :], f),
            "Gp0" => dot(operators.deta[1, :], g),
            "Tp0" => dot(operators.deta[1, :], temperature),
            "thermal_tail_length" => 1 / (PR * abs(solution.hinf)),
            "residual_inf" => solution.residual,
            "newton_iterations" => solution.iterations,
            "nearest_fold_distance_Hinf" => isempty(fold_h) ? NaN :
                                             minimum(abs.(fold_h .- solution.hinf))))
    end
    return rows
end

function save_selected_profiles(operators, solutions, output)
    targets = (-0.60, -0.532, -0.30, -0.15, -0.10, -0.08, -0.04, -0.02)
    eta_plot = vcat(collect(range(0.0, 20.0; length=401)),
                    exp.(range(log(20.1), log(2000.0); length=500)))
    rows = Dict{String,Any}[]
    for target in targets
        solution = solutions[argmin(abs.([item.hinf for item in solutions] .- target))]
        values = interpolate_fields(operators, solution, eta_plot)
        for index in eachindex(eta_plot)
            push!(rows, Dict("Hinf" => solution.hinf, "Tw" => solution.tw,
                "eta" => eta_plot[index], "H" => values[1, index],
                "F" => values[2, index], "G" => values[3, index],
                "T" => values[4, index]))
        end
    end
    write_csv_rows(output, rows;
                   columns=["Hinf", "Tw", "eta", "H", "F", "G", "T"])
end

function plot_branch(solution_sets, output)
    p1 = plot(; xlabel="Tw", ylabel="Hinf", title="Infinite-mapped branch",
              gridalpha=0.25)
    p2 = plot(; xlabel="Hinf", ylabel="Tw",
              title="Turning points in continuation coordinate", gridalpha=0.25)
    for (operators, solutions, folds) in solution_sets
        h = [item.hinf for item in solutions]
        tw = [item.tw for item in solutions]
        label = "N=$(operators.degree), L=$(operators.scale)"
        plot!(p1, tw, h; label, linewidth=1.6)
        plot!(p2, h, tw; label, linewidth=1.6)
        for fold in folds
            scatter!(p1, [fold["Tw"]], [fold["Hinf"]]; marker=:star5,
                     markersize=7, label=false)
            scatter!(p2, [fold["Hinf"]], [fold["Tw"]]; marker=:star5,
                     markersize=7, label=false)
        end
    end
    savefig(plot(p1, p2; layout=(1, 2), size=(1250, 480)), output)
end

function interpolate_solution_to_operators(source::MappedOperators,
                                           solution::MappedSolution,
                                           target::MappedOperators)
    values = barycentric_interpolate(source.x, source.weights, solution.fields,
                                     target.x)
    values[:, end] .= [solution.hinf, 0, 1, 1]
    return values
end

function mapped_similarity_spectrum(source_operators::MappedOperators,
                                    solution::MappedSolution, degree::Int)
    operators = rational_chebyshev(degree, source_operators.scale)
    fields = interpolate_solution_to_operators(source_operators, solution, operators)
    h_base, f_base, g_base, t_base = (view(fields, row, :) for row in 1:4)
    d1, d2, dx = operators.deta, operators.deta2, operators.dx
    fp, gp, tp = d1 * f_base, d1 * g_base, d1 * t_base
    nodes = length(operators.x)
    rf, rg, rh, rt = ntuple(i -> ((i - 1) * nodes + 1):(i * nodes), 4)
    identity = Matrix{Float64}(I, nodes, nodes)
    A, B = zeros(4nodes, 4nodes), zeros(4nodes, 4nodes)
    A[rf, rf] .= d2 - Diagonal(h_base) * d1 - Diagonal(2f_base)
    A[rf, rg] .= Diagonal(2 .* (g_base .- 1))
    A[rf, rh] .= -Diagonal(fp)
    A[rf, rt] .= -identity
    A[rg, rf] .= Diagonal(2 .* (1 .- g_base))
    A[rg, rg] .= d2 - Diagonal(h_base) * d1 - Diagonal(2f_base)
    A[rg, rh] .= -Diagonal(gp)
    A[rh, rf] .= 2identity
    A[rh, rh] .= d1
    A[rt, rh] .= -Diagonal(tp)
    A[rt, rt] .= d2 / PR - Diagonal(h_base) * d1
    B[rf, rf] .= identity
    B[rg, rg] .= identity
    B[rt, rt] .= identity
    function set_bc(row, column)
        A[row, :] .= 0
        B[row, :] .= 0
        A[row, column] = 1
    end
    set_bc(first(rf), first(rf)); set_bc(last(rf), last(rf))
    set_bc(first(rg), first(rg)); set_bc(last(rg), last(rg))
    set_bc(first(rh), first(rh))
    far_h = last(rh)
    A[far_h, :] .= 0; B[far_h, :] .= 0
    A[far_h, rh] .= dx[end, :]
    set_bc(first(rt), first(rt)); set_bc(last(rt), last(rt))
    values = eigen(A, B).values
    mask = isfinite.(real.(values)) .& isfinite.(imag.(values)) .& (abs.(values) .< 1e4)
    values = ComplexF64.(values[mask])
    return values[sortperm(eachindex(values);
                           by=i -> (-real(values[i]), abs(imag(values[i]))))]
end

function temporal_sample_rows(operators, solutions, folds, temporal_degrees)
    reference_fold_h = Dict(1 => -0.5327616566962796,
                            2 => -0.11330573101124149)
    targets = Tuple{String,Float64}[("principal", -0.60), ("middle", -0.30),
        ("middle_near_second", -0.15), ("upper_after_second", -0.10),
        ("upper", -0.08), ("upper", -0.05), ("near_zero_limit", -0.02)]
    for fold in folds[1:min(3, end)]
        index = Int(fold["fold_index"])
        push!(targets, ("fold_$index", get(reference_fold_h, index, fold["Hinf"])))
    end
    rows = Dict{String,Any}[]
    solution_h = [item.hinf for item in solutions]
    for (label, target) in targets
        nearest = solutions[argmin(abs.(solution_h .- target))]
        initial = vcat(BoussinesqSaddleNode.pack_fields(nearest.fields), nearest.tw)
        exact = mapped_newton(target, initial, operators; tolerance=3e-10)
        for degree in temporal_degrees
            values = mapped_similarity_spectrum(operators, exact, degree)
            near_zero = values[argmin(abs.(values))]
            positive = count(>(1e-7), real.(values))
            for (rank, value) in enumerate(values[1:min(12, end)])
                push!(rows, Dict("sample" => label, "Hinf" => target,
                    "Tw" => exact.tw, "temporal_degree" => degree,
                    "rank_by_real_part" => rank, "eigen_real" => real(value),
                    "eigen_imag" => imag(value), "positive_real_count" => positive,
                    "closest_to_zero_real" => real(near_zero),
                    "closest_to_zero_imag" => imag(near_zero)))
            end
        end
    end
    return rows
end

function main(arguments=ARGS)
    defaults = Dict{String,Any}("degree" => [100], "scale" => [8.0],
        "h-stop" => -0.02, "tolerance" => 2e-10,
        "temporal-degrees" => Int[], "output-dir" => DEFAULT_OUTPUT)
    args = parse_cli(arguments, defaults)
    output = abspath(args["output-dir"])
    mkpath(output)
    solution_sets = Any[]
    fold_summary = Dict{String,Any}[]
    for degree in args["degree"], scale in args["scale"]
        operators, solutions = trace_mapped_branch(degree, scale, args["h-stop"],
                                                   args["tolerance"])
        folds = locate_folds(solutions)
        push!(solution_sets, (operators, solutions, folds))
        tag = replace("N$(degree)_L$(@sprintf("%g", scale))", "." => "p")
        write_csv_rows(joinpath(output, "branch_$tag.csv"),
                       branch_rows(operators, solutions, folds))
        save_selected_profiles(operators, solutions,
                               joinpath(output, "profiles_$tag.csv"))
        for fold in folds
            push!(fold_summary, merge(Dict("degree" => degree,
                "mapping_scale" => scale), fold))
        end
    end
    write_csv_rows(joinpath(output, "fold_convergence.csv"), fold_summary)
    plot_branch(solution_sets, joinpath(output, "infinite_mapped_branch_convergence.png"))
    write_json(joinpath(output, "run_summary.json"), Dict(
        "mapping" => "eta=L(1-x)/(1+x)", "degrees" => args["degree"],
        "scales" => args["scale"], "h_stop" => args["h-stop"],
        "folds" => fold_summary))
    if !isempty(args["temporal-degrees"])
        operators, solutions, folds = last(solution_sets)
        rows = temporal_sample_rows(operators, solutions, folds,
                                    args["temporal-degrees"])
        write_csv_rows(joinpath(output, "mapped_temporal_samples.csv"), rows)
    end
    JSON3.pretty(stdout, fold_summary; allow_inf=true)
    println()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
