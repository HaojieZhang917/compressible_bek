module BlackburnNeutralCurveRunner

using BSplineKit
using Dates
using DelimitedFiles
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "BlackburnStability.jl"))
include(joinpath(@__DIR__, "CRD_STA.jl"))
include(joinpath(@__DIR__, "LopezBaseflow.jl"))
include(joinpath(@__DIR__, "NeutralContinuation.jl"))

using .LopezBaseflow
using .NeutralContinuation

export CurveConfig, compute_neutral_curve, run_case_with_retries,
       run_standard_batch,
       run_parallel_standard_batch, validate_curve_file,
       validate_standard_batch, case_tag, curve_path, branch_curve_path,
       load_curve, write_curve, validate_curve_data

default_output_dir() = joinpath(pwd(), "neutral_curve_batch")

Base.@kwdef struct CurveConfig
    Tw::Float64
    omega::Float64 = 0.0
    R_initial::Float64 = 500.0
    beta_initial::Float64 = 0.04
    alpha_target::Float64 = 0.1
    num_modes::Int = 1
    model::Symbol = :blackburn
    Mr::Float64 = 0.3
    Ro::Float64 = -1.0
    N_cheb::Int = 69
    property_perturbations::Bool = true
    base_property_variation::Bool = true
    beta_step::Float64 = 8.0e-4
    beta_scan_step::Float64 = 5.0e-4
    beta_bounds::Tuple{Float64,Float64} = (1.0e-3, 0.20)
    min_beta_step::Float64 = 5.0e-5
    step_recovery_successes::Int = 4
    step_growth_factor::Float64 = 2.0
    neutral_tol::Float64 = 1.0e-7
    beta_tol::Float64 = 1.0e-9
    R_tol::Float64 = 1.0e-4
    corrector_R_step::Float64 = 1.0
    max_scan_steps::Int = 80
    max_prediction_step::Float64 = 30.0
    min_mode_overlap::Float64 = 0.60
    allow_mode_switch::Bool = true
    max_curve_points::Int = 500
    min_valid_points::Int = 8
    minimum_complete_points::Int = 50
    minimum_complete_beta::Float64 = 0.08
    output_dir::String = default_output_dir()
    keep_logs::Bool = false
end

function checked_config(config::CurveConfig)
    config.model in (:blackburn, :compressible) || throw(ArgumentError(
        "model must be :blackburn or :compressible",
    ))
    config.N_cheb >= 20 || throw(ArgumentError("N_cheb must be at least 20"))
    config.num_modes in (1, 2) || throw(ArgumentError("num_modes must be 1 or 2"))
    config.R_initial > 0 || throw(ArgumentError("R_initial must be positive"))
    config.beta_step != 0 || throw(ArgumentError("beta_step must be nonzero"))
    config.beta_scan_step > 0 || throw(ArgumentError("beta_scan_step must be positive"))
    0 < config.min_beta_step <= abs(config.beta_step) || throw(ArgumentError(
        "min_beta_step must lie in (0,abs(beta_step)]",
    ))
    config.step_recovery_successes >= 1 || throw(ArgumentError(
        "step_recovery_successes must be positive",
    ))
    config.step_growth_factor > 1 || throw(ArgumentError(
        "step_growth_factor must exceed one",
    ))
    config.neutral_tol > 0 || throw(ArgumentError("neutral_tol must be positive"))
    config.beta_tol > 0 || throw(ArgumentError("beta_tol must be positive"))
    config.R_tol > 0 || throw(ArgumentError("R_tol must be positive"))
    config.corrector_R_step > 0 || throw(ArgumentError(
        "corrector_R_step must be positive",
    ))
    config.max_scan_steps >= 1 || throw(ArgumentError(
        "max_scan_steps must be positive",
    ))
    config.max_curve_points >= 1 || throw(ArgumentError(
        "max_curve_points must be positive",
    ))
    config.min_valid_points >= 1 || throw(ArgumentError(
        "min_valid_points must be positive",
    ))
    config.minimum_complete_points >= config.min_valid_points || throw(ArgumentError(
        "minimum_complete_points must not be smaller than min_valid_points",
    ))
    first(config.beta_bounds) < config.minimum_complete_beta < last(config.beta_bounds) ||
        throw(ArgumentError(
        "minimum_complete_beta must lie inside beta_bounds",
    ))
    0 <= config.min_mode_overlap <= 1 || throw(ArgumentError(
        "min_mode_overlap must lie in [0,1]",
    ))
    first(config.beta_bounds) < config.beta_initial < last(config.beta_bounds) ||
        throw(ArgumentError("beta_initial must lie inside beta_bounds"))
    return config
end

function config_namedtuple(config::CurveConfig)
    names = fieldnames(CurveConfig)
    return NamedTuple{names}(Tuple(getfield(config,name) for name in names))
end

function with_config(config::CurveConfig; kwargs...)
    return CurveConfig(; merge(config_namedtuple(config), (; kwargs...))...)
end

function number_tag(value::Real; digits::Int=6)
    text = @sprintf("%.*f", digits, Float64(value))
    text = replace(text, r"0+$" => "")
    text = endswith(text, ".") ? text[1:end-1] : text
    return occursin('.',text) ? text : text * ".0"
end

function case_tag(config::CurveConfig)
    Tw = number_tag(config.Tw; digits=4)
    if config.model === :blackburn
        return "ome=$(number_tag(config.omega))_Tw=$(Tw)_model=blackburn"
    end
    perturbation = config.property_perturbations ? "on" : "off"
    base = config.base_property_variation ? "variable" : "frozen"
    return "ome=$(number_tag(config.omega))_Tw=$(Tw)_" *
           "model=compressible_Mr=$(number_tag(config.Mr))_" *
           "propPert=$(perturbation)_baseProp=$(base)"
end

curve_path(config::CurveConfig) = joinpath(config.output_dir, case_tag(config) * ".dat")
branch_curve_path(config::CurveConfig,branch) = joinpath(
    config.output_dir,case_tag(config) * "_branch=$(branch).dat",
)
log_path(config::CurveConfig) = joinpath(config.output_dir, case_tag(config) * ".log")
function status_path(config::CurveConfig)
    worker = strip(get(ENV,"NEUTRAL_WORKER",""))
    suffix = isempty(worker) ? "" : "_" * replace(worker,r"[^A-Za-z0-9_-]" => "_")
    return joinpath(config.output_dir,"batch_status$(suffix).tsv")
end

function append_log(path, message)
    open(path, "a") do io
        println(io, Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"), '\t', message)
        flush(io)
    end
end

function sanitized(message)
    return replace(String(message), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function append_status(config::CurveConfig, status, attempt, points, message)
    path = status_path(config)
    new_file = !isfile(path)
    open(path, "a") do io
        new_file && println(io, "time\tcase\tstatus\tattempt\tpoints\tmessage")
        fields = (
            Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            case_tag(config), status, attempt, points, sanitized(message),
        )
        println(io, join(fields, '\t'))
        flush(io)
    end
end

function completed_cases(path)
    result = Dict{String,String}()
    isfile(path) || return result
    for (line_number, line) in enumerate(eachline(path))
        line_number == 1 && startswith(line, "time\t") && continue
        fields = split(line, '\t'; limit=6)
        length(fields) >= 3 || continue
        result[fields[2]] = fields[3]
    end
    return result
end

function completed_cases_in_directory(output_dir)
    result = Dict{String,String}()
    isdir(output_dir) || return result
    for path in sort(filter(
        path -> startswith(basename(path),"batch_status") && endswith(path,".tsv"),
        readdir(output_dir; join=true),
    ))
        merge!(result,completed_cases(path))
    end
    return result
end

function sample_baseflow_profile(z, values, points)
    interpolation = BSplineKit.interpolate(
        Float64.(vec(z)), Float64.(vec(values)), BSplineOrder(4),
    )
    return Float64.(interpolation.(Float64.(vec(points))))
end

function prepare_solver(config::CurveConfig)
    D, D2, x = CRD_BF.Cheb(config.N_cheb)
    if config.model === :blackburn
        z, H0, F0, G0, T0, _, _, _, info = LopezBaseflow.get_baseflow(config.Tw)
        F = sample_baseflow_profile(z,F0,x)
        G = sample_baseflow_profile(z,G0,x)
        H = sample_baseflow_profile(z,H0,x)
        T = sample_baseflow_profile(z,T0,x)
        solve_at = function (R, beta, target, nvalues; v0=nothing)
            return BlackburnStability.eigsol_blackburn(
                F,G,H,T,R,config.omega,beta,config.N_cheb,D,D2,
                target,min(nvalues,2); initial_vector=v0,
            )
        end
        return (solve_at=solve_at, baseflow_info=info, D=D, D2=D2, x=x)
    end

    gamma = 1.4
    sigma = 0.72
    Co = 2 - config.Ro - config.Ro^2
    u0,v0,w0,f,q,_,_,_ = baseflow_var(config.N_cheb,config.Ro,Co)
    H0,T0 = T_ca(config.Mr,f,q,w0,gamma,config.Tw)
    F,G,H,T,rho,_, = interp(
        u0,v0,H0,T0,x,config.N_cheb,"sim",
    )
    lambda = -(2/3) .* T
    kappa = T ./ sigma
    solve_at = function (R, beta, target, nvalues; v0=nothing)
        R > 0 || throw(ArgumentError("R must be positive"))
        Ma = config.Mr / R
        result = solve_spatial_mode(
            F,G,H,rho,lambda,kappa,T,sigma,gamma,R,Ma,
            config.omega,beta,config.N_cheb,config.Ro,Co,D,D2;
            target=target,neigs=min(nvalues,2),initial_vector=v0,
            maxit=800,tol=1.0e-10,balance=true,regularization=0.0,
            property_perturbations=config.property_perturbations,
            base_property_variation=config.base_property_variation,
        )
        return result.values,result.vectors
    end
    return (solve_at=solve_at, baseflow_info=nothing, D=D, D2=D2, x=x)
end

function continued_state(
    solve_at, R, beta, seed_values, seed_vectors, active_index, config,
)
    values,vectors,overlaps = NeutralContinuation.continue_modes(
        solve_at,R,beta,seed_values,seed_vectors;
        min_overlap=config.min_mode_overlap,
    )
    return (
        R=Float64(R), beta=Float64(beta), values=values, vectors=vectors,
        overlaps=overlaps, active_index=active_index,
        residual=imag(values[active_index]),
    )
end

select_initial_mode(values,target) = argmin(abs.(values .- target))

function refine_beta_bracket(solve_at, left, right, config)
    left.beta <= right.beta || ((left,right) = (right,left))
    left.residual * right.residual <= 0 || error("Invalid beta bracket")
    for iteration in 1:40
        best = abs(left.residual) <= abs(right.residual) ? left : right
        abs(best.residual) <= config.neutral_tol && return best

        denominator = right.residual - left.residual
        beta = abs(denominator) > eps(Float64) ?
            right.beta - right.residual * (right.beta-left.beta) / denominator :
            (left.beta + right.beta)/2
        guard = 0.1 * (right.beta-left.beta)
        if !(left.beta + guard < beta < right.beta - guard)
            beta = (left.beta + right.beta)/2
        end
        seed = abs(beta-left.beta) <= abs(right.beta-beta) ? left : right
        state = continued_state(
            solve_at,left.R,beta,seed.values,seed.vectors,
            left.active_index,config,
        )
        abs(state.residual) <= config.neutral_tol && return state
        if left.residual * state.residual <= 0
            right = state
        else
            left = state
        end
        if right.beta-left.beta <= config.beta_tol
            best = abs(left.residual) <= abs(right.residual) ? left : right
            abs(best.residual) <= config.neutral_tol || error(
                "Initial beta refinement stopped at residual=$(best.residual)",
            )
            return best
        end
    end
    best = abs(left.residual) <= abs(right.residual) ? left : right
    abs(best.residual) <= config.neutral_tol || error(
        "Initial beta refinement failed at residual=$(best.residual)",
    )
    return best
end

function find_initial_neutral(solve_at, config::CurveConfig, log_file)
    values,vectors = solve_at(
        config.R_initial,config.beta_initial,
        config.alpha_target,config.num_modes,
    )
    active_index = select_initial_mode(values,config.alpha_target)
    center = (
        R=config.R_initial,beta=config.beta_initial,
        values=values,vectors=vectors,overlaps=ones(length(values)),
        active_index=active_index,residual=imag(values[active_index]),
    )
    append_log(
        log_file,
        "initial beta=$(center.beta) R=$(center.R) values=$(center.values) " *
        "residual=$(center.residual)",
    )
    abs(center.residual) <= config.neutral_tol && return center

    preferred = center.residual > 0 ? 1 : -1
    lower_beta,upper_beta = config.beta_bounds
    for direction in (preferred,-preferred)
        previous = center
        max_steps = ceil(Int,(upper_beta-lower_beta)/config.beta_scan_step)
        for scan_index in 1:max_steps
            beta = center.beta + direction * scan_index * config.beta_scan_step
            lower_beta < beta < upper_beta || break
            state = continued_state(
                solve_at,config.R_initial,beta,previous.values,
                previous.vectors,active_index,config,
            )
            append_log(
                log_file,
                "initial_scan beta=$beta residual=$(state.residual) " *
                "overlap=$(state.overlaps)",
            )
            if previous.residual * state.residual <= 0
                root = refine_beta_bracket(solve_at,previous,state,config)
                append_log(
                    log_file,
                    "initial_neutral beta=$(root.beta) R=$(root.R) " *
                    "residual=$(root.residual)",
                )
                return root
            end
            previous = state
        end
    end
    error(
        "No initial neutral crossing was found at R=$(config.R_initial) " *
        "inside beta_bounds=$(config.beta_bounds)",
    )
end

function point_row(config, R, beta, values)
    first_value = values[1]
    last_value = values[end]
    return Float64[
        config.omega,R,beta,real(first_value),imag(first_value),
        real(last_value),imag(last_value),
    ]
end

function write_curve(path, config, data)
    mkpath(dirname(path))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        println(io, "Variables=\"omega\" \"R\" \"beta\" \"alpha_r_1\" \"alpha_i_1\" \"alpha_r_2\" \"alpha_i_2\"")
        println(io, "Zone T=\"$(case_tag(config))\"")
        writedlm(io,data)
    end
    mv(temporary,path; force=true)
    return path
end

function load_curve(path)
    lines = readlines(path)
    length(lines) >= 3 || return zeros(Float64,0,7)
    rows = Vector{Vector{Float64}}()
    for line in lines[3:end]
        isempty(strip(line)) && continue
        fields = split(strip(line))
        length(fields) == 7 || return zeros(Float64,0,7)
        push!(rows,parse.(Float64,fields))
    end
    isempty(rows) && return zeros(Float64,0,7)
    return reduce(vcat,permutedims.(rows))
end

function validate_curve_data(data, config::CurveConfig)
    issues = String[]
    size(data,2) == 7 || push!(issues,"expected 7 columns")
    size(data,1) >= config.min_valid_points || push!(
        issues,"only $(size(data,1)) points; expected at least $(config.min_valid_points)",
    )
    all(isfinite,data) || push!(issues,"non-finite values are present")
    if size(data,1) > 0 && size(data,2) == 7
        all((data[:,2] .> 0) .& (data[:,2] .< 700.001)) || push!(
            issues,"R left the configured physical interval",
        )
        beta_differences = diff(data[:,3])
        (all(beta_differences .> 0) || all(beta_differences .< 0)) || push!(
            issues,"beta is not strictly monotone",
        )
        residuals = min.(abs.(data[:,5]),abs.(data[:,7]))
        residual_limit = max(10config.neutral_tol,1.0e-6)
        maximum(residuals) <= residual_limit || push!(
            issues,"maximum neutral residual $(maximum(residuals)) exceeds $residual_limit",
        )
    end
    return (ok=isempty(issues),issues=issues)
end

function validate_curve_file(path, config::CurveConfig)
    isfile(path) || return (ok=false,issues=["file does not exist"],data=zeros(0,7))
    data = try
        load_curve(path)
    catch exception
        return (
            ok=false,issues=["failed to parse: $(sprint(showerror,exception))"],
            data=zeros(0,7),
        )
    end
    result = validate_curve_data(data,config)
    return (ok=result.ok,issues=result.issues,data=data)
end

function file_completion_issues(data,config::CurveConfig)
    issues = String[]
    size(data,1) >= config.minimum_complete_points || push!(
        issues,
        "only $(size(data,1)) points; a complete curve requires at least " *
        "$(config.minimum_complete_points)",
    )
    if size(data,1) > 0 && size(data,2) >= 3
        maximum(data[:,3]) >= config.minimum_complete_beta || push!(
            issues,
            "maximum beta=$(maximum(data[:,3])) is below " *
            "$(config.minimum_complete_beta)",
        )
    end
    return issues
end

function validate_standard_batch(;
    output_dir=default_output_dir(),
    N_cheb=69,cleanup_logs=false,
)
    configs = filter(
        !is_excluded_case,
        standard_configs(; output_dir=output_dir,N_cheb=N_cheb),
    )
    invalid = Dict{String,Vector{String}}()
    summary_path = joinpath(output_dir,"batch_validation.tsv")
    open(summary_path,"w") do io
        println(
            io,
            "case\tstatus\tpoints\tR_min\tR_max\tbeta_min\tbeta_max\t" *
            "max_neutral_residual\tbeta_direction\tissues",
        )
        for config in configs
            validation = validate_curve_file(curve_path(config),config)
            issues = vcat(
                validation.issues,
                file_completion_issues(validation.data,config),
            )
            isempty(issues) || (invalid[case_tag(config)] = issues)
            data = validation.data
            if size(data,1) == 0
                metrics = (0,NaN,NaN,NaN,NaN,NaN,"none")
            else
                beta_differences = diff(data[:,3])
                direction = all(beta_differences .> 0) ? "increasing" :
                    all(beta_differences .< 0) ? "decreasing" : "mixed"
                residuals = min.(abs.(data[:,5]),abs.(data[:,7]))
                metrics = (
                    size(data,1),minimum(data[:,2]),maximum(data[:,2]),
                    minimum(data[:,3]),maximum(data[:,3]),maximum(residuals),
                    direction,
                )
            end
            fields = (
                case_tag(config),isempty(issues) ? "ok" : "invalid",
                metrics...,isempty(issues) ? "" : join(issues,"; "),
            )
            println(io,join(fields,'\t'))
        end
    end
    if isempty(invalid) && cleanup_logs
        for path in readdir(output_dir; join=true)
            endswith(path,".log") && rm(path; force=true)
        end
    end
    return (
        ok=isempty(invalid),invalid=invalid,summary_path=summary_path,
        required_cases=length(configs),
    )
end

function find_alternate_mode_root(
    solve_at,beta,R_guess,values,vectors,active_index,config,on_evaluation,
)
    length(values) >= 2 || return nothing
    candidates = NamedTuple[]
    for candidate_index in eachindex(values)
        candidate_index == active_index && continue
        root = try
            NeutralContinuation.find_neutral_R(
                solve_at,beta,R_guess,values,vectors,candidate_index;
                R_step=min(config.corrector_R_step,0.5),
                preferred_direction=0,
                max_scan_steps=max(config.max_scan_steps,160),max_refine=40,
                neutral_tol=config.neutral_tol,R_tol=config.R_tol,
                R_bounds=(1.0e-6,700.0),max_R_deviation=150.0,
                min_overlap=config.min_mode_overlap,
                on_evaluation=on_evaluation,
            )
        catch exception
            exception isa InterruptException && rethrow()
            nothing
        end
        root === nothing || push!(candidates,root)
    end
    isempty(candidates) && return nothing
    return candidates[argmin([abs(root.R-R_guess) for root in candidates])]
end

function compute_neutral_curve(config::CurveConfig)
    checked_config(config)
    mkpath(config.output_dir)
    output_file = curve_path(config)
    trace_file = log_path(config)
    isfile(trace_file) && rm(trace_file; force=true)
    append_log(trace_file,"starting $(case_tag(config)) N=$(config.N_cheb)")

    prepared = prepare_solver(config)
    solve_at = prepared.solve_at
    initial = find_initial_neutral(solve_at,config,trace_file)
    values,vectors = initial.values,initial.vectors
    active_index = initial.active_index
    data = reshape(point_row(config,initial.R,initial.beta,values),1,:)
    write_curve(output_file,config,data)

    beta_direction = sign(config.beta_step)
    nominal_beta_step = abs(config.beta_step)
    current_beta_step = nominal_beta_step
    successful_steps_at_reduced_step = 0
    attempted_beta = initial.beta + beta_direction * current_beta_step
    stop_reason = :unknown
    while size(data,1) < config.max_curve_points
        if !(first(config.beta_bounds) < attempted_beta < last(config.beta_bounds))
            stop_reason = :beta_limit
            break
        end
        previous_R = data[end,2]
        predicted_delta_R = 0.0
        if size(data,1) >= 2
            delta_beta = data[end,3]-data[end-1,3]
            if abs(delta_beta) > eps(Float64)
                predicted_delta_R = (
                    (data[end,2]-data[end-1,2]) / delta_beta *
                    (attempted_beta-data[end,3])
                )
            end
        end
        predicted_delta_R = clamp(
            predicted_delta_R,-config.max_prediction_step,config.max_prediction_step,
        )
        R_guess = clamp(previous_R+predicted_delta_R,1.0e-5,699.999)
        preferred_direction = Int(sign(predicted_delta_R))
        evaluation_count = Ref(0)
        on_evaluation = state -> begin
            evaluation_count[] += 1
            append_log(
                trace_file,
                "corrector beta=$attempted_beta R=$(state.R) " *
                "residual=$(state.residual) overlap=$(state.overlaps)",
            )
        end

        root = try
            NeutralContinuation.find_neutral_R(
                solve_at,attempted_beta,R_guess,values,vectors,active_index;
                R_step=config.corrector_R_step,
                preferred_direction=preferred_direction,
                max_scan_steps=config.max_scan_steps,max_refine=40,
                neutral_tol=config.neutral_tol,R_tol=config.R_tol,
                R_bounds=(1.0e-6,700.0),
                max_R_deviation=max(40.0,2abs(predicted_delta_R)+10.0),
                min_overlap=config.min_mode_overlap,
                on_evaluation=on_evaluation,
            )
        catch exception
            exception isa InterruptException && rethrow()
            smaller_step = NeutralContinuation.refined_beta_step(
                current_beta_step,config.min_beta_step,
            )
            if smaller_step !== nothing
                append_log(
                    trace_file,
                    "retry failed_beta=$attempted_beta old_step=$current_beta_step " *
                    "new_step=$smaller_step reason=$(sprint(showerror,exception))",
                )
                current_beta_step = smaller_step
                successful_steps_at_reduced_step = 0
                attempted_beta = data[end,3] + beta_direction * current_beta_step
                continue
            end
            alternate_root = config.allow_mode_switch ?
                find_alternate_mode_root(
                    solve_at,attempted_beta,previous_R,values,vectors,
                    active_index,config,on_evaluation,
                ) :
                nothing
            if alternate_root !== nothing
                append_log(
                    trace_file,
                    "mode_switch beta=$attempted_beta from=$active_index " *
                    "to=$(alternate_root.active_index) R=$(alternate_root.R) " *
                    "residual=$(alternate_root.residual)",
                )
                successful_steps_at_reduced_step = 0
                alternate_root
            else
                stop_reason = :endpoint_no_root
                append_log(
                    trace_file,
                    "stop reason=$stop_reason failed_beta=$attempted_beta " *
                    "last_beta=$(data[end,3]) last_R=$(data[end,2]) " *
                    "message=$(sprint(showerror,exception))",
                )
                break
            end
        end

        abs(root.residual) <= config.neutral_tol || error(
            "Corrector returned residual=$(root.residual)",
        )
        values,vectors = root.values,root.vectors
        active_index = root.active_index
        data = vcat(data,permutedims(point_row(
            config,root.R,attempted_beta,values,
        )))
        write_curve(output_file,config,data)
        append_log(
            trace_file,
            "accepted beta=$attempted_beta R=$(root.R) " *
            "residual=$(root.residual) evaluations=$(evaluation_count[])",
        )

        if current_beta_step < nominal_beta_step
            successful_steps_at_reduced_step += 1
            recovered_step = NeutralContinuation.recovered_beta_step(
                current_beta_step,nominal_beta_step,
                successful_steps_at_reduced_step,config.step_recovery_successes;
                factor=config.step_growth_factor,
            )
            if recovered_step > current_beta_step
                old_step = current_beta_step
                current_beta_step = recovered_step
                successful_steps_at_reduced_step = 0
                append_log(
                    trace_file,
                    "step_recovery old_step=$old_step new_step=$current_beta_step",
                )
            end
        else
            successful_steps_at_reduced_step = 0
        end

        if root.R > 500 && size(data,1) > 30
            stop_reason = :R_limit
            break
        end
        attempted_beta += beta_direction * current_beta_step
    end
    stop_reason === :unknown && (stop_reason = :point_limit)
    write_curve(output_file,config,data)
    validation = validate_curve_data(data,config)
    append_log(
        trace_file,
        "finished reason=$stop_reason points=$(size(data,1)) valid=$(validation.ok) " *
        "issues=$(join(validation.issues,"; "))",
    )
    return (
        data=data,path=output_file,log_path=trace_file,
        stop_reason=stop_reason,validation=validation,
    )
end

function retry_configs(config::CurveConfig)
    robust = with_config(
        config;
        num_modes=max(config.num_modes,2),
        min_beta_step=min(config.min_beta_step,abs(config.beta_step)/32),
        corrector_R_step=min(config.corrector_R_step,0.5),
        max_scan_steps=max(config.max_scan_steps,160),
        min_mode_overlap=min(config.min_mode_overlap,0.50),
    )
    return (config,robust)
end

function curve_completion_issues(result,config::CurveConfig)
    issues = copy(result.validation.issues)
    if result.stop_reason === :endpoint_no_root &&
       (size(result.data,1) < config.minimum_complete_points ||
        maximum(result.data[:,3]) < config.minimum_complete_beta)
        push!(
            issues,
            "premature endpoint at beta=$(result.data[end,3]) with " *
            "$(size(result.data,1)) points",
        )
    end
    return issues
end

function typeI_fallback_config(config::CurveConfig)
    return with_config(
        config;
        beta_initial=0.12,alpha_target=0.65,
        beta_step=-abs(config.beta_step),num_modes=2,
        min_beta_step=min(config.min_beta_step,abs(config.beta_step)/32),
        corrector_R_step=min(config.corrector_R_step,0.5),
        max_scan_steps=max(config.max_scan_steps,160),
        min_mode_overlap=min(config.min_mode_overlap,0.50),
    )
end

function run_case_with_retries(config::CurveConfig)
    last_exception = nothing
    last_incomplete_result = nothing
    for (attempt, attempt_config) in enumerate(retry_configs(config))
        append_status(attempt_config,"running",attempt,0,"starting")
        result = try
            compute_neutral_curve(attempt_config)
        catch exception
            exception isa InterruptException && rethrow()
            last_exception = exception
            append_status(
                attempt_config,"retry",attempt,0,sprint(showerror,exception),
            )
            continue
        end
        completion_issues = curve_completion_issues(result,attempt_config)
        if isempty(completion_issues)
            append_status(
                attempt_config,"ok",attempt,size(result.data,1),
                "stop_reason=$(result.stop_reason)",
            )
            attempt_config.keep_logs || rm(result.log_path; force=true)
            return result
        end
        message = join(completion_issues,"; ")
        last_incomplete_result = result
        append_status(
            attempt_config,"retry",attempt,size(result.data,1),message,
        )
        premature_endpoint = result.stop_reason === :endpoint_no_root &&
            (size(result.data,1) < attempt_config.minimum_complete_points ||
             maximum(result.data[:,3]) < attempt_config.minimum_complete_beta)
        premature_endpoint && break
    end

    typeII_path = nothing
    if last_incomplete_result !== nothing
        typeII_path = branch_curve_path(config,:typeII)
        cp(last_incomplete_result.path,typeII_path; force=true)
    end

    typeI_config = typeI_fallback_config(config)
    typeI_message = typeII_path === nothing ?
        "base-branch initialization failed; trying Type-I fallback" :
        "saved disconnected Type-II branch to $(basename(typeII_path))"
    append_status(typeI_config,"running_typeI",3,0,typeI_message)
    typeI_result = try
        compute_neutral_curve(typeI_config)
    catch exception
        exception isa InterruptException && rethrow()
        last_exception = exception
        nothing
    end
    if typeI_result !== nothing
        typeI_issues = curve_completion_issues(typeI_result,typeI_config)
        if isempty(typeI_issues)
            branch_message = typeII_path === nothing ? "" :
                "; disconnected_typeII=$(basename(typeII_path))"
            append_status(
                typeI_config,"ok",3,size(typeI_result.data,1),
                "stop_reason=$(typeI_result.stop_reason)$(branch_message)",
            )
            typeI_config.keep_logs || rm(typeI_result.log_path; force=true)
            return typeII_path === nothing ? typeI_result :
                merge(typeI_result,(typeII_path=typeII_path,))
        end
        append_status(
            typeI_config,"retry",3,size(typeI_result.data,1),
            join(typeI_issues,"; "),
        )
    end
    message = last_exception === nothing ? "curve validation failed" :
        sprint(showerror,last_exception)
    append_status(config,"failed",3,0,message)
    error("Failed $(case_tag(config)): $message")
end

function standard_configs(;
    output_dir=default_output_dir(),
    N_cheb=69, keep_logs=false,
)
    temperatures = [1.06 + 0.02index for index in 0:7]
    configs = CurveConfig[]
    for Tw in temperatures
        common = (
            Tw=Float64(round(Tw; digits=8)),omega=0.0,R_initial=500.0,
            beta_initial=0.04,alpha_target=0.1,num_modes=1,Mr=0.3,
            Ro=-1.0,N_cheb=N_cheb,beta_step=8.0e-4,
            neutral_tol=1.0e-7,output_dir=String(output_dir),
            keep_logs=keep_logs,
        )
        push!(configs,CurveConfig(; common...,model=:blackburn))
        push!(configs,CurveConfig(
            ; common...,model=:compressible,
            property_perturbations=true,base_property_variation=true,
        ))
        push!(configs,CurveConfig(
            ; common...,model=:compressible,
            property_perturbations=false,base_property_variation=true,
        ))
        push!(configs,CurveConfig(
            ; common...,model=:compressible,
            property_perturbations=false,base_property_variation=false,
        ))
    end
    return configs
end

is_excluded_case(::CurveConfig) = false

function configs_for_worker(configs,worker::Symbol)
    worker === :all && return configs
    return filter(configs) do config
        if worker === :blackburn
            return config.model === :blackburn
        elseif config.model !== :compressible
            return false
        elseif worker === :compressible_full
            return config.property_perturbations &&
                   config.base_property_variation
        elseif worker === :compressible_no_perturb
            return !config.property_perturbations &&
                   config.base_property_variation
        elseif worker === :compressible_frozen
            return !config.property_perturbations &&
                   !config.base_property_variation
        end
        return false
    end
end

function run_standard_batch(;
    output_dir=default_output_dir(),
    N_cheb=69,keep_logs=false,resume=true,continue_on_error=true,
    worker::Symbol=:all,prefer_typeI=false,
)
    mkpath(output_dir)
    configs = standard_configs(
        ; output_dir=output_dir,N_cheb=N_cheb,keep_logs=keep_logs,
    )
    worker in (:all,:blackburn,:compressible_full,:compressible_no_perturb,
               :compressible_frozen) || throw(ArgumentError("unknown worker=$worker"))
    configs = configs_for_worker(configs,worker)
    finished = completed_cases_in_directory(output_dir)
    results = Dict{String,Any}()
    for (index,config) in enumerate(configs)
        key = case_tag(config)
        if is_excluded_case(config)
            message = "excluded by user request"
            append_status(config,"skipped",0,0,message)
            results[key] = (status=:excluded,message=message)
            println("[$index/$(length(configs))] excluded $key")
            flush(stdout)
            continue
        end
        if resume && get(finished,key,"") == "ok"
            validation = validate_curve_file(curve_path(config),config)
            completion_issues = file_completion_issues(validation.data,config)
            if validation.ok && isempty(completion_issues)
                println("[$index/$(length(configs))] skip completed $key")
                results[key] = (status=:skipped,path=curve_path(config))
                continue
            end
            append_status(
                config,"retry",0,size(validation.data,1),
                "completed file failed revalidation: " *
                join(vcat(validation.issues,completion_issues),"; "),
            )
        end

        println("[$index/$(length(configs))] computing $key")
        flush(stdout)
        try
            solve_config = prefer_typeI ? typeI_fallback_config(config) : config
            result = run_case_with_retries(solve_config)
            results[key] = (status=:ok,path=result.path,points=size(result.data,1))
            println(
                "[$index/$(length(configs))] completed $key " *
                "points=$(size(result.data,1)) stop=$(result.stop_reason)",
            )
        catch exception
            exception isa InterruptException && rethrow()
            results[key] = (status=:failed,message=sprint(showerror,exception))
            println(stderr,"[$index/$(length(configs))] FAILED $key: $(sprint(showerror,exception))")
            continue_on_error || rethrow()
        end
        flush(stdout)
        flush(stderr)
    end
    return results
end

function run_parallel_standard_batch(;
    output_dir=default_output_dir(),
    N_cheb=69,keep_logs=false,resume=true,
)
    mkpath(output_dir)
    worker_names = (
        :blackburn,:compressible_full,
        :compressible_no_perturb,:compressible_frozen,
    )
    project_dir = dirname(Base.active_project())
    script_path = abspath(@__FILE__)
    processes = Any[]
    streams = IO[]
    worker_logs = String[]

    for worker in worker_names
        worker_log = joinpath(output_dir,"worker_$(worker).log")
        stream = open(worker_log,"w")
        environment = copy(ENV)
        environment["NEUTRAL_WORKER"] = String(worker)
        environment["NEUTRAL_PARALLEL"] = "false"
        environment["NEUTRAL_OUTPUT_DIR"] = String(output_dir)
        environment["NEUTRAL_N_CHEB"] = string(N_cheb)
        environment["NEUTRAL_KEEP_LOGS"] = keep_logs ? "true" : "false"
        environment["NEUTRAL_RESUME"] = resume ? "true" : "false"
        command = `$(Base.julia_cmd()) --project=$(project_dir) $(script_path)`
        process = run(
            pipeline(setenv(command,environment); stdout=stream,stderr=stream);
            wait=false,
        )
        push!(processes,process)
        push!(streams,stream)
        push!(worker_logs,worker_log)
        println("launched worker=$worker pid=$(getpid(process))")
    end
    flush(stdout)

    for (worker,process,stream) in zip(worker_names,processes,streams)
        wait(process)
        close(stream)
        println("worker=$worker exitcode=$(process.exitcode)")
        flush(stdout)
    end

    configs = standard_configs(
        ; output_dir=output_dir,N_cheb=N_cheb,keep_logs=keep_logs,
    )
    configs = filter(!is_excluded_case,configs)
    invalid = Dict{String,Vector{String}}()
    for config in configs
        validation = validate_curve_file(curve_path(config),config)
        issues = vcat(
            validation.issues,file_completion_issues(validation.data,config),
        )
        isempty(issues) || (invalid[case_tag(config)] = issues)
    end
    process_failures = [
        String(worker) for (worker,process) in zip(worker_names,processes)
        if process.exitcode != 0
    ]
    if isempty(invalid) && isempty(process_failures)
        keep_logs || foreach(path -> rm(path; force=true),worker_logs)
        println(
            "parallel batch validated all $(length(configs)) required curves",
        )
        return (ok=true,invalid=invalid,process_failures=process_failures)
    end
    error(
        "Parallel batch validation failed: process_failures=$process_failures, " *
        "invalid=$(collect(keys(invalid)))",
    )
end

function main(args=ARGS)
    output_dir = get(ENV,"NEUTRAL_OUTPUT_DIR",default_output_dir())
    N_cheb = parse(Int,get(ENV,"NEUTRAL_N_CHEB","69"))
    keep_logs = lowercase(get(ENV,"NEUTRAL_KEEP_LOGS","false")) in ("1","true","yes")
    resume = !(lowercase(get(ENV,"NEUTRAL_RESUME","true")) in ("0","false","no"))
    worker = Symbol(lowercase(get(ENV,"NEUTRAL_WORKER","all")))
    parallel = lowercase(get(ENV,"NEUTRAL_PARALLEL","true")) in ("1","true","yes")
    prefer_typeI = lowercase(get(ENV,"NEUTRAL_PREFER_TYPEI","false")) in
        ("1","true","yes")
    validate_only = lowercase(get(ENV,"NEUTRAL_VALIDATE_ONLY","false")) in
        ("1","true","yes")
    if validate_only
        result = validate_standard_batch(
            ; output_dir=output_dir,N_cheb=N_cheb,cleanup_logs=!keep_logs,
        )
        result.ok || error("batch validation failed: $(result.invalid)")
        println(
            "validated $(result.required_cases) required curves; " *
            "summary=$(result.summary_path)",
        )
        return result
    end
    if parallel && worker === :all
        run_parallel_standard_batch(
            ; output_dir=output_dir,N_cheb=N_cheb,
            keep_logs=keep_logs,resume=resume,
        )
    else
        run_standard_batch(
            ; output_dir=output_dir,N_cheb=N_cheb,
            keep_logs=keep_logs,resume=resume,worker=worker,
            prefer_typeI=prefer_typeI,
        )
    end
end

end # module BlackburnNeutralCurveRunner

if abspath(PROGRAM_FILE) == @__FILE__
    using .BlackburnNeutralCurveRunner
    BlackburnNeutralCurveRunner.main()
end
