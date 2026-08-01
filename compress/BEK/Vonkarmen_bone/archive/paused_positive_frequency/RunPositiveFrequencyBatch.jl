using Printf

const WALL_TEMPERATURES = [1.00, 1.04, 1.08, 1.12, 1.16, 1.20]
const MODELS = [:blackburn, :compressible]

function tw_text(Tw)
    text = replace(@sprintf("%.4f", Tw), r"0+$" => "")
    endswith(text, ".") && (text *= "0")
    return text
end

function curve_path(output_dir, model, omega, Tw)
    prefix = @sprintf("ome=%s_Tw=%s_model=%s", omega, tw_text(Tw), model)
    filename = model === :blackburn ?
        "$(prefix).dat" :
        "$(prefix)_Mr=0.3_propPert=on_baseProp=variable.dat"
    return joinpath(output_dir, filename)
end

function completed_curve(path)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) >= 5 || return false
    R_values = Float64[]
    for line in lines[3:end]
        fields = split(strip(line))
        length(fields) == 7 || return false
        R = tryparse(Float64, fields[2])
        R === nothing && return false
        push!(R_values, R)
    end
    any(
        R_values[index] < R_values[index - 1] &&
        R_values[index] < R_values[index + 1]
        for index in 2:length(R_values)-1
    ) || return false
    return true
end

function launch_case(output_dir, log_dir, model, omega, N_cheb, Tw)
    log_path = joinpath(
        log_dir,
        @sprintf("model=%s_omega=%.3f_Tw=%.2f.log", model, omega, Tw),
    )
    stream = open(log_path, "w")
    environment = copy(ENV)
    environment["POSITIVE_FREQUENCY_MODEL"] = string(model)
    environment["POSITIVE_FREQUENCY_OMEGA"] = string(omega)
    environment["POSITIVE_FREQUENCY_N_CHEB"] = string(N_cheb)
    environment["POSITIVE_FREQUENCY_OUTPUT_DIR"] = output_dir
    command = `$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(
        @__DIR__, "ComputePositiveFrequencyNeutralCurves.jl",
    )) $Tw`
    process = run(
        pipeline(setenv(command, environment); stdout=stream, stderr=stream);
        wait=false,
    )
    println("launched model=$(model) Tw=$(Tw) pid=$(getpid(process))")
    flush(stdout)
    return (
        process=process,
        stream=stream,
        log_path=log_path,
        model=model,
        Tw=Tw,
    )
end

function finish_case(job)
    wait(job.process)
    close(job.stream)
    println(
        "finished model=$(job.model) Tw=$(job.Tw) " *
        "exitcode=$(job.process.exitcode) log=$(job.log_path)",
    )
    flush(stdout)
    return job.process.exitcode == 0
end

function main()
    omega = parse(Float64, get(ENV, "POSITIVE_FREQUENCY_OMEGA", "0.008"))
    N_cheb = parse(Int, get(ENV, "POSITIVE_FREQUENCY_N_CHEB", "69"))
    output_dir = get(
        ENV,
        "POSITIVE_FREQUENCY_OUTPUT_DIR",
        joinpath(
            @__DIR__,
            "positive_frequency_neutral_curve_batch",
            @sprintf("omega=%.3f", omega),
        ),
    )
    max_workers = parse(Int, get(ENV, "POSITIVE_FREQUENCY_WORKERS", "4"))
    max_workers >= 1 || error("POSITIVE_FREQUENCY_WORKERS must be positive")
    log_dir = joinpath(output_dir, "process_logs")
    mkpath(output_dir)
    mkpath(log_dir)

    cases = [
        (model=model, Tw=Tw)
        for model in MODELS for Tw in WALL_TEMPERATURES
    ]
    active = NamedTuple[]
    failures = NamedTuple[]
    skipped = 0

    for case in cases
        path = curve_path(output_dir, case.model, omega, case.Tw)
        if completed_curve(path)
            println("skip completed model=$(case.model) Tw=$(case.Tw)")
            skipped += 1
            continue
        end
        while length(active) >= max_workers
            job = popfirst!(active)
            finish_case(job) || push!(failures, job)
        end
        push!(active, launch_case(
            output_dir, log_dir, case.model, omega, N_cheb, case.Tw,
        ))
    end

    for job in active
        finish_case(job) || push!(failures, job)
    end

    println(
        "positive-frequency batch complete: total=$(length(cases)) " *
        "skipped=$(skipped) failures=$(length(failures))",
    )
    isempty(failures) || error(
        "failed cases: " *
        join(
            ["model=$(job.model),Tw=$(job.Tw)" for job in failures],
            "; ",
        ),
    )
end

main()
