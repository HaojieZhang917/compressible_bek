using Printf

const EXTENDED_TEMPERATURES = (1.40, 1.60, 1.80)
const EXTENDED_MODELS = (:blackburn, :compressible)

function tw_text(Tw::Real)
    text = replace(@sprintf("%.4f", Tw), r"0+$" => "")
    return endswith(text, ".") ? text * "0" : text
end

function output_path(model::Symbol, Tw::Real)
    filename = model === :blackburn ?
        "ome=0.0_Tw=$(tw_text(Tw))_model=blackburn.dat" :
        "ome=0.0_Tw=$(tw_text(Tw))_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat"
    directory = model === :blackburn ?
        "blackburn_neutral_curve_batch" : "neutral_curve_batch"
    return joinpath(@__DIR__, directory, filename)
end

function launch_case(log_dir, model::Symbol, Tw::Real, N_cheb::Integer)
    log_path = joinpath(
        log_dir,
        @sprintf("model=%s_Tw=%.2f.log", String(model), Tw),
    )
    stream = open(log_path, "w")
    environment = copy(ENV)
    environment["ZERO_FREQUENCY_MODEL"] = String(model)
    environment["ZERO_FREQUENCY_N_CHEB"] = string(N_cheb)
    command = `$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(
        @__DIR__, "ComputeExtendedZeroFrequencyNeutralCurves.jl",
    )) $Tw`
    process = run(
        pipeline(setenv(command, environment); stdout=stream, stderr=stream);
        wait=false,
    )
    @printf("launched model=%s Tw=%.2f pid=%d\n",
        String(model), Tw, getpid(process))
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
    @printf("finished model=%s Tw=%.2f exitcode=%d log=%s\n",
        String(job.model), job.Tw, job.process.exitcode, job.log_path)
    flush(stdout)
    return job.process.exitcode == 0
end

function main()
    N_cheb = parse(Int, get(ENV, "ZERO_FREQUENCY_N_CHEB", "69"))
    max_workers = parse(Int, get(ENV, "ZERO_FREQUENCY_WORKERS", "3"))
    max_workers >= 1 || error("ZERO_FREQUENCY_WORKERS must be positive")
    log_dir = joinpath(@__DIR__, "zero_frequency_extended_process_logs")
    mkpath(log_dir)

    cases = [
        (model=model, Tw=Tw)
        for Tw in EXTENDED_TEMPERATURES for model in EXTENDED_MODELS
    ]
    active = NamedTuple[]
    failures = NamedTuple[]
    for case in cases
        while length(active) >= max_workers
            job = popfirst!(active)
            finish_case(job) || push!(failures, job)
        end
        push!(active, launch_case(log_dir, case.model, case.Tw, N_cheb))
    end
    for job in active
        finish_case(job) || push!(failures, job)
    end

    isempty(failures) || error(
        "failed cases: " * join(
            ["model=$(job.model),Tw=$(job.Tw)" for job in failures], "; ",
        ),
    )
    missing = [
        output_path(case.model, case.Tw) for case in cases
        if !isfile(output_path(case.model, case.Tw))
    ]
    isempty(missing) || error("missing outputs: $(join(missing, ", "))")
    println("extended zero-frequency batch finished $(length(cases)) cases")
end

main()
