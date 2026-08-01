using Printf

const GRID_SIZES = [49, 59, 69, 79, 89]
const WALL_TEMPERATURES = [1.12, 1.13, 1.14, 1.16]

function curve_path(output_root, N_cheb, Tw)
    Tw_text = replace(@sprintf("%.4f", Tw), r"0+$" => "")
    endswith(Tw_text, ".") && (Tw_text *= "0")
    filename = (
        "ome=0.0_Tw=$(Tw_text)_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat"
    )
    return joinpath(output_root, "N=$(N_cheb)", filename)
end

function completed_curve(path)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) >= 42 || return false
    fields = split(strip(lines[end]))
    length(fields) == 7 || return false
    beta = tryparse(Float64, fields[3])
    return beta !== nothing && beta >= 0.0465
end

function launch_case(output_root, log_dir, N_cheb, Tw)
    log_path = joinpath(
        log_dir,
        @sprintf("N=%d_Tw=%.2f.log", N_cheb, Tw),
    )
    stream = open(log_path, "w")
    environment = copy(ENV)
    environment["TYPEII_GRID_TW"] = string(Tw)
    environment["TYPEII_GRID_N_CHEB"] = string(N_cheb)
    environment["TYPEII_GRID_OUTPUT_ROOT"] = output_root
    command = `$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(
        @__DIR__, "ComputeTypeIIGridConvergence.jl",
    ))`
    process = run(
        pipeline(setenv(command, environment); stdout=stream, stderr=stream);
        wait=false,
    )
    println("launched N=$(N_cheb) Tw=$(Tw) pid=$(getpid(process))")
    flush(stdout)
    return (
        process=process,
        stream=stream,
        log_path=log_path,
        N_cheb=N_cheb,
        Tw=Tw,
    )
end

function finish_case(job)
    wait(job.process)
    close(job.stream)
    println(
        "finished N=$(job.N_cheb) Tw=$(job.Tw) " *
        "exitcode=$(job.process.exitcode) log=$(job.log_path)",
    )
    flush(stdout)
    return job.process.exitcode == 0
end

function main()
    output_root = get(
        ENV,
        "TYPEII_GRID_OUTPUT_ROOT",
        joinpath(@__DIR__, "typeII_grid_convergence_results", "raw"),
    )
    max_workers = parse(Int, get(ENV, "TYPEII_GRID_WORKERS", "4"))
    max_workers >= 1 || error("TYPEII_GRID_WORKERS must be positive")
    log_dir = joinpath(dirname(output_root), "process_logs")
    mkpath(output_root)
    mkpath(log_dir)

    cases = [(N_cheb=N, Tw=Tw) for N in GRID_SIZES for Tw in WALL_TEMPERATURES]
    active = NamedTuple[]
    failures = NamedTuple[]
    skipped = 0

    for case in cases
        path = curve_path(output_root, case.N_cheb, case.Tw)
        if completed_curve(path)
            println("skip completed N=$(case.N_cheb) Tw=$(case.Tw)")
            skipped += 1
            continue
        end
        while length(active) >= max_workers
            job = popfirst!(active)
            finish_case(job) || push!(failures, job)
        end
        push!(active, launch_case(
            output_root, log_dir, case.N_cheb, case.Tw,
        ))
    end

    for job in active
        finish_case(job) || push!(failures, job)
    end

    println(
        "grid batch complete: total=$(length(cases)) skipped=$(skipped) " *
        "failures=$(length(failures))",
    )
    isempty(failures) || error(
        "failed cases: " *
        join(["N=$(job.N_cheb),Tw=$(job.Tw)" for job in failures], "; "),
    )
end

main()
