const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(WORKSPACE_ROOT,"NeutralCurveRunner.jl"))

using .NeutralCurveRunner

function main()
    config = CurveConfig(
        Tw=1.0,model=:lopez,N_cheb=69,num_modes=2,
        output_dir=joinpath(WORKSPACE_ROOT,"neutral_curve_batch"),
    )
    solver = NeutralCurveRunner.prepare_solver(config).solve_at
    benchmarks = (
        (name="Type-I",R=285.36,beta=0.07759,target=0.38),
        (name="Type-II",R=440.88,beta=0.04672,target=0.14),
    )
    for point in benchmarks
        values,_ = solver(point.R,point.beta,point.target,2)
        println(
            "$(point.name) R=$(point.R) beta=$(point.beta) " *
            "eigenvalues=$(values)",
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
