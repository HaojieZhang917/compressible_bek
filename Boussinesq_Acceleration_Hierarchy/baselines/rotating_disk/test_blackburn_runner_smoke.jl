include(joinpath(@__DIR__, "BlackburnNeutralCurveRunner.jl"))

using .BlackburnNeutralCurveRunner

config = CurveConfig(
    Tw=1.0,
    N_cheb=39,
    model=:blackburn,
)
prepared = BlackburnNeutralCurveRunner.prepare_solver(config)
values, vectors = prepared.solve_at(500.0, 0.04, 0.1, 2)

@assert length(values) == 2
@assert size(vectors, 2) == 2
@assert all(isfinite, real.(values))
@assert all(isfinite, imag.(values))

println("Blackburn runner smoke test values: ", values)
