include(joinpath(@__DIR__, "test_bek2_configuration.jl"))

function report(label, result)
    println(label,
            " alpha=", result.alpha,
            " Cr_raw=", result.Cr,
            " Cr_Thomas=", result.Cr_thomas,
            " ||Q|-|K||/|K|=", abs(abs(result.Q) - abs(result.K)) / abs(result.K),
            " pairing_error=", result.pairing_error,
            " direct_residual=", result.direct_residual,
            " adjoint_residual=", result.adjoint_residual)
end

neutral_R = 285.36
neutral_beta = 0.07759
neutral = solve_point(R=neutral_R,
                      n=neutral_beta * neutral_R,
                      N=99,
                      shift=0.38482 + 0im)
report("Malik-neutral benchmark", neutral)

for N in (79, 99, 129)
    result = solve_point(R=470.0, n=30.0, N=N, shift=0.3356 - 0.0676im)
    report("R470 grid N=$(N)", result)
end
