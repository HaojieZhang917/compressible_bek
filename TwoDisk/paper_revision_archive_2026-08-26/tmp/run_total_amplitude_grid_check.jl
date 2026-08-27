include(joinpath(@__DIR__, "..", "configuration_effect_validation.jl"))

const SQRT_TWO_PI = sqrt(2pi)
const CAVITY_NEUTRAL = 287.82160889457725
const SINGLE_NEUTRAL = 313.043829663044

for N in (79, 99, 129)
    cavity_ctx = cavity_context(N)
    cavity = cavity_point(CAVITY_NEUTRAL, 0.44968 + 0im, nothing, cavity_ctx)
    single = solve_point(
        R=SINGLE_NEUTRAL, n=30.0, N=N, shift=0.48858 + 0im
    )
    println(
        "N=", N,
        " cavity_alpha=", repr(cavity.alpha),
        " cavity_Cr=", cavity.Cr / SQRT_TWO_PI,
        " cavity_pair=", cavity.pairing_error,
        " single_alpha=", repr(single.alpha),
        " single_Cr=", single.Cr_thomas,
        " single_pair=", single.pairing_error,
    )
end
