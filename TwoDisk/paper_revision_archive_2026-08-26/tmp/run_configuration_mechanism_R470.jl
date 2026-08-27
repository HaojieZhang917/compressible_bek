include(joinpath(@__DIR__, "..", "configuration_effect_validation.jl"))

ctx = cavity_context(99)
cavity = cavity_point(470.0, 0.178 - 0.007im, nothing, ctx)
single = solve_point(R=470.0, n=30.0, N=99, shift=0.280 - 0.062im)

eta_cavity = vec(ctx.z) .* sqrt(ctx.cfg.Re_h)
eta_single = vec(single.grid.z)

function peak_location(eta, values)
    ids = findall(i -> isfinite(eta[i]), eachindex(eta))
    jlocal = argmax(abs.(values[ids]))
    j = ids[jlocal]
    return eta[j], abs(values[j])
end

function cancellation(a, b)
    return abs(a + b) / (abs(a) + abs(b))
end

Feta_cavity = (ctx.D * ctx.F)[1] / sqrt(ctx.cfg.Re_h)
Geta_cavity = (ctx.D * ctx.G)[1] / sqrt(ctx.cfg.Re_h)
Feta_single = (single.grid.D * single.baseflow.F)[1]
Geta_single = (single.grid.D * single.baseflow.G)[1]

metrics = [
    "R=470, n=30, omega_bar=0, a_s=0, N=99",
    "",
    "Base-flow wall gradients in the local viscous coordinate",
    "cavity F_eta(0)=$(Feta_cavity)",
    "cavity G_eta(0)=$(Geta_cavity)",
    "single F_eta(0)=$(Feta_single)",
    "single G_eta(0)=$(Geta_single)",
    "",
    "Eigenfunction peak locations",
    "cavity direct |v| peak=$(peak_location(eta_cavity, cavity.direct_fields[2]))",
    "single direct |v| peak=$(peak_location(eta_single, single.direct_fields[2]))",
    "cavity adjoint |u| peak=$(peak_location(eta_cavity, cavity.adjoint_fields[1]))",
    "single adjoint |u| peak=$(peak_location(eta_single, single.adjoint_fields[1]))",
    "",
    "Boundary projection split",
    "cavity |BC_u/Q|=$(abs(cavity.BC_u/cavity.Q))",
    "cavity |BC_v/Q|=$(abs(cavity.BC_v/cavity.Q))",
    "cavity |(BC_u+BC_v)/Q|=$(abs(cavity.BC/cavity.Q))",
    "cavity BC cancellation index=$(cancellation(cavity.BC_u,cavity.BC_v))",
    "single Thomas |BC_u/Q|=$(abs(single.BC_u_thomas/single.Q))",
    "single Thomas |BC_v/Q|=$(abs(single.BC_v_thomas/single.Q))",
    "single Thomas |(BC_u+BC_v)/Q|=$(abs(single.BC_thomas/single.Q))",
    "single BC cancellation index=$(cancellation(single.BC_u,single.BC_v))",
]

output = joinpath(OUTPUT_DIR, "mechanism_metrics_R470.txt")
open(output, "w") do io
    for line in metrics
        println(io, line)
    end
end
foreach(println, metrics)
println("Wrote ", output)
