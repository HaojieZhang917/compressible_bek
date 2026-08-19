using LinearAlgebra
using BSplineKit
using DelimitedFiles
using PyCall
using NonlinearEigenproblems

include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")

const Re_h = 2000
const Ro = -1.0
const a_s = 0.0
const mode = 1
const R = 330
const beta = 0.0503
const omega_bar = 8.0
const omega = omega_bar / R

# The base-flow boundary-value problem is solved once on its native fine grid.
# Only the Chebyshev resolution of the stability problem is varied below.
u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(Re_h, Ro, a_s, mode)

function direct_eigenvalue(N_z, target)
    N_cheb = N_z - 1
    D, D2, z = CRC_BF.Cheb(N_cheb, mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, N_cheb, mode)
    cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, R, N_cheb, D, D2, Re_h)
    L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(cof, D, D2, beta, omega, R)
    L0, L1, L2 = CRC_STA.boudary_condition(
        L0_raw, L1_raw, L2_raw, N_cheb, mode
    )
    problem = PEP([L0, L1, L2])
    values, _ = iar(
        problem;
        σ = target,
        neigs = 1,
        maxit = 500,
        tol = 1e-14,
    )
    return values[1]
end

N_values = [40, 60, 80, 100, 150, 200, 250, 300]
eigenvalues = ComplexF64[]
target = 0.46 - 0.014im

for N_z in N_values
    alpha = direct_eigenvalue(N_z, target)
    push!(eigenvalues, alpha)
    global target = alpha
    println("N_z=$(N_z), alpha=$(repr(alpha))")
    flush(stdout)
end

alpha_ref = eigenvalues[end]
relative_errors = abs.((eigenvalues .- alpha_ref) ./ alpha_ref)
output = hcat(N_values, real.(eigenvalues), imag.(eigenvalues), relative_errors)
writedlm("grid_convergence_results.dat", output)

println("reference alpha=$(repr(alpha_ref))")
for (N_z, alpha, error) in zip(N_values, eigenvalues, relative_errors)
    println(
        "RESULT N_z=$(N_z) alpha_r=$(real(alpha)) " *
        "alpha_i=$(imag(alpha)) relative_error=$(error)"
    )
end
