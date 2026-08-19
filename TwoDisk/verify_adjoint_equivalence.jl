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
const N_cheb = 129
const target = 0.46 - 0.014im

function cheb_quad(N::Int)
    x = [cos(pi * j / N) for j in 0:N]
    weights = zeros(N + 1)
    for j in 0:N
        series_sum = 0.0
        for k in 1:floor(Int, N / 2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2k * j * pi / N)
            series_sum += 2k == N ? 0.5 * term : term
        end
        endpoint_factor = (j == 0 || j == N) ? 1.0 : 2.0
        weights[j + 1] = endpoint_factor / N * (1.0 + series_sum)
    end
    jacobian = Diagonal(fill(0.5, N + 1))
    mass = kron(I(4), Diagonal(jacobian * weights))
    # Match the cavity-mode degrees of freedom removed by boudary_condition.
    constrained = (1, N + 1, N + 2, 2N + 2, 2N + 3, 3N + 3)
    free = setdiff(1:size(mass, 1), constrained)
    return mass[free, free]
end

function full_mode(vector)
    return CRC_STA.eig_full(reshape(vector, :, 1), N_cheb, 1)
end

function normalize_modes(direct, adjoint)
    direct_scale = maximum(abs.(full_mode(direct)[2]))
    adjoint_scale = maximum(abs.(full_mode(adjoint)[1]))
    return direct / direct_scale, adjoint / adjoint_scale
end

u0, v0, w0, _, _, _ = CRC_BF.BaseFlow(Re_h, Ro, a_s, mode)
D, D2, z = CRC_BF.Cheb(N_cheb, mode)
F, G, H = CRC_BF.interp(u0, v0, w0, z, N_cheb, mode)
cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, R, N_cheb, D, D2, Re_h)
L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(cof, D, D2, beta, omega, R)
L0, L1, L2 = CRC_STA.boudary_condition(L0_raw, L1_raw, L2_raw, N_cheb, mode)

direct_values, direct_vectors = iar(
    PEP([L0, L1, L2]);
    σ = target,
    neigs = 1,
    maxit = 500,
    tol = 1e-14,
)
alpha = direct_values[1]
direct = direct_vectors[:, 1]

# The transpose implementation used in the original notebook.
T0_raw = transpose(cof.D1) + im * R * beta * transpose(cof.B) -
    im * omega * transpose(cof.Ta) - beta^2 * R^2 * transpose(cof.Vyy) -
    transpose(cof.dC) - im * beta * transpose(cof.dVyz) +
    transpose(cof.d2Vzz) -
    (transpose(cof.C) + im * beta * R * transpose(cof.Vyz) -
     2transpose(cof.dVzz)) * kron(I(4), D) +
    transpose(cof.Vzz) * kron(I(4), D2)
T1_raw = im * transpose(cof.A) - beta * R * transpose(cof.Vxy) -
    im * transpose(cof.dVxz) - im * transpose(cof.Vxz) * kron(I(4), D)
T2_raw = -transpose(cof.Vxx)
T0, T1, T2 = CRC_STA.boudary_condition(
    T0_raw, T1_raw, T2_raw, N_cheb, mode
)
transpose_values, transpose_vectors = iar(
    PEP([T0, T1, T2]);
    σ = alpha,
    neigs = 1,
    maxit = 500,
    tol = 1e-14,
)
alpha_transpose = transpose_values[1]
adjoint_transpose = transpose_vectors[:, 1]

# Hermitian coefficients are the elementwise conjugates of the transpose
# coefficients. The paired eigenvalue and vector are therefore conjugated.
H0, H1, H2 = conj.(T0), conj.(T1), conj.(T2)
hermitian_values, hermitian_vectors = iar(
    PEP([H0, H1, H2]);
    σ = conj(alpha),
    neigs = 1,
    maxit = 500,
    tol = 1e-14,
)
alpha_hermitian = hermitian_values[1]
adjoint_hermitian = hermitian_vectors[:, 1]

direct, adjoint_transpose = normalize_modes(direct, adjoint_transpose)
_, adjoint_hermitian = normalize_modes(direct, adjoint_hermitian)

# Align the arbitrary global phase before comparing complex vectors.
phase = dot(conj.(adjoint_transpose), adjoint_hermitian)
adjoint_hermitian_aligned = adjoint_hermitian * exp(-im * angle(phase))
expected_hermitian = conj.(adjoint_transpose)
phase_expected = dot(expected_hermitian, adjoint_hermitian)
adjoint_hermitian_conjugate_aligned =
    adjoint_hermitian * exp(-im * angle(phase_expected))

mass = cheb_quad(N_cheb)
Q_transpose = transpose(adjoint_transpose) * mass *
    (L1 + (alpha + alpha_transpose) * L2) * direct
Q_hermitian = adjoint_hermitian' * mass *
    (L1 + (alpha + conj(alpha_hermitian)) * L2) * direct

vel_transpose = full_mode(adjoint_transpose)
vel_hermitian = full_mode(adjoint_hermitian)
height = 1 / Re_h
localization = 1 / 2
roughness_transform = height * exp(-alpha^2 / (4localization)) *
    sqrt(pi / localization)
u_wall = -(D * F)[1] * roughness_transform
v_wall = -(D * G)[1] * roughness_transform
BC_transpose = ((D * vel_transpose[1])[1] * u_wall +
    (D * vel_transpose[2])[1] * v_wall) / (R * sqrt(Re_h))
BC_hermitian = conj((D * vel_hermitian[1])[1]) * u_wall /
    (R * sqrt(Re_h)) + conj((D * vel_hermitian[2])[1]) * v_wall /
    (R * sqrt(Re_h))
Cr_transpose = abs(-im * BC_transpose / Q_transpose)
Cr_hermitian = abs(-im * BC_hermitian / Q_hermitian)

profile_magnitude_error = maximum(
    abs.(abs.(adjoint_hermitian) .- abs.(adjoint_transpose))
)
conjugate_vector_error = norm(
    adjoint_hermitian_conjugate_aligned - expected_hermitian
) / norm(expected_hermitian)

open("adjoint_equivalence_results.txt", "w") do io
    println(io, "alpha_direct = ", repr(alpha))
    println(io, "alpha_transpose = ", repr(alpha_transpose))
    println(io, "alpha_hermitian = ", repr(alpha_hermitian))
    println(io, "transpose_eigenvalue_error = ", abs(alpha_transpose - alpha))
    println(io, "hermitian_eigenvalue_error = ", abs(alpha_hermitian - conj(alpha)))
    println(io, "normalized_magnitude_profile_error = ", profile_magnitude_error)
    println(io, "conjugate_vector_error_after_phase_alignment = ", conjugate_vector_error)
    println(io, "Q_transpose = ", repr(Q_transpose))
    println(io, "Q_hermitian = ", repr(Q_hermitian))
    println(io, "BC_transpose = ", repr(BC_transpose))
    println(io, "BC_hermitian = ", repr(BC_hermitian))
    println(io, "Cr_transpose = ", Cr_transpose)
    println(io, "Cr_hermitian = ", Cr_hermitian)
    println(io, "Cr_absolute_difference = ", abs(Cr_transpose - Cr_hermitian))
end

print(read("adjoint_equivalence_results.txt", String))
