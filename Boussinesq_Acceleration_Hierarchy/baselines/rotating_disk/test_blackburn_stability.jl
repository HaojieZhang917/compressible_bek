using Test
using LinearAlgebra

include(joinpath(@__DIR__, "CRD_STA.jl"))
include(joinpath(
    @__DIR__, "..", "..", "..", "RotatingDiskFlow", "src",
    "LopezStability.jl",
))
include(joinpath(@__DIR__, "BlackburnStability.jl"))

const BS = BlackburnStability
const LS = LopezStability
const ComplexScalar = ComplexF64

function split_fields(q, n)
    return (
        q[1:n], q[n+1:2n], q[2n+1:3n],
        q[3n+1:4n], q[4n+1:5n],
    )
end

function nonlinear_fourier_residual(
    epsilon, q, F, GL, H, T, Pr, D, D2, R, alpha, beta, omega_i,
)
    n = length(F)
    u, v, w, theta, p = split_fields(q, n)
    Q = 1 .- GL
    U = F .+ epsilon .* u
    V = Q .+ epsilon .* v
    W = H ./ R .+ epsilon .* w
    temperature = T .+ epsilon .* theta
    chi = 2 .- temperature

    Ur = F ./ R .+ epsilon .* (im * alpha .* u)
    Utheta = epsilon .* (im * beta .* u)
    Uz = D * F .+ epsilon .* (D * u)
    Vr = Q ./ R .+ epsilon .* (im * alpha .* v)
    Vtheta = epsilon .* (im * beta .* v)
    Vz = D * Q .+ epsilon .* (D * v)
    Wr = epsilon .* (im * alpha .* w)
    Wtheta = epsilon .* (im * beta .* w)
    Wz = (D * H) ./ R .+ epsilon .* (D * w)

    acceleration_r = U .* Ur .+ V .* Utheta .+ W .* Uz .- V.^2 ./ R
    acceleration_theta = U .* Vr .+ V .* Vtheta .+ W .* Vz .+ U .* V ./ R
    acceleration_z = U .* Wr .+ V .* Wtheta .+ W .* Wz
    base_acceleration_z = H .* (D * H) ./ R^2
    k2 = alpha^2 + beta^2
    lap_u = D2 * F .+ epsilon .* (D2 * u .- k2 .* u)
    lap_v = D2 * Q .+ epsilon .* (D2 * v .- k2 .* v)
    lap_w = (D2 * H) ./ R .+ epsilon .* (D2 * w .- k2 .* w)

    residual_u = -im * omega_i * epsilon .* chi .* u .+ chi .* acceleration_r .+
                 epsilon .* (im * alpha .* p) .- lap_u ./ R
    residual_v = -im * omega_i * epsilon .* chi .* v .+ chi .* acceleration_theta .+
                 epsilon .* (im * beta .* p) .- lap_v ./ R
    # Malik's local parallel truncation drops the O(R^-2) density feedback
    # from the base axial acceleration while retaining all O(R^-1) terms.
    residual_w = -im * omega_i * epsilon .* chi .* w .+ chi .* acceleration_z .+
                 epsilon .* theta .* base_acceleration_z .+
                 epsilon .* (D * p) .- lap_w ./ R

    temperature_z = D * T .+ epsilon .* (D * theta)
    thermal_advection = U .* (epsilon .* im * alpha .* theta) .+
                        V .* (epsilon .* im * beta .* theta) .+
                        W .* temperature_z
    lap_temperature = D2 * T .+ epsilon .* (D2 * theta .- k2 .* theta)
    residual_temperature = -im * omega_i * epsilon .* theta .+
        thermal_advection .- lap_temperature ./ (Pr * R)
    residual_continuity = Ur .+ U ./ R .+ epsilon .* (im * beta .* v) .+ D * W

    return vcat(
        residual_u, residual_v, residual_w,
        residual_temperature, residual_continuity,
    )
end

function finite_difference_jacobian(
    F, GL, H, T, Pr, D, D2, R, alpha, beta, omega_i; step=2.0e-6,
)
    nstate = 5length(F)
    jacobian = zeros(ComplexScalar, nstate, nstate)
    direction = zeros(ComplexScalar, nstate)
    for column in 1:nstate
        fill!(direction, 0)
        direction[column] = 1
        plus = nonlinear_fourier_residual(
            step, direction, F, GL, H, T, Pr, D, D2,
            R, alpha, beta, omega_i,
        )
        minus = nonlinear_fourier_residual(
            -step, direction, F, GL, H, T, Pr, D, D2,
            R, alpha, beta, omega_i,
        )
        jacobian[:, column] .= (plus .- minus) ./ (2step)
    end
    return jacobian
end

relative_error(left, right) = norm(left - right) / max(norm(right), eps(Float64))

@testset "Blackburn operator linearization" begin
    n = 6
    D, D2, z_matrix = CRD_BF.Cheb(n - 1; domain=:finite, ymax=3.0)
    z = vec(z_matrix)
    F = 0.18 .* exp.(-0.7 .* z) .- 0.025 .* exp.(-2.0 .* z)
    GL = 1 .- exp.(-0.9 .* z)
    H = -0.72 .* (1 .- exp.(-0.6 .* z))
    T = 1 .+ 0.12 .* exp.(-0.8 .* z)
    Pr = 0.72
    R = 285.36
    beta = 0.07759
    omega_i = 0.031

    L0, L1, L2, metadata = BS.blackburn_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_i; frame=:inertial,
    )
    @test metadata.acceleration_model == :blackburn2021
    J0 = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, 0.0, beta, omega_i,
    )
    Jplus = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, 1.0, beta, omega_i,
    )
    Jminus = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, -1.0, beta, omega_i,
    )
    @test relative_error(L0, J0) < 2.0e-9
    @test relative_error(L1, (Jplus - Jminus) ./ 2) < 2.0e-9
    @test relative_error(L2, (Jplus + Jminus) ./ 2 - J0) < 2.0e-9

    L0_lopez, L1_lopez, L2_lopez, _ = LS.lopez_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_i; frame=:inertial,
    )
    chi = 2 .- T
    expected_frequency_correction = zeros(ComplexScalar, 5n, 5n)
    for offset in (0, n, 2n)
        rows = offset+1:offset+n
        expected_frequency_correction[rows, rows] .=
            -im * omega_i .* Matrix(Diagonal(ComplexF64.(chi .- 1)))
    end
    @test relative_error(
        L0 - L0_lopez,
        expected_frequency_correction,
    ) < 2.0e-14
    @test L1 == L1_lopez
    @test L2 == L2_lopez

    u_rows = 1:n
    v_rows = n+1:2n
    w_rows = 2n+1:3n
    theta_columns = 3n+1:4n
    @test L0[u_rows, theta_columns] == Matrix(Diagonal(ComplexF64.(-(
        F.^2 .+ H .* (D * F) .- (1 .- GL).^2
    ) ./ R)))
    @test L0[v_rows, theta_columns] == Matrix(Diagonal(ComplexF64.(-(
        2 .* F .* (1 .- GL) .+ H .* (D * (1 .- GL))
    ) ./ R)))
    @test all(iszero, L0[w_rows, theta_columns])

    alpha = 0.38482 + 0.002im
    omega_disk = -0.014 + 0.003im
    L0d, L1d, L2d, _ = BS.blackburn_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_disk; frame=:disk,
    )
    L0i, L1i, L2i, _ = BS.blackburn_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_disk + beta;
        frame=:inertial,
    )
    @test L0d == L0i
    @test L1d == L1i
    @test L2d == L2i

    A, B, _ = BS.blackburn_temporal_matrices(
        F, GL, H, T, Pr, D, D2, R, alpha, beta; frame=:disk,
    )
    @test relative_error(A - omega_disk .* B, L0d + alpha .* L1d + alpha^2 .* L2d) < 5.0e-15
    chi_mass = im .* Matrix(Diagonal(ComplexF64.(chi)))
    @test B[1:n, 1:n] == chi_mass
    @test B[n+1:2n, n+1:2n] == chi_mass
    @test B[2n+1:3n, 2n+1:3n] == chi_mass
    @test B[3n+1:4n, 3n+1:4n] == im .* Matrix{ComplexF64}(I, n, n)
    @test all(iszero, B[4n+1:5n, :])
    @test all(iszero, B[:, 4n+1:5n])

    T_isothermal = ones(n)
    blackburn_iso = BS.blackburn_spatial_matrices(
        F, GL, H, T_isothermal, Pr, D, D2, R, beta, omega_i;
        frame=:inertial,
    )
    lopez_iso = LS.lopez_spatial_matrices(
        F, GL, H, T_isothermal, Pr, D, D2, R, beta, omega_i;
        frame=:inertial,
    )
    @test blackburn_iso[1] == lopez_iso[1]
    @test blackburn_iso[2] == lopez_iso[2]
    @test blackburn_iso[3] == lopez_iso[3]

    reduced0, reduced1, reduced2, keep = BS.apply_homogeneous_boundaries(
        (L0, L1, L2), n,
    )
    @test size(reduced0) == (5n - 9, 5n - 9)
    @test size(reduced1) == size(reduced0)
    @test size(reduced2) == size(reduced0)
    @test length(keep) == 5n - 9
    @test 4n + 1 in keep
    @test !(5n in keep)
end
