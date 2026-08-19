using LinearAlgebra
using NonlinearEigenproblems

include(joinpath(@__DIR__, "..", "..", "TwoDiskReceptivity", "src", "TwoDiskReceptivity.jl"))
using .TwoDiskReceptivity

function mapped_weights(N::Int)
    xi = [-cos(pi * j / N) for j in 0:N]
    weights = zeros(N + 1)
    for j in 0:N
        s = 0.0
        for k in 1:floor(Int, N / 2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2k * j * pi / N)
            s += 2k == N ? 0.5 * term : term
        end
        endpoint = (j == 0 || j == N) ? 1.0 : 2.0
        weights[j + 1] = endpoint / N * (1.0 + s)
    end

    a, b, c = 2.0, 0.6, 0.5
    jacobian = similar(xi)
    for i in eachindex(xi)
        x = xi[i]
        den = 1 - b*x - (1-b)*(x^3 + c*(1-x^2))
        jacobian[i] = 2a * (b + 3*(1-b)*x^2 - 2c*(1-b)*x) / den^2
    end
    # The rational map sends xi=1 to infinity. All velocity degrees of freedom
    # at this endpoint and the last pressure degree of freedom are constrained.
    jacobian[end] = 0.0

    return jacobian .* weights
end

function mapped_mass_matrix(N::Int)
    z_weights = mapped_weights(N)
    full = kron(I(4), Diagonal(z_weights))
    constrained = (1, N + 1, N + 2, 2N + 2, 2N + 3, 3N + 3, 4N + 4)
    free = setdiff(1:4*(N+1), constrained)
    return full[free, free]
end

function solve_point(; R=470.0, n=30.0, N=99, shift=0.35-0.01im)
    params = FlowParameters(1000, -1.0, 0.0, R, n, 0.0, N, 2)
    baseflow, grid = solve_baseflow(params)
    F, G, H = baseflow.F, baseflow.G, baseflow.H
    D, D2 = grid.D, grid.D2
    beta = n / R
    omega = 0.0
    Ro = -1.0
    Co = 2 - Ro - Ro^2
    cof = assemble_coeff_matrix_BEK2(F, G, H, R, N, D, D2, Ro, Co)
    L0, L1, L2 = assemble_direct_matrices(cof, D, D2, beta, omega, R)
    L0, L1, L2 = apply_boundary_conditions!(L0, L1, L2, N)
    A0, A1, A2 = assemble_adjoint_matrices(cof, D, D2, beta, omega, R)
    A0, A1, A2 = apply_boundary_conditions!(A0, A1, A2, N)

    alpha, x = solve_eigenvalue_problem(L0, L1, L2, ComplexF64(shift), 1)
    alpha_a, y = solve_eigenvalue_problem(A0, A1, A2, ComplexF64(alpha[1]), 1)
    xn, yn, velocity, velocity_a = normalise_eigenvectors!(x, y, N)
    W = mapped_mass_matrix(N)
    Q = transpose(yn[:, 1]) * W * (L1 + (alpha[1] + alpha_a[1]) * L2) * xn[:, 1]

    u, v, w, p = velocity
    ua, va, wa, pa = velocity_a
    velocity_pair = u .* ua .+ v .* va .+ w .* wa
    K_integrand = vec(F) .* velocity_pair .+ p .* ua .+ u .* pa .-
                  (im / R) * (alpha[1] + alpha_a[1]) .* velocity_pair
    K = sum(mapped_weights(N) .* K_integrand)

    ls = 0.5
    hhat = exp(-alpha[1]^2 / (4ls)) * sqrt(pi / ls)
    u_wall = -(D * F)[1] * hhat
    v_wall = -(D * G)[1] * hhat
    BC_u = (D * velocity_a[1])[1] * u_wall / R
    BC_v = (D * velocity_a[2])[1] * v_wall / R
    BC = BC_u + BC_v
    Cr = abs(-im * BC / Q)
    Cr_K = abs(BC / K)
    transform_factor = sqrt(pi / ls)
    Cr_thomas = Cr / transform_factor
    BC_thomas = BC / transform_factor
    pairing_error = abs(alpha_a[1] - alpha[1])
    direct_residual = norm((L0 + alpha[1] * L1 + alpha[1]^2 * L2) * xn[:, 1]) /
                      (norm(L0) + abs(alpha[1]) * norm(L1) + abs(alpha[1])^2 * norm(L2)) /
                      norm(xn[:, 1])
    adjoint_residual = norm((A0 + alpha_a[1] * A1 + alpha_a[1]^2 * A2) * yn[:, 1]) /
                       (norm(A0) + abs(alpha_a[1]) * norm(A1) + abs(alpha_a[1])^2 * norm(A2)) /
                       norm(yn[:, 1])
    return (; alpha=alpha[1], alpha_a=alpha_a[1], Cr, Cr_K, Cr_thomas,
            Q, K, BC, BC_u, BC_v, BC_thomas,
            BC_u_thomas=BC_u/transform_factor,
            BC_v_thomas=BC_v/transform_factor,
            pairing_error, direct_residual,
            adjoint_residual, wall_shear=(D*F)[1] + im*(D*G)[1],
            baseflow=(F=F, G=G, H=H), grid, direct_fields=velocity,
            adjoint_fields=velocity_a)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for shift in (0.18-0.01im, 0.25-0.01im, 0.35-0.01im, 0.45-0.01im, 0.55-0.01im)
        try
            result = solve_point(shift=shift)
            println("shift=", shift,
                    " alpha=", result.alpha,
                    " alpha_a=", result.alpha_a,
                    " pairerr=", result.pairing_error,
                    " Cr_Q=", result.Cr,
                    " Cr_K=", result.Cr_K,
                    " |Q|=", abs(result.Q),
                    " |K|=", abs(result.K),
                    " |BC|=", abs(result.BC))
        catch err
            println("shift=", shift, " failed: ", sprint(showerror, err))
        end
    end

    for N in (79, 99, 129)
        result = solve_point(N=N, shift=0.35-0.06im)
        println("grid N=", N,
                " alpha=", result.alpha,
                " Cr_Q=", result.Cr,
                " Cr_Q_over_pi=", result.Cr / pi,
                " pairerr=", result.pairing_error)
    end
end
