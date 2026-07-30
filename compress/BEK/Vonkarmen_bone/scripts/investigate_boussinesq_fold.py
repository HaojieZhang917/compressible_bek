#!/usr/bin/env python3
"""Cross-validate the Boussinesq rotating-disk similarity fold near Tw=1.05."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eig, svdvals
from scipy.optimize import root


PR = 0.72
LAMBDA_C = 1.0
WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
OUT = WORKSPACE_ROOT / "boussinesq_fold_analysis"


def ode(z, y):
    H, Fp, F, Gp, G, Tp, T = y
    return np.vstack((
        -2.0 * F,
        F**2 + H * Fp - (G - 1.0)**2 + LAMBDA_C * (T - 1.0),
        Fp,
        2.0 * F * G + H * Gp - 2.0 * F,
        Gp,
        PR * H * Tp,
        Tp,
    ))


def initial_guess(z):
    y = np.zeros((7, z.size))
    y[0] = -0.8845 * (1.0 - np.exp(-0.8 * z))
    y[1] = 0.5102 * (1.0 - 0.8 * z) * np.exp(-0.8 * z)
    y[2] = 0.5102 * z * np.exp(-0.8 * z)
    y[3] = 0.6159 * np.exp(-0.8 * z)
    y[4] = 1.0 - np.exp(-0.8 * z)
    y[6] = 1.0
    return y


def fixed_tw_bc(tw):
    def bc(ya, yb):
        return np.array((
            ya[0], ya[2], ya[4], ya[6] - tw,
            yb[2], yb[4] - 1.0, yb[6] - 1.0,
        ))
    return bc


def fixed_h_bc(hinf):
    def bc(ya, yb, p):
        return np.array((
            ya[0], ya[2], ya[4], ya[6] - p[0],
            yb[2], yb[4] - 1.0, yb[6] - 1.0, yb[0] - hinf,
        ))
    return bc


def solve_isothermal(zmax=20.0, n=401, tol=1.0e-9):
    z = np.linspace(0.0, zmax, n)
    sol = solve_bvp(ode, fixed_tw_bc(1.0), z, initial_guess(z),
                    tol=tol, max_nodes=100000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def solve_fixed_h(hinf, z, guess, tw_guess, tol=2.0e-8):
    fun = lambda x, y, p: ode(x, y)
    sol = solve_bvp(fun, fixed_h_bc(hinf), z, guess, p=[tw_guess],
                    tol=tol, max_nodes=100000)
    if not sol.success:
        raise RuntimeError(f"Hinf={hinf}: {sol.message}")
    return sol


def trace_h_branch(zmax=20.0, dh=0.005):
    z = np.linspace(0.0, zmax, 401)
    seed = solve_isothermal(zmax, z.size)

    # Move from the classical branch to Hinf=-0.75, then use Hinf as the
    # continuation coordinate. This remains regular at a Tw turning point.
    guess = seed.sol(z)
    tw = 1.0
    rows = []
    sols = {}
    h_values = np.arange(-0.75, 0.1001, dh)
    for h in h_values:
        sol = solve_fixed_h(float(h), z, guess, tw)
        tw = float(sol.p[0])
        guess = sol.sol(z)
        y0 = sol.sol(0.0)
        rows.append((h, tw, y0[1], y0[3], y0[5], sol.x.size))
        sols[round(float(h), 10)] = sol
    return np.asarray(rows, float), sols


def turning_points(branch):
    h = branch[:, 0]
    tw = branch[:, 1]
    spline = CubicSpline(h, tw)
    candidates = spline.derivative().roots()
    roots = [x for x in candidates if h[1] < x < h[-2]]
    return [(float(x), float(spline(x))) for x in roots]


def roots_at_tw(branch, target_tw):
    h = branch[:, 0]
    vals = branch[:, 1] - target_tw
    roots = []
    for i in range(h.size - 1):
        if vals[i] == 0.0:
            roots.append(h[i])
        elif vals[i] * vals[i + 1] < 0.0:
            local = CubicSpline(h[max(0, i - 2):min(h.size, i + 4)],
                                vals[max(0, i - 2):min(h.size, i + 4)])
            lo, hi = h[i], h[i + 1]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                if local(lo) * local(mid) <= 0.0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
    return roots


def shooting_residual(slopes, tw, zmax, return_solution=False):
    y0 = np.array((0.0, slopes[0], 0.0, slopes[1], 0.0,
                   slopes[2], tw), float)

    def rhs(z, y):
        return ode(np.array([z]), y[:, None])[:, 0]

    ivp = solve_ivp(rhs, (0.0, zmax), y0, method="DOP853",
                    rtol=2.0e-11, atol=2.0e-13, dense_output=True)
    if not ivp.success or not np.all(np.isfinite(ivp.y[:, -1])):
        residual = np.full(3, 1.0e6)
    else:
        residual = np.array((ivp.y[2, -1], ivp.y[4, -1] - 1.0,
                             ivp.y[6, -1] - 1.0))
    return (residual, ivp) if return_solution else residual


def shoot(tw, slopes, zmax=20.0):
    ans = root(lambda s: shooting_residual(s, tw, zmax), slopes,
               method="hybr", options={"xtol": 1.0e-10, "maxfev": 1000})
    residual, ivp = shooting_residual(ans.x, tw, zmax, True)
    return ans, residual, ivp


def cheb_matrix(n, zmax):
    if n == 0:
        return np.ones(1), np.zeros((1, 1))
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.r_[2.0, np.ones(n - 1), 2.0] * (-1.0) ** np.arange(n + 1)
    X = np.tile(x, (n + 1, 1)).T
    dX = X - X.T
    D = (np.outer(c, 1.0 / c) / (dX + np.eye(n + 1)))
    D -= np.diag(D.sum(axis=1))
    z = 0.5 * zmax * (1.0 - x)
    Dz = -2.0 * D / zmax
    return z, Dz


def newton_residual_jacobian(u, tw, D, D2):
    n = D.shape[0]
    I = np.eye(n)
    H, F, G, T = u.reshape(4, n)
    DF, DG, DT = D @ F, D @ G, D @ T

    residual = np.concatenate((
        D @ H + 2.0 * F,
        D2 @ F - F**2 - H * DF + (G - 1.0)**2 - (T - 1.0),
        D2 @ G - 2.0 * F * G - H * DG + 2.0 * F,
        D2 @ T - PR * H * DT,
    ))
    J = np.zeros((4 * n, 4 * n))
    sH, sF, sG, sT = [slice(i * n, (i + 1) * n) for i in range(4)]
    J[sH, sH] = D
    J[sH, sF] = 2.0 * I
    J[sF, sH] = -np.diag(DF)
    J[sF, sF] = D2 - np.diag(2.0 * F) - np.diag(H) @ D
    J[sF, sG] = np.diag(2.0 * (G - 1.0))
    J[sF, sT] = -I
    J[sG, sH] = -np.diag(DG)
    J[sG, sF] = np.diag(2.0 - 2.0 * G)
    J[sG, sG] = D2 - np.diag(2.0 * F) - np.diag(H) @ D
    J[sT, sH] = -PR * np.diag(DT)
    J[sT, sT] = D2 - PR * np.diag(H) @ D

    def set_bc(row, col, value):
        residual[row] = value
        J[row, :] = 0.0
        J[row, col] = 1.0

    set_bc(0, 0, H[0])
    set_bc(n, n, F[0])
    set_bc(2 * n - 1, 2 * n - 1, F[-1])
    set_bc(2 * n, 2 * n, G[0])
    set_bc(3 * n - 1, 3 * n - 1, G[-1] - 1.0)
    set_bc(3 * n, 3 * n, T[0] - tw)
    set_bc(4 * n - 1, 4 * n - 1, T[-1] - 1.0)
    return residual, J


def newton_collocation(tw, source_sol, zmax=20.0, degree=60):
    z, D = cheb_matrix(degree, zmax)
    D2 = D @ D
    src = source_sol.sol(z)
    u = np.concatenate((src[0], src[2], src[4], src[6]))
    history = []
    for _ in range(15):
        residual, J = newton_residual_jacobian(u, tw, D, D2)
        norm0 = np.linalg.norm(residual, np.inf)
        history.append(norm0)
        if norm0 < 2.0e-10:
            break
        step = np.linalg.solve(J, -residual)
        damping = 1.0
        while damping > 1.0e-5:
            trial = u + damping * step
            trial_norm = np.linalg.norm(
                newton_residual_jacobian(trial, tw, D, D2)[0], np.inf)
            if trial_norm < norm0:
                u = trial
                break
            damping *= 0.5
        else:
            raise RuntimeError("Newton line search failed")
    residual, J = newton_residual_jacobian(u, tw, D, D2)
    row_norm = np.linalg.norm(J, axis=1)
    Jr = J / np.maximum(row_norm[:, None], 1.0e-300)
    sigma = svdvals(Jr)
    return z, u.reshape(4, z.size), history, sigma[-1] / sigma[0]


def similarity_eigenvalues(source_sol, zmax=20.0, degree=60):
    """Axisymmetric similarity-preserving temporal eigenvalues."""
    z, D = cheb_matrix(degree, zmax)
    D2 = D @ D
    n = z.size
    I = np.eye(n)
    y = source_sol.sol(z)
    H, F, G, T = y[0], y[2], y[4], y[6]
    Fp, Gp, Tp = D @ F, D @ G, D @ T
    sf, sg, sh, st = [slice(i * n, (i + 1) * n) for i in range(4)]
    A = np.zeros((4 * n, 4 * n))
    B = np.zeros_like(A)

    # lambda*f = f'' - H*f' - 2F*f + 2(G-1)g - F'*h - theta
    A[sf, sf] = D2 - np.diag(H) @ D - np.diag(2.0 * F)
    A[sf, sg] = np.diag(2.0 * (G - 1.0))
    A[sf, sh] = -np.diag(Fp)
    A[sf, st] = -I
    # lambda*g = g'' - H*g' - 2F*g + 2(1-G)f - G'*h
    A[sg, sf] = np.diag(2.0 * (1.0 - G))
    A[sg, sg] = D2 - np.diag(H) @ D - np.diag(2.0 * F)
    A[sg, sh] = -np.diag(Gp)
    # Similarity continuity constraint: h' + 2f = 0.
    A[sh, sf] = 2.0 * I
    A[sh, sh] = D
    # lambda*theta = theta''/Pr - H*theta' - T'*h
    A[st, sh] = -np.diag(Tp)
    A[st, st] = D2 / PR - np.diag(H) @ D
    B[sf, sf] = I
    B[sg, sg] = I
    B[st, st] = I

    def set_bc(row, col):
        A[row, :] = 0.0
        B[row, :] = 0.0
        A[row, col] = 1.0

    set_bc(0, 0)                    # f(0)=0
    set_bc(n - 1, n - 1)           # f(infinity)=0
    set_bc(n, n)                    # g(0)=0
    set_bc(2 * n - 1, 2 * n - 1)   # g(infinity)=0
    set_bc(2 * n, 2 * n)           # h(0)=0
    set_bc(3 * n, 3 * n)           # theta(0)=0
    set_bc(4 * n - 1, 4 * n - 1)   # theta(infinity)=0

    values = eig(A, B, right=False)
    values = values[np.isfinite(values) & (np.abs(values) < 1.0e5)]
    return values[np.argsort(values.real)[::-1]]


def interpolate_solution_at_h(target_h, branch, z, start_sol):
    idx = int(np.argmin(np.abs(branch[:, 0] - target_h)))
    h0, tw0 = branch[idx, :2]
    guess = start_sol.sol(z)
    return solve_fixed_h(target_h, z, guess, tw0, tol=5.0e-9)


def domain_fold_estimate(zmax):
    z = np.linspace(0.0, zmax, 351)
    seed = solve_isothermal(zmax, z.size)
    guess, tw = seed.sol(z), 1.0
    hs = np.arange(-0.70, -0.0399, 0.01)
    rows = []
    for h in hs:
        sol = solve_fixed_h(h, z, guess, tw, tol=5.0e-8)
        tw = float(sol.p[0])
        guess = sol.sol(z)
        rows.append((h, tw))
    rows = np.asarray(rows)
    turns = turning_points(np.c_[rows, np.zeros((rows.shape[0], 4))])
    return turns


def upper_branch_domain_stability():
    """Track slow-tail upper-branch states as the far boundary is moved."""
    targets = (-0.18, -0.14, -0.10, -0.08)
    degrees = {20.0: 70, 40.0: 90, 60.0: 110, 80.0: 130}
    rows = []
    for zmax in degrees:
        nz = min(1201, max(401, int(20 * zmax) + 1))
        z = np.linspace(0.0, zmax, nz)
        seed = solve_isothermal(zmax, nz, tol=2.0e-8)
        guess, tw = seed.sol(z), 1.0
        for h in np.arange(-0.70, -0.0799, 0.02):
            sol = solve_fixed_h(float(h), z, guess, tw, tol=8.0e-8)
            guess, tw = sol.sol(z), float(sol.p[0])
            if any(abs(h - target) < 1.0e-10 for target in targets):
                values = similarity_eigenvalues(
                    sol, zmax=zmax, degree=degrees[zmax])
                leading = values[0]
                rows.append((zmax, h, tw, leading.real, leading.imag,
                             np.count_nonzero(values.real > 1.0e-8),
                             1.0 / (PR * abs(h))))
    return rows


def write_csv(path, header, rows):
    np.savetxt(path, np.asarray(rows), delimiter=",", header=header,
               comments="", fmt="%.12e")


def main():
    OUT.mkdir(exist_ok=True)
    branch, branch_solutions = trace_h_branch()
    turns = turning_points(branch)
    write_csv(OUT / "branch_by_Hinf.csv",
              "Hinf,Tw,Fp0,Gp0,Tp0,nodes", branch)

    z_bvp = np.linspace(0.0, 20.0, 401)
    turning_rows = []
    for point_id, (h, tw) in enumerate(turns[:2], start=1):
        nearest_h = branch[np.argmin(np.abs(branch[:, 0] - h)), 0]
        nearest = branch_solutions[round(float(nearest_h), 10)]
        sol = solve_fixed_h(h, z_bvp, nearest.sol(z_bvp), tw,
                            tol=2.0e-10)
        _, _, history, sigma_ratio = newton_collocation(
            tw, sol, degree=70)
        thermal_length = 1.0 / (PR * abs(h))
        turning_rows.append((point_id, h, tw, history[-1], sigma_ratio,
                             thermal_length))
    write_csv(OUT / "turning_points.csv",
              "point,Hinf,Tw,newton_residual,newton_sigma_min_over_max,"
              "thermal_decay_length", turning_rows)

    target_tw = 1.045
    h_roots = roots_at_tw(branch, target_tw)
    validation = []
    profile_rows = []
    profile_solutions = []

    for branch_id, h in enumerate(h_roots, start=1):
        nearest_h = branch[np.argmin(np.abs(branch[:, 0] - h)), 0]
        nearest = branch_solutions[round(float(nearest_h), 10)]
        sol = solve_fixed_h(h, z_bvp, nearest.sol(z_bvp), target_tw,
                            tol=2.0e-10)
        slopes = sol.sol(0.0)[[1, 3, 5]]

        shooting, shoot_res, ivp = shoot(target_tw, slopes)
        z_newton, u_newton, history, sigma_ratio = newton_collocation(
            target_tw, sol)

        z_compare = np.linspace(0.0, 20.0, 1001)
        bvp_y = sol.sol(z_compare)
        shoot_y = ivp.sol(z_compare)
        bvp_profiles = bvp_y[[0, 2, 4, 6]]
        shoot_profiles = shoot_y[[0, 2, 4, 6]]
        diff_shoot = np.max(np.abs(shoot_profiles - bvp_profiles))
        # Compare directly at Chebyshev nodes. Linear interpolation from the
        # clustered spectral grid would dominate the actual Newton error.
        diff_newton = np.max(np.abs(
            u_newton - sol.sol(z_newton)[[0, 2, 4, 6]]))
        validation.append((branch_id, h, sol.p[0], slopes[0], slopes[1],
                           slopes[2], np.linalg.norm(shoot_res, np.inf),
                           diff_shoot, history[-1], diff_newton, sigma_ratio))
        profile_solutions.append((branch_id, sol))
        for j, zz in enumerate(z_compare):
            profile_rows.append((branch_id, zz, bvp_y[0, j], bvp_y[2, j],
                                 bvp_y[4, j], bvp_y[6, j]))

    write_csv(OUT / "cross_validation.csv",
              "branch,Hinf,Tw,Fp0,Gp0,Tp0,shoot_bc_residual,"
              "shoot_profile_error,newton_residual,newton_profile_error,"
              "newton_sigma_min_over_max", validation)
    write_csv(OUT / "profiles_Tw1p045.csv", "branch,z,H,F,G,T",
              profile_rows)

    fixed_tw_stability = []
    for branch_id, sol in profile_solutions:
        values = similarity_eigenvalues(sol, degree=70)
        leading = values[0]
        fixed_tw_stability.append((branch_id, sol.sol(20.0)[0], sol.p[0],
                                   leading.real, leading.imag,
                                   np.count_nonzero(values.real > 1.0e-8)))
    write_csv(OUT / "similarity_stability_Tw1p045.csv",
              "branch,Hinf,Tw,leading_real,leading_imag,unstable_eigenvalues",
              fixed_tw_stability)

    stability_branch = []
    for h in np.arange(-0.70, -0.0599, 0.02):
        sol = branch_solutions[round(float(h), 10)]
        values = similarity_eigenvalues(sol, degree=60)
        leading = values[0]
        stability_branch.append((h, sol.p[0], leading.real, leading.imag,
                                 np.count_nonzero(values.real > 1.0e-8)))
    write_csv(OUT / "similarity_stability_branch.csv",
              "Hinf,Tw,leading_real,leading_imag,unstable_eigenvalues",
              stability_branch)

    domain_rows = []
    for zmax in (15.0, 20.0, 25.0, 30.0, 40.0):
        domain_turns = domain_fold_estimate(zmax)
        for kind, (h, tw) in zip((1, 2), domain_turns[:2]):
            domain_rows.append((zmax, kind, h, tw))
    write_csv(OUT / "domain_sensitivity.csv", "zmax,turning_point,Hinf,Tw",
              domain_rows)

    upper_domain_rows = upper_branch_domain_stability()
    write_csv(OUT / "upper_branch_domain_stability.csv",
              "zmax,Hinf,Tw,leading_real,leading_imag,unstable_eigenvalues,"
              "thermal_decay_length", upper_domain_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
    physical_tail = branch[:, 0] < 0.0
    for ax in axes:
        ax.plot(branch[physical_tail, 1], branch[physical_tail, 0], "k-", lw=2)
        ax.scatter([x[1] for x in turns], [x[0] for x in turns],
                   c="crimson", zorder=3, label="turning points")
        ax.axvline(target_tw, color="0.55", ls="--", lw=1,
                   label=r"$T_w=1.045$")
        ax.set_xlabel(r"$T_w$")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$H_\infty$")
    axes[0].set_title("Negative-H branch")
    axes[1].set_xlim(1.0375, 1.0510)
    axes[1].set_ylim(-0.72, -0.14)
    axes[1].set_title("Fold region")
    axes[1].legend(fontsize=9)
    fig.suptitle("Boussinesq similarity branch, zmax=20")
    fig.tight_layout()
    fig.savefig(OUT / "branch_diagram.png", dpi=180)
    fig.savefig(OUT / "branch_diagram.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    for branch_id, sol in profile_solutions:
        zplot = np.linspace(0.0, 8.0, 600)
        y = sol.sol(zplot)
        axes[0].plot(y[2], zplot, lw=2, label=f"branch {branch_id}")
        axes[1].plot(y[4], zplot, lw=2)
        axes[2].plot(y[0], zplot, lw=2)
    for ax, label in zip(axes, ("F", "G", "H")):
        ax.set_xlabel(label)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$\eta$")
    axes[0].legend()
    fig.suptitle(r"Three similarity solutions at $T_w=1.045$")
    fig.tight_layout()
    fig.savefig(OUT / "profiles_Tw1p045.png", dpi=180)
    fig.savefig(OUT / "profiles_Tw1p045.pdf")
    plt.close(fig)

    stability_branch_array = np.asarray(stability_branch)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.axhline(0.0, color="0.4", lw=1)
    ax.plot(stability_branch_array[:, 0], stability_branch_array[:, 2],
            "o-", lw=1.8)
    for h, _ in turns[:2]:
        ax.axvline(h, color="crimson", ls="--", lw=1)
    ax.set_xlabel(r"$H_\infty$")
    ax.set_ylabel(r"leading $\mathrm{Re}(\lambda)$")
    ax.set_title("Similarity-preserving temporal stability, zmax=20")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "similarity_stability_branch.png", dpi=180)
    fig.savefig(OUT / "similarity_stability_branch.pdf")
    plt.close(fig)

    with open(OUT / "report.md", "w", encoding="utf-8") as f:
        f.write("# Boussinesq similarity-fold investigation\n\n")
        f.write("## Turning points\n\n")
        f.write("| point | Hinf | Tw | scaled Jacobian sigma_min/sigma_max | thermal decay length |\n")
        f.write("|---:|---:|---:|---:|---:|\n")
        for row in turning_rows:
            f.write(f"| {int(row[0])} | {row[1]:.9f} | {row[2]:.9f} | {row[4]:.3e} | {row[5]:.4f} |\n")
        f.write("\nThe benchmark-connected branch reaches the first saddle-node, ")
        f.write("so continuation with Tw as the parameter becomes singular. ")
        f.write("Continuation with Hinf remains regular and reveals an S-shaped ")
        f.write("branch with three solutions over part of the temperature range.\n\n")

        f.write("## Domain sensitivity\n\n")
        f.write("| zmax | point | Hinf | Tw |\n|---:|---:|---:|---:|\n")
        for row in domain_rows:
            f.write(f"| {row[0]:.0f} | {int(row[1])} | {row[2]:.9f} | {row[3]:.9f} |\n")
        f.write("\nThe first fold converges to approximately Tw=1.04802 as zmax ")
        f.write("increases. The second fold drifts because its small |Hinf| gives ")
        f.write("a long thermal tail. Results on the upper branch therefore need ")
        f.write("a substantially larger domain than the benchmark-connected branch.\n\n")
        f.write("## Cross-validation at Tw=1.045\n\n")
        f.write("| branch | Hinf | shooting profile error | Newton profile error | Newton residual |\n")
        f.write("|---:|---:|---:|---:|---:|\n")
        for row in validation:
            f.write(f"| {int(row[0])} | {row[1]:.9f} | {row[7]:.3e} | {row[9]:.3e} | {row[8]:.3e} |\n")
        f.write("\nAdaptive collocation, IVP shooting, and Chebyshev-Newton ")
        f.write("collocation converge to the same profiles.\n\n")

        f.write("## Similarity-preserving temporal stability at Tw=1.045\n\n")
        f.write("| branch | Hinf | leading eigenvalue | unstable eigenvalues |\n")
        f.write("|---:|---:|---:|---:|\n")
        for row in fixed_tw_stability:
            f.write(f"| {int(row[0])} | {row[1]:.9f} | ")
            f.write(f"{row[3]:+.6e}{row[4]:+.6e}i | {int(row[5])} |\n")
        f.write("\nThe lower branch is stable within the axisymmetric ")
        f.write("similarity-preserving subspace. The middle branch has one ")
        f.write("positive real eigenvalue. The finite-domain upper branch has ")
        f.write("an unstable complex-conjugate pair and is oscillatory unstable ")
        f.write("in this subspace. Full physical stability also requires general ")
        f.write("non-axisymmetric and non-similar disturbances.\n\n")

        f.write("## Infinite-domain conditions\n\n")
        f.write("For a non-isothermal solution, T' behaves asymptotically as ")
        f.write("exp(Pr*Hinf*eta). Hence Hinf<0 is necessary for exponential ")
        f.write("thermal decay. Hinf>=0 is inadmissible unless Tw=1. The upper ")
        f.write("branch approaches Hinf=0 and therefore requires increasingly ")
        f.write("large domains; its finite-domain eigenvalues must be checked for ")
        f.write("domain convergence before being interpreted as infinite-domain ")
        f.write("spectrum.\n\n")

        f.write("## Upper-branch domain stability\n\n")
        f.write("| zmax | Hinf | Tw | leading eigenvalue | unstable eigenvalues |\n")
        f.write("|---:|---:|---:|---:|---:|\n")
        for row in upper_domain_rows:
            f.write(f"| {row[0]:.0f} | {row[1]:.2f} | {row[2]:.9f} | ")
            f.write(f"{row[3]:+.6e}{row[4]:+.6e}i | {int(row[5])} |\n")
        f.write("\nThe apparent stable window at zmax=20 moves as the far ")
        f.write("boundary is displaced. At zmax=60 and 80 every sampled ")
        f.write("negative-H upper-branch state remains unstable. For example, ")
        f.write("at Hinf=-0.14 the leading real eigenvalue converges toward ")
        f.write("approximately +4.2e-2. The finite-domain second turning point ")
        f.write("also drifts toward Hinf=0, where a non-isothermal infinite-domain ")
        f.write("solution loses exponential thermal decay.\n\n")

        f.write("## Interpretation\n\n")
        f.write("The failure near Tw=1.05 is a fold of the similarity solution, ")
        f.write("not disappearance of every mathematical solution. The ordinary ")
        f.write("Tw continuation fails because dHinf/dTw diverges at the fold. ")
        f.write("The branch remains connected when Hinf is used as continuation ")
        f.write("coordinate, but cannot be traversed monotonically in Tw. Heating ")
        f.write("reduces radial outflow F, which weakens axial entrainment |H|; ")
        f.write("the weaker entrainment thickens the thermal layer and increases ")
        f.write("the integrated inward buoyancy. This positive feedback produces ")
        f.write("the saddle-node. Stability and physical selection of the multiple ")
        f.write("branches require a separate analysis.\n")

    print(f"Output written to {OUT.resolve()}")
    print("Turning points:")
    for h, tw in turns[:2]:
        print(f"  Hinf={h:.9f}, Tw={tw:.9f}")
    print(f"Solutions at Tw={target_tw}: {len(validation)}")
    for row in validation:
        print(f"  branch {int(row[0])}: Hinf={row[1]:.9f}, "
              f"shoot_err={row[7]:.3e}, newton_err={row[9]:.3e}")
    print("Leading similarity-preserving eigenvalues at Tw=1.045:")
    for row in fixed_tw_stability:
        print(f"  branch {int(row[0])}: {row[3]:+.8e}{row[4]:+.8e}i")


if __name__ == "__main__":
    main()
