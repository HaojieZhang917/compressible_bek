#!/usr/bin/env python3
"""Dynamical-systems analysis of the three Re_h=1000 similarity states.

Two distinct notions are kept separate.

1. Spatial dynamics: the steady seven-state BVP is an autonomous ODE in z.
   Its frozen Jacobian and the phase-space trajectory diagnose saddle
   directions and spatially damped/oscillatory structure.
2. Temporal dynamics inside the axisymmetric similarity subspace: the
   unsteady similarity equations are linearised and discretised spectrally.
   The radial pressure gradient is a Lagrange multiplier that preserves the
   zero-net-radial-flow constraint implied by H(0)=H(1)=0.

The temporal calculation is not the full three-dimensional hydrodynamic
stability problem: non-axisymmetric and radially varying perturbations are
excluded.  It nevertheless determines the stability of each steady branch
to similarity-preserving perturbations and distinguishes monotonic from
oscillatory branch instabilities.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Chebyshev, chebyshev
from scipy import linalg
from scipy.interpolate import CubicSpline

from two_disk_boussinesq_singularity import (
    SolverConfig,
    TwoDiskBoussinesqSolver,
)


SCRIPT_DIR = Path(__file__).resolve().parent
ROTOR_STATOR_ROOT = SCRIPT_DIR.parent
DEFAULT_INPUT = (
    ROTOR_STATOR_ROOT
    / "data"
    / "boussinesq_singularity_results"
    / "three_solutions_Tw1.160.csv"
)
DEFAULT_OUTPUT = ROTOR_STATOR_ROOT / "data" / "three_solution_dynamics"


BRANCH_LABELS = {
    1: "upper / low-Pi",
    2: "middle",
    3: "principal / isothermal-connected",
}


@dataclass(frozen=True)
class BaseProfile:
    branch: int
    pressure: float
    z: np.ndarray
    h: np.ndarray
    f: np.ndarray
    g: np.ndarray
    temp: np.ndarray


@dataclass
class TemporalResult:
    n: int
    z: np.ndarray
    eigvals: np.ndarray
    eigvecs: np.ndarray
    f_map: np.ndarray
    h_map: np.ndarray
    g_map: np.ndarray
    temp_map: np.ndarray
    n_f_reduced: int
    n_g: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
    )
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--re-h", type=float, default=1000.0)
    parser.add_argument("--pr", type=float, default=0.72)
    parser.add_argument(
        "--collocation-orders",
        type=int,
        nargs="+",
        default=[61, 81, 101, 121],
    )
    parser.add_argument("--branch-scan-order", type=int, default=81)
    parser.add_argument(
        "--skip-branch-scan",
        action="store_true",
        help="Skip the leading temporal eigenvalue scan along the S-curve.",
    )
    return parser.parse_args()


def load_profiles(path: Path) -> list[BaseProfile]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    profiles: list[BaseProfile] = []
    for branch_value in np.unique(data["branch"]):
        branch = int(round(float(branch_value)))
        mask = np.isclose(data["branch"], branch_value)
        order = np.argsort(data["z"][mask])
        profiles.append(
            BaseProfile(
                branch=branch,
                pressure=float(data["pressure_gradient"][mask][0]),
                z=np.asarray(data["z"][mask][order]),
                h=np.asarray(data["H"][mask][order]),
                f=np.asarray(data["F"][mask][order]),
                g=np.asarray(data["G"][mask][order]),
                temp=np.asarray(data["T"][mask][order]),
            )
        )
    return sorted(profiles, key=lambda item: item.branch)


def chebyshev_operators(n: int) -> tuple[np.ndarray, ...]:
    """Return ascending nodes, D1, D2, integration matrix, and weights."""
    x = -np.cos(np.pi * np.arange(n) / (n - 1))
    z = 0.5 * (x + 1.0)
    vandermonde = chebyshev.chebvander(x, n - 1)
    inverse = linalg.inv(vandermonde)

    derivative_1 = np.empty_like(vandermonde)
    derivative_2 = np.empty_like(vandermonde)
    integral = np.empty_like(vandermonde)
    for degree in range(n):
        basis = Chebyshev.basis(degree)
        derivative_1[:, degree] = 2.0 * basis.deriv(1)(x)
        derivative_2[:, degree] = 4.0 * basis.deriv(2)(x)
        primitive = basis.integ(1)
        integral[:, degree] = 0.5 * (primitive(x) - primitive(-1.0))

    d1 = derivative_1 @ inverse
    d2 = derivative_2 @ inverse
    integration = integral @ inverse
    weights = integration[-1].copy()
    return z, d1, d2, integration, weights


def interpolate_profile(profile: BaseProfile, z: np.ndarray) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for name in ("h", "f", "g", "temp"):
        spline = CubicSpline(profile.z, getattr(profile, name))
        result[name] = spline(z)
        result[f"{name}_z"] = spline(z, 1)
        result[f"{name}_zz"] = spline(z, 2)
    return result


def temporal_operator(
    profile: BaseProfile,
    n: int,
    re_h: float,
    pr: float,
) -> tuple[np.ndarray, TemporalResult]:
    z, d1, d2, integration, weights = chebyshev_operators(n)
    base = interpolate_profile(profile, z)
    interior = np.arange(1, n - 1)
    m = interior.size
    sqrt_re = np.sqrt(re_h)

    extension = np.zeros((n, m))
    extension[interior, np.arange(m)] = 1.0
    weights_i = weights[interior]
    q = linalg.null_space(weights_i[None, :])
    if q.shape != (m, m - 1):
        raise RuntimeError("Unexpected dimension of the zero-integral basis")

    pressure_projection = np.eye(m) - np.outer(
        np.ones(m), weights_i / np.sum(weights_i)
    )
    f_map = extension @ q
    g_map = extension.copy()
    temp_map = extension.copy()
    h_map = -2.0 * sqrt_re * integration @ f_map

    f_i = base["f"][interior]
    g_i = base["g"][interior]
    h_i = base["h"][interior]
    f_z_i = base["f_z"][interior]
    g_z_i = base["g_z"][interior]
    temp_z_i = base["temp_z"][interior]

    d1_f = (d1 @ f_map)[interior]
    d2_f = (d2 @ f_map)[interior]
    d1_g = (d1 @ g_map)[interior]
    d2_g = (d2 @ g_map)[interior]
    d1_temp = (d1 @ temp_map)[interior]
    d2_temp = (d2 @ temp_map)[interior]
    f_i_map = f_map[interior]
    h_i_map = h_map[interior]

    radial_f = (
        d2_f / re_h
        - h_i[:, None] * d1_f / sqrt_re
        - f_z_i[:, None] * h_i_map / sqrt_re
        - 2.0 * f_i[:, None] * f_i_map
    )
    radial_g = 2.0 * np.diag(g_i)
    radial_temp = -np.eye(m)

    azimuthal_f = (
        -g_z_i[:, None] * h_i_map / sqrt_re
        - 2.0 * g_i[:, None] * f_i_map
    )
    azimuthal_g = (
        d2_g / re_h
        - h_i[:, None] * d1_g / sqrt_re
        - 2.0 * np.diag(f_i)
    )

    thermal_f = -temp_z_i[:, None] * h_i_map / sqrt_re
    thermal_temp = (
        d2_temp / (re_h * pr)
        - h_i[:, None] * d1_temp / sqrt_re
    )

    n_a = m - 1
    zero_a_temp = np.zeros((n_a, m))
    zero_g_temp = np.zeros((m, m))
    operator = np.block(
        [
            [
                q.T @ pressure_projection @ radial_f,
                q.T @ pressure_projection @ radial_g,
                q.T @ pressure_projection @ radial_temp,
            ],
            [azimuthal_f, azimuthal_g, zero_g_temp],
            [thermal_f, np.zeros((m, m)), thermal_temp],
        ]
    )

    eigvals, eigvecs = linalg.eig(operator)
    order = np.lexsort((np.abs(eigvals.imag), -eigvals.real))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    return operator, TemporalResult(
        n=n,
        z=z,
        eigvals=eigvals,
        eigvecs=eigvecs,
        f_map=f_map,
        h_map=h_map,
        g_map=g_map,
        temp_map=temp_map,
        n_f_reduced=n_a,
        n_g=m,
    )


def spatial_jacobian(
    state: dict[str, float], re_h: float, pr: float
) -> np.ndarray:
    sqrt_re = np.sqrt(re_h)
    h = state["h"]
    f = state["f"]
    fp = state["f_z"]
    g = state["g"]
    gp = state["g_z"]
    tp = state["temp_z"]
    matrix = np.zeros((7, 7))
    # State ordering: H, F, F_z, G, G_z, T, T_z.
    matrix[0, 1] = -2.0 * sqrt_re
    matrix[1, 2] = 1.0
    matrix[2, 0] = sqrt_re * fp
    matrix[2, 1] = 2.0 * re_h * f
    matrix[2, 2] = sqrt_re * h
    matrix[2, 3] = -2.0 * re_h * g
    matrix[2, 5] = re_h
    matrix[3, 4] = 1.0
    matrix[4, 0] = sqrt_re * gp
    matrix[4, 1] = 2.0 * re_h * g
    matrix[4, 3] = 2.0 * re_h * f
    matrix[4, 4] = sqrt_re * h
    matrix[5, 6] = 1.0
    matrix[6, 0] = pr * sqrt_re * tp
    matrix[6, 6] = pr * sqrt_re * h
    return matrix


def count_sign_changes(values: np.ndarray, tolerance: float = 2.0e-4) -> int:
    filtered = values[np.abs(values) > tolerance]
    if filtered.size < 2:
        return 0
    return int(np.count_nonzero(np.diff(np.sign(filtered))))


def count_extrema(values: np.ndarray, z: np.ndarray, tolerance: float = 2.0e-3) -> int:
    derivative = CubicSpline(z, values)(z, 1)
    filtered = derivative[np.abs(derivative) > tolerance]
    if filtered.size < 2:
        return 0
    return int(np.count_nonzero(np.diff(np.sign(filtered))))


def branch_diagnostics(
    profile: BaseProfile, re_h: float, pr: float
) -> tuple[dict[str, float | int | str], np.ndarray]:
    z = np.linspace(0.0, 1.0, 1001)
    base = interpolate_profile(profile, z)
    core = (z >= 0.30) & (z <= 0.70)
    inviscid_residual = base["g"] ** 2 - (base["temp"] - 1.0) - profile.pressure
    radial_residual = (
        base["f_zz"] / re_h
        - base["h"] * base["f_z"] / np.sqrt(re_h)
        + base["g"] ** 2
        - base["f"] ** 2
        - (base["temp"] - 1.0)
        - profile.pressure
    )
    azimuthal_residual = (
        base["g_zz"] / re_h
        - base["h"] * base["g_z"] / np.sqrt(re_h)
        - 2.0 * base["f"] * base["g"]
    )
    thermal_residual = (
        base["temp_zz"] / (re_h * pr)
        - base["h"] * base["temp_z"] / np.sqrt(re_h)
    )
    mid_index = int(np.argmin(np.abs(z - 0.5)))
    state = {name: float(values[mid_index]) for name, values in base.items()}
    spatial_eigs = linalg.eigvals(spatial_jacobian(state, re_h, pr))
    complex_spatial = spatial_eigs[np.abs(spatial_eigs.imag) > 1.0e-8]
    min_wavelength = (
        float(np.min(2.0 * np.pi / np.abs(complex_spatial.imag)))
        if complex_spatial.size
        else float("nan")
    )
    row: dict[str, float | int | str] = {
        "branch": profile.branch,
        "label": BRANCH_LABELS[profile.branch],
        "pressure_gradient": profile.pressure,
        "H_mid": state["h"],
        "F_mid": state["f"],
        "G_mid": state["g"],
        "T_mid": state["temp"],
        "core_max_abs_F": float(np.max(np.abs(base["f"][core]))),
        "core_G_range": float(np.ptp(base["g"][core])),
        "core_rms_Gz": float(np.sqrt(np.mean(base["g_z"][core] ** 2))),
        "core_rms_inviscid_residual": float(
            np.sqrt(np.mean(inviscid_residual[core] ** 2))
        ),
        "mass_integral_F": float(np.trapezoid(base["f"], z)),
        "steady_radial_residual_max": float(
            np.max(np.abs(radial_residual[1:-1]))
        ),
        "steady_azimuthal_residual_max": float(
            np.max(np.abs(azimuthal_residual[1:-1]))
        ),
        "steady_thermal_residual_max": float(
            np.max(np.abs(thermal_residual[1:-1]))
        ),
        "F_interior_sign_changes": count_sign_changes(base["f"][1:-1]),
        "G_interior_extrema": count_extrema(base["g"][1:-1], z[1:-1]),
        "mid_spatial_positive_real_count": int(
            np.count_nonzero(spatial_eigs.real > 1.0e-8)
        ),
        "mid_spatial_negative_real_count": int(
            np.count_nonzero(spatial_eigs.real < -1.0e-8)
        ),
        "mid_spatial_neutral_count": int(
            np.count_nonzero(np.abs(spatial_eigs.real) <= 1.0e-8)
        ),
        "shortest_mid_spatial_wavelength": min_wavelength,
    }
    return row, spatial_eigs


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_temporal_spectra(
    final_results: dict[int, TemporalResult], output: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7), sharex=True, sharey=True)
    colors = {1: "#2563eb", 2: "#d97706", 3: "#7e22ce"}
    for branch, axis in zip(sorted(final_results), axes):
        eigvals = final_results[branch].eigvals
        visible = (eigvals.real > -2.5) & (np.abs(eigvals.imag) < 8.0)
        axis.scatter(
            eigvals.real[visible],
            eigvals.imag[visible],
            s=20,
            color=colors[branch],
            alpha=0.75,
        )
        axis.axvline(0.0, color="black", linewidth=1.0)
        axis.axhline(0.0, color="0.75", linewidth=0.8)
        leading = eigvals[0]
        axis.scatter(
            [leading.real], [leading.imag], s=80, marker="*", color="#dc2626", zorder=5
        )
        axis.set_title(f"Branch {branch}: {BRANCH_LABELS[branch]}\n"
                       rf"$\lambda_{{max}}={leading.real:.4f}{leading.imag:+.4f}i$")
        axis.grid(alpha=0.25)
        axis.set_xlabel(r"$\mathrm{Re}(\lambda)/\Omega$")
    axes[0].set_ylabel(r"$\mathrm{Im}(\lambda)/\Omega$")
    fig.suptitle(
        "Temporal spectra in the axisymmetric similarity subspace",
        fontsize=15,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_dominant_modes(
    final_results: dict[int, TemporalResult], output: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7), sharey=True)
    for branch, axis in zip(sorted(final_results), axes):
        result = final_results[branch]
        vector = result.eigvecs[:, 0]
        n_a = result.n_f_reduced
        n_g = result.n_g
        a = vector[:n_a]
        b = vector[n_a : n_a + n_g]
        c = vector[n_a + n_g :]
        f_mode = result.f_map @ a
        h_mode = result.h_map @ a
        g_mode = result.g_map @ b
        temp_mode = result.temp_map @ c
        scale = max(
            np.max(np.abs(f_mode)),
            np.max(np.abs(h_mode)),
            np.max(np.abs(g_mode)),
            np.max(np.abs(temp_mode)),
        )
        axis.plot(result.z, np.abs(f_mode) / scale, label=r"$|f|$")
        axis.plot(result.z, np.abs(g_mode) / scale, label=r"$|g|$")
        axis.plot(result.z, np.abs(h_mode) / scale, label=r"$|h|$")
        axis.plot(result.z, np.abs(temp_mode) / scale, label=r"$|\vartheta|$")
        axis.set_title(f"Branch {branch}: {BRANCH_LABELS[branch]}")
        axis.set_xlabel("z/h")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("normalised modal amplitude")
    axes[-1].legend(loc="upper right", fontsize=9)
    fig.suptitle("Dominant similarity-preserving temporal mode", fontsize=15)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_spatial_dynamics(profiles: list[BaseProfile], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.7))
    colors = {1: "#2563eb", 2: "#d97706", 3: "#7e22ce"}
    for profile in profiles:
        z = np.linspace(0.0, 1.0, 1001)
        base = interpolate_profile(profile, z)
        label = f"B{profile.branch}: Pi={profile.pressure:.5f}"
        axes[0].plot(base["f"], base["f_z"], color=colors[profile.branch], label=label)
        axes[1].plot(base["g"], base["g_z"], color=colors[profile.branch])
        balance = base["g"] ** 2 - (base["temp"] - 1.0) - profile.pressure
        axes[2].plot(z, balance, color=colors[profile.branch])
    axes[0].set_xlabel("F")
    axes[0].set_ylabel(r"$F_z$")
    axes[0].set_title("Radial phase projection")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("G")
    axes[1].set_ylabel(r"$G_z$")
    axes[1].set_title("Azimuthal phase projection")
    axes[2].axhline(0.0, color="black", linewidth=0.9)
    axes[2].set_xlabel("z/h")
    axes[2].set_ylabel(r"$G^2-(T-1)-\Pi$")
    axes[2].set_title("Departure from inviscid-core balance")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle("Spatial-dynamical structure of the three steady trajectories", fontsize=15)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def scan_temporal_growth_along_branch(
    script_dir: Path,
    re_h: float,
    pr: float,
    collocation_order: int,
) -> list[dict]:
    solver = TwoDiskBoussinesqSolver(
        SolverConfig(
            re_h=re_h,
            pr=pr,
            tol=2.0e-7,
            initial_nodes=401,
            model="traditional_centrifugal",
        ),
        script_dir / "baseflow_Res1000.npz",
    )
    isothermal = solver.solve_isothermal()
    fold_path = script_dir / "boussinesq_singularity_results" / (
        f"folds_Re{re_h:g}_traditional_centrifugal.json"
    )
    folds = json.loads(fold_path.read_text(encoding="utf-8"))
    fold_pressures = np.array([item["pressure_gradient"] for item in folds])
    initial_pressure = float(isothermal.p[0])
    coarse_high = np.arange(initial_pressure, 0.015, -0.002)
    moderate = np.arange(coarse_high[-1] - 0.001, 0.011, -0.001)
    fine = np.arange(moderate[-1] - 0.0005, -0.0011, -0.0005)
    pressure_values = np.sort(
        np.unique(np.concatenate((coarse_high, moderate, fine, fold_pressures)))
    )[::-1]
    branch = solver.continue_pressure(
        pressure_values,
        isothermal,
        profile_nodes=1001,
    )

    rows: list[dict] = []
    for index, pressure in enumerate(branch.columns["pressure_gradient"]):
        values = branch.profiles[index]
        profile = BaseProfile(
            branch=0,
            pressure=float(pressure),
            z=branch.z_profile,
            h=values[0],
            f=values[1],
            g=values[3],
            temp=values[5],
        )
        _operator, temporal = temporal_operator(
            profile, collocation_order, re_h, pr
        )
        leading = temporal.eigvals[0]
        rows.append(
            {
                "pressure_gradient": float(pressure),
                "Tw": float(branch.columns["Tw"][index]),
                "leading_real": float(leading.real),
                "leading_imag": float(leading.imag),
                "positive_real_count": int(
                    np.count_nonzero(temporal.eigvals.real > 1.0e-7)
                ),
                "fold_sample": int(
                    np.any(np.isclose(pressure, fold_pressures, rtol=0.0, atol=1.0e-12))
                ),
            }
        )
    return rows


def plot_branch_growth(rows: list[dict], output: Path) -> None:
    pressure = np.asarray([row["pressure_gradient"] for row in rows])
    tw = np.asarray([row["Tw"] for row in rows])
    growth = np.asarray([row["leading_real"] for row in rows])
    unstable = growth > 0.0
    folds = np.asarray([bool(row["fold_sample"]) for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    axes[0].plot(tw, pressure, color="0.65", linewidth=1.2)
    axes[0].scatter(tw[~unstable], pressure[~unstable], color="#15803d", s=30,
                    label="similarity-stable")
    axes[0].scatter(tw[unstable], pressure[unstable], color="#dc2626", s=30,
                    label="one real unstable mode")
    axes[0].scatter(tw[folds], pressure[folds], color="#111827", marker="*", s=130,
                    label="fold samples", zorder=5)
    axes[0].set_xlabel(r"$T_w$")
    axes[0].set_ylabel(r"$\Pi$")
    axes[0].set_title("Stability exchange on the S-curve")
    axes[0].legend(fontsize=8)

    axes[1].plot(pressure, growth, color="#1d4ed8", marker="o", markersize=3)
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[1].scatter(pressure[folds], growth[folds], color="#111827", marker="*", s=130,
                    zorder=5)
    axes[1].set_xlabel(r"$\Pi$")
    axes[1].set_ylabel(r"$\max\,\mathrm{Re}(\lambda)/\Omega$")
    axes[1].set_title("Real eigenvalue crosses zero at each fold")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles = load_profiles(input_path)

    diagnostic_rows: list[dict] = []
    spatial_rows: list[dict] = []
    for profile in profiles:
        diagnostic, spatial_eigs = branch_diagnostics(profile, args.re_h, args.pr)
        diagnostic_rows.append(diagnostic)
        for rank, eigval in enumerate(
            sorted(spatial_eigs, key=lambda value: (-value.real, abs(value.imag))), start=1
        ):
            spatial_rows.append(
                {
                    "branch": profile.branch,
                    "label": BRANCH_LABELS[profile.branch],
                    "rank": rank,
                    "real": float(eigval.real),
                    "imag": float(eigval.imag),
                    "spatial_wavelength": (
                        float(2.0 * np.pi / abs(eigval.imag))
                        if abs(eigval.imag) > 1.0e-10
                        else float("nan")
                    ),
                }
            )

    temporal_rows: list[dict] = []
    final_results: dict[int, TemporalResult] = {}
    convergence: dict[str, dict] = {}
    for profile in profiles:
        convergence[str(profile.branch)] = {}
        for n in args.collocation_orders:
            _operator, result = temporal_operator(profile, n, args.re_h, args.pr)
            leading = result.eigvals[0]
            convergence[str(profile.branch)][str(n)] = {
                "leading_real": float(leading.real),
                "leading_imag": float(leading.imag),
                "positive_real_count": int(np.count_nonzero(result.eigvals.real > 1.0e-7)),
            }
            for rank, eigval in enumerate(result.eigvals[:20], start=1):
                temporal_rows.append(
                    {
                        "branch": profile.branch,
                        "label": BRANCH_LABELS[profile.branch],
                        "collocation_points": n,
                        "rank": rank,
                        "real": float(eigval.real),
                        "imag": float(eigval.imag),
                        "oscillation_frequency_over_Omega": float(abs(eigval.imag)),
                        "period_Omega_t": (
                            float(2.0 * np.pi / abs(eigval.imag))
                            if abs(eigval.imag) > 1.0e-10
                            else float("nan")
                        ),
                    }
                )
            if n == args.collocation_orders[-1]:
                final_results[profile.branch] = result

    for row in diagnostic_rows:
        result = final_results[int(row["branch"])]
        leading = result.eigvals[0]
        vector = result.eigvecs[:, 0]
        n_a = result.n_f_reduced
        n_g = result.n_g
        a = vector[:n_a]
        b = vector[n_a : n_a + n_g]
        c = vector[n_a + n_g :]
        mode_components = {
            "F": result.f_map @ a,
            "H": result.h_map @ a,
            "G": result.g_map @ b,
            "T": result.temp_map @ c,
        }
        component_maxima = {
            name: float(np.max(np.abs(values)))
            for name, values in mode_components.items()
        }
        mode_scale = max(component_maxima.values())
        total_amplitude_sq = sum(
            np.abs(values / mode_scale) ** 2 for values in mode_components.values()
        )
        mode_centroid = float(
            np.trapezoid(result.z * total_amplitude_sq, result.z)
            / np.trapezoid(total_amplitude_sq, result.z)
        )
        oscillatory = result.eigvals[np.abs(result.eigvals.imag) > 1.0e-7]
        least_damped_oscillatory = oscillatory[
            np.argmax(oscillatory.real)
        ]
        row["leading_temporal_real"] = float(leading.real)
        row["leading_temporal_imag"] = float(leading.imag)
        row["temporal_positive_real_count"] = int(
            np.count_nonzero(result.eigvals.real > 1.0e-7)
        )
        for name, maximum in component_maxima.items():
            row[f"dominant_mode_{name}_relative_max"] = maximum / mode_scale
        row["dominant_mode_z_centroid"] = mode_centroid
        row["least_damped_oscillatory_real"] = float(
            least_damped_oscillatory.real
        )
        row["least_damped_oscillatory_imag_abs"] = float(
            abs(least_damped_oscillatory.imag)
        )
        row["least_damped_oscillatory_period"] = float(
            2.0 * np.pi / abs(least_damped_oscillatory.imag)
        )
        row["least_damped_oscillatory_quality_factor"] = float(
            abs(least_damped_oscillatory.imag)
            / (2.0 * abs(least_damped_oscillatory.real))
        )
        row["similarity_temporal_class"] = (
            "unstable-oscillatory"
            if leading.real > 1.0e-7 and abs(leading.imag) > 1.0e-7
            else "unstable-monotonic"
            if leading.real > 1.0e-7
            else "stable-with-damped-oscillations"
            if np.any(np.abs(result.eigvals.imag) > 1.0e-7)
            else "stable-monotonic"
        )

    write_csv(output_dir / "branch_diagnostics.csv", diagnostic_rows)
    write_csv(output_dir / "midplane_spatial_eigenvalues.csv", spatial_rows)
    write_csv(output_dir / "temporal_eigenvalues_convergence.csv", temporal_rows)
    with (output_dir / "temporal_convergence.json").open("w", encoding="utf-8") as stream:
        json.dump(convergence, stream, indent=2, ensure_ascii=False)

    plot_temporal_spectra(final_results, output_dir / "temporal_spectra.png")
    plot_dominant_modes(final_results, output_dir / "dominant_temporal_modes.png")
    plot_spatial_dynamics(profiles, output_dir / "spatial_phase_portraits.png")

    if not args.skip_branch_scan:
        branch_growth_rows = scan_temporal_growth_along_branch(
            Path(__file__).resolve().parent,
            args.re_h,
            args.pr,
            args.branch_scan_order,
        )
        write_csv(output_dir / "branch_temporal_growth.csv", branch_growth_rows)
        plot_branch_growth(
            branch_growth_rows,
            output_dir / "branch_stability_exchange.png",
        )

    print(json.dumps(diagnostic_rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
