#!/usr/bin/env python3
"""Validate and compare the von Karman and rotor--stator saddle nodes.

All original result files are read-only.  Derived spectra, normal-form checks,
mode shapes, tables, and figures are written to
``dynamical_singularity_comparison``.
"""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg


SCRIPT_DIR = Path(__file__).resolve().parent
TOOLKIT = SCRIPT_DIR.parent
TWODISK = SCRIPT_DIR
VONKARMAN = TOOLKIT / "compress" / "BEK" / "Vonkarmen_bone"
OUT = TWODISK / "dynamical_singularity_comparison"
sys.path.insert(0, str(VONKARMAN / "scripts"))
sys.path.insert(0, str(TWODISK))

import investigate_boussinesq_fold as vk  # noqa: E402
from analyze_three_solution_dynamics import temporal_operator  # noqa: E402
from two_disk_boussinesq_singularity import (  # noqa: E402
    SolverConfig,
    TwoDiskBoussinesqSolver,
)


@dataclass(frozen=True)
class VKProfile:
    branch: int
    pressure: float
    z: np.ndarray
    h: np.ndarray
    f: np.ndarray
    g: np.ndarray
    temp: np.ndarray


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv_dicts(path: Path) -> list[dict]:
    def parse(value: str):
        try:
            return float(value)
        except ValueError:
            return value

    with path.open(newline="", encoding="utf-8") as stream:
        return [
            {key: parse(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


def first_fold_by_domain() -> dict[float, tuple[float, float]]:
    rows = read_csv_dicts(
        VONKARMAN / "boussinesq_domain_branches" / "turning_points_by_zmax.csv"
    )
    return {
        row["zmax"]: (row["Hinf"], row["Tw"])
        for row in rows
        if int(round(row["turning_point"])) == 1
    }


def continue_vk(
    zmax: float, target_h_values: np.ndarray, tol: float = 2.0e-8
) -> dict[float, object]:
    target_h_values = np.asarray(sorted(set(float(x) for x in target_h_values)))
    h_max = float(np.max(target_h_values))
    n = max(401, int(round(12.0 * zmax)) + 1)
    z = np.linspace(0.0, zmax, n)
    seed = vk.solve_isothermal(zmax, n, tol=2.0e-9)
    guess = seed.sol(z)
    tw = 1.0
    stepping = np.arange(-0.75, h_max + 0.0101, 0.01)
    h_values = np.sort(np.unique(np.concatenate((stepping, target_h_values))))
    h_values = h_values[(h_values >= -0.7500001) & (h_values <= h_max + 1.0e-12)]
    solutions: dict[float, object] = {}
    for h_inf in h_values:
        sol = vk.solve_fixed_h(float(h_inf), z, guess, tw, tol=tol)
        guess = sol.sol(z)
        tw = float(sol.p[0])
        if np.any(np.isclose(h_inf, target_h_values, rtol=0.0, atol=2.0e-12)):
            solutions[round(float(h_inf), 12)] = sol
    return solutions


def vk_similarity_spectrum(source_sol, zmax: float, degree: int):
    """Generalised temporal spectrum plus eigenvectors in [f,g,h,theta]."""
    z, d1 = vk.cheb_matrix(degree, zmax)
    d2 = d1 @ d1
    n = z.size
    identity = np.eye(n)
    y = source_sol.sol(z)
    h_base, f_base, g_base, temp_base = y[0], y[2], y[4], y[6]
    fp = d1 @ f_base
    gp = d1 @ g_base
    tp = d1 @ temp_base
    sf, sg, sh, st = [slice(i * n, (i + 1) * n) for i in range(4)]
    matrix_a = np.zeros((4 * n, 4 * n))
    matrix_b = np.zeros_like(matrix_a)

    matrix_a[sf, sf] = d2 - np.diag(h_base) @ d1 - np.diag(2.0 * f_base)
    matrix_a[sf, sg] = np.diag(2.0 * (g_base - 1.0))
    matrix_a[sf, sh] = -np.diag(fp)
    matrix_a[sf, st] = -identity
    matrix_a[sg, sf] = np.diag(2.0 * (1.0 - g_base))
    matrix_a[sg, sg] = d2 - np.diag(h_base) @ d1 - np.diag(2.0 * f_base)
    matrix_a[sg, sh] = -np.diag(gp)
    matrix_a[sh, sf] = 2.0 * identity
    matrix_a[sh, sh] = d1
    matrix_a[st, sh] = -np.diag(tp)
    matrix_a[st, st] = d2 / vk.PR - np.diag(h_base) @ d1
    matrix_b[sf, sf] = identity
    matrix_b[sg, sg] = identity
    matrix_b[st, st] = identity

    def set_bc(row: int, column: int) -> None:
        matrix_a[row, :] = 0.0
        matrix_b[row, :] = 0.0
        matrix_a[row, column] = 1.0

    set_bc(0, 0)
    set_bc(n - 1, n - 1)
    set_bc(n, n)
    set_bc(2 * n - 1, 2 * n - 1)
    set_bc(2 * n, 2 * n)
    set_bc(3 * n, 3 * n)
    set_bc(4 * n - 1, 4 * n - 1)

    values, vectors = linalg.eig(matrix_a, matrix_b)
    keep = np.isfinite(values) & (np.abs(values) < 1.0e5)
    values = values[keep]
    vectors = vectors[:, keep]
    order = np.lexsort((np.abs(values.imag), -values.real))
    return z, values[order], vectors[:, order]


def leading_vk(source_sol, zmax: float, degree: int) -> complex:
    return complex(vk_similarity_spectrum(source_sol, zmax, degree)[1][0])


def through_origin_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    slope = float(np.dot(x, y) / np.dot(x, x))
    residual = y - slope * x
    r2_origin = float(1.0 - np.dot(residual, residual) / np.dot(y, y))
    return slope, r2_origin


def vk_fold_scan(zmax: float, h_fold: float, tw_fold: float):
    offsets = np.array(
        [
            -0.12, -0.10, -0.08, -0.06, -0.04, -0.03, -0.02,
            -0.015, -0.010, -0.0075, -0.005, -0.0025, 0.0,
            0.0025, 0.005, 0.0075, 0.010, 0.015, 0.020, 0.030,
            0.040, 0.060, 0.080,
        ]
    )
    targets = h_fold + offsets
    solutions = continue_vk(zmax, targets, tol=8.0e-9)
    rows: list[dict] = []
    for h_inf in targets:
        sol = solutions[round(float(h_inf), 12)]
        leading = leading_vk(sol, zmax, degree=130)
        rows.append(
            {
                "Hinf": float(h_inf),
                "Tw": float(sol.p[0]),
                "leading_real": float(leading.real),
                "leading_imag": float(leading.imag),
                "unstable_count_class": 0 if leading.real < 0.0 else 1,
                "fold_sample": int(abs(h_inf - h_fold) < 1.0e-12),
            }
        )

    fit_rows: list[dict] = []
    for side_name, side_mask in (
        ("benchmark_connected_stable", np.asarray([row["Hinf"] < h_fold for row in rows])),
        ("post_fold_unstable", np.asarray([row["Hinf"] > h_fold for row in rows])),
    ):
        mu = np.asarray([tw_fold - row["Tw"] for row in rows])
        lam2 = np.asarray([row["leading_real"] ** 2 for row in rows])
        # The square-root law is a local asymptotic result.  Restrict the fit
        # to the closest resolved points so the two sides approach the same
        # normal-form coefficient rather than fitting higher-order curvature.
        fit_cap = 2.0e-5
        use = side_mask & (mu > 1.0e-9) & (mu < fit_cap)
        slope, r2 = through_origin_fit(mu[use], lam2[use])
        fit_rows.append(
            {
                "model": "von_Karman",
                "fold": "principal",
                "side": side_name,
                "points": int(np.count_nonzero(use)),
                "temperature_distance_max": fit_cap,
                "lambda_squared_slope": slope,
                "r2_through_origin": r2,
            }
        )
    return rows, fit_rows, solutions[round(float(h_fold), 12)]


def steady_jacobian_checks(source_sol, zmax: float, tw_fold: float, degree: int):
    z, fields, history, _ratio = vk.newton_collocation(
        tw_fold, source_sol, zmax=zmax, degree=degree
    )
    d1 = vk.cheb_matrix(degree, zmax)[1]
    d2 = d1 @ d1
    u = fields.reshape(-1)
    residual, jacobian = vk.newton_residual_jacobian(u, tw_fold, d1, d2)
    row_norm = np.linalg.norm(jacobian, axis=1)
    scale = 1.0 / np.maximum(row_norm, 1.0e-300)
    scaled_jacobian = scale[:, None] * jacobian
    left, singular, right_h = linalg.svd(scaled_jacobian)
    psi = left[:, -1]
    phi = right_h[-1]

    delta_tw = 1.0e-6
    r_plus = vk.newton_residual_jacobian(u, tw_fold + delta_tw, d1, d2)[0]
    r_minus = vk.newton_residual_jacobian(u, tw_fold - delta_tw, d1, d2)[0]
    r_parameter = scale * (r_plus - r_minus) / (2.0 * delta_tw)
    transversality = float(np.dot(psi, r_parameter))

    coefficient_rows: list[dict] = []
    # Differentiate the Jacobian in the null-vector direction instead of
    # second-differencing the residual.  The latter subtracts the O(1e-11)
    # collocation residual from O(epsilon**2) quantities and loses digits.
    # Here the residual is quadratic, so the Jacobian derivative is exact up
    # to roundoff and should be independent of epsilon.
    for epsilon in (1.0e-3, 1.0e-4, 1.0e-5):
        jac_plus = vk.newton_residual_jacobian(
            u + epsilon * phi, tw_fold, d1, d2
        )[1]
        jac_minus = vk.newton_residual_jacobian(
            u - epsilon * phi, tw_fold, d1, d2
        )[1]
        second = scale * (((jac_plus - jac_minus) / (2.0 * epsilon)) @ phi)
        quadratic = float(0.5 * np.dot(psi, second))
        coefficient_rows.append(
            {
                "epsilon": epsilon,
                "transversality_coefficient": transversality,
                "quadratic_coefficient": quadratic,
            }
        )

    checks = {
        "degree": degree,
        "newton_residual": float(np.linalg.norm(residual, np.inf)),
        "newton_iterations": len(history),
        "sigma_min_over_max": float(singular[-1] / singular[0]),
        "sigma_second_min_over_max": float(singular[-2] / singular[0]),
        "nullity_gap_sigma_second_over_sigma_min": float(singular[-2] / singular[-1]),
        "transversality_coefficient": transversality,
    }
    return z, fields, phi.reshape(4, z.size), checks, coefficient_rows


def temporal_mode_metrics(z: np.ndarray, vector: np.ndarray) -> tuple[dict, np.ndarray]:
    n = z.size
    modes = vector.reshape(4, n)  # f, g, h, theta
    maxima = np.max(np.abs(modes), axis=1)
    scale = float(np.max(maxima))
    normalised = modes / scale
    amplitude_sq = np.sum(np.abs(normalised) ** 2, axis=0)
    centroid = float(np.trapezoid(z * amplitude_sq, z) / np.trapezoid(amplitude_sq, z))
    active_sq = np.sum(np.abs(normalised[[0, 1, 3]]) ** 2, axis=0)
    active_centroid = float(
        np.trapezoid(z * active_sq, z) / np.trapezoid(active_sq, z)
    )
    peak_coordinates = z[np.argmax(np.abs(normalised), axis=1)]
    metrics = {
        "F_relative_max": float(maxima[0] / scale),
        "G_relative_max": float(maxima[1] / scale),
        "H_relative_max": float(maxima[2] / scale),
        "T_relative_max": float(maxima[3] / scale),
        "coordinate_centroid": centroid,
        "active_fields_centroid": active_centroid,
        "F_peak_coordinate": float(peak_coordinates[0]),
        "G_peak_coordinate": float(peak_coordinates[1]),
        "H_peak_coordinate": float(peak_coordinates[2]),
        "T_peak_coordinate": float(peak_coordinates[3]),
    }
    return metrics, normalised


def scan_rotor_stator_near_folds(collocation_order: int = 101) -> list[dict]:
    """Densely sample both rotor--stator folds without changing source data."""
    folds = json.loads(
        (TWODISK / "boussinesq_singularity_results" /
         "folds_Re1000_traditional_centrifugal.json").read_text(encoding="utf-8")
    )
    offsets = {
        "maximum": np.array(
            [-0.0012, -0.0009, -0.0007, -0.0005, -0.00035, -0.00024,
             -0.00016, -0.00008, 0.0, 0.00008, 0.00016, 0.00024,
             0.00035, 0.0005, 0.0007, 0.0009, 0.0012]
        ),
        "minimum": np.array(
            [-0.00070, -0.00055, -0.00040, -0.00030, -0.00022,
             -0.00015, -0.00010, -0.00005, -0.00003, -0.00002,
             -0.00001, 0.0, 0.00001, 0.00002, 0.00003, 0.00005,
             0.00010, 0.00015, 0.00022, 0.00030, 0.00040, 0.00055,
             0.00070]
        ),
    }
    target_meta: list[tuple[float, str, float]] = []
    for fold in folds:
        pc = float(fold["pressure_gradient"])
        for delta in offsets[fold["kind"]]:
            target_meta.append((pc + float(delta), fold["kind"], float(delta)))
    target_pressures = np.asarray([item[0] for item in target_meta])

    solver = TwoDiskBoussinesqSolver(
        SolverConfig(
            re_h=1000.0,
            pr=0.72,
            tol=2.0e-7,
            initial_nodes=401,
            model="traditional_centrifugal",
        ),
        TWODISK / "baseflow_Res1000.npz",
    )
    isothermal = solver.solve_isothermal()
    initial_pressure = float(isothermal.p[0])
    coarse = np.arange(initial_pressure, 0.015, -0.002)
    moderate = np.arange(coarse[-1] - 0.001, 0.011, -0.001)
    fine = np.arange(moderate[-1] - 0.0005, -0.0011, -0.0005)
    path = np.sort(
        np.unique(np.concatenate((coarse, moderate, fine, target_pressures)))
    )[::-1]
    branch = solver.continue_pressure(path, isothermal, profile_nodes=1001)

    rows: list[dict] = []
    branch_pressures = np.asarray(branch.columns["pressure_gradient"])
    for pressure, fold_kind, offset in target_meta:
        index = int(np.argmin(np.abs(branch_pressures - pressure)))
        if abs(branch_pressures[index] - pressure) > 2.0e-12:
            raise RuntimeError("Requested rotor--stator pressure was not continued")
        values = branch.profiles[index]
        profile = VKProfile(
            branch=0,
            pressure=float(pressure),
            z=branch.z_profile,
            h=values[0],
            f=values[1],
            g=values[3],
            temp=values[5],
        )
        _operator, temporal = temporal_operator(
            profile, collocation_order, 1000.0, 0.72
        )
        leading = temporal.eigvals[0]
        rows.append(
            {
                "fold": fold_kind,
                "pressure_gradient": float(pressure),
                "pressure_offset": offset,
                "Tw": float(branch.columns["Tw"][index]),
                "leading_real": float(leading.real),
                "leading_imag": float(leading.imag),
                "positive_real_count": int(
                    np.count_nonzero(temporal.eigvals.real > 1.0e-7)
                ),
                "fold_sample": int(abs(offset) < 1.0e-15),
                "collocation_order": collocation_order,
            }
        )
    return sorted(rows, key=lambda row: (row["fold"], -row["pressure_gradient"]))


def two_disk_fold_fits(growth: list[dict]) -> list[dict]:
    folds = json.loads(
        (TWODISK / "boussinesq_singularity_results" /
         "folds_Re1000_traditional_centrifugal.json").read_text(encoding="utf-8")
    )
    rows: list[dict] = []
    for fold in folds:
        twc = float(fold["Tw"])
        pic = float(fold["pressure_gradient"])
        is_maximum = fold["kind"] == "maximum"
        fit_cap = 1.0e-4 if is_maximum else 5.0e-5
        for side_name, side_test in (
            ("higher_Pi", lambda pi: pi > pic),
            ("lower_Pi", lambda pi: pi < pic),
        ):
            selected = []
            for item in growth:
                if item["fold"] != fold["kind"]:
                    continue
                mu = (twc - item["Tw"]) if is_maximum else (item["Tw"] - twc)
                if side_test(item["pressure_gradient"]) and 1.0e-9 < mu < fit_cap:
                    selected.append((mu, item["leading_real"] ** 2))
            values = np.asarray(selected)
            slope, r2 = through_origin_fit(values[:, 0], values[:, 1])
            rows.append(
                {
                    "model": "rotor_stator",
                    "fold": fold["kind"],
                    "side": side_name,
                    "points": values.shape[0],
                    "temperature_distance_max": fit_cap,
                    "lambda_squared_slope": slope,
                    "r2_through_origin": r2,
                }
            )
    return rows


def plot_normal_form_scaling(
    vk_rows: list[dict], td_rows: list[dict], output: Path
) -> None:
    folds = json.loads(
        (TWODISK / "boussinesq_singularity_results" /
         "folds_Re1000_traditional_centrifugal.json").read_text(encoding="utf-8")
    )
    fold_map = {item["kind"]: item for item in folds}
    panels = [
        ("von_Karman", vk_rows, "principal"),
        ("rotor_stator", td_rows, "minimum"),
        ("rotor_stator", td_rows, "maximum"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.6))
    for axis, (model, rows, fold_kind) in zip(axes, panels):
        if model == "von_Karman":
            twc = max(row["Tw"] for row in rows)
            selected = rows
            mu = np.asarray([twc - row["Tw"] for row in selected])
            fit_cap = 2.0e-5
            title = "von Karman fold"
        else:
            selected = [row for row in rows if row["fold"] == fold_kind]
            twc = float(fold_map[fold_kind]["Tw"])
            fit_cap = 5.0e-5 if fold_kind == "minimum" else 1.0e-4
            if fold_kind == "minimum":
                mu = np.asarray([row["Tw"] - twc for row in selected])
            else:
                mu = np.asarray([twc - row["Tw"] for row in selected])
            title = f"rotor-stator {fold_kind} fold"
        lam2 = np.asarray([row["leading_real"] ** 2 for row in selected])
        stable = np.asarray([row["leading_real"] <= 0.0 for row in selected])
        use = (mu >= -1.0e-8) & (mu < 1.5e-3)
        axis.scatter(mu[use & stable], lam2[use & stable], color="#15803d", label="stable side")
        axis.scatter(mu[use & ~stable], lam2[use & ~stable], color="#dc2626", label="one-real-mode unstable")
        for mask, colour in ((stable, "#15803d"), (~stable, "#dc2626")):
            fit = use & mask & (mu > 1.0e-9) & (mu < fit_cap)
            slope, _r2 = through_origin_fit(mu[fit], lam2[fit])
            x_line = np.linspace(0.0, fit_cap, 80)
            axis.plot(x_line, slope * x_line, "--", color=colour, linewidth=1.2)
        axis.set_xlabel(r"distance in $T_w$ from fold")
        axis.set_ylabel(r"$[\mathrm{Re}(\lambda)]^2$")
        axis.set_title(title)
        axis.set_xlim(-0.02 * fit_cap, 1.25 * fit_cap)
        visible = (mu >= -1.0e-9) & (mu <= 1.25 * fit_cap)
        y_visible = float(np.max(lam2[visible]))
        axis.set_ylim(-0.03 * y_visible, 1.15 * y_visible)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_vk_fold(rows: list[dict], tw_fold: float, h_fold: float, output: Path):
    h = np.asarray([row["Hinf"] for row in rows])
    tw = np.asarray([row["Tw"] for row in rows])
    growth = np.asarray([row["leading_real"] for row in rows])
    stable = growth < 0.0
    mu = tw_fold - tw

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.7))
    axes[0].plot(tw, h, color="0.65", linewidth=1.2)
    axes[0].scatter(tw[stable], h[stable], color="#15803d", label="stable")
    axes[0].scatter(tw[~stable], h[~stable], color="#dc2626", label="one real unstable mode")
    axes[0].scatter([tw_fold], [h_fold], marker="*", s=140, color="#111827", zorder=5)
    axes[0].set_xlabel(r"$T_w$")
    axes[0].set_ylabel(r"$H_\infty$")
    axes[0].set_title("Fold and stability exchange")
    axes[0].legend(fontsize=8)

    axes[1].plot(h, growth, "o-", color="#1d4ed8", markersize=4)
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[1].axvline(h_fold, color="0.35", linestyle="--")
    axes[1].set_xlabel(r"$H_\infty$")
    axes[1].set_ylabel(r"$\max\,\mathrm{Re}(\lambda)/\Omega$")
    axes[1].set_title("One real eigenvalue crosses zero")

    axes[2].scatter(mu[stable], growth[stable] ** 2, color="#15803d", label="stable side")
    axes[2].scatter(mu[~stable], growth[~stable] ** 2, color="#dc2626", label="unstable side")
    axes[2].set_xlabel(r"$T_{w,c}-T_w$")
    axes[2].set_ylabel(r"$[\mathrm{Re}(\lambda)/\Omega]^2$")
    axes[2].set_title("Saddle-node square-root scaling")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_zero_modes(
    z_vk: np.ndarray,
    vk_mode: np.ndarray,
    z_td: np.ndarray,
    td_mode: np.ndarray,
    output: Path,
):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for values, label in zip(vk_mode, ("F", "G", "H", "temperature")):
        axes[0].plot(z_vk, np.abs(values), label=label)
    axes[0].set_xlabel(r"$\eta$")
    axes[0].set_title("von Karman principal-fold zero mode")
    axes[0].set_xlim(0.0, min(20.0, z_vk[-1]))
    axes[0].set_ylabel("normalised modal amplitude")
    axes[0].legend(fontsize=8)

    for values, label in zip(td_mode, ("F", "G", "H", "temperature")):
        axes[1].plot(z_td, np.abs(values), label=label)
    axes[1].set_xlabel("z/h")
    axes[1].set_title("rotor-stator principal-fold zero mode")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def rotor_stator_principal_fold_mode():
    fold_data = json.loads(
        (TWODISK / "boussinesq_singularity_results" /
         "folds_Re1000_traditional_centrifugal.json").read_text(encoding="utf-8")
    )
    principal = next(item for item in fold_data if item["kind"] == "maximum")
    solver = TwoDiskBoussinesqSolver(
        SolverConfig(re_h=1000.0, pr=0.72, tol=2.0e-7, initial_nodes=401),
        TWODISK / "baseflow_Res1000.npz",
    )
    isothermal = solver.solve_isothermal()
    p0 = float(isothermal.p[0])
    coarse = np.arange(p0, 0.015, -0.002)
    moderate = np.arange(coarse[-1] - 0.001, principal["pressure_gradient"], -0.001)
    pressures = np.sort(
        np.unique(np.concatenate((coarse, moderate, [principal["pressure_gradient"]])))
    )[::-1]
    branch = solver.continue_pressure(pressures, isothermal, profile_nodes=1001)
    values = branch.profiles[-1]
    profile = VKProfile(
        branch=0,
        pressure=float(principal["pressure_gradient"]),
        z=branch.z_profile,
        h=values[0],
        f=values[1],
        g=values[3],
        temp=values[5],
    )
    _operator, result = temporal_operator(profile, 121, 1000.0, 0.72)
    vector = result.eigvecs[:, 0]
    n_a, n_g = result.n_f_reduced, result.n_g
    a = vector[:n_a]
    b = vector[n_a : n_a + n_g]
    c = vector[n_a + n_g :]
    modes = np.vstack(
        (
            result.f_map @ a,
            result.g_map @ b,
            result.h_map @ a,
            result.temp_map @ c,
        )
    )
    metrics, normalised = temporal_mode_metrics(result.z, modes.reshape(-1))
    metrics["leading_real"] = float(result.eigvals[0].real)
    metrics["leading_imag"] = float(result.eigvals[0].imag)
    metrics["Tw"] = float(branch.columns["Tw"][-1])
    metrics["pressure_gradient"] = float(principal["pressure_gradient"])
    return result.z, normalised, metrics


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fold_domains = first_fold_by_domain()
    zmax = 60.0
    h_fold, tw_fold = fold_domains[zmax]

    near_rows, vk_fit_rows, fold_solution = vk_fold_scan(zmax, h_fold, tw_fold)
    write_csv(OUT / "vonkarman_near_fold_temporal.csv", near_rows)
    plot_vk_fold(near_rows, tw_fold, h_fold, OUT / "vonkarman_fold_dynamics.png")

    convergence_rows: list[dict] = []
    for domain in (20.0, 40.0, 60.0, 80.0):
        h_c, tw_c = fold_domains[domain]
        solution = continue_vk(domain, np.array([h_c]), tol=1.0e-8)[round(h_c, 12)]
        degrees = {
            20.0: (60, 80, 100),
            40.0: (80, 110, 140),
            60.0: (100, 130, 160),
            80.0: (120, 150, 180),
        }[domain]
        for degree in degrees:
            leading = leading_vk(solution, domain, degree)
            convergence_rows.append(
                {
                    "zmax": domain,
                    "degree": degree,
                    "Hinf_fold": h_c,
                    "Tw_fold": tw_c,
                    "leading_real": leading.real,
                    "leading_imag": leading.imag,
                }
            )
    write_csv(OUT / "vonkarman_fold_spectrum_convergence.csv", convergence_rows)

    z_steady, _steady_fields, steady_null, jacobian_checks, coefficient_rows = (
        steady_jacobian_checks(fold_solution, zmax, tw_fold, degree=120)
    )
    write_csv(OUT / "vonkarman_nondegeneracy_coefficients.csv", coefficient_rows)
    (OUT / "vonkarman_jacobian_checks.json").write_text(
        json.dumps(jacobian_checks, indent=2), encoding="utf-8"
    )

    z_vk, vk_values, vk_vectors = vk_similarity_spectrum(
        fold_solution, zmax, degree=120
    )
    vk_leading = vk_values[0]
    vk_mode_metrics, vk_mode = temporal_mode_metrics(z_vk, vk_vectors[:, 0])
    vk_mode_metrics.update(
        {
            "leading_real": float(vk_leading.real),
            "leading_imag": float(vk_leading.imag),
            "Tw": tw_fold,
            "Hinf": h_fold,
        }
    )
    steady_mode_reordered = steady_null[[1, 2, 0, 3]]
    steady_mode_reordered /= np.max(np.abs(steady_mode_reordered))
    temporal_real = np.real(vk_mode)
    correlation = abs(
        np.vdot(steady_mode_reordered.reshape(-1), temporal_real.reshape(-1))
    ) / (
        np.linalg.norm(steady_mode_reordered) * np.linalg.norm(temporal_real)
    )
    vk_mode_metrics["steady_temporal_zero_mode_correlation"] = float(correlation)
    (OUT / "vonkarman_zero_mode_metrics.json").write_text(
        json.dumps(vk_mode_metrics, indent=2), encoding="utf-8"
    )
    zero_rows = []
    for index, coordinate in enumerate(z_vk):
        zero_rows.append(
            {
                "eta": float(coordinate),
                "F_mode_abs": float(abs(vk_mode[0, index])),
                "G_mode_abs": float(abs(vk_mode[1, index])),
                "H_mode_abs": float(abs(vk_mode[2, index])),
                "T_mode_abs": float(abs(vk_mode[3, index])),
            }
        )
    write_csv(OUT / "vonkarman_principal_fold_zero_mode.csv", zero_rows)

    td_z, td_mode, td_metrics = rotor_stator_principal_fold_mode()
    (OUT / "rotor_stator_zero_mode_metrics.json").write_text(
        json.dumps(td_metrics, indent=2), encoding="utf-8"
    )
    td_zero_rows = []
    for index, coordinate in enumerate(td_z):
        td_zero_rows.append(
            {
                "z_over_h": float(coordinate),
                "F_mode_abs": float(abs(td_mode[0, index])),
                "G_mode_abs": float(abs(td_mode[1, index])),
                "H_mode_abs": float(abs(td_mode[2, index])),
                "T_mode_abs": float(abs(td_mode[3, index])),
            }
        )
    write_csv(OUT / "rotor_stator_principal_fold_zero_mode.csv", td_zero_rows)
    plot_zero_modes(
        z_vk, vk_mode, td_z, td_mode, OUT / "principal_fold_zero_modes.png"
    )

    td_near_rows = scan_rotor_stator_near_folds(collocation_order=101)
    write_csv(OUT / "rotor_stator_near_fold_temporal.csv", td_near_rows)
    td_fit_rows = two_disk_fold_fits(td_near_rows)
    all_fit_rows = vk_fit_rows + td_fit_rows
    write_csv(OUT / "normal_form_square_root_fits.csv", all_fit_rows)
    plot_normal_form_scaling(
        near_rows,
        td_near_rows,
        OUT / "both_models_square_root_scaling.png",
    )

    comparison = [
        {
            "model": "von_Karman",
            "domain": "semi_infinite_approximated_by_eta_max_60",
            "principal_fold_control": "Hinf",
            "principal_fold_coordinate": h_fold,
            "principal_fold_temperature": tw_fold,
            "fold_temporal_real": float(vk_leading.real),
            "fold_temporal_imag": float(vk_leading.imag),
            "zero_mode_centroid": vk_mode_metrics["coordinate_centroid"],
            "active_fields_centroid": vk_mode_metrics["active_fields_centroid"],
            "global_pressure_parameter": 0,
        },
        {
            "model": "rotor_stator",
            "domain": "finite_gap_Re_h_1000",
            "principal_fold_control": "Pi",
            "principal_fold_coordinate": td_metrics["pressure_gradient"],
            "principal_fold_temperature": td_metrics["Tw"],
            "fold_temporal_real": td_metrics["leading_real"],
            "fold_temporal_imag": td_metrics["leading_imag"],
            "zero_mode_centroid": td_metrics["coordinate_centroid"],
            "active_fields_centroid": td_metrics["active_fields_centroid"],
            "global_pressure_parameter": 1,
        },
    ]
    write_csv(OUT / "model_dynamical_comparison.csv", comparison)

    print(json.dumps(
        {
            "vonkarman_fold": {
                "Hinf": h_fold,
                "Tw": tw_fold,
                "lambda": [vk_leading.real, vk_leading.imag],
                "mode_metrics": vk_mode_metrics,
                "jacobian": jacobian_checks,
                "nondegeneracy": coefficient_rows,
            },
            "rotor_stator_fold": td_metrics,
            "normal_form_fits": all_fit_rows,
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
