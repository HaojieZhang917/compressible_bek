#!/usr/bin/env python3
"""Continuation analysis of thermally stratified rotor--stator similarity flow.

The isothermal equations and normalisation follow ``BaseFlow_cavity.jl`` and
the attached two-disk manuscript.  The thermal extension retains the unknown
constant radial pressure gradient instead of differentiating the radial
momentum equation.  This is essential: differentiating first would remove the
constant part of the centrifugal-buoyancy term and can obscure folds.

The primary model is the same centrifugal-only, traditional Boussinesq model
used in the existing half-space comparison in this repository.  With
theta=T-1 and beta*T_ref=1,

    H'   = -2 sqrt(Re_h) F,
    F''  = sqrt(Re_h) H F' + Re_h (F^2-G^2+theta+Pi),
    G''  = sqrt(Re_h) H G' + 2 Re_h F G,
    T''  = Pr sqrt(Re_h) H T'.

Here G is the absolute azimuthal velocity divided by Omega*r and Pi is the
unknown dimensionless radial-pressure-gradient constant.  Boundary conditions
are H=F=0 at both disks, (G,T)=(1,Tw) at the rotor z=0, and (G,T)=(0,1) at the
stator z=1.

A secondary ``soong_rotating_forces`` option is included only as a model-choice
sensitivity check.  It retains density changes in both centrifugal and
Coriolis forces after conversion from Soong's rotating-frame variables.  It is
not the Blackburn et al. (2021) canonical unsteady model.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


@dataclass(frozen=True)
class SolverConfig:
    re_h: float = 1000.0
    pr: float = 0.72
    tol: float = 2.0e-7
    initial_nodes: int = 401
    max_nodes: int = 120_000
    model: str = "traditional_centrifugal"


@dataclass
class BranchData:
    columns: dict[str, np.ndarray]
    z_profile: np.ndarray
    profiles: np.ndarray


class TwoDiskBoussinesqSolver:
    def __init__(self, config: SolverConfig, reference_npz: Path):
        self.config = config
        self.reference_npz = reference_npz
        if config.model not in {
            "traditional_centrifugal",
            "soong_rotating_forces",
        }:
            raise ValueError(f"Unknown model: {config.model}")

    @property
    def sqrt_re(self) -> float:
        return float(np.sqrt(self.config.re_h))

    def _derivatives(self, y: np.ndarray, pressure: float) -> np.ndarray:
        h, f, fp, g, gp, temp, tp = y
        theta = temp - 1.0
        re_h = self.config.re_h
        sqrt_re = self.sqrt_re

        if self.config.model == "traditional_centrifugal":
            radial_buoyancy = theta
            azimuthal_factor = g
        else:
            radial_buoyancy = theta * (2.0 * g - 1.0)
            azimuthal_factor = g - theta

        return np.vstack(
            (
                -2.0 * sqrt_re * f,
                fp,
                sqrt_re * h * fp
                + re_h * (f * f - g * g + radial_buoyancy + pressure),
                gp,
                sqrt_re * h * gp + 2.0 * re_h * f * azimuthal_factor,
                tp,
                self.config.pr * sqrt_re * h * tp,
            )
        )

    def _reference_guess(self) -> tuple[np.ndarray, np.ndarray]:
        ref = np.load(self.reference_npz)
        z_ref = np.asarray(ref["x"]).ravel()
        z = np.linspace(0.0, 1.0, self.config.initial_nodes)
        y = np.zeros((7, z.size))
        y[0] = np.interp(z, z_ref, np.asarray(ref["w0"]).ravel())
        y[1] = np.interp(z, z_ref, np.asarray(ref["u0"]).ravel())
        y[3] = np.interp(z, z_ref, np.asarray(ref["v0"]).ravel())
        y[2] = np.gradient(y[1], z)
        y[4] = np.gradient(y[3], z)
        y[5] = 1.0
        return z, y

    def solve_isothermal(self):
        z, y = self._reference_guess()

        def ode(_z, state, parameter):
            return self._derivatives(state, float(parameter[0]))

        def bc(ya, yb, _parameter):
            return np.array(
                [
                    ya[0],
                    ya[1],
                    ya[3] - 1.0,
                    ya[5] - 1.0,
                    yb[0],
                    yb[1],
                    yb[3],
                    yb[5] - 1.0,
                ]
            )

        solution = solve_bvp(
            ode,
            bc,
            z,
            y,
            p=np.array([0.1]),
            tol=self.config.tol,
            max_nodes=self.config.max_nodes,
        )
        self._require_success(solution, "isothermal validation")
        return solution

    def validation_errors(self, solution) -> dict[str, float | bool]:
        errors: dict[str, float | bool] = {
            "reference_re_h": 1000.0,
            "reference_comparison_applicable": bool(
                np.isclose(self.config.re_h, 1000.0)
            ),
        }
        if errors["reference_comparison_applicable"]:
            ref = np.load(self.reference_npz)
            z_ref = np.asarray(ref["x"]).ravel()
            computed = solution.sol(z_ref)
            mapping = {"H": (0, "w0"), "F": (1, "u0"), "G": (3, "v0")}
            for name, (index, key) in mapping.items():
                delta = computed[index] - np.asarray(ref[key]).ravel()
                errors[f"{name}_max_abs"] = float(np.max(np.abs(delta)))
                errors[f"{name}_rms"] = float(np.sqrt(np.mean(delta * delta)))
        errors["pressure_gradient"] = float(solution.p[0])
        errors["max_rms_residual"] = float(np.max(solution.rms_residuals))
        errors["nodes"] = int(solution.x.size)
        return errors

    def continue_pressure(
        self,
        pressure_values: Iterable[float],
        initial_solution,
        profile_nodes: int = 401,
    ) -> BranchData:
        z = initial_solution.x
        y = initial_solution.y
        tw_parameter = np.array([1.0])
        z_profile = np.linspace(0.0, 1.0, profile_nodes)

        rows: list[dict[str, float]] = []
        profiles: list[np.ndarray] = []

        def bc(ya, yb, parameter):
            return np.array(
                [
                    ya[0],
                    ya[1],
                    ya[3] - 1.0,
                    ya[5] - parameter[0],
                    yb[0],
                    yb[1],
                    yb[3],
                    yb[5] - 1.0,
                ]
            )

        for pressure in pressure_values:
            def ode(_z, state, _parameter, pressure=pressure):
                return self._derivatives(state, float(pressure))

            solution = solve_bvp(
                ode,
                bc,
                z,
                y,
                p=tw_parameter,
                tol=self.config.tol,
                max_nodes=self.config.max_nodes,
            )
            self._require_success(
                solution,
                f"pressure continuation at Pi={pressure:.9g}",
            )
            z = solution.x
            y = solution.y
            tw_parameter = solution.p.copy()
            dense = solution.sol(z_profile)
            rows.append(
                {
                    "pressure_gradient": float(pressure),
                    "Tw": float(solution.p[0]),
                    "thermal_rossby": float(solution.p[0] - 1.0),
                    "rotor_heat_flux_Tz": float(solution.y[6, 0]),
                    "G_mid": float(solution.sol(0.5)[3]),
                    "H_min": float(np.min(dense[0])),
                    "F_min": float(np.min(dense[1])),
                    "F_max": float(np.max(dense[1])),
                    "nodes": int(solution.x.size),
                    "max_rms_residual": float(np.max(solution.rms_residuals)),
                }
            )
            profiles.append(dense)

        columns = {
            key: np.asarray([row[key] for row in rows]) for key in rows[0]
        }
        return BranchData(columns, z_profile, np.asarray(profiles))

    def solve_fixed_temperature(self, tw: float, seed_pressure: float, seed_y: np.ndarray):
        z = np.linspace(0.0, 1.0, seed_y.shape[1])

        def ode(_z, state, parameter):
            return self._derivatives(state, float(parameter[0]))

        def bc(ya, yb, _parameter):
            return np.array(
                [
                    ya[0],
                    ya[1],
                    ya[3] - 1.0,
                    ya[5] - tw,
                    yb[0],
                    yb[1],
                    yb[3],
                    yb[5] - 1.0,
                ]
            )

        solution = solve_bvp(
            ode,
            bc,
            z,
            seed_y,
            p=np.array([seed_pressure]),
            tol=self.config.tol,
            max_nodes=self.config.max_nodes,
        )
        self._require_success(solution, f"fixed Tw={tw:g}")
        return solution

    def continue_temperature(
        self,
        temperature_values: Iterable[float],
        initial_solution,
        profile_nodes: int = 401,
    ) -> BranchData:
        z = initial_solution.x
        y = initial_solution.y
        pressure_parameter = initial_solution.p.copy()
        z_profile = np.linspace(0.0, 1.0, profile_nodes)
        rows: list[dict[str, float]] = []
        profiles: list[np.ndarray] = []

        def ode(_z, state, parameter):
            return self._derivatives(state, float(parameter[0]))

        for tw in temperature_values:
            def bc(ya, yb, _parameter, tw=tw):
                return np.array(
                    [
                        ya[0],
                        ya[1],
                        ya[3] - 1.0,
                        ya[5] - tw,
                        yb[0],
                        yb[1],
                        yb[3],
                        yb[5] - 1.0,
                    ]
                )

            solution = solve_bvp(
                ode,
                bc,
                z,
                y,
                p=pressure_parameter,
                tol=self.config.tol,
                max_nodes=self.config.max_nodes,
            )
            self._require_success(solution, f"temperature continuation at Tw={tw:g}")
            z = solution.x
            y = solution.y
            pressure_parameter = solution.p.copy()
            dense = solution.sol(z_profile)
            rows.append(
                {
                    "pressure_gradient": float(solution.p[0]),
                    "Tw": float(tw),
                    "thermal_rossby": float(tw - 1.0),
                    "rotor_heat_flux_Tz": float(solution.y[6, 0]),
                    "G_mid": float(solution.sol(0.5)[3]),
                    "H_min": float(np.min(dense[0])),
                    "F_min": float(np.min(dense[1])),
                    "F_max": float(np.max(dense[1])),
                    "nodes": int(solution.x.size),
                    "max_rms_residual": float(np.max(solution.rms_residuals)),
                }
            )
            profiles.append(dense)

        columns = {
            key: np.asarray([row[key] for row in rows]) for key in rows[0]
        }
        return BranchData(columns, z_profile, np.asarray(profiles))

    @staticmethod
    def _require_success(solution, context: str) -> None:
        if solution.status != 0:
            raise RuntimeError(
                f"BVP failure during {context}: {solution.message}; "
                f"nodes={solution.x.size}"
            )


def pressure_grid(initial_pressure: float, fine_step: float) -> np.ndarray:
    """Monotone grid with extra resolution around the observed folds."""
    coarse_high = np.arange(initial_pressure, 0.015, -0.001)
    fine = np.arange(coarse_high[-1] - fine_step, -0.006 - fine_step / 2, -fine_step)
    coarse_low = np.arange(fine[-1] - 0.0005, -0.0501, -0.0005)
    return np.concatenate((coarse_high, fine, coarse_low))


def turning_points(branch: BranchData) -> list[dict[str, float | str]]:
    pressure = branch.columns["pressure_gradient"]
    tw = branch.columns["Tw"]
    order = np.argsort(pressure)
    p_sorted = pressure[order]
    tw_sorted = tw[order]
    spline = CubicSpline(p_sorted, tw_sorted)
    roots = spline.derivative().roots()
    roots = roots[(roots > p_sorted[0]) & (roots < p_sorted[-1])]
    roots.sort()

    results: list[dict[str, float | str]] = []
    span = p_sorted[-1] - p_sorted[0]
    eps = min(2.0e-5, 1.0e-3 * span)
    for root in roots:
        left = float(spline.derivative()(root - eps))
        right = float(spline.derivative()(root + eps))
        if left * right >= 0.0:
            continue
        kind = "minimum" if left < 0.0 < right else "maximum"
        result: dict[str, float | str] = {
            "kind": kind,
            "pressure_gradient": float(root),
            "Tw": float(spline(root)),
            "thermal_rossby": float(spline(root) - 1.0),
        }
        for key in ("rotor_heat_flux_Tz", "G_mid", "H_min", "F_min", "F_max"):
            result[key] = float(
                CubicSpline(p_sorted, branch.columns[key][order])(root)
            )
        results.append(result)
    return results


def roots_at_temperature(branch: BranchData, target_tw: float) -> list[float]:
    pressure = branch.columns["pressure_gradient"]
    tw = branch.columns["Tw"]
    order = np.argsort(pressure)
    p = pressure[order]
    values = tw[order] - target_tw
    spline = CubicSpline(p, tw[order])
    roots: list[float] = []
    for left, right, f_left, f_right in zip(p[:-1], p[1:], values[:-1], values[1:]):
        if f_left == 0.0:
            roots.append(float(left))
        elif f_left * f_right < 0.0:
            roots.append(float(brentq(lambda value: spline(value) - target_tw, left, right)))
    return sorted(set(np.round(roots, 12)))


def write_csv(path: Path, columns: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for index in range(len(columns[names[0]])):
            writer.writerow([columns[name][index] for name in names])


def plot_branch(
    branch: BranchData,
    folds: list[dict],
    re_h: float,
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), constrained_layout=True)
    pressure = branch.columns["pressure_gradient"]
    tw = branch.columns["Tw"]
    axes[0].plot(pressure, tw, color="#1f4e79", lw=2.0)
    axes[1].plot(pressure, tw, color="#1f4e79", lw=2.0)
    for fold in folds:
        for axis in axes:
            axis.scatter(
                fold["pressure_gradient"],
                fold["Tw"],
                color="#c43c35",
                zorder=4,
            )
        offset = (-92, -30) if fold["kind"] == "minimum" else (18, 12)
        axes[1].annotate(
            f"{fold['kind']}\n$T_w$={fold['Tw']:.6f}",
            (fold["pressure_gradient"], fold["Tw"]),
            xytext=offset,
            textcoords="offset points",
            fontsize=9,
            arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.8},
        )

    axes[0].set_xlabel(r"radial-pressure parameter $\Pi$")
    axes[1].set_xlabel(r"radial-pressure parameter $\Pi$")
    axes[1].set_xlim(-0.002, 0.022)
    axes[1].set_ylim(1.145, 1.175)
    axes[1].set_title("fold region")
    for axis in axes:
        axis.set_ylabel(r"rotor temperature ratio $T_w$")
        axis.grid(alpha=0.25)
    fig.suptitle(rf"Traditional Boussinesq two-disk branch ($Re_h={re_h:g}$)")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_three_solutions(
    solutions: list,
    pressures: list[float],
    tw: float,
    path: Path,
) -> None:
    z = np.linspace(0.0, 1.0, 1001)
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0), constrained_layout=True)
    variables = [(1, "F"), (3, "G"), (0, "H"), (5, "T")]
    colors = ["#1f4e79", "#d98e04", "#8b2c6f"]
    for solution, pressure, color in zip(solutions, pressures, colors):
        values = solution.sol(z)
        label = rf"$\Pi={pressure:.6f}$"
        for axis, (index, name) in zip(axes.flat, variables):
            axis.plot(z, values[index], lw=1.8, color=color, label=label)
            axis.set_xlabel("z/h")
            axis.set_ylabel(name)
            axis.grid(alpha=0.22)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(rf"Three steady similarity solutions at $T_w={tw:.3f}$")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_temperature_continuation(branch: BranchData, re_h: float, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    tw = branch.columns["Tw"]
    axes[0].plot(tw, branch.columns["pressure_gradient"], lw=2.0, color="#375a7f")
    axes[1].plot(tw, branch.columns["G_mid"], lw=2.0, color="#8b2c6f")
    axes[0].set_ylabel(r"radial-pressure parameter $\Pi$")
    axes[1].set_ylabel(r"mid-gap swirl $G(1/2)$")
    for axis in axes:
        axis.set_xlabel(r"rotor temperature ratio $T_w$")
        axis.grid(alpha=0.25)
    fig.suptitle(rf"Soong rotating-force Boussinesq check ($Re_h={re_h:g}$)")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_analysis(args: argparse.Namespace) -> dict:
    script_dir = Path(__file__).resolve().parent
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_npz = script_dir / "baseflow_Res1000.npz"

    config = SolverConfig(
        re_h=args.re_h,
        pr=args.pr,
        tol=args.tol,
        initial_nodes=args.initial_nodes,
        model=args.model,
    )
    solver = TwoDiskBoussinesqSolver(config, reference_npz)
    isothermal = solver.solve_isothermal()
    validation = solver.validation_errors(isothermal)
    if args.model == "traditional_centrifugal":
        grid = pressure_grid(float(isothermal.p[0]), args.fine_pressure_step)
        branch = solver.continue_pressure(grid, isothermal)
        folds = turning_points(branch)
    else:
        temperatures = np.arange(1.0, args.alternative_max_temperature + 0.005, 0.01)
        branch = solver.continue_temperature(temperatures, isothermal)
        folds = []

    write_csv(output_dir / f"branch_Re{args.re_h:g}_{args.model}.csv", branch.columns)
    with (output_dir / f"folds_Re{args.re_h:g}_{args.model}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(folds, handle, indent=2, ensure_ascii=False)
    with (output_dir / f"validation_Re{args.re_h:g}_{args.model}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(validation, handle, indent=2, ensure_ascii=False)

    if args.model == "traditional_centrifugal":
        plot_branch(
            branch,
            folds,
            args.re_h,
            output_dir / f"traditional_branch_Re{args.re_h:g}.png",
        )
        target_tw = args.profile_temperature
        root_pressures = roots_at_temperature(branch, target_tw)
        solutions = []
        for pressure in root_pressures:
            nearest = int(
                np.argmin(np.abs(branch.columns["pressure_gradient"] - pressure))
            )
            solution = solver.solve_fixed_temperature(
                target_tw,
                pressure,
                branch.profiles[nearest],
            )
            solutions.append(solution)
        if solutions:
            plot_three_solutions(
                solutions,
                root_pressures,
                target_tw,
                output_dir / f"three_solutions_Tw{target_tw:.3f}.png",
            )
            profile_rows = []
            z = np.linspace(0.0, 1.0, 1001)
            for branch_index, (solution, pressure) in enumerate(
                zip(solutions, root_pressures), start=1
            ):
                values = solution.sol(z)
                for j, location in enumerate(z):
                    profile_rows.append(
                        {
                            "branch": branch_index,
                            "pressure_gradient": pressure,
                            "z": location,
                            "H": values[0, j],
                            "F": values[1, j],
                            "G": values[3, j],
                            "T": values[5, j],
                        }
                    )
            profile_columns = {
                key: np.asarray([row[key] for row in profile_rows])
                for key in profile_rows[0]
            }
            write_csv(output_dir / f"three_solutions_Tw{target_tw:.3f}.csv", profile_columns)
    else:
        plot_temperature_continuation(
            branch,
            args.re_h,
            output_dir / f"soong_temperature_continuation_Re{args.re_h:g}.png",
        )
    return {
        "config": vars(args),
        "validation": validation,
        "folds": folds,
        "branch_points": int(len(branch.columns["Tw"])),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--re-h", type=float, default=1000.0)
    parser.add_argument("--pr", type=float, default=0.72)
    parser.add_argument("--tol", type=float, default=2.0e-7)
    parser.add_argument("--initial-nodes", type=int, default=401)
    parser.add_argument("--fine-pressure-step", type=float, default=5.0e-5)
    parser.add_argument("--profile-temperature", type=float, default=1.16)
    parser.add_argument("--alternative-max-temperature", type=float, default=2.0)
    parser.add_argument(
        "--model",
        choices=["traditional_centrifugal", "soong_rotating_forces"],
        default="traditional_centrifugal",
    )
    parser.add_argument(
        "--output-dir",
        default="boussinesq_singularity_results",
    )
    return parser.parse_args()


if __name__ == "__main__":
    summary = run_analysis(parse_args())
    print(json.dumps(summary, indent=2, ensure_ascii=False))
