#!/usr/bin/env python3
"""Rational-Chebyshev analysis of the heated von Karman upper branch.

The semi-infinite coordinate is compactified with

    eta = L (1-x)/(1+x),  x in [-1, 1],

where x=1 is the disk and x=-1 is infinity.  The script uses H_infinity as
the continuation coordinate, so both temperature folds remain traversable.
Original finite-domain data are read-only; derived output defaults to
``von_karman/data/infinite_mapping`` in this consolidated project.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg
from scipy.interpolate import BarycentricInterpolator, CubicSpline


SCRIPT_DIR = Path(__file__).resolve().parent
VON_KARMAN_ROOT = SCRIPT_DIR.parent
OUT = VON_KARMAN_ROOT / "data" / "infinite_mapping"
sys.path.insert(0, str(SCRIPT_DIR))

import investigate_boussinesq_fold as vk  # noqa: E402


@dataclass
class MappedOperators:
    degree: int
    scale: float
    x: np.ndarray
    eta: np.ndarray
    dx: np.ndarray
    deta: np.ndarray
    deta2: np.ndarray


@dataclass
class MappedSolution:
    h_inf: float
    tw: float
    fields: np.ndarray  # [H,F,G,T] at mapped nodes
    residual_inf: float
    iterations: int


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def rational_chebyshev(degree: int, scale: float) -> MappedOperators:
    """Lobatto matrices ordered from disk (index 0) to infinity."""
    j = np.arange(degree + 1)
    x = np.cos(np.pi * j / degree)
    c = np.r_[2.0, np.ones(degree - 1), 2.0] * (-1.0) ** j
    xx = np.tile(x, (degree + 1, 1)).T
    delta = xx - xx.T
    dx = np.outer(c, 1.0 / c) / (delta + np.eye(degree + 1))
    dx -= np.diag(dx.sum(axis=1))

    eta = np.empty_like(x)
    eta[:-1] = scale * (1.0 - x[:-1]) / (1.0 + x[:-1])
    eta[-1] = np.inf
    dx_deta = -(1.0 + x) ** 2 / (2.0 * scale)
    deta = np.diag(dx_deta) @ dx
    deta2 = deta @ deta
    return MappedOperators(degree, scale, x, eta, dx, deta, deta2)


def mapped_residual_jacobian(
    state: np.ndarray,
    h_target: float,
    operators: MappedOperators,
) -> tuple[np.ndarray, np.ndarray]:
    """Fixed-Hinf residual with Tw as an unknown continuation parameter."""
    n = operators.x.size
    u = state[:-1]
    tw = float(state[-1])
    residual, jacobian = vk.newton_residual_jacobian(
        u, tw, operators.deta, operators.deta2
    )
    h = u[:n]

    # At rational infinity, exponentially decaying fields are flat in x.
    # The transformed continuity equation degenerates to the already imposed
    # F(infinity)=0 condition, so replace that duplicate row by regularity.
    far_continuity_row = n - 1
    residual[far_continuity_row] = operators.dx[-1] @ h
    jacobian[far_continuity_row, :] = 0.0
    jacobian[far_continuity_row, :n] = operators.dx[-1]

    augmented_residual = np.r_[residual, h[-1] - h_target]
    augmented_jacobian = np.zeros((4 * n + 1, 4 * n + 1))
    augmented_jacobian[:-1, :-1] = jacobian
    augmented_jacobian[3 * n, -1] = -1.0  # T(0)-Tw=0
    augmented_jacobian[-1, n - 1] = 1.0   # H(infinity)=H_target
    return augmented_residual, augmented_jacobian


def mapped_newton(
    h_target: float,
    initial_state: np.ndarray,
    operators: MappedOperators,
    tolerance: float = 2.0e-10,
    max_iterations: int = 18,
) -> MappedSolution:
    state = initial_state.copy()
    history: list[float] = []
    for _iteration in range(max_iterations):
        residual, jacobian = mapped_residual_jacobian(
            state, h_target, operators
        )
        norm0 = float(np.linalg.norm(residual, np.inf))
        history.append(norm0)
        if norm0 < tolerance:
            break
        row_norm = np.linalg.norm(jacobian, axis=1)
        row_scale = 1.0 / np.maximum(row_norm, 1.0e-300)
        try:
            step = linalg.solve(
                row_scale[:, None] * jacobian,
                -row_scale * residual,
                assume_a="gen",
            )
        except linalg.LinAlgError as error:
            raise RuntimeError(
                f"Mapped Newton singular at Hinf={h_target:.8g}"
            ) from error

        damping = 1.0
        accepted = False
        while damping >= 2.0e-6:
            trial = state + damping * step
            trial_residual = mapped_residual_jacobian(
                trial, h_target, operators
            )[0]
            if np.linalg.norm(trial_residual, np.inf) < norm0:
                state = trial
                accepted = True
                break
            damping *= 0.5
        if not accepted:
            raise RuntimeError(
                f"Mapped Newton line search failed at Hinf={h_target:.8g}; "
                f"residual={norm0:.3e}"
            )
    else:
        raise RuntimeError(
            f"Mapped Newton did not converge at Hinf={h_target:.8g}; "
            f"residual={history[-1]:.3e}"
        )

    final_residual = mapped_residual_jacobian(
        state, h_target, operators
    )[0]
    n = operators.x.size
    return MappedSolution(
        h_inf=h_target,
        tw=float(state[-1]),
        fields=state[:-1].reshape(4, n),
        residual_inf=float(np.linalg.norm(final_residual, np.inf)),
        iterations=len(history),
    )


def finite_domain_initial_state(
    h_start: float, operators: MappedOperators
) -> np.ndarray:
    zmax = 60.0
    z = np.linspace(0.0, zmax, 721)
    isothermal = vk.solve_isothermal(zmax, z.size, tol=2.0e-9)
    source = vk.solve_fixed_h(
        h_start,
        z,
        isothermal.sol(z),
        1.0,
        tol=2.0e-9,
    )
    sample_eta = np.minimum(operators.eta, zmax)
    source_values = source.sol(sample_eta)[[0, 2, 4, 6]]
    far = operators.eta >= zmax
    source_values[0, far] = h_start
    source_values[1, far] = 0.0
    source_values[2, far] = 1.0
    source_values[3, far] = 1.0
    return np.r_[source_values.reshape(-1), float(source.p[0])]


def default_h_values(h_stop: float) -> np.ndarray:
    pieces = (
        np.arange(-0.75, -0.6001, 0.01),
        np.arange(-0.60, -0.4699, 0.002),
        np.arange(-0.47, -0.1599, 0.005),
        np.arange(-0.16, -0.0599, 0.001),
        np.arange(-0.06, h_stop + 1.0e-10, 0.002),
    )
    values = np.sort(np.unique(np.round(np.concatenate(pieces), 12)))
    return values[(values >= -0.7500001) & (values <= h_stop + 1.0e-10)]


def trace_mapped_branch(
    degree: int,
    scale: float,
    h_stop: float,
    tolerance: float,
) -> tuple[MappedOperators, list[MappedSolution]]:
    operators = rational_chebyshev(degree, scale)
    h_values = default_h_values(h_stop)
    state = finite_domain_initial_state(float(h_values[0]), operators)
    solutions: list[MappedSolution] = []
    for index, h_inf in enumerate(h_values):
        solution = mapped_newton(
            float(h_inf), state, operators, tolerance=tolerance
        )
        solutions.append(solution)
        state = np.r_[solution.fields.reshape(-1), solution.tw]
        if index % 40 == 0:
            print(
                f"N={degree} L={scale:g}: {index + 1}/{h_values.size}, "
                f"Hinf={h_inf:.4f}, Tw={solution.tw:.8f}, "
                f"res={solution.residual_inf:.2e}",
                flush=True,
            )
    return operators, solutions


def locate_folds(solutions: list[MappedSolution]) -> list[dict]:
    h = np.asarray([item.h_inf for item in solutions])
    tw = np.asarray([item.tw for item in solutions])
    spline = CubicSpline(h, tw)
    roots = spline.derivative().roots()
    roots = roots[(roots > h[2]) & (roots < h[-3])]
    rows: list[dict] = []
    for index, root in enumerate(roots, start=1):
        second = float(spline.derivative(2)(root))
        rows.append(
            {
                "fold_index": index,
                "kind": "minimum" if second > 0.0 else "maximum",
                "Hinf": float(root),
                "Tw": float(spline(root)),
                "d2Tw_dHinf2": second,
            }
        )
    return rows


def interpolate_fields(
    operators: MappedOperators,
    solution: MappedSolution,
    eta_query: np.ndarray,
) -> np.ndarray:
    x_query = (operators.scale - eta_query) / (
        operators.scale + eta_query
    )
    values = []
    # BarycentricInterpolator accepts either node ordering.
    for field in solution.fields:
        interpolator = BarycentricInterpolator(operators.x, field)
        values.append(interpolator(x_query))
    return np.asarray(values)


def branch_rows(
    operators: MappedOperators,
    solutions: list[MappedSolution],
    folds: list[dict],
) -> list[dict]:
    fold_h = np.asarray([row["Hinf"] for row in folds])
    rows: list[dict] = []
    for solution in solutions:
        h, f, g, temp = solution.fields
        fp0 = float(operators.deta[0] @ f)
        gp0 = float(operators.deta[0] @ g)
        tp0 = float(operators.deta[0] @ temp)
        rows.append(
            {
                "degree": operators.degree,
                "mapping_scale": operators.scale,
                "Hinf": solution.h_inf,
                "Tw": solution.tw,
                "Fp0": fp0,
                "Gp0": gp0,
                "Tp0": tp0,
                "thermal_tail_length": 1.0 / (vk.PR * abs(solution.h_inf)),
                "residual_inf": solution.residual_inf,
                "newton_iterations": solution.iterations,
                "nearest_fold_distance_Hinf": (
                    float(np.min(np.abs(fold_h - solution.h_inf)))
                    if fold_h.size
                    else float("nan")
                ),
            }
        )
    return rows


def save_selected_profiles(
    operators: MappedOperators,
    solutions: list[MappedSolution],
    output: Path,
) -> None:
    targets = (-0.60, -0.532, -0.30, -0.15, -0.10, -0.08, -0.04, -0.02)
    eta_plot = np.r_[
        np.linspace(0.0, 20.0, 401),
        np.geomspace(20.1, 2000.0, 500),
    ]
    rows: list[dict] = []
    for target in targets:
        solution = min(solutions, key=lambda item: abs(item.h_inf - target))
        values = interpolate_fields(operators, solution, eta_plot)
        for index, eta in enumerate(eta_plot):
            rows.append(
                {
                    "Hinf": solution.h_inf,
                    "Tw": solution.tw,
                    "eta": float(eta),
                    "H": float(values[0, index]),
                    "F": float(values[1, index]),
                    "G": float(values[2, index]),
                    "T": float(values[3, index]),
                }
            )
    write_csv(output, rows)


def plot_branch(
    solution_sets: list[tuple[MappedOperators, list[MappedSolution], list[dict]]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for operators, solutions, folds in solution_sets:
        h = np.asarray([item.h_inf for item in solutions])
        tw = np.asarray([item.tw for item in solutions])
        label = f"N={operators.degree}, L={operators.scale:g}"
        axes[0].plot(tw, h, linewidth=1.6, label=label)
        axes[1].plot(h, tw, linewidth=1.6, label=label)
        for fold in folds:
            axes[0].scatter(fold["Tw"], fold["Hinf"], marker="*", s=90)
            axes[1].scatter(fold["Hinf"], fold["Tw"], marker="*", s=90)
    axes[0].set_xlabel(r"$T_w$")
    axes[0].set_ylabel(r"$H_\infty$")
    axes[0].set_title("Infinite-mapped branch")
    axes[1].set_xlabel(r"$H_\infty$")
    axes[1].set_ylabel(r"$T_w$")
    axes[1].set_title("Turning points in continuation coordinate")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def interpolate_solution_to_operators(
    source_operators: MappedOperators,
    solution: MappedSolution,
    target_operators: MappedOperators,
) -> np.ndarray:
    values = []
    for field in solution.fields:
        interpolator = BarycentricInterpolator(source_operators.x, field)
        values.append(interpolator(target_operators.x))
    result = np.asarray(values)
    result[0, -1] = solution.h_inf
    result[1, -1] = 0.0
    result[2, -1] = 1.0
    result[3, -1] = 1.0
    return result


def mapped_similarity_spectrum(
    source_operators: MappedOperators,
    solution: MappedSolution,
    degree: int,
) -> np.ndarray:
    """Similarity-preserving temporal spectrum on the rational grid."""
    operators = rational_chebyshev(degree, source_operators.scale)
    fields = interpolate_solution_to_operators(
        source_operators, solution, operators
    )
    h_base, f_base, g_base, temp_base = fields
    d1, d2, dx = operators.deta, operators.deta2, operators.dx
    fp, gp, tp = d1 @ f_base, d1 @ g_base, d1 @ temp_base
    n = operators.x.size
    identity = np.eye(n)
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
    # h approaches a freely selected constant at infinity; regularity of the
    # mapped perturbation replaces the degenerate far continuity row.
    far_h_row = 3 * n - 1
    matrix_a[far_h_row, :] = 0.0
    matrix_b[far_h_row, :] = 0.0
    matrix_a[far_h_row, sh] = dx[-1]
    set_bc(3 * n, 3 * n)
    set_bc(4 * n - 1, 4 * n - 1)

    values = linalg.eig(matrix_a, matrix_b, right=False)
    values = values[np.isfinite(values) & (np.abs(values) < 1.0e4)]
    order = np.lexsort((np.abs(values.imag), -values.real))
    return values[order]


def temporal_sample_rows(
    operators: MappedOperators,
    solutions: list[MappedSolution],
    folds: list[dict],
    temporal_degrees: list[int],
) -> list[dict]:
    # Strict-residual N=50--70, L=8--12 mapped runs give more reliable fold
    # coordinates than the very high-order double-precision spline alone.
    reference_fold_h = {
        1: -0.5327616566962796,
        2: -0.11330573101124149,
    }
    targets = [
        ("principal", -0.60),
        ("middle", -0.30),
        ("middle_near_second", -0.15),
        ("upper_after_second", -0.10),
        ("upper", -0.08),
        ("upper", -0.05),
        ("near_zero_limit", -0.02),
    ]
    for fold in folds[:3]:
        fold_index = int(fold["fold_index"])
        target = reference_fold_h.get(fold_index, float(fold["Hinf"]))
        targets.append((f"fold_{fold_index}", target))

    rows: list[dict] = []
    for label, target in targets:
        nearest = min(solutions, key=lambda item: abs(item.h_inf - target))
        initial = np.r_[nearest.fields.reshape(-1), nearest.tw]
        exact = mapped_newton(
            target, initial, operators, tolerance=3.0e-10
        )
        for degree in temporal_degrees:
            values = mapped_similarity_spectrum(operators, exact, degree)
            near_zero_index = int(np.argmin(np.abs(values)))
            positive = int(np.count_nonzero(values.real > 1.0e-7))
            for rank, eigenvalue in enumerate(values[:12], start=1):
                rows.append(
                    {
                        "sample": label,
                        "Hinf": target,
                        "Tw": exact.tw,
                        "temporal_degree": degree,
                        "rank_by_real_part": rank,
                        "eigen_real": float(eigenvalue.real),
                        "eigen_imag": float(eigenvalue.imag),
                        "positive_real_count": positive,
                        "closest_to_zero_real": float(values[near_zero_index].real),
                        "closest_to_zero_imag": float(values[near_zero_index].imag),
                    }
                )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, nargs="+", default=[100])
    parser.add_argument("--scale", type=float, nargs="+", default=[8.0])
    parser.add_argument("--h-stop", type=float, default=-0.02)
    parser.add_argument("--tolerance", type=float, default=2.0e-10)
    parser.add_argument(
        "--temporal-degrees", type=int, nargs="+", default=[]
    )
    parser.add_argument("--output-dir", type=Path, default=OUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    solution_sets = []
    fold_summary: list[dict] = []
    for degree in args.degree:
        for scale in args.scale:
            operators, solutions = trace_mapped_branch(
                degree, scale, args.h_stop, args.tolerance
            )
            folds = locate_folds(solutions)
            solution_sets.append((operators, solutions, folds))
            rows = branch_rows(operators, solutions, folds)
            tag = f"N{degree}_L{scale:g}".replace(".", "p")
            write_csv(args.output_dir / f"branch_{tag}.csv", rows)
            save_selected_profiles(
                operators,
                solutions,
                args.output_dir / f"profiles_{tag}.csv",
            )
            for fold in folds:
                fold_summary.append(
                    {
                        "degree": degree,
                        "mapping_scale": scale,
                        **fold,
                    }
                )
    write_csv(args.output_dir / "fold_convergence.csv", fold_summary)
    plot_branch(
        solution_sets,
        args.output_dir / "infinite_mapped_branch_convergence.png",
    )
    (args.output_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "mapping": "eta=L(1-x)/(1+x)",
                "degrees": args.degree,
                "scales": args.scale,
                "h_stop": args.h_stop,
                "folds": fold_summary,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    if args.temporal_degrees:
        operators, solutions, folds = solution_sets[-1]
        temporal_rows = temporal_sample_rows(
            operators, solutions, folds, args.temporal_degrees
        )
        write_csv(args.output_dir / "mapped_temporal_samples.csv", temporal_rows)
    print(json.dumps(fold_summary, indent=2))


if __name__ == "__main__":
    main()
