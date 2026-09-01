#!/usr/bin/env python3
"""Compare Blackburn and fully compressible neutral curves.

The critical-point extraction follows the final two slides of Boussinesq.pptx:
interior minima of R(beta) are refined with a three-point quadratic fit and
classified as Type-I or Type-II using beta=0.055.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20)
TYPE_BETA_SPLIT = 0.055


@dataclass(frozen=True)
class CriticalPoint:
    mode: str
    R: float
    beta: float
    R_discrete: float
    beta_discrete: float
    fitted: bool
    segment: int
    index: int


def number_tag(value: float) -> str:
    text = f"{value:.4f}".rstrip("0")
    if text.endswith("."):
        return text + "0"
    return text


def load_curve(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if (
            not line
            or line.lower().startswith("variables")
            or line.lower().startswith("zone")
            or line.startswith("#")
        ):
            continue
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(f"Expected seven columns in {path}: {line}")
        rows.append([float(field) for field in fields])
    if not rows:
        raise ValueError(f"No curve rows found in {path}")
    return np.asarray(rows, dtype=float)


def blackburn_paths(directory: Path, Tw: float) -> list[Path]:
    stem = f"ome=0.0_Tw={number_tag(Tw)}_model=blackburn"
    paths = [directory / f"{stem}.dat"]
    branch = directory / f"{stem}_branch=typeII.dat"
    if branch.exists():
        paths.append(branch)
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
    return paths


def compressible_path(directory: Path, Tw: float) -> Path:
    return directory / (
        f"ome=0.0_Tw={number_tag(Tw)}_model=compressible_Mr=0.3_"
        "propPert=on_baseProp=variable.dat"
    )


def validate_segment(data: np.ndarray, path: Path) -> dict[str, float | int]:
    if data.shape[1] != 7:
        raise ValueError(f"Expected seven columns in {path}")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite values in {path}")
    beta_diff = np.diff(data[:, 2])
    if not (np.all(beta_diff > 0) or np.all(beta_diff < 0)):
        raise ValueError(f"beta is not strictly monotone in {path}")
    residual = np.minimum(np.abs(data[:, 4]), np.abs(data[:, 6]))
    maximum_residual = float(np.max(residual))
    if maximum_residual > 1.0e-6:
        raise ValueError(
            f"Neutral residual {maximum_residual:.3e} exceeds 1e-6 in {path}"
        )
    return {
        "points": int(data.shape[0]),
        "R_min": float(np.min(data[:, 1])),
        "R_max": float(np.max(data[:, 1])),
        "beta_min": float(np.min(data[:, 2])),
        "beta_max": float(np.max(data[:, 2])),
        "max_residual": maximum_residual,
    }


def refine_minimum(data: np.ndarray, index: int) -> tuple[float, float, bool]:
    beta0 = data[index, 2]
    x = data[index - 1 : index + 2, 2] - beta0
    y = data[index - 1 : index + 2, 1]
    a, b, c = np.linalg.solve(
        np.column_stack((x**2, x, np.ones(3))), y
    )
    if not np.all(np.isfinite((a, b, c))) or a <= 0:
        return float(data[index, 1]), float(beta0), False
    vertex = -b / (2 * a)
    if vertex < np.min(x) or vertex > np.max(x):
        return float(data[index, 1]), float(beta0), False
    R_vertex = a * vertex**2 + b * vertex + c
    if not np.isfinite(R_vertex):
        return float(data[index, 1]), float(beta0), False
    return float(R_vertex), float(beta0 + vertex), True


def critical_points(segments: list[np.ndarray]) -> dict[str, CriticalPoint]:
    candidates: list[CriticalPoint] = []
    for segment_index, data in enumerate(segments, start=1):
        for index in range(1, data.shape[0] - 1):
            if data[index, 1] >= data[index - 1, 1]:
                continue
            if data[index, 1] >= data[index + 1, 1]:
                continue
            R, beta, fitted = refine_minimum(data, index)
            mode = "Type-I" if beta >= TYPE_BETA_SPLIT else "Type-II"
            candidates.append(
                CriticalPoint(
                    mode=mode,
                    R=R,
                    beta=beta,
                    R_discrete=float(data[index, 1]),
                    beta_discrete=float(data[index, 2]),
                    fitted=fitted,
                    segment=segment_index,
                    index=index + 1,
                )
            )
    selected: dict[str, CriticalPoint] = {}
    for mode in ("Type-I", "Type-II"):
        choices = [point for point in candidates if point.mode == mode]
        if choices:
            selected[mode] = min(choices, key=lambda point: point.R)
    return selected


def write_combined_tecplot(
    path: Path,
    cases: dict[tuple[float, str], list[np.ndarray]],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write('TITLE = "Blackburn and compressible neutral curves"\n')
        stream.write(
            'VARIABLES = "omega" "R" "beta" "alpha_r_1" '
            '"alpha_i_1" "alpha_r_2" "alpha_i_2"\n'
        )
        for (Tw, model), segments in cases.items():
            for segment_index, data in enumerate(segments, start=1):
                stream.write(
                    f'ZONE T="Tw={Tw:.2f} {model} segment={segment_index}", '
                    f"I={data.shape[0]}, DATAPACKING=POINT\n"
                )
                np.savetxt(stream, data, fmt="%.12e")


def plot_curves(
    path: Path,
    cases: dict[tuple[float, str], list[np.ndarray]],
    critical: dict[tuple[float, str], dict[str, CriticalPoint]],
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), sharex=True, sharey=True)
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        for model, color, label in (
            ("compressible", "black", "Compressible"),
            ("blackburn", "#d62728", "Blackburn"),
        ):
            for segment_index, data in enumerate(cases[(Tw, model)]):
                axis.plot(
                    data[:, 1],
                    data[:, 2],
                    color=color,
                    linewidth=1.6,
                    label=label if segment_index == 0 else None,
                )
            for mode, marker in (("Type-I", "o"), ("Type-II", "s")):
                point = critical[(Tw, model)].get(mode)
                if point is not None:
                    axis.scatter(
                        point.R,
                        point.beta,
                        s=22,
                        marker=marker,
                        facecolor=color,
                        edgecolor="white",
                        linewidth=0.45,
                        zorder=4,
                    )
        axis.text(
            0.04,
            0.93,
            rf"$T_w={Tw:.2f}$",
            transform=axis.transAxes,
            fontsize=11,
            va="top",
        )
        axis.grid(alpha=0.22, linewidth=0.6)
        axis.set_xlim(200, 525)
        axis.set_ylim(0.015, 0.135)
    axes[0, 0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.82))
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$R$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\beta$")
    fig.suptitle("Neutral curves: Blackburn vs fully compressible model", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_errors(path: Path, rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3), sharex=True)
    for mode, marker, color in (
        ("Type-I", "o", "#1f77b4"),
        ("Type-II", "s", "#d62728"),
    ):
        available = [
            row
            for row in rows
            if row["mode"] == mode and row["status"] == "available"
        ]
        temperatures = [float(row["Tw"]) for row in available]
        axes[0].plot(
            temperatures,
            [float(row["R_error_vs_compressible_percent"]) for row in available],
            marker=marker,
            color=color,
            label=mode,
        )
        axes[1].plot(
            temperatures,
            [
                float(row["beta_error_vs_compressible_percent"])
                for row in available
            ],
            marker=marker,
            color=color,
            label=mode,
        )
    axes[0].set_ylabel(r"$100(R_{c,B}-R_{c,C})/R_{c,C}$ [%]")
    axes[1].set_ylabel(r"$100(\beta_{c,B}-\beta_{c,C})/\beta_{c,C}$ [%]")
    for axis in axes:
        axis.axhline(0.0, color="0.35", linewidth=0.8)
        axis.grid(alpha=0.25, linewidth=0.6)
        axis.set_xlabel(r"$T_w$")
        axis.set_xticks(TEMPERATURES)
    axes[0].legend(frameon=False)
    fig.suptitle("Blackburn critical-parameter errors relative to compressible model")
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_error_tecplot(path: Path, rows: list[dict[str, object]]) -> None:
    variables = (
        "Tw",
        "R_c_blackburn",
        "R_c_compressible",
        "R_error_vs_compressible_percent",
        "PPT_R_difference_compressible_vs_blackburn_percent",
        "beta_c_blackburn",
        "beta_c_compressible",
        "beta_error_vs_compressible_percent",
        "PPT_beta_difference_compressible_vs_blackburn_percent",
    )
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write('TITLE = "Blackburn critical-point errors"\n')
        stream.write(
            "VARIABLES = "
            + " ".join(f'"{variable}"' for variable in variables)
            + "\n"
        )
        for mode in ("Type-I", "Type-II"):
            available = [
                row
                for row in rows
                if row["mode"] == mode and row["status"] == "available"
            ]
            stream.write(
                f'ZONE T="{mode}", I={len(available)}, DATAPACKING=POINT\n'
            )
            for row in available:
                stream.write(
                    " ".join(f"{float(row[variable]):.12e}" for variable in variables)
                    + "\n"
                )


def write_report(
    path: Path,
    rows: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
    isothermal_check: dict[str, float] | None,
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Blackburn neutral-curve comparison\n\n")
        stream.write(
            "The six wall temperatures follow slides 13-14 of "
            "`E:/Zhj/Boussinesq.pptx`: "
            "`1.00, 1.04, 1.08, 1.12, 1.16, 1.20`.\n\n"
        )
        stream.write(
            "Critical points are interior minima of `R(beta)`, refined by a "
            "three-point quadratic fit. `beta_c >= 0.055` is classified as "
            "Type-I; smaller `beta_c` is Type-II.\n\n"
        )
        stream.write(
            "The manuscript-oriented errors use the fully compressible result "
            "as reference: `100*(Blackburn-compressible)/compressible`. The TSV "
            "also records the original PPT convention "
            "`100*(compressible-Blackburn)/Blackburn`.\n\n"
        )
        stream.write("## Critical points and errors\n\n")
        stream.write(
            "| Tw | Mode | Blackburn Rc | Compressible Rc | Rc error | "
            "Blackburn beta | Compressible beta | beta error | Status |\n"
        )
        stream.write("|---:|---|---:|---:|---:|---:|---:|---:|---|\n")
        for row in rows:
            if row["status"] != "available":
                stream.write(
                    f"| {row['Tw']:.2f} | {row['mode']} | "
                    f"{row.get('R_c_blackburn', '')} |  |  | "
                    f"{row.get('beta_c_blackburn', '')} |  |  | "
                    f"{row['status']} |\n"
                )
                continue
            stream.write(
                f"| {row['Tw']:.2f} | {row['mode']} | "
                f"{row['R_c_blackburn']:.6f} | "
                f"{row['R_c_compressible']:.6f} | "
                f"{row['R_error_vs_compressible_percent']:.3f}% | "
                f"{row['beta_c_blackburn']:.8f} | "
                f"{row['beta_c_compressible']:.8f} | "
                f"{row['beta_error_vs_compressible_percent']:.3f}% | "
                "available |\n"
            )
        stream.write("\n## Curve validation\n\n")
        stream.write("| Tw | Model | Segment | Points | beta range | R range | max residual |\n")
        stream.write("|---:|---|---:|---:|---:|---:|---:|\n")
        for row in validation_rows:
            stream.write(
                f"| {row['Tw']:.2f} | {row['model']} | {row['segment']} | "
                f"{row['points']} | {row['beta_min']:.6f}-"
                f"{row['beta_max']:.6f} | {row['R_min']:.3f}-"
                f"{row['R_max']:.3f} | {row['max_residual']:.3e} |\n"
            )
        if isothermal_check is not None:
            stream.write("\n## Isothermal implementation check\n\n")
            stream.write(
                "At `Tw=1`, Blackburn must reduce exactly to Lopez. On the "
                "common beta interval, interpolation gives "
                f"`max |R_Blackburn-R_Lopez| = "
                f"{isothermal_check['max_abs_R_difference']:.3e}` and "
                f"`RMS difference = {isothermal_check['rms_R_difference']:.3e}`.\n"
            )


def main() -> None:
    workspace = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--blackburn-dir",
        type=Path,
        default=workspace / "blackburn_neutral_curve_batch",
    )
    parser.add_argument(
        "--compressible-dir",
        type=Path,
        default=workspace / "neutral_curve_batch",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "blackburn_neutral_results",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cases: dict[tuple[float, str], list[np.ndarray]] = {}
    critical: dict[tuple[float, str], dict[str, CriticalPoint]] = {}
    validation_rows: list[dict[str, object]] = []

    for Tw in TEMPERATURES:
        paths_by_model = {
            "blackburn": blackburn_paths(args.blackburn_dir, Tw),
            "compressible": [compressible_path(args.compressible_dir, Tw)],
        }
        for model, paths in paths_by_model.items():
            segments = [load_curve(path) for path in paths]
            cases[(Tw, model)] = segments
            critical[(Tw, model)] = critical_points(segments)
            for segment_index, (path, data) in enumerate(
                zip(paths, segments), start=1
            ):
                metrics = validate_segment(data, path)
                validation_rows.append(
                    {
                        "Tw": Tw,
                        "model": model,
                        "segment": segment_index,
                        "source_file": path.name,
                        **metrics,
                    }
                )

    rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for mode in ("Type-I", "Type-II"):
            blackburn = critical[(Tw, "blackburn")].get(mode)
            compressible = critical[(Tw, "compressible")].get(mode)
            row: dict[str, object] = {"Tw": Tw, "mode": mode}
            if blackburn is None:
                row["status"] = "blackburn_critical_missing"
                rows.append(row)
                continue
            row.update(
                {
                    "R_c_blackburn": blackburn.R,
                    "beta_c_blackburn": blackburn.beta,
                    "R_c_blackburn_discrete": blackburn.R_discrete,
                    "beta_c_blackburn_discrete": blackburn.beta_discrete,
                    "blackburn_quadratic_fit": blackburn.fitted,
                }
            )
            if compressible is None:
                row["status"] = "compressible_critical_missing"
                rows.append(row)
                continue
            delta_R = blackburn.R - compressible.R
            delta_beta = blackburn.beta - compressible.beta
            row.update(
                {
                    "status": "available",
                    "R_c_compressible": compressible.R,
                    "beta_c_compressible": compressible.beta,
                    "R_c_compressible_discrete": compressible.R_discrete,
                    "beta_c_compressible_discrete": compressible.beta_discrete,
                    "compressible_quadratic_fit": compressible.fitted,
                    "delta_R_blackburn_minus_compressible": delta_R,
                    "R_error_vs_compressible_percent": 100
                    * delta_R
                    / compressible.R,
                    "PPT_R_difference_compressible_vs_blackburn_percent": 100
                    * (compressible.R - blackburn.R)
                    / blackburn.R,
                    "delta_beta_blackburn_minus_compressible": delta_beta,
                    "beta_error_vs_compressible_percent": 100
                    * delta_beta
                    / compressible.beta,
                    "PPT_beta_difference_compressible_vs_blackburn_percent": 100
                    * (compressible.beta - blackburn.beta)
                    / blackburn.beta,
                }
            )
            rows.append(row)

    columns = [
        "Tw",
        "mode",
        "status",
        "R_c_blackburn",
        "beta_c_blackburn",
        "R_c_blackburn_discrete",
        "beta_c_blackburn_discrete",
        "blackburn_quadratic_fit",
        "R_c_compressible",
        "beta_c_compressible",
        "R_c_compressible_discrete",
        "beta_c_compressible_discrete",
        "compressible_quadratic_fit",
        "delta_R_blackburn_minus_compressible",
        "R_error_vs_compressible_percent",
        "PPT_R_difference_compressible_vs_blackburn_percent",
        "delta_beta_blackburn_minus_compressible",
        "beta_error_vs_compressible_percent",
        "PPT_beta_difference_compressible_vs_blackburn_percent",
    ]
    tsv_path = args.output_dir / "blackburn_vs_compressible_critical_errors.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    validation_path = args.output_dir / "neutral_curve_validation.tsv"
    with validation_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "Tw",
                "model",
                "segment",
                "source_file",
                "points",
                "R_min",
                "R_max",
                "beta_min",
                "beta_max",
                "max_residual",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(validation_rows)

    write_combined_tecplot(
        args.output_dir / "blackburn_vs_compressible_neutral_curves.dat",
        cases,
    )
    plot_curves(
        args.output_dir / "blackburn_vs_compressible_neutral_curves.png",
        cases,
        critical,
    )
    plot_errors(
        args.output_dir / "blackburn_vs_compressible_critical_errors.png",
        rows,
    )
    write_error_tecplot(
        args.output_dir / "blackburn_vs_compressible_critical_errors.dat",
        rows,
    )

    isothermal_check: dict[str, float] | None = None
    lopez_path = args.compressible_dir / "ome=0.0_Tw=1.0_model=lopez.dat"
    if lopez_path.exists():
        lopez = load_curve(lopez_path)
        blackburn = cases[(1.0, "blackburn")][0]
        lower = max(np.min(lopez[:, 2]), np.min(blackburn[:, 2]))
        upper = min(np.max(lopez[:, 2]), np.max(blackburn[:, 2]))
        sample_beta = blackburn[
            (blackburn[:, 2] >= lower) & (blackburn[:, 2] <= upper), 2
        ]
        order = np.argsort(lopez[:, 2])
        lopez_R = np.interp(sample_beta, lopez[order, 2], lopez[order, 1])
        blackburn_order = np.argsort(blackburn[:, 2])
        blackburn_R = np.interp(
            sample_beta,
            blackburn[blackburn_order, 2],
            blackburn[blackburn_order, 1],
        )
        difference = blackburn_R - lopez_R
        isothermal_check = {
            "max_abs_R_difference": float(np.max(np.abs(difference))),
            "rms_R_difference": float(np.sqrt(np.mean(difference**2))),
        }

    report_path = args.output_dir / "blackburn_vs_compressible_summary.md"
    write_report(report_path, rows, validation_rows, isothermal_check)

    print(f"Critical-error table: {tsv_path}")
    print(f"Validation table: {validation_path}")
    print(f"Summary: {report_path}")


if __name__ == "__main__":
    main()
