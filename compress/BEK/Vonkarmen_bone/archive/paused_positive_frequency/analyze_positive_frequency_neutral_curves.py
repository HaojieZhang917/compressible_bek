#!/usr/bin/env python3
"""Analyze omega>0 Type-II neutral curves for Blackburn and compressible models."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20)


@dataclass(frozen=True)
class CriticalPoint:
    R: float
    beta: float
    R_discrete: float
    beta_discrete: float
    fitted: bool


def number_tag(value: float) -> str:
    text = f"{value:.4f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def curve_path(directory: Path, model: str, omega: float, Tw: float) -> Path:
    stem = (
        f"ome={number_tag(omega)}_Tw={number_tag(Tw)}_model={model}"
    )
    if model == "compressible":
        stem += "_Mr=0.3_propPert=on_baseProp=variable"
    return directory / f"{stem}.dat"


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


def typeII_segment(path: Path) -> tuple[np.ndarray, float | None]:
    data = load_curve(path)
    log_path = path.with_suffix(".log")
    switch_beta: float | None = None
    if log_path.exists():
        match = re.search(
            r"mode_switch beta=([-+0-9.eE]+)",
            log_path.read_text(encoding="utf-8", errors="replace"),
        )
        if match is not None:
            switch_beta = float(match.group(1))
            data = data[data[:, 2] < switch_beta]
    if data.shape[0] < 3:
        raise ValueError(f"Fewer than three Type-II points in {path}")
    return data, switch_beta


def validate_curve(data: np.ndarray, path: Path) -> dict[str, float | int]:
    if data.shape[1] != 7 or not np.all(np.isfinite(data)):
        raise ValueError(f"Invalid curve data in {path}")
    beta_difference = np.diff(data[:, 2])
    if not np.all(beta_difference > 0):
        raise ValueError(f"beta is not strictly increasing in {path}")
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


def critical_point(data: np.ndarray) -> CriticalPoint:
    index = int(np.argmin(data[:, 1]))
    if index == 0 or index == data.shape[0] - 1:
        return CriticalPoint(
            R=float(data[index, 1]),
            beta=float(data[index, 2]),
            R_discrete=float(data[index, 1]),
            beta_discrete=float(data[index, 2]),
            fitted=False,
        )
    beta0 = data[index, 2]
    x = data[index - 1 : index + 2, 2] - beta0
    y = data[index - 1 : index + 2, 1]
    a, b, c = np.linalg.solve(
        np.column_stack((x**2, x, np.ones(3))), y
    )
    fitted = bool(
        np.all(np.isfinite((a, b, c)))
        and a > 0
        and np.min(x) <= -b / (2 * a) <= np.max(x)
    )
    if fitted:
        vertex = -b / (2 * a)
        R = a * vertex**2 + b * vertex + c
        beta = beta0 + vertex
    else:
        R = data[index, 1]
        beta = beta0
    return CriticalPoint(
        R=float(R),
        beta=float(beta),
        R_discrete=float(data[index, 1]),
        beta_discrete=float(beta0),
        fitted=fitted,
    )


def curve_error(blackburn: np.ndarray, compressible: np.ndarray) -> dict[str, float]:
    lower = max(float(np.min(blackburn[:, 2])), float(np.min(compressible[:, 2])))
    upper = min(float(np.max(blackburn[:, 2])), float(np.max(compressible[:, 2])))
    if upper <= lower:
        return {
            "common_beta_min": float("nan"),
            "common_beta_max": float("nan"),
            "curve_R_error_mean_percent": float("nan"),
            "curve_R_error_rms_percent": float("nan"),
            "curve_R_error_max_abs_percent": float("nan"),
        }
    beta = np.linspace(lower, upper, 401)
    R_blackburn = np.interp(beta, blackburn[:, 2], blackburn[:, 1])
    R_compressible = np.interp(beta, compressible[:, 2], compressible[:, 1])
    relative = 100.0 * (R_blackburn - R_compressible) / R_compressible
    return {
        "common_beta_min": lower,
        "common_beta_max": upper,
        "curve_R_error_mean_percent": float(np.mean(relative)),
        "curve_R_error_rms_percent": float(np.sqrt(np.mean(relative**2))),
        "curve_R_error_max_abs_percent": float(np.max(np.abs(relative))),
    }


def plot_curves(
    path: Path,
    cases: dict[tuple[float, str], np.ndarray],
    critical: dict[tuple[float, str], CriticalPoint],
    omega: float,
) -> None:
    figure, axes = plt.subplots(
        2, 3, figsize=(13.0, 7.2), sharex=True, sharey=False
    )
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        for model, color, label in (
            ("compressible", "black", "Fully compressible"),
            ("blackburn", "#d62728", "Blackburn"),
        ):
            data = cases[(Tw, model)]
            point = critical[(Tw, model)]
            axis.plot(data[:, 1], data[:, 2], color=color, linewidth=1.7, label=label)
            axis.scatter(
                point.R,
                point.beta,
                s=28,
                marker="s",
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=4,
            )
        axis.text(0.04, 0.94, rf"$T_w={Tw:.2f}$", transform=axis.transAxes, va="top")
        axis.grid(alpha=0.22, linewidth=0.6)
        beta_values = np.concatenate(
            (
                cases[(Tw, "blackburn")][:, 2],
                cases[(Tw, "compressible")][:, 2],
            )
        )
        beta_span = float(np.max(beta_values) - np.min(beta_values))
        beta_padding = max(0.12 * beta_span, 5.0e-4)
        axis.set_ylim(
            float(np.min(beta_values) - beta_padding),
            float(np.max(beta_values) + beta_padding),
        )
    axes[0, 0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.82))
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$R$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\beta$")
    figure.tight_layout()
    figure.subplots_adjust(top=0.90)
    figure.suptitle(rf"Type-II neutral lobe at $\omega={omega:.3f}$", y=0.975)
    figure.savefig(path, dpi=240)
    plt.close(figure)


def plot_errors(path: Path, rows: list[dict[str, float]], omega: float) -> None:
    Tw = [row["Tw"] for row in rows]
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), sharex=True)
    series = (
        ("R_error_percent", r"$\Delta R_c/R_{c,C}$ [%]"),
        ("beta_error_percent", r"$\Delta\beta_c/\beta_{c,C}$ [%]"),
        ("curve_R_error_rms_percent", r"curve RMS $R$ error [%]"),
    )
    for axis, (key, ylabel) in zip(axes, series):
        axis.plot(Tw, [row[key] for row in rows], "o-", color="#d62728")
        axis.axhline(0.0, color="0.35", linewidth=0.8)
        axis.axhline(5.0, color="0.55", linewidth=0.7, linestyle="--")
        axis.axhline(-5.0, color="0.55", linewidth=0.7, linestyle="--")
        axis.set_xlabel(r"$T_w$")
        axis.set_ylabel(ylabel)
        axis.set_xticks(TEMPERATURES)
        axis.grid(alpha=0.23, linewidth=0.6)
    axes[2].text(
        0.97,
        0.08,
        r"No common $\beta$ interval for $T_w\geq1.12$",
        transform=axes[2].transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="0.35",
    )
    figure.suptitle(
        rf"Blackburn errors relative to fully compressible model, $\omega={omega:.3f}$"
    )
    figure.tight_layout()
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def threshold_summary(rows: list[dict[str, float]], threshold: float) -> str:
    acceptable = [
        row
        for row in rows
        if abs(row["R_error_percent"]) <= threshold
        and abs(row["beta_error_percent"]) <= threshold
        and row["curve_R_error_rms_percent"] <= threshold
    ]
    if not acceptable:
        return "none"
    return ", ".join(f"{row['Tw']:.2f}" for row in acceptable)


def read_omega_zero_reference(path: Path) -> dict[tuple[float, str], dict[str, str]]:
    if not path.exists():
        return {}
    reference: dict[tuple[float, str], dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row.get("mode") == "Type-II":
                reference[(float(row["Tw"]), "comparison")] = row
    return reference


def frequency_shift_rows(
    rows: list[dict[str, float]],
    omega_zero: dict[tuple[float, str], dict[str, str]],
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for row in rows:
        old = omega_zero.get((row["Tw"], "comparison"))
        if old is None:
            continue
        old_blackburn = float(old["R_c_blackburn"])
        shifted: dict[str, object] = {
            "Tw": row["Tw"],
            "R_c_blackburn_omega0": old_blackburn,
            "R_c_blackburn_omega_positive": row["R_c_blackburn"],
            "blackburn_Rc_change_percent": 100.0
            * (row["R_c_blackburn"] - old_blackburn)
            / old_blackburn,
            "compressible_omega0_status": old["status"],
            "R_c_compressible_omega_positive": row["R_c_compressible"],
        }
        if old["status"] == "available":
            old_compressible = float(old["R_c_compressible"])
            shifted.update(
                {
                    "R_c_compressible_omega0": old_compressible,
                    "compressible_Rc_change_percent": 100.0
                    * (row["R_c_compressible"] - old_compressible)
                    / old_compressible,
                }
            )
        result.append(shifted)
    return result


def write_summary(
    path: Path,
    rows: list[dict[str, float]],
    omega: float,
    omega_zero: dict[tuple[float, str], dict[str, str]],
    frequency_rows: list[dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"# Positive-frequency Type-II comparison (omega={omega:.3f})\n\n")
        stream.write(
            "The selected frequency is Lingwood's modest positive-frequency "
            "case. Curves are stopped immediately before the Type-I/Type-II "
            "neutral branch point, so all reported minima belong to Type II.\n\n"
        )
        stream.write(
            "Errors use the fully compressible model as reference: "
            "`100*(Blackburn-compressible)/compressible`.\n\n"
        )
        stream.write(
            "| Tw | Blackburn Rc | Compressible Rc | Rc error | Blackburn beta_c | "
            "Compressible beta_c | beta error | curve RMS R error |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            curve_error_text = (
                f"{row['curve_R_error_rms_percent']:.3f}%"
                if np.isfinite(row["curve_R_error_rms_percent"])
                else "no common beta interval"
            )
            stream.write(
                f"| {row['Tw']:.2f} | {row['R_c_blackburn']:.6f} | "
                f"{row['R_c_compressible']:.6f} | {row['R_error_percent']:.3f}% | "
                f"{row['beta_c_blackburn']:.8f} | "
                f"{row['beta_c_compressible']:.8f} | "
                f"{row['beta_error_percent']:.3f}% | "
                f"{curve_error_text} |\n"
            )
        stream.write(
            "\nQuadratic fits are used for smooth interior minima. At "
            "`Tw>=1.12`, the Type-II minimum is the neutral-branch endpoint, "
            "so its last resolved point is reported without extrapolation.\n"
        )
        stream.write("\n## Operational applicability bands\n\n")
        stream.write(
            "A temperature passes a band only if the absolute errors in both "
            "critical parameters and the curve-wide RMS R error are within it.\n\n"
        )
        stream.write(f"- 2% band: {threshold_summary(rows, 2.0)}\n")
        stream.write(f"- 5% band: {threshold_summary(rows, 5.0)}\n")
        if omega_zero and frequency_rows:
            stream.write("\n## Comparison with omega=0 results\n\n")
            stream.write(
                "| Tw | Blackburn Rc change | Compressible Rc change | "
                "omega=0 compressible status | omega>0 status |\n"
            )
            stream.write("|---:|---:|---:|---|---|\n")
            for row in frequency_rows:
                compressible_change = (
                    f"{row['compressible_Rc_change_percent']:.3f}%"
                    if "compressible_Rc_change_percent" in row
                    else "not available"
                )
                stream.write(
                    f"| {row['Tw']:.2f} | "
                    f"{row['blackburn_Rc_change_percent']:.3f}% | "
                    f"{compressible_change} | "
                    f"{row['compressible_omega0_status']} | available |\n"
                )


def main() -> None:
    workspace = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=float, default=0.008)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=workspace / "positive_frequency_neutral_curve_batch" / "omega=0.008",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=workspace / "positive_frequency_neutral_results" / "omega=0.008",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cases: dict[tuple[float, str], np.ndarray] = {}
    critical: dict[tuple[float, str], CriticalPoint] = {}
    validation_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for model in ("blackburn", "compressible"):
            path = curve_path(args.input_dir, model, args.omega, Tw)
            data, switch_beta = typeII_segment(path)
            metrics = validate_curve(data, path)
            cases[(Tw, model)] = data
            critical[(Tw, model)] = critical_point(data)
            validation_rows.append(
                {
                    "Tw": Tw,
                    "model": model,
                    "source_file": path.name,
                    "first_mode_switch_beta": switch_beta,
                    **metrics,
                }
            )

    rows: list[dict[str, float]] = []
    for Tw in TEMPERATURES:
        blackburn = critical[(Tw, "blackburn")]
        compressible = critical[(Tw, "compressible")]
        errors = curve_error(
            cases[(Tw, "blackburn")], cases[(Tw, "compressible")]
        )
        rows.append(
            {
                "Tw": Tw,
                "R_c_blackburn": blackburn.R,
                "R_c_compressible": compressible.R,
                "R_error_percent": 100.0
                * (blackburn.R - compressible.R)
                / compressible.R,
                "beta_c_blackburn": blackburn.beta,
                "beta_c_compressible": compressible.beta,
                "beta_error_percent": 100.0
                * (blackburn.beta - compressible.beta)
                / compressible.beta,
                "R_c_blackburn_discrete": blackburn.R_discrete,
                "R_c_compressible_discrete": compressible.R_discrete,
                "beta_c_blackburn_discrete": blackburn.beta_discrete,
                "beta_c_compressible_discrete": compressible.beta_discrete,
                "blackburn_quadratic_fit": blackburn.fitted,
                "compressible_quadratic_fit": compressible.fitted,
                **errors,
            }
        )

    error_path = args.output_dir / "positive_frequency_typeII_errors.tsv"
    with error_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    validation_path = args.output_dir / "positive_frequency_curve_validation.tsv"
    with validation_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(validation_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(validation_rows)

    plot_curves(
        args.output_dir / "positive_frequency_neutral_curves.png",
        cases,
        critical,
        args.omega,
    )
    plot_errors(
        args.output_dir / "positive_frequency_typeII_errors.png",
        rows,
        args.omega,
    )
    omega_zero = read_omega_zero_reference(
        workspace
        / "blackburn_neutral_results"
        / "blackburn_vs_compressible_critical_errors.tsv"
    )
    frequency_rows = frequency_shift_rows(rows, omega_zero)
    frequency_path = (
        args.output_dir / "positive_vs_zero_frequency_typeII_shift.tsv"
    )
    if frequency_rows:
        frequency_columns = sorted(
            {key for row in frequency_rows for key in row},
            key=lambda key: (
                key != "Tw",
                key,
            ),
        )
        with frequency_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=frequency_columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(frequency_rows)
    summary_path = args.output_dir / "positive_frequency_summary.md"
    write_summary(
        summary_path, rows, args.omega, omega_zero, frequency_rows
    )
    print(f"Critical errors: {error_path}")
    print(f"Validation: {validation_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
