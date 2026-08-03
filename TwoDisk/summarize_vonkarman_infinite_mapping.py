#!/usr/bin/env python3
"""Create compact tables and figures for the mapped von Karman branch.

The script is deliberately post-processing only: it reads the original
finite-domain branch and the rational-Chebyshev calculations, then writes
derived summaries beside the mapped results.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline


HERE = Path(__file__).resolve().parent
TOOLKIT = HERE.parent
RESULT = HERE / "vonkarman_infinite_mapping"
MAPPED = RESULT / "convergence_full"
STABILITY = RESULT / "stability_N110_L8"
STRICT_LOW = RESULT / "convergence_strict_second_fold"
STRICT_HIGH = RESULT / "convergence_tight_high_order"
FINITE = (
    TOOLKIT
    / "compress"
    / "BEK"
    / "Vonkarmen_bone"
    / "boussinesq_domain_branches"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fold_summary(mapped_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    result = []
    descriptions = {
        1: ("confirmed", "principal-branch termination"),
        2: ("confirmed", "second saddle-node on returned branch"),
        3: ("candidate", "near-Hinf=0 singular-tail region"),
    }
    for index in (1, 2, 3):
        rows = [row for row in mapped_rows if int(row["fold_index"]) == index]
        h = np.asarray([float(row["Hinf"]) for row in rows])
        tw = np.asarray([float(row["Tw"]) for row in rows])
        status, meaning = descriptions[index]
        result.append(
            {
                "fold_index": index,
                "kind": rows[0]["kind"],
                "status": status,
                "Hinf_median": float(np.median(h)),
                "Hinf_min": float(h.min()),
                "Hinf_max": float(h.max()),
                "Tw_median": float(np.median(tw)),
                "Tw_min": float(tw.min()),
                "Tw_max": float(tw.max()),
                "trusted_grid_count": len(rows),
                "interpretation": meaning,
            }
        )
    return result


def strict_first_two_folds() -> list[dict[str, str]]:
    """Extract folds from well-conditioned strict-residual branch runs."""
    output: list[dict[str, str]] = []
    for directory in (STRICT_LOW, STRICT_HIGH):
        for path in directory.glob("branch_N*_L*.csv"):
            match = re.fullmatch(r"branch_N(\d+)_L([0-9p]+)\.csv", path.name)
            if match is None:
                continue
            degree = int(match.group(1))
            scale = float(match.group(2).replace("p", "."))
            if degree not in (50, 60, 70) or scale not in (8.0, 12.0):
                continue
            rows = read_csv(path)
            h = np.asarray([float(row["Hinf"]) for row in rows])
            tw = np.asarray([float(row["Tw"]) for row in rows])
            spline = CubicSpline(h, tw)
            roots = spline.derivative().roots()
            roots = roots[(roots > h[2]) & (roots < h[-3])]
            for index, root in enumerate(roots[:2], start=1):
                second = float(spline.derivative(2)(root))
                output.append(
                    {
                        "degree": str(degree),
                        "mapping_scale": str(scale),
                        "fold_index": str(index),
                        "kind": "minimum" if second > 0.0 else "maximum",
                        "Hinf": str(float(root)),
                        "Tw": str(float(spline(root))),
                        "d2Tw_dHinf2": str(second),
                    }
                )
    return output


def stability_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if int(row["temporal_degree"]) == 90:
            grouped[(row["sample"], float(row["Hinf"]))].append(row)

    output = []
    for (sample, h_inf), group in sorted(grouped.items(), key=lambda item: item[0][1]):
        group.sort(key=lambda row: int(row["rank_by_real_part"]))
        first = group[0]
        entry: dict[str, object] = {
            "sample": sample,
            "Hinf": h_inf,
            "Tw": float(first["Tw"]),
            "temporal_degree": 90,
            "positive_real_count": int(first["positive_real_count"]),
            "robust_positive_count_above_1e-5": sum(
                float(row["eigen_real"]) > 1.0e-5 for row in group
            ),
            "closest_to_zero_real": float(first["closest_to_zero_real"]),
            "closest_to_zero_imag": float(first["closest_to_zero_imag"]),
            "fold_neutral_mode_within_1e-5": (
                sample.startswith("fold_")
                and abs(float(first["closest_to_zero_real"])) < 1.0e-5
                and abs(float(first["closest_to_zero_imag"])) < 1.0e-5
            ),
        }
        for rank in range(3):
            entry[f"lambda_{rank + 1}_real"] = float(group[rank]["eigen_real"])
            entry[f"lambda_{rank + 1}_imag"] = float(group[rank]["eigen_imag"])
        output.append(entry)
    return output


def main() -> None:
    strict_folds = strict_first_two_folds()
    loose_candidate = [
        row
        for row in read_csv(MAPPED / "fold_convergence.csv")
        if int(row["degree"]) <= 110 and int(row["fold_index"]) == 3
    ]
    mapped_folds = strict_folds + loose_candidate
    folds = fold_summary(mapped_folds)
    write_csv(RESULT / "confirmed_fold_summary.csv", folds)

    temporal = read_csv(STABILITY / "mapped_temporal_samples.csv")
    spectra = stability_summary(temporal)
    write_csv(RESULT / "upper_branch_stability_summary.csv", spectra)

    finite = read_csv(FINITE / "turning_points_by_zmax.csv")
    branch = read_csv(MAPPED / "branch_N110_L8.csv")
    profiles = read_csv(MAPPED / "profiles_N110_L8.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))

    ax = axes[0, 0]
    for index, marker, label in ((1, "o", "finite-domain fold 1"), (2, "s", "finite-domain fold 2")):
        rows = [row for row in finite if int(row["turning_point"]) == index]
        zmax = np.asarray([float(row["zmax"]) for row in rows])
        tw = np.asarray([float(row["Tw"]) for row in rows])
        ax.plot(zmax, tw, marker=marker, linewidth=1.5, label=label)
    for index, color in ((1, "C0"), (2, "C1")):
        row = folds[index - 1]
        ax.axhline(
            float(row["Tw_median"]),
            color=color,
            linestyle="--",
            alpha=0.8,
            label=f"mapped fold {index}",
        )
    ax.set_xlabel(r"finite truncation $\eta_{\max}$")
    ax.set_ylabel(r"turning-point $T_w$")
    ax.set_title("Finite truncation versus true-infinity mapping")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    h = np.asarray([float(row["Hinf"]) for row in branch])
    tw = np.asarray([float(row["Tw"]) for row in branch])
    ax.plot(h, tw, color="0.25", linewidth=1.8)
    for row, color in zip(folds, ("C3", "C0", "C2")):
        ax.scatter(row["Hinf_median"], row["Tw_median"], s=85, marker="*", color=color, zorder=3)
        suffix = "" if row["status"] == "confirmed" else " (candidate)"
        ax.annotate(
            f"F{row['fold_index']}{suffix}",
            (float(row["Hinf_median"]), float(row["Tw_median"])),
            xytext=(5, 7),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel(r"$H_\infty$")
    ax.set_ylabel(r"$T_w$")
    ax.set_title("Mapped semi-infinite branch")
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    selected = (-0.10, -0.08, -0.04, -0.02)
    profile_groups: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in profiles:
        profile_groups[float(row["Hinf"])].append(row)
    for target in selected:
        key = min(profile_groups, key=lambda value: abs(value - target))
        rows = profile_groups[key]
        eta = np.asarray([float(row["eta"]) for row in rows])
        temp = np.asarray([float(row["T"]) for row in rows])
        tw0 = float(rows[0]["Tw"])
        departure = np.abs((temp - 1.0) / (tw0 - 1.0))
        mask = (eta <= 300.0) & (departure > 1.0e-9)
        tail = 1.0 / (0.72 * abs(key))
        ax.semilogy(eta[mask], departure[mask], linewidth=1.6, label=fr"$H_\infty={key:g}$, $\ell_T={tail:.1f}$")
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$|(T-1)/(T_w-1)|$")
    ax.set_title(r"Thermal tail diverges as $H_\infty$ approaches $0^-$")
    ax.set_ylim(1.0e-8, 1.4)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    sh = np.asarray([float(row["Hinf"]) for row in spectra])
    order = np.argsort(sh)
    sh = sh[order]
    for rank, marker in ((1, "o"), (2, "s"), (3, "^")):
        eig = np.asarray([float(row[f"lambda_{rank}_real"]) for row in spectra])[order]
        ax.plot(sh, eig, marker=marker, linewidth=1.25, label=fr"$\Re(\lambda_{rank})$")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel(r"$H_\infty$")
    ax.set_ylabel(r"temporal growth rate $\Re(\lambda)$")
    ax.set_title("Similarity-subspace temporal spectrum, N=90")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(RESULT / "upper_branch_infinite_mapping_summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
