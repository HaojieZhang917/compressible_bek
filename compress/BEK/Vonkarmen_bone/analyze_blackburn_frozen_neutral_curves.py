#!/usr/bin/env python3
"""Compare Blackburn with two frozen-property compressible operators."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_blackburn_neutral_curves import (
    CriticalPoint,
    blackburn_paths,
    critical_points,
    load_curve,
    number_tag,
    validate_segment,
)


TEMPERATURES = (1.04, 1.08, 1.12, 1.16, 1.20)
REFERENCE_CASES = (
    "perturbations_frozen",
    "all_properties_frozen",
)
CASE_LABELS = {
    "blackburn": "Blackburn",
    "perturbations_frozen": "Compressible: perturbation properties frozen",
    "all_properties_frozen": "Compressible: all transport properties frozen",
}


def compressible_path(directory: Path, Tw: float, case: str) -> Path:
    base = (
        f"ome=0.0_Tw={number_tag(Tw)}_model=compressible_Mr=0.3_"
        "propPert=off_baseProp="
    )
    if case == "perturbations_frozen":
        return directory / f"{base}variable.dat"
    if case == "all_properties_frozen":
        return directory / f"{base}frozen.dat"
    raise ValueError(f"Unknown frozen-property case: {case}")


def error_row(
    Tw: float,
    mode: str,
    reference_case: str,
    blackburn: CriticalPoint,
    reference: CriticalPoint,
) -> dict[str, object]:
    delta_R = blackburn.R - reference.R
    delta_beta = blackburn.beta - reference.beta
    return {
        "Tw": Tw,
        "mode": mode,
        "reference_case": reference_case,
        "R_c_blackburn": blackburn.R,
        "beta_c_blackburn": blackburn.beta,
        "R_c_reference": reference.R,
        "beta_c_reference": reference.beta,
        "delta_R_blackburn_minus_reference": delta_R,
        "R_error_vs_reference_percent": 100 * delta_R / reference.R,
        "PPT_R_difference_reference_vs_blackburn_percent": (
            100 * (reference.R - blackburn.R) / blackburn.R
        ),
        "delta_beta_blackburn_minus_reference": delta_beta,
        "beta_error_vs_reference_percent": 100 * delta_beta / reference.beta,
        "PPT_beta_difference_reference_vs_blackburn_percent": (
            100 * (reference.beta - blackburn.beta) / blackburn.beta
        ),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = (
        "Tw",
        "mode",
        "reference_case",
        "R_c_blackburn",
        "beta_c_blackburn",
        "R_c_reference",
        "beta_c_reference",
        "delta_R_blackburn_minus_reference",
        "R_error_vs_reference_percent",
        "PPT_R_difference_reference_vs_blackburn_percent",
        "delta_beta_blackburn_minus_reference",
        "beta_error_vs_reference_percent",
        "PPT_beta_difference_reference_vs_blackburn_percent",
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_tecplot_errors(path: Path, rows: list[dict[str, object]]) -> None:
    variables = (
        "Tw",
        "R_c_blackburn",
        "R_c_reference",
        "R_error_vs_reference_percent",
        "beta_c_blackburn",
        "beta_c_reference",
        "beta_error_vs_reference_percent",
    )
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write('TITLE = "Blackburn versus frozen-property errors"\n')
        stream.write(
            "VARIABLES = "
            + " ".join(f'"{variable}"' for variable in variables)
            + "\n"
        )
        for reference_case in REFERENCE_CASES:
            for mode in ("Type-I", "Type-II"):
                group = [
                    row
                    for row in rows
                    if row["reference_case"] == reference_case
                    and row["mode"] == mode
                ]
                stream.write(
                    f'ZONE T="{reference_case} {mode}", '
                    f"I={len(group)}, DATAPACKING=POINT\n"
                )
                for row in group:
                    stream.write(
                        " ".join(
                            f"{float(row[variable]):.12e}"
                            for variable in variables
                        )
                        + "\n"
                    )


def write_tecplot_curves(
    path: Path,
    cases: dict[tuple[float, str], list[np.ndarray]],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(
            'TITLE = "Blackburn and frozen-property compressible neutral curves"\n'
        )
        stream.write(
            'VARIABLES = "omega" "R" "beta" "alpha_r_1" '
            '"alpha_i_1" "alpha_r_2" "alpha_i_2"\n'
        )
        for Tw in TEMPERATURES:
            for case in ("blackburn", *REFERENCE_CASES):
                for segment_index, data in enumerate(
                    cases[(Tw, case)], start=1
                ):
                    stream.write(
                        f'ZONE T="Tw={Tw:.2f} {case} segment={segment_index}", '
                        f"I={data.shape[0]}, DATAPACKING=POINT\n"
                    )
                    np.savetxt(stream, data, fmt="%.12e")


def plot_curves(
    path: Path,
    cases: dict[tuple[float, str], list[np.ndarray]],
    critical: dict[tuple[float, str], dict[str, CriticalPoint]],
) -> None:
    styles = {
        "blackburn": ("#d62728", "-", 1.8),
        "perturbations_frozen": ("#1f77b4", "--", 1.5),
        "all_properties_frozen": ("black", "-", 1.5),
    }
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), sharex=True, sharey=True)
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        for case in ("all_properties_frozen", "perturbations_frozen", "blackburn"):
            color, linestyle, linewidth = styles[case]
            for segment_index, data in enumerate(cases[(Tw, case)]):
                axis.plot(
                    data[:, 1],
                    data[:, 2],
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    label=CASE_LABELS[case] if segment_index == 0 else None,
                )
            for mode, marker in (("Type-I", "o"), ("Type-II", "s")):
                point = critical[(Tw, case)].get(mode)
                if point is not None:
                    axis.scatter(
                        point.R,
                        point.beta,
                        s=20,
                        marker=marker,
                        facecolor=color,
                        edgecolor="white",
                        linewidth=0.4,
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
        axis.set_xlim(200, 525)
        axis.set_ylim(0.025, 0.13)
        axis.grid(alpha=0.22, linewidth=0.6)
    legend_axis = axes.flat[-1]
    legend_axis.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    legend_axis.legend(handles, labels, loc="center", frameon=False, fontsize=11)
    legend_axis.text(
        0.5,
        0.25,
        "Markers: circle = Type-I, square = Type-II",
        ha="center",
        transform=legend_axis.transAxes,
        fontsize=10,
    )
    for axis in axes[-1, :2]:
        axis.set_xlabel(r"$R$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\beta$")
    fig.suptitle(
        "Neutral curves: Blackburn vs frozen-property compressible operators",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_errors(path: Path, rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), sharex=True)
    colors = {"Type-I": "#1f77b4", "Type-II": "#d62728"}
    markers = {"Type-I": "o", "Type-II": "s"}
    linestyles = {
        "perturbations_frozen": "--",
        "all_properties_frozen": "-",
    }
    for reference_case in REFERENCE_CASES:
        for mode in ("Type-I", "Type-II"):
            group = [
                row
                for row in rows
                if row["reference_case"] == reference_case
                and row["mode"] == mode
            ]
            label = (
                f"{mode}, "
                + (
                    "perturbation properties frozen"
                    if reference_case == "perturbations_frozen"
                    else "all transport properties frozen"
                )
            )
            axes[0].plot(
                [float(row["Tw"]) for row in group],
                [float(row["R_error_vs_reference_percent"]) for row in group],
                color=colors[mode],
                marker=markers[mode],
                linestyle=linestyles[reference_case],
                label=label,
            )
            axes[1].plot(
                [float(row["Tw"]) for row in group],
                [
                    float(row["beta_error_vs_reference_percent"])
                    for row in group
                ],
                color=colors[mode],
                marker=markers[mode],
                linestyle=linestyles[reference_case],
                label=label,
            )
    axes[0].set_ylabel(r"$100(R_{c,B}-R_{c,ref})/R_{c,ref}$ [%]")
    axes[1].set_ylabel(
        r"$100(\beta_{c,B}-\beta_{c,ref})/\beta_{c,ref}$ [%]"
    )
    for axis in axes:
        axis.axhline(0.0, color="0.35", linewidth=0.8)
        axis.grid(alpha=0.25, linewidth=0.6)
        axis.set_xlabel(r"$T_w$")
        axis.set_xticks(TEMPERATURES)
    axes[0].legend(frameon=False, fontsize=8.5)
    fig.suptitle("Blackburn errors relative to frozen-property operators")
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def format_comparison_table(
    rows: list[dict[str, object]], reference_case: str
) -> str:
    lines = [
        "| Tw | Mode | reference Rc | Blackburn Rc | Rc difference | "
        "reference beta | Blackburn beta | beta difference |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        if row["reference_case"] != reference_case:
            continue
        lines.append(
            f"| {row['Tw']:.2f} | {row['mode']} | "
            f"{row['R_c_reference']:.6f} | {row['R_c_blackburn']:.6f} | "
            f"{row['R_error_vs_reference_percent']:+.3f}% | "
            f"{row['beta_c_reference']:.8f} | "
            f"{row['beta_c_blackburn']:.8f} | "
            f"{row['beta_error_vs_reference_percent']:+.3f}% |"
        )
    return "\n".join(lines)


def write_report(
    path: Path,
    rows: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Blackburn 与冻结物性可压缩算子的中性曲线比较\n\n")
        stream.write(
            "共同温度范围为 `Tw=1.04, 1.08, 1.12, 1.16, 1.20`。"
            "`neutral_curve_batch` 中没有 `Tw=1.00` 的冻结物性结果。\n\n"
        )
        stream.write(
            "误差定义为 `100*(Blackburn-reference)/reference`。"
            "所有临界点均为 `R(beta)` 的内部极小值，并采用三点二次拟合。\n\n"
        )
        stream.write(
            "这里的物性开关只修改稳定性算子：`property_perturbations=false` "
            "关闭 Chapman 律的扰动物性项；`base_property_variation=false` "
            "进一步把算子中的平均黏度、第二黏度和热导率冻结在远场值。"
            "既有可压缩基本流不会被重新求解，密度、压力和能量耦合仍然保留。\n\n"
        )
        stream.write("## 扰动物性与平均输运物性全部冻结\n\n")
        stream.write(
            format_comparison_table(rows, "all_properties_frozen") + "\n\n"
        )
        stream.write("## 仅冻结扰动物性\n\n")
        stream.write(
            format_comparison_table(rows, "perturbations_frozen") + "\n\n"
        )
        stream.write("## 验证\n\n")
        stream.write(
            f"- 共验证 {len(validation_rows)} 条模型/温度曲线；"
            f"最大中性残差为 "
            f"`{max(float(row['max_residual']) for row in validation_rows):.3e}`。\n"
        )
        stream.write(
            "- 五个温度下，Blackburn 和完全冻结可压缩算子均存在 "
            "Type-I 与 Type-II 内部极小值。\n"
        )


def main() -> None:
    workspace = Path(__file__).resolve().parent
    blackburn_dir = workspace / "blackburn_neutral_curve_batch"
    compressible_dir = workspace / "neutral_curve_batch"
    output_dir = workspace / "blackburn_frozen_neutral_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    cases: dict[tuple[float, str], list[np.ndarray]] = {}
    critical: dict[tuple[float, str], dict[str, CriticalPoint]] = {}
    validation_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        paths_by_case = {
            "blackburn": blackburn_paths(blackburn_dir, Tw),
            **{
                case: [compressible_path(compressible_dir, Tw, case)]
                for case in REFERENCE_CASES
            },
        }
        for case, paths in paths_by_case.items():
            segments = [load_curve(path) for path in paths]
            cases[(Tw, case)] = segments
            critical[(Tw, case)] = critical_points(segments)
            for segment_index, (path, data) in enumerate(
                zip(paths, segments), start=1
            ):
                validation_rows.append(
                    {
                        "Tw": Tw,
                        "case": case,
                        "segment": segment_index,
                        "source_file": path.name,
                        **validate_segment(data, path),
                    }
                )

    rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for mode in ("Type-I", "Type-II"):
            blackburn = critical[(Tw, "blackburn")][mode]
            for reference_case in REFERENCE_CASES:
                reference = critical[(Tw, reference_case)][mode]
                rows.append(
                    error_row(
                        Tw, mode, reference_case, blackburn, reference
                    )
                )

    write_tsv(
        output_dir / "blackburn_vs_frozen_critical_errors.tsv",
        rows,
    )
    write_tecplot_errors(
        output_dir / "blackburn_vs_frozen_critical_errors.dat",
        rows,
    )
    write_tecplot_curves(
        output_dir / "blackburn_vs_frozen_neutral_curves.dat",
        cases,
    )
    plot_curves(
        output_dir / "blackburn_vs_frozen_neutral_curves.png",
        cases,
        critical,
    )
    plot_errors(
        output_dir / "blackburn_vs_frozen_critical_errors.png",
        rows,
    )

    validation_path = output_dir / "neutral_curve_validation.tsv"
    with validation_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "Tw",
                "case",
                "segment",
                "source_file",
                "points",
                "R_min",
                "R_max",
                "beta_min",
                "beta_max",
                "max_residual",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(validation_rows)
    write_report(
        output_dir / "blackburn_vs_frozen_summary.md",
        rows,
        validation_rows,
    )
    print(f"Results: {output_dir}")


if __name__ == "__main__":
    main()
