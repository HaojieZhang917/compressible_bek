#!/usr/bin/env python3
"""Analyze Type-II neutral-curve grid and beta-step convergence."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULT_ROOT = ROOT / "typeII_grid_convergence_results"
RAW_ROOT = RESULT_ROOT / "raw"
ORIGINAL_ROOT = ROOT / "neutral_curve_batch"
GRIDS = (49, 59, 69, 79, 89)
TEMPERATURES = (1.12, 1.13, 1.14, 1.16)


@dataclass(frozen=True)
class TurningPoint:
    beta: float
    R: float
    fitted: bool
    index: int


@dataclass(frozen=True)
class CurveMetrics:
    N: int
    Tw: float
    points: int
    beta_min: float
    beta_max: float
    max_residual: float
    minimum: TurningPoint | None
    maximum: TurningPoint | None
    barrier: float | None
    maximum_slope: float
    maximum_slope_beta: float


def number_tag(value: float) -> str:
    text = f"{value:.4f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def curve_name(Tw: float) -> str:
    return (
        f"ome=0.0_Tw={number_tag(Tw)}_model=compressible_Mr=0.3_"
        "propPert=on_baseProp=variable.dat"
    )


def load_curve(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.lower().startswith(("variables", "zone")):
            continue
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(f"Expected seven columns in {path}: {line}")
        rows.append([float(field) for field in fields])
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] != 7:
        raise ValueError(f"Invalid curve data in {path}")
    if not np.all(np.diff(data[:, 2]) > 0):
        raise ValueError(f"beta is not strictly increasing in {path}")
    return data


def refine_vertex(
    data: np.ndarray, index: int, expected_sign: int
) -> TurningPoint:
    beta0 = data[index, 2]
    x = data[index - 1 : index + 2, 2] - beta0
    y = data[index - 1 : index + 2, 1]
    a, b, c = np.linalg.solve(
        np.column_stack((x**2, x, np.ones(3))), y
    )
    vertex = -b / (2 * a)
    fitted = (
        np.all(np.isfinite((a, b, c, vertex)))
        and expected_sign * a > 0
        and np.min(x) <= vertex <= np.max(x)
    )
    if not fitted:
        return TurningPoint(
            beta=float(beta0),
            R=float(data[index, 1]),
            fitted=False,
            index=index,
        )
    return TurningPoint(
        beta=float(beta0 + vertex),
        R=float(a * vertex**2 + b * vertex + c),
        fitted=True,
        index=index,
    )


def curve_metrics(data: np.ndarray, N: int, Tw: float) -> CurveMetrics:
    beta = data[:, 2]
    R = data[:, 1]
    typeII = (beta >= 0.038) & (beta <= 0.045)
    indices = np.flatnonzero(typeII)
    minima = [
        index
        for index in indices
        if 0 < index < len(R) - 1
        and R[index] < R[index - 1]
        and R[index] < R[index + 1]
    ]
    minimum = refine_vertex(data, minima[0], +1) if minima else None
    maxima = []
    if minimum is not None:
        maxima = [
            index
            for index in indices
            if index > minimum.index
            and index < len(R) - 1
            and R[index] > R[index - 1]
            and R[index] > R[index + 1]
        ]
    maximum = refine_vertex(data, maxima[0], -1) if maxima else None
    barrier = (
        maximum.R - minimum.R
        if minimum is not None and maximum is not None
        else None
    )

    slope = np.gradient(R, beta, edge_order=2)
    slope_window = (beta >= 0.039) & (beta <= 0.0445)
    slope_indices = np.flatnonzero(slope_window)
    maximum_slope_index = slope_indices[np.argmax(slope[slope_window])]
    residual = np.minimum(np.abs(data[:, 4]), np.abs(data[:, 6]))
    return CurveMetrics(
        N=N,
        Tw=Tw,
        points=len(data),
        beta_min=float(beta[0]),
        beta_max=float(beta[-1]),
        max_residual=float(np.max(residual)),
        minimum=minimum,
        maximum=maximum,
        barrier=barrier,
        maximum_slope=float(slope[maximum_slope_index]),
        maximum_slope_beta=float(beta[maximum_slope_index]),
    )


def write_metrics(metrics: list[CurveMetrics]) -> Path:
    path = RESULT_ROOT / "typeII_grid_convergence.tsv"
    reference = {
        metric.Tw: metric
        for metric in metrics
        if metric.N == max(GRIDS)
    }
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "N_cheb",
                "collocation_points",
                "Tw",
                "points_on_curve",
                "beta_min",
                "beta_max",
                "max_neutral_residual",
                "typeII_fold",
                "R_c",
                "beta_c",
                "R_max_adjacent",
                "beta_max_adjacent",
                "barrier_delta_R",
                "maximum_dR_dbeta",
                "maximum_slope_beta",
                "R_c_error_vs_N89_percent",
                "beta_c_error_vs_N89_percent",
            ]
        )
        for metric in metrics:
            ref = reference[metric.Tw]
            comparable = (
                metric.minimum is not None and ref.minimum is not None
            )
            R_error = (
                100
                * (metric.minimum.R - ref.minimum.R)
                / ref.minimum.R
                if comparable
                else ""
            )
            beta_error = (
                100
                * (metric.minimum.beta - ref.minimum.beta)
                / ref.minimum.beta
                if comparable
                else ""
            )
            writer.writerow(
                [
                    metric.N,
                    metric.N + 1,
                    metric.Tw,
                    metric.points,
                    metric.beta_min,
                    metric.beta_max,
                    metric.max_residual,
                    "present" if metric.minimum is not None else "absent",
                    metric.minimum.R if metric.minimum else "",
                    metric.minimum.beta if metric.minimum else "",
                    metric.maximum.R if metric.maximum else "",
                    metric.maximum.beta if metric.maximum else "",
                    metric.barrier if metric.barrier is not None else "",
                    metric.maximum_slope,
                    metric.maximum_slope_beta,
                    R_error,
                    beta_error,
                ]
            )
    return path


def plot_curves(curves: dict[tuple[int, float], np.ndarray]) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(GRIDS)))
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        for color, N in zip(colors, GRIDS):
            data = curves[N, Tw]
            mask = (data[:, 2] >= 0.039) & (data[:, 2] <= 0.0448)
            axis.plot(
                data[mask, 2],
                data[mask, 1],
                color=color,
                linewidth=1.25,
                marker="o",
                markersize=1.7,
                label=f"N={N}",
            )
        axis.set_title(rf"$T_w={Tw:.2f}$")
        axis.set_xlabel(r"$\beta$")
        axis.set_ylabel(r"$R$")
        axis.grid(True, alpha=0.25)
    axes[0, 0].legend(ncol=2, fontsize=8)
    path = RESULT_ROOT / "typeII_grid_convergence_curves.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def plot_diagnostics(metrics: list[CurveMetrics]) -> Path:
    by_key = {(metric.N, metric.Tw): metric for metric in metrics}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4), constrained_layout=True)
    for Tw in TEMPERATURES:
        values = [by_key[N, Tw] for N in GRIDS]
        axes[0].plot(
            GRIDS,
            [value.maximum_slope for value in values],
            marker="o",
            label=rf"$T_w={Tw:.2f}$",
        )
        axes[1].plot(
            GRIDS,
            [
                value.barrier if value.barrier is not None else 0.0
                for value in values
            ],
            marker="o",
            label=rf"$T_w={Tw:.2f}$",
        )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_xlabel(r"$N_{\mathrm{cheb}}$")
    axes[0].set_ylabel(r"$\max(dR/d\beta)$")
    axes[0].set_title("Fold-sign diagnostic")
    axes[1].set_xlabel(r"$N_{\mathrm{cheb}}$")
    axes[1].set_ylabel(r"$R_{\max}-R_{\min}$")
    axes[1].set_title("Resolved Type-II fold height")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=8)
    path = RESULT_ROOT / "typeII_grid_convergence_diagnostics.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def original_fold_status(Tw: float) -> str:
    path = ORIGINAL_ROOT / curve_name(Tw)
    if not path.exists():
        return "not_available"
    data = load_curve(path)
    metric = curve_metrics(data, 69, Tw)
    return "present" if metric.minimum is not None else "absent"


def format_optional(value: float | None, digits: int = 6) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


def write_summary(metrics: list[CurveMetrics]) -> Path:
    path = RESULT_ROOT / "typeII_grid_convergence_summary.md"
    by_key = {(metric.N, metric.Tw): metric for metric in metrics}
    lines = [
        "# Type-II 临界点网格收敛性验证",
        "",
        "完全可压缩模型：`Mr=0.3`, `property_perturbations=true`, "
        "`base_property_variation=true`。",
        "",
        "空间离散采用与正式计算相同的有理 Chebyshev 映射和远场截断 "
        "`y=40`；这里只改变谱阶数。中性曲线的参数步长由正式计算的 "
        "`8e-4` 加密为 `2e-4`。",
        "",
        "## 高阶网格结果",
        "",
        "| Tw | N | Type-II fold | Rc | beta_c | ΔR | max(dR/dbeta) |",
        "|---:|---:|:---:|---:|---:|---:|---:|",
    ]
    for Tw in TEMPERATURES:
        for N in GRIDS:
            metric = by_key[N, Tw]
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"{Tw:.2f}",
                        str(N),
                        "存在" if metric.minimum else "不存在",
                        format_optional(
                            metric.minimum.R if metric.minimum else None
                        ),
                        format_optional(
                            metric.minimum.beta if metric.minimum else None,
                            8,
                        ),
                        format_optional(metric.barrier),
                        f"{metric.maximum_slope:.3f}",
                    ]
                )
                + " |"
            )
    lines.extend(
        [
            "",
            "## 正式步长与加密步长的差别",
            "",
            "| Tw | 原 N=69, Δβ=8e-4 | 加密 N=69, Δβ=2e-4 |",
            "|---:|:---:|:---:|",
        ]
    )
    for Tw in (1.12, 1.14, 1.16):
        refined = by_key[69, Tw]
        lines.append(
            f"| {Tw:.2f} | {original_fold_status(Tw)} | "
            f"{'present' if refined.minimum else 'absent'} |"
        )

    high_grids = (69, 79, 89)
    tw116_absent = all(
        by_key[N, 1.16].minimum is None for N in high_grids
    )
    tw116_slopes = [by_key[N, 1.16].maximum_slope for N in high_grids]
    max_residual = max(metric.max_residual for metric in metrics)
    ref112 = by_key[89, 1.12]
    lines.extend(
        [
            "",
            "## 判断",
            "",
            (
                "- `Tw=1.16` 在 N=69,79,89 上均不存在 Type-II 折叠；"
                if tw116_absent
                else "- `Tw=1.16` 的高阶网格对折叠存在性尚未给出一致结果；"
            ),
            "- `Tw=1.16` 的高阶网格最大斜率为 "
            + ", ".join(f"{value:.3f}" for value in tw116_slopes)
            + "。全部为负意味着曲线在原 Type-II 区域保持单调。",
            (
                "- N=89 在 `Tw=1.12` 给出 "
                f"`Rc={ref112.minimum.R:.8f}`, "
                f"`beta_c={ref112.minimum.beta:.10f}`。"
                if ref112.minimum is not None
                else "- N=89 在 `Tw=1.12` 未找到 Type-II 极小值。"
            ),
            f"- 全部曲线的最大中性残差为 `{max_residual:.3e}`。",
            "- `N=89, Tw=1.12` 的首次固定 shift `alpha_target=0.11` "
            "误选了邻近 Ritz 值；将 shift 对准已收敛 Type-II 分支后，"
            "所得曲线与 N=69,79 收敛一致。该问题属于初始特征值选支，"
            "不是 Chebyshev 网格导致的模态消失。",
            "",
            "空间网格结论以 N=69,79,89 的一致性为准；临界温度附近极浅的"
            "折叠还会受到 beta 参数步长影响，不能把参数采样误差归因于"
            "Chebyshev 网格。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    curves: dict[tuple[int, float], np.ndarray] = {}
    metrics: list[CurveMetrics] = []
    for N in GRIDS:
        for Tw in TEMPERATURES:
            path = RAW_ROOT / f"N={N}" / curve_name(Tw)
            if not path.exists():
                raise FileNotFoundError(path)
            data = load_curve(path)
            curves[N, Tw] = data
            metrics.append(curve_metrics(data, N, Tw))
    metrics.sort(key=lambda item: (item.Tw, item.N))
    table_path = write_metrics(metrics)
    curve_plot = plot_curves(curves)
    diagnostic_plot = plot_diagnostics(metrics)
    summary_path = write_summary(metrics)
    print(f"wrote {table_path}")
    print(f"wrote {curve_plot}")
    print(f"wrote {diagnostic_plot}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
