#!/usr/bin/env python3
"""Analyze complete omega=0 Blackburn/compressible curves from Tw=1 to 1.8."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20, 1.40, 1.60, 1.80)
FIXED_BETA_SPLIT = 0.055
MANUAL_ENDPOINT_REASONS = {
    "ome=0.0_Tw=1.0_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.04_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.08_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.12_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.16_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.2_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.4_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.6_model=blackburn.dat": "endpoint_no_root",
    "ome=0.0_Tw=1.8_model=blackburn.dat": "R_limit",
    "ome=0.0_Tw=1.0_model=compressible_Mr=0.3_propPert=on_baseProp=variable.dat":
        "R_limit",
    "ome=0.0_Tw=1.4_model=compressible_Mr=0.3_propPert=on_baseProp=variable.dat":
        "near_degenerate_endpoint",
    "ome=0.0_Tw=1.6_model=compressible_Mr=0.3_propPert=on_baseProp=variable.dat":
        "R_limit",
    "ome=0.0_Tw=1.8_model=compressible_Mr=0.3_propPert=on_baseProp=variable.dat":
        "R_limit",
}


@dataclass(frozen=True)
class CriticalPoint:
    mode: str
    R: float
    beta: float
    fitted: bool


def number_tag(value: float) -> str:
    text = f"{value:.4f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def curve_path(workspace: Path, model: str, Tw: float) -> Path:
    tag = number_tag(Tw)
    if model == "blackburn":
        return workspace / "blackburn_neutral_curve_batch" / (
            f"ome=0.0_Tw={tag}_model=blackburn.dat"
        )
    return workspace / "neutral_curve_batch" / (
        f"ome=0.0_Tw={tag}_model=compressible_Mr=0.3_"
        "propPert=on_baseProp=variable.dat"
    )


def load_curve(path: Path) -> np.ndarray:
    data = np.loadtxt(path, skiprows=2)
    if data.ndim != 2 or data.shape[1] != 7:
        raise ValueError(f"Expected seven-column curve in {path}")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite curve data in {path}")
    if not np.all(np.diff(data[:, 2]) > 0):
        raise ValueError(f"beta is not strictly increasing in {path}")
    residual = np.minimum(np.abs(data[:, 4]), np.abs(data[:, 6]))
    if float(np.max(residual)) > 1.0e-6:
        raise ValueError(f"Neutral residual exceeds 1e-6 in {path}")
    return data


def stop_reason(path: Path) -> str:
    trace_path = path.with_suffix(".log")
    if trace_path.exists():
        for line in reversed(trace_path.read_text(encoding="utf-8").splitlines()):
            match = re.search(r"finished reason=([A-Za-z0-9_]+)", line)
            if match:
                return match.group(1)
    status_path = path.parent / "batch_status.tsv"
    if status_path.exists():
        case = path.stem
        for line in reversed(status_path.read_text(encoding="utf-8").splitlines()):
            fields = line.split("\t", 5)
            if len(fields) < 6 or fields[1] != case or fields[2] != "ok":
                continue
            match = re.search(r"stop_reason=([A-Za-z0-9_]+)", fields[5])
            if match:
                return match.group(1)
    return MANUAL_ENDPOINT_REASONS.get(path.name, "unknown")


def extrema_indices(data: np.ndarray) -> tuple[list[int], list[int]]:
    minima: list[int] = []
    maxima: list[int] = []
    for index in range(1, data.shape[0] - 1):
        left, center, right = data[index - 1 : index + 2, 1]
        if center < left and center < right:
            minima.append(index)
        if center > left and center > right:
            maxima.append(index)
    return minima, maxima


def refine_minimum(data: np.ndarray, index: int) -> tuple[float, float, bool]:
    beta0 = data[index, 2]
    x = data[index - 1 : index + 2, 2] - beta0
    y = data[index - 1 : index + 2, 1]
    a, b, c = np.linalg.solve(np.column_stack((x**2, x, np.ones(3))), y)
    if not np.all(np.isfinite((a, b, c))) or a <= 0:
        return float(data[index, 1]), float(beta0), False
    vertex = -b / (2.0 * a)
    if vertex < float(np.min(x)) or vertex > float(np.max(x)):
        return float(data[index, 1]), float(beta0), False
    return float(a * vertex**2 + b * vertex + c), float(beta0 + vertex), True


def critical_points(data: np.ndarray) -> dict[str, CriticalPoint]:
    """Classify modes by extrema topology, not by a fixed beta threshold."""
    minima, maxima = extrema_indices(data)
    candidates = [refine_minimum(data, index) for index in minima]
    candidates.sort(key=lambda point: point[1])
    if not candidates:
        return {}
    if len(candidates) >= 2:
        low, high = candidates[0], candidates[-1]
        return {
            "Type-II": CriticalPoint("Type-II", *low),
            "Type-I": CriticalPoint("Type-I", *high),
        }
    only = candidates[0]
    maxima_beta = [float(data[index, 2]) for index in maxima]
    mode = "Type-II" if any(beta > only[1] for beta in maxima_beta) else "Type-I"
    return {mode: CriticalPoint(mode, *only)}


def sampled_error(
    blackburn: np.ndarray,
    compressible: np.ndarray,
    lower: float,
    upper: float,
    samples: int = 1601,
    R_ceiling: float | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    keys = (
        "mean_signed_percent",
        "mean_abs_percent",
        "rms_percent",
        "p95_abs_percent",
        "max_abs_percent",
        "beta_at_max_abs",
        "normalized_L1_percent",
        "positive_fraction",
    )
    if upper <= lower:
        return np.array([]), np.array([]), {key: float("nan") for key in keys}
    beta = np.linspace(lower, upper, samples)
    Rb = np.interp(beta, blackburn[:, 2], blackburn[:, 1])
    Rc = np.interp(beta, compressible[:, 2], compressible[:, 1])
    if R_ceiling is not None:
        keep = (Rb <= R_ceiling) & (Rc <= R_ceiling)
        beta, Rb, Rc = beta[keep], Rb[keep], Rc[keep]
    if beta.size < 3:
        return np.array([]), np.array([]), {key: float("nan") for key in keys}
    relative = 100.0 * (Rb - Rc) / Rc
    absolute = np.abs(relative)
    maximum_index = int(np.argmax(absolute))
    metrics = {
        "mean_signed_percent": float(np.mean(relative)),
        "mean_abs_percent": float(np.mean(absolute)),
        "rms_percent": float(np.sqrt(np.mean(relative**2))),
        "p95_abs_percent": float(np.quantile(absolute, 0.95)),
        "max_abs_percent": float(absolute[maximum_index]),
        "beta_at_max_abs": float(beta[maximum_index]),
        "normalized_L1_percent": float(
            100.0
        * np.trapz(np.abs(Rb - Rc), beta)
            / np.trapz(Rc, beta)
        ),
        "positive_fraction": float(np.mean(relative > 0.0)),
    }
    return beta, relative, metrics


def curve_metrics(
    blackburn: np.ndarray, compressible: np.ndarray
) -> tuple[dict[str, float], tuple[np.ndarray, np.ndarray]]:
    common_min = max(float(blackburn[0, 2]), float(compressible[0, 2]))
    common_max = min(float(blackburn[-1, 2]), float(compressible[-1, 2]))
    union_min = min(float(blackburn[0, 2]), float(compressible[0, 2]))
    union_max = max(float(blackburn[-1, 2]), float(compressible[-1, 2]))
    regions = {
        "full": (common_min, common_max, None),
        "lower_beta": (common_min, min(common_max, FIXED_BETA_SPLIT), None),
        "upper_beta": (max(common_min, FIXED_BETA_SPLIT), common_max, None),
        "core_R_le_450": (common_min, common_max, 450.0),
    }
    output: dict[str, float] = {
        "common_beta_min": common_min,
        "common_beta_max": common_max,
        "common_beta_fraction_percent": 100.0
        * (common_max - common_min)
        / (union_max - union_min),
        "beta_min_error_percent": 100.0
        * (float(blackburn[0, 2]) - float(compressible[0, 2]))
        / float(compressible[0, 2]),
        "beta_max_error_percent": 100.0
        * (float(blackburn[-1, 2]) - float(compressible[-1, 2]))
        / float(compressible[-1, 2]),
    }
    full_profile = (np.array([]), np.array([]))
    for name, (lower, upper, ceiling) in regions.items():
        beta, relative, metrics = sampled_error(
            blackburn, compressible, lower, upper, R_ceiling=ceiling
        )
        if name == "full":
            full_profile = (beta, relative)
        for metric, value in metrics.items():
            output[f"{name}_{metric}"] = value
    return output, full_profile


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_neutral_curves(
    path: Path, curves: dict[tuple[float, str], np.ndarray]
) -> None:
    figure, axes = plt.subplots(
        3, 3, figsize=(12.6, 10.2), sharex=True, sharey=True
    )
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        for model, label, color, style in (
            ("blackburn", "Blackburn", "#1f77b4", "-"),
            ("compressible", "Compressible", "#d62728", "--"),
        ):
            data = curves[(Tw, model)]
            axis.plot(data[:, 1], data[:, 2], style, color=color, lw=1.6, label=label)
        axis.text(0.04, 0.93, rf"$T_w={Tw:.2f}$", transform=axis.transAxes, va="top")
        axis.set_xlabel(r"$R$")
        axis.grid(alpha=0.22, lw=0.6)
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\beta$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.972),
        ncol=2, frameon=False,
    )
    figure.suptitle(r"Zero-frequency neutral branches", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.935))
    figure.savefig(path, dpi=240)
    plt.close(figure)


def plot_error_profiles(
    path: Path, profiles: dict[float, tuple[np.ndarray, np.ndarray]]
) -> None:
    figure, axes = plt.subplots(3, 3, figsize=(12.6, 10.2), sharey=True)
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        beta, relative = profiles[Tw]
        axis.plot(beta, relative, color="#b22222", lw=1.6)
        axis.axhline(0.0, color="0.3", lw=0.8)
        axis.axvline(FIXED_BETA_SPLIT, color="0.5", lw=0.7, ls="--")
        axis.text(0.04, 0.93, rf"$T_w={Tw:.2f}$", transform=axis.transAxes, va="top")
        axis.set_xlabel(r"$\beta$")
        axis.grid(alpha=0.22, lw=0.6)
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$100(R_B-R_C)/R_C$ [\%]")
    figure.suptitle(r"Pointwise neutral-curve errors at $\omega=0$", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.975))
    figure.savefig(path, dpi=240)
    plt.close(figure)


def plot_summary(
    path: Path,
    critical_rows: list[dict[str, object]],
    curve_rows: list[dict[str, object]],
    density_rows: list[dict[str, object]],
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.0), sharex=True)
    available = [
        row for row in critical_rows
        if row["mode"] == "Type-I" and row["status"] == "available"
    ]
    axes[0, 0].plot(
        [row["Tw"] for row in available],
        [row["R_error_percent"] for row in available],
        "o-", label=r"$R_c$ error",
    )
    axes[0, 0].plot(
        [row["Tw"] for row in available],
        [row["beta_error_percent"] for row in available],
        "s--", label=r"$\beta_c$ error",
    )
    axes[0, 0].set_ylabel("Type-I critical error [%]")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].plot(
        [row["Tw"] for row in curve_rows],
        [row["full_rms_percent"] for row in curve_rows],
        "o-", label="Full curve",
    )
    axes[0, 1].plot(
        [row["Tw"] for row in curve_rows],
        [row["core_R_le_450_rms_percent"] for row in curve_rows],
        "s--", label=r"Core ($R\leq450$)",
    )
    axes[0, 1].set_ylabel("Pointwise RMS error [%]")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].plot(
        [row["Tw"] for row in curve_rows],
        [row["common_beta_fraction_percent"] for row in curve_rows],
        "o-", color="#2ca02c",
    )
    axes[1, 0].set_ylabel(r"Common $\beta$ coverage [%]")
    axes[1, 0].set_ylim(0, 105)

    axes[1, 1].plot(
        [row["Tw"] for row in density_rows],
        [row["linear_density_relative_error_percent"] for row in density_rows],
        "o-", color="#9467bd",
    )
    axes[1, 1].set_ylabel(r"Wall density linearization error [%]")

    for axis in axes.flat:
        axis.axhline(0.0, color="0.35", lw=0.8)
        if axis in axes[0, :]:
            axis.axhline(5.0, color="0.55", lw=0.7, ls="--")
            axis.axhline(-5.0, color="0.55", lw=0.7, ls="--")
        axis.set_xlabel(r"$T_w$")
        axis.set_xticks(TEMPERATURES)
        axis.tick_params(axis="x", rotation=35)
        axis.grid(alpha=0.22, lw=0.6)
    figure.suptitle("Extended zero-frequency model comparison", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.965))
    figure.savefig(path, dpi=240)
    plt.close(figure)


def last_discrete_within(rows: list[dict[str, object]], key: str, limit: float) -> float | None:
    values = [float(row["Tw"]) for row in rows if abs(float(row[key])) <= limit]
    return max(values) if values else None


def format_optional(value: object, digits: int = 3) -> str:
    if value == "" or value is None:
        return "—"
    return f"{float(value):.{digits}f}"


def write_report(
    path: Path,
    critical_rows: list[dict[str, object]],
    curve_rows: list[dict[str, object]],
    topology_rows: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
    density_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
) -> None:
    typeI = [row for row in critical_rows if row["mode"] == "Type-I"]
    typeI_available = [row for row in typeI if row["status"] == "available"]
    typeII = [row for row in critical_rows if row["mode"] == "Type-II"]
    common_typeII = [row for row in typeII if row["status"] == "available"]
    last_typeII_common = max(float(row["Tw"]) for row in common_typeII)
    first_topology_mismatch = min(
        float(row["Tw"]) for row in typeII if row["status"] != "available"
    )
    full_5 = last_discrete_within(curve_rows, "full_rms_percent", 5.0)
    core_5 = last_discrete_within(curve_rows, "core_R_le_450_rms_percent", 5.0)
    typeI_R_5 = last_discrete_within(typeI_available, "R_error_percent", 5.0)
    typeI_beta_5 = last_discrete_within(typeI_available, "beta_error_percent", 5.0)
    last_curve = curve_rows[-1]
    last_typeI = typeI_available[-1]
    max_residual = max(float(row["max_neutral_residual"]) for row in validation_rows)
    min_overlap = min(float(row["common_beta_fraction_percent"]) for row in curve_rows)
    max_overlap = max(float(row["common_beta_fraction_percent"]) for row in curve_rows)

    lines = [
        "# $T_w=1.0$--$1.8$ 零频中性分支综合误差分析",
        "",
        "## 1. 数据与方法",
        "",
        "比较 Blackburn 广义 Boussinesq 模型与完全可压缩模型（$M_r=0.3$、扰动物性开启、基流物性可变）的 $\\omega=0$ 中性曲线。温度点为 $1.00,1.04,1.08,1.12,1.16,1.20,1.40,1.60,1.80$。误差统一定义为 $100(Q_B-Q_C)/Q_C$。",
        "",
        "临界模态不再用固定 $\\beta$ 阈值分类，而按局部极值拓扑跟踪：若有两个局部极小值，低 $\\beta$ 为 Type-II、高 $\\beta$ 为 Type-I；若只剩一个且其后没有局部极大值，则识别为保留下来的 Type-I。这样可避免高温下 Type-I 向低 $\\beta$ 移动后被误标为 Type-II。",
        "",
        f"全部 {len(validation_rows)} 条分支均通过七列格式、有限值、$\\beta$ 单调性和中性残差检查；最大中性残差为 {max_residual:.3e}。两模型共同 $\\beta$ 区间覆盖并集的 {min_overlap:.1f}%--{max_overlap:.1f}%，所有逐点误差均只在共同区间内插值，没有外推。$T_w=1.4$ 的可压缩低-$\\beta$ 端在近简并点以 `near_degenerate_endpoint` 标记，没有将其外推到 $R=500$。",
        "",
        "## 2. 临界参数",
        "",
        "| $T_w$ | Type-I $R_{c,B}$ | Type-I $R_{c,C}$ | $R_c$ 误差 | $\\beta_c$ 误差 | Type-II 状态/误差 |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for Tw in TEMPERATURES:
        rowI = next(row for row in typeI if row["Tw"] == Tw)
        rowII = next(row for row in typeII if row["Tw"] == Tw)
        if rowII["status"] == "available":
            second = (
                f"$R_c$ {float(rowII['R_error_percent']):+.3f}%, "
                f"$\\beta_c$ {float(rowII['beta_error_percent']):+.3f}%"
            )
        else:
            second = str(rowII["status"])
        lines.append(
            f"| {Tw:.2f} | {format_optional(rowI.get('R_c_blackburn'))} | "
            f"{format_optional(rowI.get('R_c_compressible'))} | "
            f"{format_optional(rowI.get('R_error_percent'))}% | "
            f"{format_optional(rowI.get('beta_error_percent'))}% | {second} |"
        )
    lines.extend(
        [
            "",
            f"在 5% 离散判据下，Type-I $R_c$ 可支持到 $T_w={typeI_R_5:.2f}$；Type-I $\\beta_c$ 可支持到 $T_w={typeI_beta_5:.2f}$。在最高温 $T_w={float(last_typeI['Tw']):.2f}$，两者误差分别为 {float(last_typeI['R_error_percent']):+.3f}% 和 {float(last_typeI['beta_error_percent']):+.3f}%。",
            "",
            f"Type-II 只有到 $T_w={last_typeII_common:.2f}$ 才能在两模型之间定义有限的临界参数误差；从 $T_w={first_topology_mismatch:.2f}$ 起至少一方缺少对应折点，必须报告为拓扑失配而不是百分比误差。",
            "",
            "## 3. 整条曲线误差",
            "",
            "| $T_w$ | 全曲线 RMS | 核心区 RMS ($R\\le450$) | 全曲线 MAE | 95% 误差 | 最大误差 | 共同 $\\beta$ 覆盖 |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in curve_rows:
        lines.append(
            f"| {float(row['Tw']):.2f} | {float(row['full_rms_percent']):.3f}% | "
            f"{float(row['core_R_le_450_rms_percent']):.3f}% | "
            f"{float(row['full_mean_abs_percent']):.3f}% | "
            f"{float(row['full_p95_abs_percent']):.3f}% | "
            f"{float(row['full_max_abs_percent']):.3f}% | "
            f"{float(row['common_beta_fraction_percent']):.1f}% |"
        )
    lines.extend(
        [
            "",
            f"全曲线与核心区 RMS 在 5% 判据下分别只支持到 $T_w={full_5:.2f}$ 和 $T_w={core_5:.2f}$。到 $T_w={float(last_curve['Tw']):.2f}$，全曲线 RMS={float(last_curve['full_rms_percent']):.3f}%、核心区 RMS={float(last_curve['core_R_le_450_rms_percent']):.3f}%、95% 误差={float(last_curve['full_p95_abs_percent']):.3f}%。",
            "",
            "固定 $\\beta$ 的逐点误差会把曲线水平位移转化为较大的垂直 $R$ 差，尤其在接近 $R\\simeq500$ 的陡峭尾部。因此最大误差用于识别错位，适用性判断更应看核心区 RMS、临界参数和拓扑是否一致。",
            "",
            "## 4. 拓扑变化",
            "",
            "| $T_w$ | Blackburn 极小/极大值数 | Blackburn Type-II | Blackburn 终止 | 可压缩极小/极大值数 | 可压缩 Type-II | 可压缩终止 |",
            "|---:|---:|---|---|---:|---|---|",
        ]
    )
    for Tw in TEMPERATURES:
        B = next(row for row in topology_rows if row["Tw"] == Tw and row["model"] == "blackburn")
        C = next(row for row in topology_rows if row["Tw"] == Tw and row["model"] == "compressible")
        lines.append(
            f"| {Tw:.2f} | {B['local_minima_count']}/{B['local_maxima_count']} | "
            f"{B['TypeII_fold']} | {B['stop_reason']} | "
            f"{C['local_minima_count']}/{C['local_maxima_count']} | "
            f"{C['TypeII_fold']} | {C['stop_reason']} |"
        )
    lines.extend(
        [
            "",
            "Type-II 折点的有无比单个临界值误差更强：它决定模型是否预测出一整类低 $\\beta$ 中性分支。高温下若 Blackburn 仍保留而可压缩模型已消失，则广义 Boussinesq 不只是定量偏差，而是产生了额外的失稳结构。",
            "",
            "## 5. 密度线性化与高温解释",
            "",
            "| $T_w$ | Blackburn 壁面 $\\chi_w=2-T_w$ | 理想气体 $\\rho_w/\\rho_\\infty=1/T_w$ | 线性化相对误差 |",
            "|---:|---:|---:|---:|",
        ]
    )
    for row in density_rows:
        lines.append(
            f"| {float(row['Tw']):.2f} | {float(row['linear_density_ratio']):.4f} | "
            f"{float(row['ideal_gas_density_ratio']):.4f} | "
            f"{float(row['linear_density_relative_error_percent']):+.1f}% |"
        )
    lines.extend(
        [
            "",
            "该误差恰为 $-100(T_w-1)^2$%。它是 Blackburn 高温误差快速增长的明确渐近警告，但不是中性曲线误差的单一预测公式，因为两模型还同时区别于速度散度、温度扰动耦合、黏性/导热物性和能量方程。$T_w=1.8$ 时 $\\chi_w=0.2$，模型已接近 $T_w=2$ 的零线性密度极限；因此 1.4--1.8 的结果应定位为“失效诊断区”，而非 Boussinesq 定量预测区。",
            "",
            "## 6. 温度响应",
            "",
            "| 模态/量 | 共同温度区间 | Blackburn 变化 | 可压缩变化 | 变化幅值相对误差 |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in sensitivity_rows:
        lines.append(
            f"| {row['mode']} {row['quantity']} | "
            f"{float(row['Tw_start']):.2f}→{float(row['Tw_end']):.2f} | "
            f"{float(row['blackburn_change_percent']):+.3f}% | "
            f"{float(row['compressible_change_percent']):+.3f}% | "
            f"{float(row['change_magnitude_error_percent']):+.1f}% |"
        )
    lines.extend(
        [
            "",
            "若变化幅值误差在不同模态或参数之间符号不同，就不能用一个统一校正因子把 Blackburn 曲线映射到可压缩结果。论文应将这一点表述为模型选择性，而不是简单的整体高估或低估。",
            "",
            "## 7. 综合结论与论文表述",
            "",
            "1. 低温区应分别评价 Type-I 临界值、临界波数、完整曲线和 Type-II 拓扑；只看 $R_c$ 会掩盖选择性误差。",
            "2. 中温区的失效首先表现为波数位置和曲线几何偏移，随后表现为 Type-II 折点的拓扑不一致。",
            "3. 高温区中线性密度闭合本身已远离 $1/T$，Blackburn 结果主要用于展示近似如何失效；完全可压缩模型成为定量讨论所必需。",
            "4. 5% 只是报告阈值，不是控制方程给出的普适边界。所有范围结论只对应当前离散温度点，不能把相邻温度之间的真实转变位置当作已确定。",
            "",
            "## 8. 配套文件",
            "",
            "- `extended_zero_frequency_critical_errors.tsv`：拓扑跟踪后的临界参数",
            "- `extended_zero_frequency_curve_shape_errors.tsv`：全曲线、核心区和固定 $\\beta$ 分区指标",
            "- `extended_zero_frequency_topology.tsv`：全部局部极值及折点状态",
            "- `extended_zero_frequency_curve_validation.tsv`：原始曲线完整性",
            "- `extended_zero_frequency_density_indicator.tsv`：线性密度误差",
            "- `extended_zero_frequency_temperature_sensitivity.tsv`：温度响应",
            "- `extended_zero_frequency_neutral_curves.png`：九个温度的曲线叠加",
            "- `extended_zero_frequency_pointwise_errors.png`：逐点误差分布",
            "- `extended_zero_frequency_error_summary.png`：综合误差摘要",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    workspace = Path(__file__).resolve().parent
    output_dir = workspace / "zero_frequency_extended_error_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    curves: dict[tuple[float, str], np.ndarray] = {}
    critical: dict[tuple[float, str], dict[str, CriticalPoint]] = {}
    validation_rows: list[dict[str, object]] = []
    topology_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for model in ("blackburn", "compressible"):
            path = curve_path(workspace, model, Tw)
            data = load_curve(path)
            reason = stop_reason(path)
            curves[(Tw, model)] = data
            critical[(Tw, model)] = critical_points(data)
            minima, maxima = extrema_indices(data)
            residual = np.minimum(np.abs(data[:, 4]), np.abs(data[:, 6]))
            validation_rows.append(
                {
                    "Tw": Tw,
                    "model": model,
                    "source_file": path.name,
                    "points": int(data.shape[0]),
                    "beta_min": float(data[0, 2]),
                    "beta_max": float(data[-1, 2]),
                    "R_min": float(np.min(data[:, 1])),
                    "R_max": float(np.max(data[:, 1])),
                    "max_neutral_residual": float(np.max(residual)),
                    "stop_reason": reason,
                }
            )
            topology_rows.append(
                {
                    "Tw": Tw,
                    "model": model,
                    "local_minima_count": len(minima),
                    "local_maxima_count": len(maxima),
                    "local_minima_beta": ",".join(f"{data[index, 2]:.9g}" for index in minima),
                    "local_maxima_beta": ",".join(f"{data[index, 2]:.9g}" for index in maxima),
                    "TypeII_fold": "present" if "Type-II" in critical[(Tw, model)] else "absent",
                    "TypeI_minimum": "present" if "Type-I" in critical[(Tw, model)] else "absent",
                    "stop_reason": reason,
                }
            )

    critical_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for mode in ("Type-I", "Type-II"):
            B = critical[(Tw, "blackburn")].get(mode)
            C = critical[(Tw, "compressible")].get(mode)
            row: dict[str, object] = {
                "Tw": Tw,
                "mode": mode,
                "R_c_blackburn": B.R if B else "",
                "beta_c_blackburn": B.beta if B else "",
                "R_c_compressible": C.R if C else "",
                "beta_c_compressible": C.beta if C else "",
            }
            if B is not None and C is not None:
                row.update(
                    {
                        "status": "available",
                        "R_error_percent": 100.0 * (B.R - C.R) / C.R,
                        "beta_error_percent": 100.0 * (B.beta - C.beta) / C.beta,
                        "blackburn_quadratic_fit": B.fitted,
                        "compressible_quadratic_fit": C.fitted,
                    }
                )
            elif B is None and C is None:
                row["status"] = "both_missing"
            elif B is None:
                row["status"] = "blackburn_missing"
            else:
                row["status"] = "compressible_missing"
            critical_rows.append(row)

    curve_rows: list[dict[str, object]] = []
    profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for Tw in TEMPERATURES:
        metrics, profile = curve_metrics(
            curves[(Tw, "blackburn")], curves[(Tw, "compressible")]
        )
        curve_rows.append({"Tw": Tw, **metrics})
        profiles[Tw] = profile

    density_rows = [
        {
            "Tw": Tw,
            "linear_density_ratio": 2.0 - Tw,
            "ideal_gas_density_ratio": 1.0 / Tw,
            "linear_density_relative_error_percent": -100.0 * (Tw - 1.0) ** 2,
        }
        for Tw in TEMPERATURES
    ]

    sensitivity_rows: list[dict[str, object]] = []
    for mode in ("Type-I", "Type-II"):
        common = [
            Tw for Tw in TEMPERATURES
            if mode in critical[(Tw, "blackburn")]
            and mode in critical[(Tw, "compressible")]
        ]
        if len(common) < 2:
            continue
        Tw_start, Tw_end = min(common), max(common)
        for quantity, attribute in (("R_c", "R"), ("beta_c", "beta")):
            changes: dict[str, float] = {}
            values: dict[str, tuple[float, float]] = {}
            for model in ("blackburn", "compressible"):
                start = getattr(critical[(Tw_start, model)][mode], attribute)
                end = getattr(critical[(Tw_end, model)][mode], attribute)
                values[model] = (start, end)
                changes[model] = 100.0 * (end - start) / start
            reference = abs(changes["compressible"])
            magnitude_error = (
                100.0 * (abs(changes["blackburn"]) - reference) / reference
                if reference > 0 else float("nan")
            )
            sensitivity_rows.append(
                {
                    "mode": mode,
                    "quantity": quantity,
                    "Tw_start": Tw_start,
                    "Tw_end": Tw_end,
                    "blackburn_start": values["blackburn"][0],
                    "blackburn_end": values["blackburn"][1],
                    "blackburn_change_percent": changes["blackburn"],
                    "compressible_start": values["compressible"][0],
                    "compressible_end": values["compressible"][1],
                    "compressible_change_percent": changes["compressible"],
                    "change_magnitude_error_percent": magnitude_error,
                }
            )

    write_tsv(output_dir / "extended_zero_frequency_critical_errors.tsv", critical_rows)
    write_tsv(output_dir / "extended_zero_frequency_curve_shape_errors.tsv", curve_rows)
    write_tsv(output_dir / "extended_zero_frequency_topology.tsv", topology_rows)
    write_tsv(output_dir / "extended_zero_frequency_curve_validation.tsv", validation_rows)
    write_tsv(output_dir / "extended_zero_frequency_density_indicator.tsv", density_rows)
    write_tsv(
        output_dir / "extended_zero_frequency_temperature_sensitivity.tsv",
        sensitivity_rows,
    )
    plot_neutral_curves(output_dir / "extended_zero_frequency_neutral_curves.png", curves)
    plot_error_profiles(output_dir / "extended_zero_frequency_pointwise_errors.png", profiles)
    plot_summary(
        output_dir / "extended_zero_frequency_error_summary.png",
        critical_rows,
        curve_rows,
        density_rows,
    )
    write_report(
        output_dir / "extended_zero_frequency_detailed_summary.md",
        critical_rows,
        curve_rows,
        topology_rows,
        validation_rows,
        density_rows,
        sensitivity_rows,
    )
    print(f"Output: {output_dir}")


if __name__ == "__main__":
    main()
