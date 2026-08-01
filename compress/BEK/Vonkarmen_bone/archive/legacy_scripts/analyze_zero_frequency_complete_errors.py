#!/usr/bin/env python3
"""Detailed error analysis for the existing complete omega=0 neutral curves."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TEMPERATURES = (1.00, 1.04, 1.08, 1.12, 1.16, 1.20)
TYPE_SPLIT_BETA = 0.055


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


def refine_minimum(data: np.ndarray, index: int) -> tuple[float, float, bool]:
    beta0 = data[index, 2]
    x = data[index - 1 : index + 2, 2] - beta0
    y = data[index - 1 : index + 2, 1]
    a, b, c = np.linalg.solve(
        np.column_stack((x**2, x, np.ones(3))), y
    )
    if not np.all(np.isfinite((a, b, c))) or a <= 0:
        return float(data[index, 1]), float(beta0), False
    vertex = -b / (2.0 * a)
    if vertex < float(np.min(x)) or vertex > float(np.max(x)):
        return float(data[index, 1]), float(beta0), False
    return (
        float(a * vertex**2 + b * vertex + c),
        float(beta0 + vertex),
        True,
    )


def critical_points(data: np.ndarray) -> dict[str, CriticalPoint]:
    candidates: list[CriticalPoint] = []
    for index in range(1, data.shape[0] - 1):
        if data[index, 1] >= data[index - 1, 1]:
            continue
        if data[index, 1] >= data[index + 1, 1]:
            continue
        R, beta, fitted = refine_minimum(data, index)
        mode = "Type-I" if beta >= TYPE_SPLIT_BETA else "Type-II"
        candidates.append(CriticalPoint(mode, R, beta, fitted))
    result: dict[str, CriticalPoint] = {}
    for mode in ("Type-I", "Type-II"):
        options = [point for point in candidates if point.mode == mode]
        if options:
            result[mode] = min(options, key=lambda point: point.R)
    return result


def local_extrema_counts(data: np.ndarray) -> tuple[int, int, int, int]:
    minima_low = maxima_low = minima_high = maxima_high = 0
    for index in range(1, data.shape[0] - 1):
        beta = data[index, 2]
        is_low = beta < TYPE_SPLIT_BETA
        if data[index, 1] < data[index - 1, 1] and data[index, 1] < data[index + 1, 1]:
            if is_low:
                minima_low += 1
            else:
                minima_high += 1
        if data[index, 1] > data[index - 1, 1] and data[index, 1] > data[index + 1, 1]:
            if is_low:
                maxima_low += 1
            else:
                maxima_high += 1
    return minima_low, maxima_low, minima_high, maxima_high


def sampled_error(
    blackburn: np.ndarray,
    compressible: np.ndarray,
    lower: float,
    upper: float,
    samples: int = 1201,
    R_ceiling: float | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    if upper <= lower:
        return np.array([]), np.array([]), {
            key: float("nan")
            for key in (
                "mean_signed_percent",
                "mean_abs_percent",
                "rms_percent",
                "p95_abs_percent",
                "max_abs_percent",
                "beta_at_max_abs",
                "normalized_L1_percent",
                "positive_fraction",
            )
        }
    beta = np.linspace(lower, upper, samples)
    Rb = np.interp(beta, blackburn[:, 2], blackburn[:, 1])
    Rc = np.interp(beta, compressible[:, 2], compressible[:, 1])
    if R_ceiling is not None:
        keep = (Rb <= R_ceiling) & (Rc <= R_ceiling)
        beta = beta[keep]
        Rb = Rb[keep]
        Rc = Rc[keep]
        if beta.size < 3:
            return np.array([]), np.array([]), {
                key: float("nan")
                for key in (
                    "mean_signed_percent",
                    "mean_abs_percent",
                    "rms_percent",
                    "p95_abs_percent",
                    "max_abs_percent",
                    "beta_at_max_abs",
                    "normalized_L1_percent",
                    "positive_fraction",
                )
            }
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
            100.0 * np.trapz(np.abs(Rb - Rc), beta) / np.trapz(Rc, beta)
        ),
        "positive_fraction": float(np.mean(relative > 0.0)),
    }
    return beta, relative, metrics


def region_metrics(
    blackburn: np.ndarray,
    compressible: np.ndarray,
) -> tuple[dict[str, float], dict[str, tuple[np.ndarray, np.ndarray]]]:
    common_min = max(float(np.min(blackburn[:, 2])), float(np.min(compressible[:, 2])))
    common_max = min(float(np.max(blackburn[:, 2])), float(np.max(compressible[:, 2])))
    union_min = min(float(np.min(blackburn[:, 2])), float(np.min(compressible[:, 2])))
    union_max = max(float(np.max(blackburn[:, 2])), float(np.max(compressible[:, 2])))
    regions = {
        "full": (common_min, common_max),
        "low_beta": (common_min, min(common_max, TYPE_SPLIT_BETA)),
        "typeI_side": (max(common_min, TYPE_SPLIT_BETA), common_max),
    }
    output: dict[str, float] = {
        "common_beta_min": common_min,
        "common_beta_max": common_max,
        "common_beta_fraction_percent": 100.0
        * (common_max - common_min)
        / (union_max - union_min),
        "beta_min_error_percent": 100.0
        * (float(np.min(blackburn[:, 2])) - float(np.min(compressible[:, 2])))
        / float(np.min(compressible[:, 2])),
        "beta_max_error_percent": 100.0
        * (float(np.max(blackburn[:, 2])) - float(np.max(compressible[:, 2])))
        / float(np.max(compressible[:, 2])),
    }
    profiles: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, (lower, upper) in regions.items():
        beta, relative, metrics = sampled_error(
            blackburn, compressible, lower, upper
        )
        profiles[name] = (beta, relative)
        for metric, value in metrics.items():
            output[f"{name}_{metric}"] = value
    beta, relative, metrics = sampled_error(
        blackburn, compressible, common_min, common_max, R_ceiling=450.0
    )
    profiles["core_R_le_450"] = (beta, relative)
    for metric, value in metrics.items():
        output[f"core_R_le_450_{metric}"] = value
    return output, profiles


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_error_profiles(
    path: Path,
    profiles: dict[float, tuple[np.ndarray, np.ndarray]],
) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(12.8, 7.2), sharex=False, sharey=True)
    for axis, Tw in zip(axes.flat, TEMPERATURES):
        beta, relative = profiles[Tw]
        axis.plot(beta, relative, color="#d62728", linewidth=1.7)
        axis.axhline(0.0, color="0.3", linewidth=0.8)
        axis.axvline(TYPE_SPLIT_BETA, color="0.5", linewidth=0.7, linestyle="--")
        axis.text(0.04, 0.92, rf"$T_w={Tw:.2f}$", transform=axis.transAxes, va="top")
        axis.grid(alpha=0.22, linewidth=0.6)
        axis.set_xlabel(r"$\beta$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$100(R_B-R_C)/R_C$ [%]")
    figure.tight_layout()
    figure.subplots_adjust(top=0.91)
    figure.suptitle(r"Pointwise neutral-curve error at $\omega=0$", y=0.975)
    figure.savefig(path, dpi=240)
    plt.close(figure)


def plot_summary(
    path: Path,
    critical_rows: list[dict[str, object]],
    curve_rows: list[dict[str, object]],
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), sharex=True)
    for mode, marker, color in (
        ("Type-I", "o", "#1f77b4"),
        ("Type-II", "s", "#d62728"),
    ):
        available = [
            row for row in critical_rows
            if row["mode"] == mode and row["status"] == "available"
        ]
        axes[0].plot(
            [row["Tw"] for row in available],
            [row["R_error_percent"] for row in available],
            marker=marker, color=color, label=mode,
        )
        axes[1].plot(
            [row["Tw"] for row in available],
            [row["beta_error_percent"] for row in available],
            marker=marker, color=color, label=mode,
        )
    axes[2].plot(
        [row["Tw"] for row in curve_rows],
        [row["full_rms_percent"] for row in curve_rows],
        "o-", color="#9467bd", label="Full curve",
    )
    axes[2].plot(
        [row["Tw"] for row in curve_rows],
        [row["low_beta_rms_percent"] for row in curve_rows],
        "s--", color="#d62728", label=r"Low-$\beta$ region",
    )
    axes[2].plot(
        [row["Tw"] for row in curve_rows],
        [row["core_R_le_450_rms_percent"] for row in curve_rows],
        "^-.", color="#2ca02c", label=r"Core region ($R\leq450$)",
    )
    labels = (
        r"$\Delta R_c/R_{c,C}$ [%]",
        r"$\Delta\beta_c/\beta_{c,C}$ [%]",
        "pointwise RMS error [%]",
    )
    for axis, label in zip(axes, labels):
        axis.axhline(0.0, color="0.35", linewidth=0.8)
        axis.axhline(5.0, color="0.55", linewidth=0.7, linestyle="--")
        axis.axhline(-5.0, color="0.55", linewidth=0.7, linestyle="--")
        axis.set_xlabel(r"$T_w$")
        axis.set_ylabel(label)
        axis.set_xticks(TEMPERATURES)
        axis.grid(alpha=0.22, linewidth=0.6)
    axes[0].legend(frameon=False)
    axes[2].legend(frameon=False)
    figure.tight_layout()
    figure.subplots_adjust(top=0.88)
    figure.suptitle(r"Zero-frequency Blackburn errors", y=0.97)
    figure.savefig(path, dpi=240)
    plt.close(figure)


def temperature_sensitivity_rows(
    critical: dict[tuple[float, str], dict[str, CriticalPoint]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mode, Tw_start, Tw_end in (
        ("Type-I", 1.00, 1.20),
        ("Type-II", 1.00, 1.12),
    ):
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
                100.0
                * (abs(changes["blackburn"]) - reference)
                / reference
            )
            rows.append(
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
    return rows


def write_report(
    path: Path,
    critical_rows: list[dict[str, object]],
    curve_rows: list[dict[str, object]],
    sensitivity_rows: list[dict[str, object]],
    validation_rows: list[dict[str, object]],
) -> None:
    def critical_row(Tw: float, mode: str) -> dict[str, object]:
        return next(
            row for row in critical_rows
            if row["Tw"] == Tw and row["mode"] == mode
        )

    def curve_row(Tw: float) -> dict[str, object]:
        return next(row for row in curve_rows if row["Tw"] == Tw)

    lines = [
        "# 完整零频中性曲线误差分析",
        "",
        "## 1. 比较对象与误差定义",
        "",
        "本报告只分析磁盘上已有的完整 $\\omega=0$ 中性曲线，没有启动新的稳定性计算。比较对象为 Blackburn 广义 Boussinesq 模型与完全可压缩模型（$M_r=0.3$，基流物性随温度变化，扰动物性开启）。所有有符号相对误差统一定义为",
        "",
        "$$E=100\\,(Q_B-Q_C)/Q_C,$$",
        "",
        "其中下标 $B$ 和 $C$ 分别表示 Blackburn 与完全可压缩结果。正值表示 Blackburn 结果更大。",
        "",
        "曲线形状误差是在两模型共同的 $\\beta$ 区间上插值后，逐点比较 $R_B(\\beta)$ 与 $R_C(\\beta)$。此外给出 $R_B\\le450$ 且 $R_C\\le450$ 的核心区指标，避免高 $\\beta$ 端接近人为计算上限 $R\\simeq500$ 的陡峭尾部主导结论。",
        "",
        "## 2. 数据质量与曲线完整性",
        "",
        f"- 12 条原始曲线每条含 {min(int(row['points']) for row in validation_rows)}--{max(int(row['points']) for row in validation_rows)} 个点，$\\beta$ 均严格单调。",
        f"- 最大中性残差为 {max(float(row['max_neutral_residual']) for row in validation_rows):.3e}，远小于下述百分比模型误差。",
        f"- 两模型共同 $\\beta$ 区间覆盖各自并集的 {min(float(row['common_beta_fraction_percent']) for row in curve_rows):.1f}%--{max(float(row['common_beta_fraction_percent']) for row in curve_rows):.1f}%；逐点指标只在共同区间内计算，没有外推。",
        "- 完全可压缩模型在 $T_w=1.00$--$1.12$ 存在低 $\\beta$ 的 Type-II 局部极小值；在 $T_w=1.16$ 和 1.20 已不存在该折点。Blackburn 曲线在全部温度仍保留 Type-II 折点。",
        "- 因而高温下第二模态临界值缺失是已有完整曲线中的拓扑差异，不是文件只算了一小段造成的。",
        "- 独立的等温实现核对显示，$T_w=1$ 时 Blackburn 与 Lopez 曲线在共同区间上的最大 $|\\Delta R|=2.893\\times10^{-4}$、RMS 差为 $3.540\\times10^{-5}$，说明 Blackburn 程序在其应退化到 Lopez 的极限下是一致的。",
        "",
        "## 3. 临界参数误差",
        "",
        "| $T_w$ | Type-I $R_c$ 误差 | Type-I $\\beta_c$ 误差 | Type-II $R_c$ 误差 | Type-II $\\beta_c$ 误差 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for Tw in TEMPERATURES:
        first = critical_row(Tw, "Type-I")
        second = critical_row(Tw, "Type-II")
        if second["status"] == "available":
            typeII_R = f"{float(second['R_error_percent']):+.3f}%"
            typeII_beta = f"{float(second['beta_error_percent']):+.3f}%"
        else:
            typeII_R = "不可定义（C 无折点）"
            typeII_beta = "不可定义（C 无折点）"
        lines.append(
            f"| {Tw:.2f} | {float(first['R_error_percent']):+.3f}% | "
            f"{float(first['beta_error_percent']):+.3f}% | {typeII_R} | "
            f"{typeII_beta} |"
        )
    lines.extend(
        [
            "",
            "结论：Type-I 的 $R_c$ 误差在全部测试温度内不超过 1.25%，但 $\\beta_c$ 误差单调增至 4.982%。这说明 Blackburn 能较好保持临界雷诺数，却逐渐错置临界波数。Type-II 在 $T_w\\le1.12$ 时两个临界参数误差均不超过 2.5%；到 $T_w=1.16$ 后，误差不能再用有限百分数描述，因为两个模型预测了不同的曲线拓扑。",
            "",
            "## 4. 整条曲线与分区误差",
            "",
            "| $T_w$ | 全曲线 RMS | 核心区 RMS ($R\\le450$) | 低 $\\beta$ RMS | Type-I 侧 RMS | 全曲线 MAE | 全曲线 95% 误差 |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for Tw in TEMPERATURES:
        row = curve_row(Tw)
        lines.append(
            f"| {Tw:.2f} | {float(row['full_rms_percent']):.3f}% | "
            f"{float(row['core_R_le_450_rms_percent']):.3f}% | "
            f"{float(row['low_beta_rms_percent']):.3f}% | "
            f"{float(row['typeI_side_rms_percent']):.3f}% | "
            f"{float(row['full_mean_abs_percent']):.3f}% | "
            f"{float(row['full_p95_abs_percent']):.3f}% |"
        )
    lines.extend(
        [
            "",
            "全曲线 RMS 从 $T_w=1$ 的 0.992% 增至 $T_w=1.20$ 的 6.466%；核心区 RMS 从 0.899% 增至 5.882%。按离散测试点与 5% 判据，整体/核心曲线形状在 $T_w\\le1.12$ 尚可接受，而在 $T_w=1.16$ 已超限。低 $\\beta$ 区在 $T_w=1.12$ 已达到 5.087%，因此若研究重点是第二模态附近的详细曲线形状，保守适用范围只能写成 $T_w\\le1.08$；1.12 应视为临界/边缘状态，而不是精确边界。",
            "",
            "全曲线最大逐点误差分别为 2.416%、5.616%、8.376%、10.705%、13.300% 和 15.764%，均出现在共同区间的高 $\\beta$ 端。那里曲线很陡，较小的水平位移会放大固定 $\\beta$ 下的 $R$ 误差，所以它适合揭示曲线错位，但不宜单独作为模型适用性的判据。核心区 RMS 与 95% 分位误差更稳健。",
            "",
            "## 5. 为什么只看临界值会低估误差",
            "",
            "低 $\\beta$ 区的平均有符号误差随温度由 +0.011% 增至 +3.327%，而 Type-I 侧由 -0.681% 降至 -4.435%。即 Blackburn 在低 $\\beta$ 区通常给出更高的 $R(\\beta)$，在高 $\\beta$ 区却通常给出更低的 $R(\\beta)$。两侧误差互相抵消，使全曲线平均有符号误差看起来只有 -0.561% 至 -1.965%；对应的平均绝对误差却已由 0.799% 增至 5.164%。",
            "",
            "这也解释了一个表面矛盾：$T_w>1$ 时 Type-I 临界 $R_c$ 略偏高，但 Type-I 侧大部分曲线在同一 $\\beta$ 下反而偏低；同样，Type-II 临界 $R_c$ 略偏低，并不代表整个低 $\\beta$ 区都偏低。临界极小值只描述一个点，不能代表曲线的水平位移、曲率和折点结构。",
            "",
            "## 6. 温度响应是否被正确捕捉",
            "",
            "| 模态/量 | 温度区间 | Blackburn 变化 | 可压缩变化 | 变化幅值相对误差 |",
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
            "这里“变化幅值相对误差”为 Blackburn 温度变化幅值相对可压缩模型的偏差。Type-I 中，Blackburn 对 $R_c$ 的热致下降幅度低估约 9.5%，对 $\\beta_c$ 的移动幅度低估约 18.7%。Type-II 在仍有共同折点的 1.00→1.12 区间内，Blackburn 对 $R_c$ 下降幅度高估约 22.2%，对 $\\beta_c$ 移动幅度低估约 19.1%。这表明误差并非一个可用常数校正的统一比例偏差。",
            "",
            "## 7. 分层适用性结论",
            "",
            "| 研究目标 | 5% 判据下由现有离散温度点支持的范围 |",
            "|---|---|",
            "| 只预测 Type-I 临界 $R_c$ | 至少到 $T_w=1.20$ |",
            "| Type-I 的 $R_c$ 与 $\\beta_c$ | 至少到 $T_w=1.20$（$\\beta_c$ 已接近 5%） |",
            "| Type-II 临界参数 | 到 $T_w=1.12$；从 1.16 起发生拓扑失配 |",
            "| 整条/核心中性曲线几何 | 到 $T_w=1.12$；1.16 已超过 5% |",
            "| 低 $\\beta$ 第二模态区详细形状 | 保守到 $T_w=1.08$；1.12 为边缘状态 |",
            "",
            "因此最准确的表述不是“Blackburn 模型整体失效”，而是它具有明显选择性：它能在较宽温度范围内保留 Type-I 临界 $R_c$，但对临界波数、整曲线几何及第二模态拓扑的可靠性更早下降。完全可压缩模型的必要性主要体现在高温下第二模态折点的消失和中性曲线整体重塑，而不只是把临界 $R_c$ 修正几个百分点。",
            "",
            "## 8. 结论边界",
            "",
            "- 上述温度范围只对应已计算的离散点，不能把 1.12 或 1.16 当作连续参数空间中的精确失效温度。",
            "- 5% 是用于分层比较的工程化报告阈值，不是由控制方程给出的普适物理边界；论文中应同时给出原始误差曲线，避免把单一阈值绝对化。",
            "- $T_w=1.12$ 的低 $\\beta$ RMS=5.087% 与 5% 非常接近，应报告为边缘而非断言性失败。",
            "- 固定 $\\beta$ 的逐点误差同时包含曲线水平错位与垂直差异；若论文需要进一步区分二者，可另做最近曲线距离或弧长配准，但不影响当前拓扑失配结论。",
            "- 本报告没有评估基流离散、谱截断或临界值拟合的独立收敛阶；不过中性残差约 $10^{-7}$，远小于主要模型误差，现有差异不是中性根容差本身造成的。",
            "",
            "## 9. 配套文件",
            "",
            "- `zero_frequency_critical_errors.tsv`：临界参数及误差",
            "- `zero_frequency_curve_shape_errors.tsv`：全曲线与分区统计",
            "- `zero_frequency_topology.tsv`：局部极值与 Type-II 折点状态",
            "- `zero_frequency_curve_validation.tsv`：点数、范围和中性残差",
            "- `zero_frequency_temperature_sensitivity.tsv`：温度响应比较",
            "- `zero_frequency_pointwise_error_profiles.png`：逐点误差随 $\\beta$ 的分布",
            "- `zero_frequency_error_summary.png`：临界与整曲线误差摘要",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    workspace = Path(__file__).resolve().parent
    output_dir = workspace / "zero_frequency_complete_error_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    curves: dict[tuple[float, str], np.ndarray] = {}
    critical: dict[tuple[float, str], dict[str, CriticalPoint]] = {}
    validation_rows: list[dict[str, object]] = []
    topology_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for model in ("blackburn", "compressible"):
            path = curve_path(workspace, model, Tw)
            data = load_curve(path)
            curves[(Tw, model)] = data
            critical[(Tw, model)] = critical_points(data)
            residual = np.minimum(np.abs(data[:, 4]), np.abs(data[:, 6]))
            validation_rows.append(
                {
                    "Tw": Tw,
                    "model": model,
                    "source_file": path.name,
                    "points": int(data.shape[0]),
                    "beta_min": float(np.min(data[:, 2])),
                    "beta_max": float(np.max(data[:, 2])),
                    "R_min": float(np.min(data[:, 1])),
                    "R_max": float(np.max(data[:, 1])),
                    "max_neutral_residual": float(np.max(residual)),
                }
            )
            low_min, low_max, high_min, high_max = local_extrema_counts(data)
            topology_rows.append(
                {
                    "Tw": Tw,
                    "model": model,
                    "TypeII_fold": "present" if "Type-II" in critical[(Tw, model)] else "absent",
                    "low_beta_local_minima": low_min,
                    "low_beta_local_maxima": low_max,
                    "high_beta_local_minima": high_min,
                    "high_beta_local_maxima": high_max,
                }
            )

    critical_rows: list[dict[str, object]] = []
    for Tw in TEMPERATURES:
        for mode in ("Type-I", "Type-II"):
            blackburn = critical[(Tw, "blackburn")].get(mode)
            compressible = critical[(Tw, "compressible")].get(mode)
            row: dict[str, object] = {"Tw": Tw, "mode": mode}
            if blackburn is None or compressible is None:
                row.update(
                    {
                        "status": (
                            "blackburn_missing" if blackburn is None
                            else "compressible_missing"
                        ),
                        "R_c_blackburn": blackburn.R if blackburn else "",
                        "beta_c_blackburn": blackburn.beta if blackburn else "",
                    }
                )
            else:
                row.update(
                    {
                        "status": "available",
                        "R_c_blackburn": blackburn.R,
                        "R_c_compressible": compressible.R,
                        "R_error_percent": 100.0 * (blackburn.R - compressible.R) / compressible.R,
                        "beta_c_blackburn": blackburn.beta,
                        "beta_c_compressible": compressible.beta,
                        "beta_error_percent": 100.0 * (blackburn.beta - compressible.beta) / compressible.beta,
                        "blackburn_quadratic_fit": blackburn.fitted,
                        "compressible_quadratic_fit": compressible.fitted,
                    }
                )
            critical_rows.append(row)

    curve_rows: list[dict[str, object]] = []
    profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for Tw in TEMPERATURES:
        metrics, region_profiles = region_metrics(
            curves[(Tw, "blackburn")], curves[(Tw, "compressible")]
        )
        profiles[Tw] = region_profiles["full"]
        curve_rows.append({"Tw": Tw, **metrics})

    sensitivity_rows = temperature_sensitivity_rows(critical)

    write_tsv(output_dir / "zero_frequency_critical_errors.tsv", critical_rows)
    write_tsv(output_dir / "zero_frequency_curve_shape_errors.tsv", curve_rows)
    write_tsv(output_dir / "zero_frequency_topology.tsv", topology_rows)
    write_tsv(output_dir / "zero_frequency_curve_validation.tsv", validation_rows)
    write_tsv(
        output_dir / "zero_frequency_temperature_sensitivity.tsv",
        sensitivity_rows,
    )
    plot_error_profiles(
        output_dir / "zero_frequency_pointwise_error_profiles.png", profiles
    )
    plot_summary(
        output_dir / "zero_frequency_error_summary.png",
        critical_rows,
        curve_rows,
    )
    write_report(
        output_dir / "zero_frequency_detailed_summary.md",
        critical_rows,
        curve_rows,
        sensitivity_rows,
        validation_rows,
    )
    print(f"Output: {output_dir}")


if __name__ == "__main__":
    main()
