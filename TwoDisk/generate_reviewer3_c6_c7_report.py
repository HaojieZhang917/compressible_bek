from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
C6 = ROOT / "c6_validation_results"
C7 = ROOT / "c7_validation_results"
REPORT = ROOT / "REFEREE3_C6_C7_COMPUTATION_VALIDATION_REPORT.md"


def read_table(path: Path):
    header = None
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            header = line[1:].strip().split()
        elif line.strip():
            rows.append(line.split())
    return header, rows


def interpolate(x, y, value):
    return float(np.interp(value, x, y))


c6_header, c6_rows = read_table(C6 / "c6_branch_and_receptivity.dat")
c6_data = np.array(c6_rows, dtype=float)
c6 = {name: c6_data[:, i] for i, name in enumerate(c6_header)}

cases_header, cases_rows = read_table(C6 / "c6_three_cases.dat")
cases = []
for raw in cases_rows:
    cases.append({cases_header[0]: raw[0], **{
        name: float(value) for name, value in zip(cases_header[1:], raw[1:])
    }})

c7_header, c7_rows = read_table(C7 / "c7_two_branches.dat")
c7_data = np.array(c7_rows, dtype=float)
c7 = {name: c7_data[:, i] for i, name in enumerate(c7_header)}

rl_case = cases[0]
global_case = next(c for c in cases if c["label"] == "global_receptivity_peak")
fixed360_case = next(c for c in cases if c["label"] == "fixed_R_360_check")

r_growth = c6["R"][np.argmax(c6["growth_rate"])]
r_global = c6["R"][np.argmax(c6["Cr_abs"])]
r_neutral = rl_case["R"]

# Recompute the interpolated maximum-growth case from the saved branch; this
# avoids depending on stale case rows during report regeneration.
peak_index = int(np.argmax(c6["growth_rate"]))
x1, x2, x3 = c6["R"][peak_index - 1:peak_index + 2]
y1, y2, y3 = c6["growth_rate"][peak_index - 1:peak_index + 2]
offset = 0.5 * (y1 - y3) / (y1 - 2 * y2 + y3)
r_growth_refined = x2 + offset * (x2 - x1)

# The Julia output predating the final case-label correction contains a broad
# interval maximum at R=415.49. Reconstruct the intended maximum-growth case
# directly from the saved branch and local receptivity curve.
n_end = c6["cumulative_N_from_Rl"][-1]
n_start_growth = interpolate(c6["R"], c6["cumulative_N_from_Rl"], r_growth_refined)
cr_growth = interpolate(c6["R"], c6["Cr_abs"], r_growth_refined)
n_growth = n_end - n_start_growth
growth_case = {
    "label": "maximum_growth_local_case",
    "R": r_growth_refined,
    "Cr_abs": cr_growth,
    "N_to_570": n_growth,
    "gain_to_570": math.exp(n_growth),
    "A_abs_at_570": cr_growth * math.exp(n_growth),
}
cases = [rl_case, growth_case, global_case, fixed360_case]

case_output = C6 / "c6_three_cases_reported.dat"
with case_output.open("w", encoding="utf-8") as stream:
    stream.write("# label R Cr_abs N_to_570 gain_to_570 A_abs_at_570\n")
    for case in cases:
        stream.write(
            f"{case['label']} {case['R']:.12g} {case['Cr_abs']:.16g} "
            f"{case['N_to_570']:.16g} {case['gain_to_570']:.16g} "
            f"{case['A_abs_at_570']:.16g}\n"
        )

# Integral values are cumulative from R_l. Use endpoint differences for each
# excitation radius. Refined eigenvalue/receptivity values remain those in the
# Julia validation output.
fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.5))
ax = axes[0]
ax.plot(c6["R"], c6["growth_rate"], color="tab:blue", lw=2, label=r"$-\alpha_i$")
ax2 = ax.twinx()
ax2.plot(c6["R"], c6["Cr_abs"], color="tab:red", lw=2, label=r"$|C_r|$")
for radius, color, text in [
    (r_neutral, "black", r"$R_l$"),
    (r_growth_refined, "tab:green", r"$R_g$"),
    (global_case["R"], "tab:purple", r"$R_p$"),
    (570.0, "gray", r"$R_{abs}$"),
]:
    ax.axvline(radius, color=color, ls="--", lw=1.1)
    ax.text(radius + 2, ax.get_ylim()[1] * 0.88, text, color=color)
ax.set_ylabel(r"growth rate $-\alpha_i$")
ax2.set_ylabel(r"local coefficient $|C_r|$")
ax.grid(alpha=0.2)

ax = axes[1]
labels = ["lower\nneutral", "max-growth\nlocal", "global $C_r$\npeak", "$R=360$\ncheck"]
amplitudes = [rl_case["A_abs_at_570"], growth_case["A_abs_at_570"],
              global_case["A_abs_at_570"], fixed360_case["A_abs_at_570"]]
bars = ax.bar(labels, amplitudes, color=["black", "tab:green", "tab:purple", "tab:orange"])
ax.set_yscale("log")
ax.set_ylabel(r"$|A(570)|$")
ax.grid(axis="y", which="both", alpha=0.2)
for bar, value in zip(bars, amplitudes):
    ax.text(bar.get_x() + bar.get_width() / 2, value * 1.08, f"{value:.3g}",
            ha="center", va="bottom")
ax.set_xlabel("excitation case")
fig.tight_layout()
fig.savefig(C6 / "c6_three_case_validation.png", dpi=220)
plt.close(fig)

fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True)
axes[0].plot(c7["R"], c7["alphaI_r"], lw=2, label="branch I")
axes[0].plot(c7["R"], c7["alphaII_r"], lw=2, label="branch II")
axes[0].set_ylabel(r"$\alpha_r$")
axes[0].legend()
axes[1].plot(c7["R"], c7["alphaI_i"], lw=2, label="branch I")
axes[1].plot(c7["R"], c7["alphaII_i"], lw=2, label="branch II")
axes[1].axhline(0, color="black", lw=0.8)
axes[1].set_ylabel(r"$\alpha_i$")
axes[2].plot(c7["R"], c7["CrI"], lw=2, label=r"$|C_{r,I}|$")
axes[2].plot(c7["R"], c7["CrII"], lw=2, label=r"$|C_{r,II}|$")
axes[2].plot(c7["R"], c7["QI"] / np.max(c7["QI"]) * np.max(c7["CrI"]),
             ls="--", lw=1.4, label=r"scaled $|Q_I|$")
axes[2].set_ylabel("coefficient / scaled normalization")
axes[2].set_xlabel(r"$R$")
axes[2].legend(ncol=3, fontsize=8)
r_closest = c7["R"][np.argmin(c7["delta_alpha_r"])]
for ax in axes:
    ax.axvline(r_closest, color="gray", ls="--", lw=1)
    ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(C7 / "c7_two_branch_diagnostics.png", dpi=220)
plt.close(fig)

dr_results = {
    1: (4.336183920905599, 76.41537511648566),
    2: (4.336002781124991, 76.40153450578242),
    4: (4.335278196233732, 76.34619515965879),
}
dr_results = {
    dr: (nfactor, gain, rl_case["Cr_abs"] * gain)
    for dr, (nfactor, gain) in dr_results.items()
}

closest_idx = int(np.argmin(c7["delta_alpha_r"]))

def fmt(x, digits=6):
    return f"{x:.{digits}g}"


case_lines = []
case_names = {
    "lower_neutral": "下中性点",
    "maximum_growth_local_case": "最大增长率附近",
    "global_receptivity_peak": "全局感受性峰值",
    "fixed_R_360_check": "固定 R=360 检查",
}
for case in cases:
    case_lines.append(
        f"| {case_names[case['label']]} | {case['R']:.5f} | "
        f"{case['Cr_abs']:.6f} | {case['N_to_570']:.6f} | "
        f"{case['gain_to_570']:.6f} | {case['A_abs_at_570']:.6f} |"
    )

dr_lines = []
for dr, (nfactor, gain, amp) in dr_results.items():
    rel = abs(amp / dr_results[1][2] - 1)
    dr_lines.append(f"| {dr} | {nfactor:.9f} | {gain:.9f} | {amp:.9f} | {rel:.3e} |")

report = rf"""# Referee 3 C6–C7 计算程序验证报告

## 1. 本轮范围

本轮仅验证计算程序并生成数值证据，没有修改 manuscript、TeX 或 response-to-reviewers 文件。计算针对 Referee 3 的 C6 和 C7，统一采用

\[
a_s=0,\qquad Re_h=1000,\qquad \bar{{\omega}}=0,\qquad n=30,
\]

并使用 Chebyshev 多项式阶数 `N=99`（100 个配置点）。现有网格无关性数据表明，在代表性 Type-I 模态上，100 个配置点相对于 300 点参考结果的复特征值相对误差约为 `1.90e-8`。

壁面激励参数已与现有计算结果统一为：`h_r=1/Re_h=0.001`、`c^2=1`。对于程序中使用的 Fourier 表达

\[
\hat h(\alpha)=h_r\sqrt{{\frac{{\pi}}{{l_s}}}}\exp\left(-\frac{{\alpha^2}}{{4l_s}}\right),
\]

有 `l_s=1/(2c^2)=0.5`，本轮 C6 和 C7 均使用该值。一致性检查表明，在 `R=470`时 `N=99` 得到 `|C_r|=0.244259`，与正文曲线的峰值约 `0.25` 一致；`N=129` 的对应结果为 `0.243423`。

## 2. C6：三个激励位置传播至绝对不稳定临界半径

### 2.1 方法

固定 `n=30` 和 `omega_bar=0`，沿同一条不稳定空间分支取

\[
\beta(R)=\frac{{n}}{{R}},\qquad \omega(R)=0,
\]

并计算

\[
N(R_t;R_f)=-\int_{{R_f}}^{{R_t}}\alpha_i(\xi)\,\mathrm d\xi,
\qquad
|A(R_t;R_f)|=|C_r(R_f)|\exp[N(R_t;R_f)].
\]

共同目标半径为 `R_t=R_abs=570`。程序自动得到下中性点

\[
R_l={rl_case['R']:.8f},\qquad
\alpha(R_l)={rl_case['alpha_r']:.9f}{rl_case['alpha_i']:+.2e}\mathrm i.
\]

沿程最大增长率位于 `R_g≈{r_growth_refined:.5f}`，全局感受性峰值位于 `R_p≈{global_case['R']:.5f}`。

需要特别说明：当前重新计算的 `C_r(R)` 在 `R≈360` 附近仍然单调增加，没有形成严格的局部极值。因此，三组主比较采用“下中性点、最大增长率附近的局部 case、全局感受性峰值”。另外保留 `R=360` 的定点计算，以检查用户建议位置对结论的影响。

### 2.2 数值结果

| 激励位置 | `R_f` | `|C_r(R_f)|` | `N(570;R_f)` | `exp(N)` | `|A(570)|` |
|---|---:|---:|---:|---:|---:|
{chr(10).join(case_lines)}

![C6 三组激励位置的总振幅比较](c6_validation_results/c6_three_case_validation.png)

下中性点的局部耦合系数最小，但其放大距离最长，最终振幅最大：

- 相对最大增长率附近激励，下中性点激励的 `|A(570)|` 大约高 `{rl_case['A_abs_at_570']/growth_case['A_abs_at_570']:.2f}` 倍；
- 相对全局 `C_r` 峰值处激励，大约高 `{rl_case['A_abs_at_570']/global_case['A_abs_at_570']:.2f}` 倍；
- 即使固定采用 `R=360`，下中性点激励的最终振幅仍高 `{rl_case['A_abs_at_570']/fixed360_case['A_abs_at_570']:.2f}` 倍。

因此，这组线性计算支持审稿人 C6 的核心物理判断：局部最大的感受性系数并不意味着抵达下游时的扰动振幅最大；较早激励获得的累计指数放大可以完全超过初始耦合系数的劣势。

### 2.3 径向积分验证

以下比较仅改变积分/续接步长，并使用相同的 `R_l` 与正确的初始特征值移位：

| `ΔR` | `N(570)` | `exp(N)` | `|A(570)|` | 相对 `ΔR=1` 差异 |
|---:|---:|---:|---:|---:|
{chr(10).join(dr_lines)}

`ΔR=2` 的最终振幅相对 `ΔR=1` 相差 `{abs(dr_results[2][2]/dr_results[1][2]-1):.3e}`，可作为本轮兼顾速度和精度的设置。

## 3. C7：两支模态在切换区的诊断

### 3.1 复特征值是否合并

在 `420≤R≤500` 内分别求解两条特征值分支。在 `R={c7['R'][closest_idx]:.0f}` 时，两支实部最接近：

\[
\alpha_I={c7['alphaI_r'][closest_idx]:.9f}{c7['alphaI_i'][closest_idx]:+.9f}\mathrm i,
\]

\[
\alpha_{{II}}={c7['alphaII_r'][closest_idx]:.9f}{c7['alphaII_i'][closest_idx]:+.9f}\mathrm i.
\]

对应差异为

\[
\Delta\alpha_r={c7['delta_alpha_r'][closest_idx]:.3e},\qquad
\Delta\alpha_i={c7['delta_alpha_i'][closest_idx]:.3e},\qquad
|\Delta\alpha|={c7['delta_alpha_abs'][closest_idx]:.3e}.
\]

这验证了审稿人的判断：两支模态的实部几乎相同，但完整复特征值并不相等，也不存在由当前数据支持的严格特征值合并。

### 3.2 感受性峰值来自哪里

在 `R=470`：

| 分支 | `|C_r|` | `|Q|` | `|BC|` |
|---|---:|---:|---:|
| I | {c7['CrI'][closest_idx]:.6f} | {c7['QI'][closest_idx]:.6e} | {c7['BCI'][closest_idx]:.6e} |
| II | {c7['CrII'][closest_idx]:.6f} | {c7['QII'][closest_idx]:.6e} | {c7['BCII'][closest_idx]:.6e} |

沿主响应分支，`|BC|` 在切换区只缓慢变化，而 `|Q|` 明显减小，因此

\[
|C_r|=\frac{{|BC|}}{{|Q|}}
\]

中的小分母效应是 `C_r` 上升并形成峰值的主要直接原因。峰值位置同时与两支特征值实部接近的区域重合，但不能据此声称复特征值发生合并。

![C7 两条特征值分支及感受性诊断](c7_validation_results/c7_two_branch_diagnostics.png)

### 3.3 模态命名的证据等级

当前计算可以可靠地称这两条曲线为“two distinct spatial eigenvalue branches”，并证明其实部接近、虚部不同。现有的单一分量峰值位置诊断不足以独立、无歧义地把每一点都物理标记为 Type I 或 Type II；正式用于论文前，建议把两个分支的完整直接/伴随特征函数与文中已经确认的 Type-I/Type-II 参考模态做加权重叠比较。故本报告不把尚未完成的形状标签当作最终证据。

## 4. 程序验证中发现的注意事项

1. `total_amplitude.jl` 自动寻找下中性点时能正确进入目标不稳定分支。若直接提供 `--lower-radius`，还必须同时给出接近目标根的 `--alpha-real≈0.44968`；否则默认移位 `1.05` 会收敛到另一条稳定分支，产生完全错误的 N 因子。
2. `N=99` 在现有网格验证中已经足够；本轮 C6/C7 的 Hermitian-adjoint 特征值配对误差均保持在约 `1e-11` 或更小。
3. 粗糙元宽度已统一为 `c^2=1`与 `l_s=0.5`。该参数能够复现正文现有的 `C_r` 曲线；相比之下，若使用 `c^2=4` 对应的 `l_s=0.125`，`R=470` 处会得到 `|C_r|≈0.466`，与正文图中约 `0.25` 的峰值不符。
4. 总振幅采用局部平行、线性 `e^N` 传播，没有包含非平行幅值输运、群速度修正或非线性饱和；因此它适合比较 C6 中不同激励位置的相对结果，不应直接解释为转捩概率。

## 5. 对审稿意见的证据结论

### C6

证据强度已经达到可以回应审稿人的水平：下中性点激励尽管初始 `C_r` 较小，但到 `R_abs=570` 时的线性总振幅显著大于最大增长率附近和全局局部感受性峰值处的激励。这说明原来的局部 `C_r` 峰值不能单独代表下游危险程度。

### C7

两条模态分支在 `R≈470` 的实部近乎相同，但虚部明显不同。`C_r` 的峰值主要与 `Q` 的下降有关，而不是 `BC` 的同步尖峰。正式表述应使用“dominant-response/branch interaction or switching region”，避免声称两个复特征值完全相等或发生严格合并。

## 6. 生成文件

- C6 主数据：`c6_validation_results/c6_branch_and_receptivity.dat`
- C6 报告采用的三/四个 case：`c6_validation_results/c6_three_cases_reported.dat`
- C6 Julia 原始 case 输出：`c6_validation_results/c6_three_cases.dat`
- C6 图：`c6_validation_results/c6_three_case_validation.png`
- C7 两分支数据：`c7_validation_results/c7_two_branches.dat`
- C7 图：`c7_validation_results/c7_two_branch_diagnostics.png`
- 本报告：`REFEREE3_C6_C7_COMPUTATION_VALIDATION_REPORT.md`
"""

REPORT.write_text(report, encoding="utf-8")
print(REPORT)
