"""
外插中性曲线到 R=800，使用渐近模型保证 dβ/dR → 0 的物理行为。
模型: β(R) = β_∞ + A / (R - R₀)^p  (p > 0)
- 下支: β_∞ = 0 (β → 0 as R → ∞)
- 上支: β_∞ 自由拟合 (β → β_∞ as R → ∞)
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

DATA_DIR = "/home/zhj/main/code/compress/compressible_bek/TwoDisk/absolute_stability/data/"
R_TARGET = 800.0


# ============================================================
# 物理渐近模型
# ============================================================
def model_lower(R, A, R0, p):
    """下支: β(R) = A / (R - R0)^p,  A>0, R0 < min(R), p>0"""
    return A / (R - R0)**p


def model_upper(R, beta_inf, A, R0, p):
    """上支: β(R) = β_∞ + A / (R - R0)^p,  A<0 (使 β 从下方趋近 β_∞)"""
    return beta_inf + A / (R - R0)**p


def dmodel_lower(R, A, R0, p):
    """下支导数: dβ/dR = -A*p / (R-R0)^(p+1)"""
    return -A * p / (R - R0)**(p + 1)


def dmodel_upper(R, beta_inf, A, R0, p):
    """上支导数: dβ/dR = -A*p / (R-R0)^(p+1)"""
    return -A * p / (R - R0)**(p + 1)


# ============================================================
# 工具函数
# ============================================================
def split_branches(R_raw, beta_raw):
    idx_min = np.argmin(R_raw)
    beta_min = beta_raw[idx_min]
    lower = beta_raw <= beta_min
    upper = beta_raw >= beta_min
    return (R_raw[lower], beta_raw[lower],
            R_raw[upper], beta_raw[upper],
            R_raw[idx_min], beta_min)


def fit_and_extrapolate(R_data, beta_data, R_target, model, p0, bounds,
                        tail_frac=0.3, max_tail=40):
    """
    用渐近模型拟合数据尾部并外插。
    返回: R_ext, beta_ext, popt, fit_func, dbeta_dR_ext
    """
    sort_idx = np.argsort(R_data)
    R_sorted = R_data[sort_idx]
    beta_sorted = beta_data[sort_idx]

    n_tail = min(int(len(R_sorted) * tail_frac), max_tail)
    n_tail = max(n_tail, len(p0) + 3)
    R_tail = R_sorted[-n_tail:]
    beta_tail = beta_sorted[-n_tail:]

    try:
        popt, pcov = curve_fit(model, R_tail, beta_tail, p0=p0,
                               bounds=bounds, maxfev=50000, ftol=1e-14)
    except Exception as e:
        print(f"    拟合警告: {e}")
        popt, pcov = curve_fit(model, R_tail, beta_tail, p0=p0,
                               bounds=bounds, maxfev=100000, ftol=1e-14,
                               method='trf')

    fit_func = lambda R: model(R, *popt)

    # 生成外插段
    R_ext = np.linspace(R_tail[-1], R_target, 50)
    beta_ext = fit_func(R_ext)

    # 计算导数
    if model == model_lower:
        dbeta_dR = dmodel_lower(R_ext, *popt)
    else:
        dbeta_dR = dmodel_upper(R_ext, *popt)

    return R_ext, beta_ext, dbeta_dR, popt, fit_func, R_tail, beta_tail


def check_monotonicity(R_ext, beta_ext, dbeta_dR, label):
    """检查外插段的单调性和导数行为"""
    slopes = np.diff(beta_ext) / np.diff(R_ext)
    print(f"  {label}:")
    print(f"    外插段斜率范围: [{slopes[0]:.2e}, {slopes[-1]:.2e}]")
    print(f"    |斜率|递减? {'是' if abs(slopes[-1]) < abs(slopes[0]) else '否'}")
    print(f"    末端 dβ/dR = {dbeta_dR[-1]:.2e}")


# ============================================================
# 处理单个 Ts
# ============================================================
def process_one_Ts(Ts_str):
    fname = f"AS_{Ts_str}_neutral.dat"
    file_path = DATA_DIR + fname
    print(f"\n{'='*60}")
    print(f"处理 Ts = {Ts_str}")
    print(f"{'='*60}")

    data = np.loadtxt(file_path)
    R_raw, beta_raw = data[:, 0], data[:, 1]
    R_l, beta_l, R_u, beta_u, R_min, beta_min = split_branches(R_raw, beta_raw)

    print(f"  下支: {len(R_l)}点, R [{R_l.min():.1f}, {R_l.max():.1f}]")
    print(f"  上支: {len(R_u)}点, R [{R_u.min():.1f}, {R_u.max():.1f}]")

    # ---- 下支: β(R) = A / (R-R0)^p ----
    R_min_l = R_l.min()
    R0_bound_l = R_min_l * 0.9  # R0 < min(R)
    bounds_lower = ([1e-8, -np.inf, 0.05], [np.inf, R0_bound_l, 5.0])
    A_init_l = beta_l[np.argmax(R_l)] * (R_l.max() - R_min_l * 0.5)
    p0_lower = [A_init_l, R_min_l * 0.5, 1.0]

    Re_l, be_l, dbe_l, popt_l, fit_l, Rt_l, bt_l = \
        fit_and_extrapolate(R_l, beta_l, R_TARGET, model_lower, p0_lower, bounds_lower)
    beta_lower_800 = be_l[-1]
    print(f"  下支拟合: A={popt_l[0]:.4f}, R0={popt_l[1]:.2f}, p={popt_l[2]:.4f}")
    print(f"  下支 R=800: β={beta_lower_800:.8f}")
    check_monotonicity(Re_l, be_l, dbe_l, "下支")

    # 上支: β_∞ 约束在 [max, 1.2*max] 范围内（上支 β 不可能增长太多）
    beta_max_u = beta_u.max()
    R_min_u = R_u.min()
    R0_bound_u = R_min_u * 0.9
    bounds_upper = ([beta_max_u, -np.inf, -np.inf, 0.05],
                    [beta_max_u * 1.2, -1e-12, R0_bound_u, 5.0])
    beta_inf_init = beta_max_u * 1.02
    A_init_u = -(beta_inf_init - beta_u[np.argmin(R_u)]) * (R_u.max() - R_min_u * 0.5)
    p0_upper = [beta_inf_init, A_init_u, R_min_u * 0.5, 1.0]

    Re_u, be_u, dbe_u, popt_u, fit_u, Rt_u, bt_u = \
        fit_and_extrapolate(R_u, beta_u, R_TARGET, model_upper, p0_upper, bounds_upper)
    beta_upper_800 = be_u[-1]
    print(f"  上支拟合: β_∞={popt_u[0]:.6f}, A={popt_u[1]:.4f}, R0={popt_u[2]:.2f}, p={popt_u[3]:.4f}")
    print(f"  上支 R=800: β={beta_upper_800:.8f}")
    check_monotonicity(Re_u, be_u, dbe_u, "上支")

    # ---- 导出 full_curve ----
    out_path = DATA_DIR + f"AS_{Ts_str}_full_curve.dat"
    with open(out_path, "w") as f:
        f.write(f"# Neutral curve with asymptotic extrapolation to R={R_TARGET}\n")
        f.write(f"# Ts = {Ts_str}\n")
        f.write(f"# Model: beta = beta_inf + A/(R-R0)^p\n")
        f.write(f"# Lower: beta_inf=0, A={popt_l[0]:.6f}, R0={popt_l[1]:.3f}, p={popt_l[2]:.4f}\n")
        f.write(f"# Upper: beta_inf={popt_u[0]:.8f}, A={popt_u[1]:.6f}, R0={popt_u[2]:.3f}, p={popt_u[3]:.4f}\n")
        f.write(f"# R={R_TARGET} lower beta = {beta_lower_800:.8f}\n")
        f.write(f"# R={R_TARGET} upper beta = {beta_upper_800:.8f}\n")
        f.write(f"# R        beta         source\n")

        for r, b in zip(*[x[np.argsort(beta_l)] for x in (R_l, beta_l)]):
            f.write(f"{r:.8f}  {b:.8f}  lower_data\n")
        for r, b in zip(Re_l[::-1], be_l[::-1]):
            if r > R_l.max():
                f.write(f"{r:.8f}  {b:.8f}  lower_extrap\n")
        for r, b in zip(*[x[np.argsort(beta_u)] for x in (R_u, beta_u)]):
            f.write(f"{r:.8f}  {b:.8f}  upper_data\n")
        for r, b in zip(Re_u, be_u):
            if r > R_u.max():
                f.write(f"{r:.8f}  {b:.8f}  upper_extrap\n")
    print(f"  输出: {out_path}")

    return (R_l, beta_l, Re_l, be_l, dbe_l, fit_l, popt_l,
            R_u, beta_u, Re_u, be_u, dbe_u, fit_u, popt_u,
            beta_lower_800, beta_upper_800)


# ============================================================
# 批量处理
# ============================================================
results = {}
for Ts in ["-0.2", "-0.4", "0.0"]:
    results[Ts] = process_one_Ts(Ts)

# ============================================================
# 汇总
# ============================================================
print("\n" + "=" * 60)
print("外插汇总: R=800 (渐近模型)")
print("=" * 60)
print(f"{'Ts':>6s}  {'下支 β':>14s}  {'上支 β':>14s}  {'上支 β_∞':>12s}")
print("-" * 55)
for Ts in ["-0.2", "-0.4", "0.0"]:
    r = results[Ts]
    bl, bu = r[-2], r[-1]
    beta_inf_u = r[13][0]  # popt_u[0]
    print(f"{Ts:>6s}  {bl:14.8f}  {bu:14.8f}  {beta_inf_u:12.8f}")

# ============================================================
# 绘图
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
colors = {"-0.2": "C0", "-0.4": "C1", "0.0": "C2"}

ax = axes[0, 0]  # 整体曲线
ax_dr = axes[0, 1]  # 导数 dβ/dR
ax_l = axes[1, 0]   # 下支放大
ax_u = axes[1, 1]   # 上支放大

for Ts in ["-0.2", "-0.4", "0.0"]:
    (R_l, b_l, Re_l, be_l, dbe_l, fit_l, popt_l,
     R_u, b_u, Re_u, be_u, dbe_u, fit_u, popt_u,
     _, _) = results[Ts]
    c = colors[Ts]

    # (a) 整体
    ax.plot(R_l, b_l, '.', color=c, ms=2, alpha=0.4)
    ax.plot(R_u, b_u, '.', color=c, ms=2, alpha=0.4)
    ax.plot(Re_l, be_l, color=c, lw=2, label=f'Ts={Ts}')
    ax.plot(Re_u, be_u, color=c, lw=2)

    # (b) 导数
    ax_dr.plot(Re_l, np.abs(dbe_l), color=c, lw=1.5, ls='-',
               label=f'Ts={Ts} lower')
    ax_dr.plot(Re_u, np.abs(dbe_u), color=c, lw=1.5, ls='--',
               label=f'Ts={Ts} upper')

    # (c) 下支放大
    mask_l = R_l >= np.percentile(R_l, 70)
    # 计算原始数据尾部数值导数
    sort_l = np.argsort(R_l)
    Rs_l = R_l[sort_l]
    bs_l = b_l[sort_l]
    dbeta_orig_l = np.gradient(bs_l, Rs_l)
    ax_l.plot(Rs_l[-20:], np.abs(dbeta_orig_l[-20:]), '.', color=c, ms=4, alpha=0.7)
    ax_l.plot(Re_l, np.abs(dbe_l), color=c, lw=2, label=f'Ts={Ts}')

    # (d) 上支放大
    sort_u = np.argsort(R_u)
    Rs_u = R_u[sort_u]
    bs_u = b_u[sort_u]
    dbeta_orig_u = np.gradient(bs_u, Rs_u)
    ax_u.plot(Rs_u[-20:], np.abs(dbeta_orig_u[-20:]), '.', color=c, ms=4, alpha=0.7)
    ax_u.plot(Re_u, np.abs(dbe_u), color=c, lw=2, label=f'Ts={Ts}')

ax.axvline(x=R_TARGET, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('R'); ax.set_ylabel('β')
ax.set_title('Neutral curves (asymptotic extrapolation)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax_dr.set_yscale('log')
ax_dr.set_xlabel('R'); ax_dr.set_ylabel('|dβ/dR| (log scale)')
ax_dr.set_title('Derivative magnitude (must decrease)')
ax_dr.legend(fontsize=7); ax_dr.grid(True, alpha=0.3)

ax_l.set_yscale('log')
ax_l.set_xlabel('R'); ax_l.set_ylabel('|dβ/dR| (log scale)')
ax_l.set_title('Lower branch: |dβ/dR| detail')
ax_l.legend(fontsize=7); ax_l.grid(True, alpha=0.3)

ax_u.set_yscale('log')
ax_u.set_xlabel('R'); ax_u.set_ylabel('|dβ/dR| (log scale)')
ax_u.set_title('Upper branch: |dβ/dR| detail')
ax_u.legend(fontsize=7); ax_u.grid(True, alpha=0.3)

plt.tight_layout()
out_png = DATA_DIR + "all_Ts_extrapolation.png"
plt.savefig(out_png, dpi=150)
print(f"\n对比图已保存至: {out_png}")
