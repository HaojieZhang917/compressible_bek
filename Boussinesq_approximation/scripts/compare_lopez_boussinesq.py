#!/usr/bin/env python3
"""Compare traditional and Lopez-type generalized Boussinesq disk flows."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp


PR = 0.72
WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
OUT = WORKSPACE_ROOT / "lopez_boussinesq_comparison"
MAX_NODES = 200000


def initial_guess(z):
    y = np.zeros((7, z.size))
    y[0] = -0.8845 * (1.0 - np.exp(-0.8 * z))
    y[1] = 0.5102 * (1.0 - 0.8 * z) * np.exp(-0.8 * z)
    y[2] = 0.5102 * z * np.exp(-0.8 * z)
    y[3] = 0.6159 * np.exp(-0.8 * z)
    y[4] = 1.0 - np.exp(-0.8 * z)
    y[6] = 1.0
    return y


def ode(model):
    def rhs(z, y):
        H, Fp, F, Gp, G, Tp, T = y
        radial_acceleration = F**2 + H * Fp - (1.0 - G) ** 2
        azimuthal_acceleration = 2.0 * F * G + H * Gp - 2.0 * F
        theta = T - 1.0

        if model == "traditional":
            Fpp = radial_acceleration + theta
            Gpp = azimuthal_acceleration
        elif model == "lopez":
            rho_linear = 1.0 - theta
            Fpp = rho_linear * radial_acceleration
            Gpp = rho_linear * azimuthal_acceleration
        else:
            raise ValueError(model)

        return np.vstack(
            (-2.0 * F, Fpp, Fp, Gpp, Gp, PR * H * Tp, Tp)
        )

    return rhs


def bc_fixed_tw(tw):
    def bc(ya, yb):
        return np.array(
            (ya[0], ya[2], ya[4], ya[6] - tw,
             yb[2], yb[4] - 1.0, yb[6] - 1.0)
        )

    return bc


def solve_isothermal(model, z, tol=2.0e-9):
    sol = solve_bvp(
        ode(model), bc_fixed_tw(1.0), z, initial_guess(z),
        tol=tol, max_nodes=MAX_NODES,
    )
    if not sol.success:
        raise RuntimeError(f"{model} isothermal solve failed: {sol.message}")
    return sol


def continue_to_targets(
        model, targets, zmax=40.0, step=0.0025, tol=1.0e-8,
        verbose=True):
    targets = sorted(set(float(x) for x in targets))
    if targets[0] < 1.0:
        raise ValueError("This comparison traces the heated branch from Tw=1.")
    n = max(501, int(round(20.0 * zmax)) + 1)
    z = np.linspace(0.0, zmax, n)
    sol = solve_isothermal(model, z)
    current_tw = 1.0
    solutions = {1.0: sol}

    for target in targets:
        if target == 1.0:
            continue
        nsteps = int(np.ceil((target - current_tw) / step))
        for tw in np.linspace(current_tw, target, nsteps + 1)[1:]:
            trial = solve_bvp(
                ode(model), bc_fixed_tw(float(tw)), z, sol.sol(z),
                tol=tol, max_nodes=MAX_NODES,
            )
            if not trial.success:
                raise RuntimeError(
                    f"{model} continuation failed at Tw={tw:.9f}: {trial.message}"
                )
            sol = trial
        current_tw = target
        solutions[target] = sol
        if verbose:
            print(
                f"{model:11s} zmax={zmax:4.0f} Tw={target:.4f} "
                f"Hinf={sol.sol(zmax)[0]: .9f}",
                flush=True,
            )
    return z, solutions


def metrics(model, tw, sol, zmax):
    z = np.linspace(0.0, zmax, 4001)
    H, Fp, F, Gp, G, Tp, T = sol.sol(z)
    imax = int(np.argmax(F))
    if abs(tw - 1.0) > 1.0e-12:
        thermal_thickness = np.trapz((T - 1.0) / (tw - 1.0), z)
    else:
        thermal_thickness = np.nan
    return {
        "model": model,
        "Tw": tw,
        "Hinf": H[-1],
        "Fmax": F[imax],
        "z_Fmax": z[imax],
        "Fp0": Fp[0],
        "Gp0": Gp[0],
        "Tp0": Tp[0],
        "thermal_thickness": thermal_thickness,
        "rho_wall_linear": 2.0 - tw,
        "nodes": sol.x.size,
    }


def write_dict_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def profile_difference(tw, trad, lopez, zmax):
    z = np.linspace(0.0, zmax, 4001)
    yt = trad.sol(z)
    yl = lopez.sol(z)
    row = {"Tw": tw}
    for name, index in (("H", 0), ("F", 2), ("G", 4), ("T", 6)):
        scale = np.sqrt(np.trapz(yt[index] ** 2, z))
        error = np.sqrt(np.trapz((yl[index] - yt[index]) ** 2, z))
        row[f"relative_L2_{name}"] = error / max(scale, 1.0e-14)
        row[f"max_abs_{name}"] = np.max(np.abs(yl[index] - yt[index]))
    return row


def save_profiles(model, solutions, zmax):
    rows = []
    z = np.linspace(0.0, zmax, 2001)
    for tw, sol in solutions.items():
        y = sol.sol(z)
        for i, zz in enumerate(z):
            rows.append((model, tw, zz, y[0, i], y[2, i], y[4, i], y[6, i],
                         y[1, i], y[3, i], y[5, i]))
    with (OUT / f"profiles_{model}.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("model", "Tw", "z", "H", "F", "G", "T", "Fp", "Gp", "Tp"))
        writer.writerows(rows)


def make_figures(traditional, lopez, metrics_rows, differences, zmax):
    metric_lookup = {}
    for row in metrics_rows:
        metric_lookup.setdefault(row["model"], []).append(row)

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6))
    for model, color, label in (
        ("traditional", "#d1495b", "traditional centrifugal Boussinesq"),
        ("lopez", "#00777f", "Lopez generalized Boussinesq"),
    ):
        rows = metric_lookup[model]
        tw = np.array([r["Tw"] for r in rows])
        for ax, key, ylabel in zip(
            axes.flat,
            ("Hinf", "Fmax", "Fp0", "thermal_thickness"),
            (r"$H_\infty$", r"$F_{max}$", r"$F'(0)$", r"$\delta_T$"),
        ):
            ax.plot(tw, [r[key] for r in rows], "o-", ms=3.5, lw=1.8,
                    color=color, label=label)
            ax.set_xlabel(r"$T_w$")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.25)
    for ax in axes.flat:
        ax.axvline(1.048021731, color="0.3", ls="--", lw=1.0)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle("Traditional and Lopez-type Boussinesq base flows")
    fig.tight_layout()
    fig.savefig(OUT / "baseflow_metrics_comparison.png", dpi=210)
    fig.savefig(OUT / "baseflow_metrics_comparison.pdf")
    plt.close(fig)

    common = [1.02, 1.04, 1.048]
    z = np.linspace(0.0, min(15.0, zmax), 1601)
    fig, axes = plt.subplots(len(common), 4, figsize=(12.4, 8.3), sharex=True)
    for i, tw in enumerate(common):
        for model, sols, color, style in (
            ("traditional", traditional, "#d1495b", "--"),
            ("lopez", lopez, "#00777f", "-"),
        ):
            y = sols[tw].sol(z)
            for ax, index, name in zip(axes[i], (2, 4, 0, 6), ("F", "G", "H", "T")):
                ax.plot(z, y[index], style, lw=1.8, color=color, label=model)
                if i == 0:
                    ax.set_title(name)
                ax.grid(alpha=0.22)
        axes[i, 0].set_ylabel(rf"$T_w={tw:g}$")
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("z")
    fig.suptitle("Profile comparison below the traditional fold")
    fig.tight_layout()
    fig.savefig(OUT / "profiles_common_temperatures.png", dpi=210)
    fig.savefig(OUT / "profiles_common_temperatures.pdf")
    plt.close(fig)

    extended = [1.0, 1.05, 1.10, 1.20, 1.50, 1.80]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6), sharex=True)
    colors = plt.cm.plasma(np.linspace(0.05, 0.88, len(extended)))
    z = np.linspace(0.0, min(20.0, zmax), 2001)
    for color, tw in zip(colors, extended):
        y = lopez[tw].sol(z)
        for ax, index, name in zip(axes.flat, (2, 4, 0, 6), ("F", "G", "H", "T")):
            ax.plot(z, y[index], lw=1.7, color=color, label=rf"$T_w={tw:g}$")
            ax.set_ylabel(name)
            ax.grid(alpha=0.22)
    for ax in axes[-1]:
        ax.set_xlabel("z")
    axes[0, 0].legend(ncol=2, fontsize=8)
    fig.suptitle("Formal continuation of the Lopez model beyond the traditional fold")
    fig.tight_layout()
    fig.savefig(OUT / "lopez_extended_profiles.png", dpi=210)
    fig.savefig(OUT / "lopez_extended_profiles.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    tw = [r["Tw"] for r in differences]
    for key, label in (("relative_L2_F", "F"), ("relative_L2_G", "G"),
                       ("relative_L2_H", "H"), ("relative_L2_T", "T")):
        ax.semilogy(tw, [r[key] for r in differences], "o-", lw=1.7, label=label)
    ax.set_xlabel(r"$T_w$")
    ax.set_ylabel("relative L2 difference")
    ax.set_title("Lopez versus traditional profiles")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "relative_profile_differences.png", dpi=210)
    fig.savefig(OUT / "relative_profile_differences.pdf")
    plt.close(fig)


def write_report(metrics_rows, differences, domain_rows):
    by_model = {}
    for row in metrics_rows:
        by_model.setdefault(row["model"], {})[row["Tw"]] = row

    with (OUT / "report.md").open("w", encoding="utf-8") as stream:
        stream.write("# Lopez广义Boussinesq与传统离心Boussinesq基本流对比\n\n")
        stream.write("## 1. 模型定义\n\n")
        stream.write("令 $A=F^2+HF'-(1-G)^2$，$B=2FG+HG'-2F$。两种模型共有\n\n")
        stream.write("$$H'=-2F,\\qquad T''=PrHT',\\qquad Pr=0.72.$$\n\n")
        stream.write("传统固定参考离心浮力模型为\n\n")
        stream.write("$$F''=A+(T-1),\\qquad G''=B.$$\n\n")
        stream.write("Lopez惯性系广义模型在线性密度闭合 $\\rho/\\rho_0=1-(T-1)$ 下为\n\n")
        stream.write("$$F''=[1-(T-1)]A,\\qquad G''=[1-(T-1)]B.$$\n\n")
        stream.write("后者使密度扰动与完整局部惯性加速度耦合，而不是只乘固定的圆盘离心加速度。\n\n")

        stream.write("在圆盘固连旋转系中，完整惯性加速度由相对对流、Coriolis和参考离心项组成。")
        stream.write("传统模型只在最后一项保留密度扰动；Lopez模型等价于在惯性系中保留")
        stream.write("$-\\rho'(\\boldsymbol u\\cdot\\nabla)\\boldsymbol u$，转换到相似变量后使径向和周向加速度整体乘以")
        stream.write("$\\rho/\\rho_0=1-(T-1)$。因此，上述广义方程不是通过经验替换温度源项得到的。\n\n")

        stream.write("## 2. 数值方法\n\n")
        stream.write("采用自适应配置点边值法，从 $T_w=1$ 的经典von Karman解出发，以固定壁温小步延拓。")
        stream.write("主比较使用 $z_{max}=40$，并在 $z_{max}=20,40,60$ 上重复Lopez计算。")
        stream.write("边界条件为 $F=G=H=0,T=T_w$（壁面）以及 $F=0,G=1,T=1$（远场）。")
        stream.write("传统模型延拓至其收敛折点附近，Lopez模型形式延拓至 $T_w=1.8$。\n\n")

        stream.write("## 3. 主要数值结果\n\n")
        stream.write("| model | Tw | Hinf | Fmax | z(Fmax) | Fp0 | Gp0 | thermal thickness | rho_w/rho_inf (linear) |\n")
        stream.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        selected = {
            "traditional": [1.0, 1.02, 1.04, 1.048],
            "lopez": [1.0, 1.02, 1.04, 1.048, 1.05, 1.10, 1.20, 1.50, 1.80],
        }
        for model, values in selected.items():
            for tw in values:
                r = by_model[model][tw]
                thickness = "--" if not np.isfinite(r["thermal_thickness"]) else f"{r['thermal_thickness']:.8f}"
                stream.write(
                    f"| {model} | {tw:.3f} | {r['Hinf']:.8f} | {r['Fmax']:.8f} | "
                    f"{r['z_Fmax']:.5f} | {r['Fp0']:.8f} | {r['Gp0']:.8f} | "
                    f"{thickness} | {r['rho_wall_linear']:.4f} |\n"
                )

        stream.write("\n传统模型连接经典解的下分支终止于先前域收敛得到的\n\n")
        stream.write("$$T_{w,c}=1.048021731,\\qquad H_{\\infty,c}=-0.532761653.$$\n\n")
        stream.write("Lopez模型用固定壁温延拓连续通过该温度，并在形式上计算到 $T_w=1.8$。因此，")
        stream.write("$1.049$ 限制不是所有不可压缩热旋转模型的共同限制，而是传统固定参考离心闭合的分支特征。\n\n")

        stream.write("## 4. 两模型在共同温度范围内的差异\n\n")
        stream.write("| Tw | relL2(F) | relL2(G) | relL2(H) | relL2(T) |\n")
        stream.write("|---:|---:|---:|---:|---:|\n")
        for r in differences:
            stream.write(
                f"| {r['Tw']:.3f} | {r['relative_L2_F']:.6e} | "
                f"{r['relative_L2_G']:.6e} | {r['relative_L2_H']:.6e} | "
                f"{r['relative_L2_T']:.6e} |\n"
            )

        stream.write("\n在 $T_w=1.048$ 时，温度剖面的相对差异仍小于 $10^{-3}$，但 $H$ 和 $F$ 的相对差异分别约为")
        stream.write("56%和25%。传统模型中的独立 $T-1$ 源项使径向流迅速衰减，继而通过 $H'=-2F$ 削弱抽吸并加厚热尾，")
        stream.write("最终构成鞍结正反馈。Lopez模型没有该独立源项，密度扰动只修正实际局部惯性，因此基本流随壁温平滑变化。\n\n")

        stream.write("## 5. Lopez模型的远场域长检验\n\n")
        stream.write("| zmax | Tw | Hinf | Fmax | Fp0 | Gp0 |\n|---:|---:|---:|---:|---:|---:|\n")
        for r in domain_rows:
            stream.write(
                f"| {r['zmax']:.0f} | {r['Tw']:.2f} | {r['Hinf']:.9f} | "
                f"{r['Fmax']:.9f} | {r['Fp0']:.9f} | {r['Gp0']:.9f} |\n"
            )

        stream.write("\n## 6. 适用性判断\n\n")
        stream.write("1. Lopez模型仍采用 $\\nabla\\cdot u=0$、线性密度扰动和常物性，因此不是完全可压缩模型。\n")
        stream.write("2. 控制小参数仍是 $\\epsilon=\\beta|T_w-T_\\infty|$。当前尺度下理想气体有 $\\epsilon=|T_w-1|$。\n")
        stream.write("3. $T_w=1.05$ 附近 $\\epsilon\\approx0.05$，适合用来诊断传统离心闭合；")
        stream.write("$T_w=1.5$ 或 1.8 的结果只是方程的形式延拓，不能作为定量物理预测。\n")
        stream.write("4. 对理想气体，线性密度相对精确关系 $1/T$ 的误差为 $-(T-1)^2$：")
        stream.write("$T_w=1.05,1.1,1.2,1.5,1.8$ 时误差量级分别约为0.25%、1%、4%、25%、64%。")
        stream.write("线性密度在 $T_w=2$ 甚至给出壁面密度为零。\n")
        stream.write("5. 后续低马赫数渐近研究指出，对一般非保守对流加速度，Lopez的对流浮力项可能具有 $O(1)$ 相对误差。")
        stream.write("因此它是有价值的诊断闭合，但不能未经完全可压缩验证就称为更准确模型。\n")
        stream.write("6. 本报告只证明基本流分支的存在性和远场收敛，尚未证明Lopez分支的三维稳定性。\n")
        stream.write("7. 是否真实消除鞍结，应继续与低马赫数/完全可压缩、Sutherland可变物性模型在统一物理尺度下比较。\n\n")

        stream.write("## 7. 结论\n\n")
        stream.write("Lopez广义模型突破了传统模型的 $T_w\\simeq1.049$ 数学折点。")
        stream.write("这表明该折点主要与固定参考离心浮力闭合有关，而不能直接解释为Boussinesq小温差展开在4.8%温差处必然失效。")
        stream.write("不过，Lopez模型只扩展了密度与惯性/离心效应的耦合，其高温结果仍需完全可压缩模型验证。\n")

        stream.write("\n## 参考文献\n\n")
        stream.write("1. Lopez, J. M., Marques, F. & Avila, M. (2013), *The Boussinesq approximation in rapidly rotating flows*, JFM 737, 56-77, DOI: 10.1017/jfm.2013.558.\n")
        stream.write("2. Kang et al. (2021), *Perturbation analysis of baroclinic torque in low-Mach-number flows*, JFM, DOI: 10.1017/jfm.2021.896.\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    common = [1.0, 1.01, 1.02, 1.03, 1.04, 1.045, 1.048]
    generalized = common + [1.05, 1.075, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80]
    zmax = 40.0

    _, traditional = continue_to_targets("traditional", common, zmax=zmax)
    _, lopez = continue_to_targets("lopez", generalized, zmax=zmax)

    metrics_rows = []
    for model, solutions in (("traditional", traditional), ("lopez", lopez)):
        for tw, sol in solutions.items():
            metrics_rows.append(metrics(model, tw, sol, zmax))
    write_dict_csv(OUT / "baseflow_metrics.csv", metrics_rows)

    differences = [
        profile_difference(tw, traditional[tw], lopez[tw], zmax)
        for tw in common
    ]
    write_dict_csv(OUT / "profile_differences.csv", differences)
    save_profiles("traditional", traditional, zmax)
    save_profiles("lopez", lopez, zmax)

    domain_rows = []
    domain_targets = [1.05, 1.10, 1.20, 1.50]
    for length in (20.0, 40.0, 60.0):
        _, solutions = continue_to_targets(
            "lopez", domain_targets, zmax=length, step=0.005, tol=2.0e-8
        )
        for tw, sol in solutions.items():
            if tw == 1.0:
                continue
            row = metrics("lopez", tw, sol, length)
            row["zmax"] = length
            domain_rows.append(row)
    write_dict_csv(OUT / "lopez_domain_convergence.csv", domain_rows)

    make_figures(traditional, lopez, metrics_rows, differences, zmax)
    write_report(metrics_rows, differences, domain_rows)
    print(f"Results written to {OUT.resolve()}")


if __name__ == "__main__":
    main()
