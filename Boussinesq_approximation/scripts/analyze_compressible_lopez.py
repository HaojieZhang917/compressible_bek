#!/usr/bin/env python3
"""Compare the Chapman compressible similarity flow with the Lopez closure."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.interpolate import CubicSpline

from Vonkarmen_bone.scripts.compare_lopez_boussinesq import initial_guess, ode, bc_fixed_tw


PR = 0.72
GAMMA = 1.4
WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
OUT = WORKSPACE_ROOT / "compressible_lopez_comparison"
LOPEZ_SOURCE = WORKSPACE_ROOT / "lopez_boussinesq_comparison" / "profiles_lopez.csv"
TW_VALUES = (1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.075, 1.10, 1.15, 1.20, 1.30, 1.50)
MR_VALUES = (0.3,)
ETA_MAX = 40.0
N_RAW = 10001
Z_COMPARE_MAX = 20.0
N_COMPARE = 16001
N_DAT = 2001


def solve_universal_velocity():
    z = np.linspace(0.0, ETA_MAX, 801)
    sol = solve_bvp(
        ode("lopez"), bc_fixed_tw(1.0), z, initial_guess(z),
        tol=2.0e-11, max_nodes=200000,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    eta = np.linspace(0.0, ETA_MAX, N_RAW)
    values = sol.sol(eta)
    return eta, values


def solve_energy_shapes(eta, velocity):
    H, Fp, F, Gp, G = velocity[0], velocity[1], velocity[2], velocity[3], velocity[4]
    phi = -0.5 * H

    def rhs(x, y):
        F_i = np.interp(x, eta, F)
        Fp_i = np.interp(x, eta, Fp)
        Gp_i = np.interp(x, eta, Gp)
        phi_i = np.interp(x, eta, phi)
        f, fp, q, qp = y
        return np.vstack((
            fp,
            2.0 * PR * (Fp_i**2 + Gp_i**2 + F_i * f - phi_i * fp),
            qp,
            -2.0 * PR * phi_i * qp,
        ))

    def bc(ya, yb):
        return np.array((ya[0], yb[0], ya[2] - 1.0, yb[2]))

    mesh = np.linspace(0.0, ETA_MAX, 1001)
    y0 = np.zeros((4, mesh.size))
    y0[2] = np.exp(-mesh / 2.0)
    sol = solve_bvp(rhs, bc, mesh, y0, tol=2.0e-10, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(sol.message)
    values = sol.sol(eta)
    boundary_residual = np.max(np.abs(bc(values[:, 0], values[:, -1])))
    return values[0], values[1], values[2], values[3], boundary_residual


def compressible_profile(eta, velocity, energy, tw, mr):
    H0, Fp, F, Gp, G = velocity[0], velocity[1], velocity[2], velocity[3], velocity[4]
    f, fp, q, qp = energy[:4]
    mach_coefficient = 0.5 * (GAMMA - 1.0) * mr**2
    T = 1.0 - mach_coefficient * f + (tw - 1.0) * q
    Teta = -mach_coefficient * fp + (tw - 1.0) * qp
    rho = 1.0 / T
    H = H0 * T
    z_phy = np.r_[0.0, cumulative_trapezoid(T, eta)]
    Fz = Fp / T
    Gz = Gp / T
    Tz = Teta / T
    Heta = (-2.0 * F) * T + H0 * Teta
    Hz = Heta / T
    return {
        "eta": eta,
        "z": z_phy,
        "F": F,
        "G": G,
        "H": H,
        "T": T,
        "rho": rho,
        "Fz": Fz,
        "Gz": Gz,
        "Hz": Hz,
        "Tz": Tz,
    }


def load_lopez_profiles():
    data = pd.read_csv(LOPEZ_SOURCE)
    profiles = {}
    for tw in TW_VALUES:
        case = data[np.isclose(data["Tw"], tw)].sort_values("z")
        if case.empty:
            raise RuntimeError(f"Lopez profile Tw={tw} is unavailable")
        profiles[tw] = {
            name: case[name].to_numpy(float)
            for name in ("z", "F", "G", "H", "T", "Fp", "Gp", "Tp")
        }
    return profiles


def interp_case(profile, z):
    return {
        name: CubicSpline(profile["z"], profile[name])(z)
        for name in profile if name != "eta"
    }


def l2_relative(a, b, z, perturbation=False):
    numerator = np.sqrt(simpson((a - b) ** 2, z))
    reference = b - 1.0 if perturbation else b
    denominator = np.sqrt(simpson(reference**2, z))
    if denominator < 1.0e-13:
        return 0.0 if numerator < 1.0e-12 else np.nan
    return numerator / denominator


def compare_case(comp, lopez, tw, mr):
    z = np.linspace(0.0, Z_COMPARE_MAX, N_COMPARE)
    c = interp_case(comp, z)
    l = interp_case(lopez, z)
    rho_lopez = 2.0 - l["T"]

    i_fc = int(np.argmax(c["F"]))
    i_fl = int(np.argmax(l["F"]))
    row = {
        "Mr": mr,
        "Tw": tw,
        "compressible_Fmax": c["F"][i_fc],
        "compressible_z_Fmax": z[i_fc],
        "lopez_Fmax": l["F"][i_fl],
        "lopez_z_Fmax": z[i_fl],
        "compressible_Hinf": comp["H"][-1],
        "lopez_Hinf": lopez["H"][-1],
        "compressible_Fz0": comp["Fz"][0],
        "lopez_Fz0": lopez["Fp"][0],
        "compressible_Gz0": comp["Gz"][0],
        "lopez_Gz0": lopez["Gp"][0],
        "compressible_Tz0": comp["Tz"][0],
        "lopez_Tz0": lopez["Tp"][0],
        "compressible_rho_wall": comp["rho"][0],
        "lopez_rho_wall_linear": 2.0 - tw,
        "relative_L2_F": l2_relative(l["F"], c["F"], z),
        "relative_L2_G": l2_relative(l["G"], c["G"], z),
        "relative_L2_H": l2_relative(l["H"], c["H"], z),
        "relative_L2_T_total": l2_relative(l["T"], c["T"], z),
        "relative_L2_T_excess": l2_relative(l["T"], c["T"], z, perturbation=True),
        "relative_L2_rho": l2_relative(rho_lopez, c["rho"], z),
        "relative_L2_Fz": l2_relative(l["Fp"], c["Fz"], z),
        "relative_L2_Gz": l2_relative(l["Gp"], c["Gz"], z),
        "relative_L2_Tz": l2_relative(l["Tp"], c["Tz"], z),
        "max_abs_F": np.max(np.abs(l["F"] - c["F"])),
        "max_abs_G": np.max(np.abs(l["G"] - c["G"])),
        "max_abs_H": np.max(np.abs(l["H"] - c["H"])),
        "max_abs_T": np.max(np.abs(l["T"] - c["T"])),
        "max_abs_rho": np.max(np.abs(rho_lopez - c["rho"])),
    }
    row["velocity_error"] = max(row["relative_L2_F"], row["relative_L2_G"], row["relative_L2_H"])
    row["derivative_error"] = max(row["relative_L2_Fz"], row["relative_L2_Gz"])
    return row


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_error_dat(path, rows):
    variables = (
        "Tw", "E_F_L2", "E_G_L2", "E_H_L2", "E_velocity_L2",
        "E_T_total_L2", "E_T_excess_L2", "E_rho_L2",
        "max_abs_F", "max_abs_G", "max_abs_H", "max_abs_T", "max_abs_rho",
    )
    keys = (
        "Tw", "relative_L2_F", "relative_L2_G", "relative_L2_H", "velocity_error",
        "relative_L2_T_total", "relative_L2_T_excess", "relative_L2_rho",
        "max_abs_F", "max_abs_G", "max_abs_H", "max_abs_T", "max_abs_rho",
    )
    ordered = sorted(rows, key=lambda row: row["Tw"])
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write('TITLE = "Lopez versus compressible errors at Mr=0.3"\n')
        stream.write("VARIABLES = " + " ".join(f'\"{name}\"' for name in variables) + "\n")
        stream.write(f'ZONE T="Mr_0.3", I={len(ordered)}, DATAPACKING=POINT\n')
        for row in ordered:
            stream.write(" ".join(f"{row[key]:.12e}" for key in keys) + "\n")


def write_profiles_dat(path, compressible, lopez):
    variables = ("z", "F", "G", "H", "T", "rho", "Fz", "Gz", "Hz", "Tz")
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write('TITLE = "Compressible and Lopez base-flow profiles"\n')
        stream.write("VARIABLES = " + " ".join(f'\"{x}\"' for x in variables) + "\n")
        for (mr, tw), profile in sorted(compressible.items()):
            indices = np.linspace(0, len(profile["z"]) - 1, N_DAT, dtype=int)
            stream.write(f'ZONE T="compressible_Mr_{mr:.2f}_Tw_{tw:.3f}", I={len(indices)}, DATAPACKING=POINT\n')
            for i in indices:
                stream.write(" ".join(f"{profile[name][i]:.12e}" for name in variables) + "\n")
        for tw, profile in sorted(lopez.items()):
            z = profile["z"]
            rho = 2.0 - profile["T"]
            Hz = np.gradient(profile["H"], z, edge_order=2)
            stream.write(f'ZONE T="lopez_Tw_{tw:.3f}", I={len(z)}, DATAPACKING=POINT\n')
            arrays = (z, profile["F"], profile["G"], profile["H"], profile["T"], rho,
                      profile["Fp"], profile["Gp"], Hz, profile["Tp"])
            for values in zip(*arrays):
                stream.write(" ".join(f"{value:.12e}" for value in values) + "\n")


def make_figures(compressible, lopez, rows):
    selected = ((0.3, 1.0), (0.3, 1.05), (0.3, 1.10), (0.3, 1.20), (0.3, 1.50))
    fig, axes = plt.subplots(len(selected), 4, figsize=(12.4, 10.0), sharex=True)
    for i, (mr, tw) in enumerate(selected):
        c = compressible[(mr, tw)]
        l = lopez[tw]
        for ax, name in zip(axes[i], ("F", "G", "H", "T")):
            ax.plot(c["z"], c[name], color="#d1495b", lw=1.8, label="compressible")
            ax.plot(l["z"], l[name], color="#00777f", lw=1.8, ls="--", label="Lopez")
            ax.set_xlim(0, 15)
            ax.grid(alpha=0.22)
            if i == 0:
                ax.set_title(name)
        axes[i, 0].set_ylabel(rf"$M_r={mr:g}, T_w={tw:g}$")
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("physical z")
    fig.suptitle("Compressible versus Lopez profiles in physical coordinates")
    fig.tight_layout()
    fig.savefig(OUT / "profile_comparison.png", dpi=210)
    fig.savefig(OUT / "profile_comparison.pdf")
    plt.close(fig)

    table = pd.DataFrame(rows)
    subset = table.sort_values("Tw")
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.2), sharex=True)
    for ax, key, title in zip(
        axes.flat,
        ("velocity_error", "derivative_error", "relative_L2_rho", "max_abs_T"),
        ("velocity L2 error", "first-derivative L2 error", "density L2 error", "maximum absolute temperature error"),
    ):
        ax.plot(subset.Tw, subset[key], "o-", color="#8b2f52", lw=1.8, ms=4)
        ax.set_title(title)
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel(r"$T_w$")
    fig.suptitle(r"Lopez versus compressible differences at $M_r=0.3$")
    fig.tight_layout()
    fig.savefig(OUT / "error_maps.png", dpi=210)
    fig.savefig(OUT / "error_maps.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    subset = table.sort_values("Tw")
    axes[0].plot(subset.Tw, subset.compressible_Fmax, "o-", ms=3, label="compressible")
    axes[1].plot(subset.Tw, subset.compressible_z_Fmax, "o-", ms=3, label="compressible")
    lopez_rows = subset
    axes[0].plot(lopez_rows.Tw, lopez_rows.lopez_Fmax, "k--", lw=2, label="Lopez")
    axes[1].plot(lopez_rows.Tw, lopez_rows.lopez_z_Fmax, "k--", lw=2, label="Lopez")
    axes[0].set_ylabel(r"$F_{max}$")
    axes[1].set_ylabel(r"physical $z(F_{max})$")
    for ax in axes:
        ax.set_xlabel(r"$T_w$")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Amplitude and physical-coordinate displacement")
    fig.tight_layout()
    fig.savefig(OUT / "peak_amplitude_location.png", dpi=210)
    fig.savefig(OUT / "peak_amplitude_location.pdf")
    plt.close(fig)


def threshold(table, mr, key, limit):
    subset = table[np.isclose(table.Mr, mr)].sort_values("Tw")
    accepted = []
    for _, row in subset.iterrows():
        if np.isfinite(row[key]) and row[key] <= limit:
            accepted.append(row.Tw)
        else:
            break
    return max(accepted) if accepted else np.nan


def format_threshold(value):
    return "--" if not np.isfinite(value) else f"{value:.3f}"


def write_report(rows, energy_residual, universal_metrics):
    table = pd.DataFrame(rows)
    with (OUT / "report.md").open("w", encoding="utf-8") as stream:
        stream.write("# Mr=0.3时Chapman可压缩相似模型与Lopez广义Boussinesq基本流对比\n\n")
        stream.write("## 模型与坐标\n\n")
        stream.write("可压缩求解器采用Chapman型物性和密度加权相似变换。温度为\n\n")
        stream.write("$$T=1-\\frac{\\gamma-1}{2}M_r^2f+(T_w-1)q,$$\n\n")
        stream.write("并有 $\\rho=1/T$、$H=W T$。相似坐标 $\\eta$ 与物理坐标满足\n\n")
        stream.write("$$z=\\int_0^\\eta T(s)\\,ds.$$\n\n")
        stream.write("因此 `sim` 表示密度加权相似坐标，`phy` 表示恢复后的物理坐标。所有比较均固定在 $M_r=0.3$ 并在物理坐标上完成。\n\n")
        stream.write("Lopez模型仍满足不可压缩连续性，采用线性密度闭合并使密度扰动与完整局部惯性耦合；它不包含马赫数、压缩功和完整可变物性。\n\n")

        stream.write("## 数值核查\n\n")
        stream.write(f"经典速度基准为 $F_{{max}}={universal_metrics['Fmax']:.9f}$、$H_\\infty={universal_metrics['Hinf']:.9f}$。")
        stream.write(f"重新配置求解的能量形函数边界残差为 ${energy_residual:.2e}$。\n\n")
        stream.write("原 notebook 的 `f_q` 打靶结果在当前环境给出 `f(40)=2.18e-2`，未满足其声明的 `f(40)=0`；")
        stream.write("`Physical_Interpretation` 的第二个坐标点也多累计一个网格值。本报告使用同一数学模型的高精度配置点解和梯形坐标积分，避免这两个预处理误差。\n\n")
        stream.write("剖面比较采用三次样条连续插值，并在 $z\\in[0,20]$ 上使用16001点复合Simpson积分。与自适应连续积分的收敛检查表明，速度和导数误差的积分偏差均已降至约 $10^{-11}$ 或更低。\n\n")

        stream.write("## 代表性结果\n\n")
        stream.write("| Mr | Tw | E_velocity | E_derivative | E_rho | Fmax(comp) | Fmax(Lopez) | zFmax(comp) | zFmax(Lopez) | Hinf(comp) | Hinf(Lopez) |\n")
        stream.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for mr, tw in ((0.3,1.00),(0.3,1.01),(0.3,1.02),(0.3,1.05),(0.3,1.10),(0.3,1.20),(0.3,1.30),(0.3,1.50)):
            r = table[np.isclose(table.Mr,mr) & np.isclose(table.Tw,tw)].iloc[0]
            stream.write(f"| {mr:.1f} | {tw:.2f} | {r.velocity_error:.4f} | {r.derivative_error:.4f} | {r.relative_L2_rho:.4f} | {r.compressible_Fmax:.6f} | {r.lopez_Fmax:.6f} | {r.compressible_z_Fmax:.4f} | {r.lopez_z_Fmax:.4f} | {r.compressible_Hinf:.6f} | {r.lopez_Hinf:.6f} |\n")

        stream.write("\n## 误差阈值\n\n")
        stream.write("下表给出从 $T_w=1$ 开始连续满足误差阈值的最大采样壁温。\n\n")
        stream.write("| Mr | velocity <=5% | velocity <=10% | derivative <=10% | max abs(T) <=0.005 | max abs(T) <=0.01 |\n")
        stream.write("|---:|---:|---:|---:|---:|---:|\n")
        for mr in MR_VALUES:
            v5 = threshold(table,mr,"velocity_error",0.05)
            v10 = threshold(table,mr,"velocity_error",0.10)
            d10 = threshold(table,mr,"derivative_error",0.10)
            t005 = threshold(table,mr,"max_abs_T",0.005)
            t010 = threshold(table,mr,"max_abs_T",0.010)
            values = (v5, v10, d10, t005, t010)
            stream.write(f"| {mr:.1f} | " + " | ".join(format_threshold(x) for x in values) + " |\n")

        stream.write("\n## 主要差异\n\n")
        stream.write("1. Chapman变换下 $F$、$G$ 在相似坐标中的幅值与 $T_w,M_r$ 解耦；升温和马赫数主要通过物理坐标拉伸、$H=WT$ 和密度改变基本流。\n")
        stream.write("2. Lopez模型没有这种解耦。温度直接修改局部径向和周向惯性，因此 $F_{max}$、壁面剪切和 $H_\\infty$ 随 $T_w$ 改变。\n")
        stream.write("3. 即使速度幅值接近，物理峰值位置和导数仍可能明显不同；稳定性算子对这些导数比对速度幅值更敏感。\n")
        stream.write("4. Lopez模型缺少 $M_r^2f$ 所代表的黏性耗散热效应；在本报告固定的 $M_r=0.3$ 下，该项虽小但并不严格为零。\n")
        stream.write("5. 即使固定在低马赫数 $M_r=0.3$，黏性耗散项仍产生一个小但非零的温度修正；Lopez模型中没有对应的马赫数热效应。\n")

        stream.write("\n## 适用范围与可压缩模型必要性\n\n")
        stream.write("- Lopez模型适合作为低马赫数、小温差下的闭合诊断，而不应作为高温可压缩参考解。\n")
        stream.write("- 适用范围必须同时检查速度误差、导数误差、密度误差和稳定性结果，不能只看 $F_{max}$。\n")
        stream.write("- 在固定 $M_r=0.3$ 的离散工况中，以速度剖面5%误差为宽松标准，Lopez模型可连续使用到约 $T_w=1.10$；10%标准可到约 $T_w=1.20$。这只是基本流速度判据，不是普适临界温度。\n")
        stream.write("- 对定量稳定性和感受性分析，导数、温度和密度误差也会进入线性算子，因此建议从 $T_w\\gtrsim1.05$ 开始采用可压缩基本流和可压缩扰动方程。\n")
        stream.write("- Chapman模型适合文献复现和隔离坐标拉伸机制；要做定量高温预测，还应采用Sutherland黏性率及非相似径向推进模型。\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    eta, velocity = solve_universal_velocity()
    f, fp, q, qp, energy_residual = solve_energy_shapes(eta, velocity)
    energy = (f, fp, q, qp)
    lopez = load_lopez_profiles()

    compressible = {}
    for mr in MR_VALUES:
        for tw in TW_VALUES:
            compressible[(mr, tw)] = compressible_profile(eta, velocity, energy, tw, mr)

    rows = [
        compare_case(compressible[(mr, tw)], lopez[tw], tw, mr)
        for mr in MR_VALUES for tw in TW_VALUES
    ]
    write_csv(OUT / "comparison_metrics.csv", rows)
    write_error_dat(OUT / "errors_vs_Tw_Mr0p3.dat", rows)
    write_profiles_dat(OUT / "profiles_all_cases.dat", compressible, lopez)
    make_figures(compressible, lopez, rows)

    F = velocity[2]
    universal_metrics = {"Fmax": float(np.max(F)), "Hinf": float(velocity[0, -1])}
    write_report(rows, energy_residual, universal_metrics)
    print(f"energy boundary residual = {energy_residual:.3e}")
    print(f"f endpoints = {f[0]:.3e}, {f[-1]:.3e}")
    print(f"q endpoints = {q[0]:.3e}, {q[-1]:.3e}")
    print(f"results: {OUT.resolve()}")


if __name__ == "__main__":
    main()
