#!/usr/bin/env python3
"""Compare Chapman-DH and Sutherland marching base-flow results.

Defaults compare:
  Chapman data: baseflow_comparison_data/compressible_physical_raw.csv
  Sutherland data: sutherland_Tw1.5_Mr0.3_R20/
  Case: Tw=1.5, final Sutherland radius r=max, Mr=0.3.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
CHAPMAN_RAW = ROOT / "baseflow_comparison_data" / "compressible_physical_raw.csv"
SUTH_DIR = ROOT / "sutherland_Tw1.5_Mr0.3_R20"
OUT_DIR = ROOT / "chapman_sutherland_comparison"


def env_float(name: str, default: float) -> float:
    return float(os.environ.get(name, str(default)))


def load_csv(path: Path):
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8")


def rows_close(data, column: str, value: float):
    return data[np.isclose(data[column], value, rtol=0.0, atol=1e-10)]


def interp(x, y, xnew):
    order = np.argsort(x)
    xs = np.asarray(x)[order]
    ys = np.asarray(y)[order]
    keep = np.r_[True, np.diff(xs) > 0.0]
    xs = xs[keep]
    ys = ys[keep]
    return np.interp(xnew, xs, ys, left=ys[0], right=ys[-1])


def wall_derivative(z, y):
    return (y[1] - y[0]) / (z[1] - z[0])


def peak_metrics(z, F, G, H, T):
    i_f = int(np.argmax(F))
    i_g = int(np.argmin(G))
    i_h = int(np.argmin(H))
    return {
        "Fmax": float(F[i_f]),
        "z_Fmax": float(z[i_f]),
        "Gmin": float(G[i_g]),
        "z_Gmin": float(z[i_g]),
        "Hmin": float(H[i_h]),
        "z_Hmin": float(z[i_h]),
        "Tmin": float(np.min(T)),
        "Tmax": float(np.max(T)),
        "wall_dF": float(wall_derivative(z, F)),
        "wall_dG": float(wall_derivative(z, G)),
        "wall_dH": float(wall_derivative(z, H)),
        "wall_dT": float(wall_derivative(z, T)),
    }


def write_csv(path: Path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def savefig(fig, stem: str):
    png = OUT_DIR / f"{stem}.png"
    pdf = OUT_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main() -> int:
    tw = env_float("TW", 1.5)
    target_mr = env_float("MR", 0.3)
    zmax_plot = env_float("ZMAX_COMPARE", 10.0)
    zmax_common = env_float("ZMAX_COMMON", 20.0)
    n_common = int(os.environ.get("N_COMMON", "2001"))

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    chap = load_csv(CHAPMAN_RAW)
    chap_tw = rows_close(chap, "Tw", tw)
    if chap_tw.size == 0:
        raise RuntimeError(f"No Chapman data for Tw={tw}")

    suth_profiles = load_csv(SUTH_DIR / "sutherland_marching_profiles.csv")
    suth_summary = load_csv(SUTH_DIR / "sutherland_marching_summary.csv")
    final_r = float(np.max(suth_profiles["r"]))
    suth_final = rows_close(suth_profiles, "r", final_r)
    final_mr = float(suth_final["Mr"][0])

    z = np.linspace(0.0, zmax_common, n_common)
    ch_F = interp(chap_tw["z_physical"], chap_tw["F_c"], z)
    ch_G = interp(chap_tw["z_physical"], chap_tw["G_c"], z)
    ch_H = interp(chap_tw["z_physical"], chap_tw["H_c"], z)
    ch_T = interp(chap_tw["z_physical"], chap_tw["T_c"], z)

    su_F = interp(suth_final["z"], suth_final["U"], z)
    su_G = interp(suth_final["z"], suth_final["V"], z)
    su_H = interp(suth_final["z"], suth_final["W"], z)
    su_T = interp(suth_final["z"], suth_final["T"], z)

    write_csv(
        OUT_DIR / "chapman_sutherland_common_grid.csv",
        [
            "z",
            "F_chapman",
            "G_chapman",
            "H_chapman",
            "T_chapman",
            "F_sutherland",
            "G_sutherland",
            "H_sutherland",
            "T_sutherland",
            "dF_suth_minus_chap",
            "dG_suth_minus_chap",
            "dH_suth_minus_chap",
            "dT_suth_minus_chap",
        ],
        zip(
            z,
            ch_F,
            ch_G,
            ch_H,
            ch_T,
            su_F,
            su_G,
            su_H,
            su_T,
            su_F - ch_F,
            su_G - ch_G,
            su_H - ch_H,
            su_T - ch_T,
        ),
    )

    ch_m = peak_metrics(z, ch_F, ch_G, ch_H, ch_T)
    su_m = peak_metrics(z, su_F, su_G, su_H, su_T)

    metric_header = ["model", "Tw", "Mr", "r", *ch_m.keys()]
    write_csv(
        OUT_DIR / "chapman_sutherland_metrics.csv",
        metric_header,
        [
            ["Chapman-DH", tw, target_mr, "local", *ch_m.values()],
            ["Sutherland-BDF2", tw, final_mr, final_r, *su_m.values()],
        ],
    )

    mask = z <= zmax_plot
    specs = [
        ("F", ch_F, su_F),
        ("G", ch_G, su_G),
        ("H", ch_H, su_H),
        ("T", ch_T, su_T),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True)
    for ax, (name, ch, su) in zip(axes.ravel(), specs):
        ax.plot(z[mask], ch[mask], lw=2.0, label="Chapman-DH")
        ax.plot(z[mask], su[mask], lw=1.8, ls="--", label="Sutherland-BDF2")
        ax.set_title(name)
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.25)
    axes[1, 0].set_xlabel("physical z")
    axes[1, 1].set_xlabel("physical z")
    axes[0, 0].legend(frameon=False)
    fig.suptitle(f"Chapman vs Sutherland profiles, Tw={tw:g}, Mr~{target_mr:g}", y=1.02)
    fig.tight_layout()
    savefig(fig, "profiles_chapman_vs_sutherland")

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True)
    for ax, (name, ch, su) in zip(axes.ravel(), specs):
        ax.plot(z[mask], (su - ch)[mask], lw=2.0)
        ax.axhline(0.0, color="0.45", lw=0.9)
        ax.set_title(f"Sutherland - Chapman: {name}")
        ax.set_ylabel(f"d{name}")
        ax.grid(True, alpha=0.25)
    axes[1, 0].set_xlabel("physical z")
    axes[1, 1].set_xlabel("physical z")
    fig.suptitle("Profile differences on common physical grid", y=1.02)
    fig.tight_layout()
    savefig(fig, "differences_chapman_vs_sutherland")

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.1))
    axes[0].plot(suth_summary["Mr"], suth_summary["Umax"], "o-", label="Sutherland Umax")
    axes[0].axhline(ch_m["Fmax"], color="C1", ls="--", label="Chapman Fmax")
    axes[0].set_xlabel("local Mr")
    axes[0].set_ylabel("radial peak")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].plot(suth_summary["Mr"], suth_summary["Wmin"], "o-", label="Sutherland Wmin")
    axes[1].axhline(ch_m["Hmin"], color="C1", ls="--", label="Chapman Hmin")
    axes[1].set_xlabel("local Mr")
    axes[1].set_ylabel("axial minimum")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("Sutherland radial marching trend compared with Chapman local value", y=1.04)
    fig.tight_layout()
    savefig(fig, "marching_summary_vs_chapman")

    max_diffs = {
        "max_abs_dF": float(np.max(np.abs(su_F - ch_F))),
        "max_abs_dG": float(np.max(np.abs(su_G - ch_G))),
        "max_abs_dH": float(np.max(np.abs(su_H - ch_H))),
        "max_abs_dT": float(np.max(np.abs(su_T - ch_T))),
    }

    report = [
        "# Chapman-DH vs Sutherland-BDF2 Comparison",
        "",
        f"Case: `Tw={tw}`, Chapman local `Mr={target_mr}`, Sutherland final station `r={final_r:g}`, `Mr={final_mr:g}`.",
        "",
        "## Mathematical Difference",
        "",
        "Chapman-DH uses `mu=T`, so with ideal gas `rho=1/T` the product `rho*mu=1`. The Dorodnitsyn-Howarth transformation removes the temperature dependence from the radial and azimuthal basic-flow equations.",
        "",
        "Sutherland-BDF2 uses `mu=T^(3/2)*(1+S/Tinf)/(T+S/Tinf)`, so `rho*mu=mu/T` is not constant. Temperature and viscosity therefore remain in the momentum diffusion terms, and the basic velocity amplitudes can vary with radius and wall temperature.",
        "",
        "## Peak Metrics",
        "",
        "| model | F/U max | z peak | G/V min | H/W min | Tmax | wall dT |",
        "|:---|---:|---:|---:|---:|---:|---:|",
        f"| Chapman-DH | {ch_m['Fmax']:.6e} | {ch_m['z_Fmax']:.4f} | {ch_m['Gmin']:.6e} | {ch_m['Hmin']:.6e} | {ch_m['Tmax']:.6e} | {ch_m['wall_dT']:.6e} |",
        f"| Sutherland-BDF2 | {su_m['Fmax']:.6e} | {su_m['z_Fmax']:.4f} | {su_m['Gmin']:.6e} | {su_m['Hmin']:.6e} | {su_m['Tmax']:.6e} | {su_m['wall_dT']:.6e} |",
        "",
        "## Maximum Profile Differences",
        "",
        "| quantity | max abs difference |",
        "|:---|---:|",
    ]
    for key, val in max_diffs.items():
        report.append(f"| {key} | {val:.6e} |")
    report.extend(
        [
            "",
            "## Generated Files",
            "",
            "* `profiles_chapman_vs_sutherland.png/pdf`",
            "* `differences_chapman_vs_sutherland.png/pdf`",
            "* `marching_summary_vs_chapman.png/pdf`",
            "* `chapman_sutherland_common_grid.csv`",
            "* `chapman_sutherland_metrics.csv`",
        ]
    )
    (OUT_DIR / "chapman_sutherland_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote comparison to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
