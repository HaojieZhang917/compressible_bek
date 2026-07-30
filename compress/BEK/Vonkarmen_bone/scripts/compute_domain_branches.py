#!/usr/bin/env python3
"""Compute finite-domain Boussinesq branches for several far-field limits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

from investigate_boussinesq_fold import solve_fixed_h, solve_isothermal

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent

def trace_branch(zmax: float, hmin: float, hmax: float, dh: float):
    """Trace Tw(Hinf), which stays regular at a Tw turning point."""
    n = max(401, int(round(20.0 * zmax)) + 1)
    z = np.linspace(0.0, zmax, n)
    seed = solve_isothermal(zmax, n, tol=2.0e-9)
    guess = seed.sol(z)
    tw = 1.0

    count = int(round((hmax - hmin) / dh))
    h_values = np.linspace(hmin, hmax, count + 1)
    rows = []
    for i, h_inf in enumerate(h_values):
        sol = solve_fixed_h(float(h_inf), z, guess, tw, tol=1.0e-8)
        tw = float(sol.p[0])
        guess = sol.sol(z)
        wall = sol.sol(0.0)
        rows.append(
            (zmax, h_inf, tw, wall[1], wall[3], wall[5], sol.x.size)
        )
        if i % 25 == 0 or i == h_values.size - 1:
            print(
                f"zmax={zmax:4.1f}: {i + 1:3d}/{h_values.size}, "
                f"Hinf={h_inf:.4f}, Tw={tw:.9f}",
                flush=True,
            )
    return np.asarray(rows, dtype=float)


def locate_turning_points(rows):
    """Locate dTw/dHinf=0 from the uniformly sampled branch."""
    h_inf = rows[:, 1]
    tw = rows[:, 2]
    spline = CubicSpline(h_inf, tw)
    roots = spline.derivative().roots()
    roots = [float(h) for h in roots if h_inf[1] < h < h_inf[-2]]
    roots.sort()
    return [(h, float(spline(h))) for h in roots]


def write_csv(path: Path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--zmax", nargs="+", type=float, default=[15, 20, 25, 30, 40]
    )
    parser.add_argument("--hmin", type=float, default=-0.75)
    parser.add_argument("--hmax", type=float, default=-0.05)
    parser.add_argument("--dh", type=float, default=0.0025)
    parser.add_argument(
        "--out", type=Path,
        default=WORKSPACE_ROOT / "boussinesq_domain_branches",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    all_rows = []
    turning_rows = []
    branches = {}
    for zmax in args.zmax:
        rows = trace_branch(zmax, args.hmin, args.hmax, args.dh)
        branches[zmax] = rows
        all_rows.extend(rows.tolist())
        turns = locate_turning_points(rows)
        for point, (h_inf, tw) in enumerate(turns[:2], start=1):
            turning_rows.append((zmax, point, h_inf, tw))

        tag = f"{zmax:g}".replace(".", "p")
        write_csv(
            args.out / f"Hinf_Tw_zmax_{tag}.csv",
            ["zmax", "Hinf", "Tw", "Fp0", "Gp0", "Tp0", "nodes"],
            rows,
        )

    write_csv(
        args.out / "Hinf_Tw_all_zmax.csv",
        ["zmax", "Hinf", "Tw", "Fp0", "Gp0", "Tp0", "nodes"],
        all_rows,
    )
    write_csv(
        args.out / "turning_points_by_zmax.csv",
        ["zmax", "turning_point", "Hinf", "Tw"],
        turning_rows,
    )

    fig, ax = plt.subplots(figsize=(8.2, 5.7))
    colors = plt.cm.viridis(np.linspace(0.08, 0.88, len(branches)))
    for color, (zmax, rows) in zip(colors, branches.items()):
        ax.plot(rows[:, 1], rows[:, 2], lw=2.0, color=color, label=rf"$z_{{max}}={zmax:g}$")
        turns = [(r[2], r[3]) for r in turning_rows if r[0] == zmax]
        if turns:
            ax.scatter(
                [h for h, _ in turns], [tw for _, tw in turns],
                s=38, color=color, edgecolor="black", linewidth=0.6, zorder=4,
            )
    ax.set_xlabel(r"$H_\infty$")
    ax.set_ylabel(r"$T_w$")
    ax.set_title(r"Finite-domain Boussinesq branches: $T_w(H_\infty)$")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(args.out / "Hinf_Tw_curves.png", dpi=220)
    fig.savefig(args.out / "Hinf_Tw_curves.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.7))
    for color, (zmax, rows) in zip(colors, branches.items()):
        ax.plot(rows[:, 1], rows[:, 2], lw=2.0, color=color, label=rf"$z_{{max}}={zmax:g}$")
        turns = [(r[2], r[3]) for r in turning_rows if r[0] == zmax]
        ax.scatter(
            [h for h, _ in turns], [tw for _, tw in turns],
            s=42, color=color, edgecolor="black", linewidth=0.6, zorder=4,
        )
    ax.set_xlim(-0.58, -0.10)
    ax.set_ylim(1.025, 1.057)
    ax.set_xlabel(r"$H_\infty$")
    ax.set_ylabel(r"$T_w$")
    ax.set_title("Turning-point region")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(args.out / "Hinf_Tw_turning_region.png", dpi=220)
    fig.savefig(args.out / "Hinf_Tw_turning_region.pdf")
    plt.close(fig)

    with (args.out / "summary.md").open("w", encoding="utf-8") as stream:
        stream.write("# Domain-dependent Boussinesq turning points\n\n")
        stream.write(
            f"Branches were parameterized by Hinf on [{args.hmin}, {args.hmax}] "
            f"with dHinf={args.dh}. Turning points satisfy dTw/dHinf=0 "
            "on a cubic spline through the computed branch.\n\n"
        )
        stream.write("| zmax | point | Hinf | Tw |\n|---:|---:|---:|---:|\n")
        for zmax, point, h_inf, tw in turning_rows:
            stream.write(f"| {zmax:g} | {point} | {h_inf:.9f} | {tw:.9f} |\n")

    print(f"Wrote results to {args.out.resolve()}")


if __name__ == "__main__":
    main()
