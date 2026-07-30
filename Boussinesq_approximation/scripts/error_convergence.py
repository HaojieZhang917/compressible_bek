#!/usr/bin/env python3
"""Check convergence of profile-error integration at fixed Mr=0.3."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, simpson
from scipy.interpolate import CubicSpline, PchipInterpolator

from Vonkarmen_bone.scripts.analyze_compressible_lopez import (
    MR_VALUES,
    OUT,
    TW_VALUES,
    Z_COMPARE_MAX,
    compressible_profile,
    load_lopez_profiles,
    solve_energy_shapes,
    solve_universal_velocity,
)


GRID_SIZES = (251, 501, 1001, 2001, 4001, 8001, 16001, 32001)
MR = 0.3
QUAD_EPSABS = 1.0e-13
QUAD_EPSREL = 2.0e-12
FIELDS = ("F", "G", "H", "Fz", "Gz", "rho", "T")


def make_interpolants(comp, lopez, kind):
    constructor = CubicSpline if kind == "cubic" else PchipInterpolator
    comp_values = {
        "F": comp["F"], "G": comp["G"], "H": comp["H"],
        "Fz": comp["Fz"], "Gz": comp["Gz"],
        "rho": comp["rho"], "T": comp["T"],
    }
    lopez_values = {
        "F": lopez["F"], "G": lopez["G"], "H": lopez["H"],
        "Fz": lopez["Fp"], "Gz": lopez["Gp"],
        "rho": 2.0 - lopez["T"], "T": lopez["T"],
    }
    c = {
        name: constructor(comp["z"], values, extrapolate=False)
        for name, values in comp_values.items()
    }
    l = {
        name: constructor(lopez["z"], values, extrapolate=False)
        for name, values in lopez_values.items()
    }
    return c, l


def relative_l2_discrete(cfun, lfun, name, z, method):
    c = cfun[name](z)
    difference = lfun[name](z) - c
    reference = c - 1.0 if name == "T" else c
    integrate = np.trapz if method == "trapezoid" else simpson
    numerator = integrate(difference**2, z)
    denominator = integrate(reference**2, z)
    return np.sqrt(numerator / denominator)


def relative_l2_adaptive(cfun, lfun, name):
    offset = 1.0 if name == "T" else 0.0
    points = np.linspace(0.0, Z_COMPARE_MAX, 81)[1:-1]
    numerator = quad(
        lambda z: float((lfun[name](z) - cfun[name](z)) ** 2),
        0.0, Z_COMPARE_MAX, points=points,
        epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL, limit=800,
    )[0]
    denominator = quad(
        lambda z: float((cfun[name](z) - offset) ** 2),
        0.0, Z_COMPARE_MAX, points=points,
        epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL, limit=800,
    )[0]
    return np.sqrt(numerator / denominator)


def aggregate(field_errors):
    return {
        "velocity_error": max(field_errors[name] for name in ("F", "G", "H")),
        "derivative_error": max(field_errors[name] for name in ("Fz", "Gz")),
        "density_error": field_errors["rho"],
        "temperature_excess_error": field_errors["T"],
    }


def calculate_rows(compressible, lopez):
    rows = []
    continuous = {}
    for tw in TW_VALUES:
        interpolants = {}
        for kind in ("cubic", "pchip"):
            cfun, lfun = make_interpolants(compressible[tw], lopez[tw], kind)
            interpolants[kind] = (cfun, lfun)
            fields = {
                name: relative_l2_adaptive(cfun, lfun, name)
                for name in FIELDS
            }
            values = aggregate(fields)
            continuous[(tw, kind)] = values
            rows.append({
                "Tw": tw, "method": f"adaptive_{kind}", "N": 0,
                **fields, **values,
            })

        cfun, lfun = interpolants["cubic"]
        for n in GRID_SIZES:
            z = np.linspace(0.0, Z_COMPARE_MAX, n)
            for method in ("trapezoid", "simpson"):
                fields = {
                    name: relative_l2_discrete(cfun, lfun, name, z, method)
                    for name in FIELDS
                }
                rows.append({
                    "Tw": tw, "method": method, "N": n,
                    **fields, **aggregate(fields),
                })
    return rows, continuous


def write_csv(path, rows):
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def max_deviation(rows, continuous, method, n, key):
    selected = [row for row in rows if row["method"] == method and row["N"] == n]
    return max(abs(row[key] - continuous[(row["Tw"], "cubic")][key]) for row in selected)


def make_figure(rows, continuous, path):
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    for ax, key, title in zip(
        axes,
        ("velocity_error", "derivative_error"),
        ("velocity-error integration", "derivative-error integration"),
    ):
        for method, marker in (("trapezoid", "o"), ("simpson", "s")):
            errors = [max_deviation(rows, continuous, method, n, key) for n in GRID_SIZES]
            ax.loglog(GRID_SIZES, errors, marker + "-", ms=4, label=method)
        ax.set_xlabel("number of integration points")
        ax.set_ylabel("max absolute deviation from adaptive integral")
        ax.set_title(title)
        ax.grid(which="both", alpha=0.25)
        ax.legend()
    fig.suptitle(r"Error-integration convergence at $M_r=0.3$")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_report(path, rows, continuous):
    with path.open("w", encoding="ascii") as stream:
        stream.write("# Error-integration convergence at Mr=0.3\n\n")
        stream.write("The reference is adaptive quadrature of cubic-spline profile functions on z in [0,20].\n\n")
        stream.write("| Tw | E_velocity | E_derivative | trap N=4001 delta(Ev) | Simpson N=4001 delta(Ev) | cubic-PCHIP delta(Ev) |\n")
        stream.write("|---:|---:|---:|---:|---:|---:|\n")
        for tw in TW_VALUES:
            ref = continuous[(tw, "cubic")]
            trap = next(row for row in rows if row["Tw"] == tw and row["method"] == "trapezoid" and row["N"] == 4001)
            simp = next(row for row in rows if row["Tw"] == tw and row["method"] == "simpson" and row["N"] == 4001)
            pchip = continuous[(tw, "pchip")]
            stream.write(
                f"| {tw:.3f} | {ref['velocity_error']:.10e} | {ref['derivative_error']:.10e} | "
                f"{abs(trap['velocity_error']-ref['velocity_error']):.3e} | "
                f"{abs(simp['velocity_error']-ref['velocity_error']):.3e} | "
                f"{abs(pchip['velocity_error']-ref['velocity_error']):.3e} |\n"
            )
        stream.write("\n## Maximum integration deviations over all wall temperatures\n\n")
        stream.write("| method | N | velocity | derivative |\n|:--|---:|---:|---:|\n")
        for method in ("trapezoid", "simpson"):
            for n in GRID_SIZES:
                ev = max_deviation(rows, continuous, method, n, "velocity_error")
                ed = max_deviation(rows, continuous, method, n, "derivative_error")
                stream.write(f"| {method} | {n} | {ev:.6e} | {ed:.6e} |\n")


def main():
    if MR_VALUES != (MR,):
        raise RuntimeError("The comparison configuration must contain only Mr=0.3")
    eta, velocity = solve_universal_velocity()
    f, fp, q, qp, residual = solve_energy_shapes(eta, velocity)
    energy = (f, fp, q, qp)
    lopez = load_lopez_profiles()
    compressible = {
        tw: compressible_profile(eta, velocity, energy, tw, MR)
        for tw in TW_VALUES
    }
    rows, continuous = calculate_rows(compressible, lopez)
    write_csv(OUT / "integration_convergence.csv", rows)
    make_figure(rows, continuous, OUT / "integration_convergence.png")
    write_report(OUT / "integration_convergence.md", rows, continuous)
    print(f"energy boundary residual = {residual:.3e}")
    print(f"results: {(OUT / 'integration_convergence.md').resolve()}")


if __name__ == "__main__":
    main()
