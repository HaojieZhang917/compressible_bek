#!/usr/bin/env python3
"""Convert Lopez/traditional comparison CSV outputs to Tecplot ASCII DAT."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
ROOT = WORKSPACE_ROOT / "lopez_boussinesq_comparison"


def read_csv(path: Path):
    with path.open("r", newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def format_value(value) -> str:
    value = float(value)
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "+Inf" if value > 0 else "-Inf"
    return f"{value:.12e}"


def write_dat(path: Path, title: str, variables, zones):
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(f'TITLE = "{title}"\n')
        stream.write("VARIABLES = " + " ".join(f'\"{name}\"' for name in variables) + "\n")
        for zone_name, rows in zones:
            stream.write(
                f'ZONE T="{zone_name}", I={len(rows)}, DATAPACKING=POINT\n'
            )
            for row in rows:
                stream.write(" ".join(format_value(row[name]) for name in variables) + "\n")


def convert_profiles(csv_name: str, dat_name: str):
    rows = read_csv(ROOT / csv_name)
    variables = ("z", "H", "F", "G", "T", "Fp", "Gp", "Tp")
    grouped = defaultdict(list)
    for row in rows:
        key = (row["model"], number(row["Tw"]))
        grouped[key].append({name: number(row[name]) for name in variables})
    zones = []
    for (model, tw), values in sorted(grouped.items(), key=lambda item: item[0][1]):
        values.sort(key=lambda row: row["z"])
        zones.append((f"{model}_Tw_{tw:.4f}", values))
    write_dat(ROOT / dat_name, f"{dat_name[:-4]} base-flow profiles", variables, zones)


def convert_metrics():
    rows = read_csv(ROOT / "baseflow_metrics.csv")
    variables = (
        "Tw", "Hinf", "Fmax", "z_Fmax", "Fp0", "Gp0", "Tp0",
        "thermal_thickness", "rho_wall_linear", "nodes",
    )
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["model"]].append({name: number(row[name]) for name in variables})
    zones = []
    for model, values in sorted(grouped.items()):
        values.sort(key=lambda row: row["Tw"])
        zones.append((model, values))
    write_dat(
        ROOT / "baseflow_metrics.dat",
        "Traditional and Lopez Boussinesq base-flow metrics",
        variables,
        zones,
    )


def convert_differences():
    rows = read_csv(ROOT / "profile_differences.csv")
    variables = tuple(rows[0].keys())
    values = [{name: number(row[name]) for name in variables} for row in rows]
    write_dat(
        ROOT / "profile_differences.dat",
        "Lopez minus traditional profile differences",
        variables,
        [("profile_differences", values)],
    )


def convert_domain_convergence():
    rows = read_csv(ROOT / "lopez_domain_convergence.csv")
    variables = (
        "Tw", "Hinf", "Fmax", "z_Fmax", "Fp0", "Gp0", "Tp0",
        "thermal_thickness", "rho_wall_linear", "nodes",
    )
    grouped = defaultdict(list)
    for row in rows:
        grouped[number(row["zmax"])].append(
            {name: number(row[name]) for name in variables}
        )
    zones = []
    for zmax, values in sorted(grouped.items()):
        values.sort(key=lambda row: row["Tw"])
        zones.append((f"zmax_{zmax:g}", values))
    write_dat(
        ROOT / "lopez_domain_convergence.dat",
        "Lopez generalized Boussinesq domain convergence",
        variables,
        zones,
    )


def combine_profile_dat():
    rows = read_csv(ROOT / "profiles_traditional.csv")
    rows.extend(read_csv(ROOT / "profiles_lopez.csv"))
    variables = ("z", "H", "F", "G", "T", "Fp", "Gp", "Tp")
    grouped = defaultdict(list)
    for row in rows:
        key = (row["model"], number(row["Tw"]))
        grouped[key].append({name: number(row[name]) for name in variables})
    zones = []
    for (model, tw), values in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
        values.sort(key=lambda row: row["z"])
        zones.append((f"{model}_Tw_{tw:.4f}", values))
    write_dat(
        ROOT / "profiles_all_models.dat",
        "Traditional and Lopez Boussinesq base-flow profiles",
        variables,
        zones,
    )


def main():
    convert_profiles("profiles_traditional.csv", "profiles_traditional.dat")
    convert_profiles("profiles_lopez.csv", "profiles_lopez.dat")
    combine_profile_dat()
    convert_metrics()
    convert_differences()
    convert_domain_convergence()
    for path in sorted(ROOT.glob("*.dat")):
        print(f"{path}  ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
