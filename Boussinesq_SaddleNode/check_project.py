#!/usr/bin/env python3
"""Read-only structural check for the consolidated saddle-node project."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent

EXPECTED = {
    "von Karman finite solver": ROOT
    / "von_karman/scripts/investigate_boussinesq_fold.jl",
    "von Karman mapped solver": ROOT
    / "von_karman/scripts/analyze_vonkarman_infinite_mapping.jl",
    "mapped fold summary": ROOT
    / "von_karman/data/infinite_mapping/confirmed_fold_summary.csv",
    "mapped temporal summary": ROOT
    / "von_karman/data/infinite_mapping/upper_branch_stability_summary.csv",
    "rotor-stator solver": ROOT
    / "rotor_stator/scripts/two_disk_boussinesq_singularity.jl",
    "rotor-stator reference": ROOT
    / "rotor_stator/reference/baseflow_Res1000.npz",
    "rotor-stator folds": ROOT
    / "rotor_stator/data/boussinesq_singularity_results"
    / "folds_Re1000_traditional_centrifugal.json",
    "cross-model solver": ROOT
    / "cross_model/scripts/compare_saddle_node_dynamics.jl",
    "cross-model report": ROOT
    / "cross_model/data/dynamical_singularity_comparison"
    / "BEGINNER_DYNAMICAL_SYSTEMS_COMPARISON_REPORT.md",
}


def csv_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as stream:
        return sum(1 for _row in csv.DictReader(stream))


def main() -> None:
    missing = [label for label, path in EXPECTED.items() if not path.is_file()]
    if missing:
        raise SystemExit("Missing project assets: " + ", ".join(missing))

    fold_csv = EXPECTED["mapped fold summary"]
    if csv_rows(fold_csv) < 2:
        raise SystemExit("Mapped fold summary does not contain two folds")

    with EXPECTED["rotor-stator folds"].open(encoding="utf-8") as stream:
        rotor_folds = json.load(stream)
    if len(rotor_folds) != 2:
        raise SystemExit("Expected two Re_h=1000 rotor-stator folds")

    file_count = sum(1 for path in ROOT.rglob("*") if path.is_file())
    size_bytes = sum(path.stat().st_size for path in ROOT.rglob("*") if path.is_file())
    print(f"Project check passed: {file_count} files, {size_bytes / 2**20:.2f} MiB")
    print(f"Mapped fold rows: {csv_rows(fold_csv)}")
    print(f"Rotor-stator folds: {len(rotor_folds)}")


if __name__ == "__main__":
    main()
