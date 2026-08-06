#!/usr/bin/env python3
"""Verify Tecplot DAT headers, zones, row counts, and numeric columns."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ROOT = SCRIPT_DIR.parent / "tecplot"


def verify_dat(path: Path) -> tuple[int, int, int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    variables_line = next(
        (line for line in lines if line.strip().upper().startswith("VARIABLES")),
        None,
    )
    if variables_line is None:
        raise ValueError(f"Missing VARIABLES line: {path}")
    variable_count = len(re.findall(r'"(?:[^"\\]|\\.)*"', variables_line))
    if variable_count < 1:
        raise ValueError(f"No variables found: {path}")

    index = 0
    zones = 0
    total_rows = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line.upper().startswith("ZONE "):
            index += 1
            continue
        match = re.search(r"\bI\s*=\s*(\d+)", line, re.IGNORECASE)
        if match is None:
            raise ValueError(f"ZONE without I= count in {path}: {line}")
        expected = int(match.group(1))
        zones += 1
        index += 1
        while index < len(lines) and lines[index].strip().upper().startswith(
            "AUXDATA"
        ):
            index += 1
        for _ in range(expected):
            if index >= len(lines):
                raise ValueError(f"Unexpected EOF in {path}")
            data_line = lines[index].strip()
            if not data_line or data_line.upper().startswith("ZONE "):
                raise ValueError(f"Missing data row in {path} near line {index + 1}")
            tokens = data_line.split()
            if len(tokens) != variable_count:
                raise ValueError(
                    f"Column mismatch in {path} line {index + 1}: "
                    f"expected {variable_count}, got {len(tokens)}"
                )
            for token in tokens:
                float(token)
            index += 1
            total_rows += 1
    if zones < 1:
        raise ValueError(f"No zones found: {path}")
    return zones, total_rows, variable_count


def read_dat_values(path: Path) -> tuple[list[str], np.ndarray]:
    """Return variable names and all point-packed rows from a DAT file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    variables_line = next(
        line for line in lines if line.strip().upper().startswith("VARIABLES")
    )
    variables = re.findall(r'"((?:[^"\\]|\\.)*)"', variables_line)
    rows: list[list[float]] = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line.upper().startswith("ZONE "):
            index += 1
            continue
        expected = int(re.search(r"\bI\s*=\s*(\d+)", line, re.IGNORECASE).group(1))
        index += 1
        while index < len(lines) and lines[index].strip().upper().startswith(
            "AUXDATA"
        ):
            index += 1
        for _ in range(expected):
            rows.append([float(value) for value in lines[index].split()])
            index += 1
    return variables, np.asarray(rows, dtype=float)


def verify_source_values(root: Path) -> None:
    """Check representative CSV and NPZ values against generated DAT files."""
    rotor_root = root.parent
    csv_relative = Path(
        "boussinesq_singularity_convergence/three_solutions_Tw1.160.csv"
    )
    csv_path = rotor_root / "data" / csv_relative
    dat_path = root / "data" / csv_relative.with_suffix(".dat")
    variables, dat_values = read_dat_values(dat_path)
    with csv_path.open(encoding="utf-8-sig", newline="") as handle:
        records = [
            {key.strip(): value for key, value in record.items()}
            for record in csv.DictReader(handle)
        ]
    csv_values = np.asarray(
        [[float(record[name]) for name in variables] for record in records],
        dtype=float,
    )
    if csv_values.shape != dat_values.shape or not np.allclose(
        csv_values, dat_values, rtol=2e-14, atol=2e-14, equal_nan=True
    ):
        raise SystemExit(f"CSV value mismatch: {csv_path} -> {dat_path}")

    npz_path = rotor_root / "reference" / "baseflow_Res1000.npz"
    npz_dat_path = root / "reference" / "baseflow_Res1000.dat"
    variables, dat_values = read_dat_values(npz_dat_path)
    with np.load(npz_path) as source:
        npz_values = np.column_stack([source[name] for name in variables])
    if npz_values.shape != dat_values.shape or not np.allclose(
        npz_values, dat_values, rtol=2e-14, atol=2e-14, equal_nan=True
    ):
        raise SystemExit(f"NPZ value mismatch: {npz_path} -> {npz_dat_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tecplot-root", type=Path, default=DEFAULT_ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.tecplot_root.resolve()
    manifest = json.loads(
        (root / "conversion_manifest.json").read_text(encoding="utf-8")
    )
    checked_rows = 0
    for item in manifest:
        path = root / item["output"]
        zones, rows, variable_count = verify_dat(path)
        if zones != int(item["zones"]):
            raise SystemExit(
                f"Zone mismatch for {path}: manifest={item['zones']}, file={zones}"
            )
        if rows != int(item["rows"]):
            raise SystemExit(
                f"Row mismatch for {path}: manifest={item['rows']}, file={rows}"
            )
        if variable_count != len(item["numeric_variables"]):
            raise SystemExit(
                f"Variable mismatch for {path}: manifest="
                f"{len(item['numeric_variables'])}, file={variable_count}"
            )
        checked_rows += rows
    verify_source_values(root)
    print(
        f"Tecplot verification passed: {len(manifest)} files, "
        f"{checked_rows} rows; representative CSV/NPZ values match"
    )


if __name__ == "__main__":
    main()
