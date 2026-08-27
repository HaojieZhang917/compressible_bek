from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np


def numeric_rows(path: Path, columns: int) -> np.ndarray:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if len(values) != columns:
            continue
        try:
            rows.append([float(value) for value in values])
        except ValueError:
            continue
    return np.asarray(rows, dtype=float)


def temporal_keys(directory: Path) -> set[tuple[float, float]]:
    keys: set[tuple[float, float]] = set()
    for path in directory.glob("R=*.dat"):
        for line in path.read_text(encoding="utf-8").splitlines()[3:]:
            values = line.split()
            if len(values) < 5:
                continue
            try:
                radius, _, beta, _, growth = map(float, values[:5])
            except ValueError:
                continue
            if growth > 0:
                keys.add((round(radius, 10), round(beta, 10)))
    return keys


def percentile_summary(values: np.ndarray) -> str:
    if values.size == 0:
        return "none"
    return ", ".join(
        f"p{percentile}={np.percentile(values, percentile):.6e}"
        for percentile in (50, 90, 99, 100)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_directory", type=Path)
    parser.add_argument("temporal_directory", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--suspects", type=Path)
    parser.add_argument("--suspect-count", type=int, default=40)
    args = parser.parse_args()

    envelope = numeric_rows(args.result_directory / "spatial_growth_envelope.dat", 11)
    samples = numeric_rows(args.result_directory / "spatial_growth_samples.dat", 12)
    if envelope.size == 0 or samples.size == 0:
        raise RuntimeError("Merged result files contain no numeric rows")

    result_keys = {
        (round(row[0], 10), round(row[1], 10)) for row in envelope
    }
    expected_keys = temporal_keys(args.temporal_directory)
    duplicate_count = len(envelope) - len(result_keys)
    missing = expected_keys - result_keys
    unexpected = result_keys - expected_keys
    finite = bool(np.isfinite(envelope).all() and np.isfinite(samples).all())

    growth_identity = np.max(np.abs(envelope[:, 2] + envelope[:, 7]))
    omega_bar_identity = np.max(np.abs(envelope[:, 5] - envelope[:, 0] * envelope[:, 3]))
    residuals = envelope[:, 10]

    by_key = {
        (round(row[0], 10), round(row[1], 10)): row for row in envelope
    }
    radii = sorted({key[0] for key in by_key})
    betas = sorted({key[1] for key in by_key})
    radius_step = float(np.median(np.diff(radii)))
    beta_step = float(np.median(np.diff(betas)))
    neighbor_growth_jumps: list[float] = []
    neighbor_frequency_jumps: list[float] = []
    curvature_rows: list[tuple[float, ...]] = []
    for key, row in by_key.items():
        for neighbor_key in (
            (round(key[0] + radius_step, 10), key[1]),
            (key[0], round(key[1] + beta_step, 10)),
        ):
            neighbor = by_key.get(neighbor_key)
            if neighbor is None:
                continue
            neighbor_growth_jumps.append(abs(row[2] - neighbor[2]))
            neighbor_frequency_jumps.append(abs(row[3] - neighbor[3]))

        growth_curvatures: list[float] = []
        frequency_curvatures: list[float] = []
        for lower_key, upper_key in (
            (
                (round(key[0] - radius_step, 10), key[1]),
                (round(key[0] + radius_step, 10), key[1]),
            ),
            (
                (key[0], round(key[1] - beta_step, 10)),
                (key[0], round(key[1] + beta_step, 10)),
            ),
        ):
            lower = by_key.get(lower_key)
            upper = by_key.get(upper_key)
            if lower is None or upper is None:
                continue
            growth_curvatures.append(abs(row[2] - 0.5 * (lower[2] + upper[2])))
            frequency_curvatures.append(abs(row[3] - 0.5 * (lower[3] + upper[3])))
        if growth_curvatures:
            curvature_rows.append(
                (
                    key[0],
                    key[1],
                    row[2],
                    row[3],
                    max(growth_curvatures),
                    max(frequency_curvatures),
                    row[4],
                    row[10],
                )
            )

    max_index = int(np.argmax(envelope[:, 2]))
    max_row = envelope[max_index]
    sample_status_failures = int(np.count_nonzero(samples[:, 11] != 1))
    passed = (
        duplicate_count == 0
        and not missing
        and not unexpected
        and finite
        and sample_status_failures == 0
        and float(np.max(residuals)) < 1e-6
        and growth_identity < 1e-9
        and omega_bar_identity < 1e-8
    )

    lines = [
        "Zarf spatial-envelope validation",
        f"status: {'PASS' if passed else 'FAIL'}",
        f"sample rows: {len(samples)}",
        f"envelope rows: {len(envelope)}",
        f"expected temporal-Zarf points: {len(expected_keys)}",
        f"duplicate envelope keys: {duplicate_count}",
        f"missing keys: {len(missing)}",
        f"unexpected keys: {len(unexpected)}",
        f"non-finite values present: {not finite}",
        f"failed sample statuses: {sample_status_failures}",
        f"residuals: {percentile_summary(residuals)}",
        f"max |growth + alpha_i|: {growth_identity:.6e}",
        f"max |omega_bar - R*omega|: {omega_bar_identity:.6e}",
        "neighbor growth jumps: "
        + percentile_summary(np.asarray(neighbor_growth_jumps)),
        "neighbor omega jumps: "
        + percentile_summary(np.asarray(neighbor_frequency_jumps)),
        (
            "maximum growth: "
            f"R={max_row[0]:.6f}, beta={max_row[1]:.6f}, "
            f"growth={max_row[2]:.12e}, omega={max_row[3]:.12e}, "
            f"alpha={max_row[6]:.12e}{max_row[7]:+.12e}i, "
            f"residual={max_row[10]:.6e}"
        ),
    ]
    report = "\n".join(lines) + "\n"
    print(report, end="")
    report_path = args.report or args.result_directory / "validation_report.txt"
    report_path.write_text(report, encoding="utf-8")
    suspect_path = args.suspects or args.result_directory / "suspicious_points.dat"
    curvature_rows.sort(key=lambda row: (row[4], row[5]), reverse=True)
    with suspect_path.open("w", encoding="utf-8") as stream:
        stream.write(
            'VARIABLES="R","beta","growth_rate","omega_opt",'
            '"growth_curvature","omega_curvature","branch_id","residual"\n'
        )
        stream.write(
            f'ZONE T="largest local curvature", I={min(args.suspect_count, len(curvature_rows))}, F=POINT\n'
        )
        for row in curvature_rows[: args.suspect_count]:
            stream.write(" ".join(f"{value:.12e}" for value in row) + "\n")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
