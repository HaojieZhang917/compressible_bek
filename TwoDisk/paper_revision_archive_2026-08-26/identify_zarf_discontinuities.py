from __future__ import annotations

import argparse
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("envelope", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--percentile", type=float, default=99.0)
    parser.add_argument("--rings", type=int, default=1)
    args = parser.parse_args()

    envelope = numeric_rows(args.envelope, 11)
    by_key = {
        (round(row[0], 10), round(row[1], 10)): row for row in envelope
    }
    radii = sorted({key[0] for key in by_key})
    betas = sorted({key[1] for key in by_key})
    radius_step = float(np.median(np.diff(radii)))
    beta_step = float(np.median(np.diff(betas)))

    metrics: dict[tuple[float, float], np.ndarray] = {}
    for key, row in by_key.items():
        values = np.zeros(5)
        for neighbor_key in (
            (round(key[0] - radius_step, 10), key[1]),
            (round(key[0] + radius_step, 10), key[1]),
            (key[0], round(key[1] - beta_step, 10)),
            (key[0], round(key[1] + beta_step, 10)),
        ):
            neighbor = by_key.get(neighbor_key)
            if neighbor is None:
                continue
            values[0] = max(values[0], abs(row[2] - neighbor[2]))
            values[1] = max(values[1], abs(row[3] - neighbor[3]))
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
            midpoint = 0.5 * (lower + upper)
            values[2] = max(values[2], abs(row[2] - midpoint[2]))
            values[3] = max(values[3], abs(row[3] - midpoint[3]))
            alpha = complex(row[6], row[7])
            alpha_midpoint = 0.5 * (
                complex(lower[6], lower[7]) + complex(upper[6], upper[7])
            )
            values[4] = max(values[4], abs(alpha - alpha_midpoint))
        metrics[key] = values

    metric_array = np.asarray(list(metrics.values()))
    percentile_thresholds = np.percentile(metric_array, args.percentile, axis=0)
    absolute_floors = np.asarray((0.01, 0.01, 0.005, 0.005, 0.02))
    thresholds = np.maximum(percentile_thresholds, absolute_floors)
    flagged = {
        key for key, values in metrics.items() if np.any(values >= thresholds)
    }
    expanded = set(flagged)
    frontier = set(flagged)
    for _ in range(args.rings):
        next_frontier: set[tuple[float, float]] = set()
        for key in frontier:
            for neighbor_key in (
                (round(key[0] - radius_step, 10), key[1]),
                (round(key[0] + radius_step, 10), key[1]),
                (key[0], round(key[1] - beta_step, 10)),
                (key[0], round(key[1] + beta_step, 10)),
            ):
                if neighbor_key in by_key and neighbor_key not in expanded:
                    next_frontier.add(neighbor_key)
        expanded.update(next_frontier)
        frontier = next_frontier

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        for radius, beta in sorted(expanded):
            stream.write(f"{radius:.12e} {beta:.12e}\n")

    labels = (
        "growth edge jump",
        "omega edge jump",
        "growth curvature",
        "omega curvature",
        "alpha curvature",
    )
    lines = [
        f"envelope points: {len(envelope)}",
        f"raw flagged points: {len(flagged)}",
        f"expanded target points: {len(expanded)}",
        f"neighbor rings: {args.rings}",
    ]
    for label, threshold, maximum in zip(labels, thresholds, np.max(metric_array, axis=0)):
        lines.append(f"{label}: threshold={threshold:.12e}, maximum={maximum:.12e}")
    report = "\n".join(lines) + "\n"
    print(report, end="")
    report_path = args.report or args.output.with_suffix(".report.txt")
    report_path.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
