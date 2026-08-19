from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np
from repair_zarf_spatial_growth_fixed_beta import (
    detect_flags,
    format_row,
    numeric_rows,
    point_key,
)


def select_forward_targets(
    envelope: np.ndarray,
    flagged: set[tuple[float, float]],
    r_min: float,
    r_max: float,
    beta_min: float,
    beta_max: float,
) -> tuple[set[tuple[float, float]], dict[float, float]]:
    starts: dict[float, float] = {}
    for radius, beta in flagged:
        if not (r_min <= radius <= r_max and beta_min <= beta <= beta_max):
            continue
        beta_key = round(beta, 10)
        starts[beta_key] = min(starts.get(beta_key, radius), radius)

    targets: set[tuple[float, float]] = set()
    for row in envelope:
        radius = float(row[0])
        beta_key = round(float(row[1]), 10)
        start_radius = starts.get(beta_key)
        if start_radius is None:
            continue
        if r_min <= radius <= r_max and radius >= start_radius:
            targets.add(point_key(radius, beta_key))
    return targets, starts


def extrapolate_value(
    history: list[np.ndarray],
    target_radius: float,
    column: int,
    stencil_size: int,
) -> float:
    source_rows = history[-stencil_size:]
    source_radii = np.asarray([row[0] for row in source_rows], dtype=float)
    source_values = np.asarray([row[column] for row in source_rows], dtype=float)
    if len(source_rows) < 2:
        raise RuntimeError("Forward continuation requires at least two source points")
    shifted_radii = source_radii - target_radius
    slope, intercept = np.polyfit(shifted_radii, source_values, 1)
    return float(intercept)


def forward_continue(
    envelope: np.ndarray,
    targets: set[tuple[float, float]],
    stencil_size: int,
) -> tuple[np.ndarray, list[dict[str, float | str]]]:
    rows_by_beta: dict[float, list[tuple[int, np.ndarray]]] = {}
    for row_index, row in enumerate(envelope):
        beta_key = round(float(row[1]), 10)
        rows_by_beta.setdefault(beta_key, []).append((row_index, row))

    repaired = envelope.copy()
    repair_log: list[dict[str, float | str]] = []
    for beta_key, beta_rows in rows_by_beta.items():
        history: list[np.ndarray] = []
        for row_index, original_row in sorted(beta_rows, key=lambda item: item[1][0]):
            radius = float(original_row[0])
            target_key = point_key(radius, beta_key)
            if target_key not in targets:
                current_row = original_row.copy()
                history.append(current_row)
                continue
            if len(history) < 2:
                repair_log.append(
                    {
                        "R": radius,
                        "beta": beta_key,
                        "method": "unresolved-insufficient-forward-history",
                        "source_R_min": float("nan"),
                        "source_R_max": float("nan"),
                        "original_growth": original_row[2],
                        "repaired_growth": float("nan"),
                        "original_omega": original_row[3],
                        "repaired_omega": float("nan"),
                    }
                )
                history.append(original_row.copy())
                continue

            current_row = original_row.copy()
            growth = extrapolate_value(history, radius, 2, stencil_size)
            frequency = extrapolate_value(history, radius, 3, stencil_size)
            alpha_real = extrapolate_value(history, radius, 6, stencil_size)
            current_row[2] = growth
            current_row[3] = frequency
            current_row[4] = -1.0
            current_row[5] = radius * frequency
            current_row[6] = alpha_real
            current_row[7] = -growth
            current_row[10] = -1.0
            repaired[row_index] = current_row
            repair_log.append(
                {
                    "R": radius,
                    "beta": beta_key,
                    "method": "forward-fixed-beta-linear",
                    "source_R_min": history[-stencil_size][0],
                    "source_R_max": history[-1][0],
                    "original_growth": original_row[2],
                    "repaired_growth": growth,
                    "original_omega": original_row[3],
                    "repaired_omega": frequency,
                }
            )
            history.append(current_row)
    return repaired, repair_log


def write_envelope(
    source_path: Path,
    output_path: Path,
    repaired: np.ndarray,
    target_count: int,
    unresolved_count: int,
) -> None:
    source_lines = source_path.read_text(encoding="utf-8").splitlines()
    zone_index = next(
        index for index, line in enumerate(source_lines) if line.startswith("ZONE")
    )
    metadata = [
        'DATASETAUXDATA repair_method="forward fixed-beta linear continuation"',
        f'DATASETAUXDATA repair_targets="{target_count}"',
        f'DATASETAUXDATA unresolved_points="{unresolved_count}"',
        'DATASETAUXDATA original_values_preserved="yes; repaired nodes have branch_id=-1 and residual=-1"',
    ]
    output_lines = source_lines[:zone_index] + metadata + source_lines[zone_index : zone_index + 1]
    output_lines.extend(format_row(row) for row in repaired)
    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def patch_fetriangle(
    source_path: Path,
    output_path: Path,
    repaired_by_key: dict[tuple[float, float], np.ndarray],
    target_count: int,
    unresolved_count: int,
) -> int:
    source_lines = source_path.read_text(encoding="utf-8").splitlines()
    zone_index = next(
        index
        for index, line in enumerate(source_lines)
        if "ZONETYPE=FETRIANGLE" in line
    )
    node_count_match = re.search(r"\bN=(\d+)", source_lines[zone_index])
    if node_count_match is None:
        raise RuntimeError("Could not read FE-triangle node count")
    node_count = int(node_count_match.group(1))
    node_start = zone_index + 1
    node_end = node_start + node_count
    patched_lines = list(source_lines)
    replaced_nodes = 0
    for line_index in range(node_start, node_end):
        values = patched_lines[line_index].split()
        if len(values) != 11:
            continue
        try:
            node_values = np.asarray([float(value) for value in values], dtype=float)
        except ValueError:
            continue
        repaired_row = repaired_by_key.get(point_key(node_values[0], node_values[1]))
        if repaired_row is None:
            continue
        patched_lines[line_index] = format_row(repaired_row)
        replaced_nodes += 1

    metadata = [
        'DATASETAUXDATA repair_method="forward fixed-beta linear continuation on envelope nodes"',
        f'DATASETAUXDATA repair_targets="{target_count}"',
        f'DATASETAUXDATA unresolved_points="{unresolved_count}"',
        f'DATASETAUXDATA replaced_FE_nodes="{replaced_nodes}"',
    ]
    patched_lines = patched_lines[:zone_index] + metadata + patched_lines[zone_index:]
    output_path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")
    return replaced_nodes


def write_log(path: Path, entries: list[dict[str, float | str]]) -> None:
    columns = (
        "R",
        "beta",
        "method",
        "source_R_min",
        "source_R_max",
        "original_growth",
        "repaired_growth",
        "original_omega",
        "repaired_omega",
    )
    lines = [" ".join(columns)]
    for entry in entries:
        lines.append(
            f"{float(entry['R']):.12e} {float(entry['beta']):.12e} {entry['method']} "
            f"{float(entry['source_R_min']):.12e} {float(entry['source_R_max']):.12e} "
            f"{float(entry['original_growth']):.12e} {float(entry['repaired_growth']):.12e} "
            f"{float(entry['original_omega']):.12e} {float(entry['repaired_omega']):.12e}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_directory", type=Path)
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--r-min", type=float, default=350.0)
    parser.add_argument("--r-max", type=float, default=500.0)
    parser.add_argument("--beta-min", type=float, default=0.10)
    parser.add_argument("--beta-max", type=float, default=0.20)
    parser.add_argument("--stencil-size", type=int, default=3)
    args = parser.parse_args()
    if args.stencil_size < 2:
        raise ValueError("stencil-size must be at least 2")

    result_directory = args.result_directory.resolve()
    output_directory = (
        args.output_directory.resolve()
        if args.output_directory is not None
        else result_directory / "forward_continuation_first_mode_R"
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    envelope_path = result_directory / "spatial_growth_envelope.dat"
    fetriangle_path = result_directory / "spatial_growth_cloud_fetriangle.dat"
    samples_path = result_directory / "spatial_growth_samples.dat"
    envelope = numeric_rows(envelope_path, 11)
    if envelope.size == 0:
        raise RuntimeError("No envelope rows found")

    flagged, thresholds = detect_flags(envelope)
    targets, starts = select_forward_targets(
        envelope,
        flagged,
        args.r_min,
        args.r_max,
        args.beta_min,
        args.beta_max,
    )
    repaired, repair_log = forward_continue(envelope, targets, args.stencil_size)
    repaired_count = sum(
        1
        for entry in repair_log
        if entry["method"] == "forward-fixed-beta-linear"
    )
    unresolved_count = len(repair_log) - repaired_count
    repaired_by_key = {
        point_key(row[0], row[1]): row for row in repaired
    }

    output_envelope = output_directory / "spatial_growth_envelope.dat"
    output_fetriangle = output_directory / "spatial_growth_cloud_fetriangle.dat"
    output_log = output_directory / "forward_interpolated_points.dat"
    output_targets = output_directory / "forward_targets.dat"
    write_envelope(
        envelope_path,
        output_envelope,
        repaired,
        len(targets),
        unresolved_count,
    )
    replaced_nodes = patch_fetriangle(
        fetriangle_path,
        output_fetriangle,
        {key: repaired_by_key[key] for key in targets if key in repaired_by_key},
        len(targets),
        unresolved_count,
    )
    write_log(output_log, repair_log)
    output_targets.write_text(
        "\n".join(f"{radius:.12e} {beta:.12e}" for radius, beta in sorted(targets))
        + "\n",
        encoding="utf-8",
    )
    if samples_path.exists():
        shutil.copy2(samples_path, output_directory / "spatial_growth_samples.dat")

    start_lines = [
        f"beta={beta:.12e} start_R={radius:.12e}"
        for beta, radius in sorted(starts.items())
    ]
    report_lines = [
        "Forward continuation spatial-growth repair",
        f"source: {result_directory}",
        f"window: R=[{args.r_min}, {args.r_max}], beta=[{args.beta_min}, {args.beta_max}]",
        f"continuation beta lines: {len(starts)}",
        f"forward target points: {len(targets)}",
        f"repaired points: {repaired_count}",
        f"unresolved points: {unresolved_count}",
        f"replaced FE nodes: {replaced_nodes}",
        f"stencil size: {args.stencil_size}",
        "thresholds: " + " ".join(f"{value:.12e}" for value in thresholds),
        "continuation starts:",
        *start_lines,
        f"envelope: {output_envelope}",
        f"FE triangle: {output_fetriangle}",
        f"repair log: {output_log}",
    ]
    (output_directory / "forward_repair_report.txt").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )
    print("\n".join(report_lines))


if __name__ == "__main__":
    main()
