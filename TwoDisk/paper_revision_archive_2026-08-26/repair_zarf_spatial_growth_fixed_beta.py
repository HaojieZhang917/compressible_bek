from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


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


def point_key(radius: float, beta: float) -> tuple[float, float]:
    return round(float(radius), 10), round(float(beta), 10)


def detect_flags(envelope: np.ndarray) -> tuple[set[tuple[float, float]], np.ndarray]:
    by_key = {
        point_key(row[0], row[1]): row for row in envelope
    }
    radii = sorted({key[0] for key in by_key})
    betas = sorted({key[1] for key in by_key})
    if len(radii) < 2 or len(betas) < 2:
        raise RuntimeError("The envelope does not contain a two-dimensional grid")
    radius_step = float(np.median(np.diff(radii)))
    beta_step = float(np.median(np.diff(betas)))

    metrics: dict[tuple[float, float], np.ndarray] = {}
    for key, row in by_key.items():
        values = np.zeros(5)
        neighbor_keys = (
            (round(key[0] - radius_step, 10), key[1]),
            (round(key[0] + radius_step, 10), key[1]),
            (key[0], round(key[1] - beta_step, 10)),
            (key[0], round(key[1] + beta_step, 10)),
        )
        for neighbor_key in neighbor_keys:
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
    percentile_thresholds = np.percentile(metric_array, 99.0, axis=0)
    absolute_floors = np.asarray((0.01, 0.01, 0.005, 0.005, 0.02))
    thresholds = np.maximum(percentile_thresholds, absolute_floors)
    flagged = {
        key
        for key, values in metrics.items()
        if np.any(values >= thresholds)
    }
    return flagged, thresholds


def local_pchip(
    valid_rows: list[np.ndarray],
    target_radius: float,
    column: int,
    stencil_size: int,
    allow_edge_extrapolation: bool,
) -> tuple[float, str, float, float] | None:
    sorted_rows = sorted(valid_rows, key=lambda row: row[0])
    lower_rows = [row for row in sorted_rows if row[0] < target_radius]
    upper_rows = [row for row in sorted_rows if row[0] > target_radius]
    has_two_sides = bool(lower_rows and upper_rows)
    if not has_two_sides and not allow_edge_extrapolation:
        return None

    if has_two_sides:
        selected_rows = lower_rows[-stencil_size:] + upper_rows[:stencil_size]
        method = "fixed-beta-local-pchip"
    elif lower_rows:
        selected_rows = lower_rows[-max(2, stencil_size):]
        method = "fixed-beta-one-sided-pchip"
    else:
        selected_rows = upper_rows[:max(2, stencil_size)]
        method = "fixed-beta-one-sided-pchip"

    selected_rows = sorted(selected_rows, key=lambda row: row[0])
    if len(selected_rows) < 2:
        return None
    source_radii = np.asarray([row[0] for row in selected_rows], dtype=float)
    source_values = np.asarray([row[column] for row in selected_rows], dtype=float)
    interpolator = PchipInterpolator(source_radii, source_values, extrapolate=True)
    interpolated_value = float(interpolator(target_radius))
    lower_radius = float(lower_rows[-1][0]) if lower_rows else float("nan")
    upper_radius = float(upper_rows[0][0]) if upper_rows else float("nan")
    return interpolated_value, method, lower_radius, upper_radius


def repair_envelope(
    envelope: np.ndarray,
    targets: set[tuple[float, float]],
    stencil_size: int,
    allow_edge_extrapolation: bool,
) -> tuple[np.ndarray, list[dict[str, float | str]]]:
    by_key = {
        point_key(row[0], row[1]): row for row in envelope
    }
    repaired = envelope.copy()
    repair_log: list[dict[str, float | str]] = []
    for target_key in sorted(targets):
        original_row = by_key.get(target_key)
        if original_row is None:
            continue
        beta_rows = [
            row
            for key, row in by_key.items()
            if key[1] == target_key[1] and key not in targets
        ]
        growth_value = local_pchip(
            beta_rows,
            target_key[0],
            2,
            stencil_size,
            allow_edge_extrapolation,
        )
        frequency_value = local_pchip(
            beta_rows,
            target_key[0],
            3,
            stencil_size,
            allow_edge_extrapolation,
        )
        alpha_real_value = local_pchip(
            beta_rows,
            target_key[0],
            6,
            stencil_size,
            allow_edge_extrapolation,
        )
        if growth_value is None or frequency_value is None or alpha_real_value is None:
            repair_log.append(
                {
                    "R": target_key[0],
                    "beta": target_key[1],
                    "method": "unresolved-no-two-sided-stencil",
                    "original_growth": original_row[2],
                    "repaired_growth": float("nan"),
                    "original_omega": original_row[3],
                    "repaired_omega": float("nan"),
                    "lower_R": float("nan"),
                    "upper_R": float("nan"),
                }
            )
            continue

        growth, growth_method, lower_radius, upper_radius = growth_value
        frequency, frequency_method, _, _ = frequency_value
        alpha_real, alpha_method, _, _ = alpha_real_value
        row_index = next(
            index
            for index, row in enumerate(repaired)
            if point_key(row[0], row[1]) == target_key
        )
        repaired[row_index, 2] = growth
        repaired[row_index, 3] = frequency
        repaired[row_index, 4] = -1.0
        repaired[row_index, 5] = target_key[0] * frequency
        repaired[row_index, 6] = alpha_real
        repaired[row_index, 7] = -growth
        repaired[row_index, 10] = -1.0
        method = growth_method
        if frequency_method != method or alpha_method != method:
            method = "fixed-beta-pchip"
        repair_log.append(
            {
                "R": target_key[0],
                "beta": target_key[1],
                "method": method,
                "original_growth": original_row[2],
                "repaired_growth": growth,
                "original_omega": original_row[3],
                "repaired_omega": frequency,
                "lower_R": lower_radius,
                "upper_R": upper_radius,
            }
        )
    return repaired, repair_log


def format_row(row: np.ndarray) -> str:
    return " ".join(f"{value:.12e}" for value in row)


def write_envelope(
    source_path: Path,
    output_path: Path,
    envelope: np.ndarray,
    repaired_count: int,
    unresolved_count: int,
) -> None:
    source_lines = source_path.read_text(encoding="utf-8").splitlines()
    zone_index = next(
        index for index, line in enumerate(source_lines) if line.startswith("ZONE")
    )
    metadata = [
        'DATASETAUXDATA repair_method="fixed-beta-local-PCHIP"',
        f'DATASETAUXDATA repaired_points="{repaired_count}"',
        f'DATASETAUXDATA unresolved_points="{unresolved_count}"',
        'DATASETAUXDATA original_values_preserved="yes; repaired nodes have branch_id=-1 and residual=-1"',
    ]
    output_lines = source_lines[:zone_index] + metadata + source_lines[zone_index : zone_index + 1]
    output_lines.extend(format_row(row) for row in envelope)
    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def patch_fetriangle(
    source_path: Path,
    output_path: Path,
    repaired_by_key: dict[tuple[float, float], np.ndarray],
    repaired_count: int,
    unresolved_count: int,
) -> None:
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
        'DATASETAUXDATA repair_method="fixed-beta-local-PCHIP on envelope nodes"',
        f'DATASETAUXDATA repaired_envelope_points="{repaired_count}"',
        f'DATASETAUXDATA unresolved_envelope_points="{unresolved_count}"',
        f'DATASETAUXDATA replaced_FE_nodes="{replaced_nodes}"',
    ]
    patched_lines = patched_lines[:zone_index] + metadata + patched_lines[zone_index:]
    output_path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")


def write_repair_log(path: Path, repair_log: list[dict[str, float | str]]) -> None:
    columns = (
        "R",
        "beta",
        "method",
        "lower_R",
        "upper_R",
        "original_growth",
        "repaired_growth",
        "original_omega",
        "repaired_omega",
    )
    lines = [" ".join(columns)]
    for entry in repair_log:
        lines.append(
            f"{float(entry['R']):.12e} {float(entry['beta']):.12e} "
            f"{entry['method']} {float(entry['lower_R']):.12e} "
            f"{float(entry['upper_R']):.12e} {float(entry['original_growth']):.12e} "
            f"{float(entry['repaired_growth']):.12e} {float(entry['original_omega']):.12e} "
            f"{float(entry['repaired_omega']):.12e}"
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
    parser.add_argument("--allow-edge-extrapolation", action="store_true")
    parser.add_argument("--max-iterations", type=int, default=5)
    args = parser.parse_args()
    if args.stencil_size < 2:
        raise ValueError("stencil-size must be at least 2")

    result_directory = args.result_directory.resolve()
    output_directory = (
        args.output_directory.resolve()
        if args.output_directory is not None
        else result_directory / "interpolated_first_mode_R"
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    envelope_path = result_directory / "spatial_growth_envelope.dat"
    fetriangle_path = result_directory / "spatial_growth_cloud_fetriangle.dat"
    samples_path = result_directory / "spatial_growth_samples.dat"
    envelope = numeric_rows(envelope_path, 11)
    if envelope.size == 0:
        raise RuntimeError("No envelope rows found")

    flagged, initial_thresholds = detect_flags(envelope)
    targets = {
        key
        for key in flagged
        if args.r_min <= key[0] <= args.r_max
        and args.beta_min <= key[1] <= args.beta_max
    }
    initial_target_count = len(targets)
    iteration_count = 0
    final_thresholds = initial_thresholds
    for iteration_index in range(args.max_iterations + 1):
        iteration_count = iteration_index + 1
        repaired, repair_log = repair_envelope(
            envelope,
            targets,
            args.stencil_size,
            args.allow_edge_extrapolation,
        )
        post_flags, final_thresholds = detect_flags(repaired)
        post_targets = {
            key
            for key in post_flags
            if args.r_min <= key[0] <= args.r_max
            and args.beta_min <= key[1] <= args.beta_max
        }
        new_targets = post_targets - targets
        if not new_targets or iteration_index == args.max_iterations:
            break
        targets.update(new_targets)
    repaired_count = sum(
        1
        for entry in repair_log
        if entry["method"] != "unresolved-no-two-sided-stencil"
    )
    unresolved_count = len(repair_log) - repaired_count
    repaired_by_key = {
        point_key(row[0], row[1]): row for row in repaired
    }

    output_envelope = output_directory / "spatial_growth_envelope.dat"
    output_fetriangle = output_directory / "spatial_growth_cloud_fetriangle.dat"
    output_log = output_directory / "interpolated_points.dat"
    output_targets = output_directory / "raw_flagged_first_mode_points.dat"
    write_envelope(
        envelope_path,
        output_envelope,
        repaired,
        repaired_count,
        unresolved_count,
    )
    patch_fetriangle(
        fetriangle_path,
        output_fetriangle,
        {
            key: repaired_by_key[key]
            for key in targets
            if any(
                point_key(entry["R"], entry["beta"]) == key
                and entry["method"] != "unresolved-no-two-sided-stencil"
                for entry in repair_log
            )
        },
        repaired_count,
        unresolved_count,
    )
    write_repair_log(output_log, repair_log)
    output_targets.write_text(
        "\n".join(f"{radius:.12e} {beta:.12e}" for radius, beta in sorted(targets))
        + "\n",
        encoding="utf-8",
    )
    if samples_path.exists():
        shutil.copy2(samples_path, output_directory / "spatial_growth_samples.dat")

    summary_lines = [
        "Fixed-beta spatial-growth repair",
        f"source: {result_directory}",
        f"window: R=[{args.r_min}, {args.r_max}], beta=[{args.beta_min}, {args.beta_max}]",
        f"initial flagged points in window: {initial_target_count}",
        f"repair target points in window: {len(targets)}",
        f"repaired points: {repaired_count}",
        f"unresolved points: {unresolved_count}",
        f"stencil size: {args.stencil_size}",
        f"edge extrapolation enabled: {args.allow_edge_extrapolation}",
        f"iterations: {iteration_count}",
        f"post-repair flagged points in window: {len(post_targets)}",
        f"post-repair flagged points already repaired: {len(post_targets & targets)}",
        "initial thresholds: " + " ".join(f"{value:.12e}" for value in initial_thresholds),
        "final thresholds: " + " ".join(f"{value:.12e}" for value in final_thresholds),
        f"envelope: {output_envelope}",
        f"FE triangle: {output_fetriangle}",
        f"repair log: {output_log}",
    ]
    (output_directory / "repair_report.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
