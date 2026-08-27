from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr

from repair_zarf_spatial_growth_fixed_beta import (
    detect_flags,
    format_row,
    numeric_rows,
    point_key,
)


def select_forward_mask(
    envelope: np.ndarray,
    flagged: set[tuple[float, float]],
    r_min: float,
    r_max: float,
    beta_min: float,
    beta_max: float,
    growth_jump_threshold: float,
) -> tuple[set[tuple[float, float]], dict[float, float]]:
    start_by_beta: dict[float, float] = {}
    for radius, beta in flagged:
        if not (r_min <= radius <= r_max and beta_min <= beta <= beta_max):
            continue
        beta_key = round(beta, 10)
        start_by_beta[beta_key] = min(start_by_beta.get(beta_key, radius), radius)

    rows_by_beta: dict[float, list[np.ndarray]] = {}
    for row in envelope:
        radius = float(row[0])
        beta_key = round(float(row[1]), 10)
        if r_min <= radius <= r_max and beta_min <= beta_key <= beta_max:
            rows_by_beta.setdefault(beta_key, []).append(row)
    for beta_key, rows in rows_by_beta.items():
        ordered_rows = sorted(rows, key=lambda row: row[0])
        for lower_row, upper_row in zip(ordered_rows, ordered_rows[1:]):
            if abs(upper_row[2] - lower_row[2]) < growth_jump_threshold:
                continue
            start_by_beta[beta_key] = min(
                start_by_beta.get(beta_key, lower_row[0]),
                float(lower_row[0]),
            )

    targets: set[tuple[float, float]] = set()
    for row in envelope:
        radius = float(row[0])
        beta_key = round(float(row[1]), 10)
        start_radius = start_by_beta.get(beta_key)
        if start_radius is None:
            continue
        if r_min <= radius <= r_max and radius >= start_radius:
            targets.add(point_key(radius, beta_key))
    return targets, start_by_beta


def build_constraints(
    envelope: np.ndarray,
    targets: set[tuple[float, float]],
    curvature_weight: float,
    edge_weight: float,
) -> tuple[dict[tuple[float, float], int], list[tuple[list[tuple[tuple[float, float], float]], float]]]:
    by_key = {
        point_key(row[0], row[1]): row for row in envelope
    }
    radii = sorted({key[0] for key in by_key})
    betas = sorted({key[1] for key in by_key})
    radius_step = float(np.median(np.diff(radii)))
    beta_step = float(np.median(np.diff(betas)))
    unknown_keys = sorted(targets)
    unknown_index = {key: index for index, key in enumerate(unknown_keys)}
    constraints: list[tuple[list[tuple[tuple[float, float], float]], float]] = []

    def add_constraint(coefficients: list[tuple[tuple[float, float], float]], weight: float) -> None:
        if any(key in unknown_index for key, _ in coefficients):
            constraints.append((coefficients, weight))

    for key in by_key:
        radius, beta = key
        right_key = point_key(radius + radius_step, beta)
        upper_key = point_key(radius, beta + beta_step)
        if right_key in by_key:
            add_constraint([(key, -1.0), (right_key, 1.0)], edge_weight)
        if upper_key in by_key:
            add_constraint([(key, -1.0), (upper_key, 1.0)], edge_weight)

        left_key = point_key(radius - radius_step, beta)
        if left_key in by_key and right_key in by_key:
            add_constraint(
                [(left_key, 1.0), (key, -2.0), (right_key, 1.0)],
                curvature_weight,
            )
        lower_key = point_key(radius, beta - beta_step)
        if lower_key in by_key and upper_key in by_key:
            add_constraint(
                [(lower_key, 1.0), (key, -2.0), (upper_key, 1.0)],
                curvature_weight,
            )

    return unknown_index, constraints


def solve_field(
    envelope_by_key: dict[tuple[float, float], np.ndarray],
    unknown_index: dict[tuple[float, float], int],
    constraints: list[tuple[list[tuple[tuple[float, float], float]], float]],
    column: int,
) -> tuple[np.ndarray, float]:
    row_indices: list[int] = []
    column_indices: list[int] = []
    matrix_values: list[float] = []
    right_hand_side: list[float] = []
    equation_index = 0
    for coefficients, weight in constraints:
        known_contribution = 0.0
        for key, coefficient in coefficients:
            if key in unknown_index:
                row_indices.append(equation_index)
                column_indices.append(unknown_index[key])
                matrix_values.append(coefficient * np.sqrt(weight))
            else:
                known_contribution += coefficient * envelope_by_key[key][column]
        right_hand_side.append(-known_contribution * np.sqrt(weight))
        equation_index += 1
    matrix = coo_matrix(
        (matrix_values, (row_indices, column_indices)),
        shape=(equation_index, len(unknown_index)),
    ).tocsr()
    solution = lsmr(matrix, np.asarray(right_hand_side), atol=1e-12, btol=1e-12, maxiter=20000)
    return solution[0], float(solution[3])


def write_envelope(
    source_path: Path,
    output_path: Path,
    repaired: np.ndarray,
    target_count: int,
    curvature_weight: float,
    edge_weight: float,
) -> None:
    source_lines = source_path.read_text(encoding="utf-8").splitlines()
    zone_index = next(index for index, line in enumerate(source_lines) if line.startswith("ZONE"))
    metadata = [
        'DATASETAUXDATA repair_method="two-dimensional curvature-constrained inpainting"',
        f'DATASETAUXDATA repair_targets="{target_count}"',
        f'DATASETAUXDATA curvature_weight="{curvature_weight:.6e}"',
        f'DATASETAUXDATA edge_weight="{edge_weight:.6e}"',
        'DATASETAUXDATA original_values_preserved="yes; repaired nodes have branch_id=-2 and residual=-2"',
    ]
    output_lines = source_lines[:zone_index] + metadata + source_lines[zone_index : zone_index + 1]
    output_lines.extend(format_row(row) for row in repaired)
    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def patch_fetriangle(
    source_path: Path,
    output_path: Path,
    repaired_by_key: dict[tuple[float, float], np.ndarray],
    target_count: int,
    curvature_weight: float,
    edge_weight: float,
) -> int:
    source_lines = source_path.read_text(encoding="utf-8").splitlines()
    zone_index = next(index for index, line in enumerate(source_lines) if "ZONETYPE=FETRIANGLE" in line)
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
        'DATASETAUXDATA repair_method="two-dimensional curvature-constrained inpainting on envelope nodes"',
        f'DATASETAUXDATA repair_targets="{target_count}"',
        f'DATASETAUXDATA curvature_weight="{curvature_weight:.6e}"',
        f'DATASETAUXDATA edge_weight="{edge_weight:.6e}"',
        f'DATASETAUXDATA replaced_FE_nodes="{replaced_nodes}"',
    ]
    patched_lines = patched_lines[:zone_index] + metadata + patched_lines[zone_index:]
    output_path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")
    return replaced_nodes


def write_repair_log(
    path: Path,
    original_by_key: dict[tuple[float, float], np.ndarray],
    repaired_by_key: dict[tuple[float, float], np.ndarray],
    targets: set[tuple[float, float]],
) -> None:
    lines = [
        "R beta original_growth repaired_growth growth_delta original_omega repaired_omega omega_delta original_alpha_r repaired_alpha_r alpha_r_delta"
    ]
    for key in sorted(targets):
        original = original_by_key[key]
        repaired = repaired_by_key[key]
        lines.append(
            f"{key[0]:.12e} {key[1]:.12e} {original[2]:.12e} {repaired[2]:.12e} "
            f"{repaired[2] - original[2]:.12e} {original[3]:.12e} {repaired[3]:.12e} "
            f"{repaired[3] - original[3]:.12e} {original[6]:.12e} {repaired[6]:.12e} "
            f"{repaired[6] - original[6]:.12e}"
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
    parser.add_argument("--curvature-weight", type=float, default=1.0)
    parser.add_argument("--edge-weight", type=float, default=0.05)
    parser.add_argument("--growth-jump-threshold", type=float, default=0.015)
    args = parser.parse_args()

    result_directory = args.result_directory.resolve()
    output_directory = (
        args.output_directory.resolve()
        if args.output_directory is not None
        else result_directory / "inpainted_2d_first_mode_R"
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    envelope_path = result_directory / "spatial_growth_envelope.dat"
    fetriangle_path = result_directory / "spatial_growth_cloud_fetriangle.dat"
    samples_path = result_directory / "spatial_growth_samples.dat"
    envelope = numeric_rows(envelope_path, 11)
    if envelope.size == 0:
        raise RuntimeError("No envelope rows found")

    flagged, thresholds = detect_flags(envelope)
    targets, starts = select_forward_mask(
        envelope,
        flagged,
        args.r_min,
        args.r_max,
        args.beta_min,
        args.beta_max,
        args.growth_jump_threshold,
    )
    envelope_by_key = {point_key(row[0], row[1]): row for row in envelope}
    unknown_index, constraints = build_constraints(
        envelope,
        targets,
        args.curvature_weight,
        args.edge_weight,
    )
    repaired = envelope.copy()
    solver_residuals: dict[str, float] = {}
    for column, label in ((2, "growth"), (3, "omega"), (6, "alpha_r")):
        values, solver_residual = solve_field(
            envelope_by_key,
            unknown_index,
            constraints,
            column,
        )
        solver_residuals[label] = solver_residual
        for key, index in unknown_index.items():
            row_index = next(
                row_index
                for row_index, row in enumerate(repaired)
                if point_key(row[0], row[1]) == key
            )
            repaired[row_index, column] = values[index]
            repaired[row_index, 4] = -2.0
            repaired[row_index, 10] = -2.0

    for row_index, row in enumerate(repaired):
        key = point_key(row[0], row[1])
        if key not in targets:
            continue
        repaired[row_index, 5] = row[0] * row[3]
        repaired[row_index, 7] = -row[2]

    repaired_by_key = {point_key(row[0], row[1]): row for row in repaired}
    output_envelope = output_directory / "spatial_growth_envelope.dat"
    output_fetriangle = output_directory / "spatial_growth_cloud_fetriangle.dat"
    output_targets = output_directory / "inpaint_targets.dat"
    output_log = output_directory / "inpainted_points.dat"
    write_envelope(
        envelope_path,
        output_envelope,
        repaired,
        len(targets),
        args.curvature_weight,
        args.edge_weight,
    )
    replaced_nodes = patch_fetriangle(
        fetriangle_path,
        output_fetriangle,
        {key: repaired_by_key[key] for key in targets},
        len(targets),
        args.curvature_weight,
        args.edge_weight,
    )
    output_targets.write_text(
        "\n".join(f"{radius:.12e} {beta:.12e}" for radius, beta in sorted(targets))
        + "\n",
        encoding="utf-8",
    )
    write_repair_log(output_log, envelope_by_key, repaired_by_key, targets)
    if samples_path.exists():
        shutil.copy2(samples_path, output_directory / "spatial_growth_samples.dat")

    start_lines = [f"beta={beta:.12e} start_R={radius:.12e}" for beta, radius in sorted(starts.items())]
    report_lines = [
        "Two-dimensional spatial-growth inpainting",
        f"source: {result_directory}",
        f"window: R=[{args.r_min}, {args.r_max}], beta=[{args.beta_min}, {args.beta_max}]",
        f"continuation beta lines: {len(starts)}",
        f"target points: {len(targets)}",
        f"replaced FE nodes: {replaced_nodes}",
        f"constraints: {len(constraints)}",
        f"growth jump threshold: {args.growth_jump_threshold:.6e}",
        f"solver residual growth={solver_residuals['growth']:.6e}",
        f"solver residual omega={solver_residuals['omega']:.6e}",
        f"solver residual alpha_r={solver_residuals['alpha_r']:.6e}",
        "thresholds: " + " ".join(f"{value:.12e}" for value in thresholds),
        "continuation starts:",
        *start_lines,
        f"envelope: {output_envelope}",
        f"FE triangle: {output_fetriangle}",
        f"repair log: {output_log}",
    ]
    (output_directory / "inpaint_report.txt").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )
    print("\n".join(report_lines))


if __name__ == "__main__":
    main()
