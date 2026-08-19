from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, PchipInterpolator
from scipy.spatial import Delaunay


VARIABLES = (
    "R",
    "beta",
    "growth_rate_max",
    "omega_opt",
    "branch_id",
    "omega_bar_opt",
    "alpha_r",
    "alpha_i",
    "alpha_temporal",
    "omega_i_temporal",
    "residual",
)


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


def unique_branch(points: np.ndarray) -> np.ndarray:
    radii = sorted(set(points[:, 0]))
    return np.asarray(
        [(radius, np.mean(points[points[:, 0] == radius, 1])) for radius in radii]
    )


def smooth_boundary(envelope: np.ndarray, zarf: np.ndarray, count: int) -> np.ndarray:
    radius_values = sorted(set(envelope[:, 0]))
    bounds = np.asarray(
        [
            (
                radius,
                np.min(envelope[envelope[:, 0] == radius, 1]),
                np.max(envelope[envelope[:, 0] == radius, 1]),
            )
            for radius in radius_values
        ]
    )
    clipped_radius = np.clip(zarf[:, 0], bounds[0, 0], bounds[-1, 0])
    lower_edge = PchipInterpolator(bounds[:, 0], bounds[:, 1])(clipped_radius)
    upper_edge = PchipInterpolator(bounds[:, 0], bounds[:, 2])(clipped_radius)
    beta_step = float(np.median(np.diff(sorted(set(envelope[:, 1])))))
    edge_tolerance = 4.0 * beta_step
    upper = unique_branch(zarf[np.abs(zarf[:, 1] - upper_edge) <= edge_tolerance])
    lower = unique_branch(zarf[np.abs(zarf[:, 1] - lower_edge) <= edge_tolerance])
    radius_min = max(np.min(upper[:, 0]), np.min(lower[:, 0]))
    radius_max = min(np.max(upper[:, 0]), np.max(lower[:, 0]))
    dense_radius = np.linspace(radius_min, radius_max, count)
    upper_beta = PchipInterpolator(upper[:, 0], upper[:, 1])(dense_radius)
    lower_beta = PchipInterpolator(lower[:, 0], lower[:, 1])(dense_radius)
    upper_dense = np.column_stack((dense_radius, upper_beta))
    lower_dense = np.column_stack((dense_radius, lower_beta))
    return np.vstack((upper_dense, lower_dense[::-1]))


def inside_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    inside = np.zeros(len(points), dtype=bool)
    x0, y0 = polygon[-1]
    for x1, y1 in polygon:
        crossing = (y1 > y) != (y0 > y)
        denominator = y0 - y1
        denominator = denominator if abs(denominator) > np.finfo(float).eps else np.finfo(float).eps
        x_cross = (x0 - x1) * (y - y1) / denominator + x1
        inside ^= crossing & (x < x_cross)
        x0, y0 = x1, y1
    return inside


def interpolate_boundary(envelope: np.ndarray, boundary: np.ndarray) -> np.ndarray:
    radii = sorted(set(envelope[:, 0]))
    betas = sorted(set(envelope[:, 1]))
    radius_scale = float(np.median(np.diff(radii)))
    beta_scale = float(np.median(np.diff(betas)))
    source_xy = np.column_stack((envelope[:, 0] / radius_scale, envelope[:, 1] / beta_scale))
    target_xy = np.column_stack((boundary[:, 0] / radius_scale, boundary[:, 1] / beta_scale))
    values = np.empty((len(boundary), envelope.shape[1]))
    values[:, :2] = boundary
    for column in range(2, envelope.shape[1]):
        linear = LinearNDInterpolator(source_xy, envelope[:, column])
        nearest = NearestNDInterpolator(source_xy, envelope[:, column])
        interpolated = np.asarray(linear(target_xy), dtype=float)
        missing = ~np.isfinite(interpolated)
        if np.any(missing):
            interpolated[missing] = nearest(target_xy[missing])
        values[:, column] = interpolated
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("envelope", type=Path)
    parser.add_argument("zarf_boundary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mass-flux", type=float, required=True)
    parser.add_argument("--boundary-points", type=int, default=500)
    parser.add_argument("--max-scaled-edge", type=float, default=4.0)
    args = parser.parse_args()

    envelope = numeric_rows(args.envelope, 11)
    zarf = numeric_rows(args.zarf_boundary, 2)
    boundary = smooth_boundary(envelope, zarf, args.boundary_points)
    boundary_values = interpolate_boundary(envelope, boundary)

    combined = np.vstack((envelope, boundary_values))
    coordinates = np.round(combined[:, :2], decimals=12)
    _, unique_indices = np.unique(coordinates, axis=0, return_index=True)
    nodes = combined[np.sort(unique_indices)]

    radii = sorted(set(envelope[:, 0]))
    betas = sorted(set(envelope[:, 1]))
    radius_scale = float(np.median(np.diff(radii)))
    beta_scale = float(np.median(np.diff(betas)))
    scaled = np.column_stack((nodes[:, 0] / radius_scale, nodes[:, 1] / beta_scale))
    triangulation = Delaunay(scaled)
    triangles = triangulation.simplices
    centroids = np.mean(nodes[triangles, :2], axis=1)
    keep = inside_polygon(centroids, boundary)
    triangle_points = scaled[triangles]
    edge_lengths = np.stack(
        (
            np.linalg.norm(triangle_points[:, 0] - triangle_points[:, 1], axis=1),
            np.linalg.norm(triangle_points[:, 1] - triangle_points[:, 2], axis=1),
            np.linalg.norm(triangle_points[:, 2] - triangle_points[:, 0], axis=1),
        ),
        axis=1,
    )
    keep &= np.max(edge_lengths, axis=1) <= args.max_scaled_edge
    triangles = triangles[keep]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        stream.write('TITLE="Maximum spatial growth inside temporal Zarf boundary"\n')
        stream.write("VARIABLES=" + ",".join(f'\"{name}\"' for name in VARIABLES) + "\n")
        stream.write(f'DATASETAUXDATA mass_flux="{args.mass_flux:.1f}"\n')
        stream.write('DATASETAUXDATA N_cheb="99"\n')
        stream.write('DATASETAUXDATA boundary="temporal Zarf boundary"\n')
        stream.write(
            'DATASETAUXDATA note="Boundary-node values are interpolated from interior spatial-growth data; they are not forced to zero"\n'
        )
        stream.write(
            f'ZONE T="a_s={args.mass_flux:.1f} spatial growth", N={len(nodes)}, E={len(triangles)}, '
            'ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n'
        )
        for row in nodes:
            stream.write(" ".join(f"{value:.12e}" for value in row) + "\n")
        for triangle in triangles:
            stream.write(" ".join(str(index + 1) for index in triangle) + "\n")

        closed_boundary = np.vstack((boundary_values, boundary_values[0]))
        stream.write(
            f'ZONE T="temporal Zarf boundary", I={len(closed_boundary)}, '
            'DATAPACKING=POINT\n'
        )
        for row in closed_boundary:
            stream.write(" ".join(f"{value:.12e}" for value in row) + "\n")

    print(f"output: {args.output}")
    print(f"raw envelope nodes: {len(envelope)}")
    print(f"total nodes with smooth boundary: {len(nodes)}")
    print(f"triangles: {len(triangles)}")
    print(f"boundary points: {len(boundary)}")


if __name__ == "__main__":
    main()
