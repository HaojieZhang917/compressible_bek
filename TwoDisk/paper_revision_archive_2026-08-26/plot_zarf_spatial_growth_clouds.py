from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import Delaunay


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


def smooth_boundary(envelope: np.ndarray, zarf: np.ndarray, count: int = 800) -> np.ndarray:
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
    return np.vstack(
        (
            np.column_stack((dense_radius, upper_beta)),
            np.column_stack((dense_radius[::-1], lower_beta[::-1])),
        )
    )


def inside_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    inside = np.zeros(len(points), dtype=bool)
    x0, y0 = polygon[-1]
    for x1, y1 in polygon:
        crossing = (y1 > y) != (y0 > y)
        denominator = y0 - y1
        if abs(denominator) <= np.finfo(float).eps:
            denominator = np.finfo(float).eps
        x_cross = (x0 - x1) * (y - y1) / denominator + x1
        inside ^= crossing & (x < x_cross)
        x0, y0 = x1, y1
    return inside


def clipped_triangulation(envelope: np.ndarray, boundary: np.ndarray) -> mtri.Triangulation:
    radii = np.unique(envelope[:, 0])
    betas = np.unique(envelope[:, 1])
    radius_scale = float(np.median(np.diff(radii)))
    beta_scale = float(np.median(np.diff(betas)))
    scaled = np.column_stack((envelope[:, 0] / radius_scale, envelope[:, 1] / beta_scale))
    triangles = Delaunay(scaled).simplices
    centroids = np.mean(envelope[triangles, :2], axis=1)
    triangle_points = scaled[triangles]
    edge_lengths = np.stack(
        (
            np.linalg.norm(triangle_points[:, 0] - triangle_points[:, 1], axis=1),
            np.linalg.norm(triangle_points[:, 1] - triangle_points[:, 2], axis=1),
            np.linalg.norm(triangle_points[:, 2] - triangle_points[:, 0], axis=1),
        ),
        axis=1,
    )
    keep = inside_polygon(centroids, boundary) & (np.max(edge_lengths, axis=1) <= 4.0)
    return mtri.Triangulation(envelope[:, 0], envelope[:, 1], triangles[keep])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--case",
        nargs=3,
        action="append",
        metavar=("MASS_FLUX", "ENVELOPE", "BOUNDARY"),
        required=True,
    )
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    args = parser.parse_args()

    cases = []
    for mass_flux, envelope_path, boundary_path in args.case:
        envelope = numeric_rows(Path(envelope_path), 11)
        boundary = smooth_boundary(envelope, numeric_rows(Path(boundary_path), 2))
        cases.append((float(mass_flux), envelope, boundary))

    growth = np.concatenate([case[1][:, 2] for case in cases])
    vmin = float(np.min(growth)) if args.vmin is None else args.vmin
    vmax = float(np.max(growth)) if args.vmax is None else args.vmax
    levels = np.linspace(vmin, vmax, 41)
    figure, axes = plt.subplots(
        1,
        len(cases),
        figsize=(4.1 * len(cases), 4.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    contour = None
    for axis, (mass_flux, envelope, boundary) in zip(axes, cases):
        triangulation = clipped_triangulation(envelope, boundary)
        contour = axis.tricontourf(
            triangulation,
            envelope[:, 2],
            levels=levels,
            cmap="viridis",
            extend="both",
        )
        axis.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0)
        axis.set_title(rf"$a_s={mass_flux:.1f}$")
        axis.set_xlabel(r"$R$")
        axis.tick_params(direction="in")
        axis.set_facecolor("white")
    axes[0].set_ylabel(r"$\beta$")
    if contour is not None:
        figure.colorbar(contour, ax=axes, label=r"$-\alpha_i^{\max}$", shrink=0.92)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220)
    plt.close(figure)
    print(f"output: {args.output}")
    print(f"color range: {vmin:.12e} {vmax:.12e}")


if __name__ == "__main__":
    main()
