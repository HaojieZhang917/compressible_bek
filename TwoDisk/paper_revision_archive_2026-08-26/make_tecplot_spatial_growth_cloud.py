from __future__ import annotations

import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = (
    ROOT
    / "zarf_spatial_growth_results"
    / "Ts_p0p0_N99"
    / "spatial_growth_envelope.dat"
)
OUTPUT = (
    ROOT
    / "zarf_spatial_growth_results"
    / "Ts_p0p0_N99"
    / "spatial_growth_cloud_fequad.dat"
)


def read_rows(path: Path) -> list[tuple[float, ...]]:
    rows: list[tuple[float, ...]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if len(values) not in (10, 11):
            continue
        try:
            parsed = tuple(map(float, values))
        except ValueError:
            continue
        rows.append(parsed if len(parsed) == 11 else parsed[:4] + (0.0,) + parsed[4:])
    if not rows:
        raise RuntimeError(f"No numeric rows found in {path}")
    return rows


def key(radius: float, beta: float) -> tuple[float, float]:
    return round(radius, 10), round(beta, 10)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    rows = read_rows(args.source)
    rows.sort(key=lambda row: (row[0], row[1]))
    radii = sorted({row[0] for row in rows})
    betas = sorted({row[1] for row in rows})
    node_index = {key(row[0], row[1]): index + 1 for index, row in enumerate(rows)}

    elements: list[tuple[int, int, int, int]] = []
    for radius_index in range(len(radii) - 1):
        for beta_index in range(len(betas) - 1):
            corners = (
                key(radii[radius_index], betas[beta_index]),
                key(radii[radius_index], betas[beta_index + 1]),
                key(radii[radius_index + 1], betas[beta_index + 1]),
                key(radii[radius_index + 1], betas[beta_index]),
            )
            if all(corner in node_index for corner in corners):
                elements.append(tuple(node_index[corner] for corner in corners))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        stream.write('TITLE="Spatial growth-rate cloud inside the temporal Zarf boundary"\n')
        stream.write(
            'VARIABLES="R","beta","growth_rate_max","omega_opt","branch_id","omega_bar_opt",'
            '"alpha_r","alpha_i","alpha_temporal","omega_i_temporal","residual"\n'
        )
        stream.write('DATASETAUXDATA mass_flux = "0.0"\n')
        stream.write('DATASETAUXDATA N_cheb = "99"\n')
        stream.write(
            'DATASETAUXDATA description = '
            '"Maximum direct spatial growth among the available temporal-frequency samples"\n'
        )
        stream.write(
            f'ZONE T="a_s=0.0", N={len(rows)}, E={len(elements)}, '
            'ZONETYPE=FEQUADRILATERAL, DATAPACKING=POINT\n'
        )
        for row in rows:
            stream.write(" ".join(f"{value:.12e}" for value in row) + "\n")
        for element in elements:
            stream.write(" ".join(map(str, element)) + "\n")

    print(f"output: {args.output}")
    print(f"nodes: {len(rows)}")
    print(f"quadrilateral elements: {len(elements)}")
    print("contour variable: growth_rate_max")


if __name__ == "__main__":
    main()
