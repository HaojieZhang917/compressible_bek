from __future__ import annotations

import argparse
from pathlib import Path


def numeric_rows(path: Path, columns: int) -> list[tuple[float, ...]]:
    rows: list[tuple[float, ...]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if len(values) != columns:
            continue
        try:
            rows.append(tuple(map(float, values)))
        except ValueError:
            continue
    return rows


def write_samples(
    path: Path,
    rows: list[tuple[float, ...]],
    mass_flux: float,
    frequency_samples: int,
) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write('TITLE="Frequency-resolved spatial growth inside temporal Zarf"\n')
        stream.write(
            'VARIABLES="R","beta","omega","omega_bar","branch_id",'
            '"alpha_temporal","omega_i_temporal","alpha_r","alpha_i",'
            '"growth_rate","residual","status"\n'
        )
        stream.write(f'DATASETAUXDATA mass_flux="{mass_flux:.1f}"\n')
        stream.write('DATASETAUXDATA N_cheb="99"\n')
        stream.write(f'DATASETAUXDATA frequency_samples="{frequency_samples}"\n')
        stream.write(f'ZONE T="a_s={mass_flux:.1f}", I={len(rows)}, F=POINT\n')
        for row in rows:
            values = [f"{value:.12e}" for value in row]
            values[4] = str(round(row[4]))
            values[-1] = str(round(row[-1]))
            stream.write(" ".join(values) + "\n")


def write_envelope(
    path: Path,
    rows: list[tuple[float, ...]],
    mass_flux: float,
    frequency_samples: int,
) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write('TITLE="Maximum local spatial growth inside temporal Zarf"\n')
        stream.write(
            'VARIABLES="R","beta","growth_rate_max","omega_opt","branch_id",'
            '"omega_bar_opt","alpha_r","alpha_i","alpha_temporal",'
            '"omega_i_temporal","residual"\n'
        )
        stream.write(f'DATASETAUXDATA mass_flux="{mass_flux:.1f}"\n')
        stream.write('DATASETAUXDATA N_cheb="99"\n')
        stream.write(f'DATASETAUXDATA frequency_samples="{frequency_samples}"\n')
        stream.write(f'ZONE T="a_s={mass_flux:.1f}", I={len(rows)}, F=POINT\n')
        for row in rows:
            values = [f"{value:.12e}" for value in row]
            values[4] = str(round(row[4]))
            stream.write(" ".join(values) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("chunk_root", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--case-directory", default="Ts_p0p0_N99")
    parser.add_argument("--mass-flux", type=float, default=0.0)
    parser.add_argument("--frequency-samples", type=int, default=13)
    args = parser.parse_args()

    sample_paths = sorted(
        args.chunk_root.glob(f"*/{args.case_directory}/spatial_growth_samples.dat")
    )
    envelope_paths = sorted(
        args.chunk_root.glob(f"*/{args.case_directory}/spatial_growth_envelope.dat")
    )
    if not sample_paths or not envelope_paths:
        raise RuntimeError(f"No completed chunk files found below {args.chunk_root}")

    samples = [row for path in sample_paths for row in numeric_rows(path, 12)]
    envelopes = [row for path in envelope_paths for row in numeric_rows(path, 11)]
    samples.sort(key=lambda row: (row[0], row[1], row[2], row[4]))
    envelopes.sort(key=lambda row: (row[0], row[1]))

    seen: set[tuple[float, float]] = set()
    for row in envelopes:
        key = round(row[0], 10), round(row[1], 10)
        if key in seen:
            raise RuntimeError(f"Duplicate envelope point: {key}")
        seen.add(key)

    args.output_directory.mkdir(parents=True, exist_ok=True)
    write_samples(
        args.output_directory / "spatial_growth_samples.dat",
        samples,
        args.mass_flux,
        args.frequency_samples,
    )
    write_envelope(
        args.output_directory / "spatial_growth_envelope.dat",
        envelopes,
        args.mass_flux,
        args.frequency_samples,
    )
    print(f"chunks: {len(envelope_paths)}")
    print(f"samples: {len(samples)}")
    print(f"envelope points: {len(envelopes)}")


if __name__ == "__main__":
    main()
