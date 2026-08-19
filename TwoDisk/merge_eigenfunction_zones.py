"""Merge the R1-C12 eigenfunction profiles into one Tecplot data file."""

from __future__ import annotations

import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "eigenfunction_TypeI_TypeII_mass_flux.dat"

CASES = [
    {
        "mode": "Type I",
        "tag": "Type1",
        "a_s": a_s,
        "R": 330.0,
        "n": 16,
        "omega_bar": 8.0,
    }
    for a_s in (-0.4, 0.0, 0.4)
] + [
    {
        "mode": "Type II",
        "tag": "Type2",
        "a_s": a_s,
        "R": 100.0,
        "n": -7,
        "omega_bar": 8.0,
    }
    for a_s in (-0.4, 0.0, 0.4)
]


def julia_complex(text: str) -> complex:
    """Parse a scalar written by Julia, for example ``1.0 - 2.0im``."""
    compact = re.sub(r"\s+", "", text.strip()).replace("im", "j")
    return complex(compact)


def read_profile(path: Path) -> list[tuple[complex, ...]]:
    rows: list[tuple[complex, ...]] = []
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            columns = line.rstrip().split("\t")
            if len(columns) != 5:
                raise ValueError(
                    f"{path.name}:{line_number}: expected 5 columns, "
                    f"found {len(columns)}"
                )
            rows.append(tuple(julia_complex(value) for value in columns))
    if not rows:
        raise ValueError(f"{path.name}: empty input file")
    return rows


def source_path(case: dict) -> Path:
    # The source names intentionally have no .dat suffix because the decimal
    # part of a_s was interpreted as the file extension by the notebook.
    return ROOT / f"eigfunction_{case['tag']}_as={case['a_s']:.1f}"


def main() -> None:
    profiles = []
    for case in CASES:
        path = source_path(case)
        if not path.exists():
            raise FileNotFoundError(path)
        rows = read_profile(path)
        z = [row[0].real for row in rows]
        if max(abs(row[0].imag) for row in rows) > 1.0e-13:
            raise ValueError(f"{path.name}: z contains a non-zero imaginary part")
        profiles.append((case, path, rows))

    with OUTPUT.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(
            'TITLE = "Type-I and Type-II direct eigenfunction amplitudes '
            'under wall mass flux"\n'
        )
        stream.write(
            'VARIABLES = "z" "abs_u" "abs_v" "abs_w" "abs_p" '
            '"abs_velocity"\n'
        )

        for case, path, rows in profiles:
            stream.write(
                f'ZONE T="{case["mode"]}, a_s={case["a_s"]:+.1f}", '
                f'I={len(rows)}, F=POINT\n'
            )
            stream.write(f'AUXDATA SourceFile="{path.name}"\n')
            stream.write(f'AUXDATA Mode="{case["mode"]}"\n')
            stream.write(f'AUXDATA a_s="{case["a_s"]:+.1f}"\n')
            stream.write(f'AUXDATA R="{case["R"]:.1f}"\n')
            stream.write(f'AUXDATA n="{case["n"]}"\n')
            stream.write(f'AUXDATA omega_bar="{case["omega_bar"]:.1f}"\n')
            stream.write('AUXDATA Normalization="max(abs(v))=1"\n')

            for z_value, u, v, w, p in rows:
                velocity = math.sqrt(abs(u) ** 2 + abs(v) ** 2 + abs(w) ** 2)
                stream.write(
                    f"{z_value.real:.15e} {abs(u):.15e} {abs(v):.15e} "
                    f"{abs(w):.15e} {abs(p):.15e} {velocity:.15e}\n"
                )

    print(f"Wrote {OUTPUT}")
    print(f"Zones: {len(profiles)}")
    for case, _, rows in profiles:
        max_v = max(abs(row[2]) for row in rows)
        print(
            f"{case['mode']:7s} a_s={case['a_s']:+.1f}: "
            f"points={len(rows)}, max|v|={max_v:.12g}"
        )


if __name__ == "__main__":
    main()
