from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CAVITY_SOURCE = (
    ROOT
    / "configuration_total_amplitude_standard_results"
    / "total_amplitude_comparison_N99_standard.dat"
)
DISK_SOURCE = ROOT / "vonkarman_total_amplitude_results" / "vonkarman_A_N99.dat"
OUTPUT_ROOT = ROOT / "vonkarman_total_amplitude_results"


def read_zone(path: Path, zone_name: str) -> np.ndarray:
    rows: list[list[float]] = []
    active = False
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if text.startswith("ZONE"):
            active = f'T="{zone_name}"' in text
            continue
        if not active or not text:
            continue
        try:
            rows.append([float(value) for value in text.split()])
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"Zone {zone_name!r} not found in {path}")
    return np.asarray(rows)


def crossover_radii(cavity: np.ndarray, disk: np.ndarray) -> list[float]:
    radii = disk[:, 0]
    cavity_amplitude = np.interp(radii, cavity[:, 0], cavity[:, 7])
    difference = disk[:, 8] - cavity_amplitude
    crossings: list[float] = []
    for index in range(1, len(radii)):
        if difference[index - 1] * difference[index] > 0:
            continue
        fraction = -difference[index - 1] / (difference[index] - difference[index - 1])
        crossings.append(radii[index - 1] + fraction * (radii[index] - radii[index - 1]))
    return crossings


def main() -> None:
    cavity = read_zone(CAVITY_SOURCE, "rotor-stator cavity")
    disk = read_zone(DISK_SOURCE, "isolated von Karman disk")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.7,
        }
    )
    figure, axis = plt.subplots(figsize=(6.4, 3.6))
    axis.semilogy(cavity[:, 0], cavity[:, 7], color="black", label="Rotor--stator")
    axis.semilogy(
        disk[:, 0], disk[:, 8], color="#b9474d", linestyle="--", label="Isolated disk"
    )
    axis.plot(cavity[0, 0], cavity[0, 7], "o", color="black", markersize=3.5)
    axis.plot(disk[0, 0], disk[0, 8], "o", color="#b9474d", markersize=3.5)
    axis.set_xlim(280.0, 500.0)
    axis.set_ylim(0.1, 500.0)
    axis.set_xlabel(r"$R$")
    axis.set_ylabel(r"$|A(R;R_l)|$")
    axis.tick_params(direction="in", top=True, right=True)
    axis.legend(frameon=False, loc="upper left")
    figure.tight_layout()

    for extension in ("eps", "pdf", "png"):
        figure.savefig(
            OUTPUT_ROOT / f"fig21_amplitude_comparison.{extension}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)

    crossings = crossover_radii(cavity, disk)
    summary = OUTPUT_ROOT / "fig21_amplitude_comparison_summary.txt"
    summary.write_text(
        "\n".join(
            (
                f"Cavity lower neutral R = {cavity[0, 0]:.12f}",
                f"Cavity initial A = {cavity[0, 7]:.12e}",
                f"Cavity A(500) = {cavity[-1, 7]:.12e}",
                f"Disk lower neutral R = {disk[0, 0]:.12f}",
                f"Disk initial A = {disk[0, 8]:.12e}",
                f"Disk A(500) = {disk[-1, 8]:.12e}",
                "Crossover R = " + ", ".join(f"{radius:.12f}" for radius in crossings),
                f"A(500) ratio disk/cavity = {disk[-1, 8] / cavity[-1, 7]:.12e}",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
