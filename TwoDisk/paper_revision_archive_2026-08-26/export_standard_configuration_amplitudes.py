from __future__ import annotations

from pathlib import Path
from typing import Final

import matplotlib.pyplot as plt
import numpy as np


ROOT: Final = Path(__file__).resolve().parent
SOURCE: Final = (
    ROOT / "configuration_total_amplitude_results" / "total_amplitude_comparison_N99.dat"
)
OUTPUT: Final = ROOT / "configuration_total_amplitude_standard_results"
FOURIER_FACTOR: Final = float(np.sqrt(2.0 * np.pi))


def read_zones(path: Path) -> dict[str, np.ndarray]:
    zones: dict[str, list[list[float]]] = {}
    current_name: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if text.startswith("ZONE"):
            current_name = text.split('T="', 1)[1].split('"', 1)[0]
            zones[current_name] = []
            continue
        if current_name is None or not text:
            continue
        try:
            zones[current_name].append([float(value) for value in text.split()])
        except ValueError:
            continue
    return {name: np.asarray(rows) for name, rows in zones.items()}


def standardize(data: np.ndarray) -> np.ndarray:
    result = data.copy()
    result[:, 6] *= FOURIER_FACTOR
    result[:, 7] *= FOURIER_FACTOR
    return result


def write_zone(handle, name: str, data: np.ndarray) -> None:
    handle.write(f'ZONE T="{name}", I={len(data)}, F=POINT\n')
    for row in data:
        handle.write(" ".join(f"{value:.12e}" for value in row) + "\n")


def write_data(cavity: np.ndarray, disk: np.ndarray) -> tuple[Path, Path]:
    comparison_path = OUTPUT / "total_amplitude_comparison_N99_standard.dat"
    disk_path = OUTPUT / "vonkarman_total_amplitude_N99_standard.dat"
    header = (
        'TITLE="Total amplitude with the standard Gaussian Fourier convention"\n'
        'VARIABLES="R","alpha_r","alpha_i","growth_rate","N_factor",'
        '"gain","Cr_initial","A_abs"\n'
    )
    with comparison_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        write_zone(handle, "rotor-stator cavity", cavity)
        write_zone(handle, "isolated rotating disk", disk)
    with disk_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        write_zone(handle, "isolated rotating disk", disk)
    return comparison_path, disk_path


def plot_curves(cavity: np.ndarray, disk: np.ndarray) -> Path:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.7,
        }
    )
    figure, (axis_n, axis_a) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
    for data, color, linestyle, label in (
        (cavity, "black", "-", "Rotor--stator"),
        (disk, "#b9474d", "--", "Isolated disk"),
    ):
        axis_n.plot(data[:, 0], data[:, 4], color=color, linestyle=linestyle, label=label)
        axis_a.semilogy(data[:, 0], data[:, 7], color=color, linestyle=linestyle)
        axis_n.plot(data[0, 0], data[0, 4], "o", color=color, markersize=3.5)
        axis_a.plot(data[0, 0], data[0, 7], "o", color=color, markersize=3.5)

    axis_n.set_ylabel(r"$N(R;R_l)$")
    axis_n.set_ylim(0.0, 8.1)
    axis_n.legend(frameon=False, loc="upper left")
    axis_n.text(0.97, 0.92, "(a)", ha="right", transform=axis_n.transAxes)

    axis_a.set_xlabel(r"$R$")
    axis_a.set_ylabel(r"$|A(R;R_l)|$")
    axis_a.set_xlim(280.0, 505.0)
    axis_a.set_ylim(0.1, 500.0)
    axis_a.text(0.03, 0.92, "(b)", transform=axis_a.transAxes)

    for axis in (axis_n, axis_a):
        axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout(h_pad=0.18)
    for extension in ("png", "pdf"):
        figure.savefig(OUTPUT / f"total_amplitude_comparison_N99_standard.{extension}",
                       dpi=300, bbox_inches="tight")
    plt.close(figure)
    return OUTPUT / "total_amplitude_comparison_N99_standard.png"


def write_summary(cavity: np.ndarray, disk: np.ndarray) -> Path:
    summary_path = OUTPUT / "total_amplitude_comparison_N99_standard_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("Parameters: a_s=0, n=30, omega_bar=0, c^2=1, N=99\n")
        handle.write("Fourier convention: hhat(alpha) = integral h(s) exp(-i alpha s) ds\n")
        handle.write("Inverse convention: h(s) = (1/(2*pi)) integral hhat(alpha) exp(i alpha s) d alpha\n\n")
        for label, data in (("Cavity", cavity), ("Single disk", disk)):
            handle.write(f"{label} lower neutral R = {data[0, 0]:.12f}\n")
            handle.write(f"{label} |Cr(R_l)| = {data[0, 6]:.12e}\n")
            handle.write(f"{label} N(500) = {data[-1, 4]:.12e}\n")
            handle.write(f"{label} |A(500)| = {data[-1, 7]:.12e}\n\n")
        handle.write(f"Target amplitude ratio disk/cavity = {disk[-1, 7] / cavity[-1, 7]:.12e}\n")
        handle.write("The N factors and the amplitude ratio are unchanged by the Fourier convention.\n")
    return summary_path


def main() -> None:
    zones = read_zones(SOURCE)
    cavity = standardize(zones["rotor-stator cavity"])
    disk = standardize(zones["isolated rotating disk"])
    OUTPUT.mkdir(exist_ok=True)
    comparison_path, disk_path = write_data(cavity, disk)
    figure_path = plot_curves(cavity, disk)
    summary_path = write_summary(cavity, disk)
    print(f"comparison data: {comparison_path}")
    print(f"single-disk data: {disk_path}")
    print(f"figure: {figure_path}")
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
