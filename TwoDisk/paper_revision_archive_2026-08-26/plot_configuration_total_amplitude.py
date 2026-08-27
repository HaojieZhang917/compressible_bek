from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "configuration_total_amplitude_results" / "total_amplitude_comparison_N99.dat"
OUTPUT = ROOT / "configuration_total_amplitude_results"


def read_zones(path: Path) -> dict[str, np.ndarray]:
    zones: dict[str, list[list[float]]] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("ZONE"):
            current = stripped.split('T="', 1)[1].split('"', 1)[0]
            zones[current] = []
            continue
        if current is None or not stripped:
            continue
        try:
            zones[current].append([float(value) for value in stripped.split()])
        except ValueError:
            continue
    return {name: np.asarray(rows) for name, rows in zones.items()}


def main() -> None:
    zones = read_zones(DATA)
    styles = {
        "rotor-stator cavity": ("black", "-", "Rotor--stator"),
        "isolated rotating disk": ("#c44e52", "--", "Isolated disk"),
    }

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
        }
    )
    fig, (ax_n, ax_a) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True)

    for name, data in zones.items():
        color, linestyle, label = styles[name]
        radius, n_factor, amplitude = data[:, 0], data[:, 4], data[:, 7]
        ax_n.plot(radius, n_factor, color=color, linestyle=linestyle, label=label)
        ax_a.plot(radius, amplitude, color=color, linestyle=linestyle, label=label)
        ax_n.plot(radius[0], n_factor[0], "o", color=color, markersize=3.5)
        ax_a.plot(radius[0], amplitude[0], "o", color=color, markersize=3.5)

    for axis in (ax_n, ax_a):
        axis.axvline(500.0, color="0.4", linestyle=":", linewidth=0.9)
        axis.tick_params(direction="in", top=True, right=True)
        axis.grid(False)

    ax_n.set_ylabel(r"$N(R;R_l)$")
    ax_n.set_ylim(0.0, 8.1)
    ax_n.text(0.98, 0.92, "(a)", ha="right", transform=ax_n.transAxes)
    ax_n.legend(frameon=False, loc="upper left")

    ax_a.axhline(1.0, color="0.55", linestyle=":", linewidth=0.9)
    ax_a.set_yscale("log")
    ax_a.set_xlim(280.0, 505.0)
    ax_a.set_ylim(0.04, 300.0)
    ax_a.set_xlabel(r"$R$")
    ax_a.set_ylabel(r"$|A(R;R_l)|$")
    ax_a.text(0.02, 0.92, "(b)", transform=ax_a.transAxes)

    fig.tight_layout(h_pad=0.15)
    for suffix in ("png", "pdf"):
        fig.savefig(
            OUTPUT / f"total_amplitude_comparison_N99.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
