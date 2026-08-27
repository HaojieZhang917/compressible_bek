from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_DIRECTORY = Path(__file__).resolve().parent / "c_6_validation_results"
OUTPUT_DIRECTORY = Path(r"E:\Zhj\Latex\twodisk")

CASES = (
    (
        "c6_lower_neutral_A_N.dat",
        r"$R_f=R_l=287.82$",
        "black",
        "-",
    ),
    (
        "c6_maximum_growth_Rg_A_N.dat",
        r"$R_f=R_g=372.99$",
        "#c44e52",
        "--",
    ),
    (
        "c6_global_Cr_peak_A_N.dat",
        r"$R_f=R_p=468.91$",
        "#4c72b0",
        "-.",
    ),
)


def read_tecplot_points(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields:
            continue
        try:
            rows.append([float(value) for value in fields])
        except ValueError:
            continue
    return np.asarray(rows)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
        }
    )

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True)
    ax_n, ax_a = axes

    for filename, label, color, linestyle in CASES:
        data = read_tecplot_points(DATA_DIRECTORY / filename)
        radius, n_factor, amplitude = data[:, 0], data[:, 1], data[:, 2]
        ax_n.plot(radius, n_factor, color=color, linestyle=linestyle, label=label)
        ax_a.plot(radius, amplitude, color=color, linestyle=linestyle, label=label)
        ax_n.plot(radius[0], n_factor[0], marker="o", color=color, markersize=3.5)
        ax_a.plot(radius[0], amplitude[0], marker="o", color=color, markersize=3.5)

    growth_peak = 372.988887508
    absolute_radius = 570.0
    for axis in axes:
        axis.axvline(growth_peak, color="0.55", linewidth=0.9, linestyle=":")
        axis.axvline(absolute_radius, color="0.35", linewidth=0.9, linestyle="--")
        axis.grid(False)
        axis.tick_params(direction="in", top=True, right=True)

    ax_n.text(growth_peak + 2.5, 4.48, r"$R_g$", color="0.35", va="top")
    ax_n.text(absolute_radius - 3.0, 4.48, r"$R_{\mathrm{abs}}$", color="0.25", ha="right", va="top")
    ax_n.set_ylabel(r"$N(R;R_f)$")
    ax_n.set_ylim(0.0, 4.6)
    ax_n.text(0.02, 0.92, "(a)", transform=ax_n.transAxes)
    ax_n.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.10, 0.98))

    ax_a.axhline(1.0, color="0.55", linewidth=0.9, linestyle=":")
    ax_a.text(282.0, 1.08, r"$|A|=1$", color="0.35", va="bottom")
    ax_a.set_yscale("log")
    ax_a.set_ylim(0.1, 15.0)
    ax_a.set_xlim(280.0, 578.0)
    ax_a.set_xlabel(r"$R$")
    ax_a.set_ylabel(r"$|A(R;R_f)|$")
    ax_a.text(0.02, 0.92, "(b)", transform=ax_a.transAxes)

    fig.tight_layout(h_pad=0.15)
    for suffix in ("pdf", "png"):
        fig.savefig(
            OUTPUT_DIRECTORY / f"C6_amplitude_evolution.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
