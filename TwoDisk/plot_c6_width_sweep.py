from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "c6_width_sweep_results" / "c6_lower_neutral_Cr_vs_c_N99.dat"
OUTPUT = Path("/mnt/e/Zhj/Latex/twodisk/c_vary.eps")


def load_curve() -> tuple[list[float], list[float]]:
    widths = []
    coefficients = []
    for line in SOURCE.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if len(values) != 9:
            continue
        try:
            row = [float(value) for value in values]
        except ValueError:
            continue
        widths.append(row[0])
        coefficients.append(row[6])
    return widths, coefficients


def main() -> None:
    widths, coefficients = load_curve()
    peak_index = max(range(len(coefficients)), key=coefficients.__getitem__)

    plt.rcParams.update({"font.family": "serif", "font.size": 10})
    figure, axis = plt.subplots(figsize=(4.4, 3.2))
    axis.plot(widths, coefficients, color="#1f4e79", linewidth=1.8)
    axis.scatter(
        [widths[peak_index]],
        [coefficients[peak_index]],
        color="#b33a3a",
        edgecolors="none",
        s=24,
        zorder=3,
    )
    axis.set_xlabel(r"Gaussian width $c$")
    axis.set_ylabel(r"Local receptivity coefficient $|C_r|$")
    axis.set_xlim(min(widths), max(widths))
    axis.set_ylim(0.0, 1.08 * max(coefficients))
    axis.grid(True, color="#d0d0d0", linewidth=0.5, alpha=0.7)
    figure.tight_layout()
    figure.savefig(OUTPUT, format="eps", dpi=600)
    plt.close(figure)


if __name__ == "__main__":
    main()
