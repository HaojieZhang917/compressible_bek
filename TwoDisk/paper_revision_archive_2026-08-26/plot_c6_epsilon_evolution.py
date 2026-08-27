from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "c_6_epsilon_results" / "c6_lower_neutral_epsilon_N99_dr2.dat"
OUTPUT_PATH = ROOT / "c_6_epsilon_results" / "c6_lower_neutral_epsilon_N99_dr2.png"


def load_data(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments=("T", "V", "D", "Z"))


def main() -> None:
    data = load_data(DATA_PATH)
    radius = data[:, 0]
    epsilon = data[:, 3]

    figure, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    axis.semilogy(radius, epsilon, color="#007f7b", linewidth=2.2)
    for threshold, color in ((0.01, "#6b7280"), (0.05, "#d97706"), (0.1, "#b91c1c")):
        axis.axhline(threshold, color=color, linestyle="--", linewidth=1.1)
        axis.text(radius[-1] + 1.5, threshold, f"{threshold:g}", color=color,
                  va="center", fontsize=9)

    axis.set_xlim(radius[0], radius[-1] + 18)
    axis.set_xlabel(r"$R$")
    axis.set_ylabel(r"$\epsilon(R)$")
    axis.grid(True, which="both", alpha=0.25)
    figure.savefig(OUTPUT_PATH, dpi=220)


if __name__ == "__main__":
    main()
