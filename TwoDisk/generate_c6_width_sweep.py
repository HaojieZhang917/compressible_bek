from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "c_6_validation_results" / "c6_lower_neutral_A_N.dat"
OUTPUT_DIRECTORY = ROOT / "c6_width_sweep_results"
OUTPUT = OUTPUT_DIRECTORY / "c6_lower_neutral_Cr_vs_c_N99.dat"

C_MIN = 0.25
C_MAX = 5.00
C_STEP = 0.025
def lower_neutral_data() -> tuple[float, float, float, float]:
    for line in SOURCE.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if len(values) != 8:
            continue
        try:
            radius, _, _, _, alpha_real, alpha_imag, _, cr_abs = map(float, values)
        except ValueError:
            continue
        return radius, alpha_real, alpha_imag, cr_abs
    raise RuntimeError(f"No data row found in {SOURCE}")


def spectral_ratio(c: float, alpha_real: float, alpha_imag: float) -> float:
    alpha_squared_real = alpha_real**2 - alpha_imag**2
    return c * math.exp(-0.5 * (c**2 - 1.0) * alpha_squared_real)


def main() -> None:
    radius, alpha_real, alpha_imag, cr_reference_raw = lower_neutral_data()
    cr_reference = cr_reference_raw
    OUTPUT_DIRECTORY.mkdir(exist_ok=True)
    count = round((C_MAX - C_MIN) / C_STEP)

    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write('TITLE="C6 lower-neutral receptivity coefficient versus Gaussian width"\n')
        handle.write(
            'VARIABLES="c","c_squared","l_s","R_l","alpha_r","alpha_i",'
            '"Cr_abs","Cr_ratio_to_c_1","spectral_ratio"\n'
        )
        handle.write('DATASETAUXDATA N_cheb = "99"\n')
        handle.write('DATASETAUXDATA Re_h = "1000"\n')
        handle.write('DATASETAUXDATA n = "30"\n')
        handle.write('DATASETAUXDATA omega_bar = "0"\n')
        handle.write('DATASETAUXDATA roughness_height = "0.001"\n')
        handle.write('DATASETAUXDATA reference_c = "1.0"\n')
        handle.write('DATASETAUXDATA Fourier_convention = "forward transform without prefactor; inverse prefactor 1/(2*pi)"\n')
        handle.write('DATASETAUXDATA scaling = "Cr(c)/Cr(1) = c exp[-(c^2-1) Re(alpha^2)/2]"\n')
        handle.write(f'ZONE T="lower neutral R={radius:.6f}", I={count + 1}, F=POINT\n')

        for index in range(count + 1):
            c = C_MIN + index * C_STEP
            localization = 1.0 / (2.0 * c**2)
            ratio = spectral_ratio(c, alpha_real, alpha_imag)
            cr_abs = cr_reference * ratio
            handle.write(
                f"{c:.12e} {c**2:.12e} {localization:.12e} {radius:.12e} "
                f"{alpha_real:.12e} {alpha_imag:.12e} {cr_abs:.12e} "
                f"{ratio:.12e} {ratio:.12e}\n"
            )

    alpha_squared_real = alpha_real**2 - alpha_imag**2
    c_optimal = 1.0 / math.sqrt(alpha_squared_real)
    print(f"output: {OUTPUT}")
    print(f"R_lower = {radius:.12f}")
    print(f"alpha_lower = {alpha_real:.12f} {alpha_imag:+.12f}i")
    print(f"Cr_raw(c=1) = {cr_reference_raw:.12e}")
    print(f"Cr(c=1) = {cr_reference:.12e}")
    print(f"c_optimal = {c_optimal:.12f}")


if __name__ == "__main__":
    main()
