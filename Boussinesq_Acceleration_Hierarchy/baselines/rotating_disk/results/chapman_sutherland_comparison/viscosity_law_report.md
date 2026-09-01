# Chapman vs Sutherland Viscosity Laws

Temperature range: `T/T_inf = 1 ... 2`.
Sutherland constant: `S=114 K`, `T_inf=273 K`, so `S/T_inf=0.417582`.

Definitions:

* Chapman: `mu = T`.
* Sutherland: `mu = T^(3/2) * (1 + S/T_inf) / (T + S/T_inf)`.
* Ideal gas density: `rho = 1/T`.

The key momentum coefficient is `rho*mu`.  Chapman gives `rho*mu=1`, while Sutherland gives a temperature-dependent coefficient.

## Sample Values

| T | mu Chapman | mu Sutherland | mu diff % | rho mu Chapman | rho mu Sutherland | rho mu diff % |
|:---:|---:|---:|---:|---:|---:|---:|
| 1.000 | 1.000000e+00 | 1.000000e+00 | 0.000 | 1.000000e+00 | 1.000000e+00 | 0.000 |
| 1.100 | 1.100000e+00 | 1.077668e+00 | -2.030 | 1.000000e+00 | 9.796983e-01 | -2.030 |
| 1.500 | 1.500000e+00 | 1.358098e+00 | -9.460 | 1.000000e+00 | 9.053988e-01 | -9.460 |
| 2.000 | 2.000000e+00 | 1.658487e+00 | -17.076 | 1.000000e+00 | 8.292434e-01 | -17.076 |

Generated files:

* `viscosity_laws_chapman_sutherland.png/pdf`
* `viscosity_laws_chapman_sutherland.csv`
* `viscosity_law_samples.csv`
