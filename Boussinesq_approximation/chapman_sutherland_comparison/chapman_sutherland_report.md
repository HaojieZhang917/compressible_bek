# Chapman-DH vs Sutherland-BDF2 Comparison

Case: `Tw=1.5`, Chapman local `Mr=0.3`, Sutherland final station `r=20`, `Mr=0.3`.

## Mathematical Difference

Chapman-DH uses `mu=T`, so with ideal gas `rho=1/T` the product `rho*mu=1`. The Dorodnitsyn-Howarth transformation removes the temperature dependence from the radial and azimuthal basic-flow equations.

Sutherland-BDF2 uses `mu=T^(3/2)*(1+S/Tinf)/(T+S/Tinf)`, so `rho*mu=mu/T` is not constant. Temperature and viscosity therefore remain in the momentum diffusion terms, and the basic velocity amplitudes can vary with radius and wall temperature.

## Peak Metrics

| model | F/U max | z peak | G/V min | H/W min | Tmax | wall dT |
|:---|---:|---:|---:|---:|---:|---:|
| Chapman-DH | 1.807657e-01 | 1.3200 | -9.999999e-01 | -8.901750e-01 | 1.500000e+00 | -5.201553e-02 |
| Sutherland-BDF2 | 1.749419e-01 | 1.2000 | -1.000000e+00 | -8.426818e-01 | 1.500000e+00 | -1.165266e-01 |

## Maximum Profile Differences

| quantity | max abs difference |
|:---|---:|
| max_abs_dF | 1.106244e-02 |
| max_abs_dG | 2.897570e-02 |
| max_abs_dH | 6.214984e-02 |
| max_abs_dT | 1.304627e-02 |

## Generated Files

* `profiles_chapman_vs_sutherland.png/pdf`
* `differences_chapman_vs_sutherland.png/pdf`
* `marching_summary_vs_chapman.png/pdf`
* `chapman_sutherland_common_grid.csv`
* `chapman_sutherland_metrics.csv`
