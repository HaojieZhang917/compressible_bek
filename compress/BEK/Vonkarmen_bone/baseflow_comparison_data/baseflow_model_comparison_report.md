# Base-flow comparison report

## Scope

This report was generated from the scripts in `Bone.ipynb`/`CRD_STA.jl` logic, but uses a SciPy reimplementation to avoid the long Julia precompile path.

* Wall-temperature sweep: `1.000 <= Tw <= 2.000`.
* Compressible settings: `Ro=-1.0`, `Mr=0.3`, `gamma=1.4`, `sigma=0.72`.
* Boussinesq settings: `Pr=0.72`, ideal-gas centrifugal buoyancy coefficient `lambda_c=1`.
* Common comparison coordinate: physical `z in [0, 20.0]`, `N=2001`.

Because the current compressible notebook uses `Ro=-1`, its azimuthal velocity approaches `G_c(infinity)=-1`. The comparison therefore uses `G_compare=-G_i` for the Boussinesq solution while still writing the raw `G_i` columns to the CSV files.

## Main Differences

The Boussinesq model keeps density constant everywhere except the centrifugal buoyancy term. In the present nondimensionalization, the radial-momentum feedback is proportional to `-(T-1)`, so changing wall temperature directly changes the velocity field.

The compressible basic-flow path currently used by `Bone.ipynb` solves the isothermal BEK velocity field first, then constructs temperature, density, axial velocity and physical coordinate through `T_ca` and `z=int T deta`. Therefore the wall temperature mainly changes density, axial scaling and coordinate stretching; the raw tangential/radial similarity velocity is not recomputed from a temperature-coupled momentum system in this script.

## Difference Metrics

| Tw | max abs dF | max abs dG compare | max abs dH | max abs dT | wall density change | rho linearization error |
|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 2.0112e-03 | 2.5416e-03 | 1.4242e-03 | 3.2083e-03 | 0.000 | 0.000 |
| 1.10 | 2.9520e-02 | 6.9566e-02 | 2.0536e-01 | 1.0642e-02 | 0.091 | 0.010 |
| 1.20 | 5.7809e-02 | 1.2546e-01 | 3.2935e-01 | 3.1057e-02 | 0.167 | 0.040 |
| 1.50 | 1.3625e-01 | 2.5728e-01 | 5.6400e-01 | 1.4482e-01 | 0.333 | 0.250 |
| 2.00 | 2.4765e-01 | 4.0906e-01 | 8.0163e-01 | 4.4755e-01 | 0.500 | 1.000 |

Here `wall density change = |1/Tw - 1|`, the ideal-gas relative density change at the wall. `rho linearization error` compares exact ideal-gas `rho_w/rho_inf=1/Tw` with the Boussinesq linearization `1-(Tw-1)`.

## Derivatives And Inflection Points

| Tw | variable | Bouss wall d1 | Comp wall d1 | Bouss max abs d1 | Comp max abs d1 | Bouss inflection z | Comp inflection z |
|---:|:---|---:|---:|---:|---:|:---|:---|
| 1.00 | F | 5.1019e-01 | 2.0647e-01 | 5.1019e-01 | 4.9434e-01 | 1.815385 | - |
| 1.00 | G_compare | -6.1596e-01 | -2.4640e-01 | 6.1596e-01 | 6.1568e-01 | - | - |
| 1.00 | H | -6.5847e-05 | 2.4003e-03 | 3.6152e-01 | 3.6162e-01 | 0.924556 | 0.928937 |
| 1.00 | T | 0.0000e+00 | 3.4046e-03 | 7.1054e-15 | 8.1508e-03 | - | - |
| 1.10 | F | 5.7217e-01 | 1.5989e-01 | 5.7217e-01 | 4.5126e-01 | 1.884926 | - |
| 1.10 | G_compare | -6.5940e-01 | -1.9037e-01 | 6.5940e-01 | 5.6003e-01 | - | - |
| 1.10 | H | -7.2432e-05 | 2.4486e-03 | 4.1690e-01 | 3.5489e-01 | 0.959515 | 0.985564 |
| 1.10 | T | -3.6036e-02 | -7.5249e-03 | 3.6036e-02 | 2.8436e-02 | - | 1.070038 |
| 1.20 | F | 6.3301e-01 | 1.2107e-01 | 6.3301e-01 | 4.1498e-01 | 1.919692 | - |
| 1.20 | G_compare | -6.9442e-01 | -1.4368e-01 | 6.9442e-01 | 5.1367e-01 | - | - |
| 1.20 | H | -7.9025e-05 | 2.4971e-03 | 4.6813e-01 | 3.4916e-01 | 0.972689 | 1.041806 |
| 1.20 | T | -7.6775e-02 | -1.3342e-02 | 7.6775e-02 | 5.3964e-02 | - | - |
| 1.50 | F | 8.0689e-01 | 1.0288e-01 | 8.0689e-01 | 3.3021e-01 | 1.889394 | - |
| 1.50 | G_compare | -7.7487e-01 | -1.2307e-01 | 7.7487e-01 | 4.1214e-01 | - | - |
| 1.50 | H | -9.8825e-05 | 1.4130e-03 | 6.0043e-01 | 3.3623e-01 | 0.961316 | 1.211288 |
| 1.50 | T | -2.1745e-01 | -3.1150e-02 | 2.1745e-01 | 1.1407e-01 | - | - |
| 2.00 | F | 1.0737e+00 | 1.0200e-01 | 1.0737e+00 | 2.5026e-01 | 1.760359 | - |
| 2.00 | G_compare | -8.6956e-01 | -1.2313e-01 | 8.6956e-01 | 3.1074e-01 | - | - |
| 2.00 | H | -1.3185e-04 | -2.8199e-06 | 7.7569e-01 | 3.2303e-01 | 0.907683 | 1.495866;7.078451 |
| 2.00 | T | -4.9243e-01 | -6.4044e-02 | 4.9243e-01 | 1.8125e-01 | - | - |

The inflection locations are detected from finite-difference second derivatives on the common uniform physical grid. Second-derivative extrema are written to `feature_summary.csv`, with their maximum taken from the interior grid (`z>=0.2`) to avoid endpoint artifacts. These values should be used as diagnostics of structural change, not as theorem-level locations.

## Applicability Range

For an ideal gas with `T*=T/T_inf`, the Boussinesq density linearization is `rho/rho_inf ~= 1-(T-1)`. Its formal small parameter is therefore `epsilon=|Tw-1|` near the wall.

A conservative interpretation is:

* `|Tw-1| <= 0.05`: usually a safe Boussinesq range for quantitative stability comparisons.
* `0.05 < |Tw-1| <= 0.10`: useful for trend studies, but basic-flow and derivative differences should be checked case by case.
* `|Tw-1| > 0.10-0.20`: no longer a clean Boussinesq limit for an ideal gas, especially in a rotating flow where the buoyancy term enters the radial momentum balance directly.
* `Tw=2` gives `epsilon=1` and a wall density ratio `rho_w/rho_inf=0.5`; this is outside the traditional Boussinesq asymptotic range.

Using the current data-driven maximum relative velocity difference metric:
* the 5% criterion is satisfied up to approximately `Tw=1.00` in this sweep;
* the 10% criterion is satisfied up to approximately `Tw=1.00` in this sweep.

For the paper logic, this suggests treating Boussinesq as the low-temperature-difference baseline and then using the compressible model once density variation, coordinate stretching, derivative changes or inflection-point migration become non-negligible.

Important caveat: the data-driven velocity-difference range above is a comparison against the current `Bone.ipynb` compressible basic-flow construction. Since that construction does not recompute the radial/tangential velocity field with temperature-coupled density and viscosity in the momentum equations, it should not be quoted as a universal physical failure point of the Boussinesq approximation.

## Output Files

* `incompressible_boussinesq_physical_raw.csv`: raw Boussinesq profiles in physical/similarity `z`.
* `compressible_physical_raw.csv`: raw compressible profiles in physical `z=int T deta`.
* `physical_common_grid_interpolated.csv`: both models interpolated to the same physical grid, with differences.
* `cheb_grid_interpolated.csv`: both models interpolated to the notebook Chebyshev physical grid.
* `derivatives_common_grid.csv`: first and second derivatives on the common physical grid.
* `feature_summary.csv`: wall/edge values, derivative extrema and inflection-point summaries.
* `summary.csv`: one-line difference metrics for each `Tw`.
