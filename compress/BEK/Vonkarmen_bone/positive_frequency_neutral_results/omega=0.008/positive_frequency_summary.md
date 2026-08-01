# Positive-frequency Type-II comparison (omega=0.008)

The selected frequency is Lingwood's modest positive-frequency case. Curves are stopped immediately before the Type-I/Type-II neutral branch point, so all reported minima belong to Type II.

Errors use the fully compressible model as reference: `100*(Blackburn-compressible)/compressible`.

| Tw | Blackburn Rc | Compressible Rc | Rc error | Blackburn beta_c | Compressible beta_c | beta error | curve RMS R error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 296.940474 | 297.975422 | -0.347% | 0.03409902 | 0.03378400 | 0.932% | 0.501% |
| 1.04 | 292.956224 | 295.390189 | -0.824% | 0.03304893 | 0.03220748 | 2.613% | 1.467% |
| 1.08 | 288.346830 | 292.946596 | -1.570% | 0.03204713 | 0.03075064 | 4.216% | 2.071% |
| 1.12 | 339.736998 | 344.702388 | -1.440% | 0.03622437 | 0.03419015 | 5.950% | no common beta interval |
| 1.16 | 330.894829 | 346.341855 | -4.460% | 0.03511117 | 0.03261511 | 7.653% | no common beta interval |
| 1.20 | 319.561837 | 340.546967 | -6.162% | 0.03406515 | 0.03119155 | 9.213% | no common beta interval |

Quadratic fits are used for smooth interior minima. At `Tw>=1.12`, the Type-II minimum is the neutral-branch endpoint, so its last resolved point is reported without extrapolation.

## Operational applicability bands

A temperature passes a band only if the absolute errors in both critical parameters and the curve-wide RMS R error are within it.

- 2% band: 1.00
- 5% band: 1.00, 1.04, 1.08

## Comparison with omega=0 results

| Tw | Blackburn Rc change | Compressible Rc change | omega=0 compressible status | omega>0 status |
|---:|---:|---:|---|---|
| 1.00 | -34.223% | -34.496% | available | available |
| 1.04 | -33.868% | -33.829% | available | available |
| 1.08 | -33.410% | -33.141% | available | available |
| 1.12 | -19.359% | -19.844% | available | available |
| 1.16 | -18.883% | not available | compressible_critical_missing | available |
| 1.20 | -18.675% | not available | compressible_critical_missing | available |
