# Blackburn neutral-curve comparison

The six wall temperatures follow slides 13-14 of `E:/Zhj/Boussinesq.pptx`: `1.00, 1.04, 1.08, 1.12, 1.16, 1.20`.

Critical points are interior minima of `R(beta)`, refined by a three-point quadratic fit. `beta_c >= 0.055` is classified as Type-I; smaller `beta_c` is Type-II.

The manuscript-oriented errors use the fully compressible result as reference: `100*(Blackburn-compressible)/compressible`. The TSV also records the original PPT convention `100*(compressible-Blackburn)/Blackburn`.

## Critical points and errors

| Tw | Mode | Blackburn Rc | Compressible Rc | Rc error | Blackburn beta | Compressible beta | beta error | Status |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| 1.00 | Type-I | 286.049260 | 286.707530 | -0.230% | 0.07746975 | 0.07698150 | 0.634% | available |
| 1.00 | Type-II | 451.436008 | 454.894741 | -0.760% | 0.04638605 | 0.04608354 | 0.656% | available |
| 1.04 | Type-I | 282.246154 | 280.623680 | 0.578% | 0.07498795 | 0.07363198 | 1.842% | available |
| 1.04 | Type-II | 442.989748 | 446.401034 | -0.764% | 0.04527182 | 0.04455882 | 1.600% | available |
| 1.08 | Type-I | 277.988764 | 275.046188 | 1.070% | 0.07255928 | 0.07053982 | 2.863% | available |
| 1.08 | Type-II | 433.018902 | 438.158465 | -1.173% | 0.04416463 | 0.04322116 | 2.183% | available |
| 1.12 | Type-I | 273.290502 | 269.916828 | 1.250% | 0.07018870 | 0.06767163 | 3.720% | available |
| 1.12 | Type-II | 421.296549 | 430.039357 | -2.033% | 0.04314149 | 0.04209713 | 2.481% | available |
| 1.16 | Type-I | 268.172853 | 265.187457 | 1.126% | 0.06787978 | 0.06500274 | 4.426% | available |
| 1.16 | Type-II | 407.9216445644329 |  |  | 0.04219327265004532 |  |  | compressible_critical_missing |
| 1.20 | Type-I | 262.665172 | 260.813082 | 0.710% | 0.06563477 | 0.06252001 | 4.982% | available |
| 1.20 | Type-II | 392.9425052058034 |  |  | 0.041310938541705686 |  |  | compressible_critical_missing |

## Curve validation

| Tw | Model | Segment | Points | beta range | R range | max residual |
|---:|---|---:|---:|---:|---:|---:|
| 1.00 | blackburn | 1 | 110 | 0.040072-0.127272 | 286.053-509.279 | 9.007e-08 |
| 1.00 | compressible | 1 | 108 | 0.040079-0.125679 | 286.709-502.075 | 9.939e-08 |
| 1.04 | blackburn | 1 | 107 | 0.038469-0.123269 | 282.254-505.527 | 9.727e-08 |
| 1.04 | compressible | 1 | 104 | 0.038800-0.121200 | 280.638-508.579 | 9.903e-08 |
| 1.08 | blackburn | 1 | 105 | 0.036861-0.120061 | 277.999-510.995 | 9.815e-08 |
| 1.08 | compressible | 1 | 101 | 0.036352-0.116352 | 275.051-506.369 | 9.973e-08 |
| 1.12 | blackburn | 1 | 102 | 0.035258-0.116058 | 273.299-503.314 | 9.963e-08 |
| 1.12 | compressible | 1 | 97 | 0.034750-0.111550 | 269.919-500.212 | 9.944e-08 |
| 1.16 | blackburn | 1 | 100 | 0.033670-0.112870 | 268.177-505.662 | 9.943e-08 |
| 1.16 | compressible | 1 | 94 | 0.033286-0.107686 | 265.198-504.775 | 9.985e-08 |
| 1.20 | blackburn | 1 | 98 | 0.032104-0.109704 | 262.666-506.547 | 9.974e-08 |
| 1.20 | compressible | 1 | 91 | 0.031942-0.103942 | 260.818-507.999 | 9.903e-08 |

## Isothermal implementation check

At `Tw=1`, Blackburn must reduce exactly to Lopez. On the common beta interval, interpolation gives `max |R_Blackburn-R_Lopez| = 2.893e-04` and `RMS difference = 3.540e-05`.
