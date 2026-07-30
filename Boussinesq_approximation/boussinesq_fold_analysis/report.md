# Boussinesq similarity-fold investigation

## Turning points

| point | Hinf | Tw | scaled Jacobian sigma_min/sigma_max | thermal decay length |
|---:|---:|---:|---:|---:|
| 1 | -0.517317582 | 1.049031885 | 7.502e-15 | 2.6848 |
| 2 | -0.242026525 | 1.040603691 | 1.709e-09 | 5.7386 |

The benchmark-connected branch reaches the first saddle-node, so continuation with Tw as the parameter becomes singular. Continuation with Hinf remains regular and reveals an S-shaped branch with three solutions over part of the temperature range.

## Domain sensitivity

| zmax | point | Hinf | Tw |
|---:|---:|---:|---:|
| 15 | 1 | -0.442813114 | 1.054134347 |
| 15 | 2 | -0.339887319 | 1.053715196 |
| 20 | 1 | -0.517317537 | 1.049031885 |
| 20 | 2 | -0.242034777 | 1.040603680 |
| 25 | 1 | -0.529394885 | 1.048200289 |
| 25 | 2 | -0.201471929 | 1.034311918 |
| 30 | 1 | -0.532044425 | 1.048052969 |
| 30 | 2 | -0.174290096 | 1.031099891 |
| 40 | 1 | -0.532733249 | 1.048022630 |
| 40 | 2 | -0.138583758 | 1.028595662 |

The first fold converges to approximately Tw=1.04802 as zmax increases. The second fold drifts because its small |Hinf| gives a long thermal tail. Results on the upper branch therefore need a substantially larger domain than the benchmark-connected branch.

## Cross-validation at Tw=1.045

| branch | Hinf | shooting profile error | Newton profile error | Newton residual |
|---:|---:|---:|---:|---:|
| 1 | -0.640480291 | 2.065e-10 | 0.000e+00 | 1.015e-10 |
| 2 | -0.362899966 | 7.563e-10 | 7.164e-10 | 1.904e-12 |
| 3 | -0.191295921 | 5.822e-08 | 5.821e-08 | 1.223e-12 |

Adaptive collocation, IVP shooting, and Chebyshev-Newton collocation converge to the same profiles.

## Similarity-preserving temporal stability at Tw=1.045

| branch | Hinf | leading eigenvalue | unstable eigenvalues |
|---:|---:|---:|---:|
| 1 | -0.640480291 | -4.250403e-02+0.000000e+00i | 0 |
| 2 | -0.362899966 | +3.407980e-02+0.000000e+00i | 1 |
| 3 | -0.191295921 | +2.504955e-02-4.019586e-02i | 2 |

The lower branch is stable within the axisymmetric similarity-preserving subspace. The middle branch has one positive real eigenvalue. The finite-domain upper branch has an unstable complex-conjugate pair and is oscillatory unstable in this subspace. Full physical stability also requires general non-axisymmetric and non-similar disturbances.

## Infinite-domain conditions

For a non-isothermal solution, T' behaves asymptotically as exp(Pr*Hinf*eta). Hence Hinf<0 is necessary for exponential thermal decay. Hinf>=0 is inadmissible unless Tw=1. The upper branch approaches Hinf=0 and therefore requires increasingly large domains; its finite-domain eigenvalues must be checked for domain convergence before being interpreted as infinite-domain spectrum.

## Upper-branch domain stability

| zmax | Hinf | Tw | leading eigenvalue | unstable eigenvalues |
|---:|---:|---:|---:|---:|
| 20 | -0.18 | 1.048745224 | +2.592129e-02+5.127082e-02i | 2 |
| 20 | -0.14 | 1.085556150 | +2.255151e-02+1.067368e-01i | 2 |
| 20 | -0.10 | 1.149441124 | +4.512073e-03+1.705858e-01i | 2 |
| 20 | -0.08 | 1.183209164 | -1.123522e-02+2.016079e-01i | 0 |
| 40 | -0.18 | 1.030507273 | +4.656025e-02+0.000000e+00i | 1 |
| 40 | -0.14 | 1.028601080 | +4.462302e-02+0.000000e+00i | 1 |
| 40 | -0.10 | 1.053129457 | +1.512527e-02-3.782771e-02i | 2 |
| 40 | -0.08 | 1.088825859 | -4.399758e-03-7.847383e-02i | 0 |
| 60 | -0.18 | 1.031473944 | +4.438545e-02+0.000000e+00i | 1 |
| 60 | -0.14 | 1.029694394 | +4.266915e-02+0.000000e+00i | 1 |
| 60 | -0.10 | 1.028660669 | +3.930751e-02+0.000000e+00i | 2 |
| 60 | -0.08 | 1.034924164 | +2.346209e-02+0.000000e+00i | 2 |
| 80 | -0.18 | 1.031696385 | +4.414782e-02+0.000000e+00i | 1 |
| 80 | -0.14 | 1.030270107 | +4.226741e-02+0.000000e+00i | 1 |
| 80 | -0.10 | 1.029704651 | +3.909863e-02+0.000000e+00i | 1 |
| 80 | -0.08 | 1.029880313 | +3.689845e-02+0.000000e+00i | 2 |

The apparent stable window at zmax=20 moves as the far boundary is displaced. At zmax=60 and 80 every sampled negative-H upper-branch state remains unstable. For example, at Hinf=-0.14 the leading real eigenvalue converges toward approximately +4.2e-2. The finite-domain second turning point also drifts toward Hinf=0, where a non-isothermal infinite-domain solution loses exponential thermal decay.

## Interpretation

The failure near Tw=1.05 is a fold of the similarity solution, not disappearance of every mathematical solution. The ordinary Tw continuation fails because dHinf/dTw diverges at the fold. The branch remains connected when Hinf is used as continuation coordinate, but cannot be traversed monotonically in Tw. Heating reduces radial outflow F, which weakens axial entrainment |H|; the weaker entrainment thickens the thermal layer and increases the integrated inward buoyancy. This positive feedback produces the saddle-node. Stability and physical selection of the multiple branches require a separate analysis.
