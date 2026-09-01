# Current research status

Audit date: 2026-08-31.

Status labels:

- **Completed baseline**: code, numerical output and at least one independent convergence or integrity check are present.
- **Partial**: useful evidence or infrastructure exists, but the proposed density--acceleration hierarchy claim is not yet established.
- **Not started**: no directly usable implementation or result was found.

## Stage-by-stage audit

| Research stage | Status | Existing evidence | Missing work |
|---|---|---|---|
| 1. Verify the traditional-model saddle-node | **Completed baseline** | Rational-Chebyshev semi-infinite mapping confirms the first von Karman fold at `Tw=1.048021731`; a second returned-branch fold is also confirmed. Finite-domain checks, Jacobian singularity, zero modes, square-root scaling and temporal spectra are available. | Recompute only if the governing equations or nondimensionalisation are changed. |
| 2. Determine whether the fold is closure-induced | **Partial, strong evidence** | Traditional and Lopez/acceleration-consistent steady models have been compared. The Lopez branch passes the first traditional fold. A Chapman-type compressible comparison and one Sutherland marching case exist. | Perform a matched, self-consistent compressible/low-Mach comparison near the fold under the exact new parameterisation. Do not infer the conclusion only from formal high-temperature Lopez continuation. |
| 3. Derive density--acceleration scaling | **Not started** | Existing documentation derives the Lopez stability equations and explains rotating-frame terms qualitatively. | Start from dimensional equations; fix the `Ro` definition; distinguish streamline curvature from frame centrifugal acceleration; apply pressure-gradient/Leray projection; derive candidate `B_j`. |
| 4. Build a nested model hierarchy | **Not started** | Traditional, Lopez, Blackburn and compressible models exist as separate implementations. Property switches provide diagnostic comparisons. | Replace discrete model labels by continuous, term-resolved coefficients within one parent formulation and derive formal truncation orders. |
| 5. Compute fold structural sensitivity | **Partial** | Fold Jacobians, null modes, non-degeneracy coefficients and local normal-form evidence exist. | Add left-null-vector projections for each density--acceleration term and predict `dTw_c/deta_j`; validate against finite coefficient changes. |
| 6. Test scaling over parameter space | **Not started** | Temperature sweeps and a rotor--stator Reynolds-number scan exist. | Vary `Ro` with self-consistent base flows, calculate `Tw_c(Ro)`, and test collapse against the proposed `B_j`. Current rotating-disc comparisons mainly use `Ro=-1`. |
| 7. Compute modal sensitivity | **Partial infrastructure** | Blackburn spatial-stability operator, continuation runner, neutral curves, benchmarks, grid studies and tests exist. | Implement the adjoint Blackburn eigenproblem and term-resolved eigenvalue sensitivity. Existing direct neutral curves are not an adjoint sensitivity analysis. |
| 8. Construct an error estimator | **Not started** | Empirical critical-point, curve-shape, topology and density-linearisation errors have been tabulated. | Combine equation-level magnitude with fold/modal amplification, calibrate an a posteriori indicator and verify predicted versus observed errors. Do not call it a rigorous bound unless one is proved. |
| 9. Demonstrate generality | **Partial** | Traditional-model folds are available for both the semi-infinite disc and rotor--stator similarity systems. | Apply the same term-resolved hierarchy and error indicator to both geometries and check reference-frame consistency. |

## Confirmed baseline findings

### Traditional von Karman model

- First fold: `Hinf approximately -0.53276166`, `Tw approximately 1.048021731`.
- Second confirmed returned-branch fold: `Hinf approximately -0.11331`, `Tw approximately 1.03014466`.
- The isothermal-connected branch is stable within the axisymmetric similarity subspace before the first fold; the returned branch is unstable.
- A third near-zero-`Hinf` turning structure remains exploratory because of the singularly long thermal tail.

### Rotor--stator traditional model

- At `Re_h=1000`, folds occur near `Tw=1.155517549` and `1.167676484`.
- Within the axisymmetric similarity subspace, the two outer branches are stable and the middle branch has one real unstable direction.
- The Reynolds scan found two folds at `Re_h=800` and none in the scanned interval at `Re_h=400`; this is not yet a continuous fold boundary in Reynolds number.

### Acceleration-consistent and compressible comparisons

- The Lopez steady base flow continues through the traditional first fold, supporting the hypothesis that the early termination is linked to selective fixed-centrifugal density loading.
- The formal Lopez continuation to large `Tw` is a model diagnostic, not a validated high-temperature prediction.
- Blackburn versus fully compressible zero-frequency curves at `Mr=0.3` show selective accuracy: Type-I critical Reynolds numbers remain close over a wider range than critical wavenumbers, full-curve geometry or Type-II topology.
- The existing analysis reports Type-II agreement through the sampled point `Tw=1.12` and a topology mismatch beginning at the sampled point `Tw=1.16`; these are discrete observations, not an exact transition temperature.
- A refined Type-II study shows that `Delta beta=8e-4` misses the shallow fold at `Tw=1.14`, whereas `Delta beta=2e-4` resolves it consistently for high spectral orders. No fold was found at `Tw=1.16` for `N=69,79,89`.

## Verified numerical infrastructure

- `baselines/saddle_node/check_project.py` checks the reorganised project and mapped fold data.
- `baselines/rotating_disk/test_blackburn_stability.jl` contains 29 operator-linearisation tests.
- `baselines/rotating_disk/test_blackburn_runner_smoke.jl` exercises the runner on a small eigenvalue calculation.
- Neutral-curve continuation includes eigenvector-overlap mode tracking, step refinement and residual checks.
- The current Type-II result set includes spectral-order and beta-step convergence evidence.

## Claims that are not yet supported

The current files do **not** yet establish that:

1. the user's `1:Ro:Ro^2` scaling holds under a fully specified Rossby-number convention;
2. a single density--acceleration number predicts model error across `Ro`;
3. the first fold location is controlled by one identified omitted acceleration term;
4. an adjoint modal-sensitivity calculation predicts neutral-curve errors;
5. an error-controlled hierarchy works across both single-disc and rotor--stator geometries;
6. any similarity-subspace stability result implies full three-dimensional stability.

