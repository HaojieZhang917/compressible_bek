# Current research status

Audit date: 2026-09-03.

The append-only chronological evidence record is maintained in
`docs/RESEARCH_PROGRESS_LOG.md`.  New numerical conclusions must be logged
there before downstream use.

## Full-BEK traditional availability baseline (2026-09-02)

- The full heating branch of the traditional centrifugal-Boussinesq BEK
  model has been classified over `-1 <= Ro <= 1` using the regular forcing
  coordinate `B=Lambda_cf*(Tw-1)` and normalized temperature response.
- The isothermal-connected branch ends at a finite fold only in the narrow
  von Karman interval through approximately `Ro=-0.9938`; the fold pair is
  absent by `Ro=-0.9936`.
- For `-0.9936 <= Ro < 0`, the limiting mechanism is `Hinf -> 0-` and an
  infinite thermal tail.  The wall-temperature interval collapses to the
  isothermal point as `Ro -> 0-`.
- `Ro=0` is a non-uniform coefficient/scaling limit of the traditional
  differential-rotation formulation and is not assigned a non-isothermal
  basic flow.
- For `0 < Ro < 1`, the heating branch is verified without a fold through
  `Tw=2`; this is formal solvability, not a physical high-temperature
  Boussinesq accuracy claim.  At `Ro=1` the traditional temperature term is
  passive because `Lambda_cf=0`.
- The authoritative index is
  `work/results/full_bek_traditional_availability/production_v2_N120_a2.0_b0.6_c0.5/final_applicability.csv`.
  Detailed interpretation and downstream-use rules are in
  `docs/FULL_BEK_TRADITIONAL_RESULTS.md`.

Status labels:

- **Completed baseline**: code, numerical output and at least one independent convergence or integrity check are present.
- **Partial**: useful evidence or infrastructure exists, but the proposed density--acceleration hierarchy claim is not yet established.
- **Not started**: no directly usable implementation or result was found.

## Stage-by-stage audit

| Research stage | Status | Existing evidence | Missing work |
|---|---|---|---|
| 1. Verify the traditional-model saddle-node | **Completed baseline** | Rational-Chebyshev semi-infinite mapping confirms the first von Karman fold at `Tw=1.048021731`; a second returned-branch fold is also confirmed. Finite-domain checks, Jacobian singularity, zero modes, square-root scaling and temporal spectra are available. | Recompute only if the governing equations or nondimensionalisation are changed. |
| 2. Determine whether the fold is closure-induced | **Partial, strong evidence** | The infinite-mapped acceleration-consistent `Ro=-1` branch now continues from `Tw=1` to `1.99` without a fixed-`Tw` fold; pseudo-arclength and Jacobian checks distinguish its smooth `Hinf` minimum from the traditional saddle-node. A Chapman-type compressible comparison and one Sutherland marching case also exist. | Perform a matched, self-consistent compressible/low-Mach comparison near the fold under the exact new parameterisation. Do not infer quantitative high-temperature accuracy from the formal continuation. |
| 3. Derive the thermal BEK formulation | **Partial** | The dimensional disk-frame coefficient and the steady acceleration-consistent nonzero-`Ro` residual have been derived and reduce correctly at `Ro=-1`. | Complete the manuscript-grade arbitrary-frame derivation, the non-uniform `Ro=0` limit and the time-dependent canonical mass matrix. |
| 4. Build self-consistent endpoint models | **Partial** | `work/src/BEKConsistent.jl` provides a pure-Julia steady acceleration-consistent residual and continuation on the same infinite map as the traditional endpoint. The `Ro=-1` branch has been verified. | Unify both endpoints behind one tested interface and validate the consistent residual away from `Ro=-1` before a full-family sweep. |
| 5. Compute fold structural sensitivity | **Partial** | Fold Jacobians, null modes, non-degeneracy coefficients and local normal-form evidence exist. | Add left-null-vector projections for each density--acceleration term and predict `dTw_c/deta_j`; validate against finite coefficient changes. |
| 6. Map BEK bifurcations over differential rotation | **Partial, two-model topology and fold--tail continuation complete** | The common eight-`Ro` topology table is complete. The consistent `Ro=-0.5` fold has been continued with a bordered full-state method and compared with finite-`Hinf` tail proxies. The resolved fold arm approaches the delocalisation layer; apparent finite cusps retreat toward `Hinf=0` with increasing `N`. | A tail-adapted singular-limit corrector is needed only if a manuscript-grade connection coordinate is required. Do not quote the fixed-`N` cusp-like reversals as physical cusps. |
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
- The new `N=120` rational-infinite-map computation continues the
  acceleration-consistent von Karman branch to `Tw=1.99`.  It contains a
  smooth observable minimum at approximately
  `(Tw,Hinf)=(1.64705,-0.912081)`, but no fixed-`Tw` saddle-node: the
  Jacobian remains nonsingular and pseudo-arclength `Tw` remains strictly
  increasing.  This single-parameter slice cannot establish a cusp.
- Representative consistent-model slices show a broader topology
  reorganisation rather than universal fold removal.  The `Ro=-0.5` slice
  has a converged fold at `Tw=1.721704197`, whereas `Ro=-0.25,-0.1` terminate
  through infinite thermal tails and the sampled positive-`Ro` slices remain
  smooth to `Tw=1.99`.
- The unified eight-point comparison is
  `traditional = [fold,tail,tail,tail,tail,smooth,smooth,passive]` versus
  `consistent = [smooth,smooth,fold,tail,tail,smooth,smooth,smooth]` in the
  stated `Ro` order.  Thus the closure reorganises, rather than uniformly
  removes, branch boundaries.
- A bordered two-parameter continuation from the consistent `Ro=-0.5` fold
  gives a converged finite-core arm through `Hinf=-0.02`.  Fixed-resolution
  cusp-like turns do not converge: their distance from `Hinf=0` decreases
  strongly between `N=80,120,160`.  The supported interpretation is approach
  to the thermal-tail delocalisation boundary, not a verified finite cusp.
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

1. a generic `1:Ro:Ro^2` density--acceleration scaling applies to Lingwood's BEK differential-rotation `Ro`;
2. the traditional von Karman fold persists across the BEK family or terminates at a cusp;
3. the first fold location is controlled by one identified omitted acceleration term;
4. an adjoint modal-sensitivity calculation predicts neutral-curve errors;
5. an error-controlled hierarchy works across both single-disc and rotor--stator geometries;
6. any similarity-subspace stability result implies full three-dimensional stability.

## Agreed research framing as of 2026-09-01

- `Ro` denotes the BEK differential-rotation family parameter, not an independently introduced generic convective Rossby number.
- The primary physical comparison is between self-consistent traditional and canonical Boussinesq endpoint models.
- `eta_j` is a diagnostic coordinate for derivatives and homotopy tracking. Intermediate values are not proposed physical closures and will not be fitted for accuracy.
- The main hypothesis is critical amplification: a small density--acceleration closure residual may create a large fold or modal error when the steady or stability operator is nearly singular.
- The first production objective is a converged traditional/canonical fold map in `(Ro,Tw)`, followed by adjoint mechanism identification. A full neutral-curve sweep is deferred until this base-flow gate succeeds.

## Gate 1 implementation update (2026-09-01)

- Lingwood (1997) equations (3.11)--(3.14) have been transcribed into the
  standalone Julia module `work/src/BEKIsothermal.jl` using `U=-F`, `W=-H`,
  `V=G` so that `Ro=-1`, `Co=2` exactly recovers the existing von Karman
  isothermal residual.
- `work/scripts/gate1_bek_isothermal.jl` defines the first `Ro` scan and a
  von Karman profile regression against `baselines/saddle_node/src/VonKarman.jl`.
- The implementation is deliberately isothermal: no temperature, `Tw`,
  density law, or Boussinesq acceleration coefficient has been introduced.
- The Windows Julia 1.12.6 runtime is available and the isothermal scan now
  runs successfully. All sampled `Ro` points converged with residuals below
  `2.2e-10`; the `Ro=-1` residual is `4.1e-11`. This validates the mapped
  isothermal operator, but not a fold locus. Degree/mapping convergence is
  still required before any `(Ro,Tw)` saddle-node continuation.

- The traditional thermal BEK extension is now implemented in
  `work/src/BEKThermal.jl`. At `Ro=-1` with the adopted Pier map and `N=120`,
  `H_infinity` continuation reproduces the established first fold at
  `Tw=1.048021731462` (baseline `1.0480217312207`). This validates the
  thermal solver at the von Karman endpoint; the cross-`Ro` fold map remains
  to be computed.
- Cross-validation of the second fold shows that the old `zmax=20` point
  (`Hinf≈-0.242`, `Tw≈1.0406`) is finite-domain sensitive. The present
  infinite-endpoint calculation gives `Hinf≈-0.11332`, `Tw≈1.0301455`, in
  agreement with the trusted prior infinite-mapping value
  `Hinf≈-0.1133057`, `Tw≈1.03014466`.
- Targeted cusp bracketing shows continuous fold-pair approach. At
  `Ro=-0.9940`, the folds are `(-0.46827,1.0556843)` and
  `(-0.41816,1.0556248)` in `(Hinf,Tw)`; at `Ro=-0.9938` they are
  `(-0.44874,1.05603292)` and `(-0.44194,1.05603270)`. No two extrema were
  detected at `Ro=-0.9936` through `-0.9930` in the current scan. This
  brackets a cusp near `Ro≈-0.9938`; its precise coordinates still require a
  bordered fold solve.
- A first traditional-model `Ro` scan is available under
  `work/results/traditional_bek_folds/`. It uses `N=120`, Pier parameters
  `(2,0.6,0.5)`, `Delta Hinf=0.01` and detects the two `Ro=-1` folds plus one
  turning point at `Ro=0.75`. This is a preliminary trajectory only: the
  `Ro=1` branch stopped at its first continuation step and the coarse scan
  must not be interpreted as a complete fold-locus map.

## Traditional BEK solvability study update (2026-09-01)

- The corrected traditional coefficient `Lambda_cf=Co^2/(4Ro)` was used to
  map the isothermal-connected thermal branch over negative `Ro`.
- The trusted `Ro=-1` folds persist at approximately
  `(-0.53276,1.04802)` and `(-0.11331,1.03014)` in `(Hinf,Tw)`. The returned
  branch then approaches a long-tail state with `Hinf -> 0-`.
- At `Ro=-0.999,-0.998,-0.994`, the finite fold pair moves together; at
  `Ro=-0.994` the resolved pair is approximately
  `(-0.4683,1.05502)` and `(-0.4181,1.05496)`. Fixed-`Hinf` tracking becomes
  ill-conditioned near the merger, so the cusp coordinate is not yet claimed.
- For `Ro=-0.95,-0.90,-0.75,-0.50`, no finite fold was resolved on the
  connected branch. At `Hinf=-0.02`, the wall temperatures are approximately
  `1.87,1.90,1.86,1.674`, respectively. For `Ro=-0.50`, `N=80,120,160`
  agree to about `3e-6` at this diagnostic point.
- The thermal far-field mode has decay length
  `ell_T=1/(Pr*(-Hinf))`; it diverges as `Hinf -> 0-`. Thus finite folds and
  long-tail termination are locally different mechanisms, but the results
  support a connected evolution of the steady-solution boundary. A bordered
  cusp solve and explicit `Hinf=0` asymptotics remain required for a proof.
- A full-state weighted pseudo-arclength continuation on the fixed
  `Hinf=-0.02` section follows the high-temperature sheet from `Ro=-0.5` to
  `Ro=-0.99996896`, with `Tw=1.04153646` and residual `3.2e-9` at the last
  resolved point. This supports, but does not prove, connection to the
  von-Karman fold-controlled region; the final approach is limited by
  endpoint conditioning.
