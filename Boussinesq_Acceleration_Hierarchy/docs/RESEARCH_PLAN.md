# Research execution plan

## Research question

Determine how the BEK differential-rotation Rossby parameter changes the error introduced by Boussinesq density--acceleration closures, and how that error is amplified into changes in base-flow saddle-nodes, branch topology and stability predictions.

The intended causal chain is

```text
BEK differential rotation
    -> acceleration and base-flow structure
    -> traditional versus canonical closure residual
    -> fold/cusp displacement and critical amplification
    -> similarity-subspace and Type-I/Type-II stability errors.
```

`Ro` is the physical BEK family parameter. Term coefficients `eta_j` are diagnostic sensitivity/homotopy coordinates, not physical parameters and not empirical coefficients to be optimised.

## Gate 1: dimensional derivation and BEK reduction

Derive the dimensional rotating-frame momentum equation and separate local, convective/streamline-curvature, Coriolis, frame-centrifugal, Euler and imposed-body accelerations. Reproduce Lingwood's BEK definitions and obtain the similarity equations with

```text
Ro=-1: von Karman,
Ro=0:  Ekman,
Ro=1:  Bodewadt,
Co=2-Ro-Ro^2.
```

Derive the traditional and canonical thermal equations from the same parent formulation. Identify which gradient accelerations may be absorbed into the reference pressure before multiplying by density fluctuation, and verify reference-frame consistency.

Deliverables:

- checked notation, dimensional table and sign convention;
- derivation of the steady BEK ODEs and boundary conditions;
- term-by-term map from dimensional acceleration to the reduced equations;
- exact endpoint-model definitions and the permitted diagnostic `eta_j` paths;
- derivation and reference-frame checks recorded in `docs/` before numerical claims.

Stop condition: do not implement a parameter sweep until the three BEK limits and the traditional/canonical endpoint map are unambiguous.

## Gate 2: verified self-consistent steady models

Implement reusable Julia code in `work/src/` for one parent steady residual with physically derived endpoint configurations. Use `eta_j` only for derivative checks and diagnostic continuation.

Verification sequence:

1. recover the isothermal BEK profiles across representative `Ro` values;
2. recover the von Karman, Ekman and Bodewadt limits;
3. reproduce the preserved traditional von Karman fold at `Ro=-1`;
4. confirm boundary residual, spectral-order and far-field convergence;
5. compare endpoint profiles at `Tw=1` and selected small temperature differences.

Deliverables:

- Julia module and unit/regression tests;
- an isothermal BEK benchmark table;
- endpoint base-flow profiles and residual audits;
- a documented production-sweep configuration under `docs/`;
- all new outputs under a new directory in `work/results/`.

## Gate 3: BEK bifurcation map

Use pseudo-arclength continuation to compute self-consistent traditional and canonical base-flow branches as `Ro` and `Tw` vary. Construct the fold loci

```text
Tw_c^traditional(Ro),
Tw_c^canonical(Ro),
```

within the numerically converged and Boussinesq-relevant range. Search for fold-pair creation/annihilation and cusp points rather than assuming that the von Karman fold persists throughout the family.

Deliverables:

- converged branch and fold data in the `(Ro,Tw)` plane;
- fold/cusp detection and non-degeneracy checks;
- profiles and balance diagnostics on each relevant branch;
- explicit separation of self-consistent endpoint comparisons from frozen or term-toggle diagnostics.

Falsification check: if endpoint fold loci are indistinguishable throughout the valid range, revise the closure-induced-bifurcation hypothesis rather than extending to expensive stability sweeps.

## Gate 4: fold mechanism and adjoint structural sensitivity

At representative folds, compute right and left null vectors and predict the displacement caused by each density--acceleration contribution:

```text
dTw_c/deta_j = -<psi, partial_eta_j R>/<psi, partial_Tw R>.
```

Validate each derivative against small finite changes in `eta_j`. Use continuous `eta_j` paths only where needed to determine whether a fold shifts, exits the model-valid range or annihilates with another fold at a cusp.

Deliverables:

- term-by-term fold sensitivities versus `Ro`;
- finite-change linearity ranges and error tables;
- physical-space maps of the forcing terms and left/right null modes;
- an identified mechanism for any change in fold topology.

Stop condition: do not describe a single acceleration term as causal unless its adjoint prediction, finite-change response and endpoint comparison agree.

## Gate 5: matched low-Mach reference

For selected small-temperature-difference cases, compare the two Boussinesq endpoints with a self-consistent low-Mach compressible reference under matched geometry, rotation and property assumptions.

Deliverables:

- matched reference cases near and away from the predicted fold region;
- separation of density--acceleration closure, density-law, property and compressibility errors;
- evidence for or against describing the traditional fold as closure-induced.

Do not infer accuracy from the canonical derivation alone, and do not treat formal large-`Tw` continuation as quantitative Boussinesq validation.

## Gate 6: similarity-subspace temporal stability

Compute the temporal spectrum along representative branches and folds. Verify whether a single real eigenvalue crosses zero, whether stability exchanges at the folds, and how the recovery time changes with `Ro` and closure.

Deliverables:

- spectral-order convergence and constraint checks;
- leading eigenvalues and modes along the selected branches;
- a critical-slowing-down map;
- language explicitly restricted to axisymmetric similarity-subspace stability.

## Gate 7: local three-dimensional modal sensitivity

Only after the base-flow mechanism is established, add the same physically derived closure endpoints and diagnostic derivatives to the Blackburn/BEK perturbation operator. Implement the adjoint spatial eigenproblem and compare first-order eigenvalue shifts with direct endpoint or small-`eta_j` calculations for Type-I and Type-II modes.

Deliverables:

- direct/adjoint biorthogonality and original-pencil residual tests;
- term-resolved `delta alpha` predictions;
- selected neutral curves and critical-parameter validation across `Ro`;
- a two-mode or pseudospectral treatment near modal degeneracy if required;
- language restricted to local parallel-flow stability.

## Gate 8: observable error indicator

Test whether equation-level closure residuals combined with fold or modal amplification predict errors in observables such as `Tw_c`, `alpha`, `R_c` and `beta_c`.

Deliverables:

- a verified first-order observable-error estimate;
- a topology-warning criterion near folds, cusps or mode interaction;
- a map of where simplified closure results are insensitive or critically amplified;
- an honest statement of whether the indicator is empirical, asymptotic or a proved bound.

Do not present intermediate `eta_j` values as calibrated physical models. The practical comparison remains between physically derived endpoint closures.

## Gate 9: generality and manuscript scope

Apply the established mechanism to the rotor--stator system only after the BEK result is complete. Check whether confinement changes the fold mechanism without confusing similarity-subspace stability with full three-dimensional stability.

For a focused first manuscript, prioritise:

1. BEK derivation and endpoint models;
2. traditional/canonical fold and cusp map;
3. adjoint fold mechanism;
4. selected low-Mach validation.

The full Type-I/Type-II survey and rotor--stator generalisation may become separate studies if including them obscures the main closure-induced critical-amplification result.
