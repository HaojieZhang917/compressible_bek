# Research execution plan

## Gate 1: exact scaling

Derive the dimensional rotating-frame momentum equation, define `Ro`, and separate local, convective, streamline-curvature, Coriolis, frame-centrifugal and imposed-body accelerations. Determine which projected density-force contributions have candidate orders `epsilon`, `epsilon Ro` and `epsilon Ro^2`.

Deliverables:

- notation and dimensional table;
- checked nondimensional derivation;
- candidate density--acceleration numbers;
- reference-frame transformation check.

## Gate 2: term-resolved steady model

Implement a single parent steady equation with continuous coefficients for every density--acceleration contribution. Recover the traditional and acceleration-consistent equations as exact coefficient choices.

Deliverables:

- Julia module in `work/src/`;
- unit tests for limiting models;
- matched base-flow profiles at `Tw=1` and selected small temperature differences.

## Gate 3: fold mechanism

Continue the first fold as a function of each model coefficient and `Ro`. Compute right and left fold null vectors and predict fold displacement using adjoint projection.

Deliverables:

- `Tw_c(Ro, eta_j)` data;
- term-by-term fold sensitivities;
- finite-change validation;
- test of density--acceleration-number collapse.

Stop and revise the scaling if fold shifts do not collapse or if the compressible reference develops the same fold.

## Gate 4: term-resolved stability and modal sensitivity

Add the same coefficients to the Blackburn perturbation operator. Implement the adjoint spatial eigenproblem and compare first-order eigenvalue shifts with direct term-removal calculations for Type-I and Type-II modes.

Deliverables:

- direct/adjoint biorthogonality tests;
- term-resolved `delta alpha` predictions;
- neutral-curve and critical-parameter validation;
- two-mode or pseudospectral treatment near degeneracy if required.

## Gate 5: error-controlled hierarchy

Combine equation-level magnitude with fold/modal amplification. Select the simplest retained-term set for a prescribed observable and tolerance, then validate it against the parent and compressible models.

Deliverables:

- an a posteriori error indicator;
- model-selection map in `(epsilon, Ro)`;
- topology-warning criterion near folds or mode interaction;
- honest separation of density-law, low-Mach, property and truncation errors.

## Gate 6: generality

Apply the same hierarchy to the rotor--stator system and repeat a physical problem in different rotating frames. Only after these checks should the framework be described as unified.

