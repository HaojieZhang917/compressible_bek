# Model definitions and notation

This document records the model choices that must be fixed before new calculations are used in a manuscript claim. It is a working specification, not evidence that the proposed equations or results have already been verified.

## 1. Physical problem and Rossby-number convention

The primary geometry is the Bodewadt--Ekman--von Karman (BEK) family: an infinite plane or disc rotating at dimensional angular velocity `Omega_D`, below fluid in rigid-body rotation at `Omega_F`. The two rotations have the same axis. Their difference is

```text
DeltaOmega = Omega_F - Omega_D.
```

The study adopts Lingwood's BEK differential-rotation Rossby parameter

```text
Ro = DeltaOmega / Omega,
```

with the associated Coriolis parameter

```text
Co = 2 - Ro - Ro^2.
```

The system rotation scale `Omega` and the sign convention must be rederived and checked against the dimensional boundary conditions before implementation. The three canonical limits are

| `Ro` | BEK limit | Physical interpretation |
|---:|---|---|
| `-1` | von Karman | rotating disc beneath otherwise stationary fluid |
| `0` | Ekman | small differential rotation relative to system rotation |
| `1` | Bodewadt | stationary plane beneath rotating fluid |

This `Ro` is not to be identified without qualification with the generic convective-to-Coriolis ratio `U/(2 Omega L)`. Consequently, no `1:Ro:Ro^2` density--acceleration hierarchy may be claimed from dimensional reasoning alone. Every coefficient must be traced through the BEK nondimensionalisation.

Rigid-body rotation of both the wall and far field preserves the radial similarity form. The primary steady base-flow problem should therefore remain an ODE boundary-value problem in the wall-normal similarity coordinate.

## 2. Physical model endpoints

All comparisons must be constructed from one parent dimensional equation, one nondimensionalisation and identical boundary conditions. The first required endpoints are:

1. **Traditional centrifugal Boussinesq model.** This is the preserved baseline model that produces the known von Karman and rotor--stator saddle-nodes.
2. **Canonical acceleration-consistent Boussinesq model.** Density fluctuation is coupled to all acceleration terms required by the frame-consistent Blackburn formulation. The exact steady BEK reduction must be derived rather than inferred from model names or existing code switches.
3. **Matched low-Mach compressible reference.** Selected small-temperature-difference cases are required before either Boussinesq endpoint is called more accurate. Density-law, property, compressibility and closure errors must be reported separately where possible.

The Lopez model is retained as a useful historical and diagnostic comparison, but it is not automatically identical to the canonical Blackburn endpoint because its local-time mass matrix and frame formulation differ.

## 3. Term-resolved diagnostic coefficients

Let `eta_j` multiply the density coupling to the `j`th acceleration contribution after the dimensional and BEK reductions have been checked. Candidate labels include local-time, convective/streamline-curvature, Coriolis and frame-centrifugal contributions.

The coefficients have the following interpretation:

| Value | Meaning |
|---:|---|
| `eta_j=0` | the selected density--acceleration coupling is omitted |
| `eta_j=1` | the coupling is retained at its derived coefficient |
| `0<eta_j<1` | diagnostic homotopy or finite perturbation used for sensitivity and bifurcation tracking |

Intermediate values are not new material or flow parameters and should not be optimised as empirical fit coefficients. Their allowed uses are:

- evaluate derivatives such as `dTw_c/deta_j` or `dalpha/deta_j`;
- validate adjoint predictions with small finite coefficient changes;
- continue between physically derived endpoint models to determine how a fold moves, escapes the validity range or annihilates at a cusp;
- diagnose interactions between acceleration contributions while all other conditions are frozen.

Any term-toggle result must be labelled a diagnostic rather than a fully self-consistent model validation. Independently toggling terms may also break frame covariance or another structural property; this must be checked and disclosed.

## 4. Primary observables and stability language

The primary base-flow observables are the fold locus `Tw_c(Ro)`, branch topology, wall shear, heat transfer, axial entrainment and the smallest Jacobian singular/eigenvalue measures.

At a fold, the intended structural sensitivity is

```text
dTw_c/deta_j = -<psi, partial_eta_j R>/<psi, partial_Tw R>,
```

where `psi` is the left null vector of the steady residual Jacobian. The formula and normalisation must be derived for the discretised constrained problem and verified by finite changes.

Two stability problems must remain distinct:

- **Similarity-subspace temporal stability:** axisymmetric disturbances that preserve the similarity form. This may establish a real zero mode and stability exchange at a saddle-node, but it is not full three-dimensional stability.
- **Local three-dimensional BEK stability:** Type-I/Type-II spatial or temporal modes under the local parallel-flow approximation. This is not finite-radius global stability.

## 5. Central hypothesis and falsification criteria

The working hypothesis is that BEK differential rotation changes the spatial structure and relative importance of density--acceleration couplings, and that small closure residuals can be strongly amplified near a saddle-node, producing large errors in fold and modal observables.

The hypothesis must be revised if one or more of the following occurs:

- traditional and canonical endpoint models have indistinguishable converged fold loci throughout the Boussinesq-valid range;
- apparent folds disappear under domain, resolution or continuation checks;
- term sensitivities fail finite-change validation;
- matched compressible calculations develop the same fold and sensitivity, indicating that the phenomenon is not attributable mainly to the Boussinesq closure;
- model differences are dominated by density-law or property errors rather than density--acceleration coupling.

## 6. Planned parameter record

Before the first production sweep, record in a dedicated file under `docs/`:

- the exact dimensional and nondimensional equations and sign conventions;
- the sampled `Ro`, `Tw`, `Pr`, domain and spectral resolutions;
- the endpoint and diagnostic `eta` paths;
- continuation tolerances and fold/cusp detection criteria;
- all external base-flow or compressible source-data dependencies;
- the output directory under `work/results/`.
