# Full-BEK traditional-Boussinesq availability sweep

## Objective

Construct the base-flow data layer for all subsequent studies by determining,
for the full Lingwood BEK family `-1 <= Ro <= 1`, the wall-temperature interval
on the branch connected to the isothermal solution.  Every limiting point is
to be classified as one of:

1. a finite saddle-node/fold;
2. loss of the exponentially decaying thermal far-field mode as
   `Hinf -> 0-` (infinite thermal tail);
3. a non-uniform coefficient/scaling limit;
4. a passive-temperature endpoint with no traditional centrifugal feedback;
5. unresolved numerical or branch-coordinate failure.

This is a mathematical solvability/applicability map for the traditional
centrifugal-Boussinesq similarity model.  It is not, by itself, an accuracy
range for the Boussinesq approximation relative to a compressible reference.

## Model and regular thermal parameter

The velocity convention is that of `work/src/BEKThermal.jl`:

```text
H' + 2 F = 0,
F'' + Ro*(F^2 + H*F' - (G^2-1)) - Co*(G-1) + B*Theta = 0,
G'' + Ro*(2 F G + H*G') + Co*F = 0,
Theta'' - Pr*H*Theta' = 0,
Co = 2 - Ro - Ro^2.
```

The normalized temperature response satisfies

```text
Theta(0)=1, Theta(infinity)=0,
T = 1 + (Tw-1)*Theta.
```

For the disk-frame traditional centrifugal closure,

```text
Lambda_cf(Ro) = Co^2/(4 Ro),
B = Lambda_cf(Ro)*(Tw-1).
```

`B` is used internally because it remains finite in the collocation residual.
For `Ro != 0,1`, the reported wall temperature is recovered as

```text
Tw = 1 + B/Lambda_cf.
```

At `Ro -> 0`, fixed non-isothermal `Tw` gives unbounded `B`; this is a
non-uniform Ekman limit of the traditional differential-rotation scaling, not
a saddle-node.  The two one-sided limits are therefore computed and reported,
while `Ro=0` is classified separately.  At `Ro=1`, `Co=Lambda_cf=0`; the
traditional disk-frame centrifugal temperature term vanishes and temperature
is passive in the momentum equations.  This endpoint has no model-induced
finite wall-temperature limit, although the physical small-density-variation
requirement remains.

The relation `Lambda_cf=Co^2/(4Ro)` follows from Lingwood's velocity and
viscous scaling: the retained frame centrifugal acceleration has coefficient
`(Omega_D/Omega)^2/Ro = Co^2/(4Ro)` after division by the BEK radial momentum
scale.  This explains rather than removes the non-uniform `Ro=0` limit.

## Primary sweep

- `Pr = 0.72`.
- Pier semi-infinite map with production values
  `N=120, (a,b,c)=(2,0.6,0.5)`.
- Nonlinear residual tolerance: `1e-9` or tighter for production points.
- Negative-Ro refinement retains the established near-von-Karman cusp grid.
- Full-family grid:

```text
[-1, -0.999, -0.998, -0.995, -0.994, -0.9938, -0.99,
 -0.95, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,
 -0.2, -0.1, -0.05, -0.02, -0.01,
  0,
  0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
  0.7, 0.75, 0.8, 0.9, 0.95, 0.98, 0.99, 1]
```

Each regular branch starts at `B=0` and is continued primarily in `Hinf`
toward `Hinf=0-`.  Extrema of `B(Hinf)` are candidate folds and must be
cross-checked with the fixed-`B` Jacobian.  Failure of fixed-`Hinf` Newton
continuation is classified as unresolved unless pseudo-arclength or a bordered
solve identifies the limiting mechanism.

## Classification and wall-temperature limit

- `fold`: the first finite extremum encountered from the isothermal-connected
  branch, with a small fixed-`B` Jacobian singular-value ratio.
- `thermal_tail`: no earlier fold, `Hinf -> 0-`, thermal length
  `ell_T=1/(Pr*(-Hinf)) -> infinity`, and the limiting `B`/`Tw` is estimated
  from several small negative `Hinf` values.
- `coefficient_singularity`: `Ro=0`; the differential-rotation-scaled
  traditional model has no uniform finite-`Tw-1` thermal limit.
- `passive_temperature`: `Ro=1`; `Lambda_cf=0`, so the traditional model
  supplies no momentum-based upper wall-temperature limit.
- `other_fold_or_branch_event`: a verified Jacobian singularity not captured
  by a simple `Tw(Hinf)` extremum.
- `unresolved`: continuation or conditioning stops before one of the above is
  verified.

For heating (`Tw >= 1`), the applicability interval reported for the
isothermal-connected branch is `1 <= Tw <= Tw_limit` when a finite limit is
identified.  Multiple returned branches are stored but do not extend the
isothermal-connected operating interval past its first fold.

## Verification

1. Reproduce both infinite-mapping folds at `Ro=-1`.
2. Repeat representative mechanisms at `N=80,120,160`.
3. Repeat selected cases with at least two mapping parameter sets.
4. Require fold candidates to have a fixed-`B` Jacobian singular-value check.
5. Require tail classifications to remain fold-free under refined sampling and
   to show increasing `ell_T` as `Hinf -> 0-`.
6. Store all new results under
   `work/results/full_bek_traditional_availability/`; do not overwrite
   preserved baselines.

## Source dependencies

- Lingwood (1997), equations (3.5)-(3.18), for the BEK scaling, `Ro`, `Co`,
  similarity equations and boundary conditions.
- Existing converged `Ro=-1` infinite-mapping folds under
  `baselines/saddle_node/` and `work/results/traditional_bek_cusp/` for
  regression only.

