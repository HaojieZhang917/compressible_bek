# Consistent-model fold and thermal-tail continuation

## Objective

Starting from the verified acceleration-consistent fold at

```text
Ro=-0.5, Tw=1.721704197, Hinf=-0.29758392,
```

continue the finite-fold locus in the two-parameter `(Ro,Tw)` plane in both
directions and independently continue the thermal-tail boundary
`Hinf=0`.  The primary question is whether the finite-fold locus ends through
a cusp-type degeneracy or approaches the delocalisation boundary.

## Model and discretisation

- Fully self-consistent steady acceleration-density model in
  `work/src/BEKConsistent.jl`; this is not a frozen-base-flow or term-toggle
  diagnostic.
- `Pr=0.72`, `gamma=1`.
- Pier semi-infinite map with
  `N=120`, `(a,b,c)=(2,0.6,0.5)` and Newton tolerance `1e-9`.
- No temperature-dependent viscosity, low-Mach correction, stability
  calculation or additional property parameter is introduced.

## Fold definition and continuation

At a fold, solve the augmented system

```text
R(u,Tw;Ro) = 0,
J(u,Tw;Ro) v = 0,
v'v = 1,
```

where `u` contains all four collocation fields and `J=R_u` is the fixed-`Tw`
Jacobian.  A full-state pseudo-arclength condition closes the system when
`Ro` is released.  The fold corrector uses the analytic directional Hessian
needed for the `d(Jv)/du` block and the analytic residual derivative with
respect to `Ro`.

A loss of a fixed-`Ro` fold by itself will not be labelled a cusp.  A cusp
claim additionally requires loss of fold non-degeneracy or a resolved meeting
of the two fold arms in the parameter plane.

## Thermal-tail definition and continuation

The exact `Hinf=0` state is non-localised on the semi-infinite domain.  It is
therefore approached using several fixed negative levels

```text
Hinf = -0.02, -0.01, -0.005, -0.0025
```

continued in `Ro`, followed by an `Hinf -> 0-` extrapolation of `Tw` at common
`Ro`.  Agreement across these levels and divergence of the thermal length are
required.  The extrapolated boundary is a delocalisation limit, not a regular
finite-domain solution curve.

## Planned outputs

- Shared continuation code: `work/src/BEKConsistent.jl`.
- Executable study: `work/scripts/continue_consistent_fold_tail.jl`.
- New results: `work/results/consistent_fold_tail_continuation/`.
- A parameter-plane figure containing the fold locus and tail boundary.
- Endpoint diagnostics stating whether the evidence supports a cusp,
  fold--tail connection, or remains unresolved.

All completed numerical conclusions must be appended to
`docs/RESEARCH_PROGRESS_LOG.md`.

## Results

### Fold seed and resolved primary arm

The bordered fold corrector gives at the production discretisation

```text
Ro   = -0.5
Tw   = 1.721704195473
Hinf = -0.297584911495
```

The augmented residual is below `1e-9`, and the row-scaled fixed-`Tw`
Jacobian singular value is at roundoff.  The analytic directional Hessian
and `Ro` derivative used by the corrector were independently checked against
finite differences with relative errors `7.5e-8` and `7.5e-10`.

The primary fold arm is converged while it remains separated from the
thermal-tail singular limit.  Representative intersections are

| `Hinf` | `Ro`, `N=120` | `Ro`, `N=160` | `Tw`, `N=120` | `Tw`, `N=160` |
|---:|---:|---:|---:|---:|
| -0.10 | -0.36821936 | -0.36822351 | 1.55281242 | 1.55281799 |
| -0.05 | -0.34300303 | -0.34300432 | 1.52191682 | 1.52191850 |
| -0.02 | -0.32938114 | -0.32932181 | 1.50531834 | 1.50524625 |

### Apparent cusp test

At fixed `N=120`, the fold non-degeneracy diagnostic has two zero brackets,
near

```text
(Ro,Tw,Hinf) = (-0.31698, 1.49020, -0.00680)
                (-0.41765, 1.60635, -0.00293).
```

They are **not converged finite cusp coordinates**.  The first apparent
primary-arm reversal migrates as

| `N` | `Ro` | `Tw` | `Hinf` |
|---:|---:|---:|---:|
| 80 | -0.32515 | 1.50016 | -1.7335e-2 |
| 120 | -0.31671 | 1.48987 | -6.9130e-3 |
| 160 | -0.31171 | 1.48374 | -3.8180e-3 |

The secondary near-tail structures migrate even more strongly.  Increasing
resolution pushes these reversals toward `Hinf=0`, exactly where the thermal
decay length diverges.  They are therefore classified as resolution-limited
tail structures, not as verified finite-core cusps.

### Relation to the thermal-tail boundary

The independently extrapolated low-temperature tail boundary uses
`Hinf=-0.02,-0.01,-0.005,-0.0025`.  On the production parameter-plane plot,
the outgoing fold arm approaches and then lies within about `1e-4` in wall
temperature of this boundary; the closest production comparison is

```text
Hinf = -1.4753e-3,
fold Tw - extrapolated tail Tw = -7.93e-5.
```

The reliable conclusion is thus:

- no resolution-independent finite cusp has been found;
- the converged finite-core fold arm moves toward `Hinf=0`;
- the apparent fold turns retreat into the singular tail layer as `N`
  increases;
- current evidence favours fold termination through far-field
  delocalisation rather than a finite-core cusp.

This is strong numerical evidence for a fold--tail connection, but the exact
`Hinf=0` connection coordinate is a singular limit and has not been solved as
a regular boundary-value problem.  A quadratic estimate from the last
`N=160` primary-arm points gives `(Ro,Tw)≈(-0.3072,1.4783)`, but it is retained
only as an exploratory location, not a manuscript-grade critical coordinate.

## Produced evidence

- `work/results/consistent_fold_tail_continuation/production_N120/`
- `work/results/consistent_fold_tail_continuation/convergence_N80/`
- `work/results/consistent_fold_tail_continuation/convergence_N160/`
- `work/results/consistent_fold_tail_continuation/primary_fold_level_convergence.csv`
- `work/results/consistent_fold_tail_continuation/apparent_cusp_migration.csv`
- `work/results/consistent_fold_tail_continuation/fold_tail_limit_estimate.txt`

The main parameter-plane plot is
`production_N120/fold_tail_parameter_plane.png`.  Cusp-bracket plots are
diagnostic only and must not be used to claim physical cusp coordinates.
