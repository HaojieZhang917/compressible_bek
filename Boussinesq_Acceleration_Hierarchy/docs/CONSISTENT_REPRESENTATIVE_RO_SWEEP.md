# Representative-Ro acceleration-consistent BEK sweep

## Objective

Before constructing a full two-parameter diagram, classify the heating branch
of the steady acceleration-consistent BEK model at

```text
Ro = -1, -0.75, -0.5, -0.25, -0.1, 0.25, 0.5, 1.
```

For each Rossby number the study asks only:

1. Is the branch connected to the isothermal solution single-valued in `Tw`?
2. Does the fixed-`Tw` steady Jacobian approach singularity?
3. Does `Hinf` approach zero, giving an infinite thermal tail?
4. Is there a verified branch termination before the positive-linear-density
   boundary `Tw=2`?

The result will be compared with the already completed disk-frame traditional
topology map.  It is a representative slice study, not yet a full
two-parameter fold/cusp continuation.

## Model and numerical definition

- Steady acceleration-consistent Boussinesq residual in
  `work/src/BEKConsistent.jl`.
- `Pr=0.72`, `gamma=beta_inf*T_inf=1`.
- Pier rational semi-infinite map with
  `N=120`, `(a,b,c)=(2,0.6,0.5)`.
- Branch starts from the independently converged isothermal BEK solution at
  each `Ro`.
- Fixed-`Tw` continuation uses adaptive steps no larger than `0.01` from
  `Tw=1` toward `Tw=1.99`.
- Newton tolerance is `1e-9` for the multi-`Ro` production sweep; stored
  solutions must have residual below this threshold.
- The fixed-`Tw` row-scaled Jacobian singular-value ratio is sampled along
  every branch.  A fold is not assigned from a projected `Hinf` extremum
  alone.
- If fixed-`Tw` continuation stops, the step is refined and the endpoint is
  classified only after Jacobian, `Hinf`, thermal length and, where needed,
  pseudo-arclength checks.
- `1 <= Tw <= 1.2` is the primary small-temperature interpretation range.
  Results for `1.2 < Tw < 2` are formal closure diagnostics.

## Classification rules

- `smooth_to_density_boundary`: fixed-`Tw` branch reaches `Tw=1.99`, no
  Jacobian collapse, and `Hinf` stays bounded away from zero.
- `thermal_tail`: `Hinf -> 0-`, the thermal length
  `1/(Pr*(-Hinf))` diverges and no earlier fold is found.
- `fold`: fixed-`Tw` Jacobian approaches a simple zero singular direction and
  pseudo-arclength shows `dTw/ds=0`.
- `other_termination`: a reproducible endpoint that is neither of the above.
- `unresolved`: numerical continuation stops before a mechanism is verified.

Single-valued means single-valued only on the isothermal-connected branch over
the sampled interval.  This sweep does not exclude disconnected solutions.
A cusp cannot be claimed from an individual fixed-`Ro` slice.

## Outputs

- Executable Julia study:
  `work/scripts/consistent_representative_ro_sweep.jl`.
- Results:
  `work/results/consistent_representative_ro_sweep/`.
- Final conclusions must be appended to
  `docs/RESEARCH_PROGRESS_LOG.md` before downstream use.

## Results

The representative production sweep gives

| `Ro` | consistent-model mechanism | endpoint / verified range | endpoint `Hinf` |
|---:|---|---:|---:|
| -1.00 | smooth to density boundary | verified to `Tw=1.99` | -0.89671408 |
| -0.75 | smooth to density boundary | verified to `Tw=1.99` | -0.84813769 |
| -0.50 | finite fold | `Tw_c=1.721704197` | -0.29758392 |
| -0.25 | infinite thermal tail | `Tw_tail≈1.405123` | 0 |
| -0.10 | infinite thermal tail | `Tw_tail≈1.182254` | 0 |
| +0.25 | smooth to density boundary | verified to `Tw=1.99` | -2.10851591 |
| +0.50 | smooth to density boundary | verified to `Tw=1.99` | -1.71026682 |
| +1.00 | smooth to density boundary | verified to `Tw=1.99` | -1.45829570 |

Here `density boundary` means the imposed linear ideal-gas law approaches
zero wall density as `Tw -> 2`; it is not a physical Boussinesq-accuracy
boundary.

### Fold verification at `Ro=-0.5`

Fixed-`Tw` continuation approaches a Jacobian singularity while
`Hinf≈-0.30`, well before a thermal-tail limit.  Full-state
pseudo-arclength continuation crosses the maximum in `Tw` and produces
negative `dTw` steps on the returned branch.  The interpolated fold is

```text
Tw_c = 1.721704197,
Hinf_c = -0.29758392.
```

Resolution results are

| `N` | `Tw_c` | `Hinf_c` |
|---:|---:|---:|
| 80 | 1.7217041979 | -0.29758138 |
| 120 | 1.7217041973 | -0.29758392 |
| 160 | 1.7217041945 | -0.29758340 |

The fold temperature is stable to about `4e-9` over this test.

### Thermal-tail verification

At `Ro=-0.25` and `-0.1`, `Hinf` approaches zero from below and the thermal
decay length grows to several hundred on the last resolved states.  No
earlier projected turning point is accompanied by a fixed-`Tw` fold.  The
quadratic `Tw(Hinf)` extrapolations are less precise than the fold coordinate:

| `Ro` | `N=80` | `N=120` | `N=160` |
|---:|---:|---:|---:|
| -0.25 | 1.405003 | 1.405123 | 1.404882 |
| -0.10 | 1.182239 | 1.182254 | 1.182281 |

These values classify the mechanism; a dedicated tail asymptotic calculation
is required before quoting a higher-precision endpoint.

## Interpretation

The consistent closure does not simply remove the traditional-model
singularities throughout the BEK family.  On the sampled slices it reorganises
them:

- the traditional von Karman fold is absent at `Ro=-1` in the consistent
  model;
- a different consistent-model fold exists at `Ro=-0.5`, where the
  traditional production map instead terminates through a thermal tail;
- tail termination persists at `Ro=-0.25` and `-0.1`;
- the sampled positive-`Ro` branches remain smooth through `Tw=1.99`, but
  `Ro=1` is no longer passive because the consistent model retains density
  coupling to the actual inertial acceleration.

Thus the main topology claim should be closure-dependent relocation and
conversion of branch boundaries, not universal removal of all folds.  The
transitions between the smooth, fold and tail slices now define the targets
for the two-model comparison and subsequent two-parameter continuation.

