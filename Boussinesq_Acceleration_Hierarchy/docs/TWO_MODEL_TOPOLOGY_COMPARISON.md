# Unified traditional/consistent topology comparison

## Scope

This study compares the isothermal-connected heating branches of the
disk-frame traditional and acceleration-consistent steady BEK models at the
same eight Rossby numbers

```text
Ro = -1, -0.75, -0.5, -0.25, -0.1, 0.25, 0.5, 1.
```

It closes the one-parameter comparison before any two-parameter continuation.
It does not add further isolated `Ro` slices and does not address stability,
low-Mach corrections or additional material-property parameters.

## Common numerical definition

- `Pr=0.72`, `gamma=1` for the consistent model and the corresponding
  ideal-gas linear coefficient for the traditional model.
- Pier semi-infinite rational map with
  `N=120`, `(a,b,c)=(2,0.6,0.5)`.
- Newton tolerance `1e-9`.
- Fixed-`Tw` continuation starts at the independently converged isothermal
  solution, uses an adaptive step no larger than `0.01`, and stops at a
  verified boundary or `Tw=1.99`.
- A finite fold requires both a fixed-`Tw` Jacobian singularity and a
  pseudo-arclength crossing with a sign change of `dTw/ds`.
- A thermal tail requires `Hinf -> 0-` and divergence of
  `-1/(Pr*Hinf)` without an earlier verified fold.
- `passive` is reserved for the traditional `Ro=1` endpoint, where the
  traditional thermal forcing coefficient vanishes exactly.
- `smooth` means only that the branch remains regular through `Tw=1.99`;
  this imposed density boundary is not a physical Boussinesq-validity limit.

## Reported quantities

Each model/`Ro` row records

```text
topology, Tw_c, Hinf_c, sigma_min(J_scaled), sigma_min/sigma_max, dTw/ds.
```

`J_scaled` is the fixed-`Tw` Jacobian after unit row-norm scaling.  This makes
the diagnostic reproducible under the fixed discretisation, but it remains a
discretisation diagnostic rather than a dimensional physical quantity.

For a fold, all critical quantities are evaluated at a fixed-`Hinf` refined
fold and `dTw/ds` is a centred full-state arclength derivative.  For a tail,
`Tw_c` is a quadratic `Tw(Hinf)` extrapolation to `Hinf=0`; the singular value
is necessarily reported at the last resolved finite-tail state and
`dTw/ds` is not a fold diagnostic, so it is recorded as not applicable.  For
smooth/passive rows, the values refer to `Tw=1.99`.

## Source and outputs

- Traditional residual: `work/src/BEKThermal.jl`.
- Consistent residual: `work/src/BEKConsistent.jl`.
- Executable study: `work/scripts/compare_traditional_consistent_topology.jl`.
- New results: `work/results/two_model_topology_comparison/`.

The completed conclusions must be appended to
`docs/RESEARCH_PROGRESS_LOG.md` before manuscript use.

## Results

The unified classification is

| `Ro` | traditional | consistent | closure-induced change |
|---:|---|---|---|
| -1.00 | fold | smooth | traditional fold removed |
| -0.75 | tail | smooth | delocalisation boundary removed from the tested interval |
| -0.50 | tail | fold | finite fold created/relocated |
| -0.25 | tail | tail | mechanism retained |
| -0.10 | tail | tail | mechanism retained |
| +0.25 | smooth | smooth | no topology change through `Tw=1.99` |
| +0.50 | smooth | smooth | no topology change through `Tw=1.99` |
| +1.00 | passive | smooth | thermal--momentum coupling activated |

The two verified folds are

| model | `Ro` | `Tw_c` | `Hinf_c` | scaled `sigma_min` | numerical `dTw/ds` |
|---|---:|---:|---:|---:|---:|
| traditional | -1.0 | 1.0480217310 | -0.53275559 | `2.94e-9` | `-4.81e-6` |
| consistent | -0.5 | 1.7217041955 | -0.29758769 | `6.85e-11` | `2.18e-7` |

The small nonzero arclength derivatives are centred numerical estimates at
the fixed-`Hinf` refined vertices; the pseudo-arclength branches contain
opposite signs on the two sides.  The traditional fold temperature agrees
with the preserved baseline `1.048021731` and the consistent value agrees
with the independent `N=80,120,160` result.

For tail rows the complete endpoint and last-resolved Jacobian diagnostics
are in `work/results/two_model_topology_comparison/topology_table.csv`.
Those singular values must not be interpreted as fold singular values:
their location is a finite negative `Hinf`, while `Tw_c` is extrapolated to
the non-localised `Hinf=0` limit.

The comparison establishes that the closure does not monotonically
regularise the BEK family.  It reorganises where finite-core folding,
far-field delocalisation, smooth continuation and passive thermal response
occur.
