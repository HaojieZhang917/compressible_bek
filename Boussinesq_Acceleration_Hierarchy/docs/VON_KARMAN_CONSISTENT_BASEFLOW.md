# Acceleration-consistent von Karman base-flow study

## Objective

Compute the heated von Karman similarity branch for the steady,
acceleration-consistent Boussinesq endpoint and determine whether the branch
connected to the isothermal solution contains a wall-temperature saddle-node.
The calculation also distinguishes an extremum of the diagnostic observable
`Hinf(Tw)` from a true fold in the control parameter `Tw`.

A cusp is a codimension-two event and cannot be established on the single
slice `Ro=-1`, `gamma=1`.  This study therefore reports whether a fold exists
on that slice; a later two-parameter continuation in `(Ro,Tw)` or
`(gamma,Tw)` is required to locate or exclude a cusp of the full family.

## Model

Let

```text
gamma = beta_inf*T_inf,
chi(T) = 1 - gamma*(T-1).
```

For the ideal-gas linearization used here, `gamma=1`.  At `Ro=-1`, the
inertial azimuthal velocity is proportional to `1-G`, and the steady
acceleration-consistent similarity equations are

```text
H' + 2F = 0,
F'' + chi(T)*(-F^2 - H*F' + (1-G)^2) = 0,
G'' + chi(T)*(-2F*G - H*G' + 2F) = 0,
T'' - Pr*H*T' = 0.
```

Boundary conditions are

```text
H(0)=F(0)=G(0)=0,  T(0)=Tw,
F(infinity)=0,     G(infinity)=1, T(infinity)=1.
```

The more general nonzero-`Ro` residual implemented for regression is

```text
s = Co/2, f = s+Ro, Co=2-Ro-Ro^2,
A_r = Ro*(F^2+H*F') - (s+Ro*G)^2/Ro,
A_theta = Ro*(2F*G+H*G') + Co*F,
F'' + chi*A_r + f^2/Ro = 0,
G'' + chi*A_theta = 0.
```

At `T=1`, these equations reduce exactly to the isothermal BEK equations.
This is a fully self-consistent steady endpoint, not a frozen-base-flow or
term-toggle diagnostic.  It remains a linear-density, constant-property
Boussinesq model rather than a compressible model.

## Numerical definition

- `Ro=-1`, `Pr=0.72`, `gamma=1`.
- Pier rational semi-infinite map with
  `N=120`, `(a,b,c)=(2,0.6,0.5)`.
- Fixed-wall-temperature continuation from `Tw=1` to `Tw=1.99`.
- Fine physical-range sampling for `1 <= Tw <= 1.2`; continuation above this
  range is a formal closure diagnostic only.
- Newton residual tolerance `1e-10`.
- A row-scaled fixed-`Tw` Jacobian singular-value ratio is evaluated along the
  branch.  A true wall-temperature fold requires this Jacobian to become
  singular and a pseudo-arclength tangent to satisfy `dTw/ds=0`.
- Selected profiles: `Tw=1,1.05,1.1,1.2,1.5,1.8,1.95,1.99`.

## Source and output dependencies

- Isothermal operators: `work/src/BEKIsothermal.jl`.
- New consistent residual and continuation:
  `work/src/BEKConsistent.jl`.
- Executable study:
  `work/scripts/vonkarman_consistent_baseflow.jl`.
- Outputs:
  `work/results/vonkarman_consistent_baseflow/`.
- Preserved finite-domain Lopez values under `baselines/rotating_disk/` are
  used only as an independent small-temperature regression check and are not
  overwritten.

## Results

The isothermal-connected branch was continued without a fixed-`Tw` failure
from `Tw=1` through `Tw=1.99`.  Selected production values are

| `Tw` | `Hinf` | `Fmax` | `F'(0)` | `G'(0)` |
|---:|---:|---:|---:|---:|
| 1.00 | -0.8844741102 | 0.1807006018 | 0.5102326189 | 0.6159220144 |
| 1.05 | -0.8879786240 | 0.1797738824 | 0.4968556672 | 0.6038746123 |
| 1.20 | -0.8974590313 | 0.1763750272 | 0.4548232078 | 0.5669994273 |
| 1.50 | -0.9101888204 | 0.1671767264 | 0.3603475692 | 0.4897062720 |
| 1.80 | -0.9095332845 | 0.1533597104 | 0.2462495563 | 0.4070412060 |
| 1.99 | -0.8967140780 | 0.1413379678 | 0.1578824434 | 0.3512316373 |

All production residuals are below `1e-10`.  The values at `Tw=1.2` and
`Tw=1.5` reproduce the preserved finite-domain Lopez results to the displayed
digits.

`Hinf(Tw)` has one smooth minimum.  A quadratic interpolation through the
three nearest production samples gives

```text
Tw = 1.64704657,
Hinf = -0.9120805930.
```

This point is not a wall-temperature saddle-node: the branch remains a
single-valued function of `Tw`, the fixed-`Tw` Jacobian remains nonsingular,
and the pseudo-arclength continuation crosses the `Hinf` minimum with `Tw`
still increasing.  Over the full pseudo-arclength trace the smallest forward
wall-temperature increment is `6.84098e-3`; no `dTw/ds=0` point is detected.

The row-scaled fixed-`Tw` singular-value ratio decreases smoothly from
`7.70e-7` at `Tw=1` to `5.77e-7` at `Tw=1.99`.  This mild change is unlike the
collapse to approximately machine precision at the traditional-model fold.
The absolute ratio is resolution dependent, so fold rejection is based on
its smooth within-resolution behaviour, the successful fixed-`Tw`
continuation and the independent pseudo-arclength tangent check together.

Resolution checks give

| `Tw` | quantity | `N=80` | `N=120` | `N=160` |
|---:|---|---:|---:|---:|
| 1.2 | `Hinf` | -0.8974590309 | -0.8974590314 | -0.8974590314 |
| 1.2 | `F'(0)` | 0.4548232078 | 0.4548232078 | 0.4548232078 |
| 1.8 | `Hinf` | -0.9095332838 | -0.9095332845 | -0.9095332845 |
| 1.8 | `F'(0)` | 0.2462495563 | 0.2462495563 | 0.2462495563 |

## Bifurcation conclusion

No saddle-node is found on the heated, isothermal-connected
acceleration-consistent von Karman branch for
`1 <= Tw <= 1.99`.  The interval above the small-temperature Boussinesq range
is only a formal model diagnostic; `Tw=2` is the zero-wall-density boundary
of the adopted linear ideal-gas law.  Disconnected solution components have
not been excluded by continuation from the isothermal state.

There is no cusp on this one-parameter slice.  A cusp claim requires
two-parameter fold continuation.  Since this `Ro=-1`, `gamma=1` slice contains
no fold, it supplies no cusp candidate by itself.
