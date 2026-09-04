# Traditional Boussinesq BEK solvability study

## Question

Determine whether the finite-temperature saddle nodes near the von Karman
limit and the `Hinf -> 0-` long-tail limit away from that endpoint are
separate phenomena or different features of one connected steady-solution
surface.

## Model and numerical specification

The study uses the traditional, disk-frame centrifugal Boussinesq reduction

```text
Lambda_cf(Ro) = Co^2/(4 Ro),  Co = 2 - Ro - Ro^2,
```

in the radial BEK residual. The convention is valid only for `Ro < 0` in
this sweep; the chosen differential-rotation scaling is singular at `Ro=0`.
All production points use the Pier semi-infinite mapping with
`N=120, (a,b,c)=(2,0.6,0.5), Pr=0.72`.

For each selected negative `Ro`, the branch connected to the isothermal
solution is continued at fixed `Hinf`. Each converged state records `Tw`,
the fixed-`Tw` Jacobian singular-value indicator, and the exact far-field
thermal e-fold length

```text
ell_T = 1/(Pr*(-Hinf)),  Hinf < 0.
```

The latter follows by linearizing `T''-Pr*H*T'=0` at infinity. A failure of
fixed-`Hinf` continuation is reported as a coordinate/conditioning limit,
not by itself as a proof of nonexistence. Fold claims require extrema in the
continued branch together with the fixed-`Tw` Jacobian diagnostic.

### Near-von-Karman high-tail connection test

To separate a genuine branch termination from a poor continuation
coordinate, a second calculation holds `Hinf=-0.02` fixed and starts from
the converged long-tail state at `Ro=-0.995`. It continues in `Ro` with a
full-state pseudo-arclength constraint: the collocation fields, `Tw`, and
`Ro` all enter the tangent. The field block is normalized by the square root
of its number of degrees of freedom, so it cannot dominate the scalar
parameters purely because of vector dimension. This is a topological
diagnostic for the semi-infinite collocation equations; it does not replace
a bordered fold/cusp calculation, and it is not a physical
Boussinesq-validity test.

## Interpretation criteria

* A saddle node is a finite-`Hinf` extremum of `Tw` with a near-singular
  fixed-`Tw` Jacobian.
* A long-tail endpoint is approached when `Hinf -> 0-` and `ell_T` grows
  without a finite-`Hinf` `Tw` extremum.
* The two are globally connected only if they occur on the same continued
  branch; they remain locally distinct because the former is a Jacobian
  singularity whereas the latter loses the exponentially decaying thermal
  far-field mode.

## Computed evidence

The production `N=120` scan over negative `Ro` found the following branch
pattern. The endpoint `Hinf=-0.02` is used only as a common diagnostic point,
not as a finite-domain truncation.

| `Ro` | status of isothermal-connected branch | `Tw` at `Hinf=-0.02` or stop | interpretation |
|---:|---|---:|---|
| -1.000 | two established folds, then returned branch | `1.03188` | finite folds plus tail approach |
| -0.999 | two folds | branch tracking becomes ill-conditioned after the second fold | near-cusp projection |
| -0.998 | two folds | branch tracking becomes ill-conditioned after the second fold | near-cusp projection |
| -0.994 | close fold pair | branch tracking becomes ill-conditioned | near-merger region |
| -0.990 | no resolved finite fold | fixed-`Hinf` continuation stops near `Hinf=-0.087`, `Tw=1.312` | continuation conditioning, not proof of nonexistence |
| -0.950 | no resolved finite fold | `Tw=1.87` at `Hinf=-0.02` | long-tail branch |
| -0.900 | no resolved finite fold | `Tw=1.90` at `Hinf=-0.02` | long-tail branch |
| -0.750 | no resolved finite fold | `Tw=1.86` at `Hinf=-0.02` | long-tail branch |
| -0.500 | no resolved finite fold | `Tw=1.674` at `Hinf=-0.02` | long-tail branch |

Higher-order spot checks at fixed `Hinf=-0.02` for `Ro=-0.5` gave
`Tw=1.67412` for `N=80,120,160`, demonstrating that the long-tail branch
value is spectrally converged at this diagnostic distance from zero. At
`Ro=-0.5`, continuation to `Hinf=-0.0001` gave `Tw≈1.68921` and
`ell_T≈1.39e4`; the slow convergence in this final limit is expected because
the thermal decay length diverges.

The far-field rates quantify the mechanism change. At the first and second
von Karman folds, `(alpha_T,alpha_V)` are approximately `(0.384,0.533)` and
`(0.0816,0.113)`, respectively: thermal and hydrodynamic tails are comparable.
At `Hinf=-0.02`, `Ro=-0.95` gives `(0.0144,0.281)`, while `Ro=-0.5` gives
`(0.0144,0.796)`. Thus the tail away from the von Karman endpoint is strongly
thermal and is no longer a joint velocity-temperature slow mode.

For `Ro=-1`, the two trusted infinite-mapping folds remain
`(-0.53276,1.04802)` and `(-0.11331,1.03014)`. Along the returned branch,
`Hinf=-0.05` gives `ell_T≈27.8`; this branch continues toward the long-tail
regime after the finite folds. This is direct evidence that the two
termination mechanisms occur on the same steady solution branch at the
von Karman endpoint, while the finite fold pair becomes difficult to follow
as `Ro` moves away from `-1`.

The full-state fixed-`Hinf=-0.02` pseudo-arclength test provides an additional
cross-section of the same surface (`work/results/traditional_bek_solvability/
long_tail_ro_fullstate_toward_vk.csv`). Starting from the converged high-tail
branch, it reaches `Ro=-0.99996896` with `Tw=1.04153646` and residual
`3.2×10^-9`; the path is monotone in `Tw` over the last resolved segment and
remains converged while `Ro` approaches `-1`. The final step is limited by
Newton conditioning near the endpoint, not by a detected residual blow-up or
loss of the solution. This supports a connected high-tail sheet approaching
the von Karman fold-controlled region, but it does not prove that the chosen
fixed-`Hinf` slice intersects the exact returned branch at `Ro=-1`.

## Applicability conclusion

The traditional closure is not uniformly applicable over the BEK family at
arbitrary `Tw`. Its practical domain is controlled by two boundaries:

1. near `Ro=-1`, a finite-temperature saddle-node boundary;
2. away from `Ro=-1`, a singular long-tail boundary where `Hinf -> 0-` and
   the thermal far-field mode loses exponential decay.

The calculations support a connected evolution from the fold-controlled
region to the tail-controlled region, but they do not yet prove a smooth
codimension-two cusp curve. Such a claim requires a bordered fold solve and
an explicit asymptotic continuation at `Hinf=0`. The current evidence is
therefore a mechanism map and an applicability diagnosis, not a universal
existence theorem.

The values near `Tw=1.7-1.9` are formal solvability diagnostics, not
quantitative Boussinesq predictions: with temperature nondimensionalized by
the ambient value, the linear-density approximation requires `|Tw-1| << 1`.
