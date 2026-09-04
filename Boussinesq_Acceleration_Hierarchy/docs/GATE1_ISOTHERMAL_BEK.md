# Gate 1: isothermal BEK extension

## Scope

This gate implements the steady, isothermal Bodewadt--Ekman--von Karman
(BEK) similarity equations on a rational Chebyshev semi-infinite mapping. No
temperature equation or wall-temperature continuation is included yet. The
purpose of this stage is to establish the `Ro`-dependent base-flow operator and
its von Karman regression before adding the thermal Boussinesq closure.

## Source equations and convention

Lingwood (1997), equations (3.11)--(3.14), uses the dimensional velocity
scaling in which the radial and axial similarity functions are `U` and `W`.
The implementation uses the existing von Karman convention

```text
U = -F,  W = -H,  V = G.
```

With this convention the isothermal equations are

```text
H' + 2 F = 0,
F'' + Ro*(F^2 + H*F' - (G^2 - 1)) - Co*(G - 1) = 0,
G'' + Ro*(2 F G + H*G') + Co*F = 0,
Co = 2 - Ro - Ro^2.
```

The boundary conditions are `H(0)=F(0)=G(0)=0`, `F(infinity)=0`,
`G(infinity)=1`, and `H'(infinity)=0`. The last condition replaces the
redundant far-field continuity row in the collocation system.

At `Ro=-1`, `Co=2`, the radial equation becomes

```text
F'' - F^2 - H*F' + (G - 1)^2 = 0,
```

and the azimuthal equation becomes

```text
G'' - 2 F G - H*G' + 2 F = 0,
```

which are the current isothermal von Karman equations in
`baselines/saddle_node/src/VonKarman.jl`.

## Numerical record

- Mapping: Pier (2003) rational map
  `z=a*(1+xi')/(1-xi')`, `xi'=b*xi+(1-b)*(xi^3+c*(1-xi^2))`,
  with `a=2`, `b=0.6`, `c=0.5`.
- Collocation degree: `N=120` for the `Ro` scan and subsequent calculations.
- Sampled `Ro`: `[-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1]`.
- Nonlinear tolerance: `2e-10` in the supplied Gate 1 script.
- Output: `work/results/gate1_isothermal/`.
- Code: `work/src/BEKIsothermal.jl` and
  `work/scripts/gate1_bek_isothermal.jl`.

The scan is a base-flow verification only. It does not identify a thermal
saddle-node because `Tw` and temperature are intentionally absent at this
gate. Saddle-node evolution in `(Ro,Tw)` requires the next thermal extension.

## Thermal validation update

The traditional centrifugal Boussinesq radial term and the isothermal energy
equation are implemented in `work/src/BEKThermal.jl`. At `Ro=-1`, continuation
in `H(infinity)` with `N=120`, `a=2`, `b=0.6`, `c=0.5` gives the local fold fit

```text
Hinf = -0.532760583932
Tw   = 1.048021731462
```

The established baseline is `Tw=1.0480217312207`; the difference is about
`2.4e-10`. This is a solver validation at the von Karman endpoint, not yet a
BEK `Ro`-dependent fold map.

## Second-fold cross-validation

The legacy `branch_diagram.png` in the external `Vonkarmen_bone` directory
uses a finite `zmax=20` domain. Its second turning point (`Hinf` about
`-0.242`, `Tw` about `1.0406`) is domain-sensitive: the accompanying report
shows it drifting toward `Hinf` about `-0.113` as `zmax` is increased.

Using the present infinite-endpoint implementation, a dense `N=80` continuation
gives `Hinf=-0.114` and `Tw=1.03014422`. A local `N=120` continuation seeded
from the prior infinite-domain profile gives

```text
Hinf = -0.113323594317
Tw   = 1.030145508131
```

The trusted prior infinite-mapping summary is
`Hinf=-0.113305731011` and `Tw=1.030144663889`. The remaining difference is
consistent with the relaxed local Newton tolerance and seed interpolation;
the agreement in `Tw` is about `8.4e-7`. The finite `zmax=20` second fold must
therefore not be used as the infinite-domain result.

## Verification requirements

1. Run the script with Julia 1.10 or newer.
2. Check every reported collocation residual is below `2e-10`.
3. Repeat `Ro=-1` at increased degree and confirm the fields and residual
   converge to the existing von Karman reference.
4. Repeat selected `Ro` values with a second mapping scale before using any
   fold or topology claim.

The map places `z=infinity` at an exact collocation endpoint; the finite
von Karman solution is used only as a Newton seed and is never imposed as a
far-field boundary. The first run completed with residuals below `2.2e-10`
for all sampled `Ro`; at `Ro=-1` the residual was `4.1e-11` and the computed
`H(infinity)` was `-0.88447411`. The coarse nearest-node profile check gave a
maximum difference `4.1e-3`; this is an interpolation/seed diagnostic, not a
converged error estimate, and must be replaced by a degree/mapping convergence
comparison before making a precision claim. No baseline files are overwritten.

## Source dependency

R. J. Lingwood, “Absolute instability of the Ekman layer and related rotating
flows,” *Journal of Fluid Mechanics* **331** (1997), 405--428, equations
(3.11)--(3.16). The attached PDF is the local source used for transcription.
