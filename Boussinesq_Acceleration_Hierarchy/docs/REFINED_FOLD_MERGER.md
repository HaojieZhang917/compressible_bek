# Refined fold-merger study

The corrected traditional thermal BEK coefficient
`Lambda_cf=Co^2/(4Ro)` is used with `N=120`, Pier parameters
`(a,b,c)=(2,0.6,0.5)` and `Pr=0.72`. The fold merger is resolved by fixed-
`Hinf` continuation with an initial step `5e-4` (down to `2e-5` on failure).
The scan uses `Delta Ro=0.001` and `0.0001` over `[-1,-0.99]`.

Resolved pairs include `Ro=-0.999` at `Hinf=-0.5291503,-0.2093312`,
`Ro=-0.998` at `Hinf=-0.5243172,-0.2590517`, and `Ro=-0.994` at
`Hinf=-0.4683198,-0.4180556` (`Delta Hinf=0.0502642`). At intermediate
points the single-branch scan detects only one extremum, indicating branch
tracking loss near the merger rather than physical reappearance. A bordered
fold solve is required for the final cusp coordinate.
