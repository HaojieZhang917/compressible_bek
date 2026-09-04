# Systematic thermal availability study

This study tests whether the traditional Boussinesq BEK endpoint changes
continuation mechanism as `Ro` varies. It uses the corrected coefficient
`Lambda_cf=Co^2/(4Ro)`, `Co=2-Ro-Ro^2`, `N=60` for exploratory coverage,
Pier parameters `(a,b,c)=(2,0.6,0.5)`, and `Pr=0.72`. Each branch starts from
the `Tw=1` isothermal solution and is continued at fixed `Hinf` toward
`Hinf=-0.001`. Fold extrema and the final continuation status are recorded.
The lower order is for mechanism mapping; endpoint values are to be checked
at `N=120` before quantitative claims.
