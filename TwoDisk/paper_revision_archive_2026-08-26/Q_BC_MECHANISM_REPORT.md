# Diagnostic report: mathematical origin of the non-monotonic receptivity coefficient

## Scope

This is a stand-alone numerical and mathematical diagnostic for

$$
C_r=-\mathrm{i}\frac{\mathrm{B.C.}}{Q}.
$$

It has not been inserted into the manuscript or the consolidated response
letter.

## 1. Exact extremum condition

For nonzero $Q$ and $\mathrm{B.C.}$,

$$
\frac{\mathrm{d}}{\mathrm{d}R}\ln|C_r|
=
\frac{\mathrm{d}}{\mathrm{d}R}\ln|\mathrm{B.C.}|
-
\frac{\mathrm{d}}{\mathrm{d}R}\ln|Q|.
$$

Therefore, a stationary point of $|C_r|$ satisfies

$$
\frac{\mathrm{d}}{\mathrm{d}R}\ln|\mathrm{B.C.}|
=
\frac{\mathrm{d}}{\mathrm{d}R}\ln|Q|.
$$

It does not generally coincide with either a maximum of
$|\mathrm{B.C.}|$ or a minimum of $|Q|$.

For the continuously tracked unstable branch:

| Quantity | Extremum radius | Extremum value |
|---|---:|---:|
| $|\mathrm{B.C.}|$ | $R\simeq389.01$ (maximum) | $1.01641\times10^{-4}$ |
| $|C_r|$ | $R\simeq468.93$ (maximum) | $0.244415$ |
| $|Q|$ | $R\simeq470.00$ (minimum) | $3.57755\times10^{-4}$ |

At the $C_r$ maximum, both numerator and denominator are still decreasing,
with nearly equal logarithmic slopes:

$$
\frac{\mathrm{d}}{\mathrm{d}R}\ln|\mathrm{B.C.}|
\simeq
\frac{\mathrm{d}}{\mathrm{d}R}\ln|Q|
\simeq-1.4\times10^{-3}.
$$

The $C_r$ maximum consequently occurs approximately one radius unit upstream
of the $Q$ minimum. At the $Q$ minimum, the numerator is already decreasing,
so $C_r$ has already begun to fall.

## 2. Which term controls the rise and fall?

The radial evolution has three regimes:

1. From the lower neutral point to $R\simeq389$, $|\mathrm{B.C.}|$ increases
   by approximately $25.9\%$ while $|Q|$ decreases by approximately $9.3\%$.
   Both effects increase $|C_r|$.
2. From $R\simeq389$ to $R\simeq469$, $|\mathrm{B.C.}|$ decreases by
   approximately $13.9\%$, but $|Q|$ decreases by approximately $29.9\%$.
   The denominator decrease therefore overcomes the numerator decrease and
   produces the renewed rise of $|C_r|$.
3. From $R=470$ to $R=500$, $|Q|$ increases by approximately $47.2\%$ and
   $|\mathrm{B.C.}|$ decreases by approximately $3.7\%$, causing
   $|C_r|$ to decrease by approximately $34.6\%$.

Thus, the downstream enhancement and turnover are predominantly
denominator-controlled. The numerator modulates the response and shifts the
$C_r$ maximum upstream of the $Q$ minimum.

The previously described local $C_r$ maximum near the maximum-growth
location is not present on the recomputed continuous branch. The curve has a
shoulder, where its positive slope becomes small, followed by renewed growth.
This shoulder is produced by partial cancellation between the numerator and
denominator logarithmic slopes.

## 3. What causes the minimum of Q?

For the polynomial eigenvalue problem

$$
P(\alpha)\phi
=(L_0+\alpha L_1+\alpha^2L_2)\phi=0,
$$

the denominator is

$$
Q
=
\left\langle\phi^\dagger,
P_\alpha(\alpha)\phi\right\rangle
=
\left\langle\phi^\dagger,
(L_1+2\alpha L_2)\phi\right\rangle.
$$

Splitting $Q=Q_1+Q_2$ with

$$
Q_1=\langle\phi^\dagger,L_1\phi\rangle,
\qquad
Q_2=\langle\phi^\dagger,2\alpha L_2\phi\rangle,
$$

gives, at $R=470$,

$$
|Q_1|=3.6140\times10^{-4},\qquad
|Q_2|=5.8956\times10^{-6},
$$

and

$$
\frac{|Q|}{|Q_1|+|Q_2|}=0.974.
$$

Hence, the $Q$ valley is not caused by cancellation between the linear and
quadratic eigenvalue terms; the quadratic contribution is only about
$1.6\%$ of the linear contribution.

A scale-independent alignment diagnostic is

$$
\chi_Q
=
\frac{|Q|}
{\|\phi^\dagger\|_M
 \|P_\alpha(\alpha)\phi\|_M}.
$$

For the unstable branch, $\chi_Q$ decreases from $0.02270$ at $R=420$ to a
minimum of $0.01710$ at $R=470$, and subsequently recovers. Its reciprocal
increases from $44.0$ to $58.5$. This demonstrates that the $Q$ valley is not
merely an arbitrary eigenfunction-scaling effect. It corresponds to weakened
biorthogonal alignment and moderately increased eigenvalue sensitivity.

Within the adopted local $L_2$ metric, the dominant $Q_1$ term is itself the
residual of complex contributions from the different equation/state blocks.
Their cancellation strengthens near $R=470$. This is a plausible algebraic
source of the reduced biorthogonal pairing. Because this block decomposition
depends on the chosen metric and variable scaling, it should be presented as
a representation-dependent diagnostic, not a universal physical invariant.

The complete complex eigenvalues of the two tracked branches remain separated
by at least approximately $3.10\times10^{-2}$ over $420\leq R\leq500$.
Therefore, the reduced $Q$ is not an exceptional point, eigenvalue
coalescence, or pole singularity.

## 4. What happens in the boundary projection?

The boundary term can be decomposed into its radial- and azimuthal-velocity
contributions:

$$
\mathrm{B.C.}=\mathrm{B.C.}_u+\mathrm{B.C.}_v.
$$

At $R=470$,

$$
|\mathrm{B.C.}_u|=2.4004\times10^{-4},\qquad
|\mathrm{B.C.}_v|=1.5267\times10^{-4},
$$

while

$$
|\mathrm{B.C.}|=8.7384\times10^{-5}.
$$

Their phase difference is approximately $179.6^\circ$, giving the
cancellation index

$$
\eta_{\mathrm{BC}}
=
\frac{|\mathrm{B.C.}_u+\mathrm{B.C.}_v|}
{|\mathrm{B.C.}_u|+|\mathrm{B.C.}_v|}
=0.223.
$$

Thus, the numerator contains strong destructive interference. Nevertheless,
its total magnitude varies smoothly and decreases throughout the
$R=420$--$500$ interval. This cancellation suppresses the response and helps
shift the $C_r$ maximum upstream, but it does not generate the downstream
rise or sharp turnover.

## 5. Numerical robustness and limitations

- A refined radial scan confirms the ordering
  $R_{C_r,\max}<R_{Q,\min}$.
- At $N=99$, the refined values are approximately
  $R_{C_r,\max}=468.93$ and $R_{Q,\min}=470.03$.
- At $N=129$, the corresponding values are approximately
  $468.71$ and $470.18$; the peak $C_r$ differs by approximately $0.31\%$.
- A continuation calculation using multiple candidate roots and
  eigenfunction overlap follows two distinct branches throughout
  $420\leq R\leq500$ and reproduces the extrema above. The Hermitian-adjoint
  eigenvalue pairing errors remain below $8\times10^{-12}$.

The earlier fixed-shift two-root script solved each radius independently and
could jump to the wrong root near the end of the interval. It should not be
used as the final evidence for a Type-I/Type-II mode conversion. The
continuation result is sufficient for the $Q$--B.C. mechanism, but a strict
physical Type-I/Type-II identity assignment still requires full-mode overlap
with unambiguous reference modes.

## 6. Defensible interpretation

The strongest conclusion currently supported is:

> The non-monotonic receptivity coefficient is governed by the competition
> between the logarithmic variations of the forcing projection and the
> polynomial-eigenvalue normalization denominator. The downstream enhancement
> is predominantly associated with a finite reduction of the direct-adjoint
> biorthogonal pairing represented by $Q$, whereas destructive interference
> between the radial and azimuthal boundary projections suppresses and shifts
> the response peak. The behavior is a finite non-normal conditioning effect,
> not an eigenvalue coalescence.

The coincidence of this conditioning feature with the two-branch
near-real-part crossing may be reported as a correlation. A causal statement
that mode conversion produces the $Q$ minimum would require additional
mode-overlap/rotation evidence.
