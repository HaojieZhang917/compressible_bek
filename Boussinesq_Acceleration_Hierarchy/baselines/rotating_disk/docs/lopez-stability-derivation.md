# Lopez generalized-Boussinesq stability equations for a rotating disk

## 1. Scope and conventions

This document applies the inertial-frame model proposed by Lopez, Marques & Avila
(2013) to the rotating-disk boundary layer used in this workspace. Lopez et al.
studied rotating Taylor-Couette flow, not the rotating disk itself; the equations
below are therefore a derivation for the present geometry, not equations copied
from their paper.

Temperature is nondimensionalized by the free-stream temperature. For the
ideal-gas linear Boussinesq closure used here,

$$
\chi(T)=\frac{\rho}{\rho_\infty}\simeq 1-(T-1)=2-T,
\qquad \delta\chi=-\vartheta,
$$

where $T(z)$ is the basic temperature and $\vartheta$ is its infinitesimal
disturbance. There is no independent `Lambda` parameter.

The Lopez basic-flow file stores the azimuthal deficit $G_L$, with
$G_L(0)=0$ and $G_L(\infty)=1$. Define the inertial azimuthal velocity profile

$$
Q(z)=1-G_L(z),\qquad Q(0)=1,\quad Q(\infty)=0.
$$

If a disk-fixed code uses a relative azimuthal profile, its input is
$G_R=Q-1=-G_L$. Feeding the raw Lopez $G_L$ into the existing stability code
without this conversion reverses the physical convention.

## 2. Governing approximation

Lopez, Marques & Avila (2013) start from the inertial-frame Boussinesq equation
and retain the density fluctuation in the nonlinear advection acceleration. With
gravity omitted for the rotating-disk problem, their equation (2.12) is

$$
\rho_0(\partial_t+\boldsymbol V\cdot\nabla)\boldsymbol V
=-\nabla p^*+\mu\nabla^2\boldsymbol V
-\rho'(\boldsymbol V\cdot\nabla)\boldsymbol V.
$$

Moving the last term to the left gives

$$
\rho_0\partial_t\boldsymbol V
+(\rho_0+\rho')(\boldsymbol V\cdot\nabla)\boldsymbol V
=-\nabla p^*+\mu\nabla^2\boldsymbol V.
$$

This step fixes which accelerations are density weighted: the nonlinear
advection is weighted by the complete linearized density, whereas the local
time derivative is weighted by the reference density. For

$$
\frac{\rho'}{\rho_0}=-\beta_T(T_d-T_\infty),
\qquad T=\frac{T_d}{T_\infty},
$$

the ideal-gas reference value $\beta_TT_\infty=1$ gives

$$
\chi(T)=1+\frac{\rho'}{\rho_0}=2-T,
\qquad \delta\chi=-\vartheta.
$$

After division by the reference density, the dimensionless system used here is

$$
\nabla\cdot\boldsymbol V=0,
$$

$$
\partial_t\boldsymbol V
+\chi(T)(\boldsymbol V\cdot\nabla)\boldsymbol V
=-\nabla P+\frac{1}{R}\nabla^2\boldsymbol V,
$$

$$
\partial_tT+\boldsymbol V\cdot\nabla T
=\frac{1}{Pr\,R}\nabla^2T.
$$

The coefficient $1/R$ follows the local Reynolds-number scaling already used
by the legacy `archive/legacy_scripts/Stability.jl`. Most importantly, $\chi$ multiplies the nonlinear advection
term, but not the local time derivative, pressure gradient, or viscous term.
This is the specific 2013 Lopez model. A later canonical generalized
Boussinesq model that multiplies every local acceleration is a different model
and has a different mass matrix.

The momentum equation is non-conservative. Consequently no $\chi'$ terms are
created by differentiating the momentum flux.

This is not the same as replacing the complete material derivative by
$\chi(\partial_t+\boldsymbol V\cdot\nabla)\boldsymbol V$. Such a replacement
would put $\chi$ in the temporal mass matrix and would define a different
generalized-Boussinesq model.

## 3. Basic flow

At the local radius $r=R$, use

$$
\overline{\boldsymbol V}
=\left(\frac rR F(z),\frac rR Q(z),\frac1R H(z)\right).
$$

Define the basic inertial accelerations

$$
A_r=F^2+HF'-Q^2,
$$

$$
A_\theta=2FQ+HQ',
$$

$$
A_z=\frac{HH'}{R^2}.
$$

Then the radial, azimuthal, continuity, and thermal equations are

$$
F''=\chi A_r,
\qquad
Q''=\chi A_\theta,
\qquad
H'=-2F,
\qquad
T''=Pr\,H T'.
$$

Since $Q=1-G_L$, these equations are identical to the implemented form

$$
F''=\chi\left[F^2+HF'-(1-G_L)^2\right],
$$

$$
G_L''=\chi\left[2FG_L+HG_L'-2F\right].
$$

The axial momentum balance determines the basic pressure and is not needed to
compute $F,Q,H,T$.

## 4. Linearization before the normal-mode substitution

Write

$$
\boldsymbol V=\overline{\boldsymbol V}+\epsilon\boldsymbol u,
\quad T=\overline T+\epsilon\vartheta,
\quad P=\overline P+\epsilon p,
$$

with $\boldsymbol u=(u,v,w)$. Since

$$
\chi(\overline T+\epsilon\vartheta)
=\overline\chi-\epsilon\vartheta+O(\epsilon^2),
$$

the first-order momentum equation is

$$
\partial_t\boldsymbol u
+\overline\chi\left[
(\overline{\boldsymbol V}\cdot\nabla)\boldsymbol u
+(\boldsymbol u\cdot\nabla)\overline{\boldsymbol V}
\right]
-\vartheta(\overline{\boldsymbol V}\cdot\nabla)
\overline{\boldsymbol V}
=-\nabla p+\frac1R\nabla^2\boldsymbol u.
$$

The last term on the left is the Lopez thermal-inertial feedback. It is the
basic local acceleration, not a constant centrifugal source.

The other two linearized equations are

$$
\nabla\cdot\boldsymbol u=0,
$$

$$
\partial_t\vartheta
+\overline{\boldsymbol V}\cdot\nabla\vartheta
+\boldsymbol u\cdot\nabla\overline T
=\frac1{Pr\,R}\nabla^2\vartheta.
$$

## 5. Local normal modes

Use an inertial-frame normal mode

$$
\{u,v,w,\vartheta,p\}(r,\theta,z,t)
=\{\hat u,\hat v,\hat w,\hat\vartheta,\hat p\}(z)
e^{i[\alpha(r-R)+\beta R\theta-\omega_I t]}.
$$

Thus $R^{-1}\partial_\theta\mapsto i\beta$. Let

$$
D=\frac{d}{dz},\qquad k^2=\alpha^2+\beta^2,
\qquad K(z)=\alpha F+\beta Q,
\qquad \Delta_k=D^2-k^2.
$$

The local boundary-layer form of the linearized equations is

$$
\begin{aligned}
0={}&\left[-i\omega_I+i\chi K+\frac{\chi H}{R}D
+\frac{\chi F}{R}-\frac1R\Delta_k\right]u
-\frac{2\chi Q}{R}v+\chi F'w
-\frac{A_r}{R}\vartheta+i\alpha p,
\\
0={}&\frac{2\chi Q}{R}u
+\left[-i\omega_I+i\chi K+\frac{\chi H}{R}D
+\frac{\chi F}{R}-\frac1R\Delta_k\right]v
+\chi Q'w-\frac{A_\theta}{R}\vartheta+i\beta p,
\\
0={}&\left[-i\omega_I+i\chi K+\frac{\chi H}{R}D
+\frac{\chi H'}{R}-\frac1R\Delta_k\right]w
-\frac{HH'}{R^2}\vartheta+Dp,
\\
0={}&T'w+\left[-i\omega_I+iK+\frac{H}{R}D
-\frac1{Pr\,R}\Delta_k\right]\vartheta,
\\
0={}&\left(i\alpha+\frac1R\right)u+i\beta v+Dw.
\end{aligned}
$$

In terms of the stored $G_L$, the two main thermal-inertial entries are

$$
L_{u\vartheta}=-\frac1R
\left[F^2+HF'-(1-G_L)^2\right],
$$

$$
L_{v\vartheta}=\frac1R
\left[2FG_L+HG_L'-2F\right].
$$

The second sign follows from $A_\theta=-(2FG_L+HG_L'-2F)$.

The displayed equations retain the curvature terms used by the current local
stability formulation but use a scalar local viscous Laplacian. If all
$O(R^{-2})$ cylindrical viscous terms are retained, replace $\Delta_k$ in the
radial and azimuthal equations by the vector-Laplacian blocks

$$
(\Delta_k+i\alpha/R-1/R^2)u-2i\beta v/R,
$$

$$
(\Delta_k+i\alpha/R-1/R^2)v+2i\beta u/R,
$$

respectively. These terms are absent from the legacy
`archive/legacy_scripts/Stability.jl` local
boundary-layer operator and should not be added to only one model.

## 6. Block operator

For $q=(u,v,w,\vartheta,p)^T$, define

$$
\mathcal L_u=-i\omega_I+i\chi K+\frac{\chi H}{R}D
+\frac{\chi F}{R}-\frac1R\Delta_k,
$$

$$
\mathcal L_w=-i\omega_I+i\chi K+\frac{\chi H}{R}D
+\frac{\chi H'}{R}-\frac1R\Delta_k,
$$

$$
\mathcal L_T=-i\omega_I+iK+\frac{H}{R}D
-\frac1{Pr\,R}\Delta_k.
$$

Then

$$
\mathcal L q=
\begin{bmatrix}
\mathcal L_u&-2\chi Q/R&\chi F'&-A_r/R&i\alpha\\
2\chi Q/R&\mathcal L_u&\chi Q'&-A_\theta/R&i\beta\\
0&0&\mathcal L_w&-HH'/R^2&D\\
0&0&T'&\mathcal L_T&0\\
i\alpha+1/R&i\beta&D&0&0
\end{bmatrix}q=0.
$$

The signs of the three thermal-inertial entries follow directly from
$\delta\chi=-\vartheta$:

$$
\delta\{\chi(\boldsymbol V\cdot\nabla)\boldsymbol V\}
=\chi\,\delta\{(\boldsymbol V\cdot\nabla)\boldsymbol V\}
-\vartheta\,\overline{(\boldsymbol V\cdot\nabla)\boldsymbol V}.
$$

Therefore the temperature column in the three momentum rows is
$(-A_r/R,-A_\theta/R,-HH'/R^2)^T$. No temperature-feedback term occurs in
continuity, pressure, or viscous diffusion.

## 7. Temporal generalized eigenvalue problem

For real $\alpha,\beta$, separate the frequency term:

$$
\mathcal L=\mathcal L_s-i\omega_I\mathcal M,
$$

$$
\mathcal M=\operatorname{diag}(I,I,I,I,0).
$$

Therefore

$$
\boxed{\mathcal L_s q=i\omega_I\mathcal M q.}
$$

Equivalently, in the `A*q = omega*B*q` convention,

$$
A=\mathcal L_s,\qquad B=i\mathcal M.
$$

For an inertial-frame eigenvalue, $A$ is obtained by setting $\omega_I=0$ in
$L_0+\alpha L_1+\alpha^2L_2$. For a disk-fixed eigenvalue, the same construction
sets $\omega_R=0$, which means $\omega_I=\beta$ inside $A$. In both cases

$$
B=i\,\operatorname{diag}(I,I,I,I,0),
$$

and the pressure/continuity block remains a descriptor constraint.

With the chosen factor $e^{-i\omega t}$, $\operatorname{Im}(\omega)>0$ is
temporally unstable. The singular pressure block of $\mathcal M$ produces
constraint/infinite eigenvalues that must be filtered.

The identity velocity blocks in `Ta` are correct for Lopez (2013). Replacing
them by $\chi I$ would implement a different generalized model in which the
local time acceleration is density weighted.

## 8. Spatial quadratic eigenvalue problem

For prescribed real $\omega_I,\beta$, collect powers of the complex radial
wavenumber:

$$
\boxed{(L_0+\alpha L_1+\alpha^2L_2)q=0.}
$$

Let $I$ be the collocation identity and let diagonal functions denote diagonal
matrices. Define the alpha-independent operators

$$
M_0=-i\omega_I I+i\beta\,\operatorname{diag}(\chi Q)
+\operatorname{diag}\!\left(\frac{\chi H}{R}\right)D
-\frac1R(D^2-\beta^2I),
$$

$$
E_0=-i\omega_I I+i\beta\,\operatorname{diag}(Q)
+\operatorname{diag}\!\left(\frac{H}{R}\right)D
-\frac1{PrR}(D^2-\beta^2I).
$$

The complete zero-order matrix is

$$
L_0=
\begin{bmatrix}
M_0+\chi F/R&-2\chi Q/R&\chi F'&-A_r/R&0\\
2\chi Q/R&M_0+\chi F/R&\chi Q'&-A_\theta/R&i\beta I\\
0&0&M_0+\chi H'/R&-HH'/R^2&D\\
0&0&T'&E_0&0\\
I/R&i\beta I&D&0&0
\end{bmatrix}.
$$

Every product such as $\chi F$ in this expression denotes a diagonal matrix.
This displayed $L_0$ is exactly the block layout assembled by
`lopez_spatial_matrices`.

In the code-compatible local viscous approximation,

$$
L_2=\operatorname{diag}
\left(\frac1R,\frac1R,\frac1R,\frac1{PrR},0\right),
$$

and

$$
L_1=
\begin{bmatrix}
i\chi F&0&0&0&i\\
0&i\chi F&0&0&0\\
0&0&i\chi F&0&0\\
0&0&0&iF&0\\
i&0&0&0&0
\end{bmatrix}.
$$

Thus the only pressure entries are $L_{1,up}=iI$, $L_{0,vp}=i\beta I$, and
$L_{0,wp}=D$. The only continuity entries are $L_{1,cu}=iI$,
$L_{0,cu}=I/R$, $L_{0,cv}=i\beta I$, and $L_{0,cw}=D$.

$L_0$ contains the $D,D^2,\beta,\omega_I$ terms displayed in section 5.
The standard companion linearization, with $s=\alpha q$, is

$$
\begin{bmatrix}
0&I\\-L_0&-L_1
\end{bmatrix}
\begin{bmatrix}q\\s\end{bmatrix}
=\alpha
\begin{bmatrix}
I&0\\0&L_2
\end{bmatrix}
\begin{bmatrix}q\\s\end{bmatrix}.
$$

Because the pressure/continuity part makes $L_2$ singular, this descriptor
problem also contains infinite eigenvalues. Physical spatial modes must satisfy
the original quadratic residual after reconstruction.

## 9. Boundary conditions

For a fixed-temperature no-slip disk and a decaying far field,

$$
u=v=w=\vartheta=0\quad\text{at }z=0,
$$

$$
u,v,w,\vartheta\rightarrow0\quad\text{as }z\rightarrow\infty.
$$

Pressure disturbance decays in the far field, so

$$
p\rightarrow0\quad\text{as }z\rightarrow\infty.
$$

No pressure value is prescribed at the wall: $p(0)$ is determined by the
axial momentum equation and continuity. A separate pressure gauge is needed
only in a zero-horizontal-wavenumber formulation that does not already impose
far-field pressure decay.

In the present collocation implementation, eliminate the eight endpoint
degrees of freedom for $u,v,w,\vartheta$ and the far-field pressure degree of
freedom. The wall-pressure degree of freedom is retained. An equivalent
tau-row replacement is also valid.

## 10. Relation to the disk-fixed frequency

The Lopez closure is naturally defined in the inertial frame. If the normal
mode is reported in the disk-fixed frame, the real frequency is shifted by the
azimuthal frame rate. In the present local scaling the disk angular speed is
$1/R$, the integer azimuthal mode is $n=\beta R$, and
$\theta_D=\theta_I-t/R$. Hence

$$
e^{i(n\theta_D-\omega_Dt)}
=e^{i[n\theta_I-(\omega_D+n/R)t]},
$$

so that

$$
\omega_I=\omega_R+\beta,
$$

where the subscript $R$ denotes the disk-fixed value used by the calling
notebook. The imaginary part, and therefore the growth rate and neutral
condition, is unchanged. For the temperature equation, whose advection is not
density weighted,

$$
-i\omega_I+i\beta Q
=-i\omega_R+i\beta(Q-1),
$$

which is exactly the classical disk-fixed transport written with the relative
azimuthal profile $G_R=Q-1=-G_L$. The momentum transport is different because
its advection is multiplied by $\chi$:

$$
-i\omega_I+i\beta\chi Q
=-i\omega_R+i\beta\chi G_R+i\beta(\chi-1).
$$

The final term is required by the Lopez model whenever $T\ne1$. In the code it
is included automatically by retaining the inertial profile $Q$ and applying
the unweighted frequency shift $\omega_I=\omega_R+\beta$. Multiplying a
classical disk-fixed advection operator by $\chi$ would omit this term.

One must not obtain a Lopez operator by simply multiplying all terms of a
rotating-frame traditional-Boussinesq matrix by $\chi$. The 2013 closure
weights inertial-frame advection only and is not term-by-term frame invariant.

## 11. Required changes relative to the legacy `archive/legacy_scripts/Stability.jl`

The present `Spatial_mode_BEK` is a traditional fixed-reference centrifugal
model. A Lopez builder needs the following changes.

1. Convert the base azimuthal profile to $Q=1-G_L$ and assemble the operator in
   the inertial frame. Do not feed $G_L$ directly into the old `G` blocks.
2. Multiply every velocity advection, curvature, and base-shear coefficient by
   $\chi=2-T$: the momentum parts of the `A`, `B`, `C`, and `D1` blocks all
   change. Pressure, diffusion, continuity, and the temperature equation do
   not receive this multiplier.
3. Replace the constant traditional radial feedback `D1[1,4] = -1/R` by
   $-A_r/R$.
4. Add the missing azimuthal feedback $-A_\theta/R$ and axial feedback
   $-HH'/R^2$.
5. Keep `Ta=diag(I,I,I,I,0)` for the Lopez (2013) model.
6. Use $Pr$ consistently in the thermal diffusion block.
7. Impose $p(\infty)=0$, but do not prescribe $p(0)$.

The compressible matrices in `CRD_STA.jl` are useful as an organizational
example (five variables, singular mass matrix, temporal and spatial operator
separation), but their density equation, equation of state, variable-property
diffusion, and pressure-temperature coupling must not be copied into this
incompressible Lopez system.

## 12. Verification sequence

Before tracing a neutral curve, verify the implementation in this order.

1. Check the basic-flow residuals $F''-\chi A_r$,
   $Q''-\chi A_\theta$, $H'+2F$, and $T''-PrHT'$ on the stability grid.
2. At $T_w=1$, recover the classical hydrodynamic benchmark. Temperature is
   one-way coupled when $T'=0$, so retaining thermal feedback must not move the
   hydrodynamic eigenvalues.
3. Compare every analytic matrix column against a finite-difference
   Jacobian-vector product of the nonlinear Lopez residual.
4. Verify temporal eigenvalues with increasing $N$ and $z_{max}$; discard
   infinite/constraint modes by residual and boundary-condition tests.
5. For a spatial mode, check
   $\|(L_0+\alpha L_1+\alpha^2L_2)q\|/\|q\|$ after companion reconstruction.
6. Only after these checks compare neutral curves with the traditional and
   compressible models.

The executable audit `test_lopez_stability.jl` implements these checks without
reusing the matrix block formulas. It evaluates the nonlinear cylindrical
Lopez residual at positive and negative perturbation amplitudes and constructs
every Jacobian column by centered finite differences. Evaluating that Jacobian
at $\alpha=0,+1,-1$ independently recovers $L_0,L_1,L_2$. The test also checks
the temporal/spatial operator identity, the disk/inertial frequency shift, the
isothermal classical hydrodynamic submatrix, and all boundary-condition indices.

## 13. Isothermal two-mode validation against Malik (1986)

Malik's two stationary minima are used as separate benchmarks:

| Mode | Mechanism | $R$ | $\beta$ | $\omega_R$ | $\alpha_{r,\mathrm{ref}}$ |
|---|---|---:|---:|---:|---:|
| Type I | upper-branch cross-flow | 285.36 | 0.07759 | 0 | 0.38482 |
| Type II | lower-branch viscous curvature/Coriolis | 440.88 | 0.04672 | 0 | 0.13228 |

At $T_w=1$, $T=1$, $\chi=1$, and $T'=0$. The Lopez system therefore becomes
the classical incompressible rotating-disk hydrodynamic operator. In primitive
variables its four hydrodynamic equations reduce term by term to Malik's
(2.16)-(2.19) after, as in Malik's final sixth-order system, dropping terms of
order $R^{-2}$ and smaller. The thermal equation is one-way coupled and cannot
shift either hydrodynamic mode.

Using the mapped $z_{max}=20$ grid, $p(\infty)=0$, and the continuous BVP basic
flow gives

| $N$ | Type-I $\alpha_r$ | Type-I $\alpha_i$ | Type-II $\alpha_r$ | Type-II $\alpha_i$ |
|---:|---:|---:|---:|---:|
| 39 | 0.385052489904 | $3.48665925\times10^{-4}$ | 0.131849539510 | $9.07981525\times10^{-4}$ |
| 49 | 0.385052489830 | $3.48665830\times10^{-4}$ | 0.131849549642 | $9.07976466\times10^{-4}$ |
| 69 | 0.385052489834 | $3.48665836\times10^{-4}$ | 0.131849548389 | $9.07976809\times10^{-4}$ |
| 99 | 0.385052489835 | $3.48665831\times10^{-4}$ | 0.131849548368 | $9.07976804\times10^{-4}$ |

The Type-I and Type-II real-wavenumber errors are respectively $0.060\%$ and
$0.325\%$. Original quadratic-pencil residuals are below $10^{-13}$. Changing
the finite far field from $z_{max}=20$ to 40 leaves both eigenvalues unchanged
to at least nine significant digits. Thus the remaining non-zero $\alpha_i$
is a converged systematic difference between the present primitive-variable
collocation/far-field treatment and Malik's eliminated-pressure asymptotic
conditions evaluated at the rounded published minima; it is not an eigensolver
or domain-resolution failure.

For automated regression, `validate_lopez_benchmark.jl` requires both modes to
satisfy a real-wavenumber error below $0.5\%$, $|\alpha_i|<10^{-3}$, and an
original-pencil residual below $10^{-10}$. Both modes pass. With
$\exp[i(\alpha r+\beta R\theta-\omega t)]$, positive $\alpha_i$ denotes weak
downstream spatial decay at the rounded benchmark point.

## 14. Neutral-curve tracker interface

`LopezStability.jl` provides an `eigsol`-compatible targeted solver:

```julia
include("LopezStability.jl")
using .LopezStability

eigval, eigvec = LopezStability.eigsol_lopez(
    F, G, H, T, R, omega, beta, N_cheb, D, D2, c, 1,
)
```

Here `G` is the Lopez deficit $G_L$, with $G_L(0)=0$ and
$G_L(\infty)=1$. The default return layout matches the old `eigsol`: `eigval`
is ordered by distance from the shift `c`, and each column of `eigvec` contains
the reduced collocation degrees of freedom after applying the nine boundary
conditions. For continuation, set `c=eigval[1]` at the next point and request
one eigenvalue.

To obtain fields with their boundary zeros reinserted, use

```julia
eigval, eigvec, info = LopezStability.eigsol_lopez(
    F, G, H, T, R, omega, beta, N_cheb, D, D2, c, 1;
    full_eigenvectors=true,
    return_info=true,
)

n = N_cheb + 1
u     = eigvec[1:n, 1]
v     = eigvec[n+1:2n, 1]
w     = eigvec[2n+1:3n, 1]
theta = eigvec[3n+1:4n, 1]
p     = eigvec[4n+1:5n, 1]
```

`info.residuals` contains the residuals of the original quadratic pencil, not
only the internal IAR convergence estimate. `info.partial=true` indicates that
fewer than the requested number of modes converged; the converged target mode
is still returned.

## References

1. J. M. Lopez, F. Marques & M. Avila (2013), "The Boussinesq approximation
   in rapidly rotating flows", JFM 737, 56-77,
   https://doi.org/10.1017/jfm.2013.558.
2. H. M. Blackburn et al. (2021), "On the Boussinesq approximation in
   arbitrarily accelerating frames of reference", JFM 924, R1,
   https://doi.org/10.1017/jfm.2021.640.
3. M. R. Malik (1986), "The neutral curve for stationary disturbances in
   rotating-disk flow", JFM 164, 275-287,
   https://doi.org/10.1017/S0022112086002550.
