# Stand-alone draft responses to Referee 3, Comments 6 and 7

> Draft status: the additional calculations reported below have been completed,
> but this text has not been inserted into the manuscript or the consolidated
> response letter. Manuscript changes are therefore described as planned
> revisions.

## Comment 6

> That the peak of the amplitude (which probably should be referred to as
> coupling or receptivity coefficient) appears close to the location of maximum
> growth rate (the last paragraph of the left column, page 17) was taken as a
> main result by the authors. There seems no mathematical reason for this. Is
> this purely a coincidence? This peak is unlikely to be significant because
> the mode excited has missed the exponential amplification upstream. If the
> same form of wall forcing is placed near the lower-branch neutral curve, the
> mode excited would probably acquire a much greater amplitude when it reaches
> the location of the maximum growth.
>
> With the above consideration, I would like to suggest that the authors
> calculate the initial amplitude at the lower-branch neutral curve, and follow
> the development of the mode to the critical location for absolute
> instability. The amplitude gained would give some useful idea of how likely
> the convective or absolute instability causes transition. The authors made a
> good point in the introduction on the relative role of convective and
> absolute instabilities. They should demonstrate it using concrete result.

### Response

We thank the Reviewer for this important observation. We agree that $C_r$ is
a local coupling (or receptivity) coefficient, rather than the downstream
amplitude of a single disturbance, and that there is no general mathematical
argument requiring its maximum to coincide with the maximum spatial growth
rate. Our original discussion therefore attached too much physical
significance to this local correlation.

Following the Reviewer's suggestion, we performed an additional calculation
for the representative stationary disturbance

$$
a_s=0,\qquad Re_h=1000,\qquad \bar{\omega}=0,\qquad n=30,
$$

using the same roughness profile and width as in the manuscript
($c^2=1$, $l_s=0.5$, and $h_r=1/Re_h$). For forcing centered at $R_f$, the
downstream amplification and total amplitude were evaluated as

$$
N(R;R_f)=-\int_{R_f}^{R}\alpha_i(\xi)\,\mathrm{d}\xi,
\qquad
\lvert A(R;R_f)\rvert
=\lvert C_r(R_f)\rvert\exp[N(R;R_f)].
$$

The lower-branch neutral point is $R_l=287.822$, the maximum-growth location
is $R_g=372.989$, and the independently determined onset radius of absolute
instability used in the manuscript is $R_{\mathrm{abs}}=570$. A disturbance introduced at $R_l$
accumulates $N(R_g;R_l)=1.7914$ before reaching $R_g$, giving
$\lvert A(R_g;R_l)\rvert=0.8613$. By comparison, introducing the same forcing
locally at $R_g$ gives $\lvert C_r(R_g)\rvert=0.19396$. The upstream
disturbance is therefore already approximately $4.44$ times larger when it
reaches the maximum-growth location.

Propagating the disturbances to $R_{\mathrm{abs}}=570$ gives the following
comparison:

| Forcing location | $R_f$ | $\lvert C_r(R_f)\rvert$ | $N(570;R_f)$ | $\lvert A(570;R_f)\rvert$ |
|---|---:|---:|---:|---:|
| Lower-branch neutral point | 287.822 | 0.143595 | 4.33600 | 10.9709 |
| Maximum-growth location | 372.989 | 0.193957 | 2.54460 | 2.47064 |
| Peak of the local $C_r$ curve | 468.910 | 0.244415 | 0.426114 | 0.374271 |

Although the local coupling coefficient is smallest at the lower neutral
point, its longer amplification distance gives the largest downstream
amplitude. At $R_{\mathrm{abs}}$, this amplitude is approximately $4.44$
times that obtained by forcing at the maximum-growth location and $29.3$
times that obtained by forcing at the peak of the local $C_r$ curve. The
corresponding convective amplification from $R_l$ to $R_{\mathrm{abs}}$ is
$\exp(4.336)\simeq76.4$. As an additional upstream check, forcing at $R=360$
gives $\lvert A(570)\rvert=3.61342$. The recomputed $C_r$ curve is still
increasing at $R=360$, so we do not identify this point as a strict local
maximum.

These results support the Reviewer's central point: local receptivity and
downstream amplified amplitude must be distinguished, and the local $C_r$
peak is not, by itself, the most dangerous forcing location. In the linear
estimate, the disturbance introduced at $R_l$ already reaches order-one
amplitude at approximately $R=378$, well upstream of $R_{\mathrm{abs}}$.
This result suggests that convective amplification can become dynamically
important before the onset of absolute instability. However, because this is
a locally parallel linear $e^N$ estimate and no calibrated transition
threshold or nonlinear saturation is included, the order-one value also marks
the limit beyond which the extrapolated amplitude should not be interpreted
quantitatively. The calculation therefore provides a concrete comparison but
does not, by itself, determine the actual transition radius or prove which
mechanism ultimately triggers transition.

In the revision, we plan to define $N(R;R_f)$ and $A(R;R_f)$ explicitly,
report the comparison above, refer to $C_r$ consistently as a local coupling
or receptivity coefficient, and weaken the earlier claim of a general
correspondence between the $C_r$ and growth-rate maxima. Exact page and line
numbers will be supplied after these changes are incorporated.

## Comment 7

> There appears a mathematical reason why $C_r$ attains its maximum in the
> region of "modal switching" (the last paragraph of the left column, page 17).
> The authors need to explain this further by stating clearly which two modes
> are switched over. I do not think that the modes would have exactly the same
> eigenvalues at the switching point, but probably their real parts become
> (nearly) identical.

### Response

We thank the Reviewer for requesting this clarification. The two modes
intended in our original discussion were the Type-I (crossflow/inviscid) and
Type-II (viscous) spatial modes. To test the proposed spectral interpretation,
we separately resolved the two spatial branches that had been associated with
these modal families in the original discussion over $420\leq R\leq500$.
Near $R=470$, where their real parts are closest on the computed radial grid,
the eigenvalues are

$$
\alpha_{\mathrm{I}}
=0.177999-0.007407\,\mathrm{i},
\qquad
\alpha_{\mathrm{II}}
=0.177519+0.024768\,\mathrm{i}.
$$

Consequently,

$$
\lvert\Delta\alpha_r\rvert=4.80\times10^{-4},
\qquad
\lvert\Delta\alpha_i\rvert=3.22\times10^{-2},
\qquad
\lvert\alpha_{\mathrm{I}}-\alpha_{\mathrm{II}}\rvert
=3.22\times10^{-2}.
$$

We therefore agree with the Reviewer that the real parts become nearly
identical, whereas the complete complex eigenvalues remain distinct; the
calculation provides no evidence of an exact eigenvalue coalescence.

The separated-branch calculation also shows that our original phrase "modal
switching" was too strong. Over the interval examined, the receptivity
coefficient of branch 1 remains larger than that of branch 2, including at
$R=470$, where

| Spatial branch | $\lvert C_r\rvert$ | $\lvert Q\rvert$ | $\lvert\mathrm{B.C.}\rvert$ |
|---|---:|---:|---:|
| Branch 1 (originally associated with Type I) | 0.244259 | $3.57753\times10^{-4}$ | $8.73844\times10^{-5}$ |
| Branch 2 (originally associated with Type II) | 0.172435 | $3.20296\times10^{-4}$ | $5.52302\times10^{-5}$ |

The directly demonstrated origin of the Type-I $C_r$ maximum follows from

$$
\lvert C_{r,j}\rvert
=\frac{\lvert\mathrm{B.C.}_j\rvert}{\lvert Q_j\rvert},
\qquad j=1,2.
$$

Along branch 1 from $R=420$ to $R=470$,
$\lvert Q\rvert$ decreases by approximately $26\%$, whereas
$\lvert\mathrm{B.C.}\rvert$ decreases by only approximately $12\%$;
consequently, $\lvert C_r\rvert$ increases by approximately $20\%$. The
maximum is therefore associated primarily with the reduction of the
normalization denominator. Its occurrence in the same radial range as the
near equality of the two eigenvalue real parts is a coincident spectral
feature of this case, but it does not establish either an exact crossing of
the complex eigenvalues or a change of the dominant receptivity branch.
In particular, the present calculation does not demonstrate that proximity
of the second branch causes the minimum of $\lvert Q\rvert$; that proximity
is therefore reported as a correlation rather than a causal mechanism.

In the revision, we plan to replace "modal switching" with the more precise
description "near real-part crossing of two distinct spatial branches
associated with the Type-I and Type-II modal families" and to report both
complex eigenvalues together with the
$\lvert\mathrm{B.C.}\rvert/\lvert Q\rvert$ decomposition. Exact page and line
numbers will be supplied after these changes are incorporated.
