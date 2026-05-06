# Notes on HLL / HLLC Treatment of NC Terms

This note describes the numerical constructions themselves and intentionally avoids implementation-specific option names.

## 1. Two equivalent forms

For the volume-fraction coupling, the product rule gives two equivalent forms:

$$
\partial_t \alpha + u\,\partial_x \alpha = 0,
$$

$$
\partial_t \alpha + \partial_x(\alpha u) = \alpha\,\partial_x u.
$$

With the Kapila correction, these become

$$
\partial_t \alpha + u\,\partial_x \alpha = \pm K\,\partial_x u,
$$

$$
\partial_t \alpha + \partial_x(\alpha u) = (\alpha \pm K)\,\partial_x u.
$$

Numerically, these lead to two methods:

Method 1, or alpha-interface:

$$
\text{flux side: } (U,F)=(\alpha,0),
\qquad
\text{NC side: } (U,F)=(1,\alpha).
$$

This is not a separate physical auxiliary PDE. It is simply the HLL evaluation obtained when the left-hand transport part is written with zero flux, so that the same HLL diffusion and regularization are retained.

Method 2, or u-interface:

$$
\text{flux side: } (U,F)=(\alpha,\alpha u),
\qquad
\text{NC side: } (U,F)=(1,u).
$$

For tangential hypo terms, the same Method 2 idea is used with

$$
(U,F)=(1,u_t).
$$

## 2. Which method is used

For the volume-fraction term:

- HLL uses both Method 1 and Method 2.
- HLLC uses Method 2 only.

For all other NC terms:

- The $K\,\partial_x u$ term always uses Method 2.
- The hypoelastic terms involving $\partial_x u_n$ always use Method 2.
- The hypoelastic terms involving $\partial_x u_t$ always use Method 2.

For example, in the 1D $x$-directed hypo subsystem one has terms such as

$$
\partial_t(\rho\tau_{xx})+\partial_x(\rho u\tau_{xx})
=
\rho\left(\frac{4G}{3}+\tau_{xx}\right)\partial_x u,
$$

$$
\partial_t(\rho\tau_{xy})+\partial_x(\rho u\tau_{xy})
=
\rho\left(G+\tau_{xx}\right)\partial_x v.
$$

So these terms explicitly require consistent approximations of both the normal velocity gradient $\partial_x u$ and the tangential velocity gradient $\partial_x v$.

So even when HLL uses alpha-interface for the pure $\alpha$ transport part, the $K\,\partial_x u$ and hypo terms are still built from interface-consistent velocity traces.

## 3. HLL

With wave-speed bounds $S_L<S_R$, the HLL flux is

$$
F_{\mathrm{HLL}}(U,F)=
\begin{cases}
F_L, & 0\le S_L,\\[4pt]
\dfrac{S_R F_L-S_L F_R+S_LS_R(U_R-U_L)}{S_R-S_L}, & S_L\le 0\le S_R,\\[10pt]
F_R, & S_R\le 0.
\end{cases}
$$

### HLL Method 1: alpha-interface

For the flux side,

$$
F_{\mathrm{HLL}}^{(\alpha,\;0)}=
F_{\mathrm{HLL}}(U=\alpha,F=0)=
\begin{cases}
0, & 0\le S_L,\\[4pt]
\dfrac{S_LS_R(\alpha_R-\alpha_L)}{S_R-S_L}, & S_L\le 0\le S_R,\\[10pt]
0, & S_R\le 0.
\end{cases}
$$

For the interface-consistent $\alpha$ trace on the NC side,

$$
\Psi_{\alpha,\mathrm{HLL}}=
F_{\mathrm{HLL}}(U=1,F=\alpha)=
\begin{cases}
\alpha_L, & 0\le S_L,\\[4pt]
\dfrac{S_R\alpha_L-S_L\alpha_R}{S_R-S_L}, & S_L\le 0\le S_R,\\[10pt]
\alpha_R, & S_R\le 0.
\end{cases}
$$

### HLL Method 2: u-interface

For the flux side,

$$
F_{\mathrm{HLL}}^{(\alpha,\;\alpha u)}=
F_{\mathrm{HLL}}(U=\alpha,F=\alpha u)=
\begin{cases}
\alpha_L u_L, & 0\le S_L,\\[4pt]
\dfrac{S_R\alpha_Lu_L-S_L\alpha_Ru_R+S_LS_R(\alpha_R-\alpha_L)}{S_R-S_L},
& S_L\le 0\le S_R,\\[10pt]
\alpha_R u_R, & S_R\le 0.
\end{cases}
$$

For the normal velocity trace on the NC side,

$$
\Psi_{u,\mathrm{HLL}}=
F_{\mathrm{HLL}}(U=1,F=u)=
\begin{cases}
u_L, & 0\le S_L,\\[4pt]
\dfrac{S_Ru_L-S_Lu_R}{S_R-S_L}, & S_L\le 0\le S_R,\\[10pt]
u_R, & S_R\le 0.
\end{cases}
$$

For tangential hypo terms,

$$
\Psi_{u_t,\mathrm{HLL}}=
F_{\mathrm{HLL}}(U=1,F=u_t)=
\begin{cases}
u_{t,L}, & 0\le S_L,\\[4pt]
\dfrac{S_Ru_{t,L}-S_Lu_{t,R}}{S_R-S_L}, & S_L\le 0\le S_R,\\[10pt]
u_{t,R}, & S_R\le 0.
\end{cases}
$$

## 4. HLLC

HLLC uses Method 2 only. Let

$$
S_L<S_M<S_R,
$$

and define

$$
\zeta_K=\frac{S_K-u_K}{S_K-S_M}
=\frac{\rho_K^*}{\rho_K},
\qquad K\in\{L,R\}.
$$

In this construction, HLLC is used only in Method 2.

### HLLC Method 2: u-interface

For the flux side,

$$
F_{\mathrm{HLLC}}^{(\alpha,\;\alpha u)}=
\begin{cases}
\alpha_L u_L, & 0\le S_L,\\[4pt]
\alpha_L S_M\zeta_L, & S_L\le 0\le S_M,\\[4pt]
\alpha_R S_M\zeta_R, & S_M\le 0\le S_R,\\[4pt]
\alpha_R u_R, & S_R\le 0.
\end{cases}
$$

For the normal velocity trace on the NC side,

$$
\Psi_{u,\mathrm{HLLC}}=
F_{\mathrm{HLLC}}(U=1,F=u)=
\begin{cases}
u_L, & 0\le S_L,\\[4pt]
S_M\zeta_L, & S_L\le 0\le S_M,\\[4pt]
S_M\zeta_R, & S_M\le 0\le S_R,\\[4pt]
u_R, & S_R\le 0.
\end{cases}
$$

For tangential hypo terms, the tangential trace comes from the HLLC tangential star state:

$$
\Psi_{u_t,\mathrm{HLLC}}=
\begin{cases}
u_{t,L}, & 0\le S_L,\\[4pt]
u_t^*, & S_L\le 0\le S_R,\\[4pt]
u_{t,R}, & S_R\le 0.
\end{cases}
$$

### Why the star-branch trace is $S_M\zeta_L$

Take any co-moving scalar $z$ on the left star branch. The HLLC jump condition gives

$$
S_L(z_L^*-z_L)=z_L^*S_M-z_Lu_L.
$$

Hence

$$
z_L^*=z_L\,\frac{S_L-u_L}{S_L-S_M}
=z_L\zeta_L.
$$

The left star-branch HLLC flux is then

$$
F_{\mathrm{HLLC},L}(z)
=z_Lu_L+S_L(z_L^*-z_L).
$$

Substituting $z_L^*=z_L\zeta_L$ gives

$$
F_{\mathrm{HLLC},L}(z)
=z_L\left[u_L+S_L(\zeta_L-1)\right]
=z_L S_M\zeta_L.
$$

Therefore, for the unit-variable problem $z=1$,

$$
\Psi_{u,\mathrm{HLLC},L}=S_M\zeta_L.
$$

The right branch is identical:

$$
\Psi_{u,\mathrm{HLLC},R}=S_M\zeta_R.
$$

This is the key nonstandard point in HLLC: the transport trace entering the NC terms is not $S_M$, but $S_M\zeta_K$.

## 5. ADC remark

With ADC, the non-conservative transport traces should be blended between HLL Method 2 and HLLC:

$$
\Psi^{\mathrm{ADC}}
=
\Psi^{\mathrm{HLL,\,M2}}
+\phi\left(\Psi^{\mathrm{HLLC}}-\Psi^{\mathrm{HLL,\,M2}}\right),
$$

applied to both the normal and tangential traces, with the conservative transport flux blended by the same sensor. This is the mathematically consistent pairing because both endpoints are transport-consistent velocity-trace constructions for gradient-driven NC terms, whereas Method 1 is based on an $\alpha$-interface construction and is not the correct low-order partner for terms involving velocity gradients.
