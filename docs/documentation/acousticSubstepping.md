@page acousticSubstepping Split-Explicit Acoustic Substepping (Low-Mach)

# Split-Explicit Acoustic Substepping

This page describes the `acoustic_substepping` time integrator, which relaxes the
acoustic-CFL time-step restriction at low Mach number without introducing any global
(elliptic / pressure-Poisson) solve. For the parameters that activate it, see
@ref case "Case Files".

## Motivation

An explicit compressible solver is limited by the acoustic CFL condition,

\f[ \Delta t \;\lesssim\; \frac{\Delta x}{|u| + c}, \f]

where \f$c\f$ is the sound speed. When \f$c \gg |u|\f$ (low Mach number \f$M = |u|/c \ll 1\f$),
this step is a factor of \f$\sim 1/M\f$ smaller than the time scale on which the flow
itself evolves, \f$\Delta x/|u|\f$. The solver then spends almost all of its work
resolving fast acoustic waves that carry little of the dynamics.

Pressure-projection and implicit-acoustic methods remove this restriction but require a
global elliptic solve for the pressure each step. The associated global reductions and
sparse linear solves scale poorly on GPU clusters. Acoustic substepping keeps the solver
fully explicit, with the same nearest-neighbor halo-exchange communication pattern as the
standard time stepper, and is well suited to GPU/exascale hardware.

The scheme follows the split-explicit approach of atmospheric dynamical cores
(Klemp, Wicker & Skamarock; Wicker & Skamarock, *Mon. Wea. Rev.* 2002), adapted to MFC's
Godunov finite-volume framework after Nazari & Nair, *J. Adv. Model. Earth Syst.* 2017.

## Method

Each SSP-RK3 stage splits the governing equations into a **slow** (advective) part,
advanced once at the large advective step \f$\Delta t \sim \Delta x / |u|\f$, and a
**fast** (acoustic) part subcycled on \f$n_s \approx (|u|+c)/|u|\f$ micro-steps of size
\f$\Delta\tau = \Delta t / n_s\f$:

| Term | Mode | Discretization |
|------|------|----------------|
| Momentum advection \f$\nabla\!\cdot(\rho \mathbf{u}\otimes\mathbf{u})\f$ | slow | WENO + Riemann (contact-speed dissipation) |
| Volume-fraction and viscous/source terms | slow | existing schemes |
| Mass transport \f$\nabla\!\cdot(\alpha_k\rho_k \mathbf{u})\f$ | fast | 2nd-order centered, subcycled |
| Energy transport \f$\nabla\!\cdot((\rho E + p)\mathbf{u})\f$ | fast | 2nd-order centered, subcycled |
| Pressure gradient \f$\nabla p\f$ (momentum) | fast | 2nd-order centered, subcycled |

The expensive WENO + Riemann flux is evaluated once per RK stage and held frozen during
the subcycle, while each acoustic micro-step is a low-order stencil. Because the costly
work runs at the advective rate rather than the acoustic rate, the wall-clock cost drops
by roughly \f$O(1/M)\f$ relative to the standard explicit scheme.

### Slow flux

The slow flux reuses the HLLC Riemann solver with two modifications, both active only
when `acoustic_substepping` is set: the pressure term is removed from the momentum and
energy flux (it is handled by the subcycle), and the numerical dissipation is capped to
the contact wave speed \f$s_\star\f$. The latter is essential: the standard HLLC
dissipation scales with \f$|u| \pm c\f$, which is unstable at the advective time step;
restricting it to \f$s_\star\f$ makes the convective flux \f$|u|\f$-stable. Low-Mach
accuracy of the convective flux is handled by the existing `low_Mach` correction.

### Acoustic substep

The acoustic micro-step (`s_acoustic_substep` in `src/simulation/m_acoustic_substep.fpp`)
performs a forward–backward update: mass and energy transport are advanced first, the
pressure is recomputed from the equation of state, and the momentum is then advanced with
the new pressure gradient. The frozen slow forcing is added as a \f$\Delta\tau\f$-scaled
source at each micro-step. To suppress the acoustic noise that a centered scheme would
otherwise accumulate, a grad–div divergence damping term is applied; its discrete operator
is rank-one and therefore annihilates discretely divergence-free (vortical) modes, damping
only the compressive/acoustic content. The forward sweep computes flux divergences from a
frozen snapshot of the conserved state and applies them afterward, which keeps the scheme
species-mass conservative.

### Two-tier flux at discontinuities

A purely centered acoustic substep rings at a sharp jump. To handle embedded
discontinuities the substep runs **two tiers**, selected per face each RK stage by a
discontinuity criterion:

- **Smooth tier (default).** The cheap 2nd-order centered transport above, used wherever
  the flow is smooth. This is the low-Mach-accurate, low-cost path and is bit-identical to
  the pre-shock-capture scheme on a smooth field.
- **Robust tier (flagged faces).** At a flagged face the centered mass and energy
  convective fluxes are replaced by the **full HLLC** flux (the convective Riemann flux
  `s_convective_face_flux` plus the acoustic star-pressure flux `s_acoustic_face_flux`)
  built on one consistently WENO-reconstructed left/right state. The per-component delta is
  applied through the same single-valued telescoping structure, so species mass and total
  energy remain exactly conservative. Total energy at the face is rebuilt from the
  reconstructed primitives via the stiffened-gas mixture EOS.

The robust tier replaces the **mass, energy, and momentum** convective transport with full
HLLC. Mass and energy are made live directly from the reconstructed left/right state. The
**convective momentum** flux is also made live, but through a live-minus-frozen delta: the
slow (HLLC) path already computes an exact per-face convective-momentum flux at the RK stage
state, which `s_compute_rhs` stores; the substep subtracts that stored flux and adds the
live `s_convective_face_flux` value. Subtracting the *stored* stage-entry flux (rather than a
reconstruction of it) makes the slow/frozen cancellation exact on every RK stage, which is
what keeps the scheme stable when convective momentum transport is leading order — i.e. at a
**propagating shock**. The transverse (off-normal) momentum slots of the rotated HLLC flux
are mapped back into the stored global-component flux through the same single-valued
telescoping delta, so transverse momentum is transported consistently in multiple
dimensions (validated in 2D, below). Only the acoustic star-pressure momentum, which the
slow path omits under `acoustic_substepping`, is added on top.

With this, **propagating shocks and shock–interface interactions are captured** in split
mode, not merely material contacts. The per-component delta is applied through a
single-valued telescoping face flux, so species mass and total energy remain exactly
conservative (drift \f$\sim 10^{-16}\f$). Total energy at the face is rebuilt from the
reconstructed primitives via the stiffened-gas mixture EOS.

The criterion (`s_acoustic_flag`) flags a face when either
- the largest per-cell pressure jump over the 4-point stencil exceeds \f$\sim 1\%\f$ of the
  local pressure scale (a scale-invariant, `weno_eps`-independent max-indicator test), or
- the volume-fraction variation \f$|\Delta\alpha|\f$ across the face exceeds a small
  threshold (the material-interface flag).

The flag mask is computed once from the stage-entry state and held fixed for that stage's
microsteps; the reconstruction is **conditional** (skipped entirely when no face on the
rank is flagged), so the smooth-flow cost — and the \f$O(1/M)\f$ speedup — is preserved.
Because both flux entries hardcode the direct wave-speed estimate, the robust tier requires
``wave_speeds = 'direct'`` (enforced at input checking).

### Time step

`s_compute_dt` sets \f$\Delta t\f$ from the advective CFL (using \f$|u|\f$, not
\f$|u|+c\f$) and computes \f$n_s\f$ from a domain-maximum reduction of
\f$(|u|+c)/|u|\f$. The only global collective is that existing time-step reduction.

### Temporal accuracy

The Wicker–Skamarock RK3 macro coupling is **second-order accurate in the advective
step** \f$\Delta t\f$ when the acoustic micro-step \f$\Delta\tau\f$ is held fixed (i.e.
\f$n_s\f$ scaled with \f$\Delta t\f$ so the acoustic CFL is constant). A fixed-grid
\f$\Delta t\f$-refinement of a smooth translating low-Mach vortex gives an \f$L_2\f$
convergence rate of the conserved vector of \f$\approx 2.0\f$ in this regime
(see `examples/2D_isentropic_vortex_lowmach_tconv`).

The forward–backward acoustic micro-step is itself a symplectic-Euler update, which is
first-order accurate in \f$\Delta\tau\f$. With the recommended auto substep count
(\f$n_s \approx (|u|+c)/|u|\f$ set by the Mach number), \f$\Delta\tau \propto \Delta t\f$,
so refining \f$\Delta t\f$ refines \f$\Delta\tau\f$ in lockstep and the practical
convergence rate is first order — the standard behaviour of a split-explicit /
forward-backward scheme. The acoustic substep is therefore well suited to its design
purpose (resolving the fast waves cheaply at a fixed acoustic CFL) rather than to
asymptotic \f$\Delta t \to 0\f$ refinement.

## Usage

| Parameter | Type | Description |
|-----------|------|-------------|
| `acoustic_substepping` | Logical | Enable the split-explicit low-Mach integrator |
| `n_acoustic_substeps` | Integer | Fixed substep count; `0` auto-computes it each step (recommended) |
| `acoustic_div_damp` | Real | Dimensionless grad–div damping coefficient (default `0.1`; stable for \f$\lesssim 0.5/\text{num\_dims}\f$) |

The mode requires `model_eqns = 2` (5-equation model), a CFL-based time step
(`cfl_adap_dt` or `cfl_const_dt`), and `time_stepper = 3` (SSP-RK3). These constraints are
enforced at input checking. It is incompatible with bubbles, immersed boundaries,
(hypo/hyper)elasticity, chemistry, and phase change.

```python
"model_eqns": 2,
"time_stepper": 3,
"cfl_adap_dt": "T",
"cfl_target": 0.5,
"acoustic_substepping": "T",
"n_acoustic_substeps": 0,
"acoustic_div_damp": 0.1,
```

Multiple fluids (`num_fluids > 1`) are supported: the subcycle uses the stiffened-gas
mixture equation of state and mixture velocity, and advects the volume fractions on the
slow step. For `num_fluids = 1` the path reduces exactly to the single-fluid expressions.

The mode runs on CPU and on GPU through both OpenACC and OpenMP target offload, using the
backend-agnostic `GPU_*` macros (see @ref gpuParallelization "GPU Parallelization").

## Scope and limitations

- **Smooth low-Mach flow** and **material interfaces** are in scope. On a smooth field the
  scheme is bit-identical to the standard solver on the smooth tier and recovers the full
  \f$O(1/M)\f$ per-step speedup (a smooth 2D \f$M=0.1\f$ case measures \f$\approx 5.9\times\f$
  overall, with the blended cost ratio back near the all-smooth value once the conditional
  reconstruction skips the unflagged field). At a material interface — resolved or
  1-cell-sharp — the robust tier keeps the solution non-oscillatory (bounded few-percent
  pressure ringing) and species mass exactly conservative (drift \f$\sim 10^{-16}\f$), with
  the interface slightly more diffused than full HLLC.
- **Propagating shocks and shock–interface interactions are now captured** at flagged faces
  via the live full-HLLC two-tier flux (mass, energy, and convective momentum all live; see
  *Two-tier flux at discontinuities*). Validated cases (all with stable BCs, see below):
  - A 1D periodic shock tube (`examples/1D_shock_periodic_lowmach`): stable, non-oscillatory,
    species mass conserved to \f$\sim 10^{-16}\f$, matching full HLLC at the same resolution
    to \f$\approx 3\text{–}5\%\f$ \f$L_2\f$ on the conserved fields.
  - A 2D periodic shock with transverse shear (`examples/2D_shock_transverse_lowmach`),
    which exercises the off-normal momentum path: **both** momentum components match full
    HLLC to \f$\approx 5\text{–}7\%\f$ \f$L_2\f$, the transverse velocity stays uniform in
    the homogeneous direction to \f$\sim 10^{-15}\f$ (no spurious transverse momentum from
    the rotated-flux mapping), mass drift \f$\sim 10^{-16}\f$, no overshoot.
  - A 2D periodic multifluid shock–droplet (`examples/2D_shockdroplet_lowmach`,
    `num_fluids = 2`): stable, non-oscillatory, per-species mass conserved to
    \f$\sim 10^{-16}\f$, volume fraction bounded, matching full HLLC to \f$\approx 4\text{–}7\%\f$
    \f$L_2\f$.
  Accuracy caveat: because each acoustic micro-step is first-order in \f$\Delta\tau\f$ (the
  forward–backward update is symplectic-Euler; the macro coupling is 2nd-order — see
  *Temporal accuracy*), the captured front is a few cells more dissipative than full HLLC,
  and the dissipation accumulates over a long run (the multifluid droplet, run for many more
  steps, shows the largest peak-amplitude loss). For `num_fluids > 1` the face EOS is rebuilt
  from the **cell** volume fractions (exact only for `num_fluids = 1`); this does not
  destabilize or introduce oscillations at the tested multifluid shock, but it adds
  dissipation there — reconstructing the mixture EOS coefficients at the face is a future
  accuracy refinement.
- **The centered smooth tier is unstable at extrapolation/outflow boundaries with a
  background flow.** A uniform outflow cell carrying a mean flow can go NaN, because the
  centered (non-upwind) acoustic transport has no boundary dissipation there. This is
  independent of the shock-capture machinery. Keep discontinuities away from outflow
  boundaries, or use stable BCs (periodic, or reflecting/slip walls) near them; all of the
  validated cases above are periodic for this reason. Robustifying the smooth tier at
  outflow boundaries is future work.
- The macro (WS-RK3) coupling is second-order in the advective step at fixed acoustic CFL;
  with auto substeps (\f$\Delta\tau \propto \Delta t\f$) the symplectic-Euler micro-step
  makes the practical convergence first order (see *Temporal accuracy* above). Spatial
  order is unaffected.
- Restricted to `model_eqns = 2`.

## Source files

- `src/simulation/m_acoustic_substep.fpp` — forward–backward subcycle kernel
- `src/simulation/m_time_steppers.fpp` — `s_split_explicit_rk` orchestration and the two-CFL `s_compute_dt`
- `src/simulation/m_riemann_solver_hllc.fpp` — slow-flux variant
