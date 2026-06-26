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

- Intended for **smooth** low-Mach flow. The acoustic substep is centered (non-upwinded),
  so flows with embedded shocks are out of scope.
- The macro (WS-RK3) coupling is second-order in the advective step at fixed acoustic CFL;
  with auto substeps (\f$\Delta\tau \propto \Delta t\f$) the symplectic-Euler micro-step
  makes the practical convergence first order (see *Temporal accuracy* above). Spatial
  order is unaffected.
- Restricted to `model_eqns = 2`.

## Source files

- `src/simulation/m_acoustic_substep.fpp` — forward–backward subcycle kernel
- `src/simulation/m_time_steppers.fpp` — `s_split_explicit_rk` orchestration and the two-CFL `s_compute_dt`
- `src/simulation/m_riemann_solver_hllc.fpp` — slow-flux variant
