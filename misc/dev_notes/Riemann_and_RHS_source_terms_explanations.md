# Riemann Source Variable Explanation

## 0. Purpose

The Riemann solvers in `m_riemann_solvers.fpp` produce face-centered quantities
at every cell interface. Some of these are ordinary conservative fluxes; others
carry information needed by the non-conservative (NC) source terms, which are
computed separately in `m_rhs.fpp`. This note explains which arrays carry what,
how they get from the solver to the RHS, and how the meaning changes depending
on the solver mode.

Companion document: `notes/explanations/hll_hllc_nc_terms_notes.md` (mathematical
derivations of the NC constructions).

All examples below use the x-sweep. The y/z sweeps follow the same logic with
index permutations.

## 1. Three NC Advection Modes

The volume-fraction advection equation has a non-conservative source term, and
the Riemann solver must export the face quantity needed to discretize it. The
code selects exactly one of three modes at initialization
(`m_global_parameters.fpp`):

| Flag | True when | What the solver exports in `flux_src_rs(adv*)` |
|------|-----------|------------------------------------------------|
| `adv_src_alpha_iface` | HLL Method 1 (`riemann_solver == 1, .not. hll_u_interface`) | Per-fluid interface $\Psi_{\alpha_k}$ in slots `adv%beg:adv%end` |
| `adv_src_vel_iface` | HLL Method 2, HLLC, Exact, LF | Shared face-normal $\Psi_u$ in slot `adv%beg` only |
| `adv_src_none` | HLLD (`riemann_solver == 4`) | Zeros (all NC terms handled inside the solver) |

Notes:

- "Method 1" and "Method 2" refer only to the NC advection construction
  (see companion doc §§1–2).
- `flux_src_rs` is **overloaded**: its content depends on the active mode.
  Always read it together with the mode flag and the slot index.

## 2. Face Arrays and Handoff

The Riemann solver computes face quantities in **Riemann-space** arrays (local
to `m_riemann_solvers`, indexed as `(j, k, l, component)`). These are then
copied into **physical-space** arrays that `m_rhs` reads. The copy step
("finalization") handles the index permutation between sweep-local and
physical ordering.

There are four face arrays, plus one that is local-only:

| Array (x-sweep) | Content | Finalized to | Read by `m_rhs` as |
|------------------|---------|-------------|---------------------|
| `flux_rsx_vf(j,k,l, 1:sys_size)` | Conservative face flux $\hat{F}$ | `flux_vf` | `flux_n(dir)%vf(i)` |
| `flux_src_rsx_vf(j,k,l, adv%beg:...)` | NC advection handoff (overloaded, see §1) | `flux_src_vf` | `flux_src_n(dir)%vf(...)` |
| `nc_iface_vel_rsx_vf(j,k,l, 1:num_dims)` | Face velocities for NC terms that `flux_src` cannot carry (see §3.3) | `nc_iface_vel_vf` | `nc_iface_vel_n(dir)%vf(...)` |
| `flux_gsrc_rsx_vf(j,k,l, 1:sys_size)` | Geometric source flux (axisym/3D-cyl only) | `flux_gsrc_vf` | `flux_gsrc_n(dir)%vf(i)` |
| `vel_src_rsx_vf(j,k,l, 1:num_vels)` | HLLC star-state face velocities (local only) | — | Not directly; used by viscous/capillary source builders |

In addition, `flux_src_vf(mom%beg:E)` carries viscous/capillary momentum-energy
source fluxes, but these are **not** copied from `flux_src_rsx_vf`. They are
built directly in physical space by `s_compute_viscous_source_flux` and
`s_compute_capillary_source_flux`.

## 3. Array Details

### 3.1 `flux_src_rsx_vf` overloading

| Slot | `adv_src_alpha_iface` (HLL M1) | `adv_src_vel_iface` (HLL M2 / HLLC) | `adv_src_none` (HLLD) | Finalized? |
|------|------|------|------|------|
| `adv%beg` | First $\Psi_\alpha$ | $\Psi_u$ | 0 | Always |
| `adv%beg+1:adv%end` | Remaining $\Psi_\alpha$ per fluid | Unused | 0 | Only when `adv_src_alpha_iface .or. adv_src_none` |

When `adv_src_vel_iface` is active, `flux_src_n(dir)%vf(adv%beg+1:adv%end)` are
pointer-aliased to `flux_src_n(dir)%vf(adv%beg)` in `m_rhs.fpp`, so the RHS loop
reads the same shared velocity for all alpha equations without branching.

### 3.2 `vel_src_rsx_vf`

Local to the Riemann solver — never finalized. Stores HLLC star-state face
velocities. Used by:

- Viscous/capillary energy source builders (which read `vel_src_rs` directly
  inside `m_riemann_solvers`)
- HLLC advection export: the normal component is copied to
  `flux_src_rs(adv%beg)` as the exported $\Psi_u$

### 3.3 `nc_iface_vel_rsx_vf`

A second face-velocity export channel, independent of `flux_src`. Stores
$F_\text{HLL}(U{=}1,\,F{=}u_i)$ for all `num_dims` velocity components in
physical order. Finalized to `nc_iface_vel_n(dir)%vf(...)` via its own path.

It exists because `flux_src` cannot always carry the velocity traces that the
NC source terms need:

1. **Interface-consistent hypo RHS** (`hypo_nc_interface`): the stress evolution
   equation needs velocity gradients in all directions (e.g., $\partial_x v$
   from the x-sweep). `flux_src` only exports one scalar per face, so the
   full velocity vector goes here.

2. **HLL M1 + K*div(u)** (`adv_src_alpha_iface .and. alt_soundspeed`):
   `flux_src(adv*)` is already occupied by per-fluid $\Psi_\alpha$, so the
   velocity trace needed for K*div(u) must go in a separate channel.

3. **HLLD + axisymmetric** (`hypo_nc_dual_pass .and. grid_geometry == 2`):
   the axisymmetric geometric source needs face velocities exported from the
   dual-pass solver.

For HLLC hypo ADC, the code updates all three export channels consistently:
`vel_src_rs`, `flux_src_rs(adv%beg)`, and `nc_iface_vel_rs`.

### 3.4 `flux_gsrc_rsx_vf`

Geometric source flux for cylindrical coordinates. Only active for the y-sweep
(axisymmetric, `cyl_coord`) and z-sweep (3D cylindrical, `grid_geometry == 3`).
The x-sweep buffer exists but is not exported.

## 4. Finalization Summary

| Source | Destination | Condition |
|--------|-------------|-----------|
| `flux_rs(1:sys_size)` | `flux_vf(1:sys_size)` | Always |
| `flux_src_rs(adv%beg)` | `flux_src_vf(adv%beg)` | Always |
| `flux_src_rs(adv%beg+1:adv%end)` | `flux_src_vf(adv%beg+1:adv%end)` | `adv_src_alpha_iface .or. adv_src_none` |
| `nc_iface_vel_rs` | `nc_iface_vel_vf` | `use_nc_iface_vel` |
| `flux_gsrc_rs` | `flux_gsrc_vf` | `cyl_coord` (y-sweep) or `grid_geometry == 3` (z-sweep) |
| `vel_src_rs` | — | Never finalized |

## 5. What `m_rhs` Reads

### 5.1 Advection source (`s_compute_advection_source_term`)

Three branches, selected by the `adv_src_*` flag:

| Branch | Face quantity read | RHS formula per $\alpha_k$ | K*div(u) velocity source |
|--------|-------------------|---------------------------|--------------------------|
| `adv_src_alpha_iface` | `flux_src_n(dir)%vf(j_adv)` = per-fluid $\Psi_{\alpha_k}$ | $u_\text{cell} \cdot \Delta\Psi_\alpha / \Delta x$ | `nc_iface_vel_n(dir)%vf(1)` |
| `adv_src_vel_iface` | `flux_src_n(dir)%vf(adv\%beg)` = shared $\Psi_u$ | $\alpha_k \cdot \Delta\Psi_u / \Delta x$ | Same `flux_src_n` slot (already $\Psi_u$) |
| `adv_src_none` | — | Skipped (HLLD handles internally) | — |

The K*div(u) correction (active when `alt_soundspeed = T`) adds
$\pm K \cdot \Delta\Psi_u / \Delta x$ to `adv%beg` / `adv%end`.
Its velocity source differs between the first two branches because in M1 the
`flux_src` slots are occupied by alpha traces.

### 5.2 Viscous / capillary / chemistry (`s_compute_additional_physics_rhs`)

Reads `flux_src_n(dir)%vf(mom%beg:E)`. Independent of the advection mode —
this channel is built in physical space by the viscous/capillary source builders.

Surface tension color-function transport also reads `flux_src_n(dir)%vf(adv%beg)`,
which requires $\Psi_u$ in that slot — so it is only meaningful with
`adv_src_vel_iface`.

### 5.3 Hypoelastic stress source

Three NC modes for the stress evolution equation, selected by solver:

| Solver | User parameter | Derived flag | Velocity source for $\nabla u$ |
|--------|---------------|--------------|-------------------------------|
| HLL (default) | `hypo_hll_interface_rhs = F` | `hypo_nc_finite_diff` | Cell-center FD stencil (`s_compute_hypoelastic_rhs_finite_diff_per_sweep`) |
| HLL (optional) | `hypo_hll_interface_rhs = T` | `hypo_nc_interface` | `nc_iface_vel_n` (`s_compute_hypoelastic_rhs_iface`) |
| HLLC | — | `hypo_nc_interface` | `nc_iface_vel_n` (`s_compute_hypoelastic_rhs_iface`) |
| HLLD | — | `hypo_nc_dual_pass` | Built into Riemann flux (no separate RHS routine) |

**Code organization (three shapes).** These three stress treatments also live in three
different code *shapes*: HLLC and HLL add their hypoelastic handling as inline
`if (hypoelasticity)` branches inside `s_hllc_riemann_solver` / `s_hll_riemann_solver`;
HLLD is a separate module `m_riemann_solver_hypo_hlld`, reached via the `hypo_nc_dual_pass`
path in `s_riemann_solver`. HLLD is separate because its anchored dual pass produces both the
hat_L and hat_R anchored flux sets in one fused solve (whose partial RHS are summed in
`m_rhs`), a different control flow from the single-pass inline branches the other two use.

### 5.4 `use_nc_iface_vel` — allocation condition

```
use_nc_iface_vel = hypo_nc_interface
                   .or. (hypo_nc_dual_pass .and. grid_geometry == 2)
                   .or. (adv_src_alpha_iface .and. alt_soundspeed)
```

| Case | Why `nc_iface_vel` is needed |
|------|------------------------------|
| `hypo_nc_interface` | Stress RHS needs all velocity components ($\partial_x v$, $\partial_y u$, etc.) |
| `hypo_nc_dual_pass .and. grid_geometry == 2` | HLLD axisym geometric source needs exported face velocities |
| `adv_src_alpha_iface .and. alt_soundspeed` | `flux_src` is occupied by $\Psi_\alpha$; K*div(u) needs $\Psi_u$ elsewhere |
