# Non-Newtonian (Herschel–Bulkley) Viscosity — Design

**Date:** 2026-06-08
**Status:** Approved design, pre-implementation
**Supersedes:** PR #1298 (closed — branch too stale to merge; master had since refactored the
Riemann reconstruction buffers and IB derived types out from under it). This is a fresh
implementation on current `master` that reuses the validated pieces of #1298.

## Goal

Add per-fluid Herschel–Bulkley (HB) non-Newtonian viscosity with Papanastasiou
regularization to the `simulation` target. Covers power-law, Bingham, and general HB
fluids. The effective viscosity varies spatially with the local strain rate, recomputed
every step. Backward compatible: all parameters default to Newtonian behavior.

HB effective viscosity:

```
mu_eff(gamma_dot) = (tau0 / gamma_dot) * (1 - exp(-hb_m * gamma_dot)) + K * gamma_dot^(nn - 1)
gamma_dot = sqrt(2 * D_ij * D_ij)   (second invariant of the strain-rate tensor D)
```

Special cases: `tau0=0, nn=1` → Newtonian `mu=K`; `tau0=0, nn≠1` → power-law;
`nn=1, tau0>0` → Bingham. All HB parameters are nondimensional (scaled by
`rho_ref*U_ref*L_ref`), so `1/mu_eff` equals the local effective Reynolds number, matching
how `fluid_pp(i)%Re(1)` is used for Newtonian fluids.

## Key insight that shapes the whole design

In MFC's existing machinery a fluid's contribution to the mixture viscosity is
`alpha_q / Re_q = alpha_q * mu_q`, with `mu_q = 1/Re_q` **constant** for a Newtonian fluid.
Non-Newtonian simply makes `mu_q` a function of the local shear rate, `mu_q(gamma_dot)`.
The mixture rule, the viscous stress assembly, the flux source, and the CFL form are all
**unchanged** — only the per-fluid `mu_q` value changes.

Crucially, on current `master` the viscous path **already computes the velocity gradients**
we need (`s_get_viscous` in `m_viscous.fpp` passes `dqL_prim_dx_vf/dy/dz` and `dqR_*` into
the Riemann solvers), and **already maps rotated `(j,k,l)` loop indices to physical
coordinates** via the existing `idx_right_phys` helper
(`m_riemann_solvers.fpp:~1386`). We assemble the strain-rate tensor from those
already-available, already-correctly-indexed gradients. We therefore do **not** create a
gradient-recomputation module — that module (`m_re_visc` in PR #1298) is exactly where the
rotated→physical index-mapping bug lived. By piggybacking on the solver's own gradient
indexing, that bug class cannot recur.

## Architecture (Approach 1: in-place substitution, reuse solver gradients)

Two substitution points, one reused formula module, one reused shear-rate helper.

### Reused module: `src/simulation/m_hb_function.fpp` (from #1298, already validated)
- `f_compute_hb_viscosity(tau0, K, nn, mu_min, mu_max, gamma_dot, hb_m)` — pure
  `GPU_ROUTINE(parallelism='[seq]')`. Papanastasiou-regularized HB formula with a
  `verysmall` zero-shear guard and correct analytic limits; clamps to `[mu_min, mu_max]`
  when set (sentinel-checked against `dflt_real`). No hard shear-rate clamp.
- `f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)` — pure
  `GPU_ROUTINE(parallelism='[seq]')`. `gamma_dot = sqrt(2 D:D)`; absent dimensions pass 0.

Strain-rate tensor from velocity gradients:
```
D_xx = du/dx   D_yy = dv/dy   D_zz = dw/dz
D_xy = 0.5*(du/dy + dv/dx)   D_xz = 0.5*(du/dz + dw/dx)   D_yz = 0.5*(dv/dz + dw/dy)
```

### Data model and GPU strategy
New per-fluid parameters (8), registered the current-master way:

| param | type | default | meaning |
|---|---|---|---|
| `non_newtonian` | logical | `F` | enable HB for this fluid |
| `K` | real | `dflt_real` | consistency index (required if NN) |
| `nn` | real | `dflt_real` | flow behavior index (required if NN) |
| `tau0` | real | `0` | yield stress (0 ⇒ power-law) |
| `hb_m` | real | `dflt_real` | Papanastasiou param (required if `tau0>0`) |
| `mu_min` | real | `dflt_real` (inactive sentinel) | lower viscosity clamp |
| `mu_max` | real | `dflt_real` (inactive sentinel) | upper clamp (required if NN) |
| `mu_bulk` | real | `dflt_real` (inactive sentinel) | bulk viscosity for NN |

Registration locations (current master; namelist binding is auto-generated at cmake
configure — no `m_start_up.fpp` edit):
1. `toolchain/mfc/params/definitions.py` — `_r()` per param inside the `fluid_pp` loop.
2. `toolchain/mfc/params/descriptions.py` — pattern entries.
3. `src/common/m_derived_types.fpp` — fields on `physical_parameters`.
4. `src/{pre_process,simulation,post_process}/m_global_parameters.fpp` — default init (×3).
5. `src/{pre_process,simulation,post_process}/m_mpi_proxy.fpp` — MPI broadcast (×3).
6. `toolchain/mfc/case_validator.py` — physics guards (see below).

**GPU strategy — mirror the existing `Res_gs`/`Res_viscous` pattern.** Master precomputes
those as plain module arrays so GPU kernels never dereference `fluid_pp%...`. We do the
same: at init in `s_initialize_global_parameters_module`, populate `GPU_DECLARE`'d
module-level arrays `is_non_newtonian(num_fluids)` and
`hb_tau0/hb_K/hb_nn/hb_m/hb_mu_min/hb_mu_max/hb_mu_bulk(num_fluids)`. Kernels read those.
This is the established idiom and avoids the cross-compiler derived-type-on-device
fragility #1298 fought.

**`any_non_newtonian`** — a single global logical, set at init by scanning the fluids,
`GPU_DECLARE`'d. It gates every new branch.

### Substitution point 1 — Cartesian Riemann (`m_riemann_solvers.fpp`)
Where the per-fluid loop builds `Re_L(i) = sum_q alpha_L(q)/Res_gs(i,q)` (and `Re_R`),
branch on `is_non_newtonian(q)`:
- **shear (i=1):** Newtonian → `1/Res_gs(1,q)` (unchanged); NN →
  `mu_q = f_compute_hb_viscosity(hb_*(q), gamma_dot)`. `gamma_dot_L` / `gamma_dot_R`
  computed **once per interface** from the L / R gradients (assembled from
  `dqL_prim_dx_vf(mom%beg+c-1)%sf(...)` etc.; right-state via the existing `idx_right_phys`
  mapping), reused across fluids.
- **bulk (i=2):** Newtonian → `1/Res_gs(2,q)`; NN → `hb_mu_bulk(q)` if set, else 0.

The harmonic-mean mixture, `Re_avg_rsx_vf`, and the downstream
`s_compute_cartesian_viscous_source_flux` (~line 4100) are **untouched** — they consume
`Re_L/Re_R`, so the substitution propagates automatically to HLL/HLLC/HLLD.

### Substitution point 2 — Cylindrical viscous path (`m_viscous.fpp`)
`s_compute_viscous_stress_cylindrical_boundary` builds `Re_visc(i) = 1/sum_q
alpha/Res_viscous` the same way — same substitution, `gamma_dot` from the gradients
available there. The axis-singularity geometry handling is orthogonal to the viscosity
value and stays as-is.

### IBM (`m_ibm.fpp`)
Force/torque is built by calling `s_compute_viscous_stress_tensor(..., dynamic_viscosity,
i±1, j, k)` at stencil neighbors. #1298's bug: it computed `mu` once at the center cell and
reused it for every sample. Fix: **evaluate `mu` per stencil sample from that sample's own
local shear rate.** Preferred form — evaluate `mu` inside `s_compute_viscous_stress_tensor`
from the gradients it computes for the cell (confirm in planning whether the routine
computes its gradients internally or receives them; that decides inside-routine vs
per-call-site). IBM runs outside the flux-gradient context, so this uses plain physical
`(i,j,k)` central differences — no rotation, no index-mapping risk.

### CFL / dt (`m_sim_helpers.fpp` / `m_time_steppers.fpp`)
Viscous limit is `vcfl_dt = cfl * dx^2 / max(nu)`, `nu = mu/rho`. For a non-Newtonian fluid
use its **upper clamp `mu_max`** as the bounding viscosity. Since `mu(gamma_dot) <= mu_max`
everywhere, the dt is guaranteed stable (at worst slightly conservative) and needs no
gradient plumbing into the time-stepper. This is why `mu_max` is required for all NN fluids.

## Newtonian-invariance invariant (non-negotiable)
When `any_non_newtonian = .false.`, every modified site takes the existing `Res_gs` /
`Res_viscous` path — bitwise identical to today. **No existing golden file may change.**
Verified empirically: run the current viscous regression tests without `--generate` and
confirm zero diff before generating any new goldens.

## Guards (placement principle)
- **Python `case_validator.py`** — everything inferable from the case parameters before the
  solver runs (fail fast at `./mfc.sh validate`):
  - `non_newtonian` ⇒ `viscous`
  - `non_newtonian` ⇒ `K` and `nn` set
  - `tau0 > 0` ⇒ `hb_m` set
  - `non_newtonian` ⇒ `mu_max` set
  - `mu_max > mu_min` when both set
  - `non_newtonian` incompatible with `igr`
  - `model_eqns` restricted to the values where the multi-fluid viscous path is valid
    (confirm 2/3 in planning)
- **Fortran `m_checker.fpp`** — only genuinely runtime-determined constraints. None expected
  for this feature; add one only if planning surfaces a runtime-only constraint.

## Validation anchors
- **`examples/2D_poiseuille_nn/`** — power-law Poiseuille. Closed-form `u(y)` exists; the
  case carries the analytic comparison, giving a quantitative L2 error. This is the trust
  anchor, checked **before** generating any golden.
- **`examples/2D_lid_driven_cavity_nn/`** — shear-thinning and shear-thickening cavity,
  documented against Li et al. (2015), as a qualitative example.

## Tests
`toolchain/mfc/test/cases.py` — a small, fast `stack.push("Non-Newtonian", {...})` block:
power-law shear-thinning and shear-thickening on a tiny grid, 1D/2D. Goldens generated only
after the Poiseuille analytic check passes, so they encode verified output. Existing viscous
goldens remain untouched (the invariance check above).

## Docs
- `docs/documentation/case.md` — non-Newtonian section (reuse #1298's, which was good).
- `docs/module_categories.json` — add `m_hb_function` only. **No `m_re_visc`** — that module
  does not exist in this design.

## Build / compiler matrix
Must compile and behave identically under gfortran, nvfortran, Cray ftn, Intel ifx (CI), and
AMD flang for OpenMP-offload GPU. GPU via `GPU_*` Fypp macros only. The new device routines
are pure `GPU_ROUTINE(parallelism='[seq]')`. CPU-only build must work (every `#ifdef` covers
CPU/ACC/OMP, with/without MPI).

## Out of scope (explicit)
- No `m_re_visc`-style gradient recomputation module.
- No new halo exchange / MPI buffers (Approach 1 reuses existing gradients).
- No precomputed full-domain `mu_eff` field (that was Approach 2, not chosen).

## Files touched (summary)
- New: `src/simulation/m_hb_function.fpp`; `examples/2D_poiseuille_nn/`,
  `examples/2D_lid_driven_cavity_nn/`.
- Modified: `m_riemann_solvers.fpp`, `m_viscous.fpp`, `m_ibm.fpp`, `m_sim_helpers.fpp`
  (and/or `m_time_steppers.fpp`), `m_global_parameters.fpp` (×3),
  `m_mpi_proxy.fpp` (×3), `m_derived_types.fpp`, `definitions.py`, `descriptions.py`,
  `case_validator.py`, `cases.py`, `docs/documentation/case.md`,
  `docs/module_categories.json`.
