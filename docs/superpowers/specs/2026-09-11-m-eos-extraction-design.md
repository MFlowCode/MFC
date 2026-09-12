# Extracting `m_eos` from `m_variables_conversion`

Date: 2026-09-11
Branch: `module-eos`
Status: design approved, implementation plan pending

## Context

`src/common/m_variables_conversion.fpp` is 1997 lines. Roughly a third of it is
equation-of-state machinery, most of it added by three PRs on master:

- #1762 `c9c8b0f6` — centralized the EOS expressions in the Riemann solvers
- #1808 `80100576` — cases set only the parameters their EOS reads
- #1811 `808df619` — state-dependent EOS: Mie-Gruneisen, JWL, Vinet

That work established a single abstraction: every EOS is expressed in Gamma/Pi form,
`rho e = Gamma(rho) p + Pi(rho)`, so the solvers never learned a new equation of state --
they stopped treating `gamma`/`pi_inf` as constants. The code is coherent; it is simply
living in the wrong file.

## Goals

- Move the EOS machinery into `src/common/m_eos.fpp` as a leaf module shared by all
  three executables.
- Keep the change a **pure move**: no numerical change, no kernel-count change.
- Leave the per-phase API surface as small as the callers actually require.

## Non-goals

- Changing any EOS formulation, adding a family, or altering solver behaviour.
- Moving device globals or their `GPU_DECLARE` clauses.
- Touching the `cray_inline` macro chain (see Known issues).
- The family-registry work (see Follow-ups) -- that is a separate PR.

## Approach

### Chosen: one module, thin waist

`m_eos` absorbs the whole EOS chain, so the deep call paths stay *internal* to it.
`s_phase_pressure_on_isentrope -> s_rk4 -> s_ode_slope -> s_reference_curve` is four
levels deep, but only the top call crosses a module boundary; the lower three keep
their current bare directives and need no change.

### Rejected: split curves from closure

A `m_eos_curves` + `m_eos` pair reads better and tests better, but it puts
`s_reference_curve` behind a module boundary from `s_ode_slope` -- the innermost call
of every RK4 step, 8 steps x 4 stages = 32 crossings per evaluation. A cross-module
`!$acc routine seq` call is a real device call on every backend absent LTO/IPA. The
hot inner loop is the wrong place to put an interface.

Recorded here so it is not re-proposed later.

## Module boundary

`src/common/m_eos.fpp` uses `m_derived_types`, `m_constants`, and
`m_global_parameters_common` -- nothing else. It does not need `m_mpi_proxy`,
`m_helper`, or `m_thermochem`, which `m_variables_conversion` does. It therefore sits
*below* the per-target `m_global_parameters`, the way `m_variables_conversion` already
reaches at the common globals directly for `shear_indices`.

| Tier | Routines | Visibility |
| --- | --- | --- |
| Gamma/Pi primitives | `f_pressure`, `f_bulk_modulus`, `f_relativistic_enthalpy`, `f_isentrope_exponent`, `f_isentrope_pressure`, `f_sg_thermal`, `f_c2_from_coefficients` | public |
| Family layer | `s_reference_curve`, `s_eos_coefficients`, `s_ode_slope`, `s_rk4`, `s_phase_c2`, `f_has_isentropic_reference`, `f_hugoniot_compression_limit` | private |
| Per-phase API | `f_is_state_dependent`, `s_phase_coefficients`, `s_phase_pressure_on_isentrope`, `s_phase_temperature`, `s_phase_density_on_isentrope`, `s_phase_internal_energy`, `s_phase_bulk_modulus` | public |

`s_eos_coefficients` and `f_c2_from_coefficients` are exported by
`m_variables_conversion` today but have no caller outside it. Making them private to
`m_eos` is free and shrinks the exported surface.

## What stays in `m_variables_conversion`

The mixture *closure rules* are not equations of state and stay put:

- `s_compute_mixture_coefficients`, `s_compute_mixture_coefficients_dt`
- `s_compute_energy`, `s_compute_pressure`
- `s_compute_speed_of_sound`, `s_compute_speed_of_sound_avg`,
  `s_compute_fast_magnetosonic_speed`
- all `s_convert_*` routines

`s_compute_speed_of_sound` additionally reads the `chemistry`, `relativity` and `mhd`
switches, which would drag `m_thermochem` back into the leaf module.

`f_elastic_energy` and `f_hypoelastic_energy` also stay. They sit inside the EOS region
of the file today but are hypoelasticity, not equation of state.

## Device globals: unmoved

`eoss`, `eos_coeffs`, `gammas`, `cvs`, `isentrope_n`, `isentrope_B`,
`any_state_dependent_eos` and their `$:GPU_DECLARE(create=...)` stay in
`m_global_parameters_common`. Moving them buys nothing and walks into two documented
traps: the Cray `ftn-7066 Global in accelerator routine without declare` failure, and
amdflang whole-image device codegen instability.

## Initialization split

The EOS portion of `s_initialize_variables_conversion_module` becomes
`s_initialize_eos_module`, called immediately before it from each target's start-up.
It carries `gammas`, `eoss`, `isentrope_n`, `isentrope_B`, `pi_infs`, `cvs`, `qvs`,
`qvps`, `eos_coeffs`, `any_state_dependent_eos`, and the `@:PROHIBIT` case-optimization
consistency check.

**The loop is interleaved and does not cut cleanly.** `Gs_vc(i)` (shear modulus, private
to `m_variables_conversion`) is assigned in the middle of the EOS coefficient block, and
`Gs_vc` shares the single closing
`$:GPU_UPDATE(device='[gammas, ..., Gs_vc, eoss, eos_coeffs]')`. That update must be
split in two, with `Gs_vc` and `Res_vc` staying behind. This is the one part of the move
that is not mechanical.

## Cross-module directive cleanup

`s_phase_pressure_on_isentrope` and `s_phase_temperature` are already called
cross-module on master (`m_riemann_solver_hllc`, `m_reactive_burn`,
`post_process/m_start_up`) with a bare `parallelism='[seq]'`. Both gain
`cray_inline=True` and `function_name=` to match the convention every other
cross-module device helper here follows.

This cannot perturb results on any GPU build -- see Known issues -- so the bit-for-bit
guarantee below survives it.

## Verification

A pure move must produce **bit-for-bit identical goldens**. No tolerance changes, no
regeneration. That is the acceptance test for correctness, and it is stronger than any
test that could be written for the new module.

Kernel count is unchanged, so an amdflang A/B between this branch and its merge base is
not confounded by the whole-image codegen effect that makes wall-time comparisons
meaningless across commits differing in target-region count. GPU benchmark lanes should
therefore be directly comparable.

The one genuinely unmeasured thing: nested `cray_inline=True` across a module boundary
(`s_compute_mixture_coefficients` -> `s_phase_coefficients`, both forced-inline, now in
different files). Same-module today, and no precedent elsewhere in the repo.

## Known issues (recorded, not acted on)

`cray_inline=True` expands to `!DIR$ INLINEALWAYS` **only** on `_CRAYFTN` builds with
neither offload backend. On any GPU build the macro chain at
`src/common/include/parallel_macros.fpp:86` selects the acc/omp directive instead, and
`ACC_ROUTINE` discards `function_name` entirely -- so `cray_inline=True` and a bare
`parallelism='[seq]'` emit byte-identical directives on GPU. This is a known issue and
is explicitly out of scope here; `docs/documentation/gpuParallelization.md:637` describes
it in GPU terms and is likewise left alone.

The consequence for this design is only that the module boundary is not a GPU
performance interface, which is why the rejection of the split-module alternative rests
on device-call cost on all backends rather than on a Cray-specific directive.

## Follow-ups (separate PRs)

1. **Family registry.** Adding an EOS family touches ~8 places: `m_constants`,
   `m_derived_types`, the init `select case`, `s_reference_curve`,
   `f_is_state_dependent`, `f_has_isentropic_reference`,
   `toolchain/mfc/params/definitions.py`, `toolchain/mfc/eos.py`. Since
   `toolchain/mfc/params/generators/fortran_gen.py` already generates Fortran from
   Python, one registry could emit the `eos_*` constants, the `physical_parameters`
   fields, both predicates, and the validator's parameter sets. This is the real DRY
   defect; the line count is not.

2. **Per-fluid case-optimized `eoss`.** `any_state_dependent_eos` already becomes a
   `parameter` under `--case-optimization`, but `eoss(i)` stays allocatable, so a
   two-fluid case with one Mie-Gruneisen and one stiffened gas still branches at runtime
   on the cheap fluid. Baking per-fluid `eoss` into fypp parameters would let the
   compiler drop the dead family per loop.
