# Non-Newtonian (Herschel–Bulkley) Viscosity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-fluid Herschel–Bulkley non-Newtonian viscosity (Papanastasiou-regularized) to the MFC `simulation` target, with spatially varying effective viscosity computed from the local strain rate, fully backward-compatible with Newtonian runs.

**Architecture:** A fluid's mixture contribution is `alpha_q * mu_q`; non-Newtonian makes `mu_q` a function of the local shear rate. We reuse the velocity gradients the viscous path already computes, assemble the strain-rate tensor, evaluate `mu_q(gamma_dot)` via a validated HB formula, and feed it through MFC's *existing* mixture/stress/CFL machinery. A single GPU device helper (`s_compute_mixture_inv_re`) encapsulates the per-interface mixture so the 4 Riemann + 4 cylindrical-viscous call sites share one code path (no copy-drift). When no fluid is non-Newtonian, every path is bitwise-identical to today.

**Tech Stack:** Modern Fortran 2008+ with Fypp macros; `GPU_*` macro API (OpenACC/OpenMP); Python toolchain (`./mfc.sh`); golden-file regression tests.

**Reference spec:** `docs/superpowers/specs/2026-06-08-non-newtonian-viscosity-design.md`

---

## Conventions for every task

- All commands run from repo root: `/storage/project/r-sbryngelson3-0/sbryngelson3/mfc`.
- Build: `./mfc.sh build -j 8` (all 3 targets, CPU). GPU is CI-gated; if nvfortran is available locally, also `./mfc.sh build -j 8 -t simulation --gpu acc`.
- Format before every commit: `./mfc.sh format -j 8`. Lint: `./mfc.sh precheck -j 8`.
- Branch is already `feature/non-newtonian-viscosity` on remote `fork` (sbryngelson/MFC).
- **Precision:** only `wp`/`stp` kinds, generic intrinsics, `1.0_wp` literals. **GPU:** only `GPU_*` macros. **Abort:** only `@:PROHIBIT`/`s_mpi_abort`.
- Commit messages: imperative, scoped (e.g. `feat:`, `test:`, `docs:`). Use `git commit -m "..."` (no heredocs).

---

## File structure

**New:**
- `src/simulation/m_hb_function.fpp` — HB viscosity formula, shear-rate helper, and the mixture helper. Pure GPU device routines only.
- `examples/2D_poiseuille_nn/case.py` — power-law Poiseuille (analytic anchor).
- `examples/2D_lid_driven_cavity_nn/case.py` — shear-thinning + shear-thickening cavity (qualitative example).

**Modified:**
- `toolchain/mfc/params/definitions.py`, `toolchain/mfc/params/descriptions.py`, `toolchain/mfc/case_validator.py`, `toolchain/mfc/test/cases.py`
- `src/common/m_derived_types.fpp`
- `src/{pre_process,simulation,post_process}/m_global_parameters.fpp` (×3)
- `src/{pre_process,simulation,post_process}/m_mpi_proxy.fpp` (×3)
- `src/simulation/m_riemann_solvers.fpp`, `src/simulation/m_viscous.fpp`, `src/simulation/m_ibm.fpp`, `src/simulation/m_sim_helpers.fpp`
- `docs/documentation/case.md`, `docs/module_categories.json`

---

## Phase 0: Baseline — capture Newtonian invariance reference

### Task 0: Record the current viscous golden-test state

**Files:** none modified.

- [ ] **Step 1: List the existing viscous regression tests**

Run: `./mfc.sh test -l 2>/dev/null | grep -i visc | tee /tmp/nn_baseline_visc_tests.txt`
Expected: a non-empty list of viscous test UUIDs/traces. These must NOT change at the end.

- [ ] **Step 2: Run the viscous tests on a clean build to confirm green baseline**

Run: `./mfc.sh build -j 8 && ./mfc.sh test -j 8 --only Viscous`
Expected: all pass. If any fail pre-change, STOP and report — the baseline must be green first.

- [ ] **Step 3: No commit** (read-only baseline).

---

## Phase 1: Parameters, data model, and the HB module (no behavior change)

After this phase the code builds, all existing tests pass unchanged, and the HB module exists and is unit-exercised — but nothing reads the new params yet, so Newtonian invariance is trivially preserved.

### Task 1: Add the 8 parameters to the derived type and Python definitions

**Files:**
- Modify: `src/common/m_derived_types.fpp:337-345`
- Modify: `toolchain/mfc/params/definitions.py` (fluid_pp `_r` loop)
- Modify: `toolchain/mfc/params/descriptions.py` (PATTERNS)

- [ ] **Step 1: Add fields to `physical_parameters`**

In `src/common/m_derived_types.fpp`, change the type (after `real(wp) :: G`):

```fortran
    type physical_parameters
        real(wp)               :: gamma   !< Sp. heat ratio
        real(wp)               :: pi_inf  !< Liquid stiffness
        real(wp), dimension(2) :: Re      !< Reynolds number
        real(wp)               :: cv      !< heat capacity
        real(wp)               :: qv      !< reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
        real(wp)               :: qvp     !< reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
        real(wp)               :: G
        logical                :: non_newtonian  !< Enable Herschel-Bulkley non-Newtonian viscosity
        real(wp)               :: K              !< HB consistency index
        real(wp)               :: nn             !< HB flow behavior index
        real(wp)               :: tau0           !< HB yield stress (0 => power-law)
        real(wp)               :: hb_m           !< Papanastasiou regularization parameter
        real(wp)               :: mu_min         !< Lower viscosity clamp (inactive sentinel = dflt_real)
        real(wp)               :: mu_max         !< Upper viscosity clamp (required when non_newtonian)
        real(wp)               :: mu_bulk        !< Bulk viscosity for NN (inactive sentinel = dflt_real)
    end type physical_parameters
```

- [ ] **Step 2: Register the parameters in `definitions.py`**

Locate the `fluid_pp` loop (where `{px} = f"fluid_pp({f})%"` and existing `_r(f"{px}G", REAL, ...)` lives). Add inside that loop, after the `G` registration:

```python
            _r(f"{px}non_newtonian", LOG, {"viscosity"}, math=r"\mathrm{non\text{-}Newtonian}_k")
            _r(f"{px}K", REAL, {"viscosity"}, math=r"K_k")
            _r(f"{px}nn", REAL, {"viscosity"}, math=r"n_k")
            _r(f"{px}tau0", REAL, {"viscosity"}, math=r"\tau_{0,k}")
            _r(f"{px}hb_m", REAL, {"viscosity"}, math=r"m_k")
            _r(f"{px}mu_min", REAL, {"viscosity"}, math=r"\mu_{\min,k}")
            _r(f"{px}mu_max", REAL, {"viscosity"}, math=r"\mu_{\max,k}")
            _r(f"{px}mu_bulk", REAL, {"viscosity"}, math=r"\mu_{\mathrm{bulk},k}")
```

(`LOG`, `REAL` are the `ParamType` aliases already imported in this file; match the exact alias names used for the neighboring `_r` calls — verify by reading the top of the loop.)

- [ ] **Step 3: Add description patterns in `descriptions.py`**

In the `PATTERNS` list add:

```python
    (r"fluid_pp\((\d+)\)%non_newtonian", "Enable Herschel-Bulkley non-Newtonian viscosity for fluid {0}"),
    (r"fluid_pp\((\d+)\)%K", "HB consistency index for fluid {0}"),
    (r"fluid_pp\((\d+)\)%nn", "HB flow behavior index for fluid {0}"),
    (r"fluid_pp\((\d+)\)%tau0", "HB yield stress for fluid {0}"),
    (r"fluid_pp\((\d+)\)%hb_m", "Papanastasiou regularization parameter for fluid {0}"),
    (r"fluid_pp\((\d+)\)%mu_min", "Lower viscosity clamp for fluid {0}"),
    (r"fluid_pp\((\d+)\)%mu_max", "Upper viscosity clamp for fluid {0}"),
    (r"fluid_pp\((\d+)\)%mu_bulk", "Bulk viscosity (non-Newtonian) for fluid {0}"),
```

- [ ] **Step 4: Verify the parameters are discoverable**

Run: `./mfc.sh params "non_newtonian" && ./mfc.sh params "mu_max"`
Expected: both resolve to `fluid_pp(i)%...` entries with the descriptions above.

- [ ] **Step 5: Commit**

```bash
./mfc.sh format -j 8
git add src/common/m_derived_types.fpp toolchain/mfc/params/definitions.py toolchain/mfc/params/descriptions.py
git commit -m "feat(params): register Herschel-Bulkley non-Newtonian fluid parameters"
```

### Task 2: Default-initialize the new fields in all three targets

**Files:**
- Modify: `src/pre_process/m_global_parameters.fpp`, `src/simulation/m_global_parameters.fpp`, `src/post_process/m_global_parameters.fpp` (the `do i = 1, num_fluids_max` fluid_pp default loop; in simulation it's at lines 498-499)

- [ ] **Step 1: Add defaults in each of the three files**

In each target's fluid_pp default loop, after `fluid_pp(i)%G = 0._wp`, add:

```fortran
            fluid_pp(i)%non_newtonian = .false.
            fluid_pp(i)%K = dflt_real
            fluid_pp(i)%nn = dflt_real
            fluid_pp(i)%tau0 = 0._wp
            fluid_pp(i)%hb_m = dflt_real
            fluid_pp(i)%mu_min = dflt_real
            fluid_pp(i)%mu_max = dflt_real
            fluid_pp(i)%mu_bulk = dflt_real
```

(Confirm the indentation/loop variable matches each file — pre/post may use a slightly different surrounding block; match local style.)

- [ ] **Step 2: Build**

Run: `./mfc.sh build -j 8`
Expected: all three targets compile (namelist auto-regenerates at cmake configure).

- [ ] **Step 3: Commit**

```bash
./mfc.sh format -j 8
git add src/pre_process/m_global_parameters.fpp src/simulation/m_global_parameters.fpp src/post_process/m_global_parameters.fpp
git commit -m "feat(params): default-initialize non-Newtonian fields in all targets"
```

### Task 3: MPI-broadcast the new fields in all three targets

**Files:**
- Modify: `src/pre_process/m_mpi_proxy.fpp`, `src/simulation/m_mpi_proxy.fpp`, `src/post_process/m_mpi_proxy.fpp`

- [ ] **Step 1: Extend the fluid_pp broadcast loop in each file**

Find the `#:for VAR in [ ... ]` Fypp loop that broadcasts `gamma, pi_inf, G, cv, qv, qvp` for `fluid_pp(i)`. Add the new real scalars to that list:

```fortran
        #:for VAR in [ 'gamma','pi_inf','G','cv','qv','qvp','K','nn','tau0','hb_m','mu_min','mu_max','mu_bulk' ]
            call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
```

Then broadcast the logical separately (after the `#:endfor`, inside the same `do i` loop):

```fortran
        call MPI_BCAST(fluid_pp(i)%non_newtonian, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
```

(Match the existing pattern in each file; the simulation copy is near line 173. `mpi_p` matches `wp`; the logical uses `MPI_LOGICAL`.)

- [ ] **Step 2: Build with MPI**

Run: `./mfc.sh build -j 8`
Expected: compiles. (MPI is on by default.)

- [ ] **Step 3: Commit**

```bash
./mfc.sh format -j 8
git add src/pre_process/m_mpi_proxy.fpp src/simulation/m_mpi_proxy.fpp src/post_process/m_mpi_proxy.fpp
git commit -m "feat(params): MPI-broadcast non-Newtonian fields in all targets"
```

### Task 4: Create the `m_hb_function` module (formula + shear-rate + mixture helper)

**Files:**
- Create: `src/simulation/m_hb_function.fpp`
- Modify: `docs/module_categories.json`

- [ ] **Step 1: Write the module**

Create `src/simulation/m_hb_function.fpp`:

```fortran
!>
!! @file m_hb_function.f90
!! @brief Herschel-Bulkley non-Newtonian viscosity: formula, shear rate, and mixture inverse-Re.

#:include 'macros.fpp'

module m_hb_function

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Re_size, Re_idx, sgm_eps, dflt_real, any_non_newtonian, hb_* arrays
    use m_constants          !< verysmall

    implicit none

    private; public :: f_compute_hb_viscosity, f_compute_shear_rate_from_components, &
                       s_compute_mixture_inv_re

contains

    !> Papanastasiou-regularized Herschel-Bulkley viscosity.
    pure function f_compute_hb_viscosity(tau0, K_val, nn_val, mu_min_val, mu_max_val, shear_rate, hb_m_val) result(mu)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: tau0, K_val, nn_val
        real(wp), intent(in) :: mu_min_val, mu_max_val
        real(wp), intent(in) :: shear_rate, hb_m_val
        real(wp)             :: mu
        real(wp)             :: yield_term, power_law_term, g_eff

        g_eff = max(shear_rate, verysmall)
        if (shear_rate <= verysmall) then
            yield_term = tau0*hb_m_val
        else
            yield_term = tau0*(1._wp - exp(-hb_m_val*shear_rate))/shear_rate
        end if
        power_law_term = K_val*(g_eff**(nn_val - 1._wp))

        mu = yield_term + power_law_term
        if (mu_min_val > dflt_real) mu = max(mu, mu_min_val)
        if (mu_max_val > dflt_real) mu = min(mu, mu_max_val)

    end function f_compute_hb_viscosity

    !> Shear rate gamma_dot = sqrt(2 D_ij D_ij). Absent dims pass 0.
    pure function f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz) result(shear_rate)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp)             :: shear_rate

        shear_rate = sqrt(2._wp*(D_xx*D_xx + D_yy*D_yy + D_zz*D_zz + 2._wp*(D_xy*D_xy + D_xz*D_xz + D_yz*D_yz)))

    end function f_compute_shear_rate_from_components

    !> Mixture inverse Reynolds (= 1/mu_mix) per direction (1=shear, 2=bulk) at one state.
    !! Reproduces the legacy Newtonian arithmetic exactly when any_non_newtonian is .false.
    !! @param alpha     volume fractions at this state
    !! @param shear_rate local shear rate at this state (used only for non-Newtonian shear term)
    !! @param Res       precomputed per-fluid Reynolds array (Res_gs or Res_viscous)
    !! @param Re_out    output (1:2) = 1/mu_mix for shear and bulk
    pure subroutine s_compute_mixture_inv_re(alpha, shear_rate, Res, Re_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(:), intent(in)   :: alpha
        real(wp), intent(in)                 :: shear_rate
        real(wp), dimension(:, :), intent(in) :: Res
        real(wp), dimension(2), intent(out)  :: Re_out
        integer  :: i, q, fl
        real(wp) :: mu_q

        ! Plain serial loops: this is already a seq device routine, so no GPU_LOOP
        ! directives (they emit empty on Cray/AMD and are redundant here).
        do i = 1, 2
            Re_out(i) = dflt_real
            if (Re_size(i) > 0) Re_out(i) = 0._wp
            do q = 1, Re_size(i)
                fl = Re_idx(i, q)
                if (any_non_newtonian .and. is_non_newtonian(fl)) then
                    if (i == 1) then
                        mu_q = f_compute_hb_viscosity(hb_tau0(fl), hb_K(fl), hb_nn(fl), &
                                                      hb_mu_min(fl), hb_mu_max(fl), shear_rate, hb_m_arr(fl))
                    else
                        mu_q = 0._wp
                        if (hb_mu_bulk(fl) > dflt_real) mu_q = hb_mu_bulk(fl)
                    end if
                    Re_out(i) = alpha(fl)*mu_q + Re_out(i)
                else
                    Re_out(i) = alpha(fl)/Res(i, q) + Re_out(i)
                end if
            end do
            Re_out(i) = 1._wp/max(Re_out(i), sgm_eps)
        end do

    end subroutine s_compute_mixture_inv_re

end module m_hb_function
```

(The module-level arrays `any_non_newtonian`, `is_non_newtonian`, `hb_tau0`, `hb_K`, `hb_nn`, `hb_m_arr`, `hb_mu_min`, `hb_mu_max`, `hb_mu_bulk` are added to `m_global_parameters` in Task 5 — this module will not compile until then. Implement Task 5 immediately after.)

- [ ] **Step 2: Register the module**

In `docs/module_categories.json`, add `"m_hb_function"` to the appropriate category array (alongside `m_viscous`). Do NOT add `m_re_visc`.

- [ ] **Step 3: Commit** (build deferred to end of Task 5)

```bash
git add src/simulation/m_hb_function.fpp docs/module_categories.json
git commit -m "feat(nn): add m_hb_function (HB formula, shear rate, mixture inv-Re)"
```

### Task 5: Add the global flag + GPU device arrays and populate them at init

**Files:**
- Modify: `src/simulation/m_global_parameters.fpp` (declarations near the `Re_size`/`Re_idx` block ~lines 211-216; init in `s_initialize_global_parameters_module` near the Re_idx population ~lines 952-982; deallocate in the finalizer)

- [ ] **Step 1: Declare the flag and device arrays (module scope)**

Near the existing `Re_size`/`Re_idx` declarations add:

```fortran
    logical :: any_non_newtonian   !< .true. if any fluid is non-Newtonian
    logical, allocatable, dimension(:) :: is_non_newtonian   !< per-fluid NN flag
    real(wp), allocatable, dimension(:) :: hb_tau0, hb_K, hb_nn, hb_m_arr
    real(wp), allocatable, dimension(:) :: hb_mu_min, hb_mu_max, hb_mu_bulk
    $:GPU_DECLARE(create='[any_non_newtonian, is_non_newtonian, hb_tau0, hb_K, hb_nn, hb_m_arr, hb_mu_min, hb_mu_max, hb_mu_bulk]')
```

- [ ] **Step 2: Default the flag where other scalars are defaulted**

In `s_assign_default_values_to_user_inputs` (or the equivalent default block), add:

```fortran
        any_non_newtonian = .false.
```

- [ ] **Step 3: Allocate + populate at init**

In `s_initialize_global_parameters_module`, after `Re_idx` is populated and `viscous` is known, add:

```fortran
        @:ALLOCATE(is_non_newtonian(1:num_fluids))
        @:ALLOCATE(hb_tau0(1:num_fluids), hb_K(1:num_fluids), hb_nn(1:num_fluids), hb_m_arr(1:num_fluids))
        @:ALLOCATE(hb_mu_min(1:num_fluids), hb_mu_max(1:num_fluids), hb_mu_bulk(1:num_fluids))

        any_non_newtonian = .false.
        do i = 1, num_fluids
            is_non_newtonian(i) = fluid_pp(i)%non_newtonian
            if (is_non_newtonian(i)) any_non_newtonian = .true.
            hb_tau0(i)   = fluid_pp(i)%tau0
            hb_K(i)      = fluid_pp(i)%K
            hb_nn(i)     = fluid_pp(i)%nn
            hb_m_arr(i)  = fluid_pp(i)%hb_m
            hb_mu_min(i) = fluid_pp(i)%mu_min
            hb_mu_max(i) = fluid_pp(i)%mu_max
            hb_mu_bulk(i) = fluid_pp(i)%mu_bulk
        end do
        $:GPU_UPDATE(device='[any_non_newtonian, is_non_newtonian, hb_tau0, hb_K, hb_nn, hb_m_arr, hb_mu_min, hb_mu_max, hb_mu_bulk]')
```

(Place this so `num_fluids` and `fluid_pp` are valid. Reuse the existing `i` iterator if present.)

- [ ] **Step 4: Deallocate in the finalizer**

In `s_finalize_global_parameters_module` add matching `@:DEALLOCATE(...)` for all nine arrays.

- [ ] **Step 5: Build (this is the first build of m_hb_function)**

Run: `./mfc.sh build -j 8`
Expected: all three targets compile. (`m_hb_function` is simulation-only; ensure it's only `use`d in simulation files. If pre/post try to use it, they will fail — they should not.)

- [ ] **Step 6: Confirm Newtonian invariance is still trivially intact**

Run: `./mfc.sh test -j 8 --only Viscous`
Expected: all pass, identical to Task 0 baseline (nothing reads the new arrays yet).

- [ ] **Step 7: Commit**

```bash
./mfc.sh format -j 8
git add src/simulation/m_global_parameters.fpp
git commit -m "feat(nn): add any_non_newtonian flag and GPU HB parameter arrays"
```

---

## DESIGN CORRECTION (2026-06-08): Riemann injection point

During Task 6 it was found that HLL (`riemann_solver=1`) and HLLC (`=2`) do NOT apply
viscosity at the per-fluid `Re_L`/`Re_R` sites. They store a harmonic-mean `Re_avg_rsx_vf`
and assemble the actual viscous stress in the SHARED routine
`s_compute_cartesian_viscous_source_flux` (m_riemann_solvers.fpp:4257), which already builds
the interface strain tensor `vel_grad_avg` with correct PHYSICAL indexing (`idx_right_phys`)
and reads `Re_avg_rsx_vf`. That shared routine is the correct, low-risk injection point.

**Decision (user-approved): non-Newtonian supports HLL + HLLC only.** LF (`=5`) and HLLD
(`=4`) are guarded off (`non_newtonian` requires `riemann_solver in {1,2}`), added in Task 11.

The Re_L/Re_R sites (320, 1004, 1438, 2899) are LEFT UNCHANGED. Tasks 6 and 8 below are
superseded by the corrected versions in this section.

### Task 6 (CORRECTED): Inject shear-dependent viscosity in the shared Cartesian flux routine

**Files:** Modify `src/simulation/m_riemann_solvers.fpp` only.

Thread `q_prim_vf` (full `dimension(sys_size)`) into the viscous-flux dispatcher and the
Cartesian routine so it can read per-fluid `alpha` at the interface:

1. `use m_hb_function` at module top.
2. `s_compute_viscous_source_flux` (line 102): add `q_prim_vf` as an argument
   (`type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf`); pass it on to
   `s_compute_cartesian_viscous_source_flux`. (Leave the cylindrical call as-is; Task 8
   extends it.)
3. `s_compute_cartesian_viscous_source_flux` (line 4257): add `q_prim_vf` argument.
4. Both call sites (HLL ~778, HLLC ~3273): pass `q_prim_vf` (the full field — it is in
   scope at both; the solver receives it).
5. Inside the loop of `s_compute_cartesian_viscous_source_flux`, AFTER `vel_grad_avg` is
   built (line 4320) and BEFORE the `norm_dir` Re_shear/Re_bulk assignment (4328+), add:

```fortran
                    if (any_non_newtonian) then
                        D_xx = vel_grad_avg(1, 1); D_yy = 0._wp; D_zz = 0._wp
                        D_xy = 0._wp; D_xz = 0._wp; D_yz = 0._wp
                        if (num_dims > 1) then
                            D_yy = vel_grad_avg(2, 2)
                            D_xy = 0.5_wp*(vel_grad_avg(1, 2) + vel_grad_avg(2, 1))
                        end if
                        if (num_dims > 2) then
                            D_zz = vel_grad_avg(3, 3)
                            D_xz = 0.5_wp*(vel_grad_avg(1, 3) + vel_grad_avg(3, 1))
                            D_yz = 0.5_wp*(vel_grad_avg(2, 3) + vel_grad_avg(3, 2))
                        end if
                        gamma_dot = f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
                        do fl = 1, num_fluids
                            alpha_avg(fl) = 0.5_wp*(q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(j_loop, k_loop, l_loop) + &
                                                    q_prim_vf(eqn_idx%adv%beg + fl - 1)%sf(idx_right_phys(1), &
                                                    idx_right_phys(2), idx_right_phys(3)))
                        end do
                        call s_compute_mixture_inv_re(alpha_avg, gamma_dot, Res_gs, Re_nn)
                    end if
```

   The `vel_grad_avg(2,2)` etc. reads are guarded because `vel_grad_avg` is
   `dimension(num_dims,num_dims)` in the normal build (out-of-bounds in 1D otherwise).

6. Replace EACH of the three `norm_dir` Re_shear/Re_bulk assignment blocks (4328-4346) so
   non-Newtonian uses the freshly computed mixture:

```fortran
                    if (any_non_newtonian) then
                        Re_shear = Re_nn(1)
                        Re_bulk = Re_nn(2)
                    else
                        Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                        Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                    end if
```
   (The `vel_src_at_interface` assignment in each norm_dir branch is unchanged.)

7. Declarations in `s_compute_cartesian_viscous_source_flux`: add
   `real(wp) :: gamma_dot, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz`,
   `real(wp), dimension(2) :: Re_nn`, `real(wp), dimension(num_fluids) :: alpha_avg`,
   `integer :: fl`. Add all of these to the `$:GPU_PARALLEL_LOOP(... private='[...]')` list.

8. Build (`-t simulation`), then the Newtonian invariance gate
   (`--only Viscous`, all pass, zero golden changes), then format + commit
   `feat(nn): shear-dependent viscosity in shared Cartesian viscous flux (HLL/HLLC)`.

Newtonian invariance holds by construction: when `any_non_newtonian=.false.`, the new block
is skipped and `Re_shear/Re_bulk = Re_avg_rsx_vf` exactly as before.

---

### Task 6 (ORIGINAL — SUPERSEDED, do not implement): Substitute mixture inv-Re in the Riemann solvers

**Files:**
- Modify: `src/simulation/m_riemann_solvers.fpp` — 4 sites matching `alpha_L(Re_idx(i, q))/Res_gs(i, q)` (HLL, HLLC, HLLD and variants; find them with the grep below). Add `use m_hb_function` at module top.

- [ ] **Step 1: Add the module use**

At the top of `m_riemann_solvers.fpp` module (with the other `use` statements):

```fortran
    use m_hb_function
```

- [ ] **Step 2: Locate the 4 substitution sites**

Run: `grep -n "alpha_L(Re_idx(i, q))/Res_gs(i, q)" src/simulation/m_riemann_solvers.fpp`
Expected: 4 line numbers. Each sits inside a block of the exact shape (lines shown for the HLLC copy near 1428):

```fortran
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 2
                            Re_L(i) = dflt_real
                            Re_R(i) = dflt_real
                            if (Re_size(i) > 0) Re_L(i) = 0._wp
                            if (Re_size(i) > 0) Re_R(i) = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, Re_size(i)
                                Re_L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                Re_R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re_R(i)
                            end do
                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                        end do
```

- [ ] **Step 3: Replace each of the 4 blocks with the helper call**

Replace the entire `do i = 1, 2 ... end do` block above with:

```fortran
                        gamma_dot_L = f_compute_shear_rate_from_components( &
                            dqL_prim_dx_vf(eqn_idx%mom%beg)%sf(j, k, l), &
                            dqL_prim_dy_vf(eqn_idx%mom%beg + 1)%sf(j, k, l), &
                            dqL_prim_dz_vf(eqn_idx%mom%beg + 2)%sf(j, k, l), &
                            0.5_wp*(dqL_prim_dy_vf(eqn_idx%mom%beg)%sf(j, k, l) + dqL_prim_dx_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)), &
                            0.5_wp*(dqL_prim_dz_vf(eqn_idx%mom%beg)%sf(j, k, l) + dqL_prim_dx_vf(eqn_idx%mom%beg + 2)%sf(j, k, l)), &
                            0.5_wp*(dqL_prim_dz_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) + dqL_prim_dy_vf(eqn_idx%mom%beg + 2)%sf(j, k, l)))
                        gamma_dot_R = f_compute_shear_rate_from_components( &
                            dqR_prim_dx_vf(eqn_idx%mom%beg)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)), &
                            dqR_prim_dy_vf(eqn_idx%mom%beg + 1)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)), &
                            dqR_prim_dz_vf(eqn_idx%mom%beg + 2)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)), &
                            0.5_wp*(dqR_prim_dy_vf(eqn_idx%mom%beg)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)) + dqR_prim_dx_vf(eqn_idx%mom%beg + 1)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3))), &
                            0.5_wp*(dqR_prim_dz_vf(eqn_idx%mom%beg)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)) + dqR_prim_dx_vf(eqn_idx%mom%beg + 2)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3))), &
                            0.5_wp*(dqR_prim_dz_vf(eqn_idx%mom%beg + 1)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)) + dqR_prim_dy_vf(eqn_idx%mom%beg + 2)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3))))
                        call s_compute_mixture_inv_re(alpha_L, gamma_dot_L, Res_gs, Re_L)
                        call s_compute_mixture_inv_re(alpha_R, gamma_dot_R, Res_gs, Re_R)
```

**Important guards for absent dimensions and case-optimization:** the `dq*_dy`/`dq*_dz` accesses must be compiled out in 1D/2D. Wrap the `mom%beg+1`/`mom%beg+2` gradient reads in the same `#:if not MFC_CASE_OPTIMIZATION or num_dims > 1` / `> 2` Fypp guards used at lines 1452-1465, substituting 0 for the missing component. Concretely, build `gamma_dot_L`/`gamma_dot_R` via intermediate locals `dl_* = 0._wp` set under the dimension guards, then call `f_compute_shear_rate_from_components(dl_xx, dl_yy, ...)`. (This keeps 1D/2D builds valid; the helper already treats absent components as 0.)

- [ ] **Step 4: Declare the new scalars in each solver subroutine**

Add `gamma_dot_L, gamma_dot_R` (and any `dl_*` intermediates) to the local `real(wp)` declarations of each solver subroutine, and add them to the `private='[...]'` list of the enclosing `GPU_PARALLEL_LOOP`.

- [ ] **Step 5: Build (CPU)**

Run: `./mfc.sh build -j 8 -t simulation`
Expected: compiles in 1D/2D/3D configurations. If nvfortran present: `./mfc.sh build -j 8 -t simulation --gpu acc`.

- [ ] **Step 6: Newtonian invariance check (CRITICAL)**

Run: `./mfc.sh test -j 8 --only Viscous`
Expected: ALL pass with NO golden changes. The helper reproduces the legacy arithmetic when `any_non_newtonian=.false.`. If any golden differs, the helper is not bitwise-equivalent — STOP, diff the helper arithmetic against the original loop, fix before proceeding. Do NOT `--generate`.

- [ ] **Step 7: Commit**

```bash
./mfc.sh format -j 8
git add src/simulation/m_riemann_solvers.fpp
git commit -m "feat(nn): spatially varying viscosity in Cartesian Riemann viscous flux"
```

### Task 7: Power-law Poiseuille example with analytic comparison (the correctness anchor)

**Files:**
- Create: `examples/2D_poiseuille_nn/case.py`

- [ ] **Step 1: Write the case**

Create `examples/2D_poiseuille_nn/case.py` as a 2D channel (periodic in x, no-slip walls in y) driven by a constant body force `g_x`, single power-law fluid (`tau0=0`, `nn != 1`), with `mu_max`/`mu_min` set. Choose nondimensional `K`, `nn`, force, and half-height `H` so the steady analytic profile is well resolved. The power-law Poiseuille steady velocity profile (force per unit mass `g`, density `rho`, half-height `H`, wall at `y=±H`) is:

```
u(y) = (n / (n+1)) * (rho * g / K)^(1/n) * ( H^((n+1)/n) - |y|^((n+1)/n) )
```

Include a trailing Python block (run only when invoked directly, not during MFC's import of the dict) that, given the simulation output, computes the L2 norm of `(u_numeric - u_analytic)/u_max` across the centerline column and prints it. Keep the grid modest (e.g. `m=63, n=127`) and integrate to steady state.

- [ ] **Step 2: Validate the case file parses and passes the validator**

Run: `./mfc.sh validate examples/2D_poiseuille_nn/case.py`
Expected: passes (after Task 10's validator additions; if run before Task 10, it should still pass the base validator).

- [ ] **Step 3: Run and check the analytic error**

Run: `./mfc.sh run examples/2D_poiseuille_nn/case.py -n 2`
Then run the comparison block. Expected: relative L2 error below a documented tolerance (target `< 2e-2` on this grid; record the actual value in the case docstring). If the error is large or the profile is qualitatively wrong (e.g. mirrored, wrong curvature for `n<1` vs `n>1`), STOP — this indicates a sign/index error in Task 6; debug before continuing. **This is the gate that proves the physics is right before any golden is frozen.**

- [ ] **Step 4: Commit**

```bash
./mfc.sh format -j 8
git add examples/2D_poiseuille_nn/case.py
git commit -m "test(nn): power-law Poiseuille example with analytic velocity-profile check"
```

---

## Phase 3: Cylindrical viscous path

### Task 8: Substitute mixture inv-Re in `m_viscous` cylindrical boundary

**Files:**
- Modify: `src/simulation/m_viscous.fpp` — 4 sites matching `alpha_visc(Re_idx(i, q))/Res_viscous(i, q)` (lines 152-160, 255-263, 352-360, 452-460). Add `use m_hb_function`.

- [ ] **Step 1: Add the module use** at the top of `m_viscous.fpp`:

```fortran
    use m_hb_function
```

- [ ] **Step 2: Locate the sites**

Run: `grep -n "alpha_visc(Re_idx(i, q))/Res_viscous(i, q)" src/simulation/m_viscous.fpp`
Expected: 4 line numbers, each in a block:

```fortran
                                        Re_visc(i) = dflt_real
                                        if (Re_size(i) > 0) Re_visc(i) = 0._wp
                                        do q = 1, Re_size(i)
                                            Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) + Re_visc(i)
                                        end do
                                        Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
```

- [ ] **Step 3: Replace each block**

These blocks live inside loops over `i = 1, 2` with cell indices available as `(j, k, l)` and gradient fields `grad_x_vf, grad_y_vf, grad_z_vf` (args of `s_compute_viscous_stress_cylindrical_boundary`). Compute the shear rate from those gradient fields at the local cell and call the helper. Replace the `do i = 1, 2 ... end do` structure (the four occurrences each have their own enclosing `do i`) with:

```fortran
                                    gamma_dot_c = f_compute_shear_rate_from_components( &
                                        grad_x_vf(1)%sf(j, k, l), grad_y_vf(2)%sf(j, k, l), grad_z_vf(3)%sf(j, k, l), &
                                        0.5_wp*(grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l)), &
                                        0.5_wp*(grad_z_vf(1)%sf(j, k, l) + grad_x_vf(3)%sf(j, k, l)), &
                                        0.5_wp*(grad_z_vf(2)%sf(j, k, l) + grad_y_vf(3)%sf(j, k, l)))
                                    call s_compute_mixture_inv_re(alpha_visc, gamma_dot_c, Res_viscous, Re_visc)
```

**Verify the gradient field layout first:** confirm `grad_x_vf(c)%sf` is `d(vel_c)/dx` (i.e. index `c` = velocity component, derivative direction = the `_x/_y/_z` of the field name). Read the call site of `s_compute_viscous_stress_cylindrical_boundary` and the start of `s_get_viscous` to confirm before applying. If the component/direction convention differs, adjust the `f_compute_shear_rate_from_components` argument mapping accordingly. Guard `grad_y_vf`/`grad_z_vf` accesses under the existing `n>0`/`p>0` conditions (cylindrical implies `n>0`); pass 0 for absent dims.

- [ ] **Step 4: Declare `gamma_dot_c`** as a local `real(wp)` in `s_compute_viscous_stress_cylindrical_boundary` and add to the relevant GPU loop `private` list if inside one.

- [ ] **Step 5: Build**

Run: `./mfc.sh build -j 8 -t simulation`
Expected: compiles.

- [ ] **Step 6: Newtonian invariance check**

Run: `./mfc.sh test -j 8 --only Viscous`
Expected: all pass, no golden change (cylindrical viscous tests included). Do NOT `--generate`.

- [ ] **Step 7: Commit**

```bash
./mfc.sh format -j 8
git add src/simulation/m_viscous.fpp
git commit -m "feat(nn): spatially varying viscosity in cylindrical viscous stress"
```

---

## Phase 4: IBM

### Task 9: Evaluate per-sample viscosity inside the stress-tensor routine

**Files:**
- Modify: `src/simulation/m_viscous.fpp:1265-1327` (`s_compute_viscous_stress_tensor`)

- [ ] **Step 1: Compute shear rate and per-cell mixture viscosity inside the routine**

`s_compute_viscous_stress_tensor` already builds `velocity_gradient_tensor` and has `(i,j,k)` and `q_prim_vf`. After the gradient tensor is built (after line 1300) and before the stress assembly (line 1308), insert:

```fortran
        ! Non-Newtonian: override mu with the local strain-rate-dependent mixture viscosity,
        ! so each stencil sample (i,j,k) uses its own viscosity.
        mu_eff = dynamic_viscosity
        if (any_non_newtonian) then
            gamma_dot_c = f_compute_shear_rate_from_components( &
                velocity_gradient_tensor(1, 1), velocity_gradient_tensor(2, 2), velocity_gradient_tensor(3, 3), &
                0.5_wp*(velocity_gradient_tensor(1, 2) + velocity_gradient_tensor(2, 1)), &
                0.5_wp*(velocity_gradient_tensor(1, 3) + velocity_gradient_tensor(3, 1)), &
                0.5_wp*(velocity_gradient_tensor(2, 3) + velocity_gradient_tensor(3, 2)))
            mu_eff = 0._wp
            do l = 1, num_fluids
                if (is_non_newtonian(l)) then
                    mu_eff = mu_eff + q_prim_vf(eqn_idx%adv%beg + l - 1)%sf(i, j, k)* &
                             f_compute_hb_viscosity(hb_tau0(l), hb_K(l), hb_nn(l), hb_mu_min(l), hb_mu_max(l), gamma_dot_c, hb_m_arr(l))
                else
                    mu_eff = mu_eff + q_prim_vf(eqn_idx%adv%beg + l - 1)%sf(i, j, k)/max(fluid_pp(l)%Re(1), sgm_eps)
                end if
            end do
        end if
```

Then replace `dynamic_viscosity` with `mu_eff` in the stress assembly (lines 1311 and 1317).

(`velocity_gradient_tensor(c, d)` is `d(vel_c)/dx_d` per lines 1292-1298, so the mapping above is correct: `(1,1)=du/dx`, `(2,2)=dv/dy`, `(1,2)=du/dy`, etc.)

- [ ] **Step 2: Declare locals and confirm `m_hb_function` use**

Add `real(wp) :: mu_eff, gamma_dot_c` to the routine's declarations. `m_viscous` already `use`s `m_hb_function` from Task 8.

- [ ] **Step 3: Build**

Run: `./mfc.sh build -j 8 -t simulation`
Expected: compiles. Note `s_compute_viscous_stress_tensor` is `GPU_ROUTINE seq`; the `do l = 1, num_fluids` loop and `fluid_pp` read — prefer `hb_*` arrays; `fluid_pp(l)%Re(1)` is read-only and already used on-device elsewhere in IBM (line 929-930), so acceptable, but if a GPU build fails on the derived-type read, pre-extract `fluid_pp%Re(1)` into a `GPU_DECLARE`'d `inv_re_newt(num_fluids)` array in Task 5 and use it here.

- [ ] **Step 4: Newtonian invariance check**

Run: `./mfc.sh test -j 8 --only IBM` and `./mfc.sh test -j 8 --only Viscous`
Expected: all pass, no golden change (`any_non_newtonian=.false.` ⇒ `mu_eff = dynamic_viscosity`). Do NOT `--generate`.

- [ ] **Step 5: Commit**

```bash
./mfc.sh format -j 8
git add src/simulation/m_viscous.fpp
git commit -m "feat(nn): per-stencil-sample viscosity in IBM viscous stress tensor"
```

---

## Phase 5: CFL bounding viscosity

### Task 10: Use `mu_max` as the bounding viscosity in the viscous CFL limit

**Files:**
- Modify: `src/simulation/m_sim_helpers.fpp` (`s_compute_dt_from_cfl`, the `Re_l` viscous-CFL block ~lines 191-236)
- Modify: `src/simulation/m_time_steppers.fpp` (where `Re_l` is assembled before calling `s_compute_dt_from_cfl`, ~line 658)

- [ ] **Step 1: Inspect how `Re_l(1:2)` is built before the CFL call**

Run: `grep -n "Re_l\|s_compute_dt_from_cfl" src/simulation/m_time_steppers.fpp`
Read that block. It builds `Re_l` as the per-cell mixture `1/mu` the same `alpha/Res` way.

- [ ] **Step 2: Bound the shear viscosity by `mu_max` for non-Newtonian fluids**

Where `Re_l(1)` (shear) is assembled for the CFL, replace the non-Newtonian fluids' contribution with their `mu_max` bound. Concretely, build the mixture inverse-Re using `mu_max` for NN fluids:

```fortran
            if (any_non_newtonian) then
                Re_l(1) = 0._wp
                do q = 1, Re_size(1)
                    fl = Re_idx(1, q)
                    if (is_non_newtonian(fl)) then
                        Re_l(1) = alpha(fl)*hb_mu_max(fl) + Re_l(1)
                    else
                        Re_l(1) = alpha(fl)/Res(1, q) + Re_l(1)
                    end if
                end do
                Re_l(1) = 1._wp/max(Re_l(1), sgm_eps)
            else
                ! ... existing Newtonian assembly unchanged ...
            end if
```

Since `mu(gamma_dot) <= mu_max`, `Re_l(1) >= 1/mu_max`, so `1/(rho*Re_l)` (the kinematic viscosity used in `vcfl_dt`) is bounded above — guaranteeing a stable (possibly slightly conservative) dt. Bulk uses `hb_mu_bulk` analogously if set. Reuse the local `Res`/`alpha` names already present at that call site; add `fl` to declarations.

- [ ] **Step 3: Build**

Run: `./mfc.sh build -j 8 -t simulation`
Expected: compiles.

- [ ] **Step 4: Newtonian invariance check**

Run: `./mfc.sh test -j 8 --only Viscous`
Expected: all pass, no golden change (`any_non_newtonian=.false.` branch is the untouched assembly).

- [ ] **Step 5: Commit**

```bash
./mfc.sh format -j 8
git add src/simulation/m_sim_helpers.fpp src/simulation/m_time_steppers.fpp
git commit -m "feat(nn): bound viscous CFL by mu_max for non-Newtonian fluids"
```

---

## Phase 6: Guards, examples, tests, docs

### Task 11: Add case-validator guards (Python)

**Files:**
- Modify: `toolchain/mfc/case_validator.py`

- [ ] **Step 1: Add a `check_non_newtonian` method**

Following the existing check-method pattern (and the `PHYSICS_DOCS` registration style), add per-fluid checks. For each fluid `i` with `fluid_pp(i)%non_newtonian == 'T'`:

```python
        for i in range(1, num_fluids + 1):
            if self.get(f"fluid_pp({i})%non_newtonian", "F") != "T":
                continue
            self.prohibit(self.get("viscous", "F") != "T",
                          f"fluid_pp({i})%non_newtonian requires viscous = T")
            self.prohibit(self.is_default(f"fluid_pp({i})%K"),
                          f"fluid_pp({i})%non_newtonian requires K")
            self.prohibit(self.is_default(f"fluid_pp({i})%nn"),
                          f"fluid_pp({i})%non_newtonian requires nn")
            self.prohibit(self.is_default(f"fluid_pp({i})%mu_max"),
                          f"fluid_pp({i})%non_newtonian requires mu_max")
            tau0 = self.get(f"fluid_pp({i})%tau0", 0.0)
            self.prohibit(float(tau0) > 0.0 and self.is_default(f"fluid_pp({i})%hb_m"),
                          f"fluid_pp({i})%tau0 > 0 requires hb_m")
            mu_min = self.get(f"fluid_pp({i})%mu_min", None)
            mu_max = self.get(f"fluid_pp({i})%mu_max", None)
            if mu_min is not None and mu_max is not None and not self.is_default(f"fluid_pp({i})%mu_min"):
                self.prohibit(float(mu_max) <= float(mu_min),
                              f"fluid_pp({i})%mu_max must exceed mu_min")
            self.prohibit(self.get("igr", "F") == "T",
                          f"fluid_pp({i})%non_newtonian is incompatible with igr")
            self.prohibit(int(self.get("model_eqns", 0)) not in (2, 3),
                          f"fluid_pp({i})%non_newtonian requires model_eqns 2 or 3")
```

Match the exact helper names this file uses (`self.get`, `self.prohibit`, a default-sentinel check — read the file to confirm whether it's `self.is_default(...)` or a comparison to `dflt_real`/`None`, and adapt). Register `check_non_newtonian` wherever the other `check_*` methods are dispatched, with a `PHYSICS_DOCS` entry pointing at the case.md section.

- [ ] **Step 2: Confirm both example cases validate**

Run: `./mfc.sh validate examples/2D_poiseuille_nn/case.py`
Expected: passes.
Run: a deliberately-broken probe — temporarily set `mu_max` to default in a scratch copy and confirm `./mfc.sh validate` now errors with "requires mu_max". Discard the scratch copy.

- [ ] **Step 3: Commit**

```bash
git add toolchain/mfc/case_validator.py
git commit -m "feat(nn): case-validator guards for non-Newtonian parameters"
```

### Task 12: Lid-driven cavity example (qualitative, Li et al. 2015)

**Files:**
- Create: `examples/2D_lid_driven_cavity_nn/case.py`

- [ ] **Step 1: Write the case** — a 2D lid-driven cavity with a single non-Newtonian fluid, parameterized so one run is shear-thinning (`nn<1`) and the documented configuration matches the Li et al. (2015) `Re=500` setup. Set `mu_max`/`mu_min`. Put the literature reference and expected qualitative behavior in the docstring.

- [ ] **Step 2: Validate and smoke-run**

Run: `./mfc.sh validate examples/2D_lid_driven_cavity_nn/case.py`
Run: `./mfc.sh run examples/2D_lid_driven_cavity_nn/case.py -n 2 -t pre_process simulation` for a few steps to confirm it runs without NaN.
Expected: runs; centerline velocity qualitatively matches the shear-thinning trend.

- [ ] **Step 3: Commit**

```bash
./mfc.sh format -j 8
git add examples/2D_lid_driven_cavity_nn/case.py
git commit -m "test(nn): lid-driven cavity non-Newtonian example (Li et al. 2015)"
```

### Task 13: Regression tests + golden generation

**Files:**
- Modify: `toolchain/mfc/test/cases.py`

- [ ] **Step 1: Add a Non-Newtonian test block**

Following the `stack.push("Viscous", {...})` pattern, add a small fast block on a tiny grid. Include a shear-thinning (`nn=0.5`) and a shear-thickening (`nn=1.5`) power-law variant (`tau0=0`), each with `mu_max`/`mu_min` set and `viscous='T'`, exercised in 1D and 2D via the existing `dimInfo` mechanism:

```python
    stack.push("Non-Newtonian", {
        "viscous": "T",
        "fluid_pp(1)%Re(1)": 1.0e4,
        "fluid_pp(1)%non_newtonian": "T",
        "fluid_pp(1)%tau0": 0.0,
        "fluid_pp(1)%K": 1.0e-4,
        "fluid_pp(1)%nn": 0.5,
        "fluid_pp(1)%mu_max": 1.0e-1,
        "fluid_pp(1)%mu_min": 1.0e-6,
        "fluid_pp(1)%hb_m": 1.0e3,
        "dt": 1e-8,
        "patch_icpp(1)%vel(1)": 1.0,
    })
    cases.append(define_case_d(stack, "nn=0.5", {}))
    cases.append(define_case_d(stack, "nn=1.5", {"fluid_pp(1)%nn": 1.5}))
    stack.pop()
```

(Adjust keys to whatever `BASE_CFG` requires for a minimal viscous run; mirror the neighboring Viscous block's grid/patch setup. Keep runtime small.)

- [ ] **Step 2: List the new tests**

Run: `./mfc.sh test -l | grep -i "non-newtonian"`
Expected: the new cases appear with fresh UUIDs.

- [ ] **Step 3: Generate goldens (ONLY after the Task 7 analytic check passed)**

Run: `./mfc.sh test --generate --only Non-Newtonian -j 8`
Expected: golden files created under `tests/<UUID>/`. Spot-check a golden for sane (non-NaN, non-constant) field values.

- [ ] **Step 4: Verify the goldens reproduce**

Run: `./mfc.sh test --only Non-Newtonian -j 8`
Expected: all pass.

- [ ] **Step 5: Final Newtonian-invariance check across the whole suite**

Run: `./mfc.sh test -j 8 --only Viscous` and `./mfc.sh test -j 8 --only IBM`
Expected: all pass; `git status` shows NO modifications under existing `tests/<UUID>/` dirs (only new ones added).

- [ ] **Step 6: Commit**

```bash
git add toolchain/mfc/test/cases.py tests/
git commit -m "test(nn): regression cases for shear-thinning and shear-thickening fluids"
```

### Task 14: Documentation

**Files:**
- Modify: `docs/documentation/case.md`

- [ ] **Step 1: Add the non-Newtonian section** under the fluid parameters table: the parameter table (8 params), the HB formula, the special cases (power-law / Bingham / Newtonian), the nondimensionalization note, and a pointer that `Re(1)` is the reference Reynolds number that registers the fluid as viscous. (Reuse the wording from PR #1298's `case.md` addition, which was accurate.)

- [ ] **Step 2: Run docs lint**

Run: `./mfc.sh precheck -j 8`
Expected: `lint_docs.py` passes (docs freshness satisfied).

- [ ] **Step 3: Commit**

```bash
git add docs/documentation/case.md
git commit -m "docs(nn): document Herschel-Bulkley non-Newtonian viscosity"
```

---

## Phase 7: Full verification

### Task 15: Whole-suite verification and push

- [ ] **Step 1: Format + precheck**

Run: `./mfc.sh format -j 8 && ./mfc.sh precheck -j 8`
Expected: clean.

- [ ] **Step 2: Full build, all targets, CPU**

Run: `./mfc.sh build -j 8`
Expected: pre_process, simulation, post_process all compile.

- [ ] **Step 3: GPU build if available**

Run (if nvfortran present): `./mfc.sh build -j 8 -t simulation --gpu acc`
Expected: compiles. (Cray/AMD/Intel are CI-gated.)

- [ ] **Step 4: Targeted test pass**

Run: `./mfc.sh test -j 8 --only Viscous` ; `--only IBM` ; `--only Non-Newtonian`
Expected: all pass; no diffs to pre-existing goldens.

- [ ] **Step 5: Push the branch**

```bash
git push fork feature/non-newtonian-viscosity
```

- [ ] **Step 6: Open the PR** (when ready) against `MFlowCode/MFC` from `sbryngelson:feature/non-newtonian-viscosity`, referencing the closed #1298 and this plan. (Defer to user.)

---

## Notes / risks flagged for the implementer

- **Newtonian invariance is the #1 acceptance gate.** Every phase re-runs `--only Viscous` without `--generate`. A single golden diff before Phase 6 means the helper is not bitwise-equivalent — fix, don't regenerate.
- **`s_compute_mixture_inv_re` arithmetic order** must match the original inlined loop exactly (same `+ Re_out(i)` accumulation order, same `1._wp/max(..., sgm_eps)`). This is what makes the Newtonian path bitwise-identical.
- **Gradient component/direction convention** in `m_viscous` (`grad_x_vf(c)%sf`) must be confirmed by reading `s_get_viscous` before Task 8; the Riemann `dq*_d{x,y,z}_vf(mom%beg+c-1)` convention is confirmed (m_riemann_solvers.fpp:1449-1462).
- **GPU derived-type reads:** prefer the `hb_*`/`is_non_newtonian` device arrays over `fluid_pp%...` inside kernels. If Task 9's `fluid_pp(l)%Re(1)` read fails a GPU build, pre-extract it into a device array in Task 5.
- **Case-optimization:** `num_fluids`/`num_dims` are baked into the binary under `--case-optimization`; the Fypp dimension guards in Task 6 must compile out unused-dimension gradient reads.
- **Assumed-shape device dummies:** `s_compute_mixture_inv_re` takes `alpha(:)` and `Res(:,:)` assumed-shape. Assumed-shape works well for the AMD case-opt path where `alpha_L` is `dimension(3)` vs `dimension(num_fluids)`, but some GPU compilers dislike assumed-shape descriptors in `routine seq`. If the first `--gpu acc`/`--gpu mp` build rejects it, switch the dummies to explicit shape (`alpha(1:num_fluids)`, `Res(1:2, 1:Re_size_max)`) and replicate the AMD case-opt dimension handling the original Riemann block uses.
