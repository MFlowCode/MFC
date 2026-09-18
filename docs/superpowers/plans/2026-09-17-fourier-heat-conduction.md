# Fourier Heat Conduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add $\nabla\cdot(k\nabla T)$ to MFC's mixture energy equation for single- and multi-fluid Navier-Stokes runs, in Cartesian and cylindrical geometry.

**Architecture:** Conduction enters as a face-centered central difference accumulated into `flux_src_vf(eqn_idx%E)`, the same array the viscous work term uses, so the existing `flux_src` differencing carries it into the RHS unchanged. It is called from `m_rhs.fpp` beside the existing chemistry-diffusion hook, not from inside a Riemann solver, so no Riemann solver signature changes. Temperature is the thermal-equilibrium mixture temperature, computed from primitives during conservative-to-primitive conversion and stored in the already-existing `q_T_sf`.

**Tech Stack:** Fortran 2008 preprocessed by Fypp (`.fpp`), GPU offload via the `$:GPU_*` macro layer, Python toolchain (`./mfc.sh`) for builds, case generation, and golden-file testing.

**Spec:** `docs/superpowers/specs/2026-09-17-fourier-heat-conduction-design.md`

## Global Constraints

- New branches are made on forks, never on `MFlowCode/MFC`. This work is on branch `heat` of `sbryngelson/MFC`.
- A PR made with AI tools must say so, and must follow the PR template.
- A PR that changes CFD results needs verification that the change is correct. Task 3 is that verification.
- Build and run through `./mfc.sh`. It can leave a sticky lock file in `build/`; if a build refuses to start, remove `build/lock.yaml` before retrying.
- DRY, including side-effect code. Comments as short as possible without losing value. GPU macros follow the existing patterns in the file being edited. Prefer short subroutines with separated concerns.
- Real literals use the `_wp` suffix. Never write a bare `0.0` in numeric Fortran code.
- The feature must be bit-identical to master when no fluid sets `k_therm`.

## Deviations from the spec, decided during planning

Three refinements were found while reading the code. They are improvements on what the spec describes, and each task below implements the refined version.

1. **Call site.** The spec put the conduction flux call inside the HLL/HLLC solvers next to `s_compute_viscous_source_flux`. `m_rhs.fpp:719-724` already calls `s_compute_chemistry_diffusion_flux` with exactly the arguments conduction needs, including `q_T_sf`, which the Riemann solvers do not receive. Conduction is called from there instead. No Riemann solver is touched.
2. **Temperature from primitives.** The spec's mixture temperature $T=(\rho e + p - \sum\alpha_i\rho_i q_{v,i})/(\sum\alpha_i\rho_i c_{v,i} n_i)$ is algebraically identical to $T=\big((\sum_i\alpha_i\gamma_i + 1)p + \sum_i\alpha_i\pi_{\infty,i}\big)/\sum_i\alpha_i\rho_i c_{v,i}n_i$ in MFC's stored variables, because $\rho e = \gamma_{\mathrm{mix}}p + \pi_{\infty,\mathrm{mix}} + \sum\alpha_i\rho_i q_{v,i}$. The second form uses `gamma_K` and `pi_inf_K`, which `s_convert_conservative_to_primitive_variables` has already computed in the same loop iteration. Single-fluid sanity check: $\gamma+1=\Gamma/(\Gamma-1)$ and $\pi_\infty^{\mathrm{MFC}}=\Gamma\pi_\infty/(\Gamma-1)$ give $T=(p+\pi_\infty)/((\Gamma-1)\rho c_v)$, which is exactly `f_sg_thermal`.
3. **Cylindrical is almost free.** The spec expected new code for the geometric source. `m_rhs.fpp:1885-1897` and `:1913-1925` already apply `-0.5/y_cc(k)*(flux_src(k-1) + flux_src(k))` over `i = mom%beg, E` gated on `cyl_coord` alone, so the $\tfrac{1}{r}k\,\partial_rT$ source appears as soon as the energy source flux carries conduction. Only the axis cell (`bc_y%beg == -2` or `-14`), which the generic loop skips, needs new code.

4. **The RHS gates are load-bearing for Task 3, not Task 5.** `flux_src(mom%beg:E)` reaches the RHS only through guards at `m_rhs.fpp:727`, `:1769/:1774` (x), `:1854/:1859` (y), and `:1946/:1951` (z), all written as `surface_tension .or. viscous`. Until those include `heat_conduction`, the conduction flux is computed and then discarded. They are widened in Task 3, where the verification catches their absence.

5. **No phase-change refactor.** The spec asked for one shared temperature function with `m_phase_change.fpp` switched over. `s_infinite_pt_relaxation_k` accumulates `mCP` and `mQ` in the same sequential loop that feeds its Newton solver's `gp`/`gpp`; routing it through a shared helper would split that loop and perturb solver code this change has no business touching. Instead the new function is added for the primitive-variable path and both sites carry a one-line comment noting the forms are equivalent. Flagged rather than done silently.

## File Structure

**Created**
- `src/simulation/m_conduction.fpp` — the conduction source flux. Only physics: face conductivity, face temperature gradient, accumulation into the energy source flux. Cartesian and cylindrical.
- `examples/1D_conduction_convergence/case.py` — sinusoidal-temperature, uniform-pressure, quiescent 1D case.
- `examples/1D_conduction_convergence/compare_analytic.py` — one-step consistency and spatial convergence check against the exact Laplacian.
- `examples/1D_conduction_convergence/README.md` — what the case verifies and how to run it.

**Modified**
- `src/common/m_derived_types.fpp` — `k_therm` member on `physical_parameters`.
- `src/common/m_global_parameters_common.fpp` — `heat_conduction` flag and `fluid_k_therm` array, shared by all three executables.
- `src/common/m_variables_conversion.fpp` — fill `fluid_k_therm`/`heat_conduction` at init, free them at finalize, add `f_mixture_temperature`, populate `q_T_sf`.
- `src/{simulation,pre_process,post_process}/m_global_parameters.fpp` — `k_therm` default.
- `src/simulation/m_rhs.fpp` — allocate the energy source flux, widen the four RHS gates, call the conduction flux, cylindrical axis source.
- `src/simulation/m_riemann_state.fpp` — zero the energy source flux when conduction is on.
- `src/simulation/m_sim_helpers.fpp` — thermal CFL limit and diagnostic.
- `src/simulation/m_time_steppers.fpp` — thermal `dt` candidate.
- `src/simulation/m_data_output.fpp` — TCFL in run-time info.
- `src/common/m_mpi_common.fpp` — stability-extrema reduction widened.
- `src/simulation/m_checker.fpp` — input validation.
- `toolchain/mfc/params/definitions.py`, `toolchain/mfc/params/descriptions.py` — register `k_therm`.
- `toolchain/mfc/case_validator.py` — mirror the Fortran checks.
- `toolchain/mfc/test/cases.py` — golden tests.
- `docs/documentation/case.md`, `docs/documentation/equations.md` — user docs.

---

### Task 1: `k_therm` input and the `heat_conduction` gate

Adds the input parameter and the derived flag, wired end to end but affecting no physics. Deliverable: a case file that sets `fluid_pp(1)%k_therm` runs and produces output bit-identical to the same case without it.

**Files:**
- Modify: `src/common/m_derived_types.fpp` (the `physical_parameters` type)
- Modify: `src/common/m_global_parameters_common.fpp`
- Modify: `src/common/m_variables_conversion.fpp:268-350` (init) and `:1244-1256` (finalize)
- Modify: `src/simulation/m_global_parameters.fpp:474`, `src/pre_process/m_global_parameters.fpp:425`, `src/post_process/m_global_parameters.fpp:229`
- Modify: `toolchain/mfc/params/definitions.py:129` and `:962`
- Modify: `toolchain/mfc/params/descriptions.py:351`

**Interfaces:**
- Consumes: nothing.
- Produces: `fluid_pp(i)%k_therm` (real, default `0._wp`); module variables `heat_conduction` (logical) and `fluid_k_therm(1:num_fluids)` (real array, device-resident), both exported from `m_global_parameters_common` and therefore visible anywhere `m_global_parameters` is used.

- [ ] **Step 1: Add the type member**

In `src/common/m_derived_types.fpp`, in `type physical_parameters`, immediately after the `Re` line:

```fortran
        real(wp), dimension(2) :: Re                 !< Reynolds number
        real(wp)               :: k_therm            !< Thermal conductivity (name avoids %K, the Herschel-Bulkley index)
```

- [ ] **Step 2: Default it in all three executables**

In each of `src/simulation/m_global_parameters.fpp`, `src/pre_process/m_global_parameters.fpp`, and `src/post_process/m_global_parameters.fpp`, in the `fluid_pp` default loop, immediately after the `fluid_pp(i)%Re(:) = dflt_real` line:

```fortran
            fluid_pp(i)%k_therm = 0._wp
```

Zero, not `dflt_real`: zero conductivity is the physically meaningful "off" state and makes the gate a simple positivity test.

- [ ] **Step 3: Declare the shared flag and array**

In `src/common/m_global_parameters_common.fpp`, in the "Material properties derived from fluid_pp" block, after the `eos_coeffs` declaration:

```fortran
    !> Fourier heat conduction: true when any fluid sets k_therm > 0. Derived, never read from the namelist.
    logical                             :: heat_conduction
    real(wp), allocatable, dimension(:) :: fluid_k_therm
    $:GPU_DECLARE(create='[heat_conduction, fluid_k_therm]')
```

- [ ] **Step 4: Fill them at initialization**

In `src/common/m_variables_conversion.fpp`, in `s_initialize_variables_conversion_module`, add to the allocation block next to `@:ALLOCATE(cvs (1:num_fluids))`:

```fortran
        @:ALLOCATE(fluid_k_therm(1:num_fluids))
```

Inside the `do i = 1, num_fluids` loop, next to `cvs(i) = fluid_pp(i)%cv`:

```fortran
            fluid_k_therm(i) = fluid_pp(i)%k_therm
```

Immediately after that loop ends, before the existing `$:GPU_UPDATE(device='[gammas, ...]')`:

```fortran
        heat_conduction = any(fluid_k_therm > 0._wp)
```

Then extend that same `GPU_UPDATE` list to carry the two new names:

```fortran
        $:GPU_UPDATE(device='[gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps, Gs_vc, eoss, eos_coeffs, &
                     & fluid_k_therm, heat_conduction]')
```

- [ ] **Step 5: Free them at finalization**

In `s_finalize_variables_conversion_module`, extend the existing deallocate:

```fortran
        @:DEALLOCATE(gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps, Gs_vc, eoss, fluid_k_therm)
```

- [ ] **Step 6: Register the parameter in the toolchain**

In `toolchain/mfc/params/definitions.py`, add the tag group next to the existing `"viscosity"` entry at line 129:

```python
    "heat_conduction": "Heat conduction",
```

and register the parameter in the `fluid_pp` loop, immediately after the two `Re` registrations:

```python
        _r(f"{px}k_therm", REAL, {"heat_conduction"}, math=r"\f$k_k\f$")
```

`fluid_pp` is already a namelist root for all three targets, so `NAMELIST_VARS` needs no change.

In `toolchain/mfc/params/descriptions.py`, add to the `fluid_pp patterns` list after the `cv` entry:

```python
    (r"fluid_pp\((\d+)\)%k_therm", "Thermal conductivity for fluid {0}"),
```

- [ ] **Step 7: Build**

```bash
./mfc.sh build -t pre_process simulation post_process -j 8
```

Expected: a clean build. The cmake reconfigure regenerates `generated_decls.fpp` with the new member. If the build refuses to start, delete `build/lock.yaml` and retry.

- [ ] **Step 8: Verify the parameter is accepted and inert**

```bash
./mfc.sh run examples/1D_vacuum/case.py -n 1
cp -r examples/1D_vacuum/restart_data /tmp/conduction_baseline
```

Then add `"fluid_pp(1)%k_therm": 1.0e-3,` to that case's parameter dict, rerun, and diff:

```bash
./mfc.sh run examples/1D_vacuum/case.py -n 1
diff -r /tmp/conduction_baseline examples/1D_vacuum/restart_data
```

Expected: the run succeeds (the parameter is recognized, not rejected as unknown) and the diff is empty. Revert the case file edit afterwards.

- [ ] **Step 9: Confirm no regression**

```bash
./mfc.sh test -j 8 --percent 20
```

Expected: all selected tests pass.

- [ ] **Step 10: Commit**

```bash
git add src/common/m_derived_types.fpp src/common/m_global_parameters_common.fpp \
        src/common/m_variables_conversion.fpp src/simulation/m_global_parameters.fpp \
        src/pre_process/m_global_parameters.fpp src/post_process/m_global_parameters.fpp \
        toolchain/mfc/params/definitions.py toolchain/mfc/params/descriptions.py
git commit -m "feat: add fluid_pp%k_therm input and heat_conduction gate"
```

---

### Task 2: Mixture temperature and `q_T_sf` population

Makes temperature available everywhere `q_prim` is, including ghost cells, whenever conduction is on. Still no physics change: nothing reads the field yet.

**Files:**
- Modify: `src/common/m_variables_conversion.fpp` (new function; `s_convert_conservative_to_primitive_variables` around line 662)
- Modify: `src/common/m_phase_change.fpp:278` (comment only)

**Interfaces:**
- Consumes: `heat_conduction`, `fluid_k_therm` from Task 1; the existing `cvs`, `isentrope_n`, `gammas`, `pi_infs` arrays.
- Produces:
  ```fortran
  function f_mixture_temperature(alpha_K, alpha_rho_K, pres, gamma_K, pi_inf_K) result(T)
      real(wp), dimension(num_fluids), intent(in) :: alpha_K, alpha_rho_K
      real(wp), intent(in)                        :: pres, gamma_K, pi_inf_K
      real(wp)                                    :: T
  ```
  and a `q_T_sf` that holds the mixture temperature over `idwbuff` whenever `heat_conduction` is true.

- [ ] **Step 1: Add the temperature function**

In `src/common/m_variables_conversion.fpp`, directly after `f_sg_thermal` (around line 1558), add:

```fortran
    !> Thermal-equilibrium mixture temperature for stiffened gas, from primitives.
    !! Algebraically identical to the conservative form in m_phase_change's s_infinite_pt_relaxation_k,
    !! T = (rho*e + p - sum(alpha_rho_i*qv_i)) / sum(alpha_rho_i*cv_i*n_i), because
    !! rho*e = gamma_mix*p + pi_inf_mix + sum(alpha_rho_i*qv_i) in MFC's stored variables.
    function f_mixture_temperature(alpha_K, alpha_rho_K, pres, gamma_K, pi_inf_K) result(T)

        $:GPU_ROUTINE(function_name='f_mixture_temperature', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_K, alpha_rho_K
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_K, alpha_rho_K
        #:endif
        real(wp), intent(in) :: pres, gamma_K, pi_inf_K
        real(wp)             :: T
        real(wp)             :: mCP  !< sum of alpha_rho_i*cp_i; cp_i = n_i*cv_i
        integer              :: i

        mCP = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            mCP = mCP + alpha_rho_K(i)*cvs(i)*isentrope_n(i)
        end do

        T = ((gamma_K + 1._wp)*pres + pi_inf_K)/max(mCP, sgm_eps)

    end function f_mixture_temperature
```

`alpha_K` is unused in the expression but is kept in the signature so the cylindrical and per-phase extensions in later work do not change callers. If the compiler warns about an unused dummy, drop `alpha_K` from both the signature and every call site rather than silencing the warning.

Add `f_mixture_temperature` to the module's `public` list at line 32, next to `s_phase_temperature`.

- [ ] **Step 2: Populate `q_T_sf` during conversion**

In `s_convert_conservative_to_primitive_variables`, the existing block reads:

```fortran
                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = pres

                    if (chemistry) then
                        q_T_sf%sf(j, k, l) = T
                    end if
```

Change it to:

```fortran
                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = pres

                    if (chemistry) then
                        q_T_sf%sf(j, k, l) = T
                    else if (heat_conduction) then
                        q_T_sf%sf(j, k, l) = f_mixture_temperature(alpha_K, alpha_rho_K, pres, gamma_K, pi_inf_K)
                    end if
```

`alpha_K`, `alpha_rho_K`, `gamma_K`, and `pi_inf_K` are already local to this loop iteration. Because this routine runs over `ibounds`, which the simulation calls with `idwbuff`, ghost-cell temperatures come out correct with no new halo exchange.

- [ ] **Step 3: Cross-reference the phase-change form**

In `src/common/m_phase_change.fpp`, change the comment above line 278 from `! common temperature` to:

```fortran
        ! common temperature; same closure as f_mixture_temperature in m_variables_conversion, written in conservative variables
```

- [ ] **Step 4: Leave the start-up gate alone**

`src/simulation/m_start_up.fpp:948` calls `s_compute_q_T_sf`, which lives in `m_chemistry.fpp:46` and divides by species partial densities. It is chemistry-only by construction and must not be called for conduction. No edit here: `s_convert_conservative_to_primitive_variables` runs over `idwbuff` before the first RHS evaluation, so Step 2 already fills `q_T_sf` everywhere conduction reads it. This step exists to record the decision, not to change code.

- [ ] **Step 5: Build**

```bash
./mfc.sh build -t simulation -j 8
```

Expected: clean build.

- [ ] **Step 6: Verify the temperature is right**

The Task 3 case does not exist yet, so verify against a single-fluid case whose exact temperature is known in closed form. Write `/tmp/check_T.py`:

```python
import numpy as np

# Single ideal-gas fluid, model_eqns = 2, num_fluids = 1 -> conservative record is
# [alpha_rho(1), mom_x, rho*E, alpha(1)], variable-major, C order.
NVAR, E_IDX = 4, 2
GAM, CV, NX = 1.4, 1.0, 100

q = np.fromfile("restart_data/lustre_0.dat", dtype=np.float64).reshape((NVAR, NX))
rho, mom, rhoE = q[0], q[1], q[E_IDX]
pres = (GAM - 1.0) * (rhoE - 0.5 * mom**2 / rho)          # pi_inf = 0
T_exact = pres / ((GAM - 1.0) * rho * CV)                 # f_sg_thermal, single fluid

T_num = np.fromfile("restart_data/lustre_T_0.dat", dtype=np.float64)  # written when T_wrt = T
print("max relative error:", np.abs(T_num - T_exact).max() / T_exact.max())
```

Run a 1D single-fluid case with `"fluid_pp(1)%cv": 1.0`, `"fluid_pp(1)%k_therm": 1.0e-3`, and post-process `T_wrt = T`, then run the script from that case's directory. Confirm the temperature output filename against the actual contents of `restart_data/` before trusting the comparison.

Expected: agreement to better than `1e-12` relative.

- [ ] **Step 7: Confirm no regression**

```bash
./mfc.sh test -j 8 --percent 20
```

Expected: all selected tests pass, phase-change tests included.

- [ ] **Step 8: Commit**

```bash
git add src/common/m_variables_conversion.fpp src/common/m_phase_change.fpp
git commit -m "feat: compute mixture temperature into q_T_sf when heat conduction is on"
```

---

### Task 3: Cartesian conduction flux, verified against the exact Laplacian

The core of the feature. The verification is a one-step consistency check, which isolates the new term completely: at $t=0$ with uniform pressure and zero velocity, every Euler flux vanishes, so the entire RHS is $\nabla\cdot(k\nabla T)$ and the first-step change in $\rho E$ is exactly $\Delta t\,k\nabla^2 T$ up to $O(\Delta t)$ time error and $O(\Delta x^2)$ space error.

**Files:**
- Create: `src/simulation/m_conduction.fpp`
- Create: `examples/1D_conduction_convergence/case.py`
- Create: `examples/1D_conduction_convergence/compare_analytic.py`
- Modify: `src/simulation/m_rhs.fpp:224` (allocation), `:719-724` (call site), `:727` (dispatch gate), `:1769`/`:1774`, `:1854`/`:1859`, `:1946`/`:1951` (RHS gates)
- Modify: `src/simulation/m_riemann_state.fpp:637` and its y/z counterparts (zeroing)

**Interfaces:**
- Consumes: `heat_conduction`, `fluid_k_therm` (Task 1); `q_T_sf` (Task 2).
- Produces:
  ```fortran
  subroutine s_compute_conduction_source_flux(idir, q_prim_qp, q_T_sf, flux_src_vf, irx, iry, irz)
      integer, intent(in)                                    :: idir
      type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_qp
      type(scalar_field), intent(in)                         :: q_T_sf
      type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
      type(int_bounds_info), intent(in)                      :: irx, iry, irz
  ```

- [ ] **Step 1: Write the failing verification case**

Create `examples/1D_conduction_convergence/case.py`. A sinusoidal temperature at uniform pressure on a periodic domain, with density set analytically so that $T(x) = T_0(1 + A\sin(2\pi x/L))$ exactly:

```python
#!/usr/bin/env python3
# 1D Fourier conduction verification: uniform pressure, zero velocity, periodic.
# Ideal gas (pi_inf = 0) with T = p / ((Gamma - 1) * rho * cv), so setting
#     rho(x) = RHO0 / (1 + A*sin(2*pi*x/L))
# gives exactly T(x) = T0 * (1 + A*sin(2*pi*x/L)) with T0 = p / ((Gamma-1)*RHO0*cv).
# At t = 0 velocity is zero and pressure is uniform, so every Euler flux vanishes and
#     d(rho*E)/dt = k * d2T/dx2 = -k * T0 * A * (2*pi/L)**2 * sin(2*pi*x/L)
# exactly. One time step therefore measures the conduction term in isolation.
import json
import os

NX = int(os.environ.get("NX", "100"))

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = 0.1
K_THERM = 1.0e-3
L = 1.0
DT = 1.0e-8

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": NX - 1,
            "n": 0,
            "p": 0,
            "dt": DT,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,  # periodic
            "bc_x%end": -1,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": P0,
            "patch_icpp(1)%alpha_rho(1)": f"{RHO0} / (1.0 + {AMP} * sin(2.0 * pi * x / {L}))",
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (GAM - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": CV,
            "fluid_pp(1)%k_therm": K_THERM,
            "parallel_io": "F",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
        }
    )
)
```

- [ ] **Step 2: Write the analytic comparison**

Create `examples/1D_conduction_convergence/compare_analytic.py`:

```python
#!/usr/bin/env python3
"""
Verification for examples/1D_conduction_convergence.

At t = 0 the state is quiescent (u = 0) at uniform pressure, so every Euler flux is
zero and the only RHS contribution is Fourier conduction. One time step therefore gives

    (rho*E(dt) - rho*E(0)) / dt  =  k * d2T/dx2 + O(dt) + O(dx^2)

with the exact right-hand side

    d2T/dx2 = -T0 * A * (2*pi/L)**2 * sin(2*pi*x/L).

Run at several NX and confirm the L2 error falls at second order.

    for nx in 50 100 200 400; do
        NX=$nx ./mfc.sh run examples/1D_conduction_convergence/case.py -n 1
        NX=$nx ./build/venv/bin/python3 examples/1D_conduction_convergence/compare_analytic.py
    done
"""

import os
import sys

import numpy as np

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = 0.1
K_THERM = 1.0e-3
L = 1.0
DT = 1.0e-8

NX = int(os.environ.get("NX", "100"))
HERE = os.path.dirname(os.path.abspath(__file__))
RESTART = os.path.join(HERE, "restart_data")

NVAR = 4  # alpha_rho(1), mom_x, E, alpha(1)
E_IDX = 2  # zero-based index of rho*E in the record


def read_step(step):
    path = os.path.join(RESTART, f"lustre_{step}.dat")
    if not os.path.exists(path):
        sys.exit(f"missing {path}; run the case first")
    raw = np.fromfile(path, dtype=np.float64)
    return raw.reshape((NVAR, NX))


def main():
    q0 = read_step(0)
    q1 = read_step(1)

    dx = L / NX
    x = (np.arange(NX) + 0.5) * dx

    t0 = P0 / ((GAM - 1.0) * RHO0 * CV)
    d2t_exact = -t0 * AMP * (2.0 * np.pi / L) ** 2 * np.sin(2.0 * np.pi * x / L)
    rhs_exact = K_THERM * d2t_exact

    rhs_num = (q1[E_IDX] - q0[E_IDX]) / DT

    err = np.sqrt(np.mean((rhs_num - rhs_exact) ** 2))
    scale = np.sqrt(np.mean(rhs_exact**2))
    print(f"NX={NX:5d}  L2 error={err:.6e}  relative={err / scale:.6e}")


if __name__ == "__main__":
    main()
```

The record layout matches `examples/2D_ibm_poiseuille_nn/compare_analytic.py`, which reads `restart_data/lustre_<step>.dat` as `np.fromfile(dtype=np.float64).reshape((NVAR, NY, NX))` — variable-major, C order, conservative variables. For `model_eqns = 2` with one fluid the record is `[alpha_rho(1), mom_x, rho*E, alpha(1)]`, so `NVAR = 4` and `E_IDX = 2`. Confirm `q0.shape == (4, NX)` on the first run; a wrong column would make this task's verification meaningless.

- [ ] **Step 3: Run it against the current build to see it fail**

```bash
NX=100 ./mfc.sh run examples/1D_conduction_convergence/case.py -n 1
NX=100 ./build/venv/bin/python3 examples/1D_conduction_convergence/compare_analytic.py
```

Expected: relative error of order 1 — the numerator is essentially zero because no conduction term exists yet, so the measured RHS is ~0 against a nonzero exact value. That is the red state.

- [ ] **Step 4: Write the conduction module**

Create `src/simulation/m_conduction.fpp`:

```fortran
!>
!! @file
!! @brief Contains module m_conduction

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Fourier heat conduction, div(k grad T), as a face-centered source flux on the energy equation.
!! The term has no cross-derivatives, so the direction-split face difference below is exact for it.
module m_conduction

    use m_global_parameters

    implicit none

    private; public :: s_compute_conduction_source_flux

    type(int_bounds_info) :: isc1, isc2, isc3
    $:GPU_DECLARE(create='[isc1, isc2, isc3]')
    integer, dimension(3) :: offsets_c
    $:GPU_DECLARE(create='[offsets_c]')

contains

    !> Accumulate -k*dT/dx_idir into the energy source flux at each idir-normal face.
    subroutine s_compute_conduction_source_flux(idir, q_prim_qp, q_T_sf, flux_src_vf, irx, iry, irz)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_qp
        type(scalar_field), intent(in)                         :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in)                      :: irx, iry, irz

        real(wp) :: k_face, dT_dxi, grid_spacing, alpha_face
        integer  :: x, y, z, i

        isc1 = irx; isc2 = iry; isc3 = irz
        offsets_c = 0
        offsets_c(idir) = 1

        $:GPU_UPDATE(device='[isc1, isc2, isc3, offsets_c]')

        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_face, dT_dxi, grid_spacing, alpha_face, i]')
        do z = isc3%beg, isc3%end
            do y = isc2%beg, isc2%end
                do x = isc1%beg, isc1%end
                    select case (idir)
                    case (1)
                        grid_spacing = x_cc(x + 1) - x_cc(x)
                    case (2)
                        grid_spacing = y_cc(y + 1) - y_cc(y)
                    case (3)
                        grid_spacing = z_cc(z + 1) - z_cc(z)
                    end select

                    ! Volume-fraction-weighted face conductivity. Raw cell-centered alphas over- and
                    ! undershoot near interfaces, so clamp the face average as the viscous path does.
                    k_face = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_face = 0.5_wp*(q_prim_qp(eqn_idx%adv%beg + i - 1)%sf(x, y, z) &
                                             & + q_prim_qp(eqn_idx%adv%beg + i - 1)%sf(x + offsets_c(1), &
                                                                                       & y + offsets_c(2), z + offsets_c(3)))
                        alpha_face = min(max(alpha_face, 0._wp), 1._wp)
                        k_face = k_face + alpha_face*fluid_k_therm(i)
                    end do

                    dT_dxi = (q_T_sf%sf(x + offsets_c(1), y + offsets_c(2), z + offsets_c(3)) &
                              & - q_T_sf%sf(x, y, z))/grid_spacing

                    flux_src_vf(eqn_idx%E)%sf(x, y, z) = flux_src_vf(eqn_idx%E)%sf(x, y, z) - k_face*dT_dxi
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_conduction_source_flux

end module m_conduction
```

The sign matches `m_chemistry.fpp:470`, which subtracts its energy diffusion flux (containing `lambda_Cell*dT_dxi`) from `flux_src_vf(eqn_idx%E)`.

- [ ] **Step 5: Allocate the energy source flux for conduction**

In `src/simulation/m_rhs.fpp`, change line 224 from:

```fortran
                    if (viscous .or. surface_tension) then
```

to:

```fortran
                    if (viscous .or. surface_tension .or. heat_conduction) then
```

Without this the energy slot of `flux_src_n` is never allocated in a conduction-only run.

- [ ] **Step 6: Zero it each Riemann sweep**

In `src/simulation/m_riemann_state.fpp`, `s_initialize_riemann_solver` has one `if (viscous .or. (surface_tension)) then` guard per `norm_dir` branch. Change every one of them to:

```fortran
            if (viscous .or. surface_tension .or. heat_conduction) then
```

Confirm with `grep -n "viscous .or. (surface_tension)" src/simulation/m_riemann_state.fpp` that all occurrences were changed; a missed branch leaves stale flux in that direction and shows up as a direction-dependent error.

- [ ] **Step 7: Widen the RHS gates so the flux is not discarded**

`flux_src(mom%beg:E)` reaches the RHS only through four guards, all currently written without `heat_conduction`. Change every one of them.

`src/simulation/m_rhs.fpp:727`:

```fortran
                    if (viscous .or. surface_tension .or. chem_params%diffusion .or. heat_conduction) then
```

Then, in `s_compute_additional_physics_rhs`, the three per-direction blocks at lines 1769 (x), 1854 (y), and 1946 (z). Each has an outer and an inner guard; widen both in all three:

```fortran
            if ((surface_tension .or. viscous) .or. chem_params%diffusion .or. heat_conduction) then
```

```fortran
                            if (surface_tension .or. viscous .or. heat_conduction) then
```

The inner loop runs over `i = eqn_idx%mom%beg, eqn_idx%E`. In a conduction-only run the momentum components of `flux_src` are allocated (Step 5) and zeroed (Step 6), so they contribute exactly zero and no extra branch is needed.

Verify with `grep -n "surface_tension .or. viscous\|viscous .or. surface_tension" src/simulation/m_rhs.fpp` that seven guards now mention `heat_conduction` — one dispatch plus two per direction. A missed direction shows up as a direction-dependent error in the Step 10 convergence run, which is exactly what that test is for.

- [ ] **Step 8: Call it from the RHS**

In `src/simulation/m_rhs.fpp`, after the chemistry diffusion block that ends at line 724, add:

```fortran
                    ! RHS for Fourier heat conduction
                    if (heat_conduction) then
                        call nvtxStartRange("RHS-CONDUCTION")
                        call s_compute_conduction_source_flux(id, q_prim_qp%vf, q_T_sf, flux_src_n(id)%vf, irx, iry, irz)
                        call nvtxEndRange
                    end if
```

Add `use m_conduction` to the module's `use` block, next to the other simulation-module imports.

- [ ] **Step 9: Build**

```bash
./mfc.sh build -t simulation -j 8
```

Expected: clean build. A new `.fpp` under `src/simulation/` is picked up by the CONFIGURE_DEPENDS glob in `cmake/Fypp.cmake`; if it is not, force a reconfigure with `./mfc.sh build --clean -t simulation -j 8`.

- [ ] **Step 10: Run the verification and confirm second order**

```bash
for nx in 50 100 200 400; do
    NX=$nx ./mfc.sh run examples/1D_conduction_convergence/case.py -n 1
    NX=$nx ./build/venv/bin/python3 examples/1D_conduction_convergence/compare_analytic.py
done
```

Expected: the relative error is small at NX=50 and falls by a factor near 4 for each doubling, i.e. an observed order near 2.0. Anything that plateaus means a coefficient or indexing error, not a discretization limit — do not proceed past this step until the order is right. A factor-of-2 offset in the magnitude at every resolution means a wrong constant, which the convergence rate alone will not reveal, so also check the relative error itself is at the percent level at NX=400, not merely converging.

- [ ] **Step 11: Confirm the term is inert when off**

```bash
./mfc.sh test -j 8 --percent 25
```

Expected: all selected tests pass. None of them set `k_therm`, so every one exercises the `heat_conduction == .false.` path and must be unchanged.

- [ ] **Step 12: Write the example README**

Create `examples/1D_conduction_convergence/README.md` stating what the case verifies, the exact solution used, the run commands from Step 10, and the observed orders recorded in Step 10.

- [ ] **Step 13: Commit**

```bash
git add src/simulation/m_conduction.fpp src/simulation/m_rhs.fpp src/simulation/m_riemann_state.fpp \
        examples/1D_conduction_convergence
git commit -m "feat: add Cartesian Fourier conduction flux with analytic verification"
```

---

### Task 4: Thermal time step constraint

Conduction imposes $\Delta t \le \mathrm{CFL}\,\Delta x^2\rho c_v/k$. Without it, `cfl_adap_dt` runs pick a step the conduction term is unstable at.

**Files:**
- Modify: `src/simulation/m_sim_helpers.fpp:21` (`dt_limiter_names`), `:108-180` (`s_compute_stability_from_dt`), `:183-253` (`s_compute_dt_from_cfl`)
- Modify: `src/simulation/m_time_steppers.fpp:660-745`
- Modify: `src/simulation/m_data_output.fpp:160-250`
- Modify: `src/common/m_mpi_common.fpp` (`s_mpi_reduce_stability_criteria_extrema`)

**Interfaces:**
- Consumes: `heat_conduction`, `fluid_k_therm` (Task 1).
- Produces: `max_dt` widened from `dimension(3)` to `dimension(4)` with `max_dt(4)` the thermal limit; `dt_limiter_names` gains `'TCFL'`; two signatures change:
  ```fortran
  subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, alpha, alpha_rho, j, k, l)
  subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, alpha, alpha_rho, j, k, l, icfl, vcfl, Rc, ccfl, tcfl)
  ```
  `alpha` and `alpha_rho` are `real(wp), dimension(num_fluids), intent(in)`; `tcfl` is `real(wp), intent(inout)`. Both are already computed by `s_compute_cell_state` in the caller's loop body (`m_time_steppers.fpp` around line 700 and `m_data_output.fpp` around line 210), so they are passed through, never recomputed.

- [ ] **Step 1: Widen the limiter name list**

In `src/simulation/m_sim_helpers.fpp:21`:

```fortran
    character(len=4), dimension(5), parameter :: dt_limiter_names = (/'ICFL', 'VCFL', 'CCFL', 'TCFL', 'COLL'/)
```

`'TCFL'` is inserted before `'COLL'` so it sits with the other CFL criteria; the candidate array in Step 3 must use the same order.

- [ ] **Step 2: Add the thermal `dt` candidate**

In `s_compute_dt_from_cfl`, change the declaration `real(wp), dimension(3), intent(out) :: max_dt` to `dimension(4)`, add `max_dt(4) = huge(1._wp)` beside the existing initializations, add a local `real(wp) :: tcfl_dt, k_mix, rho_cv`, and after the capillary block add:

```fortran
        ! Thermal diffusion CFL: dt <= cfl * dx^2 * rho * cv / k
        if (heat_conduction) then
            k_mix = 0._wp
            rho_cv = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                k_mix = k_mix + alpha(i)*fluid_k_therm(i)
                rho_cv = rho_cv + alpha_rho(i)*cvs(i)
            end do

            if (p > 0) then
                if (grid_geometry == 3) then
                    fltr_dtheta = f_compute_filtered_dtheta(k, l)
                    tcfl_dt = cfl_target*(min(dx(j), dy(k), fltr_dtheta)**2._wp)*rho_cv/max(k_mix, sgm_eps)
                else
                    tcfl_dt = cfl_target*(min(dx(j), dy(k), dz(l))**2._wp)*rho_cv/max(k_mix, sgm_eps)
                end if
            else if (n > 0) then
                tcfl_dt = cfl_target*(min(dx(j), dy(k))**2._wp)*rho_cv/max(k_mix, sgm_eps)
            else
                tcfl_dt = cfl_target*(dx(j)**2._wp)*rho_cv/max(k_mix, sgm_eps)
            end if
            max_dt(4) = tcfl_dt
        end if
```

`s_compute_dt_from_cfl` currently takes `(vel, c, max_dt, rho, Re_l, j, k, l)`. Add `alpha, alpha_rho` after `Re_l` as shown in the Interfaces block above, and pass the caller's existing locals at `m_time_steppers.fpp:720`:

```fortran
                    call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, alpha, alpha_rho, j, k, l)
```

Those two locals are filled by the `s_compute_cell_state` call a few lines earlier in the same loop body.

- [ ] **Step 3: Feed it into the step selection**

In `src/simulation/m_time_steppers.fpp`, declare `tcfl_dt_local` beside `ccfl_dt_local`, initialize it to `huge(1._wp)`, add it to the reduction clause of the enclosing `$:GPU_PARALLEL_LOOP`, and add after the `ccfl_dt_local` line:

```fortran
                    tcfl_dt_local = min(tcfl_dt_local, max_dt(4))
```

Widen `dt_candidates_loc`/`dt_candidates_glb` to `dimension(5)` and set, matching the `dt_limiter_names` order from Step 1:

```fortran
        dt_candidates_loc(1) = icfl_dt_local
        dt_candidates_loc(2) = vcfl_dt_local
        dt_candidates_loc(3) = ccfl_dt_local
        dt_candidates_loc(4) = tcfl_dt_local
        dt_candidates_loc(5) = coll_dt_local
```

Also widen the `max_dt` declaration at its declaration site in this file to `dimension(4)`.

- [ ] **Step 4: Report it in run-time info**

In `s_compute_stability_from_dt` (`src/simulation/m_sim_helpers.fpp:108`), add `alpha, alpha_rho` after `Re_l` and `tcfl` after `ccfl` per the Interfaces block, declare the locals `real(wp) :: k_mix, rho_cv` and `integer :: i`, and after the capillary block add:

```fortran
        ! Thermal diffusion CFL
        if (heat_conduction) then
            k_mix = 0._wp
            rho_cv = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                k_mix = k_mix + alpha(i)*fluid_k_therm(i)
                rho_cv = rho_cv + alpha_rho(i)*cvs(i)
            end do

            if (p > 0) then
                if (grid_geometry == 3) then
                    fltr_dtheta = f_compute_filtered_dtheta(k, l)
                    tcfl = dt*k_mix/(rho_cv*min(dx(j), dy(k), fltr_dtheta)**2._wp)
                else
                    tcfl = dt*k_mix/(rho_cv*min(dx(j), dy(k), dz(l))**2._wp)
                end if
            else if (n > 0) then
                tcfl = dt*k_mix/(rho_cv*min(dx(j), dy(k))**2._wp)
            else
                tcfl = dt*k_mix/(rho_cv*dx(j)**2._wp)
            end if
        end if
```

In `src/simulation/m_data_output.fpp`, mirror every `vcfl_max` construct for `tcfl_max`: the module variable at line 35, the `_loc`/`_glb` locals at line 165, the initialization at line 173, the `private` and `reduction` clauses at lines 179-180, the accumulation at line 216 (`tcfl_max_loc = max(tcfl_max_loc, merge(tcfl, 0.0_wp, heat_conduction))`), the serial fallback at line 230, and the header/value writes.

`s_mpi_reduce_stability_criteria_extrema` in `src/common/m_mpi_common.fpp` already takes ten positional arguments. Rather than adding an eleventh, convert its max-reduced scalars into one array argument and its min-reduced scalars into another, and update both call sites. This is the targeted cleanup the spec calls for in code this task has to touch anyway.

- [ ] **Step 5: Build**

```bash
./mfc.sh build -t simulation -j 8
```

Expected: clean build.

- [ ] **Step 6: Verify the limiter engages**

Add `"cfl_adap_dt": "T"`, `"cfl_target": 0.5`, and `"t_step_stop": 50` to a copy of `examples/1D_conduction_convergence/case.py` at `NX=400` with `k_therm` raised to `1.0`, so conduction dominates. Run with `run_time_info = T` and read the limiter column.

Expected: the reported limiter is `TCFL` for most steps, and the run is stable. With the thermal candidate removed (temporarily force `max_dt(4) = huge(1._wp)`), the same case should go unstable or pick a visibly larger step — confirm that contrast so the test is known to be load-bearing.

- [ ] **Step 7: Confirm no regression**

```bash
./mfc.sh test -j 8 --percent 25
```

Expected: all selected tests pass. Pay attention to the `cfl_adap_dt` cases near `toolchain/mfc/test/cases.py:3178`, which are viscous-CFL limited and must stay so.

- [ ] **Step 8: Commit**

```bash
git add src/simulation/m_sim_helpers.fpp src/simulation/m_time_steppers.fpp \
        src/simulation/m_data_output.fpp src/common/m_mpi_common.fpp
git commit -m "feat: add thermal diffusion CFL constraint and TCFL reporting"
```

---

### Task 5: Cylindrical axis source and golden tests

Away from the axis, cylindrical geometry needs no new physics: `m_rhs.fpp:1885-1897` and `:1913-1925` already apply `-0.5/y_cc(k)*(flux_src(k-1) + flux_src(k))` over `i = mom%beg, E` under `if (cyl_coord)` alone, so the $\tfrac{1}{r}k\,\partial_rT$ source appears as soon as the energy source flux carries conduction. The axis cell `k = 0`, which the generic loop skips, is the one place that needs code.

**Files:**
- Modify: `src/simulation/m_conduction.fpp` (axis routine)
- Modify: `src/simulation/m_rhs.fpp:350` (allocation), `:1815-1836` and `:1901-1911` (axis handling)
- Modify: `toolchain/mfc/test/cases.py`

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces:
  ```fortran
  subroutine s_compute_conduction_axis_source(q_prim_vf, q_T_sf, tau_Re_vf, ix, iy, iz)
      type(scalar_field), dimension(sys_size), intent(in)      :: q_prim_vf
      type(scalar_field), intent(in)                           :: q_T_sf
      type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
      type(int_bounds_info), intent(in)                        :: ix, iy, iz
  ```
  It accumulates, so it must run after `s_compute_viscous_stress_cylindrical_boundary`, which zeroes `tau_Re_vf(mom%beg:E)` at its top (`m_viscous.fpp:124-134`).

- [ ] **Step 1: Allocate `tau_Re_vf` for conduction-only runs**

`src/simulation/m_rhs.fpp:350` reads `if (viscous) then` and allocates `tau_Re_vf` plus the `dq_prim_d*_qp` gradient fields. Conduction needs only the energy component, so add a separate allocation rather than widening that guard and dragging the gradient fields in with it. Immediately after the `end if` that closes the `if (viscous)` block:

```fortran
            if (heat_conduction .and. .not. viscous) then
                @:ALLOCATE(tau_Re_vf(1:sys_size))
                @:ALLOCATE(tau_Re_vf(eqn_idx%E)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(tau_Re_vf(eqn_idx%E))
            end if
```

Mirror this in the finalizer beside the existing `@:DEALLOCATE(tau_Re_vf(eqn_idx%E)%sf)` at `m_rhs.fpp:2193`.

- [ ] **Step 2: Add the axis routine**

The generic geometric loop uses `flux_src(E)`, which holds `-k dT/dr` at faces. The axis substitution at `m_rhs.fpp:1909` uses `tau_Re_vf(E)` at cell centers in the same role, so the axis fill is the cell-centered value of the same quantity. Append to `src/simulation/m_conduction.fpp`:

```fortran
    !> Cell-centered -k*dT/dr for the axis cell of a cylindrical grid. The generic geometric source in
    !! m_rhs uses face values of flux_src(E) in this role; the axis cell has no face pair, so it reads
    !! tau_Re_vf(E) instead. Accumulates: call after s_compute_viscous_stress_cylindrical_boundary,
    !! which zeroes tau_Re_vf(mom%beg:E).
    subroutine s_compute_conduction_axis_source(q_prim_vf, q_T_sf, tau_Re_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in)      :: q_prim_vf
        type(scalar_field), intent(in)                           :: q_T_sf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
        type(int_bounds_info), intent(in)                        :: ix, iy, iz

        real(wp) :: k_cell, dT_dr, alpha_cell
        integer  :: j, k, l, i

        isc1 = ix; isc2 = iy; isc3 = iz
        $:GPU_UPDATE(device='[isc1, isc2, isc3]')

        if (.not. viscous) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = isc3%beg, isc3%end
                do k = isc2%beg, isc2%end
                    do j = isc1%beg, isc1%end
                        tau_Re_vf(eqn_idx%E)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_cell, dT_dr, alpha_cell, i]')
        do l = isc3%beg, isc3%end
            do k = isc2%beg + 1, isc2%end - 1
                do j = isc1%beg, isc1%end
                    k_cell = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_cell = min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), 0._wp), 1._wp)
                        k_cell = k_cell + alpha_cell*fluid_k_therm(i)
                    end do

                    dT_dr = (q_T_sf%sf(j, k + 1, l) - q_T_sf%sf(j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))

                    tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, l) - k_cell*dT_dr
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_conduction_axis_source
```

Add it to the module's `public` list. The `k` loop is inset by one because the central difference reads `k-1` and `k+1`; `m_rhs` only ever reads rows `-1`, `0`, and `1` of this array, all of which the inset range covers.

- [ ] **Step 3: Call it and widen the axis guards**

In `src/simulation/m_rhs.fpp`, inside `if (cyl_coord .and. ((bc_y%beg == -2) .or. (bc_y%beg == -14)))` starting at line 1815, the `if (viscous) then` block calls `s_compute_viscous_stress_cylindrical_boundary` and then applies the axis reflection. Restructure it so conduction participates:

```fortran
                if (viscous) then
                    call s_compute_viscous_stress_cylindrical_boundary(q_prim_vf, &
                        & dq_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), dq_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                        & dq_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                end if

                if (heat_conduction) then
                    call s_compute_conduction_axis_source(q_prim_vf, q_T_sf, tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                end if

                if (viscous .or. heat_conduction) then
                    $:GPU_PARALLEL_LOOP(private='[i, j, l]', collapse=2)
                    do l = 0, p
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%mom%beg, eqn_idx%E
                                rhs_vf(i)%sf(j, 0, l) = rhs_vf(i)%sf(j, 0, l) + 1._wp/(y_cc(1) - y_cc(-1))*(tau_Re_vf(i)%sf(j, &
                                       & -1, l) - tau_Re_vf(i)%sf(j, 1, l))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
```

The original code duplicates the `s_compute_viscous_stress_cylindrical_boundary` call in both arms of an `if (p > 0)` whose two branches are byte-identical (`m_rhs.fpp:1817-1824`); collapsing it to the single call above removes that duplication, which the DRY guideline in `AGENTS.md` asks for in code being edited. Diff the two branches before collapsing to confirm they really are identical.

Apply the same `viscous .or. heat_conduction` widening to the second axis block at `m_rhs.fpp:1901`:

```fortran
                    if (viscous .or. heat_conduction) then
                        $:GPU_PARALLEL_LOOP(private='[i, j, l]', collapse=2)
                        do l = 0, p
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, 0, l) = rhs_vf(i)%sf(j, 0, l) - 1._wp/y_cc(0)*tau_Re_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
```

Add `q_T_sf` to the argument list of the enclosing routine if it is not already in scope there; `m_rhs` holds it as a module-level field used by the chemistry diffusion call at line 722, so check before threading anything new.

- [ ] **Step 4: Build**

```bash
./mfc.sh build -t simulation -j 8
```

Expected: clean build.

- [ ] **Step 5: Verify the radial term against the exact Laplacian**

Create `examples/2D_axisym_conduction_convergence/` as an axisymmetric twin of the Task 3 case: `cyl_coord = T`, `bc_y%beg = -2` (axis), `n > 0`, uniform pressure, zero velocity, and a radial temperature profile whose cylindrical Laplacian is exact in closed form. Use

$$T(r) = T_0\left(1 + A\left(1 - \tfrac{r^2}{R^2}\right)\right), \qquad \nabla^2 T = \frac{1}{r}\frac{d}{dr}\left(r\frac{dT}{dr}\right) = -\frac{4AT_0}{R^2},$$

a constant, so the expected one-step change in $\rho E$ is the constant $-4kAT_0/R^2$ everywhere including the axis cell. Set density analytically as in Task 3 so that $T(r)$ comes out exactly, i.e. `patch_icpp(1)%alpha_rho(1)` as an expression in `y`.

Apply the same one-step measurement as `examples/1D_conduction_convergence/compare_analytic.py`, at several radial resolutions.

Expected: second-order convergence, and — because the exact answer is a constant — the axis cell's value must match the interior to the same order. An axis cell that is off by a factor or has the wrong sign stands out immediately against a constant field. That is why this profile was chosen over a sinusoid.

- [ ] **Step 6: Add the golden tests**

In `toolchain/mfc/test/cases.py`, following the `stack.push("Viscous", {...})` pattern at line 855, add a conduction group in the 1D Cartesian section:

```python
                stack.push("Conduction", {"fluid_pp(1)%k_therm": 1.0e-3, "fluid_pp(1)%cv": 1.0, "dt": 1e-11})
                cases.append(define_case_d(stack, "", {}))
                stack.pop()
```

and an equivalent one in the axisymmetric section, so both the Cartesian flux and the axis source are covered. Set `cv` explicitly: it defaults to zero, and a conducting fluid with `cv = 0` is rejected by the Task 6 checks.

- [ ] **Step 7: Generate the golden files**

```bash
./mfc.sh test -l | grep -i conduction
./mfc.sh test --generate -o Conduction -j 8
```

Expected: golden files appear under `tests/<hash>/` for the new traces only. Confirm with `git status` that no pre-existing golden file changed. If one did, the feature is not inert when off and Task 3 Step 11 missed it — stop and fix that before regenerating anything.

- [ ] **Step 8: Run the full suite**

```bash
./mfc.sh test -j 8
```

Expected: all tests pass, including the new conduction traces.

- [ ] **Step 9: Commit**

```bash
git add src/simulation/m_conduction.fpp src/simulation/m_rhs.fpp \
        toolchain/mfc/test/cases.py tests/ examples/2D_axisym_conduction_convergence
git commit -m "feat: cylindrical axis conduction source and conduction golden tests"
```

---

### Task 6: Input validation and documentation

**Files:**
- Modify: `src/simulation/m_checker.fpp`
- Modify: `toolchain/mfc/case_validator.py`
- Modify: `docs/documentation/case.md`, `docs/documentation/equations.md`

**Interfaces:**
- Consumes: everything above.
- Produces: `s_check_inputs_conduction`, called from `s_check_inputs`.

- [ ] **Step 1: Add the Fortran checks**

In `src/simulation/m_checker.fpp`, add:

```fortran
    !> Checks constraints on Fourier heat conduction inputs
    impure subroutine s_check_inputs_conduction

        integer :: i

        do i = 1, num_fluids
            @:PROHIBIT(fluid_pp(i)%k_therm < 0._wp, "fluid_pp(i)%k_therm must be non-negative")
            if (fluid_pp(i)%k_therm > 0._wp) then
                @:PROHIBIT(fluid_pp(i)%cv <= 0._wp, &
                           & "fluid_pp(i)%cv must be positive when fluid_pp(i)%k_therm is set: the mixture &
                           & temperature is undefined without it")
                @:PROHIBIT(fluid_pp(i)%eos /= eos_stiffened_gas .and. fluid_pp(i)%eos /= eos_ideal_gas, &
                           & "heat conduction supports only the stiffened-gas and ideal-gas equations of state")
            end if
        end do

        @:PROHIBIT(heat_conduction .and. igr, "heat conduction is not supported with igr")
        @:PROHIBIT(heat_conduction .and. chemistry, &
                   & "heat conduction is not supported with chemistry: the reacting path already carries &
                   & mixture-averaged conduction through chem_params%diffusion")

    end subroutine s_check_inputs_conduction
```

Call it unconditionally from `s_check_inputs` (`src/simulation/m_checker.fpp:23`), after `s_check_inputs_compilers`. The `k_therm < 0` check must run even when `heat_conduction` is false, so a negative value is caught rather than silently ignored.

The chemistry rule matters beyond tidiness: with both on, `s_convert_conservative_to_primitive_variables` takes the `chemistry` branch for `q_T_sf` while `m_rhs` runs both diffusion hooks, double-counting conduction. Prohibiting the combination is why Task 3 never had to guard the `flux_src(E)` allocation at `m_rhs.fpp:259` against a double allocate.

- [ ] **Step 2: Mirror the checks in the Python validator**

Add the same five rules to `toolchain/mfc/case_validator.py` so a bad case fails before a job is submitted rather than after it starts. Read `toolchain/mfc/test_case_validator.py` first: it is the unit-test file for that validator, and each new rule needs a case there asserting it rejects the bad input and accepts the good one. Run them with `./mfc.sh lint` or directly via the repo's pytest entry point.

- [ ] **Step 3: Verify each check fires**

For each of the five rules — negative `k_therm`, `k_therm` with `cv = 0`, `k_therm` on a Mie-Gruneisen fluid, `heat_conduction` with `igr`, `heat_conduction` with `chemistry` — copy `examples/1D_conduction_convergence/case.py`, break exactly that one rule, and run it.

Expected: the run aborts with the specific message above, not a generic failure and not a NaN partway through. Record all five messages. A rule that does not fire is worse than no rule, because it advertises a guarantee the code does not provide.

- [ ] **Step 4: Document the term**

In `docs/documentation/equations.md`, add the conduction term to the energy equation with the closure actually implemented: single thermal-equilibrium mixture temperature, $k=\sum_i\alpha_ik_i$, stiffened- and ideal-gas only.

In `docs/documentation/case.md`, document `fluid_pp(i)%k_therm` beside the `Re` entries: what it is, its units, that it requires `cv`, that it is independent of `viscous`, and that `igr` and non-stiffened equations of state are unsupported.

- [ ] **Step 5: Lint**

```bash
./mfc.sh lint
```

Expected: clean. The docs linters (`lint_docs.py`, `lint_param_docs.py`) check that every registered parameter is documented, so a missing `case.md` entry fails here.

- [ ] **Step 6: Full suite**

```bash
./mfc.sh test -j 8
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/simulation/m_checker.fpp toolchain/mfc/case_validator.py docs/documentation
git commit -m "feat: validate heat conduction inputs and document the term"
```

---

## Before opening the PR

- State in the PR description that it was written with Claude Code.
- Follow the PR template.
- This PR changes CFD results only for cases that set `k_therm`. Include the Task 3 Step 10 convergence table and the Task 5 Step 5 cylindrical table as the correctness verification the contributing guidelines require, and state that all pre-existing golden files are unchanged.
