# m_eos Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the equation-of-state machinery out of `src/common/m_variables_conversion.fpp` into a new leaf module `src/common/m_eos.fpp`, with zero numerical change.

**Architecture:** One module absorbs the whole EOS chain so the deep call paths (`s_phase_pressure_on_isentrope -> s_rk4 -> s_ode_slope -> s_reference_curve`) stay internal to it. `m_eos` depends only on `m_derived_types`, `m_constants` and `m_global_parameters_common`, so it sits below the per-target `m_global_parameters` and is shared by all three executables. Device globals do not move.

**Tech Stack:** Fortran 2008 + Fypp macros, CMake/Ninja, `./mfc.sh` toolchain, golden-file regression tests.

**Spec:** `docs/superpowers/specs/2026-09-11-m-eos-extraction-design.md`

## Global Constraints

- **Bit-for-bit goldens.** Every task's test run must pass with goldens unchanged. Never run `./mfc.sh test --generate`. A golden diff means the task is wrong, not that the golden is stale.
- **No numerical change.** Moving code only. Do not reorder arithmetic, rename variables, change literals, or "tidy" expressions while moving.
- **Device globals stay in `m_global_parameters_common`.** Do not move `eoss`, `eos_coeffs`, `gammas`, `cvs`, `isentrope_n`, `isentrope_B`, `qvs`, `qvps`, `any_state_dependent_eos` or their `$:GPU_DECLARE` clauses.
- **Do not touch the `cray_inline` macro chain** in `src/common/include/parallel_macros.fpp` or its description in `docs/documentation/gpuParallelization.md`. Known issue, explicitly out of scope.
- **Formatting is enforced.** Every commit must pass `./mfc.sh format` and `./mfc.sh lint`.
- **No CMake edit is needed.** `cmake/Fypp.cmake:82` globs `src/common/*.fpp` with `CONFIGURE_DEPENDS`, and CMake's Fortran scanner orders modules automatically.
- Local verification is CPU-only. GPU backends are verified by CI: `.github/workflows/test.yml`, `frontier/`, `frontier_amd/`, `bench.yml`.

---

### Task 1: Capture the baseline

No code change. This establishes the safety net every later task is measured against.

**Files:**
- Create: `/tmp/mfc-eos-baseline.txt` (scratch, not committed)

- [ ] **Step 1: Confirm a clean tree at the branch point**

```bash
git status --porcelain    # expect empty
git log -1 --oneline      # expect the design-doc commit
```

- [ ] **Step 2: Build**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -20
```
Expected: builds `pre_process`, `simulation`, `post_process` with no errors.

- [ ] **Step 3: Run the full test suite and record the result**

```bash
./mfc.sh test -j $(nproc) 2>&1 | tee /tmp/mfc-eos-baseline.txt | tail -20
```
Expected: all tests pass. Record the pass count — later tasks must match it exactly.

- [ ] **Step 4: Record the file size being reduced**

```bash
wc -l src/common/m_variables_conversion.fpp    # expect 1997
```

No commit for this task.

---

### Task 2: Create `m_eos` and move the code, with temporary re-export

The move is mechanical. `m_variables_conversion` keeps re-exporting every moved name, so no consumer changes and the build stays green. Consumers migrate in Task 3.

**Files:**
- Create: `src/common/m_eos.fpp`
- Modify: `src/common/m_variables_conversion.fpp` (delete lines 1370-1756 and 1792-1829; add `use m_eos`; adjust `public ::` list)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: module `m_eos`, public names — `f_pressure(e_int, gamma, pi_inf, qv) result(pres)`, `f_bulk_modulus(pres, gamma, pi_inf) result(blkmod)`, `f_relativistic_enthalpy(pres, rho, gamma) result(H)`, `f_isentrope_exponent(gamma) result(n)`, `f_isentrope_pressure(pi_inf, gamma) result(B)`, `f_sg_thermal(pres, rho_or_T, n, B, cv) result(T_or_rho)`, `f_is_state_dependent(i) result(yes)`, `s_phase_coefficients(alpha_rho, alpha, i, rho, gamma, pi_inf, dpi, dgamma)`, `s_phase_pressure_on_isentrope(pres, rho, xi, i, p_isen)`, `s_phase_temperature(rho, pres, i, T)`, `s_phase_density_on_isentrope(i, rho_from, p_from, p_to, rho_to, c2_to)`, `s_phase_internal_energy(pres, alpha, alpha_rho, i, e_phase)`, `s_phase_bulk_modulus(pres, alpha, alpha_rho, i, blkmod)`. All signatures are unchanged from their current definitions.

- [ ] **Step 1: Extract the two contiguous blocks to a scratch file**

Block A is `s_reference_curve` through `end subroutine s_phase_bulk_modulus`. Block B is `f_pressure` through `end function f_relativistic_enthalpy`. The hypoelastic pair between them (1758-1790) stays behind.

```bash
sed -n '1370,1756p' src/common/m_variables_conversion.fpp  > /tmp/eos_block_a.f
sed -n '1792,1829p' src/common/m_variables_conversion.fpp  > /tmp/eos_block_b.f
head -3 /tmp/eos_block_a.f   # expect the s_reference_curve doc comment
tail -3 /tmp/eos_block_b.f   # expect "end function f_relativistic_enthalpy"
```

- [ ] **Step 2: Write the new module header**

Create `src/common/m_eos.fpp` with exactly this header, then append Block A followed by Block B, then the closing `end module`:

```fortran
!>
!! @file
!! @brief Contains module m_eos

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Equations of state in Gamma/Pi form, rho e = Gamma(rho) p + Pi(rho).
!!
!! Stiffened and ideal gas keep constant coefficients, resolved once at start-up. The
!! state-dependent families (Mie-Gruneisen, JWL, Vinet) evaluate theirs per cell from a
!! reference curve. Mixture closure rules -- Wood's law, the six-equation mean, the bubbly
!! branch -- are not equations of state and live in m_variables_conversion.
!!
!! This module is a leaf: it reads the material arrays in m_global_parameters_common and
!! depends on nothing else in MFC. Adding an EOS family means one case in s_reference_curve.
module m_eos

    use m_derived_types
    use m_global_parameters_common
    use m_constants, only: eos_stiffened_gas, eos_ideal_gas, eos_mie_gruneisen, eos_jwl, eos_vinet, eos_rk4_steps, &
        & ode_isentrope, ode_reference_temperature, sgm_eps, verysmall, dflt_real

    implicit none

    private

    public :: f_pressure, f_bulk_modulus, f_relativistic_enthalpy, f_isentrope_exponent, f_isentrope_pressure, f_sg_thermal, &
        & f_hugoniot_compression_limit, f_is_state_dependent, s_phase_coefficients, s_phase_pressure_on_isentrope, &
        & s_phase_temperature, s_phase_density_on_isentrope, s_phase_internal_energy, s_phase_bulk_modulus

contains

    ! <<< Block A here, then Block B >>>

end module m_eos
```

`f_hugoniot_compression_limit` is public in this task **only** because `s_initialize_variables_conversion_module` still calls it (it is not moved until Task 4, which makes it private again). Everything else in the family layer is private from the start — `s_eos_coefficients`, `s_reference_curve`, `s_ode_slope`, `s_rk4`, `s_phase_c2`, `f_has_isentropic_reference` and `f_c2_from_coefficients` have no caller outside this module.

- [ ] **Step 3: Delete the moved blocks from `m_variables_conversion`**

Delete the later range first so the earlier line numbers stay valid:

```bash
sed -i '1792,1829d' src/common/m_variables_conversion.fpp
sed -i '1370,1756d' src/common/m_variables_conversion.fpp
wc -l src/common/m_variables_conversion.fpp    # expect 1572
```

- [ ] **Step 4: Add the `use` and keep the re-export**

In `src/common/m_variables_conversion.fpp`, add after line 11 (`use m_derived_types`):

```fortran
    use m_eos
```

Leave the existing `public ::` list untouched. Because `m_variables_conversion` declares `private` and then names the moved routines in its `public` list, those names re-export from `m_eos` through it and every existing consumer keeps compiling unchanged. Remove only `s_eos_coefficients` from that list — it is private to `m_eos` now and naming it would not compile.

- [ ] **Step 5: Build**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -20
```
Expected: clean build. A "symbol not found" error here names a routine that was moved but still referenced — add it to the `m_eos` public list rather than moving it back.

- [ ] **Step 6: Run the full test suite**

```bash
./mfc.sh test -j $(nproc) 2>&1 | tail -20
```
Expected: identical pass count to Task 1. Any golden diff means arithmetic changed during the move — revert and redo the extraction with `sed`, not by hand.

- [ ] **Step 7: Format, lint, commit**

```bash
./mfc.sh format && ./mfc.sh lint
git add src/common/m_eos.fpp src/common/m_variables_conversion.fpp
git commit -m "refactor: move the equation-of-state machinery into m_eos"
```

---

### Task 3: Migrate consumers to `use m_eos` and drop the re-export

**Files:**
- Modify: `src/common/m_phase_change.fpp` (`f_sg_thermal`)
- Modify: `src/common/m_variables_conversion.fpp` (drop moved names from `public ::`)
- Modify: `src/post_process/m_data_output.fpp` (`s_phase_internal_energy`)
- Modify: `src/post_process/m_derived_variables.fpp` (`f_isentrope_exponent`, `f_isentrope_pressure`)
- Modify: `src/post_process/m_start_up.fpp` (`s_phase_temperature`)
- Modify: `src/simulation/m_acoustic_src.fpp` (`s_phase_bulk_modulus`, `f_bulk_modulus`)
- Modify: `src/simulation/m_bubbles_EE.fpp` (`f_isentrope_exponent`, `f_isentrope_pressure`)
- Modify: `src/simulation/m_bubbles_EL.fpp` (`f_pressure`, `f_bulk_modulus`)
- Modify: `src/simulation/m_hypoelastic.fpp` (`s_phase_bulk_modulus`, `f_bulk_modulus`)
- Modify: `src/simulation/m_ibm.fpp` (`s_phase_internal_energy`)
- Modify: `src/simulation/m_igr.fpp` (`f_pressure`, `f_bulk_modulus`)
- Modify: `src/simulation/m_pressure_relaxation.fpp` (`f_is_state_dependent`, `s_phase_coefficients`, `s_phase_density_on_isentrope`, `s_phase_internal_energy`, `f_pressure`)
- Modify: `src/simulation/m_qbmm.fpp` (`f_bulk_modulus`)
- Modify: `src/simulation/m_reactive_burn.fpp` (`s_phase_temperature`, `f_pressure`)
- Modify: `src/simulation/m_rhs.fpp` (`s_phase_bulk_modulus`)
- Modify: `src/simulation/m_riemann_solver_hll.fpp` (`f_isentrope_exponent`, `f_isentrope_pressure`, `f_relativistic_enthalpy`)
- Modify: `src/simulation/m_riemann_solver_hllc.fpp` (`f_isentrope_exponent`, `f_isentrope_pressure`, `s_phase_pressure_on_isentrope`, `s_phase_internal_energy`)
- Modify: `src/simulation/m_riemann_solver_hypo_hlld.fpp` (`s_phase_bulk_modulus`)
- Modify: `src/simulation/m_start_up.fpp` (`s_phase_internal_energy`)

**Interfaces:**
- Consumes: the `m_eos` public names produced by Task 2.
- Produces: nothing new. After this task `m_variables_conversion` no longer re-exports any EOS name.

- [ ] **Step 1: Add the import to each consumer**

Two shapes exist in the tree. For a module with a bare `use m_variables_conversion` (e.g. `m_riemann_solver_hllc.fpp:13`, `m_phase_change.fpp:14`, `m_igr.fpp:13`), add a line immediately after it:

```fortran
    use m_variables_conversion
    use m_eos
```

For a module with an `only:` list (e.g. `m_pressure_relaxation.fpp:15-16`), split the list so EOS names come from `m_eos`:

```fortran
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel
    use m_eos, only: f_pressure, s_phase_internal_energy, s_phase_coefficients, s_phase_density_on_isentrope, &
        & f_is_state_dependent
```

- [ ] **Step 2: Drop the moved names from the `m_variables_conversion` public list**

Remove exactly these from the `public ::` statement, leaving the rest in place: `f_bulk_modulus`, `f_pressure`, `s_phase_internal_energy`, `f_isentrope_exponent`, `f_isentrope_pressure`, `f_sg_thermal`, `f_relativistic_enthalpy`, `s_phase_coefficients`, `s_phase_pressure_on_isentrope`, `s_phase_temperature`, `f_is_state_dependent`, `s_phase_bulk_modulus`, `s_phase_density_on_isentrope`.

Keep `gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps` in the list — those are `m_global_parameters_common` arrays re-exported for convenience, not EOS routines, and moving them is out of scope.

- [ ] **Step 3: Build**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -30
```
Expected: clean. Every error here is a consumer that still needs `use m_eos`; the message names the file and the symbol.

- [ ] **Step 4: Verify no consumer still relies on the re-export**

```bash
grep -rn "use m_eos" src/ --include=*.fpp | wc -l    # expect 19 (18 consumers + m_variables_conversion)
```

- [ ] **Step 5: Run the full test suite**

```bash
./mfc.sh test -j $(nproc) 2>&1 | tail -20
```
Expected: identical pass count to Task 1.

- [ ] **Step 6: Format, lint, commit**

```bash
./mfc.sh format && ./mfc.sh lint
git add src/
git commit -m "refactor: import the EOS operators from m_eos directly"
```

---

### Task 4: Split the initialization

This is the one part of the move that is not mechanical: `Gs_vc` is assigned in the middle of the EOS coefficient loop and shares the closing `$:GPU_UPDATE`.

**Files:**
- Modify: `src/common/m_eos.fpp` (add `s_initialize_eos_module`)
- Modify: `src/common/m_variables_conversion.fpp:271-350` (remove the EOS portion, keep `Gs_vc`)
- Modify: `src/pre_process/m_start_up.fpp`, `src/simulation/m_start_up.fpp`, `src/post_process/m_start_up.fpp` (call the new routine)

**Interfaces:**
- Consumes: the `m_eos` routines `f_isentrope_exponent`, `f_isentrope_pressure`, `f_hugoniot_compression_limit` and `f_is_state_dependent`, all of which become module-internal calls once the init routine lands here.
- Produces: `impure subroutine s_initialize_eos_module()` — no arguments. Must be called before `s_initialize_variables_conversion_module`.

- [ ] **Step 1: Add `s_initialize_eos_module` to `m_eos`**

Move lines 271-278 (the EOS `@:ALLOCATE` calls), the body of the `do i = 1, num_fluids` loop at 281-338 **excluding line 294** (`Gs_vc(i) = fluid_pp(i)%G`), and lines 339-350 (the case-optimization check and both `$:GPU_UPDATE` calls) into a new routine. Drop `Gs_vc` from the device update:

```fortran
    !> Resolve every fluid's EOS coefficients once, before any conversion runs.
    impure subroutine s_initialize_eos_module()

        integer :: i
        logical :: state_dependent  !< Whether this case's fluids need a density-dependent EOS

        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(eoss (1:num_fluids))
        @:ALLOCATE(isentrope_n (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(isentrope_B(1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps    (1:num_fluids))

        ! <<< the loop body from lines 281-338, minus line 294 >>>
        ! <<< the case-optimization block from lines 339-346 >>>

        $:GPU_UPDATE(device='[gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps, eoss, eos_coeffs]')
        #:if not MFC_CASE_OPTIMIZATION
            $:GPU_UPDATE(device='[any_state_dependent_eos]')
        #:endif

    end subroutine s_initialize_eos_module
```

Add `s_initialize_eos_module` to the `m_eos` public list, and remove `f_hugoniot_compression_limit` from it — its only caller is now inside this module.

- [ ] **Step 2: Leave `Gs_vc` behind in `m_variables_conversion`**

`@:ALLOCATE(Gs_vc (1:num_fluids))` is at line 279, outside the ranges being moved — leave it exactly where it is. What needs a new home is the single assignment at line 294, which loses its enclosing loop. Replace the deleted block with:

```fortran
        do i = 1, num_fluids
            Gs_vc(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs_vc]')
```

Delete the now-unused `state_dependent` local from `s_initialize_variables_conversion_module`'s declarations. Keep `i` — the new loop and the `Re_idx` loop below both use it.

- [ ] **Step 3: Call the new routine from all three start-ups**

In each of `src/pre_process/m_start_up.fpp`, `src/simulation/m_start_up.fpp` and `src/post_process/m_start_up.fpp`, add `use m_eos, only: s_initialize_eos_module` and insert the call on the line immediately before the existing `s_initialize_variables_conversion_module(...)` call:

```fortran
    call s_initialize_eos_module()
```

Find the existing call sites with:

```bash
grep -rn "s_initialize_variables_conversion_module" src/*/m_start_up.fpp
```

- [ ] **Step 4: Build**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -20
```
Expected: clean build.

- [ ] **Step 5: Run the full test suite**

```bash
./mfc.sh test -j $(nproc) 2>&1 | tail -20
```
Expected: identical pass count to Task 1. A failure that appears only in some cases means the new call is ordered after a consumer of `gammas`/`eos_coeffs` — move it earlier.

- [ ] **Step 6: Format, lint, commit**

```bash
./mfc.sh format && ./mfc.sh lint
git add src/
git commit -m "refactor: give m_eos its own initialization routine"
```

---

### Task 5: Mark the two cross-module device helpers

`s_phase_pressure_on_isentrope` and `s_phase_temperature` are called cross-module and carry a bare `parallelism='[seq]'`. Bring them in line with every other cross-module device helper.

**Files:**
- Modify: `src/common/m_eos.fpp` (two `$:GPU_ROUTINE` directives)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing. Signatures unchanged.

- [ ] **Step 1: Add the attributes**

Change the directive inside `s_phase_pressure_on_isentrope` from

```fortran
        $:GPU_ROUTINE(parallelism='[seq]')
```

to

```fortran
        $:GPU_ROUTINE(function_name='s_phase_pressure_on_isentrope', parallelism='[seq]', cray_inline=True)
```

and the one inside `s_phase_temperature` to

```fortran
        $:GPU_ROUTINE(function_name='s_phase_temperature', parallelism='[seq]', cray_inline=True)
```

Change nothing else. Do not touch `s_reference_curve`, `s_eos_coefficients`, `s_ode_slope`, `s_rk4` or `s_phase_c2` — they are private to `m_eos` and correctly carry the bare directive.

- [ ] **Step 2: Build and test**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -20
./mfc.sh test -j $(nproc) 2>&1 | tail -20
```
Expected: clean build, identical pass count. On a CPU build this directive expands to nothing, so goldens cannot move.

- [ ] **Step 3: Format, lint, commit**

```bash
./mfc.sh format && ./mfc.sh lint
git add src/common/m_eos.fpp
git commit -m "refactor: mark the cross-module EOS device helpers for inlining"
```

---

### Task 6: Update the documentation

**Files:**
- Modify: `docs/documentation/contributing.md:467-490` (the "How to Add an Equation of State" section)
- Modify: `docs/documentation/equations.md:28` (key source files)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

- [ ] **Step 1: Correct the location claim in `contributing.md`**

The section opens with "Every stiffened-gas expression lives in `src/common/m_variables_conversion.fpp`." Replace that sentence with:

```markdown
The equation-of-state operators live in `src/common/m_eos.fpp`; the mixture closure rules that
combine them (`s_compute_mixture_coefficients`, `s_compute_speed_of_sound` and their variants)
stay in `src/common/m_variables_conversion.fpp`. Adding a second EOS means supplying these, not
grepping for `gammas`:
```

Then add a row to the operator table, after the `f_sg_thermal` row:

```markdown
| `s_reference_curve` | the reference curve \f$p_{ref}, e_{ref}\f$ and \f$\Gamma_G\f$ of a state-dependent family - one `case` per family, and nothing else |
```

Leave the rest of the section, including the mechanical/caloric split and the stored-forms paragraph, unchanged.

- [ ] **Step 2: Correct the source-file pointer in `equations.md`**

Change line 28 from `src/common/m_variables_conversion.fpp` (EOS and variable conversion) to:

```markdown
**Key source files:** `src/simulation/m_rhs.fpp` (RHS evaluation), `src/common/m_eos.fpp` (equations of state), `src/common/m_variables_conversion.fpp` (variable conversion and mixture rules).
```

- [ ] **Step 3: Verify no other doc points at the old location**

```bash
grep -rn "m_variables_conversion" docs/documentation/
```
Expected: only references that genuinely concern conversion or mixture rules remain.

- [ ] **Step 4: Commit**

```bash
git add docs/documentation/
git commit -m "docs: point the equation-of-state guide at m_eos"
```

---

### Task 7: Final verification

**Files:** none modified.

- [ ] **Step 1: Confirm the size reduction**

```bash
wc -l src/common/m_variables_conversion.fpp src/common/m_eos.fpp
```
Expected: `m_variables_conversion` near 1500 lines (down from 1997), `m_eos` near 530. Task 2 removes 425 lines and Task 4 a further ~77.

- [ ] **Step 2: Confirm `m_eos` is a leaf**

```bash
grep -n "^ *use " src/common/m_eos.fpp
```
Expected: exactly `m_derived_types`, `m_global_parameters_common`, `m_constants`. Anything else means a dependency crept in and the module is no longer shareable across all three targets.

- [ ] **Step 3: Confirm the goldens never moved**

```bash
git diff --stat upstream/master -- tests/
```
Expected: empty. Any output here is a failed refactor.

- [ ] **Step 4: Confirm the whole diff is a move**

```bash
git diff upstream/master --stat -- src/
```
Expected: `m_eos.fpp` created, `m_variables_conversion.fpp` shrunk by a comparable amount, and small `use`-line edits in the 18 consumers plus three start-ups.

- [ ] **Step 5: Full clean build and test**

```bash
./mfc.sh build -j $(nproc) 2>&1 | tail -20
./mfc.sh test -j $(nproc) 2>&1 | tail -20
```
Expected: identical pass count to Task 1.

- [ ] **Step 6: Push and open the PR**

GPU backends are not verifiable locally. Push the branch and let CI cover them:
`test.yml` (NVHPC/CPU lanes), `frontier/` (Cray), `frontier_amd/` (amdflang), `bench.yml` (performance).

Because the kernel count is unchanged, the amdflang whole-image codegen effect does not confound the benchmark comparison — a wall-time regression in `bench.yml` is a real signal, not noise, and should be investigated before merge.

```bash
git push -u origin module-eos
```
