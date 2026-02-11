@page contributing Contributing

# Contributing to MFC

We welcome contributions of all kinds: bug fixes, new features, documentation, tests, and issue triage.
This guide covers everything you need to get started and get your changes merged.

## Getting Set Up

1. **Fork and clone**
   ```bash
   git clone https://github.com/<your-user>/MFC.git
   cd MFC
   git remote add upstream https://github.com/MFlowCode/MFC.git
   ```
2. **Build MFC** (see @ref getting-started for full details):
   ```bash
   ./mfc.sh build -j $(nproc)
   ```
3. **Run the test suite** to verify your environment:
   ```bash
   ./mfc.sh test -j $(nproc)
   ```

## Architecture Overview

Understanding MFC's structure helps you know where to make changes and what they affect.

### Three-Phase Pipeline

MFC runs simulations in three phases, each a separate Fortran executable:

1. **pre_process** — Reads case parameters, generates initial conditions (patch geometries, flow states), writes binary grid and flow data to disk.
2. **simulation** — Reads the initial state, advances the solution in time via TVD Runge-Kutta integration, writes solution snapshots at specified intervals.
3. **post_process** — Reads simulation snapshots, computes derived quantities (vorticity, Schlieren, sound speed, etc.), writes Silo/HDF5 output for visualization.

All three share code in `src/common/`. Only `simulation` is GPU-accelerated.

### Directory Layout

| Directory | Contents |
|-----------|----------|
| `src/simulation/` | Time-stepping, RHS, Riemann solvers, WENO, physics models (GPU-accelerated) |
| `src/pre_process/` | Initial condition generation |
| `src/post_process/` | Derived variable computation and formatted output |
| `src/common/` | Derived types, global parameters, MPI, precision, I/O — shared by all three executables |
| `toolchain/` | Python CLI, parameter system, case validation, build orchestration |
| `tests/` | Golden files for 500+ regression tests |
| `examples/` | Sample case files |
| `docs/` | Doxygen documentation source |

### Simulation Data Flow

Each time step, `simulation` computes the right-hand side through this pipeline:

```
q_cons_vf (conservative variables: density, momentum, energy, volume fractions)
    → convert to primitive (density, velocity, pressure)
    → WENO reconstruction (left/right states at cell faces)
    → Riemann solve (numerical fluxes)
    → flux divergence + source terms (viscous, surface tension, body forces)
    → RHS assembly
    → Runge-Kutta update → q_cons_vf (next stage/step)
```

Key data structures (defined in `src/common/m_derived_types.fpp`):
- `scalar_field` — wraps a 3D `real(stp)` array (`%sf(i,j,k)`)
- `vector_field` — array of `scalar_field` (`%vf(1:sys_size)`)
- `q_cons_vf` / `q_prim_vf` — conservative and primitive state vectors

### Build Toolchain

`./mfc.sh` is a shell wrapper that invokes the Python toolchain (`toolchain/main.py`), which orchestrates:

1. **CMake** configures the build (compiler detection, dependencies, GPU backend)
2. **Fypp** preprocesses `.fpp` files into `.f90` (expands GPU macros, code generation)
3. **Fortran compiler** builds three executables from the generated `.f90` files

See @ref parameters for the full list of ~3,400 simulation parameters. See @ref case_constraints for feature compatibility and example configurations.

## Development Workflow

| Step | Command / Action |
|------|-----------------|
| Sync your fork | `git checkout master && git pull upstream master` |
| Create a branch on your fork | `git checkout -b feature/<short-name>` |
| Code, test, document | Follow the standards below |
| Run tests | `./mfc.sh test` |
| Commit | Clear, atomic commits (see below) |
| Push to your fork | `git push origin feature/<short-name>` |
| Open a PR | From your fork to `MFlowCode/MFC:master`. Every push triggers CI -- bundle changes to avoid flooding the queue |

### Commit Messages

- Start with a concise (50 chars or fewer) summary in imperative mood: `Fix out-of-bounds in EOS module`
- Add a blank line, then a detailed explanation if needed
- Reference related issues: `Fixes #123` or `Part of #456`

## Coding Standards

MFC is written in modern Fortran 2008+ with [Fypp](https://github.com/aradi/fypp) metaprogramming.
The standards below are split into **hard rules** (enforced in CI and review) and **soft guidelines** (goals for new code).

### Hard Rules

These are enforced. CI and reviewers will flag violations.

| Element | Rule |
|---------|------|
| Formatting | Enforced automatically by pre-commit hook (`./mfc.sh format` and `./mfc.sh lint`) |
| Indentation | 2 spaces; continuation lines align beneath `&` |
| Case | Lowercase keywords and intrinsics (`do`, `end subroutine`, ...) |
| Modules | `m_<feature>` (e.g. `m_riemann_solvers`) |
| Public subroutines | `s_<verb>_<noun>` (e.g. `s_compute_flux`) |
| Public functions | `f_<verb>_<noun>` (e.g. `f_create_bbox`) |
| Variables | Every argument has explicit `intent`; use `implicit none`, `dimension`/`allocatable`/`pointer` as appropriate |
| Forbidden | `goto`, `COMMON` blocks, global `save` variables |
| Error handling | Call `s_mpi_abort(<msg>)` -- never `stop` or `error stop` |
| GPU macros | Do not use raw OpenACC/OpenMP pragmas. Use the project's Fypp GPU macros (see below) |
| Compiler support | Code must compile with GNU gfortran, NVIDIA nvfortran, Cray ftn, and Intel ifx |

### Soft Guidelines

Aim for these in new and modified code. Existing code may not meet all of them.

| Element | Guideline |
|---------|-----------|
| Routine size | Prefer subroutine ≤ 500 lines, helper ≤ 150, function ≤ 100, file ≤ 1000 |
| Arguments | Prefer ≤ 6; consider a derived-type params struct for more |
| DRY | Avoid duplicating logic; factor shared code into helpers |

## Common Pitfalls

This section documents domain-specific issues that frequently appear in MFC contributions.
Both human reviewers and AI code reviewers reference this section.

### Array Bounds and Indexing

- MFC uses **non-unity lower bounds** (e.g., `idwbuff(1)%beg:idwbuff(1)%end` with negative ghost-cell indices). Always verify loop bounds match array declarations.
- **Riemann solver indexing:** Left states at `j`, right states at `j+1`. Off-by-one here corrupts fluxes.

### Precision and Type Safety

- **`stp` vs `wp` mixing:** In mixed-precision mode, `stp` (storage) may be half-precision while `wp` (working) is double. Conversions between them must be intentional, especially in MPI pack/unpack and RHS accumulation.
- **No double-precision intrinsics:** `dsqrt`, `dexp`, `dlog`, `dble`, `dabs`, `real(8)`, `real(4)` are forbidden. Use generic intrinsics with `wp` kind.
- **MPI type matching:** `mpi_p` must match `wp`; `mpi_io_p` must match `stp`. Mismatches corrupt communicated data.

### Memory and Allocation

- **ALLOCATE/DEALLOCATE pairing:** Every `@:ALLOCATE()` must have a matching `@:DEALLOCATE()`. Missing deallocations leak GPU memory.
- **@:ACC_SETUP_VFs / @:ACC_SETUP_SFs:** Vector/scalar fields must have GPU pointer setup before use in kernels.
- **Conditional allocation:** If an array is allocated inside an `if` block, its deallocation must follow the same condition.
- **Out-of-bounds access:** Fortran is permissive with assumed-shape arrays. Check that index arithmetic stays within declared bounds.

### MPI Correctness

- **Halo exchange:** Pack/unpack offset calculations (`pack_offset`, `unpack_offset`) must be correct for both interior and periodic boundaries. Off-by-one causes data corruption.
- **GPU data coherence:** Non-RDMA MPI requires `GPU_UPDATE(host=...)` before send and `GPU_UPDATE(device=...)` after receive. Missing these causes stale data.
- **Buffer sizing:** `halo_size` depends on dimensionality and QBMM state. `v_size` must account for extra bubble variables when QBMM is active.
- **Deadlocks:** Mismatched send/recv counts or tags across MPI ranks.

### Physics and Model Consistency

- **Pressure formula** must match `model_eqns` value. Model 2/3 (multi-fluid), model 4 (bubbles), MHD, and hypoelastic each use different EOS formulations. Wrong formula = wrong physics.
- **Conservative-primitive conversion:** Density recovery, kinetic energy, and pressure each have model-specific paths. Verify the correct branch is taken.
- **Volume fractions** must sum to 1. `alpha_rho_K` must be non-negative. Species mass fractions should be clipped to [0,1].
- **Boundary conditions:** Periodic BCs must match at both ends (`bc_x%beg` and `bc_x%end`). Cylindrical coordinates have special requirements (`bc_y%beg = -14` for axis in 3D).
- **Parameter constraints:** New parameters or physics features must be validated in `toolchain/mfc/case_validator.py`. New features should add corresponding validation.

### Python Toolchain

- New parameters in `toolchain/mfc/params/definitions.py` must have correct types, constraints, and tags.
- Validation in `case_validator.py` must cover new interdependencies.
- CLI schema in `toolchain/mfc/cli/commands.py` must match argument parsing.
- Check subprocess calls for shell injection risks and missing error handling.

### Compiler Portability

- Any compiler-specific code (`#ifdef __INTEL_COMPILER` etc.) must have fallbacks for all four supported compilers.
- Fypp macros must expand correctly for both GPU and CPU builds (macros are `#ifdef`'d out for non-GPU).
- No hardcoded GPU architectures without CMake detection.

### Architecture Notes

- **`src/common/` affects all three executables** (pre_process, simulation, post_process). Changes here have wide blast radius.
- No new global state; private helpers stay inside their defining module.
- Flag modifications to public subroutine signatures, parameter defaults, or output formats.
- Avoid unnecessary host/device transfers in hot loops, redundant allocations, and algorithmic inefficiency.

## Fypp and GPU

MFC uses [Fypp](https://github.com/aradi/fypp) macros (in `src/*/include/`) to generate accelerator-specific Fortran for OpenACC and OpenMP backends.
Only `simulation` (plus its `common` dependencies) is GPU-accelerated.

- **Raw OpenACC/OpenMP pragmas are not allowed.** Use the project's Fypp GPU macros instead.
- Add `collapse(n)` when safe, declare loop-local variables with `private(...)`.
- Avoid `stop`/`error stop` inside device code.
- Keep macros simple and readable.

See @ref gpuParallelization for the full GPU macro API reference, including all parameters, restrictions, examples, and debugging tools.

## How-To Guides

Step-by-step recipes for common development tasks.

### How to Add a New Simulation Parameter

Adding a parameter touches both the Python toolchain and Fortran source. Follow these steps in order. See @ref parameters for the full list of existing parameters and @ref case_constraints for feature compatibility.

**Step 1: Register in Python** (`toolchain/mfc/params/definitions.py`)

Add a call to `_r()` inside the `_load()` function:

```python
_r("my_param", REAL, {"my_feature_tag"})
```

The arguments are: name, type (`INT`, `REAL`, `LOG`, `STR`), and a set of feature tags. You can add an explicit description with `desc="..."`, otherwise one is auto-generated from `_SIMPLE_DESCS` or `_ATTR_DESCS`.

**Step 2: Add constraints** (same file, `CONSTRAINTS` dict)

If the parameter has valid ranges or choices:

```python
CONSTRAINTS = {
    # ...
    "my_param": {"min": 0, "max": 100},
    # or: "my_param": {"choices": [1, 2, 3]},
}
```

**Step 3: Add dependencies** (same file, `DEPENDENCIES` dict)

If enabling one parameter requires or recommends others:

```python
DEPENDENCIES = {
    # ...
    "my_param": {
        "when_true": {
            "requires": ["other_param"],
            "recommends": ["optional_param"],
        }
    },
}
```

Triggers include `when_true` (logical is `T`), `when_set` (parameter is not `None`), and `when_value` (parameter equals a specific value).

**Step 4: Add physics validation** (`toolchain/mfc/case_validator.py`)

If the parameter has cross-parameter constraints that go beyond simple min/max:

```python
def check_my_feature(self):
    if self.params["my_param"] > 0 and not self.params["other_param"]:
        self.errors.append("my_param requires other_param to be set")
```

**Step 5: Declare in Fortran** (`src/<target>/m_global_parameters.fpp`)

Add the variable declaration in the appropriate target's global parameters module. Choose the target(s) where the parameter is used:

- `src/pre_process/m_global_parameters.fpp`
- `src/simulation/m_global_parameters.fpp`
- `src/post_process/m_global_parameters.fpp`

```fortran
real(wp) :: my_param    !< Description of the parameter
```

If the parameter is used in GPU kernels, add a GPU declaration:

```fortran
$:GPU_DECLARE(create='[my_param]')
```

**Step 6: Add to Fortran namelist** (`src/<target>/m_start_up.fpp`)

Add the parameter name to the `namelist /user_inputs/` declaration:

```fortran
namelist /user_inputs/ ... , my_param, ...
```

The toolchain writes the parameter to the input file and Fortran reads it via this namelist. No other I/O code is needed.

**Step 7: Use in Fortran code**

Reference `my_param` anywhere in the target's modules. It is available as a global after the namelist is read at startup.

### How to Write a GPU Parallel Loop

All GPU loops use Fypp macros. See @ref gpuParallelization for the full API.

**Simple parallel loop** (3D with collapse):

```fortran
$:GPU_PARALLEL_LOOP(collapse=3)
do l = 0, p
    do k = 0, n
        do j = 0, m
            q_sf(j, k, l) = 0._wp
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

**With private variables** (temporaries local to each thread):

```fortran
$:GPU_PARALLEL_LOOP(collapse=3, private='[rho, pres, vel]')
do l = 0, p
    do k = 0, n
        do j = 0, m
            rho = q_prim_vf(1)%sf(j, k, l)
            pres = q_prim_vf(E_idx)%sf(j, k, l)
            ! ... use rho, pres as thread-local ...
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

**With reduction:**

```fortran
$:GPU_PARALLEL_LOOP(collapse=3, &
    & reduction='[[my_sum], [my_max]]', &
    & reductionOp='[+, MAX]', &
    & copy='[my_sum, my_max]')
do l = 0, p
    do k = 0, n
        do j = 0, m
            my_sum = my_sum + q_sf(j, k, l)
            my_max = max(my_max, q_sf(j, k, l))
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

**Sequential inner loop** within a parallel region:

```fortran
$:GPU_PARALLEL_LOOP(collapse=3)
do l = 0, p
    do k = 0, n
        do j = 0, m
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                alpha(i) = q_prim_vf(advxb + i - 1)%sf(j, k, l)
            end do
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

Key rules:
- Always pair `$:GPU_PARALLEL_LOOP(...)` with `$:END_GPU_PARALLEL_LOOP()`
- Use `collapse(n)` to fuse nested loops when the loop bounds are independent
- Declare all loop-local temporaries in `private='[...]'`
- Never use `stop` or `error stop` inside a GPU loop

### How to Allocate and Manage GPU Arrays

The full lifecycle of a GPU-resident array:

**Step 1: Declare** with GPU directive for module-level variables:

```fortran
real(wp), allocatable, dimension(:,:,:) :: my_array
$:GPU_DECLARE(create='[my_array]')
```

**Step 2: Allocate** in your initialization subroutine:

```fortran
@:ALLOCATE(my_array(0:m, 0:n, 0:p))
```

`@:ALLOCATE` handles both the Fortran `allocate` and the GPU `enter data create`.

**Step 3: Setup pointer fields** (only needed for derived types with pointer components like `scalar_field`):

```fortran
@:ALLOCATE(my_field%sf(0:m, 0:n, 0:p))
@:ACC_SETUP_SFs(my_field)
```

`@:ACC_SETUP_SFs` registers the pointer with the GPU runtime (required on Cray).

**Step 4: Deallocate** in your finalization subroutine, mirroring every allocation:

```fortran
@:DEALLOCATE(my_array)
```

If an array is allocated inside an `if` block, its deallocation must follow the same condition.

### How to Add a Test Case

**Step 1: Create a case file**

Test cases are Python scripts that print a JSON dict of parameters. See `examples/` for templates:

```python
#!/usr/bin/env python3
import json

print(json.dumps({
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "m": 49,
    "n": 0,
    "p": 0,
    "dt": 1e-6,
    "t_step_start": 0,
    "t_step_stop": 100,
    "t_step_save": 100,
    "num_patches": 1,
    "model_eqns": 2,
    "num_fluids": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "riemann_solver": 1,
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
}))
```

Keep grids small and runtimes short.

**Step 2: Register as a regression test** (`toolchain/mfc/test/cases.py`)

Add your case using the `Case` dataclass and the stack pattern for parameterized variations:

```python
stack.push("my_feature", {"my_param": value})
cases.append(define_case_d(stack, '', {}))
stack.pop()
```

**Step 3: Generate golden files**

```bash
./mfc.sh test --generate -o <test_id>
```

Golden files are stored as binary snapshots in `tests/<hash>/`.

**Step 4: Run**

```bash
./mfc.sh test -j $(nproc)
```

### How to Create a New Fortran Module

**Step 1: Create the file**

Name it `src/<target>/m_<feature>.fpp`. CMake auto-discovers `.fpp` files — no build system changes needed.

**Step 2: Use this boilerplate:**

```fortran
!> @file m_my_feature.fpp
!! @brief Description of the module

#:include 'case.fpp'
#:include 'macros.fpp'

module m_my_feature

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_my_feature, &
                       s_compute_my_feature, &
                       s_finalize_my_feature

    ! Module-level data
    real(wp), allocatable, dimension(:,:,:) :: work_array

contains

    !> Initialize module data
    impure subroutine s_initialize_my_feature()
        @:ALLOCATE(work_array(0:m, 0:n, 0:p))
    end subroutine s_initialize_my_feature

    !> Core computation
    subroutine s_compute_my_feature(q_prim_vf, rhs_vf)
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        ! ...
    end subroutine s_compute_my_feature

    !> Clean up module data
    impure subroutine s_finalize_my_feature()
        @:DEALLOCATE(work_array)
    end subroutine s_finalize_my_feature

end module m_my_feature
```

Key conventions:
- `private` by default, explicitly `public` for the module API
- Initialize/finalize subroutines for allocation lifecycle
- Every `@:ALLOCATE` has a matching `@:DEALLOCATE`
- Every argument has explicit `intent`

### Working with the Precision System

MFC supports double (default), single, and mixed precision. The types are defined in `src/common/m_precision_select.f90`:

| Type | Purpose | Example |
|------|---------|---------|
| `wp` | Working precision (computation) | `real(wp) :: velocity` |
| `stp` | Storage precision (I/O, field storage) | `real(stp), pointer :: sf(:,:,:)` |
| `mpi_p` | MPI type matching `wp` | `call MPI_BCAST(var, 1, mpi_p, ...)` |
| `mpi_io_p` | MPI type matching `stp` | Used in parallel I/O |

Rules:
- Use `real(wp)` for all computational variables
- Literal constants need the `_wp` suffix: `1.0_wp`, `3.14159_wp`, `1e-6_wp`
- Use **generic** intrinsics only: `sqrt`, `abs`, `sin`, `exp`, `log`, `max`, `min`
- **Forbidden** double-precision intrinsics: `dsqrt`, `dexp`, `dlog`, `dble`, `dabs`, `real(8)`, `real(4)`
- Conversions between `stp` and `wp` must be intentional, especially in MPI pack/unpack

### How to Extend MPI Halo Exchange

Halo exchange is in `src/simulation/m_mpi_proxy.fpp` (and `src/common/m_mpi_common.fpp` for buffer allocation).

To add new data to the halo exchange:

**Step 1: Update buffer sizing** (`src/common/m_mpi_common.fpp`)

`v_size` determines how many variables are packed per cell. If your new data adds fields per cell, increase `v_size`:

```fortran
v_size = sys_size + my_extra_fields
```

**Step 2: Add pack loop** (`src/simulation/m_mpi_proxy.fpp`)

Pack your data into the send buffer using a linear index:

```fortran
$:GPU_PARALLEL_LOOP(collapse=3, private='[j,k,l,r]')
do l = 0, p
    do k = 0, n
        do j = 0, buff_size - 1
            r = j + buff_size*(k + (n + 1)*l)
            buff_send(r) = my_data%sf(j + pack_offset, k, l)
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

**Step 3: GPU data coherence**

For non-RDMA MPI, add host/device transfers around the MPI call:

```fortran
$:GPU_UPDATE(host='[buff_send]')       ! GPU → CPU before send
call MPI_SENDRECV(buff_send, ..., buff_recv, ..., ierr)
$:GPU_UPDATE(device='[buff_recv]')     ! CPU → GPU after receive
```

**Step 4: Add unpack loop** mirroring the pack loop with `unpack_offset`.

### How to Add a Post-Processing Output Variable

Post-processing derived variables live in `src/post_process/m_derived_variables.fpp`.

**Step 1: Allocate storage** in `s_initialize_derived_variables_module`:

```fortran
if (my_var_wrt) then
    allocate(my_var_sf(-offset_x%beg:m + offset_x%end, &
                       -offset_y%beg:n + offset_y%end, &
                       -offset_z%beg:p + offset_z%end))
end if
```

**Step 2: Create derivation subroutine:**

```fortran
subroutine s_derive_my_variable(q_prim_vf, q_sf)
    type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
    real(wp), dimension(-offset_x%beg:m + offset_x%end, &
                        -offset_y%beg:n + offset_y%end, &
                        -offset_z%beg:p + offset_z%end), &
        intent(inout) :: q_sf
    integer :: i, j, k

    do k = -offset_z%beg, p + offset_z%end
        do j = -offset_y%beg, n + offset_y%end
            do i = -offset_x%beg, m + offset_x%end
                q_sf(i, j, k) = ! ... compute from q_prim_vf ...
            end do
        end do
    end do
end subroutine s_derive_my_variable
```

**Step 3: Call from output** in `m_data_output.fpp`:

```fortran
if (my_var_wrt) then
    call s_derive_my_variable(q_prim_vf, q_sf)
    call s_write_variable_to_formatted_database_file(q_sf, 'my_variable', dbfile, dbroot)
end if
```

### Modifying `src/common/`

Code in `src/common/` is compiled into all three executables (pre_process, simulation, post_process). Changes here have wide blast radius.

Checklist:
- Test all three targets: `./mfc.sh test` covers this
- If adding GPU code, remember that only `simulation` is GPU-accelerated. Guard GPU macros with `#:if MFC_SIMULATION`
- Check that new `use` statements don't create circular dependencies
- New modules need `implicit none` and explicit `intent` on all arguments

### Debugging

See @ref troubleshooting for debugging workflows, profiling tools, GPU diagnostic environment variables, common build/runtime errors, and fixes.

## Testing

MFC has 500+ regression tests. See @ref testing for the full guide.

- **Add tests** for any new feature or bug fix
- Use `./mfc.sh test --generate` to create golden files for new cases
- Keep tests fast: use small grids and short runtimes
- Test with `-a` to include post-processing validation

## CI Pipeline

Every push to a PR triggers CI. Understanding the pipeline helps you fix failures quickly.

### Lint Gate (runs first, blocks all other jobs)

All four checks must pass before any builds start:

1. **Formatting** — `./mfc.sh format` (auto-handled by pre-commit hook)
2. **Spelling** — `./mfc.sh spelling`
3. **Toolchain lint** — `./mfc.sh lint` (Python code quality)
4. **Source lint** — checks for:
   - Raw `!$acc` or `!$omp` directives (must use Fypp GPU macros)
   - Double-precision intrinsics (`dsqrt`, `dexp`, `dble`, etc.)

### Build and Test Matrix

After the lint gate passes:

- **Platforms:** Ubuntu and macOS
- **Compilers:** GNU (both), Intel OneAPI (Ubuntu only)
- **Modes:** debug + release, MPI + no-MPI, double + single precision
- **HPC runners:** Phoenix (NVIDIA/nvfortran), Frontier (AMD/Cray ftn) — both OpenACC and OpenMP backends
- **Retries:** Tests retry up to 3 times before failing
- **Cleanliness check:** Compiler warnings are tracked — your PR cannot increase the warning count

### Common CI Failures

| Failure | Fix |
|---------|-----|
| Formatting check | Pre-commit hook handles this; if you bypassed it, run `./mfc.sh format` |
| Raw pragma detected | Replace `!$acc`/`!$omp` with Fypp GPU macros (see @ref gpuParallelization) |
| Double-precision intrinsic | Use generic intrinsic with `wp` kind (e.g., `sqrt` not `dsqrt`) |
| Golden file mismatch | If intentional: `./mfc.sh test --generate --only <UUID>` |
| Warnings increased | Fix the new compiler warnings before merging |

See @ref troubleshooting for detailed debugging workflows.

## Documentation

- Add or update **Doxygen docstrings** in source files for new public routines
- Update **markdown docs** under `docs/` if user-facing behavior changes
- Provide a minimal **example case** in `examples/` for new features when practical

## Submitting a Pull Request

1. **PRs come from your fork.** Do not create branches on `MFlowCode/MFC` directly. Push to your fork and open a PR from there against `MFlowCode/MFC:master`.
2. **One PR = one logical change.** Split large changes into focused PRs.
3. **Fill out the PR template.** Remove checklist items that don't apply.
4. **Link issues** with `Fixes #<id>` or `Part of #<id>`.
5. **Ensure CI passes** before requesting review. Run `./mfc.sh test` locally first. Formatting and linting are handled automatically by the pre-commit hook.
6. **Describe your testing**: what you ran, which compilers/platforms you used.

If your change touches GPU code (`src/simulation/`), see the GPU checklist in the PR template.

## Code Review and Merge

- Respond to reviewer comments promptly
- Push focused updates; each push re-runs CI, so batch your fixes
- A maintainer will merge your PR once all reviews are approved and CI is green

If your PR is large or architectural, consider opening an issue first to discuss the approach.
