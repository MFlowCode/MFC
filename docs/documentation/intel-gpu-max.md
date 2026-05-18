# Building MFC for Intel Data Center GPU Max (Ponte Vecchio)

This documents how to build and run MFC with Intel GPU Max (Xe-HPC / Ponte Vecchio)
using ifx OpenMP target offload to SPIR64, as tested on GT CRNCH RoboGator (`dash3`).

## System configuration

| Component | Version / Path |
|---|---|
| Hardware | Intel Data Center GPU Max 1100 (Ponte Vecchio, PCI 8086:0bda) |
| OS | Linux (RHEL 8 compatible, kernel 5.15) |
| Fortran compiler | ifx 2025.3.3 (`/opt/intel/oneapi/compiler/2025.3/`) |
| MKL | oneMKL 2026.0 (`/opt/intel/oneapi/mkl/2026.0/`) |
| SYCL runtime | `libsycl.so` in `/opt/intel/oneapi/compiler/2026.0/lib/` |
| GPU device | `/dev/dri/renderD128` (requires `render` group membership) |

## Environment setup

Load the required oneAPI environment before building or running:

```bash
export PATH=/opt/intel/oneapi/compiler/2025.3/bin:$PATH
export MKLROOT=/opt/intel/oneapi/mkl/2026.0
export LIBRARY_PATH=/opt/intel/oneapi/compiler/2026.0/lib:\
/opt/intel/oneapi/compiler/2025.3/lib:\
${MKLROOT}/lib:\
/opt/intel/oneapi/tbb/2022.1/lib:\
$LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/umf/1.1/lib:\
/opt/intel/oneapi/compiler/2026.0/lib:\
/opt/intel/oneapi/compiler/2025.3/lib:\
${MKLROOT}/lib:\
/opt/intel/oneapi/tbb/2022.1/lib:\
$LD_LIBRARY_PATH
export FC=/opt/intel/oneapi/compiler/2025.3/bin/ifx
```

> **Important**: `FC` must be set explicitly. Without it, CMake may cache an older
> ifx (2025.0) from a system module, which does not support `need_device_addr` in
> the MKL 2026.0 OpenMP offload Fortran module.

> **Important**: `LIBRARY_PATH` (not just `LD_LIBRARY_PATH`) must include the
> compiler 2026.0 lib directory so the linker finds `libsycl.so` at build time.

> **Important**: `/opt/intel/oneapi/umf/1.1/lib` must be in `LD_LIBRARY_PATH` at
> runtime. The Level Zero and OpenCL UR adapters in the 2026.0 compiler depend on
> `libumf.so.1`, which lives in the separate `umf/1.1` package, not in the compiler
> lib directories themselves.

## Building

```bash
./mfc.sh build -t simulation --gpu mp --no-mpi -j 8
```

- `--gpu mp`: OpenMP target offload backend (SPIR64)
- `--no-mpi`: omit for MPI-enabled runs; include for single-node testing
- `-j 8`: parallel build jobs

MFC will automatically:
1. Compile `$MKLROOT/include/mkl_dfti_omp_offload.f90` with minimal flags
   (no `-free -fpp`) via a CMake `add_custom_command` to avoid OpenMP 5.2
   clause compatibility issues with global compile flags
2. Link `-qmkl=parallel` for MKL threading + core
3. Link `libmkl_sycl_dft`, `libsycl`, `libOpenCL` for GPU FFT dispatch

## GPU FFT implementation

MFC uses oneMKL DFTI with the OpenMP 5.1 `!$omp dispatch` construct for FFT
in cylindrical geometry (the azimuthal Fourier filter in `m_fftw.fpp`).
This is activated when `__INTEL_LLVM_COMPILER` is defined (i.e., compiled with ifx).

Key verified properties (oneMKL 2026.0, ifx 2025.3.3):
- Batch R2C transform with `INPUT_DISTANCE != OUTPUT_DISTANCE` works correctly.
  MFC uses `real_size = p+1` and `cmplx_size = (p+1)/2+1` which always differ.
- `!$omp dispatch` correctly dispatches DFTI calls to device-mapped allocatables.

## Running MFC cases

Build all three targets (pre_process, simulation, post_process) before running:

```bash
./mfc.sh build --gpu mp --no-mpi -j 8
```

Then run a case normally:

```bash
./mfc.sh run examples/1D_convergence/case.py --no-build --no-mpi
```

To run individual stages directly (useful when `syscheck` blocks due to GPU access):

```bash
export MKLROOT=/opt/intel/oneapi/mkl/2026.0
export LD_LIBRARY_PATH=/opt/intel/oneapi/umf/1.1/lib:\
/opt/intel/oneapi/compiler/2026.0/lib:\
/opt/intel/oneapi/compiler/2025.3/lib:\
${MKLROOT}/lib:\
/opt/intel/oneapi/tbb/2022.1/lib:\
$LD_LIBRARY_PATH

cd examples/my_case
/path/to/build/install/<hash>/bin/pre_process
/path/to/build/install/<hash>/bin/simulation
```

The install hashes are printed by `./mfc.sh build`; look for lines like
`✓ Installed simulation`.

## Chemistry/pyrometheus cases

MFC's pyrometheus-generated thermochemistry (`m_thermochem.f90`) works with the
Intel GPU build. Pyrometheus emits `!$omp declare target` for all thermochem
routines when `directive_offload="mp"` is requested (automatically set when
`--gpu mp` is active).

Verified on ifx 2025.3.3: `m_thermochem.f90` and `m_chemistry.fpp` both compile
at O3 + SPIR64 without ICE. The `1D_reactive_shocktube` example (H2/O2/Ar, 29
reactions, 10 species) runs to completion. No Intel-specific source workarounds
are needed for chemistry beyond the general levelset fix.

To build and run a reactive chemistry case:
```bash
# Build (chemistry module generated automatically from cantera_file in case.py)
./mfc.sh run examples/1D_reactive_shocktube/case.py \
    --gpu mp --no-mpi -t pre_process simulation

# Or bypass syscheck if no GPU render-group access:
cd examples/1D_reactive_shocktube
/path/to/build/install/<pre_hash>/bin/pre_process
OMP_TARGET_OFFLOAD=DISABLED \
    /path/to/build/install/<sim_hash>/bin/simulation
```

ifx warning during chemistry compilation:
```
warning #8694: A procedure called by a procedure with the DECLARE TARGET
attribute must have the DECLARE TARGET attribute.  [GET_MIXTURE_VISCOSITY_MIXAVG]
```
This is a false positive from ifx's module-interface tracking; `m_thermochem.f90`
does declare all routines as target. The warning is harmless and the code runs
correctly.

## GPU device access

The Intel GPU requires membership in the `render` group (GID 109) to access
`/dev/dri/renderD128` via Level Zero.

Without render group access, `ZE_RESULT_ERROR_UNINITIALIZED` is returned by
Level Zero. OpenMP target offload falls back to the CPU host plugin
(correct results, but no GPU acceleration).

To diagnose GPU visibility:

```bash
ls -la /dev/dri/renderD128          # should show rw permissions for your user/group=render

# With full LD_LIBRARY_PATH set:
LD_LIBRARY_PATH=/opt/intel/oneapi/umf/1.1/lib:... \
    /opt/intel/oneapi/compiler/2026.0/bin/sycl-ls --verbose
# Look for: "[opencl:gpu]" or "[ext_oneapi_level_zero:gpu]" platforms

LIBOMPTARGET_DEBUG=1 ./simulation   # look for "Level0 NG plugin initialization"
                                    # and absence of "ZE_RESULT_ERROR_UNINITIALIZED"
```

To get GPU access:
- **Interactive shell**: request from system admin to add user to `render` group
  (`sudo usermod -a -G render $USER`, then re-login)
- **SLURM**: submit with `--gres=gpu:max_1100=1`; if Level Zero still fails,
  the SLURM epilog may not have configured device cgroup ACLs for the job user —
  contact the system admin

> **Note on `sycl-ls`**: the 2026.0 `sycl-ls` requires `libumf.so.1` from
> `/opt/intel/oneapi/umf/1.1/lib` in `LD_LIBRARY_PATH`, otherwise all adapters
> fail to load and it reports "No platforms found".

## Link flags (what MFC's CMake generates)

The full set of flags the compiler uses for the simulation target:

**Compile flags:**
```
-fiopenmp -fopenmp-targets=spir64 -free -fpp -march=native
```

**Link flags:**
```
-fiopenmp -fopenmp-targets=spir64
-qmkl=parallel
-L$MKLROOT/lib -lmkl_sycl_dft
-L/opt/intel/oneapi/compiler/2026.0/lib -lsycl -lOpenCL
```

**MKL OMP module (compiled separately, no global flags):**
```bash
ifx -fiopenmp -fopenmp-targets=spir64 \
    -c -I$MKLROOT/include \
    $MKLROOT/include/mkl_dfti_omp_offload.f90 \
    -o mkl_dfti_omp_offload.o
```

## Known issues

### `need_device_addr` compilation error
`mkl_dfti_omp_offload.f90` from MKL 2026.0 uses `need_device_addr` in
`!$omp declare variant` (OpenMP 5.2). This requires ifx **2025.3** or newer.
If CMake finds an older ifx (e.g., 2025.0 from a system module path), the
compile fails with:
```
error #5082: Syntax error, found IDENTIFIER 'NEED_DEVICE_ADDR'
```
Fix: set `FC=/opt/intel/oneapi/compiler/2025.3/bin/ifx` before building
and run `./mfc.sh clean` first so CMake re-detects the compiler.

### Two routines with ifx SPIR64 codegen bugs

**`s_apply_levelset` (`m_compute_levelset.fpp`)** — ifx SPIR64 inliner ICE:

The LLVM inliner (at O1+) pulls `!$omp declare target (seq)` geometry
routines into the `target teams loop` kernel and generates LLVM IR that
crashes the SPIR-V converter with a segfault in `llvm-spirv`. At O0 the
crash does not occur (no inlining). Two fixes combined:

1. **Split loops**: replaced the single if-else dispatch loop with one
   `GPU_PARALLEL_LOOP` per geometry type so each kernel calls exactly one
   declare-target routine. The original multi-callee dispatch also triggers
   `"Instruction does not dominate all uses!"` in llvm-link.

2. **Per-file `-fno-inline`**: in `CMakeLists.txt`, `set_source_files_properties`
   adds `-fno-inline` to `m_compute_levelset.fpp.f90` for IntelLLVM+OpenMP
   builds. This prevents the inliner from pulling declare-target routines into
   the kernel body where they crash the SPIR-V backend. The routines remain
   callable as proper device-side function calls via `!$omp declare target`.

**`s_pressure_relaxation_procedure` (`m_pressure_relaxation.fpp`)** — SPIR-V
InvalidArraySize in declare-target helpers:

`!$omp declare target (seq)` routines with `dimension(sys_size)` explicit-shape
dummy arguments trigger `"InvalidArraySize: Array size must be at least 1"` in
llvm-spirv. SPIR-V requires compile-time constant array bounds; `sys_size` is
a runtime module integer. Fixed by changing `dimension(sys_size)` →
`dimension(:)` (assumed-shape) on all helper routine interfaces. The loop now
runs on GPU for Intel.

### syscheck GPU assertion
`syscheck` runs `assert(omp_get_num_devices() > 0)` and aborts if the GPU
is not accessible. This is a runtime check, not a build issue. See GPU device
access section above.

To run a case anyway (testing code correctness on CPU fallback), invoke
`pre_process` and `simulation` directly from their install paths, bypassing
the `./mfc.sh run` wrapper that calls `syscheck` first.

### `libumf.so.1` not found at runtime
The 2026.0 Level Zero and OpenCL UR adapters link against `libumf.so.1`.
If not in `LD_LIBRARY_PATH`, all adapters fail silently and sycl-ls reports
"No platforms found". Fix:

```bash
export LD_LIBRARY_PATH=/opt/intel/oneapi/umf/1.1/lib:$LD_LIBRARY_PATH
```
