# GPU and MPI Patterns

## GPU Offloading Architecture

Only `src/simulation/` is GPU-accelerated. Pre/post_process run on CPU only.

MFC uses a **backend-agnostic GPU abstraction** via Fypp macros. The same source code
compiles to either OpenACC or OpenMP target offload depending on the build flag:

- `./mfc.sh build --gpu acc` → OpenACC backend (NVIDIA nvfortran, Cray ftn)
- `./mfc.sh build --gpu mp`  → OpenMP target offload backend (Cray ftn, AMD flang)
- `./mfc.sh build` (no --gpu) → CPU-only, GPU macros expand to plain Fortran

### Macro Layers (in src/common/include/)
- `parallel_macros.fpp` — **Use these.** Generic `GPU_*` macros that dispatch to the
  correct backend based on `MFC_OpenACC` / `MFC_OpenMP` compile definitions.
- `acc_macros.fpp` — OpenACC-specific `ACC_*` implementations (do not call directly)
- `omp_macros.fpp` — OpenMP target offload `OMP_*` implementations (do not call directly)
  - OMP macros generate **compiler-specific** directives: NVIDIA uses `target teams loop`,
    Cray uses `target teams distribute parallel do simd`, AMD uses
    `target teams distribute parallel do`
- `shared_parallel_macros.fpp` — Shared helpers (collapse, private, reduction generators)

### Key GPU Macros (always use the `GPU_*` prefix)

Full set with signatures in `parallel_macros.fpp`. The ones you reach for most:
- `$:GPU_PARALLEL_LOOP(collapse=N, private=[...], reduction=[...], reductionOp='+')`
  + `$:END_GPU_PARALLEL_LOOP()` — parallel spatial loop; by far the most common (see pattern below).
- `$:GPU_LOOP(collapse=N, ...)` — inner loop *within* a parallel region.
- `$:GPU_UPDATE(host=[...])` / `$:GPU_UPDATE(device=[...])` — device↔host copies (around MPI; see below).
- `#:call GPU_PARALLEL(...)` — block region for scalar reductions (`maxval`/`minval`).

Others in `parallel_macros.fpp`: `GPU_ENTER_DATA`/`GPU_EXIT_DATA`, `GPU_DECLARE`, `GPU_ROUTINE`,
`GPU_ATOMIC`, `GPU_WAIT`, and the block macros `GPU_DATA`, `GPU_HOST_DATA`.

Typical GPU loop pattern (used 750+ times in the codebase):
```
$:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=3)
do l = idwbuff(3)%beg, idwbuff(3)%end
    do k = idwbuff(2)%beg, idwbuff(2)%end
        do j = idwbuff(1)%beg, idwbuff(1)%end
            ! loop body
        end do
    end do
end do
$:END_GPU_PARALLEL_LOOP()
```

WARNING: Do NOT use `GPU_PARALLEL` wrapping `GPU_LOOP` for spatial loops. `GPU_LOOP`
emits empty directives on Cray and AMD compilers, causing silent serial execution.
Use `GPU_PARALLEL_LOOP` / `END_GPU_PARALLEL_LOOP` for all parallel spatial loops.

NEVER write raw `!$acc` or `!$omp` directives. Always use `GPU_*` Fypp macros.
The precheck source lint will catch raw directives and fail.

### Memory Management Macros (from macros.fpp)
- `@:ALLOCATE(var1, var2, ...)` — Fortran allocate + `GPU_ENTER_DATA(create=...)`
- `@:DEALLOCATE(var1, var2, ...)` — `GPU_EXIT_DATA(delete=...)` + Fortran deallocate
- `@:PREFER_GPU(var1, var2, ...)` — NVIDIA unified memory page placement hint
- Every `@:ALLOCATE` MUST have a matching `@:DEALLOCATE` in finalization
- Conditional allocation MUST have conditional deallocation

### GPU Field Setup (Cray-specific, from macros.fpp)
- `@:ACC_SETUP_VFs(...)` / `@:ACC_SETUP_SFs(...)` — GPU pointer setup for vector/scalar fields
- These compile only for Cray (`_CRAYFTN`); other compilers skip them

### Compiler-Backend Matrix

CI-gated compilers (must always pass): gfortran, nvfortran, Cray ftn, Intel ifx.
AMD flang is additionally supported for GPU builds but not in the CI matrix.

| Compiler        | `--gpu acc` (OpenACC) | `--gpu mp` (OpenMP)    | CPU-only |
|-----------------|----------------------|------------------------|----------|
| GNU gfortran    | No                   | Experimental (AMD GCN) | Yes      |
| NVIDIA nvfortran| Yes (primary)        | Yes                    | Yes      |
| Cray ftn (CCE)  | Yes                  | Yes (primary)          | Yes      |
| Intel ifx       | No                   | Experimental (SPIR64)  | Yes      |
| AMD flang       | No                   | Yes                    | Yes      |

## Preprocessor Defines (`#ifdef` / `#ifndef`)

Raw `#ifdef` / `#ifndef` preprocessor guards are **normal and expected** in MFC.
They are NOT the same as raw `!$acc`/`!$omp` pragmas (which are forbidden).

Use `#ifdef` for feature, target, compiler, and library gating:

### Feature gating
- `MFC_MPI` — MPI-enabled build (`--mpi` flag, default ON)
- `MFC_OpenACC` — OpenACC GPU backend (`--gpu acc`)
- `MFC_OpenMP` — OpenMP target offload backend (`--gpu mp`)
- `MFC_GPU` — Any GPU build (either OpenACC or OpenMP)
- `MFC_DEBUG` — Debug build (`--debug`)
- `MFC_SINGLE_PRECISION` — Single-precision mode (`--single`)
- `MFC_MIXED_PRECISION` — Mixed-precision mode (`--mixed`)

### Target gating (for code in `src/common/` shared across executables)
- `MFC_PRE_PROCESS` — Only in pre_process builds
- `MFC_SIMULATION` — Only in simulation builds
- `MFC_POST_PROCESS` — Only in post_process builds

### Compiler gating (for compiler-specific workarounds)
Compiler/feature macros: `_CRAYFTN`, `__NVCOMPILER_GPU_UNIFIED_MEM` (NVIDIA unified mem, GH-200 /
`--unified`), `__PGI` (legacy PGI/NVIDIA), `__INTEL_COMPILER`, `FRONTIER_UNIFIED`. Library code is
similarly gated (FFTW in `m_fftw.fpp` on `MFC_GPU`/`__PGI`; CUDA Fortran `cudafor` on
`__NVCOMPILER_GPU_UNIFIED_MEM`; SILO/HDF5 paths). Grep the relevant file for exact usage.

When adding new `#ifdef` blocks, always provide an `#else` or `#endif` path so
the code compiles in all configurations (CPU-only, GPU-ACC, GPU-OMP, with/without MPI).

## MPI

### Halo Exchange
- Pack/unpack offset calculations are error-prone — verify carefully
- Buffer sizing depends on dimensionality and QBMM state
- GPU coherence: always `GPU_UPDATE(host=...)` before MPI send,
  `GPU_UPDATE(device=...)` after MPI receive

### Error Handling
- Use `call s_mpi_abort()` for fatal errors, never `stop` or `error stop`
- MPI must be finalized before program exit
