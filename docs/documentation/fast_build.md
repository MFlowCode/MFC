# Fast Dev Builds (`--fast-build`)

Status: **prototype / work in progress.** The NVHPC path is implemented and measured;
the AMD (LLVMFlang) path described below is an analysis with a proposed change that
still needs validation on an AMD GPU + AMD-compiler machine.

## Motivation

GPU builds are slow to iterate on when you just want to add a `print` and re-run.
Two different compilers, two different bottlenecks:

- **NVHPC** — the per-iteration cost is the two-pass IPO link (`-Mextract` + `-Minline`),
  plus device codegen for every targeted compute capability. Editing one hot file and
  rebuilding takes minutes because the IPO passes re-run and device code is generated for
  every arch in `MFC_CUDA_CC`.
- **AMD / LLVMFlang** — the per-iteration cost is the **device LTO link**. The link step
  can take 20+ minutes *every* build, because the OpenMP offload device link does
  whole-program LTO regardless of what changed.

`--fast-build` is a dedicated build mode that strips the expensive, optimization-oriented
machinery that is pointless during print-debugging.

## Usage

```bash
./mfc.sh build -t simulation --gpu acc --fast-build -j 8     # NVHPC (OpenACC)
./mfc.sh build -t simulation --gpu mp  --fast-build -j 8     # AMD/Cray (OpenMP offload)
```

`--fast-build` is mutually exclusive with `--debug` / `--reldebug`. It is **not** a
correctness build: no bounds checking, no `MFC_DEBUG` asserts. Add your own `print`/`write`
statements; pair with `--debug` when you need runtime checks.

## What it does

`--fast-build` selects a new CMake build type, `Fast`, that deliberately matches **none** of
the existing conditional flag blocks in `CMakeLists.txt`:

- Not `Release`, so no IPO/LTO and no `-march=native`.
- Not `Debug`/`RelDebug`, so no `MFC_DEBUG` and no `-gpu=debug`.

It then adds a light `-O1` (via `add_compile_options`, since the `CMAKE_*_FLAGS_FAST` cache
variables do not inject flags in this codebase). Because `MFC_DEBUG` is off, device routines
contain no host-only debug aborts, so the binary compiles cleanly **without** IPO.

On NVHPC GPU builds it also restricts device codegen to a **single** compute capability —
the GPU on the build node, detected via `nvidia-smi` — overriding the multi-arch
`MFC_CUDA_CC` that the module files set. Set `MFC_FAST_ARCH=<cc>` (e.g. `MFC_FAST_ARCH=90`)
to override the detection on a login node with no visible GPU.

## NVHPC results (measured)

NVHPC 24.5, Quadro RTX 6000 (cc75), generic `simulation` build, 8 cores:

| Scenario | Release (fat 5-arch) | `--fast-build` (single-arch) |
| --- | --- | --- |
| Clean full build | 641 s | 170 s (3.8x) |
| Hot-module incremental (`m_riemann_solvers`) | 385 s | 79 s (4.9x) |

Verified: builds with no IPO (`-Mextract` absent), no `MFC_DEBUG`, single `-gpu=cc75`, and
the resulting binary runs a 1D case on the GPU to exit code 0 with finite output.

## AMD / LLVMFlang: the device-LTO link (proposed, needs validation)

The AMD GPU offload flags live in `CMakeLists.txt` (`MFC_SETUP_TARGET`):

```cmake
# compile
target_compile_options(${a_target} PRIVATE
    -fopenmp --offload-arch=gfx90a -O3
    -fopenmp-assume-threads-oversubscription
    -fopenmp-assume-teams-oversubscription)
# link
target_link_options(${a_target} PRIVATE
    -fopenmp --offload-arch=gfx90a -flto-partitions=${MFC_BUILD_JOBS})
```

The `-flto-partitions` at link is the tell: the OpenMP offload **device link runs
whole-program LTO every time**, so even a one-file edit re-LTOs all device code. Single-arch
is not a lever here (already a single `gfx90a`).

### Levers, best first

1. **JIT the device code: `-fopenmp-target-jit`.** Instead of AOT-compiling and LTO-linking
   device code into the binary at link time, embed device LLVM-IR and JIT each kernel at
   runtime on first launch. The device LTO link essentially disappears, so the link drops to
   roughly host-link time. Cost: a one-time JIT warmup on the *run* (tunable with
   `LIBOMPTARGET_JIT_OPT_LEVEL`). This is the real fix for AMD iteration.

2. **Build with a high `-j` (no code change — try this first).** `-flto-partitions` is set to
   `MFC_BUILD_JOBS`, which is your `-j`. Building with `-j 8` on a 64-core node runs the
   device LTO link only 8-way parallel. Use `-j 32`/`-j 64` to give the LTO link more
   partitions; this alone may cut the link time substantially with no toolchain change.

3. **Lower device optimization `-O3` -> `-O1`/`-O0` for dev builds.** The `-O3` drives the
   heavy LTO optimization; lowering it cuts link work (slower runtime, fine for debugging).

4. **`-fno-lto` (AOT, non-LTO device link).** Links per-translation-unit device objects
   instead of whole-program LTO. Potentially faster, but less certain across ROCm/flang
   versions — only if JIT does not pan out.

### Proposed `Fast` branch for LLVMFlang

Make the flags above build-type-aware so `--fast-build` emits, for `LLVMFlang`:

```cmake
# compile
-fopenmp --offload-arch=gfx90a -O1 -fopenmp-target-jit \
    -fopenmp-assume-threads-oversubscription -fopenmp-assume-teams-oversubscription
# link (no -flto-partitions; JIT removes the whole-program device LTO)
-fopenmp --offload-arch=gfx90a -fopenmp-target-jit
```

### How to validate on an AMD machine

On a Frontier AMD / AFAR-style node (`source ./mfc.sh load -c famd -m g` or equivalent):

1. **Baseline** — time the current link:
   `./mfc.sh build -t simulation --gpu mp -j 8` and note the link duration.
2. **Free lever** — rebuild with a high `-j` (more LTO partitions) and compare the link time:
   `./mfc.sh build -t simulation --gpu mp -j 64`.
3. **JIT lever** — once the `Fast` LLVMFlang branch is wired in, build with `--fast-build`
   and confirm: (a) the link time collapses, (b) a small case runs to exit 0 (expect a
   one-time JIT warmup on first launch). `OMP_TARGET_OFFLOAD=MANDATORY` is already set.

### Caveat

None of the AMD numbers are measured — LLVMFlang is not available on the development machine
used so far. The diagnosis follows directly from the build flags and from how LLVM OpenMP
offload works, but the exact flag spelling and ROCm/flang-version behavior must be confirmed
on real hardware before the LLVMFlang `Fast` branch is trusted.

## Implementation notes / TODO

- Implemented: `Fast` CMake build type, `fast_build` field in `MFCConfig` (auto
  `--fast-build`/`--no-fast-build`, own build slug), NVHPC single-arch autodetect, lock-file
  version bump for the new config field.
- Not yet: a `--gpu-arch` CLI flag (only the `MFC_FAST_ARCH` env escape hatch exists), the
  LLVMFlang `Fast` branch above, Cray-on-AMD validation, and `--help`/docs polish.
