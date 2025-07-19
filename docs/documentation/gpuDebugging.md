# Debugging Tools and Tips for GPUs

## Compiler agnostic tools

## OpenMP tools
```bash
OMP_DISPLAY_ENV=true | false | verbose
```
- Prints out the internal control values and environment variables at beginning of program if `true` or `verbose`
- `verbose` will also print out vendor-specific internal control values and environment variables

```bash
OMP_TARGET_OFFLOAD = MANDATORY | DISABLED | DEFAULT
```
- Quick way to turn off off-load (DISABLED) or make it abort if a GPU isn't found (MANDATORY)
- great first test: does the problem disappear when you drop back to the CPU? 

```bash
OMP_THREAD_LIMIT=<positive_integer>
```
- Sets maximum number of OpenMP threads to use in a contention group
- Might be useful in checking for issues with contention or race conditions

```bash
OMP_DISPLAY_AFFINITY=TRUE
```
- Will display affinity bindings for each OpenMP thread, containing hostname, process identifier, OS thread identifier, OpenMP thread identifier, and affinity binding.

## Cray Compiler Tools

### Cray General Options

```bash
CRAY_ACC_DEBUG: 0 (off), 1, 2, 3 (very noisy)
```
- Dumps a time-stamped log line ("ACC: …) for every allocation, data transfer, kernel launch, wait, etc. Great first stop when "nothing seems to run on the GPU.

- Outputs on STDERR by default. Can be changed by setting `CRAY_ACC_DEBUG_FILE`. 
  - Recognizes `stderr`, `stdout`, and `process`.
  - `process` automatically generates a new file based on `pid` (each MPI process will have a different file)

- While this environment variable specifies ACC, it can be used for both OpenACC and OpenMP

```bash
CRAY_ACC_FORCE_EARLY_INIT=1
```
- Force full GPU initialization at program start so you can see start-up hangs immediately
- Default behavior without environment variable is to defer initialization on first use
-  Device initialization includes initializing the GPU vendor’s low-level device runtime library (e.g., libcuda for NVIDIA GPUs) and establishing all necessary software contexts for interacting with the device

### Cray OpenACC Options

```bash
CRAY_ACC_PRESENT_DUMP_SAVE_NAMES=1
```
- Will cause acc_present_dump() to output variable names and file locations in addition to variable mappings
- Add acc_present_dump() around hotspots to help find problems with data movements
  - Helps more if adding `CRAY_ACC_DEBUG` environment variable

## NVHPC Compiler Options

### NVHPC General Options

```bash
STATIC_RANDOM_SEED=1
```
- Forces the seed returned by RANDOM_SEED to be constant, so generates same sequence of random numbers 
- Useful for testing issues with randomized data

```bash
NVCOMPILER_TERM=option[,option]
```
- `[no]debug`: Enables/disables just-in-time debugging (debugging invoked on error)
- `[no]trace`: Enables/disables stack traceback on error

### NVHPC OpenACC Options

```bash
NVCOMPILER_ACC_NOTIFY= <bitmask>
```
- Assign the environment variable to a bitmask to print out information to stderr for the following
  - kernel launches: 1
  - data transfers: 2
  - region entry/exit: 4
  - wait operation of synchronizations with the device: 8
  - device memory allocations and deallocations: 16
- 1 (kernels only) is the usual first step.3 (kernels + copies) is great for "why is it so slow?" 

```bash
NVCOMPILER_ACC_TIME=1
```
- Lightweight profiler
- prints a tidy end-of-run table with per-region and per-kernel times and bytes moved
- Do not use with CUDA profiler at the same time

```bash
NVCOMPILER_ACC_DEBUG=1
```
- Spews everything the runtime sees: host/device addresses, mapping events, present-table look-ups, etc.
- Great for "partially present" or "pointer went missing" errors. 
- [Doc for NVCOMPILER_ACC_DEBUG](https://docs.nvidia.com/hpc-sdk/archive/20.9/pdf/hpc209openacc_gs.pdf)
  - Ctrl+F for `NVCOMPILER_ACC_DEBUG`

### NVHPC OpenMP Options

```bash
LIBOMPTARGET_PROFILE=run.json
```
- Emits a Chrome-trace (JSON) timeline you can open in chrome://tracing or Speedscope
- great lightweight profiler when Nsight is over-kill.
- Granularity in µs via `LIBOMPTARGET_PROFILE_GRANULARITY` (default 500).

```bash
LIBOMPTARGET_INFO=<bitmask>
```
- Prints out different types of runtime information
- Human-readable log of data-mapping inserts/updates, kernel launches, copies, waits.
- Perfect first stop for "why is nothing copied?" 
- Flags
  - Print all data arguments upon entering an OpenMP device kernel: 0x01
  - Indicate when a mapped address already exists in the device mapping table: 0x02
  - Dump the contents of the device pointer map at kernel exit: 0x04
  - Indicate when an entry is changed in the device mapping table: 0x08
  - Print OpenMP kernel information from device plugins: 0x10
  - Indicate when data is copied to and from the device: 0x20

```bash
LIBOMPTARGET_DEBUG=1
```
- Developer-level trace (host-side)
- Much noisier than INFO
- only works if the runtime was built with -DOMPTARGET_DEBUG.

```bash
LIBOMPTARGET_JIT_OPT_LEVEL=-O{0,1,2,3}
```
- This environment variable can be used to change the optimization pipeline used to optimize the embedded device code as part of the device JIT. 
- The value is corresponds to the -O{0,1,2,3} command line argument passed to clang.

```bash
LIBOMPTARGET_JIT_SKIP_OPT=1
```
- This environment variable can be used to skip the optimization pipeline during JIT compilation.
- If set, the image will only be passed through the backend.
- The backend is invoked with the LIBOMPTARGET_JIT_OPT_LEVEL flag.

## Compiler Documentation
- [Cray & OpenMP Docs](https://cpe.ext.hpe.com/docs/24.11/cce/man7/intro_openmp.7.html#environment-variables)
- [Cray & OpenACC Docs](https://cpe.ext.hpe.com/docs/24.11/cce/man7/intro_openacc.7.html#environment-variables)
- [NVHPC & OpenACC Docs](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html?highlight=NVCOMPILER_#environment-variables)
- [NVHPC & OpenMP Docs](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html?highlight=NVCOMPILER_#id2)
- [LLVM & OpenMP Docs] (https://openmp.llvm.org/design/Runtimes.html)
    - NVHPC is built on top of LLVM
- [OpenMP Docs](https://www.openmp.org/spec-html/5.1/openmp.html)
- [OpenACC Docs](https://www.openacc.org/sites/default/files/inline-files/OpenACC.2.7.pdf)