@page troubleshooting Troubleshooting Guide

# Troubleshooting Guide

This guide covers debugging tools, common issues, and troubleshooting workflows for MFC.

## Quick Reference

### Debugging Flags

| Flag | Command | Purpose |
|------|---------|---------|
| `-v` | build, run, test | Show commands being executed |
| `-vv` | build, run, test | Full compiler/cmake output |
| `-vvv` | build | Add cmake dependency debugging |
| `-d` | all | Write debug log to file |
| `--debug` | build | Build with debug symbols |
| `--gcov` | build | Build with code coverage |
| `--no-gpu` | build | Disable GPU to isolate issues |
| `--no-mpi` | build | Disable MPI to isolate issues |

### Profiling Flags

| Flag | Command | Purpose |
|------|---------|---------|
| `--ncu` | run | NVIDIA Nsight Compute (kernel profiling) |
| `--nsys` | run | NVIDIA Nsight Systems (timeline profiling) |
| `--rcu` | run | AMD rocprof-compute |
| `--rsys` | run | AMD rocprof |

### Useful Commands

```bash
./mfc.sh validate case.py      # Check case for errors before running
./mfc.sh build --debug         # Build with debug symbols
./mfc.sh build -vv             # Build with full compiler output
./mfc.sh run case.py -v        # Run with verbose output
./mfc.sh test --only <UUID>    # Run a specific test
./mfc.sh clean                 # Clean and start fresh
```

---

## Debugging Workflow

### Build Fails?

```
1. Run with -vv to see full error output
   ./mfc.sh build -vv

2. If dependency issue, try -vvv for cmake debug
   ./mfc.sh build -vvv

3. Try disabling features to isolate:
   ./mfc.sh build --no-gpu      # Rule out GPU issues
   ./mfc.sh build --no-mpi      # Rule out MPI issues

4. Clean and rebuild:
   ./mfc.sh clean && ./mfc.sh build
```

### Run Fails?

```
1. Validate your case first:
   ./mfc.sh validate case.py

2. Run with verbose output:
   ./mfc.sh run case.py -v

3. Check for NaN or constraint errors in output

4. Try with fewer MPI ranks:
   ./mfc.sh run case.py -n 1

5. Build in debug mode for better error messages:
   ./mfc.sh build --debug
   ./mfc.sh run case.py
```

### Test Fails?

```
1. Run the specific failing test:
   ./mfc.sh test --only <UUID>

2. Run serially for clearer output:
   ./mfc.sh test --only <UUID> -j 1

3. Check test output in:
   build/tests/<UUID>/

4. If changes are intentional, regenerate golden files:
   ./mfc.sh test --generate --only <UUID>
```

---

## Build Debugging

### Verbosity Levels

MFC supports three verbosity levels for builds:

**`-v` (Level 1):** Shows build progress with file counts
```bash
./mfc.sh build -v
# Shows: [1/42] Compiling m_global_parameters.fpp.f90
```

**`-vv` (Level 2):** Full compiler commands and cmake output
```bash
./mfc.sh build -vv
# Shows: Full compiler invocations with all flags
# Useful for: Seeing exact compiler errors, checking flags
```

**`-vvv` (Level 3):** Adds cmake dependency debugging
```bash
./mfc.sh build -vvv
# Shows: cmake --debug-find output
# Useful for: Finding why cmake can't locate MPI, HDF5, etc.
```

### Debug Builds

Build with debug symbols for better error messages and debugging:

```bash
./mfc.sh build --debug
```

This sets `CMAKE_BUILD_TYPE=Debug`, which:
- Adds `-g` debug symbols
- Reduces optimization (`-O0` or `-O1`)
- Enables runtime checks in some compilers
- Makes stack traces more readable

### Code Coverage

Build with code coverage instrumentation:

```bash
./mfc.sh build --gcov
```

Useful for identifying which code paths are exercised by tests.

### Isolating Issues

When builds fail, isolate the problem by disabling features:

```bash
# Rule out GPU compiler issues
./mfc.sh build --no-gpu

# Rule out MPI issues
./mfc.sh build --no-mpi

# Rule out precision issues
./mfc.sh build --single    # Single precision
./mfc.sh build --mixed     # Mixed precision
```

### Clean Rebuild

Sometimes a clean rebuild fixes issues:

```bash
./mfc.sh clean
./mfc.sh build -j $(nproc)
```

Or manually:
```bash
rm -rf build/
./mfc.sh build -j $(nproc)
```

---

## Runtime Debugging

### Case Validation

**Always validate your case before running:**

```bash
./mfc.sh validate case.py
```

This checks for:
- Syntax errors in case.py
- Invalid parameter values
- Constraint violations (incompatible parameter combinations)
- Typos in parameter names (with "did you mean?" suggestions)

**Example output:**
```
Errors:
  [red] weno_order=7 but must be one of: [1, 3, 5]
  [red] Unknown parameter 'weno_ordr' - did you mean 'weno_order'?

Warnings:
  [yellow] bubbles=True recommends setting: nb, polytropic
```

For more detailed validation output:
```bash
./mfc.sh validate case.py -d
```

### Verbose Runs

Add verbosity to see what's happening:

```bash
# Show the exact command being run
./mfc.sh run case.py -v

# Show job script details
./mfc.sh run case.py -vv
```

### Debug Log

Write detailed debug information to a log file:

```bash
./mfc.sh run case.py -d
```

This creates a debug log with detailed internal state, useful for reporting issues.

### Running with Fewer Ranks

To isolate parallel issues, try running with fewer MPI ranks:

```bash
./mfc.sh run case.py -n 1    # Single rank
./mfc.sh run case.py -n 2    # Two ranks
```

---

## Test Debugging

### Running Specific Tests

Run only a specific test by UUID:

```bash
./mfc.sh test --only 3D_sphbubcollapse
```

Run multiple specific tests:
```bash
./mfc.sh test --only 3D_sphbubcollapse 2D_shockbubble
```

### Serial Test Execution

Run tests serially for clearer output:

```bash
./mfc.sh test -j 1
```

### Test Output Location

Test outputs are written to:
```
build/tests/<UUID>/
```

This directory contains:
- Input files generated from case.py
- Output files from each stage (pre_process, simulation, post_process)
- Any error messages or logs

### Regenerating Golden Files

If you intentionally changed physics/numerics, regenerate reference files:

```bash
# Regenerate for specific tests
./mfc.sh test --generate --only <UUID>

# Regenerate for a range of tests
./mfc.sh test --generate --from <UUID1> --to <UUID2>
```

### Understanding Test Failures

**"Golden file mismatch"** means numerical results differ from reference values.

Possible causes:
1. **Intentional changes:** Regenerate golden files with `--generate`
2. **Compiler differences:** Different compilers produce slightly different results
3. **Precision settings:** Single vs double precision
4. **Platform differences:** CPU architecture, GPU differences

---

## Performance Profiling

### NVIDIA GPU Profiling

**Nsight Compute** (kernel-level analysis):
```bash
./mfc.sh run case.py --ncu
# Or with additional ncu flags:
./mfc.sh run case.py --ncu --set full
```

**Nsight Systems** (timeline/system analysis):
```bash
./mfc.sh run case.py --nsys
# Or with additional nsys flags:
./mfc.sh run case.py --nsys -o profile_output
```

### AMD GPU Profiling

**rocprof-compute** (kernel analysis):
```bash
./mfc.sh run case.py --rcu
```

**rocprof** (system profiling):
```bash
./mfc.sh run case.py --rsys
```

### Performance vs Correctness

If results seem wrong, first verify correctness:
1. Run `./mfc.sh validate case.py`
2. Compare against a known working case
3. Try `--debug` build for better error detection

If results are correct but slow:
1. Use profiling tools above to identify bottlenecks
2. Try `--case-optimization` for production runs (10x speedup)
3. Check GPU utilization with `nvidia-smi` or `rocm-smi`

---

## Common Build Errors

### "CMake could not find MPI"

**Cause:** MPI is not installed or not in PATH.

**Fix:**
- **Ubuntu/Debian:** `sudo apt install libopenmpi-dev openmpi-bin`
- **macOS (Homebrew):** `brew install open-mpi`
- **HPC systems:** `module load openmpi` or use `./mfc.sh load`

Verify: `mpirun --version`

### "CMake could not find a Fortran compiler"

**Cause:** No Fortran compiler installed or not in PATH.

**Fix:**
- **Ubuntu/Debian:** `sudo apt install gfortran`
- **macOS (Homebrew):** `brew install gcc`
- **HPC systems:** `module load gcc` or `module load nvhpc`

### "Fypp preprocessing failed"

**Cause:** Syntax error in `.fpp` files or missing Fypp.

**Fix:**
1. Ensure Fypp is installed: `pip install fypp`
2. Check the specific error line in the output
3. Run `./mfc.sh format` to check for formatting issues

### Build fails with GPU/OpenACC errors

**Cause:** GPU compiler not properly configured.

**Fix:**
1. Verify GPU compiler is available:
   - NVIDIA: `nvfortran --version`
   - AMD: `ftn --version` or `amdflang --version`
2. Load appropriate modules: `module load nvhpc` or `module load rocm`
3. Set CUDA compute capability if needed: `export MFC_CUDA_CC=80`
4. Try building without GPU first: `./mfc.sh build --no-gpu`

### "Module not found" on HPC

**Cause:** Required modules not loaded.

**Fix:**
1. Use MFC's module loader if available for your system:
   ```bash
   source ./mfc.sh load -c <cluster> -m <mode>
   ```
2. Or manually load modules:
   ```bash
   module load gcc openmpi cmake python
   ```

---

## Common Runtime Errors

### "Case parameter constraint violations"

**Cause:** Invalid combination of simulation parameters.

**Fix:**
1. Run `./mfc.sh validate case.py` for detailed diagnostics
2. Check the specific constraint mentioned
3. Review similar examples in `examples/`

Common issues:
- Grid dimensions (`m`, `n`, `p`) not matching dimensionality
- `weno_order` too high for grid size
- Boundary conditions not matching domain setup

### "NaN detected"

**Cause:** Numerical instability.

**Fix:**
1. Reduce time step (`dt`)
2. Check initial conditions for unphysical values
3. Verify boundary conditions are appropriate
4. Try a more diffusive scheme initially
5. Check CFL condition

### GPU not detected

**Cause:** GPU drivers or runtime not configured.

**Fix:**
1. Verify GPU is visible:
   - NVIDIA: `nvidia-smi`
   - AMD: `rocm-smi`
2. Check compiler supports GPU offloading
3. Ensure correct modules are loaded
4. Set `OMP_TARGET_OFFLOAD=MANDATORY` to force GPU (fails fast if unavailable)

### MPI errors or hangs

**Cause:** MPI configuration issues or deadlocks.

**Fix:**
1. Verify MPI works: `mpirun -n 2 hostname`
2. Try with fewer ranks: `./mfc.sh run case.py -n 1`
3. Check for incompatible MPI/compiler combinations
4. On clusters, ensure you're using the correct MPI for the interconnect

---

## Getting Help

If you can't resolve an issue:

1. **Search existing issues:** [GitHub Issues](https://github.com/MFlowCode/MFC/issues)

2. **Check documentation:** [MFC Documentation](https://mflowcode.github.io/)

3. **Open a new issue** with:
   - Your OS and compiler versions
   - The exact command you ran
   - Complete error output (use `-vv` for builds)
   - Your case file (if applicable)
   - Output of `./mfc.sh validate case.py`

4. **Use the CLI help:**
   ```bash
   ./mfc.sh help debugging
   ./mfc.sh help gpu
   ./mfc.sh <command> -h
   ```
