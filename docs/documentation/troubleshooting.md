# Troubleshooting Guide

This guide covers common issues you may encounter when building, running, or testing MFC.

## Build Errors

### "CMake could not find MPI"

**Cause:** MPI is not installed or not in your PATH.

**Fix:**
- **Ubuntu/Debian:** `sudo apt install libopenmpi-dev openmpi-bin`
- **macOS (Homebrew):** `brew install open-mpi`
- **HPC systems:** `module load openmpi` (or similar)

After installing, verify with: `mpirun --version`

### "CMake could not find a Fortran compiler"

**Cause:** No Fortran compiler is installed or not in PATH.

**Fix:**
- **Ubuntu/Debian:** `sudo apt install gfortran`
- **macOS (Homebrew):** `brew install gcc`
- **HPC systems:** `module load gcc` or `module load nvhpc`

### "Fypp preprocessing failed"

**Cause:** Usually a syntax error in `.fpp` files or missing Fypp.

**Fix:**
1. Ensure Fypp is installed: `pip install fypp`
2. Check the specific error line mentioned in the output
3. If contributing code, run `./mfc.sh format` to check for issues

### Build fails with OpenACC/GPU errors

**Cause:** GPU compiler not properly configured.

**Fix:**
1. For NVIDIA GPUs, ensure NVHPC is installed and loaded: `module load nvhpc`
2. Set `CUDA_CC` environment variable if needed: `export MFC_CUDA_CC=80`
3. Try building without GPU first: `./mfc.sh build --no-gpu`

## Runtime Errors

### "Case parameter constraint violations"

**Cause:** Invalid combination of simulation parameters.

**Fix:**
1. Run `./mfc.sh validate case.py` for detailed diagnostics
2. Check the specific constraint mentioned in the error
3. Refer to `docs/documentation/case.md` for parameter documentation
4. Look at similar examples in the `examples/` directory

Common issues:
- `m`, `n`, `p` grid dimensions not matching the dimensionality
- `weno_order` not compatible with grid size
- Boundary conditions not matching the domain setup

### "Golden file mismatch" (during tests)

**Cause:** Numerical results differ from the expected reference values.

**Fix:**
1. If you intentionally changed the physics/numerics, regenerate golden files:
   ```bash
   ./mfc.sh test --generate --from <UUID> --to <UUID>
   ```
2. If unexpected, check for:
   - Compiler differences (different compilers may give slightly different results)
   - Precision settings (single vs. double)
   - Platform-specific numerical differences

### "NaN detected in the case"

**Cause:** Numerical instability in the simulation.

**Fix:**
1. Reduce the time step (`dt`)
2. Check initial conditions for discontinuities
3. Verify boundary conditions are appropriate
4. Consider using a more diffusive scheme initially

### Test timeout (exceeded 1 hour)

**Cause:** Test is hanging or taking too long.

**Fix:**
1. Check if the case is too large for the test environment
2. Verify MPI is working: `mpirun -n 2 hostname`
3. Check for deadlocks in parallel code
4. Review the case configuration for infinite loops

## Common Issues

### `./mfc.sh clean` doesn't work

**Cause:** The command was not properly registered.

**Fix:** This has been fixed in recent versions. If you have an older version:
```bash
rm -rf build/
```

### Tests pass locally but fail in CI

**Cause:** Environment differences between local and CI.

**Fix:**
1. Check compiler versions match
2. Verify precision settings (single/double/mixed)
3. Look for platform-specific code paths
4. Run with the same flags as CI: check `.github/workflows/test.yml`

### GPU not detected

**Cause:** GPU drivers or compiler not properly configured.

**Fix:**
1. Verify GPU is visible: `nvidia-smi` (NVIDIA) or `rocm-smi` (AMD)
2. Check compiler supports GPU: `nvfortran --version` or `ftn --version`
3. Ensure correct modules are loaded on HPC systems
4. Set `OMP_TARGET_OFFLOAD=MANDATORY` to force GPU usage

### "Module not found" on HPC

**Cause:** Required modules not loaded.

**Fix:**
1. Check available modules: `module avail`
2. Load required modules (example for typical setup):
   ```bash
   module load gcc openmpi cmake python
   ```
3. Use `./mfc.sh load` if available for your system

## Getting Help

If you can't resolve an issue:

1. Check existing [GitHub Issues](https://github.com/MFlowCode/MFC/issues)
2. Search the [documentation](https://mflowcode.github.io/)
3. Open a new issue with:
   - Your OS and compiler versions
   - The exact command you ran
   - The complete error output
   - Your case file (if applicable)
