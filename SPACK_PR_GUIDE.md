# Spack Package Submission Guide for MFC

## Overview
This guide walks through submitting MFC to the Spack package manager, making it installable via `spack install mfc` on HPC systems worldwide.

## Steps

### 1. Fork and Clone Spack

```bash
# Fork https://github.com/spack/spack on GitHub first
git clone https://github.com/YOUR_USERNAME/spack.git
cd spack
git remote add upstream https://github.com/spack/spack.git
source share/spack/setup-env.sh
```

### 2. Create the Package

```bash
# Create the package directory
mkdir -p var/spack/repos/builtin/packages/mfc

# Copy the package.py file from MFC repo
cp /path/to/MFC/package.py var/spack/repos/builtin/packages/mfc/
```

### 3. Test the Package Locally

```bash
# Check package syntax
spack info mfc

# Check dependency resolution
spack spec mfc
spack spec mfc+openacc %nvhpc

# Try a build (optional but recommended)
spack install --test=root mfc
```

### 4. Run Spack's Style Checks

```bash
# Run flake8 linter
spack style --tool flake8 var/spack/repos/builtin/packages/mfc/package.py

# Run black formatter
spack style --tool black var/spack/repos/builtin/packages/mfc/package.py

# Fix any issues reported
```

### 5. Create a Branch and Commit

```bash
git checkout develop
git pull upstream develop
git checkout -b mfc-package

git add var/spack/repos/builtin/packages/mfc/package.py
git commit -m "mfc: new package for exascale multiphysics flow solver"
```

### 6. Push and Create PR

```bash
git push origin mfc-package
```

Then go to https://github.com/YOUR_USERNAME/spack and create a Pull Request to `spack/spack:develop`.

---

## PR Title

```
mfc: new package for exascale multiphysics flow solver
```

## PR Description Template

```markdown
### Description

Add MFC (Multicomponent Flow Code), an exascale multiphase/multiphysics compressible flow solver.

**Key features:**
- Scales to 43K+ GPUs on leadership-class supercomputers (OLCF Frontier, LLNL El Capitan)
- High-order WENO/TENO schemes for shock and interface capturing
- Immersed boundary method for complex geometries
- GPU acceleration via OpenACC and OpenMP offloading
- Modern Fortran 2008+ with Fypp metaprogramming
- 2025 Gordon Bell Prize Finalist

**Homepage:** https://mflowcode.github.io/  
**Repository:** https://github.com/MFlowCode/MFC  
**Paper:** https://doi.org/10.48550/arXiv.2503.07953

### Checklist

- [x] Package builds successfully with `spack install mfc`
- [x] Package passes `spack style` checks
- [x] Added maintainer (@sbryngelson)
- [x] Tested on macOS/Linux
- [ ] (Optional) Tested GPU variants

### Testing

```bash
# CPU-only build
spack install mfc

# With MPI
spack install mfc+mpi

# GPU build (requires NVHPC or Cray)
spack install mfc+openacc %nvhpc
```

### Notes

- MFC requires Fypp preprocessor (py-fypp) which is already in Spack
- GPU builds require NVHPC or Cray compilers (conflicts with GCC for GPU)
- Post-processing tools depend on HDF5 and SILO
```

---

## Expected Timeline

- **PR submission:** Immediate
- **Initial review:** 1-2 weeks
- **Revisions (if needed):** 1-2 rounds
- **Merge:** 2-4 weeks total

## After Merge

Once merged, MFC will be installable on any HPC system via:

```bash
spack install mfc
spack install mfc+mpi+openacc %nvhpc  # GPU version
```

This gives MFC:
- ✅ One-line install on thousands of HPC systems
- ✅ Automatic dependency management
- ✅ Integration with facility module systems
- ✅ Visibility in Spack package search

## Contact

Questions about the package? Tag @sbryngelson in the PR or reach out via:
- GitHub: https://github.com/MFlowCode/MFC/discussions
- Email: shb@gatech.edu







