# MFC Spack Package

This document describes the Spack package for installing MFC on HPC systems and supercomputers.

## Overview

The Spack package enables installation of MFC on high-performance computing systems including leadership-class supercomputers. It handles complex dependencies, supports multiple compilers, enables GPU acceleration, and integrates with HPC module systems.

## Continuous Integration

The Spack package is automatically tested via GitHub Actions (`.github/workflows/spack.yml`) on every commit that modifies `packaging/spack/`. The CI workflow includes:

### Automated Tests
1. **Lint & Audit**: Validates package follows Spack style guidelines and conventions
2. **Concretization Tests**: Verifies dependency resolution for all variant combinations
3. **Installation Tests**: Performs actual builds with different configurations
4. **Conflict Detection**: Ensures invalid configurations are properly rejected

This ensures the package remains functional and compliant with Spack standards.

## What Gets Installed

When users run `spack install mfc`, they get:

### Binaries
- `pre_process` - Preprocessing binary for setting up initial conditions
- `simulation` - Main simulation binary for computational fluid dynamics
- `post_process` - Post-processing binary for analyzing results (optional)

### Integration
- Automatic module file generation for environment setup
- Integration with site-specific compiler and MPI configurations
- GPU support configuration when available
- Dependency management for all required libraries

## Package Structure

### Core Information

Package class: CMakePackage (uses CMake build system)
Homepage: https://mflowcode.github.io/
Source: GitHub releases and git repository
License: MIT
Maintainer: sbryngelson

### Versions

Available versions:
- master: Latest development version from master branch
- 5.1.0: Stable release with verified checksum

The package uses semantic versioning and can install any tagged release.

### Build Variants

Users can customize the installation with build variants:

#### MPI Support (default: enabled)
```
spack install mfc+mpi    # With MPI (default)
spack install mfc~mpi    # Without MPI
```
Enables parallel execution across multiple nodes.

#### OpenACC GPU Support (default: disabled)
```
spack install mfc+openacc
```
Enables GPU acceleration using OpenACC directives. Requires NVHPC or Cray compilers.

#### OpenMP GPU Support (default: disabled)
```
spack install mfc+openmp
```
Enables GPU acceleration using OpenMP target offloading.

#### Precision (default: double)
```
spack install mfc precision=single
spack install mfc precision=double
```
Selects floating-point precision for computations.

#### Post-processing (default: enabled)
```
spack install mfc+post_process    # Build post_process binary (default)
spack install mfc~post_process    # Skip post_process binary
```
Controls whether the post-processing tool is built.

#### Chemistry (default: disabled)
```
spack install mfc+chemistry
```
Enables thermochemistry support by generating the `m_thermochem.f90` module using Pyrometheus and Cantera during the build. When enabled, MFC can perform reactive flow simulations with detailed chemical kinetics. Requires `cantera+python` and vendors Pyrometheus automatically.

### Dependencies

Build-time dependencies (required during compilation):
- cmake@3.20 or newer - Build system
- py-fypp - Fortran preprocessor
- python@3.0 or newer - Build scripts

Runtime dependencies (always required):
- fftw@3.0 or newer - Fast Fourier Transform library
- lapack - Linear algebra routines

Optional dependencies (variant-controlled):
- mpi - Message Passing Interface (when +mpi)
- silo - Silo data format with HDF5 support (when +post_process)
- hdf5 - HDF5 data format (transitive dependency via Silo when +post_process)
- cantera+python - Thermochemical kinetics library (when +chemistry)
- pyrometheus - Fortran code generator for chemistry (vendored automatically when +chemistry)
- cuda - NVIDIA CUDA toolkit (when +openacc or +openmp with NVHPC)
- hip - AMD ROCm HIP (when +openacc or +openmp with Cray)

### Compiler Support

Supported compilers:
- GCC 5.0 or newer
- NVHPC 21.7 or newer (for OpenACC)
- Cray Compiler Environment (for OpenACC/OpenMP)
- Intel compilers

Explicitly unsupported:
- Apple Clang (conflicts with MFC requirements)
- GCC versions older than 5.0
- NVHPC versions older than 21.7

### Compiler Conflicts

The package defines explicit conflicts to prevent invalid configurations:
- Apple Clang is not supported due to Fortran requirements
- OpenACC requires NVHPC or Cray compilers (not available with GCC)
- Older compiler versions are rejected with helpful error messages

## Installation Process

### Build Configuration

The package generates CMake arguments based on selected variants:
- MFC_MPI: Enable/disable MPI support
- MFC_OpenACC: Enable/disable OpenACC GPU support
- MFC_OpenMP: Enable/disable OpenMP GPU support
- MFC_PRE_PROCESS: Always enabled
- MFC_SIMULATION: Always enabled
- MFC_POST_PROCESS: Controlled by post_process variant
- MFC_SINGLE_PRECISION: Set when precision=single

### Build Environment

The package configures the build environment:
- Adds py-fypp to PATH for Fortran preprocessing
- Sets up compiler wrappers for MPI when enabled
- Configures GPU toolchains when GPU variants are selected

### CMake Integration

Spack automatically handles CMake configuration:
1. Generates build directory
2. Runs CMake with package-specific arguments
3. Invokes parallel build
4. Installs binaries to Spack prefix
5. Generates module file for environment setup

## Usage Examples

### Basic Installation
```
spack install mfc
```
Installs MFC with default settings (MPI enabled, double precision, post-processing enabled).

### Custom Installation
```
spack install mfc+openacc %nvhpc
```
Installs with NVIDIA GPU support using NVHPC compiler.

### Development Version
```
spack install mfc@master
```
Installs the latest development version from the master branch.

### Full Custom Build
```
spack install mfc+mpi+openmp precision=double %gcc@12.0.0 ^openmpi@4.1.0
```
Specifies compiler version and MPI implementation explicitly.

### Load Module
```
spack load mfc
```
Loads the MFC environment including all dependencies.

## Platform Support

The package is designed for HPC systems:

### Tested Platforms
- OLCF Frontier (AMD GPUs with Cray compilers)
- OLCF Summit (NVIDIA GPUs with NVHPC)
- NERSC Perlmutter (NVIDIA GPUs)
- Generic Linux clusters

### Architecture Support
- x86_64 processors (Intel, AMD)
- ARM64 processors (on supported systems)
- NVIDIA GPUs (V100, A100, H100)
- AMD GPUs (MI100, MI250, MI300)

### Operating Systems
- Linux distributions (RHEL, CentOS, Ubuntu, SLES)
- HPC-specific OS variants

## Integration with HPC Systems

### Module System

After installation, Spack generates a module file that sets up:
- PATH to include MFC binaries
- Library paths for dependencies
- Environment variables for MPI and GPU support

Users load the module with:
```
module load mfc
```
or
```
spack load mfc
```

### Compiler Wrappers

On HPC systems, Spack uses system-provided compiler wrappers:
- cc/CC/ftn on Cray systems
- mpicc/mpicxx/mpif90 on clusters
- Automatically includes system libraries and paths

### GPU Configuration

When GPU variants are enabled:
- Selects appropriate GPU architecture flags
- Links against CUDA or HIP libraries
- Configures offloading runtime
- Sets up GPU-aware MPI if available

## Distribution

The package can be included in Spack in two ways:

### Official Spack Repository

Submit to the official Spack repository at https://github.com/spack/spack:
1. Fork the Spack repository
2. Add package.py to var/spack/repos/builtin/packages/mfc/
3. Run Spack style and audit checks
4. Submit pull request
5. Address review feedback
6. Package becomes available in official Spack

### Custom Spack Repository

Create a custom repository for immediate availability:
```
spack repo create my-repo
mkdir -p my-repo/packages/mfc
cp package.py my-repo/packages/mfc/package.py
spack repo add my-repo
spack install mfc
```
This allows distribution before official Spack acceptance.

## Advanced Features

### Concretization

Spack automatically determines the best configuration:
- Selects compatible compiler versions
- Resolves dependency conflicts
- Optimizes for the target architecture
- Reuses already-installed dependencies

### Dependency DAG

Users can view the dependency tree:
```
spack spec mfc
```
Shows all dependencies and their versions before installation.

### Binary Caches

Spack can use pre-compiled binaries:
```
spack mirror add facility https://mirror.url
spack buildcache keys --install --trust
spack install mfc
```
Downloads pre-built packages when available, falling back to source builds.

## Testing and Validation

### Package Audits
```
spack audit packages mfc
spack style --fix package.py
```
Validates package follows Spack conventions.

### Installation Testing
```
spack install --test=root mfc
```
Runs package tests after installation.

### Dependency Testing
```
spack install --test=all mfc
```
Tests MFC and all dependencies.

## Updates and Maintenance

### Adding New Versions

To add a new release:
1. Download tarball from GitHub
2. Calculate SHA256 checksum
3. Add version line to package.py
4. Test installation
5. Submit update

Example:
```python
version("5.2.0", sha256="<new_checksum>")
```

### Updating Dependencies

When MFC requirements change:
1. Update depends_on statements
2. Adjust version constraints if needed
3. Add new variants if needed
4. Test on representative systems
5. Submit update

### Deprecating Old Versions

Old versions can be deprecated:
```python
version("4.0.0", sha256="...", deprecated=True)
```

## Troubleshooting

### Build Failures

Check build logs:
```
spack install --verbose mfc
spack cd -b mfc  # Change to build directory
```

### Dependency Conflicts

View concrete spec before installation:
```
spack spec mfc
```

### Compiler Issues

Specify compiler explicitly:
```
spack install mfc %gcc@12.0.0
```

### GPU Problems

Check GPU variant configuration:
```
spack spec mfc+openacc
```

## Advantages Over Manual Installation

HPC users benefit from Spack installation:

1. Dependency Management: Automatically builds all required libraries
2. Reproducibility: Exact package versions recorded and reproducible
3. Module Integration: Seamless integration with existing module systems
4. Multi-Version Support: Multiple MFC versions can coexist
5. Compiler Management: Works with site-specific compiler installations
6. GPU Support: Automatic configuration of GPU acceleration
7. MPI Integration: Uses site-specific MPI implementations
8. Clean Uninstallation: Complete removal with dependencies

## Performance Considerations

### Optimization Flags

Spack applies appropriate optimization flags:
- Target-specific CPU optimizations
- GPU architecture flags when applicable
- Compiler-specific performance tuning

### MPI Configuration

The package respects site MPI configuration:
- Uses pre-installed MPI when available
- Builds MPI from source if needed
- Configures GPU-aware MPI on supported systems

### BLAS/LAPACK

Allows selection of optimized linear algebra libraries:
```
spack install mfc ^openblas
spack install mfc ^intel-mkl
```

## Security and Checksums

All release versions include SHA256 checksums to verify:
- Downloaded source matches expected file
- No corruption during download
- Authenticity of release

The master branch version does not have a checksum as it tracks the latest code.

## Contributing

To contribute package improvements:
1. Fork Spack repository
2. Modify package.py
3. Test changes on relevant systems
4. Run Spack style checks
5. Submit pull request to Spack
6. Work with reviewers to address feedback

## Documentation References

- Spack Documentation: https://spack.readthedocs.io/
- MFC Documentation: https://mflowcode.github.io/
- Package Development: https://spack.readthedocs.io/en/latest/packaging_guide.html

## Validation Status

The package has been tested with:
- Multiple compiler versions (GCC, NVHPC, Cray)
- Various MPI implementations (OpenMPI, MPICH, Cray MPICH)
- GPU configurations (CUDA, HIP)
- Different precision settings
- All variant combinations
- Installation on multiple HPC systems

All configurations install successfully and pass basic functionality tests.

