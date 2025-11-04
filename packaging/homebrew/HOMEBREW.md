# MFC Homebrew Formula

This document describes the Homebrew formula for installing MFC on macOS systems.

## Overview

The Homebrew formula enables one-command installation of MFC on macOS (both Intel and Apple Silicon). The formula handles all dependencies, builds the three main MFC binaries, and installs them in a standard Homebrew location.

## What Gets Installed

When users run `brew install mfc`, they get:

### Binaries
- `pre_process` - Preprocessing binary for setting up initial conditions
- `simulation` - Main simulation binary for computational fluid dynamics
- `post_process` - Post-processing binary for analyzing results
- `mfc` - Wrapper script that provides the full MFC interface

### Additional Components
- Python toolchain in `/usr/local/Cellar/mfc/VERSION/toolchain/`
- Example cases in `/usr/local/Cellar/mfc/VERSION/share/mfc/examples/`
- Documentation references and usage information

## Formula Structure

### Dependencies

Build-time dependencies (only needed during installation):
- cmake - Build system generator
- gcc - GNU Compiler Collection (provides gfortran)
- python@3.12 - Python 3.12 for build scripts

Runtime dependencies (needed to run MFC):
- boost - C++ libraries
- fftw - Fast Fourier Transform library
- hdf5 - Hierarchical Data Format 5 for data storage
- open-mpi - Message Passing Interface for parallel computing
- openblas - Optimized BLAS library

### Build Process

The formula executes the following steps during installation:

1. Sets up environment variables for compilers and libraries
2. Runs `./mfc.sh build` to compile all three binaries
3. Installs binaries to Homebrew's bin directory
4. Creates a wrapper script that sets up the environment
5. Installs Python toolchain and examples

### Environment Configuration

The formula ensures proper environment setup:
- BOOST_INCLUDE points to Homebrew's boost installation
- FC (Fortran compiler) is set to gfortran
- CC (C compiler) is set to gcc
- CXX (C++ compiler) is set to g++

### Wrapper Script

The installed `mfc` wrapper provides the complete MFC interface. It automatically configures the environment and delegates to the main `mfc.sh` script. Users can run any MFC command through this wrapper:

```
mfc build
mfc run examples/case.py
mfc test
mfc clean
```

## Installation Methods

### Standard Installation

Users install MFC with:
```
brew install mfc
```

This fetches the source tarball from GitHub releases, verifies the checksum, and builds from source.

### Development Installation

For the latest development version:
```
brew install --HEAD mfc
```

This clones from the master branch instead of using a release tarball.

## Testing

The formula includes automated tests that verify:
- All three binary files exist after installation
- The mfc wrapper script is functional
- The help command executes without errors

These tests run automatically during `brew install` and can be run manually with `brew test mfc`.

## Post-Installation

After installation completes, Homebrew displays usage information including:
- Locations of installed binaries
- Path to example cases
- Example command to run a simulation
- Link to full documentation

## Usage Examples

Once installed, users can immediately start using MFC:

```bash
# Run a test case
mfc run /usr/local/share/mfc/examples/1D_sodshocktube/case.py

# Run just preprocessing
pre_process -i input.dat

# Run simulation
simulation -i input.dat

# Post-process results
post_process -i input.dat
```

## Distribution

The formula can be distributed in two ways:

### Official Homebrew Repository
Submit a pull request to homebrew-core for inclusion in the main Homebrew repository. This requires:
- Stable release with version tag
- Verified tarball checksum
- Formula code review
- Automated testing passes

### Third-Party Tap
Create a separate tap (custom repository) for immediate availability:
```
brew tap organization/mfc
brew install organization/mfc/mfc
```

This allows distribution before official Homebrew acceptance.

## Platform Support

The formula supports:
- macOS 11 (Big Sur) and later
- Intel x86_64 processors
- Apple Silicon (ARM64) processors

Homebrew automatically selects the appropriate compiler flags and optimization settings for each architecture.

## Updates and Versioning

The formula specifies version 5.1.0 as the current release. To update:

1. Change the version number in the URL
2. Download the new tarball
3. Calculate new SHA256 checksum
4. Update the sha256 line in the formula
5. Test the installation
6. Submit updated formula

## Technical Details

### Source URL
The formula downloads from: https://github.com/MFlowCode/MFC/archive/refs/tags/v5.1.0.tar.gz

### Checksum Verification
SHA256: 4684bee6a529287f243f8929fb7edb0dfebbb04df7c1806459761c9a6c9261cf

This ensures the downloaded source matches the expected file exactly.

### Build Parallelization
The formula uses all available CPU cores for building (`ENV.make_jobs`) to minimize compilation time.

### Installation Prefix
Files install to the standard Homebrew prefix:
- Binaries: `/usr/local/bin/` (Intel) or `/opt/homebrew/bin/` (Apple Silicon)
- Data: `/usr/local/share/mfc/` or `/opt/homebrew/share/mfc/`
- Toolchain: `/usr/local/Cellar/mfc/VERSION/` or `/opt/homebrew/Cellar/mfc/VERSION/`

## Advantages Over Manual Installation

Users benefit from Homebrew installation:

1. Automatic dependency management - Homebrew installs all required libraries
2. Pre-compiled binaries - On some systems, bottles (binary packages) may be available
3. Easy updates - `brew upgrade mfc` gets the latest version
4. Clean uninstallation - `brew uninstall mfc` removes everything
5. Standard paths - Binaries are automatically in PATH
6. Version management - Multiple versions can coexist if needed

## Maintenance

The formula requires minimal maintenance once accepted:
- Update version and checksum when new releases are published
- Adjust dependencies if MFC requirements change
- Update minimum macOS version if newer features are needed
- Monitor for deprecated Homebrew APIs and update accordingly

## Validation

The formula has been tested on:
- macOS 14 (Sonoma) with Apple Silicon
- macOS 13 (Ventura) with Intel
- Fresh installations and upgrades
- Both release (5.1.0) and HEAD versions

All tests pass and the installation completes successfully on all tested platforms.

