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
- Python toolchain in `/usr/local/Cellar/mfc/VERSION/toolchain/` or `/opt/homebrew/Cellar/mfc/VERSION/toolchain/`
- Python virtual environment in `/usr/local/Cellar/mfc/VERSION/libexec/venv/` or `/opt/homebrew/Cellar/mfc/VERSION/libexec/venv/`
  - Pre-installed with Cantera 3.1.0 and MFC toolchain packages
- Example cases in `/usr/local/Cellar/mfc/VERSION/examples/` or `/opt/homebrew/Cellar/mfc/VERSION/examples/`
- Main `mfc.sh` script in `libexec/`
- Documentation references and usage information

## Formula Structure

### Dependencies

Build-time dependencies (only needed during installation):
- cmake - Build system generator
- gcc - GNU Compiler Collection (provides gfortran)

Runtime dependencies (needed to run MFC):
- boost - C++ libraries
- fftw - Fast Fourier Transform library
- hdf5 - Hierarchical Data Format 5 for data storage
- open-mpi - Message Passing Interface for parallel computing
- openblas - Optimized BLAS library
- python@3.12 - Python 3.12 (used for build scripts and runtime virtual environment)

Python package dependencies (installed in virtual environment):
- cantera==3.1.0 - Chemical kinetics library (required for MFC build and runtime)
- MFC toolchain package - Python utilities for MFC operations (installed in editable mode)

### Virtual Environment Setup

Before building, the formula creates a Python virtual environment to isolate MFC's Python dependencies:

1. Creates a new virtual environment in `libexec/venv` using Python 3.12
2. Upgrades pip, setuptools, and wheel to latest versions
3. Installs Cantera 3.1.0 from PyPI (required dependency for MFC)
4. Installs the MFC toolchain package in editable mode (`-e`) to avoid RECORD file issues
5. Creates a symlink from `build/venv` to the venv so `mfc.sh` can find it during build

This virtual environment persists after installation and is used by the wrapper script at runtime.

### Build Process

The formula executes the following steps during installation:

1. **Virtual Environment Setup** (see section above)
2. Sets `VIRTUAL_ENV` environment variable so `mfc.sh` uses the pre-configured venv
3. Runs `./mfc.sh build -t pre_process simulation post_process -j <cores>` to compile all three binaries
4. Installs binaries from `build/install/*/bin/*` to Homebrew's bin directory
5. Installs `mfc.sh` to `libexec/` for script execution
6. Installs Python toolchain directory to `prefix/toolchain/` (required for mfc.sh functionality)
7. Installs examples directory to `prefix/examples/`
8. Creates and installs the `mfc` wrapper script that handles runtime environment setup

### Environment Configuration

The formula relies on Homebrew's automatic environment setup:
- Compiler flags are set via the `gcc` dependency (Homebrew's superenv)
- Library paths are automatically configured
- The wrapper script sets BOOST_INCLUDE at runtime for user commands

### Wrapper Script

The installed `mfc` wrapper provides the complete MFC interface. Due to Homebrew's read-only Cellar structure, the wrapper implements a sophisticated workaround:

**Key Features:**
1. **Temporary Working Directory**: Creates a temporary directory since the Cellar is read-only and `mfc.sh` may need to write build artifacts
2. **Environment Setup**: Copies `mfc.sh`, toolchain, examples, and the virtual environment to the temp directory
3. **Toolchain Patching**: Dynamically patches the toolchain Python code to:
   - Use pre-installed binaries from Homebrew's bin directory instead of building new ones
   - Skip building main targets (`pre_process`, `simulation`, `post_process`, `syscheck`) since they're already installed
4. **Build Optimization**: Automatically adds `--no-build` flag to `mfc run` commands to skip unnecessary compilation
5. **Cleanup**: Automatically removes the temporary directory when the command completes

Users can run any MFC command through this wrapper:

```bash
mfc build
mfc run examples/case.py
mfc test
mfc clean
```

The wrapper ensures all MFC functionality works correctly while respecting Homebrew's installation constraints.

## Installation Methods

### Standard Installation

Users install MFC with:
```bash
brew install mfc
```

This fetches the source tarball from GitHub releases, verifies the checksum, and builds from source.

### Development Installation

For the latest development version:
```bash
brew install --HEAD mfc
```

This clones from the master branch instead of using a release tarball.

## Testing

The formula includes automated tests that verify:
- All three binary files exist after installation and are executable
- The Python toolchain directory is installed correctly
- The virtual environment exists and contains an executable Python interpreter
- The examples directory is installed
- The mfc wrapper script is functional and responds to `--help`

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
# Copy an example case to your working directory
cp $(brew --prefix mfc)/examples/1D_sodshocktube/case.py .

# Run a test case
mfc run case.py

# Or run directly from the installed examples directory
mfc run $(brew --prefix mfc)/examples/1D_sodshocktube/case.py

# Run just preprocessing
pre_process -i input.dat

# Run simulation
simulation -i input.dat

# Post-process results
post_process -i input.dat
```

Note: The `brew --prefix mfc` command returns the installation prefix (e.g., `/usr/local/Cellar/mfc/VERSION` or `/opt/homebrew/Cellar/mfc/VERSION`), making examples work on both Intel and Apple Silicon systems.

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
```bash
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

When a new MFC version is released, the formula must be updated. To update to a new version:

1. Change the version number in the URL
2. Download the new tarball
3. Calculate new SHA256 checksum
4. Update the sha256 line in the formula
5. Test the installation
6. Submit updated formula

## Technical Details

### Source URL
The formula downloads from GitHub releases using the version tag specified in the formula (e.g., v5.1.0):
```
https://github.com/MFlowCode/MFC/archive/refs/tags/vVERSION.tar.gz
```

### Checksum Verification
The formula includes a SHA256 checksum that must match the downloaded tarball. This ensures:
- The downloaded source matches the expected file exactly
- No corruption occurred during download
- Authenticity of the release

Maintainers must update this checksum when releasing new versions.

### Build Parallelization
The formula uses all available CPU cores for building (`ENV.make_jobs`) to minimize compilation time.

### Installation Prefix
Files install to the standard Homebrew prefix:
- Binaries: `/usr/local/bin/` (Intel) or `/opt/homebrew/bin/` (Apple Silicon)
- Main script: `/usr/local/Cellar/mfc/VERSION/libexec/mfc.sh` or `/opt/homebrew/Cellar/mfc/VERSION/libexec/mfc.sh`
- Toolchain: `/usr/local/Cellar/mfc/VERSION/toolchain/` or `/opt/homebrew/Cellar/mfc/VERSION/toolchain/`
- Virtual environment: `/usr/local/Cellar/mfc/VERSION/libexec/venv/` or `/opt/homebrew/Cellar/mfc/VERSION/libexec/venv/`
- Examples: `/usr/local/Cellar/mfc/VERSION/examples/` or `/opt/homebrew/Cellar/mfc/VERSION/examples/`

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
- Both stable release versions and HEAD (development) versions

All tests pass and the installation completes successfully on all tested platforms.

