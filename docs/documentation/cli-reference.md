@page cli-reference Command Line Reference

# Command Line Reference

> **Auto-generated** from `toolchain/mfc/cli/commands.py`
> 
> Regenerate with: `./mfc.sh generate`

## Overview

Welcome to the MFC master script. This tool automates and manages building, testing, running, and cleaning of MFC in various configurations on all supported platforms. The README documents this tool and its various commands in more detail. To get started, run `./mfc.sh build -h`.

## Quick Reference

| Command | Alias | Description |
|---------|-------|-------------|
| [<code>build</code>](#build) | `b` | Build MFC and its dependencies. |
| [<code>run</code>](#run) | `r` | Run a case with MFC. |
| [<code>test</code>](#test) | `t` | Run MFC's test suite. |
| [<code>clean</code>](#clean) | `c` | Clean build artifacts. |
| [<code>validate</code>](#validate) | `v` | Validate a case file without running. |
| [<code>new</code>](#new) | - | Create a new case from a template. |
| [<code>params</code>](#params) | - | Search and explore MFC case parameters. |
| [<code>packer</code>](#packer) | - | Packer utility (pack/unpack/compare). |
| [<code>completion</code>](#completion) | - | Install shell tab-completion. |
| [<code>help</code>](#help) | - | Show help on a topic. |
| [<code>generate</code>](#generate) | - | Regenerate completion scripts from CLI schema. |
| [<code>load</code>](#load) | - | Loads the MFC environment with source. |
| [<code>lint</code>](#lint) | - | Lints all code after editing. |
| [<code>format</code>](#format) | - | Formats all code after editing. |
| [<code>spelling</code>](#spelling) | - | Runs the spell checker after editing. |
| [<code>interactive</code>](#interactive) | - | Launch interactive menu-driven interface. |
| [<code>bench</code>](#bench) | - | Benchmark MFC (for CI). |
| [<code>bench_diff</code>](#bench_diff) | - | Compare MFC Benchmarks (for CI). |
| [<code>count</code>](#count) | - | Count LOC in MFC. |
| [<code>count_diff</code>](#count_diff) | - | Compare LOC between branches. |

## Commands

### build (alias: `b`)

Build MFC targets with optional GPU support and case optimization.

**Usage:** `./mfc.sh build [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `-i`, `--input` | (GPU Optimization) Build a version of MFC optimized for a case. | - |
| `--case-optimization` | (GPU Optimization) Compile MFC targets with some case parameters hard-coded (requires --input). | `false` |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Build all default targets (CPU)
./mfc.sh build

# Build with 8 parallel jobs
./mfc.sh build -j 8

# Build with GPU (OpenACC) support
./mfc.sh build --gpu

# Build only simulation target
./mfc.sh build -t simulation

```

---

### run (alias: `r`)

Run an MFC simulation case interactively or submit as a batch job.

**Usage:** `./mfc.sh run INPUT [OPTIONS]`

**Arguments:**

- `INPUT` - Input file to run.

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `-g`, `--gpus` | (Optional GPU override) List of GPU #s to use (environment default if unspecified). | - |
| `-e`, `--engine` | Job execution/submission engine choice. | `interactive` |
| `-p`, `--partition` | (Batch) Partition for job submission. | - |
| `-q`, `--quality_of_service` | (Batch) Quality of Service for job submission. | - |
| `-N`, `--nodes` | (Batch) Number of nodes. | `1` |
| `-n`, `--tasks-per-node` | Number of tasks per node. | `1` |
| `-w`, `--walltime` | (Batch) Walltime. | `01:00:00` |
| `-a`, `--account` | (Batch) Account to charge. | - |
| `-@`, `--email` | (Batch) Email for job notification. | - |
| `-#`, `--name` | (Batch) Job name. | `MFC` |
| `-s`, `--scratch` | Build from scratch. | `false` |
| `-b`, `--binary` | (Interactive) Override MPI execution binary | - |
| `--dry-run` | (Batch) Run without submitting batch file. | `false` |
| `--case-optimization` | (GPU Optimization) Compile MFC targets with some case parameters hard-coded. | `false` |
| `--no-build` | (Testing) Do not rebuild MFC. | `false` |
| `--wait` | (Batch) Wait for the job to finish. | `false` |
| `-c`, `--computer` | (Batch) Path to a custom submission file template or one of the built-in templates. | `default` |
| `-o`, `--output-summary` | Output file (YAML) for summary. | - |
| `--clean` | Clean the case before running. | `false` |
| `--ncu` | Profile with NVIDIA Nsight Compute. | - |
| `--nsys` | Profile with NVIDIA Nsight Systems. | - |
| `--rcu` | Profile with ROCM rocprof-compute. | - |
| `--rsys` | Profile with ROCM rocprof-systems. | - |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Run interactively with 1 rank
./mfc.sh run case.py

# Run with 4 MPI ranks
./mfc.sh run case.py -n 4

# Submit batch job: 2 nodes, 4 ranks/node
./mfc.sh run case.py -e batch -N 2 -n 4

# Submit with account
./mfc.sh run case.py -e batch -a myaccount

```

---

### test (alias: `t`)

Run MFC's test suite with various filtering and generation options.

**Usage:** `./mfc.sh test [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `-g`, `--gpus` | (Optional GPU override) List of GPU #s to use (environment default if unspecified). | - |
| `-l`, `--list` | List all available tests. | - |
| `-f`, `--from` | First test UUID to run. | - |
| `-t`, `--to` | Last test UUID to run. | - |
| `-o`, `--only` | Only run tests with specified properties. | `[]` |
| `-a`, `--test-all` | Run the Post Process Tests too. | `false` |
| `-%`, `--percent` | Percentage of tests to run. | `100` |
| `-m`, `--max-attempts` | Maximum number of attempts to run a test. | `1` |
| `--rdma-mpi` | Run tests with RDMA MPI enabled | `false` |
| `--no-build` | (Testing) Do not rebuild MFC. | `false` |
| `--no-examples` | Do not test example cases. | `false` |
| `--case-optimization` | (GPU Optimization) Compile MFC targets with some case parameters hard-coded. | `false` |
| `--dry-run` | Build and generate case files but do not run tests. | `false` |
| `--generate` | (Test Generation) Generate golden files. | `false` |
| `--add-new-variables` | (Test Generation) If new variables are found in D/ when running tests, add them to the golden files. | `false` |
| `--remove-old-tests` | (Test Generation) Delete tests directories that are no longer. | `false` |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Run all tests
./mfc.sh test

# Run with 4 parallel jobs
./mfc.sh test -j 4

# Run only 3D tests
./mfc.sh test --only 3D

# Regenerate golden files
./mfc.sh test --generate

```

---

### clean (alias: `c`)

Remove build artifacts and cache files.

**Usage:** `./mfc.sh clean [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Clean all build files
./mfc.sh clean

```

---

### validate (alias: `v`)

Check a case file for syntax errors and constraint violations.

**Usage:** `./mfc.sh validate INPUT [OPTIONS]`

**Arguments:**

- `INPUT` - Path to case file to validate.

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |

**Examples:**

```bash
# Check syntax and constraints
./mfc.sh validate case.py

# Validate with debug output
./mfc.sh validate case.py -d

```

---

## Utility Commands

### new

Create a new simulation case directory from a built-in or example template.

**Usage:** `./mfc.sh new [NAME] [OPTIONS]`

**Arguments:**

- `NAME` - Name/path for the new case directory.

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--template` | Template to use (e.g., 1D_minimal, 2D_minimal, 3D_minimal, or example:<name>). | `1D_minimal` |
| `-l`, `--list` | List available templates. | - |

**Examples:**

```bash
# Create with 1D_minimal template
./mfc.sh new my_case

# Create with 2D template
./mfc.sh new my_case -t 2D_minimal

# Copy from example
./mfc.sh new my_case -t example:3D_sphbubcollapse

# List available templates
./mfc.sh new --list

```

---

### params

Search, list, and get information about MFC's ~3,300 case parameters.

**Usage:** `./mfc.sh params [QUERY] [OPTIONS]`

**Arguments:**

- `QUERY` - Search query (parameter name or pattern to search for).

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--type` | Filter by type: int, real, log, str. | - |
| `-f`, `--families` | List parameter families (grouped by prefix). | `false` |
| `-c`, `--count` | Show count statistics only. | `false` |
| `-n`, `--limit` | Maximum number of results to show. | `25` |
| `-d`, `--describe` | Show parameter descriptions. | `false` |

**Examples:**

```bash
# Search for parameters containing 'weno'
./mfc.sh params weno

# Search for bubble-related parameters
./mfc.sh params bubble

# Show patch_icpp parameters
./mfc.sh params patch_icpp

# List all parameter families
./mfc.sh params -f

# Show parameter count statistics
./mfc.sh params -c

# Show only REAL type weno params
./mfc.sh params -t real weno

```

---

### packer

Pack simulation output into a single file or compare packed files.

**Usage:** `./mfc.sh packer [OPTIONS]`

**Examples:**

```bash
# Pack case output
./mfc.sh packer pack case.py

# Compare two packed files
./mfc.sh packer compare a.pack b.pack

```

**Subcommands:**

#### packer pack

Pack a case into a single file.

Arguments:
- `INPUT` - Input file of case to pack.

Options:

| Option | Description | Default |
|--------|-------------|---------|
| `-o`, `--output` | Base name of output file. | - |

#### packer compare

Compare two cases.

Arguments:
- `INPUT1` - First pack file.
- `INPUT2` - Second pack file.

Options:

| Option | Description | Default |
|--------|-------------|---------|
| `-rel`, `--reltol` | Relative tolerance. | `1e-12` |
| `-abs`, `--abstol` | Absolute tolerance. | `1e-12` |

---

### completion

Install or manage shell tab-completion for MFC commands.

**Usage:** `./mfc.sh completion [COMPLETION_ACTION] [COMPLETION_SHELL] [OPTIONS]`

**Arguments:**

- `COMPLETION_ACTION` - Action: install, uninstall, or status
- `COMPLETION_SHELL` - Shell type: bash or zsh (auto-detected if not specified)

**Examples:**

```bash
# Install tab completion for current shell
./mfc.sh completion install

# Install for bash specifically
./mfc.sh completion install bash

# Check completion installation status
./mfc.sh completion status

```

---

### help

Show help on a topic.

**Usage:** `./mfc.sh help [TOPIC] [OPTIONS]`

**Arguments:**

- `TOPIC` - Help topic: gpu, clusters, batch, debugging

---

### generate

Regenerate shell completion scripts, documentation, and JSON schema from the CLI schema.

**Usage:** `./mfc.sh generate [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--check` | Check if generated files are up to date (exit 1 if not). | `false` |
| `--json-schema` | Generate JSON Schema for IDE auto-completion of case files. | `false` |

**Examples:**

```bash
# Regenerate completion scripts
./mfc.sh generate

# Check if completions are up to date
./mfc.sh generate --check

# Generate JSON Schema for IDE support
./mfc.sh generate --json-schema

```

---

## Development Commands

### lint

Run pylint on MFC's toolchain Python code.

**Usage:** `./mfc.sh lint [OPTIONS]`

**Examples:**

```bash
# Run pylint on Python code
./mfc.sh lint

```

---

### format

Format all code in the repository.

**Usage:** `./mfc.sh format [OPTIONS]`

**Examples:**

```bash
# Format all code
./mfc.sh format

# Format with 8 parallel jobs
./mfc.sh format -j 8

```

---

### spelling

Check spelling in the codebase.

**Usage:** `./mfc.sh spelling [OPTIONS]`

**Examples:**

```bash
# Run spell checker
./mfc.sh spelling

```

---

### count

Count LOC in MFC.

**Usage:** `./mfc.sh count [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Show LOC statistics
./mfc.sh count

```

---

### count_diff

Compare LOC between branches.

**Usage:** `./mfc.sh count_diff [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

---

## CI Commands

### bench

Benchmark MFC (for CI).

**Usage:** `./mfc.sh bench [OPTIONS]`

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-t`, `--targets` | Space separated list of targets to act upon. | `pre_process, simulation, post_process` |
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `-g`, `--gpus` | (Optional GPU override) List of GPU #s to use (environment default if unspecified). | - |
| `-o`, `--output` | Path to the YAML output file to write the results to. | - |
| `-m`, `--mem` | Memory per task for benchmarking cases | `1` |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

**Examples:**

```bash
# Run benchmarks and save results
./mfc.sh bench -o results.yaml

```

---

### bench_diff

Compare MFC Benchmarks (for CI).

**Usage:** `./mfc.sh bench_diff LHS RHS [OPTIONS]`

**Arguments:**

- `LHS` - Path to a benchmark result YAML file.
- `RHS` - Path to a benchmark result YAML file.

**Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `-j`, `--jobs` | Allows for JOBS concurrent jobs. | `1` |
| `-v`, `--verbose` | Increase verbosity level. Use -v, -vv, or -vvv for more detail. | `0` |
| `-d`, `--debug-log` | Enable debug logging for troubleshooting. | - |
| `--mpi`, `--no-mpi` | Enable/disable MPI | `true` |
| `--gpu [acc/mp]`, `--no-gpu` | Enable GPU (OpenACC/OpenMP) | `no` |
| `--debug`, `--no-debug` | Enable debug mode | `false` |

---

## Other Commands

### load

Load MFC environment modules (use with source).

**Usage:** `./mfc.sh load [OPTIONS]`

**Examples:**

```bash
# Load Phoenix GPU modules
source ./mfc.sh load -c p -m g

# Load Frontier CPU modules
source ./mfc.sh load -c f -m c

```

---

### interactive

Launch an interactive menu for MFC operations.

**Usage:** `./mfc.sh interactive [OPTIONS]`

---

## Common Options

Many commands share common option sets:

### Target Selection (`-t, --targets`)

Available targets:
- `pre_process` - Pre-processor
- `simulation` - Main simulation
- `post_process` - Post-processor
- `syscheck` - System check utility
- `documentation` - Build documentation

### Build Configuration Flags

| Flag | Description |
|------|-------------|
| `--mpi` / `--no-mpi` | Enable/disable MPI support |
| `--gpu [acc/mp]` / `--no-gpu` | Enable GPU with OpenACC or OpenMP |
| `--debug` / `--no-debug` | Enable debug build |
| `--gcov` / `--no-gcov` | Enable code coverage |
| `--single` / `--no-single` | Single precision |
| `--mixed` / `--no-mixed` | Mixed precision |

### Verbosity (`-v, --verbose`)

- `-v` - Basic verbose output
- `-vv` - Show build commands
- `-vvv` - Full debug output including CMake debug
