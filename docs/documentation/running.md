@page running Running

# Running

MFC can be run using `mfc.sh`'s `run` command.
It supports interactive and batch execution.
Batch mode is designed for multi-node distributed systems (supercomputers) equipped with a scheduler such as PBS, SLURM, or LSF.
A full (and up-to-date) list of available arguments can be acquired with `./mfc.sh run -h`.

MFC supports running simulations locally (Linux, MacOS, and Windows) as well as
several supercomputer clusters, both interactively and through batch submission.

## Using the Homebrew package (macOS)

If you installed MFC via Homebrew, run cases with the `mfc` wrapper:

```bash
mfc <path/to/case.py> -n 2
```

- Use `-n X` to control the number of MPI processes (ranks).
- The Homebrew package uses a simplified syntax: just `mfc <case.py>` to run cases.
- To use developer commands (`build`, `test`, `clean`, etc.), clone the repository and use `./mfc.sh`.
- The wrapper passes through runtime flags like `-t pre_process simulation`, `-n`, and others; it always runs with preinstalled binaries.
- Examples live at `$(brew --prefix mfc)/examples/`.

> [!IMPORTANT]
> Running simulations locally should work out of the box. On supported clusters,
> you can append `-c <computer name>` on the command line to instruct the MFC toolchain
> to make use of the template file `toolchain/templates/<computer name>.mako`. You can
> browse that directory and contribute your own files. Since systems and their schedulers
> do not have a standardized syntax to request certain resources, MFC can only provide
> support for a restricted subset of common or user-contributed configuration options.
>
> Adding a new template file or modifying an existing one will most likely be required if:
> - You are on a cluster that does not have a template yet.
> - Your cluster is configured with SLURM, but interactive job launches fail when
>   using `srun`. You might need to invoke `mpirun` instead.
> - Something in the existing default or computer template file is incompatible with
>   your system or does not provide a feature you need.
>
> If `-c <computer name>` is left unspecified, it defaults to `-c default`.

Please refer to `./mfc.sh run -h` for a complete list of arguments and options, along with their defaults.

---

## Interactive Execution

To run all stages of MFC, that is [pre_process](https://github.com/MFlowCode/MFC/tree/master/src/pre_process/), [simulation](https://github.com/MFlowCode/MFC/tree/master/src/simulation/), and [post_process](https://github.com/MFlowCode/MFC/tree/master/src/post_process/) on the sample case [2D_shockbubble](https://github.com/MFlowCode/MFC/tree/master/examples/2D_shockbubble/),

```shell
./mfc.sh run examples/2D_shockbubble/case.py
```

If you want to run a subset of the available stages, you can use the `-t` argument.
To use multiple threads, use the `-n` option along with the number of threads you wish to use.
If a (re)build is required, it will be done automatically, with the number of threads
specified with the `-j` option.

For example,

- Running [pre_process](https://github.com/MFlowCode/MFC/tree/master/src/pre_process/) with 2 cores:

```shell
./mfc.sh run examples/2D_shockbubble/case.py -t pre_process -n 2
```

- Running [simulation](https://github.com/MFlowCode/MFC/tree/master/src/simulation/) and [post_process](https://github.com/MFlowCode/MFC/tree/master/src/post_process/)
using 4 cores:

```shell
./mfc.sh run examples/2D_shockbubble/case.py -t simulation post_process -n 4
```

---

## Running on GPUs

MFC supports GPU acceleration via OpenACC (default) or OpenMP offloading.
This section covers how to build and run MFC on GPU systems.

### Building with GPU Support

First, build MFC with GPU support enabled:

```shell
# Using OpenACC (default, recommended for NVIDIA)
./mfc.sh build --gpu

# Explicitly specify OpenACC
./mfc.sh build --gpu acc

# Using OpenMP offloading (alternative)
./mfc.sh build --gpu mp
```

On HPC clusters, load the appropriate modules first:

```shell
source ./mfc.sh load -c <cluster> -m g   # 'g' for GPU mode
./mfc.sh build --gpu -j $(nproc)
```

### Running on GPUs

Run simulations with GPU support:

```shell
# Basic GPU run
./mfc.sh run case.py --gpu

# Specify GPU IDs (useful for multi-GPU nodes)
./mfc.sh run case.py --gpu -g 0 1 2 3

# Run with 4 MPI ranks (typically one per GPU)
./mfc.sh run case.py -n 4 --gpu
```

### Supported Compilers

| Vendor | Compiler | OpenACC | OpenMP |
|--------|----------|---------|--------|
| NVIDIA | NVHPC (nvfortran) | Yes | Yes |
| AMD | Cray (ftn) | Yes | Yes |
| AMD | AMD (amdflang) | No | Yes |

### Environment Setup

Most HPC systems require loading GPU-specific modules:

**NVIDIA systems:**
```shell
module load nvhpc cuda
# Or use MFC's module loader:
source ./mfc.sh load -c phoenix -m g
```

**AMD systems:**
```shell
module load rocm amdflang
# Or use MFC's module loader:
source ./mfc.sh load -c frontier -m g
```

### Verifying GPU Detection

Check that GPUs are visible before running:

```shell
# NVIDIA
nvidia-smi

# AMD
rocm-smi
```

To force GPU usage (fails fast if no GPU available):
```shell
export OMP_TARGET_OFFLOAD=MANDATORY
./mfc.sh run case.py --gpu
```

### GPU Profiling

MFC integrates with vendor profiling tools for performance analysis.

#### NVIDIA GPUs

**Nsight Systems** (timeline/system-level profiling):
```shell
./mfc.sh run case.py -t simulation --nsys [nsys flags]
```
- Best for understanding execution order and timing of major subroutines (WENO, Riemann, etc.)
- Generates `.nsys-rep` files in the case directory
- Open with [NVIDIA Nsight Systems GUI](https://developer.nvidia.com/nsight-systems/get-started)
- Use few timesteps to keep report files small

**Nsight Compute** (kernel-level profiling):
```shell
./mfc.sh run case.py -t simulation --ncu [ncu flags]
```
- Detailed per-kernel metrics: cycles, memory usage, occupancy
- Significantly slower than normal execution
- Use only with few timesteps

#### AMD GPUs

**rocprof-systems** (timeline profiling):
```shell
./mfc.sh run case.py -t simulation --rsys --hip-trace [rocprof flags]
```
- Generates files viewable in [Perfetto UI](https://ui.perfetto.dev/)
- Use few timesteps for manageable file sizes

**rocprof-compute** (kernel-level profiling):
```shell
./mfc.sh run case.py -t simulation --rcu -n <name> [rocprof-compute flags]
```
- Provides rooflines, cache usage, register usage
- Runs the executable multiple times (moderately slower)

> [!NOTE]
> Place profiling arguments at the end of the command so their flags can be appended.

---

## Batch Execution

The MFC detects which scheduler your system is using and handles the creation and execution of batch scripts.
Use batch mode for running on HPC clusters with job schedulers (SLURM, PBS, LSF).

### Basic Usage

```shell
./mfc.sh run case.py -e batch -N 2 -n 4 -c <computer>
```

### Batch Options

| Option | Long Form | Description |
|--------|-----------|-------------|
| `-e batch` | `--engine batch` | Enable batch submission (required) |
| `-N` | `--nodes` | Number of nodes to request |
| `-n` | `--tasks-per-node` | MPI ranks per node |
| `-w` | `--walltime` | Maximum job time (HH:MM:SS) |
| `-a` | `--account` | Account/allocation to charge |
| `-p` | `--partition` | Queue/partition name |
| `-q` | `--qos` | Quality of service level |
| `-@` | `--email` | Email for job notifications |
| `-#` | `--name` | Job name (default: MFC) |
| `-c` | `--computer` | Submission template to use |

### Examples

**Basic 4-node job:**
```shell
./mfc.sh run case.py -e batch -N 4 -n 8 -w 02:00:00
```

**Job with account and email:**
```shell
./mfc.sh run case.py -e batch -N 2 -a myproject -@ user@example.com
```

**GPU job on Frontier:**
```shell
./mfc.sh run case.py -e batch -N 4 -n 8 -c frontier --gpu -w 01:00:00
```

**Dry run (preview script without submitting):**
```shell
./mfc.sh run case.py -e batch -N 2 --dry-run
```

**Wait for job completion:**
```shell
./mfc.sh run case.py -e batch -N 2 --wait
```

### Computer Templates

MFC includes pre-configured templates in `toolchain/templates/` for many clusters.
Use `-c <name>` to select one:

```shell
./mfc.sh run case.py -e batch -c phoenix   # Georgia Tech Phoenix
./mfc.sh run case.py -e batch -c frontier  # OLCF Frontier
./mfc.sh run case.py -e batch -c delta     # NCSA Delta
```

If no template exists for your cluster, use `-c default` and customize as needed,
or contribute a new template.

### Scheduler Notes

**SLURM systems:**
Most clusters use SLURM. MFC automatically generates appropriate `sbatch` scripts.

**LSF systems (e.g., Summit):**
IBM's JSRUN does not use the traditional node-based approach. MFC constructs equivalent resource sets for task and GPU counts.

---

<a name="restarting-cases"></a>
## Restarting Cases

When running a simulation, MFC generates a `./restart_data` folder in the case directory that contains `lustre_*.dat` files that can be used to restart a simulation from saved timesteps.
This allows a user to simulate some timestep $X$, then continue it to run to another timestep $Y$, where $Y > X$.
The user can also choose to add new patches at the intermediate timestep.

If you want to restart a simulation,

- For a simulation that uses a constant time step set up the initial case file with:
    - `t_step_start` : $t_i$
    - `t_step_stop`  : $t_f$
    - `t_step_save`  : $SF$ in which $t_i$ is the starting time, $t_f$ is the final time, and $SF$ is the saving frequency time.
    For a simulation that uses adaptive time-stepping, set up the initial case file with:
    - `n_start` : $t_i$
    - `t_stop`  : $t_f$
    - `t_save`  : $SF$ in which $t_i$ is the starting time, $t_f$ is the final time, and $SF$ is the saving frequency time.

- Run `pre_process` and `simulation` on the case.
    - `./mfc.sh run case.py -t pre_process simulation `
- As the simulation runs, Lustre files will be created for each saved timestep in `./restart_data`.
- When the simulation stops, choose any Lustre file as the restarting point (lustre_ $t_s$.dat)
- Create a new duplicate input file (e.g., `restart_case.py`), which should have:

1. For the Computational Domain Parameters
    - Have the following removed __except__ `m`, `n`, and `p`:
		- All domain/mesh information
			- `(xyz)_domain%%beg`
			- `(xyz)_domain%%end`
			- `stretch_(xyz)`
			- `a_(xyz)`
			- `(xyz)_a`
			- `(xyz)_b`
	- When using a constant time-step, alter the following:
		- `t_step_start` : $t_s$ (the point at which the simulation will restart)
		- `t_step_stop`  : $t_{f2}$ (new final simulation time, which can be the same as $t_f$)
		- `t_step_save`  : ${SF}_2$ (if interested in changing the saving frequency)

        If using a CFL-based time-step, alter the following:
		- `n_start` : $t_s$ (the save file at which the simulation will restart)
		- `t_stop`  : $t_{f2}$ (new final simulation time, which can be the same as $t_f$)
		- `t_save`  : ${SF}_2$ (if interested in changing the saving frequency)

	- Add the following:
		- `old_ic` : 'T' (to specify that we have initial conditions from previous simulations)
		- `old_grid` : 'T' (to specify that we have a grid from previous simulations)
		- `t_step_old` : $t_i$ (the time step used as the `t_step_start` of the original `case.py` file)
2. For the Simulation Algorithm Parameters
	- Substitute `num_patches` to reflect the number of ADDED patches in the `restart_case.py` file. If no patches are added, set `num_patches: 0`

3. For Patches
	- Have all information about old patches (used in the `case.py` file) removed.
		- `patch_icpp(1)%%all variables`
		- `patch_icpp(2)%%all variables`
		- `patch_icpp(num_patches)%%all variables`
	- Add information about new patches that will be introduced, if any. The parameter num_patches should reflect this addition.
		- e.g. `patch_icpp(1)%%some variables of interest`

4. For Fluid properties
	- Keep information about the fluid properties

- Run pre-process and simulation on `restart_case.py`
    - `./mfc.sh run restart_case.py -t pre_process simulation `

- Run the post_process
	- There are several ways to do this. Keep in mind that, regardless of the .py file used, the post_process command will generate output files in the [`t_step_start`, `t_step_stop`] range, with `t_step_save` as the spacing between files.
	- One way is to set `t_step_stop` to the restarting point $t_s$ in `case.py`. Then, run the commands below. The first command will run on timesteps $[t_i, t_s]$. The second command will run on $[t_s, t_{f2}]$. Therefore, the whole range $[t_i, t_{f2}]$ will be post-processed.

```shell
./mfc.sh run case.py -t post_process
./mfc.sh run restart_case.py -t post_process
```

We have provided an example, `case.py` and `restart_case.py` in `/examples/1D_vacuum_restart/`. This simulation is a duplicate of the `1D_vacuum` case. It demonstrates stopping at timestep 7000, adding a new patch, and restarting the simulation. To test this code, run:

```shell
./mfc.sh run examples/1D_vacuum_restart/case.py -t pre_process simulation
./mfc.sh run examples/1D_vacuum_restart/restart_case.py -t pre_process simulation
./mfc.sh run examples/1D_vacuum_restart/case.py -t post_process
./mfc.sh run examples/1D_vacuum_restart/restart_case.py -t post_process
```

---

## Example Runs

- Oak Ridge National Laboratory's [Summit](https://www.olcf.ornl.gov/summit/):

```shell
./mfc.sh run examples/2D_shockbubble/case.py -e batch \
               -N 2 -n 4 -t simulation -a <redacted> -c summit
```
