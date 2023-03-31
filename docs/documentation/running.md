# Running

MFC can be run using `mfc.sh`'s `run` command. It supports both interactive and
batch execution, the latter being designed for multi-socket systems, namely supercomputers,
equipped with a scheduler such as PBS, SLURM, and LSF. A full (and updated) list
of available arguments can be acquired with `./mfc.sh run -h`. 

## Interactive Execution

To run all stages of MFC, that is [pre_process](https://github.com/MFlowCode/MFC/tree/master/src/pre_process/), [simulation](https://github.com/MFlowCode/MFC/tree/master/src/simulation/), and [post_process](https://github.com/MFlowCode/MFC/tree/master/src/post_process/) on the sample case [2D_shockbubble](https://github.com/MFlowCode/MFC/tree/master/examples/2D_shockbubble/),

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py
```

If you want to run a subset of the available stages, you can use the `-t` argument.
To use multiple threads, use the `-n` option along with the number of threads you wish to use.
If a (re)build is required, it will be done automatically, with the number of threads
specified with the `-j` option.

For example,

- Running [pre_process](https://github.com/MFlowCode/MFC/tree/master/src/pre_process/) with 2 cores:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -t pre_process -n 2
```

- Running [simulation](https://github.com/MFlowCode/MFC/tree/master/src/simulation/) and [post_process](https://github.com/MFlowCode/MFC/tree/master/src/post_process/)
using 4 cores:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -t simulation post_process -n 4
```

On some computer clusters, MFC might select the wrong MPI program to execute your application
because it uses a general heuristic for its selection. Notably, `srun` is known to fail on some SLURM
systems when using GPUs or MPI implementations from different vendors, whereas `mpirun` functions properly. To override and manually specify which
MPI program you wish to run your application with, please use the `-b <program name>` option (i.e `--binary`).

Additional flags can be appended to the MPI executable call using the `-f` (i.e `--flags`) option.

Please refer to `./mfc.sh run -h` for a complete list of arguments and options, along with their defaults.

## Batch Execution

The MFC detects which scheduler your system is using and handles the creation and
execution of batch scripts. The batch engine is requested with the `-e batch` option.
Whereas the interactive engine can execute all of MFC's codes in succession, the batch engine
requires you to only specify one target with the `-t` option. The number of nodes and GPUs can, 
respectively be specified with the `-N` (i.e `--nodes`) and `-g` (i.e `--gpus-per-node`) options.

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch -N 2 -n 4 -g 4 -t simulation
```

Other useful arguments include:

- `-# <job name>` to name your job. (i.e `--name`)
- `-@ sample@example.com` to receive emails from the scheduler. (i.e `--email`)
- `-w hh:mm:ss` to specify the job's maximum allowed walltime. (i.e `--walltime`)
- `-a <account name>` to identify the account to be charged for the job. (i.e `--account`)
- `-p <partition name>` to select the job's partition. (i.e `--partition`)

Since some schedulers don't have a standardized syntax to request certain resources, MFC can only
provide support for a restricted subset of common configuration options. If MFC fails
to execute on your system, or if you wish to adjust how the program runs and resources
are requested to be allocated, you are invited to modify the template batch script for your queue system.
Upon execution of `./mfc.sh run`, MFC fills in the template with runtime parameters, to
generate the batch file it will submit. These files are located in the [templates](https://github.com/MFlowCode/MFC/tree/master/toolchain/templates/)
directory. To request GPUs, modification of the template will be required on most systems.

- Lines that begin with `#>` are ignored and won't figure in the final batch
script, not even as a comment.

- Statements of the form `${expression}` are string-replaced to provide
runtime parameters, most notably execution options. You can perform therein any Python operation recognized by the built-in `expr()` function.

As an example, one might request GPUs on a SLURM system using the following:

```
#SBATCH --gpus=v100-32:{gpus_per_node*nodes}
```

- Statements of the form `{MFC::expression}` tell MFC where to place the common code,
across all batch files, that is required for proper execution. They are not intended to be
modified by users.

**Disclaimer**: IBM's JSRUN on LSF-managed computers does not use the traditional node-based approach to
allocate resources. Therefore, the MFC constructs equivalent resource-sets in task and GPU count.

### Profiling with NVIDIA Nsight

MFC provides two different argument to facilitate profiling with NVIDIA Nsight. **Please ensure that the used argument is placed at the end so that their respective flags can be appended.**
- Nsight Systems (Nsys): `./mfc.sh run ... --nsys [nsys flags]` allows one to visualize MFC's system-wide performance with [NVIDIA Nsight Systems](https://developer.nvidia.com/nsight-systems). NSys is best for getting a general understanding of the order and execution times of major subroutines (WENO, Riemann, etc.) in MFC. When used, `--nsys` will run the simulation and generate `.nsys-rep` files in the case directory for all targets. These files can then be imported into Nsight System's GUI, which can be downloaded [here](https://developer.nvidia.com/nsight-systems/get-started#latest-Platforms). It is best to run case files with a few timesteps so that the report files remain small. Learn more about NVIDIA Nsight Systems [here](https://docs.nvidia.com/nsight-systems/UserGuide/index.html).
- Nsight Compute (NCU): `./mfc.sh run ... --ncu [ncu flags]` allows one to conduct kernel-level profiling with [NVIDIA Nsight Compute](https://developer.nvidia.com/nsight-compute). NCU provides profiling information for every subroutine called and is more detailed than NSys. When used, `--ncu` will output profiling information for all subroutines, including elapsed clock cycles, memory used, and more after the simulation is run. Please note that adding this argument will significantly slow down the simulation and should only be used on case files with a few timesteps. Learn more about NVIDIA Nsight Compute [here](https://docs.nvidia.com/nsight-compute/NsightCompute/index.html).

### Restarting Cases

When running a simulation, MFC generates a `./restart_data` folder in the case directory that contains `lustre_*.dat` files that can be used to restart a simulation from saved timesteps. This allows a user to run a simulation to some timestep $X$, then later continue it to run to another timestep $Y$, where $Y > X$. The user can also choose to add new patches at the intermediate timestep.

If you want to restart a simulation, 

- Set up the initial simulation, with:
    - `t_step_start` : $t_i$
    - `t_step_stop`  : $t_f$
    - `t_step_save`  : $SF$
in which $t_i$ is the starting time, $t_f$ is the final time, and $SF$ is the saving frequency time.
- Run pre-process and simulation on the case.
    - `./mfc.sh run case.py -t pre_process simulation `
- As the simulation runs, it will create LUSTRE files for each saved timestep in `./restart_data`.
- When the simulation stops, choose any LUSTRE file as the restarting point (lustre_ $t_s$.dat)
- Create a new duplicate input file, (ex. `restart_case.py`), on which it should:

1. For the Computational Domain Parameters 
    - Have the following removed BUT `m`, `n`, and `p`:
		- All domaing/mesh information
			- `(xyz)_domain%beg`
			- `(xyz)_domain%end`
			- `stretch_(xyz)`
			- `a_(xyz)`
			- `(xyz)_a`
			- `(xyz)_b`
	- Have the following altered:
		- `t_step_start` : $t_s$ # this is the point at which the simulation will restart
		- `t_step_stop`  : $t_{f2}$ # your new final simulation time, which can be the same as $t_f$
		- `t_step_save`  : ${SF}_2$ # if interested in changing the saving frequency 
	- Have the following ADDED:
		- `old_ic` : 'T' # to specify that we have initial conditions from previous simulations
		- `old_grid` : 'T' # to specify that we have a grid from previous simulations (maybe I do not need m, n, and p, then?)
		- `t_step_old` : $t_i$ # this is the time step used as the `t_step_start` of the original `case.py` file
2. For the Simulation Algorithm Parameters
	- Substitute `num_patches` to reflect the number of ADDED patches in the `restart_case.py` file. If no patches are added, set `num_patches: 0`

3. For Patches
	- Have all information about old patches (used in the `case.py` file) removed.
		- `patch_icpp(1)%all variables`
		- `patch_icpp(2)%all variables`
		- `patch_icpp(num_patches)%all variables`
	- Add information about new patches that will be introduced, if any. The parameter num_patches should reflect this addition.
		- e.g. `patch_icpp(1)%some variables of interest`	

4. For Fluid properties	
	- Keep information about the fluid properties

- Run pre-process and simulation on restart_case.py
    - `./mfc.sh run restart_case.py -t pre_process simulation `
	
- Run the post_process
	- There are several ways to do this. Keep in mind that, regardless of the .py file used, the post_process command will generate output files in the [`t_step_start`, `t_step_stop`] range, with `t_step_save` as the spacing between files.
	- One way is to set `t_step_stop` to the restarting point $t_s$ in `case.py`. Then, run the commands below. The first command will run on timesteps $[t_i, t_s]$. The second command will run on $[t_s, t_{f2}]$. Therefore, the whole range $[t_i, t_{f2}]$ will be post processed.

```console
$ ./mfc.sh run case.py -t post_process
$ ./mfc.sh run restart_case.py -t post_process
```	

We have provided an example `case.py` and `restart_case.py` in `/examples/1D_vacuum_restart/`. This simulation is a duplicate of the `1D_vacuum` case. It demonstrates stopping at timestep 7000, adding a new patch, and restarting the simulation. To test this code, run:

```console
$ ./mfc.sh run examples/1D_vacuum_restart/case.py -t pre_process simulation
$ ./mfc.sh run examples/1D_vacuum_restart/restart_case.py -t pre_process simulation
$ ./mfc.sh run examples/1D_vacuum_restart/case.py -t post_process
$ ./mfc.sh run examples/1D_vacuum_restart/restart_case.py -t post_process
```

### Example Runs

- Oak Ridge National Laboratory's [Summit](https://www.olcf.ornl.gov/summit/):

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch    \
               -N 2 -n 4 -g 4 -t simulation -a <redacted>
```

- University of California, San Diego's [Expanse](https://www.sdsc.edu/services/hpc/expanse/):

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -e batch -p GPU -t simulation \
               -N 2 -n 8 -g 8 -f="--gpus=v100-32:16" -b mpirun –w 00:30:00
```
