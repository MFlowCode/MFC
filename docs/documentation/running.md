# Running

MFC can be run using `mfc.sh`'s `run` command. It supports both interactive and
batch execution, the latter being designed for multi-socket systems, namely supercomputers,
equipped with a scheduler such as PBS, SLURM, and LSF. A full (and updated) list
of available arguments can be acquired with `./mfc.sh run -h`. 

## Interactive Execution

To run all stages of MFC, that is [pre_process](src/pre_process/), [simulation](src/simulation/), and [post_process](src/post_process/)   on the sample case [2D_shockbubble](examples/2D_shockbubble/),

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py
```

If you want to run a subset of the available stages, you can use the `-t` argument.
To use multiple threads, use the `-n` option along with the number of threads you wish to use.
If a (re)build is required, it will be done automatically, with the number of threads
specified with the `-j` option.

For example,

- Running [pre_process](src/pre_process/) with 2 cores:

```console
$ ./mfc.sh run examples/2D_shockbubble/case.py -t pre_process -n 2
```

- Running [simulation](src/simulation/) and [post_process](src/post_process/)
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
generate the batch file it will submit. These files are located in the [templates](templates/)
directory. To request GPUs, modification of the template will be required on most systems.

- Lines that begin with `#>` are ignored and won't figure in the final batch
script, not even as a comment.

- Statements of the form `${expression}` are string-replaced to provide
runtime parameters, most notably execution options. They reference the variables in the
same format as those under the "run" section of [defaults.yaml](defaults.yaml), replacing
`-` for `_`. You can perform therein any Python operation recognized by the built-in `expr()` function.

As an example, one might request GPUs on a SLURM system using the following:

```
#SBATCH --gpus=v100-32:{gpus_per_node*nodes}
```

- Statements of the form `{MFC::expression}` tell MFC where to place the common code,
across all batch files, that is required for proper execution. They are not intended to be
modified by users.

**Disclaimer**: IBM's JSRUN on LSF-managed computers does not use the traditional node-based approach to
allocate resources. Therefore, the MFC constructs equivalent resource-sets in task and GPU count.

## Example Runs

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
