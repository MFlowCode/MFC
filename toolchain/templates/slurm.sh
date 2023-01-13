#!/usr/bin/env bash
#>
#>           - SLURM Batch File Template -
#>
#> This file is part of the ./mfc.sh run subsystem.
#> For more information, please consult the README.
#>
#> - You are invited to modify this file to suit your
#>   needs, in order to get MFC running properly on
#>   your system.
#>
#> - Lines that begin with "#>" are ignored and won't
#>   figure in the final batch script, not even as a
#>   comment.
#>
#> - Statements of the form `${expression}` are string-
#>   -replaced by mfc.sh run to provide runtime parameters,
#>   most notably execution options. They reference the
#>   variables in the same format as those under the "run"
#>   section of [mfc.user.yaml](mfc.user.yaml), replacing
#>   "-" for "_". You can perform therein any Python operation
#>   recognized by the built-in `expr()` function.
#>
#> - Statements of the form {MFC::expression} tell MFC
#>   where to place the common code, across all batch
#>   files that is required to run MFC. They are not
#>   intended to be modified by users.
#>
#SBATCH --job-name="{name}"
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --gpus=v100-16:{(1 if gpu else 0)*tasks_per_node*nodes}
#SBATCH --time={walltime}
#SBATCH --partition="{partition}"
#SBATCH --output="{name}.out"
#SBATCH --account="{account}"
#SBATCH --error="{name}.err"
#SBATCH --mail-user="{email}"
#SBATCH --export=ALL
#SBATCH --mail-type="BEGIN, END, FAIL"
#>
#> Note: The following options aren't enabled by default.
#>       They serve as a guide to users that wish to pass
#>       more options to the batch system.
#>
#> #SBATCH --mem=...
#> #SBATCH --constraint="lustre"
#> #SBATCH --gpus-per-task={1 if gpu else 0}


#>
#> Note: If your system requires you to load environment
#>       modules inside of your batch script, please load
#>       them bellow.
#>


#>
#> Note: The MFC prologue sets up the environment required
#>       prior to execution.
#>
{MFC::PROLOGUE}


#>
#> Note: This MPI executable might not be well supported
#>       on your system - if at all. {MFC::BIN} refers to
#>       the path the MFC executable.
#>

for binpath in {MFC::BINARIES}; do

    echo -e ":) Running $binpath:"

#>
#> Note: This MPI executable might not be well supported
#>       on your system - if at all. {MFC::BIN} refers to
#>       the path the MFC executable.
#>
#>srun                                   \
#>     --nodes={nodes}                   \
#>     --ntasks-per-node {cpus_per_node} \
#>     --mpi=pmi2		       \
#>     "{MFC::BIN}"
#>
#> srun --mpi=pmix   \
#>      "{MFC::BIN}"
#>

    mpirun -np {nodes*tasks_per_node} "$binpath"

done

{MFC::EPILOGUE}

#>
#> Note: Lines after the MFC Epilogue will not be executed.
#>
