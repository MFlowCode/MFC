#!/usr/bin/env bash
#>
#>           - LSF Batch File Template -
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
#BSUB -J {name}
#BSUB -nnodes {nodes}
#BSUB -N
#BSUB -P {account}
#BSUB -W {"walltime"[:-3]}
#>
#> Note: The above expression for the walltime converts
#>       the expression "hh:mm:ss" to the appropriate
#>       format for the batch system ("hh:mm"). It is
#>       a python expression evaluated at runtime.
#>
#>
#> Note: The following options aren't enabled by default.
#>       They serve as a guide to users that wish to pass
#>       more options to the batch system.
#>



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
jsrun {'--smpiargs="-gpu"' if gpus_per_node > 0 else ''} \
      --nrs          {cpus_per_node*nodes}               \
      --cpu_per_rs   1                                   \
      --gpu_per_rs   {min(gpus_per_node, 1)}             \
      --tasks_per_rs 1                                   \
      "{MFC::BIN}"

{MFC::EPILOGUE}
#>
#> Note: Lines after the MFC Epilogue will not be executed.
#>
