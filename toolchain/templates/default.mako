<%!
import os
from mako.exceptions import RuntimeException
%>
% if os.name == 'nt':
@echo off
% else:
#!/usr/bin/env bash
% endif

<%namespace name="helpers" file="helpers.mako"/>

${helpers.template_prologue()}

<%
if engine == 'batch':
    raise RuntimeException("The default template does not support batch jobs. Please use a different template via the --computer option.")
%>

% if os.name != 'nt':
    warn "This is the$MAGENTA default$COLOR_RESET template."
    warn "It is not intended to support all systems and execution engines."
    warn "Consider using a different template via the $MAGENTA--computer$COLOR_RESET option if you encounter problems."

    % if mpi:
        # Find a suitable MPI launcher and store it in the variable "binary".
        for binary in ${binary or ''} jsrun srun mpirun mpiexec; do
            if command -v $binary > /dev/null; then
                break
            fi
        done

        if ! command -v $binary > /dev/null; then
            error ":( Could not find a suitable MPI launcher.\n"
            exit 1
        else
            ok ":) Selected MPI launcher $MAGENTA$binary$COLOR_RESET. Use$MAGENTA --binary$COLOR_RESET to override."
        fi
    % endif

    % for target in targets:
        ${helpers.run_prologue(target)}

        % if not mpi:
            (set -x; ${profiler} "${target.get_install_binpath(case)}")
        % else:
            if [ "$binary" == "jsrun" ]; then
                (set -x; ${profiler}   \
                    jsrun --nrs          ${tasks_per_node*nodes} \
                        --cpu_per_rs   1                       \
                        --gpu_per_rs   ${1 if gpu_enabled else 0}      \
                        --tasks_per_rs 1                       \
                        "${target.get_install_binpath(case)}")
            elif [ "$binary" == "srun" ]; then
                (set -x; ${profiler}                            \
                    srun --ntasks ${nodes*tasks_per_node}       \
                        "${target.get_install_binpath(case)}")
            elif [ "$binary" == "mpirun" ]; then
                (set -x; ${profiler}     \
                    $binary -np ${nodes*tasks_per_node}            \
                            "${target.get_install_binpath(case)}")
            elif [ "$binary" == "mpiexec" ]; then
                (set -x; ${profiler}                               \
                    $binary --ntasks ${nodes*tasks_per_node}       \
                            "${target.get_install_binpath(case)}")
            fi
        % endif

        ${helpers.run_epilogue(target)}

        echo
    % endfor
% else:
    % for target in targets:
        ${helpers.run_prologue(target)}

        % if not mpi:
            ${profiler} "${target.get_install_binpath(case)}.exe"
        % else:
            ${profiler}                                            \
                mpiexec -n ${nodes*tasks_per_node}                 \
                        "${target.get_install_binpath(case)}.exe"
        % endif

        ${helpers.run_epilogue(target)}
    % endfor
% endif

${helpers.template_epilogue()}
