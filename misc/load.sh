#!/usr/bin/env bash

RED="\u001b[31m"
CYAN="\u001b[36m"
GREEN="\u001b[32m"
YELLOW="\u001b[33m"
COLOR_RESET="\033[m"

echo -en "$YELLOW"
echo -e "Please source this script:"
echo -e "$ source load.sh\n"
echo -en "$COLOR_RESET"

on_error() {
    echo -en "\n$YELLOW"
    echo -e "We would welcome your contribution:"
    echo -e " - Github: https://github.com/MFlowCode/MFC"
    echo -en "$COLOR_RESET"

    exit 1
}

LLNL="$GREEN""LLNL$COLOR_RESET"
SUMMIT="$LLNL/Summit"
ASCENT="$LLNL/Ascent"
WOMBAT="$LLNL/Wombat"

XSEDE="$CYAN""XSEDE$COLOR_RESET"
BRIDGES2="$XSEDE/Bridges2"
EXPANSE="$XSEDE/Expanse"

RICHARDSON="Richardson"

C_LLNL=$GREEN"s$COLOR_RESET/$GREEN""a$COLOR_RESET/$GREEN""w$COLOR_RESET"
C_XSEDE=$CYAN"b$COLOR_RESET/$CYAN""e$COLOR_RESET"
C_OTHER="r"

echo "[MFC] Which computer would you like to load submodules for?"
echo -e " - $SUMMIT (s) - $BRIDGES2 (b) - $RICHARDSON (r)"
echo -e " - $ASCENT (a) - $EXPANSE  (e)"
echo -e " - $WOMBAT (w)"
echo -en "($C_LLNL/$C_XSEDE/$C_OTHER): "
read response

COMPUTER=""
declare -a MODULES

if [ "$response" == "s" ]; then # For Summit
    COMPUTER="$SUMMIT"
    MODULES=("lsf-tools/2.0"
             "darshan-runtime/3.3.1-lite"
             "spectrum-mpi/10.4.0.3-20210112"
             "hsi/5.0.2.p5"
             "xalt/1.2.1"
             "nvhpc/21.9"
             "cuda/11.4.2")
elif [ "$response" == "b" ]; then # Bridges2
    COMPUTER="$BRIDGES2"
    MODULES=("openmpi/4.0.5-nvhpc21.7"
             "cuda/11.1.1")
elif [ "$response" == "a" ]; then # For Ascent
    COMPUTER="$ASCENT"
    MODULES=("nvhpc/21.9"
             "spectrum-mpi"
             "cuda/11.2.0"
             "nsight-compute"
             "nsight-systems")
elif [ "$response" == "r" ]; then # Richardson
    COMPUTER="$RICHARDSON"
    MODULES=("gcc/9.3.0"
             "openmpi-2.0/gcc-9.3.0"
             "python/3.7")
elif [ "$response" == "w" ]; then # For Wombat
    COMPUTER="$WOMBAT"
    MODULES=("python/3.9.9"
             "/sw/wombat/Nvidia_HPC_SDK/modulefiles/nvhpc/22.1")
elif [ "$response" == "e" ]; then # Expanse
    COMPUTER="$EXPANSE"
    MODULES=("gpu/0.15.4"
             "python/3.8.5"
             "openmpi/4.0.5"
             "cuda/11.0.2"
             "nvhpc/22.2"
             "slurm/expanse/current")
else
    echo -e $RED"Error: Requested system is not supported"$COLOR_RESET
    on_error
fi

echo -e "[MFC] Loading modules for \u001b[35m$COMPUTER$COLOR_RESET:"

module reset

# Find length of longest module_name in $MODULES for $COMPUTER
max_module_length="0"
for module_name in ${MODULES[@]}; do
    module_length="${#module_name}"

    if [ "$module_length" -gt "$max_module_length" ]; then
        max_module_length="$module_length"
    fi
done

# Load modules for $COMPUTER with $MODULES
for module_name in ${MODULES[@]}; do
    echo -en " - Loading $CYAN$module_name$COLOR_RESET "

    # Add padding spaces
    module_length="${#module_name}"
    delta="$((max_module_length-module_length-2))"
    if [ "$delta" -ne "-2" ]; then
        printf "%0.s-" $(seq 0 $delta)
        echo -n " "
    fi

    # Load the module
    module load "$module_name" > /dev/null 2>&1

    # Handle Success / Failure
    code=$?
    if [ "$code" == "0" ]; then
        echo -e "[$GREEN""SUCCESSFUL$COLOR_RESET]"
    else
        echo -e "[$RED""FAILED$COLOR_RESET]"
        
        # Run load again to show error message
        module load "$module_name"

        on_error
    fi
done

echo -e "[MFC] You should now be able to build MFC on $COMPUTER!"
