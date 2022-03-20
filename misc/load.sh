#!/usr/bin/env bash

RED="\u001b[31m"
CYAN="\u001b[36m"
GREEN="\u001b[32m"
YELLOW="\u001b[33m"
MAGENTA="\u001b[35m"
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

echo -e "$MAGENTA[Q 1/2]$COLOR_RESET Which computer would you like to load submodules for?"
echo -e " - $SUMMIT (s) - $BRIDGES2 (b) - $RICHARDSON (r)"
echo -e " - $ASCENT (a) - $EXPANSE  (e)"
echo -e " - $WOMBAT (w)"
echo -en "($C_LLNL/$C_XSEDE/$C_OTHER): "
read u_computer

echo -e "$MAGENTA[Q 2/2]$COLOR_RESET Would you like to run solely on CPUs or GPUs as well?"
echo -e " - CPU (c) - GPU (g)"
echo -n "(c/g): " 
read u_cg

# User input to lowercase
u_computer=$(echo "$u_computer" | tr '[:upper:]' '[:lower:]')
u_cg=$(echo "$u_cg" | tr '[:upper:]' '[:lower:]')

COMPUTER=""
CG=""
declare -a MODULES

if [ "$u_cg" == "c" ]; then # for CPU
    CG="CPU"
elif [ "$u_cg" == "g" ]; then # For GPU
    CG="GPU"
else
    echo -e $RED"Error: Invalid choice \"$u_cg\" for Q2."$COLOR_RESET

    on_error
fi

if [ "$u_computer" == "s" ]; then # For Summit
    COMPUTER="$SUMMIT"

    if [ "$u_cg" == "c" ]; then
        echo -e $RED"Error: CPU modules not (yet) supported on Summit."$COLOR_RESET
        on_error
    elif [ "$u_cg" == "g" ]; then
        MODULES=("lsf-tools/2.0"
                 "darshan-runtime/3.3.1-lite"
                 "spectrum-mpi/10.4.0.3-20210112"
                 "hsi/5.0.2.p5"
                 "xalt/1.2.1"
                 "nvhpc/21.9"
                 "cuda/11.4.2")
    fi
elif [ "$u_computer" == "b" ]; then # Bridges2
    COMPUTER="$BRIDGES2"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/10.2.0" "openmpi/4.0.5-gcc10.2.0")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("nvhpc/22.1" "openmpi/4.0.5-nvhpc22.1" "cuda/11.1.1")
    fi

    MODULES=("${MODULES[@]}" "python/3.8.6")
elif [ "$u_computer" == "a" ]; then # For Ascent
    COMPUTER="$ASCENT"

    if [ "$u_cg" == "c" ]; then
        echo -e $RED"Error: CPU modules not (yet) supported on Ascent."$COLOR_RESET
        on_error
    elif [ "$u_cg" == "g" ]; then
        MODULES=("nvhpc/21.9"     "spectrum-mpi"   "cuda/11.2.0"
                 "nsight-compute" "nsight-systems")
    fi
    
elif [ "$u_computer" == "r" ]; then # Richardson
    COMPUTER="$RICHARDSON"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/9.3.0" "openmpi-2.0/gcc-9.3.0")
    elif [ "$u_cg" == "g" ]; then
        echo -e $RED"Error: GPU not supported on Richardson."$COLOR_RESET
        on_error
    fi
    
    MODULES=("${MODULES[@]}" "python/3.7")
elif [ "$u_computer" == "w" ]; then # For Wombat
    COMPUTER="$WOMBAT"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/11.1.0" "openmpi/4.0.5_gcc")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("cuda/11.5.1" "/sw/wombat/Nvidia_HPC_SDK/modulefiles/nvhpc/22.1")
    fi

    MODULES=("${MODULES[@]}" "python/3.9.9")
elif [ "$u_computer" == "e" ]; then # Expanse
    COMPUTER="$EXPANSE"

    if [ "$u_cg" == "c" ]; then
        MODULES=("cpu/0.15.4" "gcc/10.2.0" "openmpi/gcc/64/1.10.7")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("gpu/0.15.4" "openmpi/4.0.5" "cuda/11.0.2" "nvhpc/22.2")
    fi

    MODULES=("${MODULES[@]}" "python/3.8.5" "slurm/expanse/current")
else
    echo -e $RED"Error: Requested system $u_computer is not supported (yet!)"$COLOR_RESET
    
    on_error
fi

echo -e "$MAGENTA[ MFC ]$COLOR_RESET Loading modules for $COMPUTER ($CG mode):"

# Reset modules to default system configuration
module reset > /dev/null 2>&1
code="$?"

# Purge if reset is not available
if [ "$code" -ne "0" ]; then
    module purge > /dev/null 2>&1
fi

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
