#!/usr/bin/env bash

RED="\u001b[31m"
CYAN="\u001b[36m"
GREEN="\u001b[32m"
YELLOW="\u001b[33m"
MAGENTA="\u001b[35m"
COLOR_RESET="\033[m"

echo -en "$YELLOW"
echo -e "Please source this script:"
echo -e "$(whoami) $ . ./mfc.sh load\n"
echo -en "$COLOR_RESET"

on_error() {
    echo -en "\n$YELLOW"
    echo -e "We would welcome your contribution:"
    echo -e " - Github: https://github.com/MFlowCode/MFC"
    echo -en "$COLOR_RESET"
    return
}

ORNL="$GREEN""ORNL$COLOR_RESET"
SUMMIT="$ORNL/Summit"
ASCENT="$ORNL/Ascent"
WOMBAT="$ORNL/Wombat"

XSEDE="$CYAN""XSEDE$COLOR_RESET"
BRIDGES2="$XSEDE/Bridges2"
EXPANSE="$XSEDE/Expanse"

GT="$YELLOW""GATECH$COLOR_RESET"
PHOENIX="$GT/Phoenix"

CALTECH=$RED"CALTECH$COLOR_RESET"
RICHARDSON="$CALTECH/Richardson"

C_ORNL=$GREEN"s$COLOR_RESET/$GREEN""a$COLOR_RESET/$GREEN""w$COLOR_RESET"
C_XSEDE=$CYAN"b$COLOR_RESET/$CYAN""e$COLOR_RESET"
C_GT=$YELLOW"p$COLOR_RESET"
C_OTHER=$RED"r$COLOR_RESET"

# Reset u_computer & u_cg to known values since this script is run with "source"
# Therefore, values of variables defined here are kept when the script runs again.
u_computer=""
u_cg=""

# If there are command-line arguments, parse them:
while [[ $# -gt 0 ]]; do
    case $1 in
    -c|--computer)
      u_computer="$2"
      shift; shift
      ;;
    -m|--mode)
      u_cg="$2"
      shift; shift
      ;;
    -*|--*)
        echo "Unknown option $1"
        on_error
        return
        ;;
    esac
done

# Get computer (if not supplied in command-line)
if [ -v $u_computer ]; then
    echo -e "$MAGENTA[Q 1/2]$COLOR_RESET Which computer would you like to load submodules for?"
    echo -e " | $SUMMIT ($GREEN""s$COLOR_RESET) | $BRIDGES2 ($CYAN""b$COLOR_RESET) | $RICHARDSON ($RED""r$COLOR_RESET)"
    echo -e " | $ASCENT ($GREEN""a$COLOR_RESET) | $EXPANSE  ($CYAN""e$COLOR_RESET) |"
    echo -e " | $WOMBAT ($GREEN""w$COLOR_RESET) | $PHOENIX ($YELLOW""p$COLOR_RESET) |"
    echo -en "($C_ORNL/$C_XSEDE/$C_GT/$C_OTHER): "
    read u_computer
fi

# Get CPU/GPU (if not supplied in command-line)
if [ -v $u_cg ]; then
    echo -e "$MAGENTA[Q 2/2]$COLOR_RESET Would you like to run solely on CPUs or GPUs as well?"
    echo -e " - CPU (c) - GPU (g)"
    echo -n "(c/g): " 
    read u_cg
fi

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
    return
fi

if [ "$u_computer" == "s" ]; then # For Summit
    COMPUTER="$SUMMIT"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/11.1.0" "spectrum-mpi/10.4.0.3-20210112")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("spectrum-mpi/10.4.0.3-20210112" "nvhpc/21.9" "cuda/11.4.2")
    fi

    MODULES=("${MODULES[@]}" "python/3.8.10" "darshan-runtime/3.3.1-lite"
             "hsi/5.0.2.p5"  "xalt/1.2.1"    "lsf-tools/2.0")
elif [ "$u_computer" == "b" ]; then # Bridges2
    COMPUTER="$BRIDGES2"

    if [ "$u_cg" == "c" ]; then
        MODULES=("intel/2021.3.0" "intelmpi/2021.3.0-intel2021.3.0")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("nvhpc/22.1" "cuda/11.1.1")
    fi

    MODULES=("${MODULES[@]}" "openmpi/4.0.5-gcc10.2.0" "python/3.8.6")
elif [ "$u_computer" == "a" ]; then # For Ascent
    COMPUTER="$ASCENT"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/11.1.0" "spectrum-mpi")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("nvhpc/21.11"    "spectrum-mpi"   "cuda/nvhpc"
                 "nsight-compute" "nsight-systems")
    fi
elif [ "$u_computer" == "r" ]; then # Richardson
    COMPUTER="$RICHARDSON"

    if [ "$u_cg" == "c" ]; then
        MODULES=("gcc/9.3.0" "openmpi-2.0/gcc-9.3.0")
    elif [ "$u_cg" == "g" ]; then
        echo -e $RED"Error: GPU not supported on Richardson."$COLOR_RESET
        on_error
        return
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
        MODULES=("cpu/0.15.4" "intel/19.1.1.217" "intel-mpi/2019.8.254")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("gpu/0.15.4" "openmpi/4.0.5" "cuda/11.0.2" "nvhpc/22.2")
    fi

    MODULES=("${MODULES[@]}" "slurm/expanse/current")
elif [ "$u_computer" == "p" ]; then # Phoenix
    COMPUTER="$PHOENIX"

    if [ "$u_cg" == "c" ]; then
        MODULES=("intel/19.0.5" "mvapich2/2.3.2")
    elif [ "$u_cg" == "g" ]; then
        MODULES=("cuda/11.2" "nvhpc/22.1")
    fi

    MODULES=("${MODULES[@]}" "python/3.7.4")
else
    echo -e $RED"Error: Requested system $u_computer is not supported (yet!)"$COLOR_RESET
    
    on_error
    return
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
        return
    fi
done

echo -e "$MAGENTA[ MFC ]$COLOR_RESET You should now be able to build MFC on $COMPUTER ($CG mode)!"
