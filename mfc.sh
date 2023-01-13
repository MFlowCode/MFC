#!/bin/bash


MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=6

# Check whether this script was called from MFC's root directory.
if [ ! -f "$(pwd)/toolchain/mfc.py" ]; then
    echo "mfc: ERROR > You must call this script from within MFC's root folder."

    exit 1
fi

# Load utility script
source "$(pwd)/toolchain/util.sh"

# Handle upgrading from older MFC build systems
if [ -d "$(pwd)/bootstrap" ] || [ -d "$(pwd)/dependencies" ] || [ -f "$(pwd)/build/mfc.lock.yaml" ]; then
    error "You are upgrading from an older version of$MAGENTA MFC$COLOR_RESET. Please remove, if applicable, the dependencies/, bootstrap/, and build/ directories before running this command again."

    exit 1
fi

# If the user wishes to run the "load" script
if [ "$1" == "load" ]; then
    shift;

    # Reset u_computer & u_cg to known values since this script is run with "source"
    # Therefore, values of variables defined here are kept when the script runs again.
    u_computer=""; u_cg=""

    # If there are command-line arguments, parse them:
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--computer) u_computer="$2"; shift; shift;    ;;
            -m|--mode)     u_cg="$2";       shift; shift;    ;;
            -*|--*)        echo "Unknown option $1"; return; ;;
        esac
    done

    # Get computer (if not supplied in command-line)
    if [ -v $u_computer ]; then
        log   "Select a system:"
        log   "$G""ORNL$W:    Ascent     (a), Crusher (c), Summit (s), Wombat (w)"
        log   "$C""ACCESS$W:  Bridges2   (b), Expanse (e)"
        log   "$Y""GaTech$W:  Phoenix    (p)"
        log   "$R""CALTECH$W: Richardson (r)"
        log_n "($G""a$W/$G""c$W/$G""s$W/$G""w$W/$C""b$W/$C""e$CR/$Y""p$CR/$R""r$CR): "
        read u_computer
        log
    fi

    # Get CPU/GPU (if not supplied in command-line)
    if [ -v $u_cg ]; then
        log   "Select configuration:"
        log   " - CPU (c)"
        log   " - GPU (g)"
        log_n "(c/g): "
        read u_cg
        log
    fi

    # User input to lowercase
    u_computer=$(echo "$u_computer" | tr '[:upper:]' '[:lower:]')
    u_cg=$(echo "$u_cg" | tr '[:upper:]' '[:lower:]')

    COMPUTER=""; CG=""
    declare -a MODULES

    if [ "$u_cg" == "c" ]; then # for CPU
        CG="CPU"
    elif [ "$u_cg" == "g" ]; then # For GPU
        CG="GPU"
    fi

    if [ "$u_computer" == "s" ]; then # For Summit
        if [ "$u_cg" == "c" ]; then
            MODULES=("gcc/12.1.0")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("nvhpc/22.5" "cuda/nvhpc")
        fi

        MODULES=("${MODULES[@]}" "python/3.8.10" "darshan-runtime/3.3.1-lite"
                "hsi/5.0.2.p5" "xalt/1.2.1" "lsf-tools/2.0"
                "cmake/3.23.1" "ninja/1.10.2" "spectrum-mpi/10.4.0.3-20210112")
    elif [ "$u_computer" == "b" ]; then # Bridges2
        if [ "$u_cg" == "c" ]; then
            MODULES=("allocations/1.0" "gcc/10.2.0" "python/3.8.6"
                     "openmpi/4.0.5-gcc10.2.0")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("nvhpc/22.9" "openmpi/4.0.5-nvhpc22.9")
        fi

        MODULES=("${MODULES[@]}" "python/3.8.6")
    elif [ "$u_computer" == "a" ]; then # For Ascent
        if [ "$u_cg" == "c" ]; then
            MODULES=("gcc/11.1.0" "spectrum-mpi" "cuda")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("nvhpc/21.11" "spectrum-mpi" "cuda/nvhpc"
                    "nsight-compute" "nsight-systems")
        fi

        MODULES=("${MODULES[@]}" "python" "cmake/3.22.2")
    elif [ "$u_computer" == "r" ]; then # Richardson
        if [ "$u_cg" == "c" ]; then
            MODULES=("gcc/9.3.0" "openmpi-2.0/gcc-9.3.0")
        elif [ "$u_cg" == "g" ]; then
            error "GPU not supported on Richardson."

            return
        fi

        MODULES=("${MODULES[@]}" "python/3.7")
    elif [ "$u_computer" == "w" ]; then # For Wombat
        if [ "$u_cg" == "c" ]; then
            MODULES=("gcc/11.1.0" "openmpi/4.0.5_gcc")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("cuda/11.5.1" "/sw/wombat/Nvidia_HPC_SDK/modulefiles/nvhpc/22.1")
        fi

        MODULES=("${MODULES[@]}" "cmake/3.22.1" "python/3.9.9")
    elif [ "$u_computer" == "e" ]; then # Expanse
        if [ "$u_cg" == "c" ]; then
            warn "Please set CC=icc, CXX=icx, and FC=ifort."
            log

            MODULES=("cpu/0.15.4" "gcc/9.2.0" "openmpi/4.1.1" "cmake/3.18.2")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("gpu/0.15.4" "cuda/11.0.2" "nvhpc/22.2" "openmpi/4.0.5" "cmake/3.19.8")
        fi

        MODULES=("${MODULES[@]}" "python/3.8.5")
    elif [ "$u_computer" == "p" ]; then # Phoenix
        if [ "$u_cg" == "c" ]; then
            MODULES=("intel/19.0.5" "mvapich2/2.3.2")
        elif [ "$u_cg" == "g" ]; then
            MODULES=("cuda/11.2" "nvhpc/22.1")
        fi

        MODULES=("${MODULES[@]}" "python/3.7.4" "cmake/3.20.3")
    elif [ "$u_computer" == "c" ]; then # Crusher
        if [ "$u_cg" == "c" ]; then
            MODULES=()
        elif [ "$u_cg" == "g" ]; then
            MODULES=("rocm/5.1.0" "craype-accel-amd-gfx90a")
        fi

        MODULES=("${MODULES[@]}" "cmake/3.23.2" "cray-fftw/3.3.10.2" "hdf5/1.12.1" "cray-python/3.9.13.1" "ninja/1.10.2" "cray-mpich/8.1.23")
    else
        echo -e $RED"Error: Requested system $u_computer is not supported (yet!)"$COLOR_RESET

        return
    fi

    log "Loading modules for $CG mode:"

    # Reset modules to default system configuration
    if [ "$u_computer" != "p" ]; then
        module reset > /dev/null 2>&1
        code="$?"

        # Purge if reset is not available
        if [ "$code" -ne "0" ]; then
            module purge > /dev/null 2>&1
        fi
    else
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

    # Load modules ($MODULES)
    for module_name in ${MODULES[@]}; do
        log_n " - Load $CYAN$module_name$COLOR_RESET "

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
            echo -e "[$G""SUCCESS$W]"
        else
            echo -e "[$R""FAILURE$W]"

            # Run load again to show error message
            module load "$module_name"

            return
        fi
    done

    ok "All modules have been loaded."

    return
elif [ "$1" == "format" ]; then
    if ! command -v fprettify > /dev/null 2>&1; then
        pip3 install --upgrade fprettify
    fi

    fprettify src --exclude "src/*/autogen" --recursive \
        --indent 4 --c-relations --enable-replacements --enable-decl \
        --whitespace-comma 1 --whitespace-multdiv 0 --whitespace-plusminus 1 \
        --case 1 1 1 1 --strict-indent
    ret="$?"

    if [ "$ret" != "0" ]; then
        error "failed to execute fprettify."
        error "MFC has not been fprettify'ied."

        exit 1
    fi

    ok "MFC has been fprettify'ied."

    exit 0
elif [ "$1" == "cloc" ]; then
    if ! command -v cloc > /dev/null 2>&1; then
        error "cloc (github.com/AlDanial/cloc) is not installed."

        exit 1
    fi

    cloc .          --exclude-dir=build,tests,examples,.vscode,.github \
         --fullpath --not-match-d=src/*/*/autogen

    exit $?
elif [ "$1" == "docker" ]; then
    shift;

    if ! command -v docker > /dev/null 2>&1; then
        error "$MAGENTA""Docker$COLOR_RESET is not installed."

        exit 1
    fi

    log "Running in$MAGENTA Docker$COLOR_RESET mode."
    log "  - Fetching image..."
    if ! docker pull henryleberre/mfc; then
        error "Failed to fetch$MAGENTA Docker$COLOR_RESET image from$MAGENTA Docker Hub$COLOR_RESET."

        exit 1
    fi

    log "  - Starting container..."
    docker run --interactive --tty --rm                              \
               --mount type=bind,source="$(pwd)",target=/home/me/MFC \
               henryleberre/mfc
    if (($?)); then
        error "Failed to start Docker container."

        exit 1
    fi

    exit 0
fi


mkdir -p "$(pwd)/build"


# If we downloaded our own version of CMake, let us make it accessible to
# any command we run hereafter.

if [ -f "$(pwd)/build/cmake/bin/cmake" ]; then
    export PATH="$(pwd)/build/cmake/bin:$PATH"

    log "Found$MAGENTA CMake$COLOR_RESET: $MAGENTA$(pwd)/build/cmake/bin/cmake$COLOR_RESET."
fi

bShouldInstallCMake=false
if ! command -v cmake > /dev/null 2>&1; then
    # Not installed
    bShouldInstallCMake=true

    warn "$MAGENTA""CMake$COLOR_RESET is not installed."
else
    cmake_verstr=$(cmake --version | tr ' ' '\n' | sed -n 3p)
    cmake_major=$(echo $cmake_verstr | tr '.' '\n' | sed -n 1p)
    cmake_minor=$(echo $cmake_verstr | tr '.' '\n' | sed -n 2p)
    cmake_patch=$(echo $cmake_verstr | tr '.' '\n' | sed -n 3p)
    cmake_version="$(printf %05d%05d%05d $cmake_major $cmake_minor $cmake_patch)"

    MFC_CMAKE_MIN_VERSTR=$(cat CMakeLists.txt | head -n 1 | sed 's/[^0-9,.]*//g')
    MFC_CMAKE_MIN_MAJOR=$(echo $MFC_CMAKE_MIN_VERSTR | tr '.' '\n' | head -n 1)
    MFC_CMAKE_MIN_MINOR=$(echo $MFC_CMAKE_MIN_VERSTR | tr '.' '\n' | head -n 2 | tail -n 1)
    MFC_CMAKE_MIN_PATCH=0
    MFC_CMAKE_MIN_VERSION="$(printf %05d%05d%05d $MFC_CMAKE_MIN_MAJOR $MFC_CMAKE_MIN_MINOR $MFC_CMAKE_MIN_PATCH)"

    if [ "$cmake_version" -lt "$MFC_CMAKE_MIN_VERSION" ]; then
        # Out of date
        bShouldInstallCMake=true

        log "$MAGENTA CMake$YELLOW is out of date (current: $MAGENTA$cmake_verstr$YELLOW < minimum: $MAGENTA$MFC_CMAKE_MIN_VERSTR$YELLOW).$COLOR_RESET"
    fi
fi

if [ "$bShouldInstallCMake" = true ]; then
    version="3.24.2"
    arch="$(uname -m)"

    bErrorDoSelfDownload=false
    if [[ "$OSTYPE" == "darwin"* ]]; then
        error "Cannot fetch$MAGENTA CMake$COLOR_RESET for$MAGENTA macOS$COLOR_RESET."
        bErrorDoSelfDownload=true
    fi

    if ! [[ "$arch" == "x86_64" || "$arch" == "aarch64" ]]; then
        error "Cannot fetch$MAGENTA CMake$COLOR_RESET for architecture $MAGENTA$arch$COLOR_RESET."
        bErrorDoSelfDownload=true
    fi

    if [ "$bErrorDoSelfDownload" = true ]; then
        log "Please install$MAGENTA CMake$COLOR_RESET manually:"
        log " - via$MAGENTA Aptitude$COLOR_RESET:  $ sudo apt install cmake"
        log " - via$MAGENTA Homebrew$COLOR_RESET:  $ brew install cmake"
        log " - via$MAGENTA Spack$COLOR_RESET:     $ spack install cmake"
        log " - via$MAGENTA Pacman$COLOR_RESET:    $ sudo pacman -S cmake"
        log " - via$MAGENTA Lmod$COLOR_RESET:      $ module load cmake (if available)"
        log " - via$MAGENTA cmake.org$COLOR_RESET: https://cmake.org/download/"

        exit 1
    fi

    if [ -d "$(pwd)/build/cmake" ]; then
        rm -rf "$(pwd)/build/cmake"
    fi

    mkdir -p "$(pwd)/build/cmake"

    filename="cmake-$version-linux-$arch.sh"
    repository="https://github.com/Kitware/CMake"
    url="$repository/releases/download/v$version/$filename"

    log "Downloading$MAGENTA CMake v$version$COLOR_RESET for $MAGENTA$arch$COLOR_RESET from $CYAN$repository$COLOR_RESET."

    cmake_fatal_error() {
        log   "$MAGENTA""CMake$COLOR_RESET is not discoverable or is an older release, incompatible with$MAGENTA MFC$COLOR_RESET. Please download"
        log   "or install a recent version of$MAGENTA CMake$COLOR_RESET to get past this step. If you are currently on a"
        log   "managed system like a cluster, provided there is no suitable environment module, you can"
        log   "either build it from source, or get it via$MAGENTA Spack$COLOR_RESET."
        log   "- The minimum required version is currently$MAGENTA CMake$COLOR_RESET v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
        log   "- We attempted to download$MAGENTA CMake$COLOR_RESET v$version from $url."

        exit 1
    }

    if ! command -v wget > /dev/null 2>&1; then
        error "$MAGENTA""wget$COLOR_RESET is not installed but is necessary to download$MAGENTA CMake$COLOR_RESET."

        cmake_fatal_error
    fi

    if ! wget -P "$(pwd)/build/cmake" "$url"; then
        error "Failed to download a compatible version of$MAGENTA CMake$COLOR_RESET."

        cmake_fatal_error
    fi

    log "Installing$MAGENTA CMake$COLOR_RESET into $MAGENTA$(pwd)/build/cmake$COLOR_RESET."

    if ! $SHELL "$(pwd)/build/cmake/$filename" "--skip-license" "--prefix=$(pwd)/build/cmake"; then
        error "Failed to install a compatible version of CMake."

        cmake_fatal_error
    fi

    rm "$(pwd)/build/cmake/$filename"

    ok "Found$MAGENTA CMake$COLOR_RESET: $MAGENTA$(pwd)/build/cmake/bin/cmake$COLOR_RESET."
    export PATH="$(pwd)/build/cmake/bin:$PATH"
fi


if [ -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check Python is still working within the VENV
    if ! $(pwd)/build/venv/bin/python3 --version > /dev/null 2>&1; then
        # If not, delete it and install it again
        warn "$MAGENTA""Python$COLOR_RESET is no longer working inside the Virtualenv."
        warn "Deleting the Virtualenv and starting from scratch..."

        rm -r "$(pwd)/build/venv"
    fi
fi

if ! command -v pip3 > /dev/null 2>&1 && [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check whether python3 is in the $PATH / is accessible.
    if ! command -v python3 > /dev/null 2>&1; then
        error "Couldn't find$MAGENTA Python$COLOR_RESET. Please ensure it is discoverable."

        exit 1
    fi

    # CHeck Python's version for compatibility
    if ! python3 -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"; then
        error "$(python3 --version) is incompatible.$MAGENTA Python$COLOR_RESET v$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR or higher is required."

        exit 1
    fi

    get_pip_url="https://bootstrap.pypa.io/pip/get-pip.py"

    warn "$MAGENTA""Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET is not installed."
    log  "Downloading$MAGENTA Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET from $get_pip_url..."

    if ! wget -O "$(pwd)/build/get-pip.py" "$get_pip_url"; then
        error "Couldn't download get-pip.py."

        exit 1
    fi

    # Suppress PIP version warning (out of date)
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    if ! python3 "$(pwd)/build/get-pip.py" --user; then
        error "Couldn't install$MAGENTA pip$COLOR_RESET with get-pip.py"

        exit 1
    fi

    ok "Installed pip."
fi


# Create a Python virtualenv if it hasn't already been created
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    if ! python3 -m venv "$(pwd)/build/venv"; then
        error "Failed to create a$MAGENTA Python$COLOR_RESET virtual environment. Delete the build/venv folder and try again."

        exit 1
    fi

    ok "Created a$MAGENTA Python$COLOR_RESET virtual environment (venv)."
fi


# If the user wishes to enter the python venv
# mfc.sh must be source'd for this to work
if [ "$1" == "venv" ]; then
    if [[ "$VIRTUAL_ENV" != "" ]]; then
        # Already inside the venv, exit out of it

        log " > Exiting the$MAGENTA MFC Python$COLOR_RESET virtual environment."

        deactivate
    else
        # Enter the venv
        log " > Entering the$MAGENTA MFC Python$COLOR_RESET virtual environment."

        source "$(pwd)/build/venv/bin/activate"

        log " > To exit, you can do any of the following:"
        log "    - Run 'deactivate'."
        log "    - Run '. ./mfc.sh venv'."
        log "    - Close your terminal."
    fi

    return
fi


# Activate the Python venv
source "$(pwd)/build/venv/bin/activate"
ok "(venv) Entered the$MAGENTA Python$COLOR_RESET virtual environment."


# Install Python dependencies if, either:
# - This script is running for the first time
# (or)
# - The requirements.txt file has changed
if ! cmp "$(pwd)/toolchain/requirements.txt" "$(pwd)/build/requirements.txt" > /dev/null 2>&1; then
    log "(venv) (Re)Installing mfc.sh's Python dependencies (via Pip)."

    if ! PIP_DISABLE_PIP_VERSION_CHECK=1 pip3 install -r "$(pwd)/toolchain/requirements.txt"; then
        error "(venv) Installation failed."

        log   "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
        deactivate

        exit 1
    fi

    ok "(venv) Installation succeeded."

    # Save the new/current requirements.txt
    cp "$(pwd)/toolchain/requirements.txt" "$(pwd)/build/"
fi


# Run the main.py bootstrap script
python3 "$(pwd)/toolchain/mfc.py" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
log "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
deactivate


# Exit proper exit code
exit $code
