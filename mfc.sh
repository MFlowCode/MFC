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
    source "$(pwd)/toolchain/scripts/load.sh" "$@"

    return
elif [ "$1" == "format" ]; then
    "$SHELL" "$(pwd)/toolchain/scripts/format.sh"

    exit $?
fi

if [ "$1" == "docker" ]; then
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

    echo "  - Starting container..."
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

    warning "$MAGENTA""CMake$COLOR_RESET is not installed."
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

    warning "$MAGENTA""Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET is not installed."
    log     "Downloading$MAGENTA Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET from $get_pip_url..."

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
ok "Entered the$MAGENTA Python$COLOR_RESET virtual environment (venv)."


# Fetch required Python modules.
# Some modules which are now in Python's standard library
#                    are imported as backports to support Python v3.6.
declare -a REQUIRED_PYTHON_MODULES=("wheel,wheel" "argparse,argparse" "dataclasses,dataclasses" "typing,typing" "yaml,pyyaml" "rich,rich" "fypp,fypp")

for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    import_name=$(echo $module | tr ',' '\n' | head -n 1)
    install_name=$(echo $module | tr ',' '\n' | tail -n 1)

    if ! python3 -c "import $import_name" > /dev/null 2>&1; then
        if ! PIP_DISABLE_PIP_VERSION_CHECK=1 pip3 install "$install_name"; then
            error "Failed to install $import_name/$install_name through Python3's pip."
        
            exit 1
        fi

        ok "Installed$MAGENTA Python PIP$COLOR_RESET package $MAGENTA$install_name$COLOR_RESET (into venv)."
    fi
done


# Run the mfc.py bootstrap script
python3 "$(pwd)/toolchain/mfc.py" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
log "Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
deactivate


# Exit proper exit code
exit $code
