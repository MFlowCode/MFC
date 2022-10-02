#!/usr/bin/env bash


MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=6

# Check whether this script was called from MFC's root directory.
if [ ! -f "$(pwd)/toolchain/mfc.py" ]; then
    echo "[mfc.sh]: Error: You must call this script from within MFC's root folder."

    exit 1
fi

# Load utility script
source "$(pwd)/toolchain/util.sh"

# Handle upgrading from older MFC build systems
if [ -d "$(pwd)/bootstrap" ] || [ -d "$(pwd)/dependencies" ] || [ -f "$(pwd)/build/mfc.lock.yaml" ]; then
    echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: You are upgrading from an older version of MFC. Please remove, if applicable, the dependencies/, bootstrap/, and build/ directories before running this command again.$COLOR_RESET"

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
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Docker is not installed.$COLOR_RESET"

        exit 1
    fi

    echo "$CYAN[mfc.sh]$COLOR_RESET: Running in Docker mode."
    echo "  - Fetching image..."
    if ! docker pull henryleberre/mfc; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to fetch Docker image from Docker Hub.$COLOR_RESET"

        exit 1
    fi

    echo "  - Starting container..."
    docker run --interactive --tty --rm                              \
               --mount type=bind,source="$(pwd)",target=/home/me/MFC \
               henryleberre/mfc
    if (($?)); then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to start Docker container.$COLOR_RESET"

        exit 1
    fi

    exit 0
fi


mkdir -p "$(pwd)/build"


# If we downloaded our own version of CMake, let us make it accessible to
# any command we run hereafter.

if [ -f "$(pwd)/build/cmake/bin/cmake" ]; then
    export PATH="$(pwd)/build/cmake/bin:$PATH"

    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Found CMake: $(pwd)/build/cmake/bin/cmake."
fi

bShouldInstallCMake=false
if ! command -v cmake > /dev/null 2>&1; then
    # Not installed
    bShouldInstallCMake=true

    echo -e "$CYAN[mfc.sh]$COLOR_RESET:$MAGENTA CMake$YELLOW is not installed.$COLOR_RESET"
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

        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$MAGENTA CMake$YELLOW is out of date (current: $MAGENTA$cmake_verstr$YELLOW < minimum: $MAGENTA$MFC_CMAKE_MIN_VERSTR$YELLOW).$COLOR_RESET"
    fi
fi

if [ "$bShouldInstallCMake" = true ]; then
    if [ -d "$(pwd)/build/cmake" ]; then
        rm -rf "$(pwd)/build/cmake"
    fi

    mkdir -p "$(pwd)/build/cmake"

    version="3.24.2"
    arch="$(uname -m)"
    filename="cmake-$version-linux-$arch.sh"
    repository="https://github.com/Kitware/CMake"
    url="$repository/releases/download/v$version/$filename"

    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Downloading$MAGENTA CMake v$version$COLOR_RESET for $MAGENTA$arch$COLOR_RESET from $CYAN$repository$COLOR_RESET."

    if ! wget -P "$(pwd)/build/cmake" "$url"; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to download a compatible version of CMake."
        echo "CMake is not discoverable or is an older release, incompatible with MFC. Please download"
        echo "or install a recent version of CMake to get past this step. If you are currently on a"
        echo "managed system like a cluster, provided there is no suitable environment module, you can"
        echo "either build it from source, or get it via Spack."
        echo "- The minimum required version is currently CMake v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
        echo -e "- We attempted to download CMake v$version from $url.$COLOR_RESET"

        exit 1
    fi

    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Installing$MAGENTA CMake$COLOR_RESET into $MAGENTA$(pwd)/build/cmake$COLOR_RESET."

    if ! $SHELL "$(pwd)/build/cmake/$filename" "--skip-license" "--prefix=$(pwd)/build/cmake"; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to install a compatible version of CMake."
        echo "CMake is not discoverable or is an older release, incompatible with MFC. Please download"
        echo "or install a recent version of CMake to get past this step. If you are currently on a"
        echo "managed system like a cluster, provided there is no suitable environment module, you can"
        echo "either build it from source, or get it via Spack."
        echo "- The minimum required version is currently CMake v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
        echo -e "- We attempted to install CMake v$version from $url.$COLOR_RESET"

        exit 1
    fi

    rm "$(pwd)/build/cmake/$filename"

    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Found$MAGENTA CMake$COLOR_RESET: $MAGENTA$(pwd)/build/cmake/bin/cmake$COLOR_RESET."
    export PATH="$(pwd)/build/cmake/bin:$PATH"
fi


if [ -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check Python is still working within the VENV
    if ! $(pwd)/build/venv/bin/python3 --version > /dev/null 2>&1; then
        # If not, delete it and install it again
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$YELLOW WARNING: Python is no longer working inside the Virtualenv."
        echo -e "                  Deleting the Virtualenv and starting from scratch... $COLOR_RESET"

        rm -r "$(pwd)/build/venv"
    fi
fi

if ! command -v pip3 > /dev/null 2>&1 && [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check whether python3 is in the $PATH / is accessible.
    if ! command -v python3 > /dev/null 2>&1; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Couldn't find Python. Please ensure it is discoverable.$COLOR_RESET"

        exit 1
    fi

    # CHeck Python's version for compatibility
    if ! python3 -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: $(python3 --version) is incompatible. Python v$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR or higher is required.$COLOR_RESET"

        exit 1
    fi

    get_pip_url="https://bootstrap.pypa.io/pip/get-pip.py"

    echo -e "$CYAN[mfc.sh]$YELLOW Python's PIP is not installed.$COLOR_RESET"
    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Downloading Python's PIP from $get_pip_url..."

    if ! wget -O "$(pwd)/build/get-pip.py" "$get_pip_url"; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Couldn't download get-pip.py.$COLOR_RESET"

        exit 1
    fi

    # Suppress PIP version warning (out of date)
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    if ! python3 "$(pwd)/build/get-pip.py" --user; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Couldn't install pip with get-pip.py.$COLOR_RESET"

        exit 1
    fi
fi


# Create a Python virtualenv if it hasn't already been created
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    echo -e "$CYAN[mfc.sh]$COLOR_RESET: Creating a Python virtual environment (venv)..."

    if ! python3 -m venv "$(pwd)/build/venv"; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to create a Python virtual environment. Delete the build/venv folder and try again.$COLOR_RESET"

        exit 1
    fi
fi


# If the user wishes to enter the python venv
# mfc.sh must be source'd for this to work
if [ "$1" == "venv" ]; then
    if [[ "$VIRTUAL_ENV" != "" ]]; then
        # Already inside the venv, exit out of it

        echo " > Exiting the MFC Python virtual environment."

        deactivate
    else
        # Enter the venv
        echo " > Entering the MFC Python virtual environment."

        source "$(pwd)/build/venv/bin/activate"

        echo " > To exit, you can do any of the following:"
        echo "    - Run 'deactivate'."
        echo "    - Run '. ./mfc.sh venv'."
        echo "    - Close your terminal."
    fi

    return
fi


# Activate the Python venv
echo -e "$CYAN[mfc.sh]$COLOR_RESET: Entering the Python virtual environment (venv)."
source "$(pwd)/build/venv/bin/activate"


# Fetch required Python modules.
# Some modules which are now in Python's standard library
#                    are imported as backports to support Python v3.6.
declare -a REQUIRED_PYTHON_MODULES=("wheel,wheel" "argparse,argparse" "dataclasses,dataclasses" "typing,typing" "yaml,pyyaml" "rich,rich" "fypp,fypp")

for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    import_name=$(echo $module | tr ',' '\n' | head -n 1)
    install_name=$(echo $module | tr ',' '\n' | tail -n 1)

    if ! python3 -c "import $import_name" > /dev/null 2>&1; then
        echo -e "$CYAN[mfc.sh]$COLOR_RESET: Installing Python PIP package $MAGENTA$install_name$COLOR_RESET (into venv)."

        if ! pip3 install "$install_name"; then
            echo -e "$CYAN[mfc.sh]$COLOR_RESET:$RED Error: Failed to install $import_name/$install_name through Python3's pip.$COLOR_RESET"

            exit 1
        fi
    fi
done


# Run the mfc.py bootstrap script
python3 "$(pwd)/toolchain/mfc.py" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
echo -e "$CYAN[mfc.sh]$COLOR_RESET: Exiting the Python virtual environment."
deactivate


# Exit proper exit code
exit $code
