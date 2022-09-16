#!/usr/bin/env bash


MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=6

MFC_CMAKE_MIN_MAJOR=3
MFC_CMAKE_MIN_MINOR=18
MFC_CMAKE_MIN_PATCH=0
MFC_CMAKE_MIN_VERSION="$(printf %05d%05d%05d $MFC_CMAKE_MIN_MAJOR $MFC_CMAKE_MIN_MINOR $MFC_CMAKE_MIN_PATCH)"

# Check whether this script was called from MFC's root directory.
if [ ! -f "$(pwd)/toolchain/main.py" ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi

# Load utility script
source "$(pwd)/misc/util.sh"

# Handle upgrading from older MFC build systems
if [ -d "$(pwd)/bootstrap" ] || [ -d "$(pwd)/dependencies" ] || [ -f "$(pwd)/build/mfc.lock.yaml" ]; then
    echo -en "$RED"
    echo "[mfc.sh] Error: You are upgrading from an older version of MFC. Please remove, if applicable, the dependencies/, bootstrap/, and build/ directories before running this command again."
    echo -en "$COLOR_RESET"

    exit 1
fi

# If the user wishes to run the "load" script
if [ "$1" == "load" ]; then
    shift;
    source "$(pwd)/toolchain/load.sh" "$@"
    return
fi

if [ "$1" == "docker" ]; then
    shift;

    if ! command -v docker > /dev/null 2>&1; then
        echo -en "$RED"
        echo "[mfc.sh] Error: Docker is not installed."
        echo -en "$COLOR_RESET"

        exit 1
    fi

    echo "[mfc.sh] Running in Docker mode."
    echo "  - Fetching image..."
    if ! docker pull henryleberre/mfc; then
        echo -en "$RED"
        echo "[mfc.sh] Error: Failed to fetch Docker image from Docker Hub."
        echo -en "$COLOR_RESET"

        exit 1
    fi

    echo "  - Starting container..."
    docker run --interactive --tty --rm                              \
               --mount type=bind,source="$(pwd)",target=/home/me/MFC \
               henryleberre/mfc
    if (($?)); then
        echo -en "$RED"
        echo "[mfc.sh] Error: Failed to start Docker container."
        echo -en "$COLOR_RESET"

        exit 1
    fi

    exit 0
fi

mkdir -p "$(pwd)/build"


# CMake
cmake_verstr=$(cmake --version | tr ' ' '\n' | sed -n 3p)
cmake_major=$(echo $cmake_verstr | tr '.' '\n' | sed -n 1p)
cmake_minor=$(echo $cmake_verstr | tr '.' '\n' | sed -n 2p)
cmake_patch=$(echo $cmake_verstr | tr '.' '\n' | sed -n 3p)
cmake_version="$(printf %05d%05d%05d $cmake_major $cmake_minor $cmake_patch)"

if [ ! -d "$(pwd)/build/cmake/bin" ]; then
    if ! command -v cmake > /dev/null 2>&1 || [ "$cmake_version" -lt "$MFC_CMAKE_MIN_VERSION" ]; then
        # Cmake is either not installed or too old for our use

        mkdir -p "$(pwd)/build/cmake"

        version="3.24.2"
        filename="cmake-$version-linux-$(uname -m).sh"
        url="https://github.com/Kitware/CMake/releases/download/v$version/$filename"

        if ! wget -P "$(pwd)/build/cmake" "$url"; then
            echo -en "$RED"
            echo "[mfc.sh] Error: Failed to download a compatible version of CMake."
            echo "CMake is not discoverable or is an older release, incompatible with MFC. Please download"
            echo "or install a recent version of CMake to get past this step. If you are currently on a"
            echo "managed system like a cluster, provided there is no suitable environment module, you can"
            echo "either build it from source, or get it via Spack."
            echo "- The minimum required version is currently CMake v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
            echo "- We attempted to download CMake v$version from $url."
            echo -en "$COLOR_RESET"

            exit 1
        fi

        if ! $SHELL "$(pwd)/build/cmake/$filename" "--skip-license" "--prefix=$(pwd)/build/cmake"; then
            echo -en "$RED"
            echo "[mfc.sh] Error: Failed to install a compatible version of CMake."
            echo "CMake is not discoverable or is an older release, incompatible with MFC. Please download"
            echo "or install a recent version of CMake to get past this step. If you are currently on a"
            echo "managed system like a cluster, provided there is no suitable environment module, you can"
            echo "either build it from source, or get it via Spack."
            echo "- The minimum required version is currently CMake v$MFC_CMAKE_MIN_MAJOR.$MFC_CMAKE_MIN_MINOR.$MFC_CMAKE_MIN_PATCH."
            echo "- We attempted to install CMake v$version from $url."
            echo -en "$COLOR_RESET"

            exit 1
        fi

        rm "$(pwd)/build/cmake/$filename"
    fi
fi

# If we downloaded our own version of CMake, let us make it accessible to
# any command we run hereafter.

if [ -d "$(pwd)/build/cmake" ]; then
    export PATH="$(pwd)/build/cmake/bin:$PATH"
    export MANPATH="$(pwd)/build/cmake/man:$MANPATH"
fi

# Check whether python3 is in the $PATH / is accessible.
if ! command -v python3 > /dev/null 2>&1; then
    echo -en "$RED"
    echo "[mfc.sh] Error: Couldn't find Python. Please ensure it is discoverable."
    echo -en "$COLOR_RESET"

    exit 1
fi

# CHeck Python's version for compatibility
if ! python3 -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"; then
    echo -en "$RED"
    echo "[mfc.sh] Error: $(python3 --version) is incompatible. Python v$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR or higher is required."
    echo -en "$COLOR_RESET"

    exit 1
fi

if [ -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check Python is still working within the VENV
    if ! $(pwd)/build/venv/bin/python3 --version > /dev/null 2>&1; then
        # If not, delete it and install it again
        echo -en "$YELLOW"
        echo "[mfc.sh]: WARNING: Python is no longer working inside the Virtualenv."
        echo "                   Deleting the Virtualenv and starting from scratch..."
        echo -en "$COLOR_RESET"

        rm -r "$(pwd)/build/venv"
    fi
fi

if ! command -v pip3 > /dev/null 2>&1 && [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    if ! wget -O "$(pwd)/build/get-pip.py" https://bootstrap.pypa.io/pip/get-pip.py; then
        echo -en "$RED"
        echo "[mfc.sh] Error: Couldn't download get-pip.py."
        echo -en "$COLOR_RESET"

        exit 1
    fi

    # Suppress PIP version warning (out of date)
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    if ! python3 "$(pwd)/build/get-pip.py" --user; then
        echo -en "$RED"
        echo "[mfc.sh] Error: Couldn't install pip with get-pip.py."
        echo -en "$COLOR_RESET"

        exit 1
    fi
fi


# Create a Python virtualenv if it hasn't already been created
if [ ! -d "$(pwd)/build/venv" ]; then
    if ! python3 -m venv "$(pwd)/build/venv"; then
        echo -en "$RED"
        echo "[mfc.sh] Error: Failed to create a Python virtual environment. Delete the build/venv folder and try again."
        echo -en "$COLOR_RESET"

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
source "$(pwd)/build/venv/bin/activate"


# Fetch required Python modules.
# Some modules which are now in Python's standard library
#                    are imported as backports to support Python v3.6.
declare -a REQUIRED_PYTHON_MODULES=("wheel,wheel" "argparse,argparse" "dataclasses,dataclasses" "typing,typing" "yaml,pyyaml" "rich,rich" "fypp,fypp")

for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    import_name=$(echo $module | tr ',' '\n' | head -n 1)
    install_name=$(echo $module | tr ',' '\n' | tail -n 1)

    if ! python3 -c "import $import_name" > /dev/null 2>&1; then
        if ! pip3 install "$install_name"; then
            echo -en "$RED"
            echo "[mfc.sh] Error: Failed to install $import_name/$install_name through Python3's pip."
            echo -en "$COLOR_RESET"

            exit $?
        fi
    fi
done


# Run the mfc.py bootstrap script
python3 "$(pwd)/toolchain/main.py" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
deactivate


# Exit proper exit code
exit $code
