#!/usr/bin/env bash


# Script Constants (use MFC_xx prefix because this file can be source'd)
MFC_GET_PIP_PATH="$MFC_ENV_DIR/get-pip.py"

MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=6

# Check whether this script was called from MFC's root directory.
if [ ! -f "$(pwd)/toolchain/mfc/main.py" ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
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
        echo "[mfc.sh] Error: Docker is not installed."
        exit 1
    fi

    echo "[mfc.sh] Running in Docker mode."
    echo "  - Fetching image..."
    docker pull henryleberre/mfc
    if (($?)); then
        echo "[mfc.sh] Error: Failed to fetch Docker image from Docker Hub."
        exit 1
    fi

    echo "  - Starting container..."
    docker run --interactive --tty --rm                              \
               --mount type=bind,source="$(pwd)",target=/home/me/MFC \
               henryleberre/mfc
    if (($?)); then
        echo "[mfc.sh] Error: Failed to start Docker container."
        exit $?
    fi

    exit 0
fi

mkdir -p "$(pwd)/build"

# Check whether python3 is in the $PATH / is accessible.
command -v python3 > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python. Please ensure it is discoverable."
    exit 1
fi

# CHeck Python's version for compatibility
python3 -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"
if (($?)); then
    echo "[mfc.sh] Error: $(python3 --version) is incompatible. Python v$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR or higher is required."
    exit 1
fi


# (Re)-Install Pip via get-pip to make sure it is properly configured and working.
# Note: Some supercomputers require(d) this workaround to install and import python packages.
if [ ! -f "$(pwd)/build/get-pip.py" ]; then
    wget -O "$(pwd)/build/get-pip.py" https://bootstrap.pypa.io/pip/3.6/get-pip.py
    if (($?)); then
        echo "[mfc.sh] Error: Couldn't download get-pip.py."
        exit $?
    fi

    # Suppress PIP version warning (out of date)
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    python3 "$(pwd)/build/get-pip.py" --user
    
    if (($?)); then
        echo "[mfc.sh] Error: Coudln't install pip with get-pip.py."
        exit $?
    fi
fi


# Create a Python virtualenv if it hasn't already been created
bVenvIsNew=0
if [ ! -d "$(pwd)/build/venv" ]; then
    bVenvIsNew=1

    python3 -m venv "$(pwd)/build/venv"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to create a Python virtual environment."
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


# Upgrade Pip
if [ "$bVenvIsNew" == "1" ]; then
    python3 -m pip install --upgrade pip > /dev/null
    if (($?)); then
        echo "[mfc.sh] Error: Failed to update Pip."
        exit 1
    fi
fi


# Fetch required Python modules, including
# the MFC toolchain/kit and its dependencies,
# inside the Python venv.

python3 -c "import mfc" > /dev/null 2>&1
if (($?)); then
    python3 -m pip install -e "$(pwd)/toolchain/"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to install MFC's toolchain through Python3's pip."
        exit $?
    fi
fi

# Run the mfc.py bootstrap script
python3 "$(pwd)/toolchain/mfc/main.py" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
deactivate


# Exit proper exit code
exit $code
