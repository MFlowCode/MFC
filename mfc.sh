#!/usr/bin/env bash


# Script Constants (use MFC_xx prefix because this file can be source'd)
MFC_DIR="$(pwd)/build"


MFC_TOOLCHAIN_DIR="$(pwd)/toolchain"
MFC_EXEC_PATH="$MFC_TOOLCHAIN_DIR/mfc/main.py"

MFC_PYTHON_BIN="python3"
MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=6
MFC_PYTHON_PIP_BIN="$MFC_PYTHON_BIN -m pip"
MFC_PYTHON_VENV_DIR="$MFC_DIR/environment/venv"
MFC_GET_PIP_PATH="$MFC_DIR/environment/get-pip.py"

# Check whether this script was called from MFC's root directory.
if [ ! -f "$MFC_EXEC_PATH" ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi


# If the user wishes to run the "load" script
if [ "$1" == "load" ]; then
    shift;
    source "$MFC_TOOLCHAIN_DIR/load.sh" "$@"
    return
fi


# Create main subdirectories inside $MFC_DIR
mkdir -p "$MFC_DIR/mfc"
mkdir -p "$MFC_DIR/environment"
mkdir -p "$MFC_DIR/dependencies"


# Make bootstrap files executable
chmod +x "$MFC_EXEC_PATH"


# Check whether python3 is in the $PATH / is accessible.
which $MFC_PYTHON_BIN > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python. Please ensure it is discoverable."
    exit 1
fi


# Check if Python is at least minimally functionnal.
$MFC_PYTHON_BIN -c 'print("")' > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Python is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi


# CHeck Python's version for compatibility with bootstrap/*.py scripts
$MFC_PYTHON_BIN -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"
if (($?)); then
    echo "[mfc.sh] Error: $($MFC_PYTHON_BIN --version) is incompatible. Python v$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR or higher is required."
    exit 1
fi


# (Re)-Install Pip via get-pip to make sure it is properly configured and working.
# Note: Some supercomputers require(d) this workaround to install and import python packages.
if [ ! -f "$MFC_GET_PIP_PATH" ]; then
    mkdir -p "$MFC_DIR"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to create directory '$MFC_DIR'."
        exit 1
    fi

    wget -O "$MFC_GET_PIP_PATH" https://bootstrap.pypa.io/pip/3.6/get-pip.py
    if (($?)); then
        echo "[mfc.sh] Error: Couldn't download get-pip.py using wget to '$MFC_GET_PIP_PATH'."
        exit 1
    fi

    # Suppress PIP version warning (out of date)
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    $MFC_PYTHON_PIP_BIN config set global.disable-pip-version-check true || true
    $MFC_PYTHON_BIN "$MFC_GET_PIP_PATH" --user
    
    if (($?)); then
        echo "[mfc.sh] Error: Coudln't install pip with get-pip.py ($MFC_GET_PIP_PATH)."
        exit 1
    fi
fi


# Create a Python virtualenv if it hasn't already been created
bVenvIsNew=0
if [ ! -d "$MFC_PYTHON_VENV_DIR" ]; then
    bVenvIsNew=1

    $MFC_PYTHON_BIN -m venv "$MFC_PYTHON_VENV_DIR"
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
        
        source "$MFC_PYTHON_VENV_DIR/bin/activate"
        
        echo " > To exit, you can do any of the following:"
        echo "    - Run 'deactivate'."
        echo "    - Run '. ./mfc.sh venv'."
        echo "    - Close your terminal."
    fi
    
    return
fi


# Activate the Python venv
source "$MFC_PYTHON_VENV_DIR/bin/activate"


# Upgrade Pip
if [ "$bVenvIsNew" == "1" ]; then
    $MFC_PYTHON_PIP_BIN install --upgrade pip > /dev/null
    if (($?)); then
        echo "[mfc.sh] Error: Failed to update Pip."
        exit 1
    fi
fi


# Fetch required Python modules, including
# the MFC toolchain/kit and its dependencies,
# inside the Python venv.

$MFC_PYTHON_BIN -c "import mfc" > /dev/null 2>&1
if (($?)); then
    $MFC_PYTHON_PIP_BIN install -e "$MFC_TOOLCHAIN_DIR/"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to install MFC's toolchain through Python3's pip."
        exit 1
    fi
fi


# Run the mfc.py bootstrap script
$MFC_PYTHON_BIN "$MFC_EXEC_PATH" "$@"
code=$?


# Deactivate the Python virtualenv in case the user "source"'d this script
deactivate


# Exit proper exit code
exit $code
