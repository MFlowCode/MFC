#!/usr/bin/env bash

# If the user wishes to run the "load" script
if [ "$1" == "load" ]; then
    shift;
    source ./bootstrap/load.sh "$@"
    return
fi

# Script Constants
MFC_DIR="./build"
MFC_GET_PIP_PATH="$MFC_DIR/get-pip.py"
PYTHON_VENV_DIR="$MFC_DIR/.venv"

PYTHON_BIN="python3"
PYTHON_PIP_BIN="$PYTHON_BIN -m pip"
PYTHON_MIN_MAJOR=3
PYTHON_MIN_MINOR=6

EXEC_PATH="bootstrap/mfc.py"

# Check whether this script was called from MFC's root directory.
if [ ! -f "$EXEC_PATH" ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi

# Make bootstrap files executable
chmod +x "./$EXEC_PATH"

# Check whether python3 is in the $PATH / is accessible.
which $PYTHON_BIN > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python. Please ensure it is discoverable."
    exit 1
fi

# Check if Python is at least minimally functionnal.
$PYTHON_BIN -c 'print("")' > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Python is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi

# CHeck Python's version for compatibility with bootstrap/*.py scripts
$PYTHON_BIN -c "import sys; exit(int(not (sys.version_info[0]==$PYTHON_MIN_MAJOR and sys.version_info[1] >= $PYTHON_MIN_MINOR)))"
if (($?)); then
    echo "[mfc.sh] Error: $($PYTHON_BIN --version) is incompatible. Python v$PYTHON_MIN_MAJOR.$PYTHON_MIN_MINOR or higher is required."
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

    $PYTHON_BIN "$MFC_GET_PIP_PATH" --user
    if (($?)); then
        echo "[mfc.sh] Error: Coudln't install pip with get-pip.py ($MFC_GET_PIP_PATH)."
        exit 1
    fi
fi

# Create a Python virtualenv if it hasn't already been created
bVenvIsNew=0
if [ ! -d "$PYTHON_VENV_DIR" ]; then
    bVenvIsNew=1

    $PYTHON_BIN -m venv "$PYTHON_VENV_DIR"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to create a Python virtual environment."
        exit 1
    fi
fi

# Activate the Python venv
source "$PYTHON_VENV_DIR"/bin/activate

# Upgrade Pip
if [ "$bVenvIsNew" == "1" ]; then
    $PYTHON_PIP_BIN install --upgrade pip > /dev/null
    if (($?)); then
        echo "[mfc.sh] Error: Failed to update Pip."
        exit 1
    fi
fi

# Fetch required Python modules.
# Some modules which are now in Python's standard library
#                    are imported as backports to support Python v3.6.
declare -a REQUIRED_PYTHON_MODULES=("argparse,argparse" "dataclasses,dataclasses" "typing,typing" "yaml,pyyaml" "rich,rich" "fypp,fypp")

for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    import_name=$(echo $module | tr ',' '\n' | head -n 1)
    install_name=$(echo $module | tr ',' '\n' | tail -n 1)

    $PYTHON_BIN -c "import $import_name" > /dev/null 2>&1
    if (($?)); then
        $PYTHON_PIP_BIN install $install_name
        if (($?)); then
            echo "[mfc.sh] Error: Failed to install $import_name/$install_name through Python3's pip."
            exit 1
        fi
    fi
done

# Run the mfc.py bootstrap script
$PYTHON_BIN "./$EXEC_PATH" "$@"
code=$?

# Deactivate the Python virtualenv in case the user "source"'d this script
deactivate

# Exit proper exit code
exit $code
