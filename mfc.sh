#!/usr/bin/env bash


# Script Constants
MFC_DIR="./.mfc"
MFC_GET_PIP_PATH="$MFC_DIR/get-pip.py"
PYTHON_VENV_DIR="$MFC_DIR/.venv"


# Check whether this script was called from MFC's root directory.
if [ ! -f ./bootstrap/delegate.py ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi

# Make bootstrap files executable
chmod +x ./bootstrap/delegate.py
chmod +x ./bootstrap/prepare_make.sh

# Check whether python3 is in the $PATH / is accessible.
which python3 > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python3. Please ensure it is discoverable."
    exit 1
fi

# Check if Python is at least minimally functionnal.
python3 -c 'print("")' > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Python3 is present but can't execute a simple program. Please ensure that python3 is working."
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

    wget -O "$MFC_GET_PIP_PATH" https://bootstrap.pypa.io/get-pip.py
    if (($?)); then
        echo "[mfc.sh] Error: Couldn't download get-pip.py using wget to '$MFC_GET_PIP_PATH'."
        exit 1
    fi

    python3 "$MFC_GET_PIP_PATH" --user
    if (($?)); then
        echo "[mfc.sh] Error: Coudln't install pip with get-pip.py ($MFC_GET_PIP_PATH)."
        exit 1
    fi
fi

# Create a Python virtualenv if it hasn't already been created
if [ ! -d "$PYTHON_VENV_DIR" ]; then
    python3 -m venv "$PYTHON_VENV_DIR"
    if (($?)); then
        echo "[mfc.sh] Error: Failed to create a Python virtual environment."
        exit 1
    fi
fi

# Activate the Python venv
source "$PYTHON_VENV_DIR"/bin/activate

# Fetch required Python modules
declare -a REQUIRED_PYTHON_MODULES=("argparse,argparse" "dataclasses,dataclasses" "typing,typing" "yaml,pyyaml" "colorama,colorama" "fypp,fypp")

for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    import_name=$(echo $module | tr ',' '\n' | head -n 1)
    install_name=$(echo $module | tr ',' '\n' | tail -n 1)

    python3 -c "import $import_name" > /dev/null 2>&1
    if (($?)); then
        python3 -m pip install $install_name
        if (($?)); then
            echo "[mfc.sh] Error: Failed to install $import_name/$install_name through Python3's pip."
            exit 1
        fi
    fi
done

# Run the mfc.py bootstrap script
python3 ./bootstrap/delegate.py "$@"
code=$?

# Deactivate the Python virtualenv in case the user "source"'d this script
deactivate

# Exit with ./bootstrap/delegate.py's exit code
exit $code
