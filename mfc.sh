#!/usr/bin/env bash

# Check whether this script was called from MFC's root directory.
if [ ! -f ./bootstrap/mfc.py ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi

# Check whether python3 is in the $PATH / is accessible.
which python3 > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python3. Please ensure it is discoverable."
    exit 1
fi

# Check if Python is at least minimally functionnal.
python3 -c 'print("")' > /dev/null
if (($?)); then
    echo "[mfc.sh] Error: Python3 is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi

MFC_DIR="./.mfc"
MFC_GET_PIP_PATH="$MFC_DIR/get-pip.py"

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

# Optional: Use a Python virtual environment (venv)
#
#python3 -m venv ./venv
#if (($?)); then
#    echo "[mfc.sh] Error: Failed to create a Python virtual environment."
#    exit 1
#fi
#
#source ./venv/bin/activate
#if (($?)); then
#    echo "[mfc.sh] Error: Faild to activate the Python virtual environment."
#    exit 1
#fi

# Install PyYAML if it isn't installed
python3 -m pip show -q pyyaml
if (($?)); then
    python3 -m pip install pyyaml
    if (($?)); then
        echo "[mfc.sh] Error: Failed to install PyYAML through Python3's pip."
        exit 1
    fi
fi

# Install Colorama if it isn't installed
python3 -m pip show -q colorama
if (($?)); then
    python3 -m pip install colorama
    if (($?)); then
        echo "[mfc.sh] Error: Failed to install Colorama through Python3's pip."
        exit 1
    fi
fi

# Run the mfc.py bootstrap script
cd bootstrap
python3 ./mfc.py $@
code=$?

cd ..

exit $code

