#!/bin/bash

MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=11
MFC_PYTHON_MIN_STR="$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR"

is_python_compatible() {
    if ! ${1:-python3} -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"; then
        return 1
    fi

    return 0
}

assert_python_compatible() {
    if ! is_python_compatible $1; then
        error "$MAGENTA$(${1:-python3} --version)$COLOR_RESET (${1:-python3}) is out of date. Required >= $MAGENTA$MFC_PYTHON_MIN_STR$COLOR_RESET."

        exit 1
    fi
}

if [ -f "$(pwd)/build/venv/bin/activate" ]; then
    if ! is_python_compatible "$(pwd)/build/venv/bin/python3"; then
        warn "$MAGENTA""Python$COLOR_RESET is outdated inside the Virtualenv."
        warn "Deleting the Virtualenv and starting from scratch..."

        rm -r "$(pwd)/build/venv"
    fi
fi

# Only bootstrap pip if we don't already have a venv
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check whether python3 is in the $PATH / is accessible.
    if ! command -v python3 > /dev/null 2>&1; then
        error "Couldn't find$MAGENTA Python$COLOR_RESET. Please ensure it is discoverable."

        exit 1
    fi

    assert_python_compatible

    # Check if pip is already available as a Python module
    # This works on both laptops and HPC systems with module-loaded Python
    if ! python3 -c "import pip" > /dev/null 2>&1; then
        warn "$MAGENTA""Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET is not installed."
        
        # Try ensurepip first (standard library, safe)
        log "Attempting to install pip via ensurepip..."
        if python3 -m ensurepip --upgrade 2>/dev/null; then
            ok "Installed pip via ensurepip."
        else
            # Fall back to get-pip.py only if ensurepip fails
            get_pip_url="https://bootstrap.pypa.io/pip/get-pip.py"
            log "Downloading$MAGENTA Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET from $get_pip_url..."

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

            ok "Installed pip via get-pip.py."
            
            # Ensure user-site bin directory is on PATH for this session
            user_base_bin="$(python3 -m site --user-base)/bin"
            if [ -d "$user_base_bin" ]; then
                export PATH="$user_base_bin:$PATH"
            fi
        fi
    fi
fi


# Create a Python virtualenv if it hasn't already been created
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    assert_python_compatible

    if ! python3 -m venv "$(pwd)/build/venv"; then
        error "Failed to create a$MAGENTA Python$COLOR_RESET virtual environment. Delete the build/venv folder and try again."

        exit 1
    fi

    ok "Created a$MAGENTA Python$COLOR_RESET virtual environment (venv)."

    rm "$(pwd)/build/pyproject.toml" > /dev/null 2>&1 || true
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
ok "(venv) Entered the $MAGENTA$(python3 --version)$COLOR_RESET virtual environment (>= $MAGENTA$MFC_PYTHON_MIN_STR$COLOR_RESET)."


# Install Python dependencies if, either:
# - This script is running for the first time
# (or)
# - The pyproject.toml file has changed
if ! cmp "$(pwd)/toolchain/pyproject.toml" "$(pwd)/build/pyproject.toml" > /dev/null 2>&1; then
    log "(venv) (Re)Installing mfc.sh's Python dependencies (via Pip)."

    next_arg=0
    nthreads=1
    for arg in "$@"; do
        if [ "$arg" == "-j" ] || [ "$arg" == "--jobs" ]; then
            next_arg=1
            continue
        fi
        if [ "$next_arg" == 1 ]; then
            next_arg=0
            nthreads=$arg
            continue
        fi
    done

    if ! PIP_DISABLE_PIP_VERSION_CHECK=1 MAKEFLAGS=$nthreads pip3 install "$(pwd)/toolchain"; then
        error "(venv) Installation failed."

        log   "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
        deactivate

        exit 1
    fi

    ok "(venv) Installation succeeded."

    # Save the new/current pyproject.toml
    cp "$(pwd)/toolchain/pyproject.toml" "$(pwd)/build/"
fi
