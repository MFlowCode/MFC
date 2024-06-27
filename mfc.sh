#!/bin/bash

# Check whether this script was called from MFC's root directory.
if [ ! -f "$(pwd)/toolchain/util.sh" ]; then
    echo "mfc: ERROR > You must call this script from within MFC's root folder."

    exit 1
fi

# Load utility script
. "$(pwd)/toolchain/util.sh"

# Handle upgrading from older MFC build systems
if [ -d "$(pwd)/bootstrap" ] || [ -d "$(pwd)/dependencies" ] || [ -f "$(pwd)/build/mfc.lock.yaml" ]; then
    error "Please remove, if applicable, the following directories:"
    error "* dependencies/"
    error "* bootstrap/"
    error "* build/"
    error "These are incompatible with modern versions of MFC."

    exit 1
fi

# If the user wishes to run the "load" script
if [ "$1" '==' 'load' ]; then
    shift; . "$(pwd)/toolchain/bootstrap/modules.sh" $@; return
elif [ "$1" '==' "lint" ]; then
    . "$(pwd)/toolchain/bootstrap/python.sh"

    shift; . "$(pwd)/toolchain/bootstrap/lint.sh"    $@; exit 0
elif [ "$1" '==' "format" ]; then
    . "$(pwd)/toolchain/bootstrap/python.sh"

    shift; . "$(pwd)/toolchain/bootstrap/format.sh"  $@; exit 0
elif [ "$1" '==' "docker" ]; then
    shift; . "$(pwd)/toolchain/bootstrap/docker.sh"  $@; exit 0
elif [ "$1" '==' "venv" ]; then
    shift; . "$(pwd)/toolchain/bootstrap/python.sh"  $@; return
fi

mkdir -p "$(pwd)/build"

. "$(pwd)/toolchain/bootstrap/cmake.sh"
. "$(pwd)/toolchain/bootstrap/python.sh"

echo

# Run the main.py bootstrap script
python3 "$(pwd)/toolchain/main.py" "$@"
code=$?

echo

if [ $code -ne 0 ]; then
    error "mfc.py finished with a $code exit code."
fi

# Deactivate the Python virtualenv in case the user "source"'d this script
log "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
deactivate

# Exit with proper exit code
exit $code
