#!/bin/bash

. toolchain/util.sh

echo "> $1"

# Use autopep8 for Python formatting
# (black has issues with Python 3.12.5 which is common on HPC systems)
if ! autopep8 --in-place --max-line-length 200 "$1" 2>&1; then
    error "Failed to format $1 with autopep8"
    exit 1
fi

