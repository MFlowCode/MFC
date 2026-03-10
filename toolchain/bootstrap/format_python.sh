#!/bin/bash

. toolchain/util.sh

echo "> $1"

# Auto-fix lint issues, then format
ruff check --fix "$1" > /dev/null || true
if ! ruff format "$1" 2>&1; then
    error "Failed to format $1 with ruff"
    exit 1
fi

