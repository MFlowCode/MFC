#!/bin/bash
set -e
set -o pipefail

# Parse arguments
RUN_TESTS=true
for arg in "$@"; do
    case $arg in
        --no-test)
            RUN_TESTS=false
            ;;
    esac
done

log "(venv) Auto-fixing safe lint issues with$MAGENTA ruff$COLOR_RESET..."

ruff check --fix-only toolchain/ examples/*/case.py benchmarks/*/case.py

log "(venv) Running$MAGENTA ruff$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""toolchain$COLOR_RESET."

ruff check toolchain/

log "(venv) Running$MAGENTA ruff$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""examples$COLOR_RESET."

ruff check examples/*/case.py

log "(venv) Running$MAGENTA ruff$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""benchmarks$COLOR_RESET."

ruff check benchmarks/*/case.py

# Run toolchain unit tests unless --no-test is specified
if [ "$RUN_TESTS" = true ]; then
    log "(venv) Running$MAGENTA unit tests$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""toolchain$COLOR_RESET."

    cd toolchain
    python3 -m pytest . -v
    cd - > /dev/null
fi

exit 0
