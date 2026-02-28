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

# Install viz optional deps only if not already present.  On air-gapped systems or
# networks where PyPI is unreachable the install may fail; in that case we skip the
# viz-specific lint and tests rather than aborting the entire precheck.
VIZ_LINT=true
if ! python3 -c "import matplotlib, dash, textual, imageio, h5py, plotext, plotly" 2>/dev/null; then
    log "(venv) Installing$MAGENTA viz$COLOR_RESET optional dependencies for linting..."
    if ! { uv pip install -q "$(pwd)/toolchain[viz]" 2>/dev/null \
           || python3 -m pip install -q "$(pwd)/toolchain[viz]" 2>/dev/null; }; then
        log "${YELLOW}Warning:${COLOR_RESET} viz optional dependencies could not be installed (no network?). Skipping viz lint/tests."
        VIZ_LINT=false
    fi
fi

log "(venv) Running$MAGENTA pylint$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""toolchain$COLOR_RESET."

# Exclude the viz subpackage from pylint when its optional deps are unavailable,
# since pylint needs to import matplotlib/dash/etc. to analyse those modules.
PYLINT_VIZ_OPT=""
[ "$VIZ_LINT" = false ] && PYLINT_VIZ_OPT="--ignore-paths=.*/mfc/viz/.*"
# shellcheck disable=SC2086
pylint -d R1722,W0718,C0301,C0116,C0115,C0114,C0410,W0622,W0640,C0103,W1309,C0411,W1514,R0401,W0511,C0321,C3001,R0801,R0911,R0912 $PYLINT_VIZ_OPT "$(pwd)/toolchain/"

log "(venv) Running$MAGENTA pylint$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""examples$COLOR_RESET."

pylint -d C0103,C0114,C0301,R0801,C0410,W0611,W1514,E0401,C0115,C0116,C0200,W1309,W0401,E0602,R1720,W0614,E1101 $(pwd)/examples/*/case.py

log "(venv) Running$MAGENTA pylint$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""benchmarks$COLOR_RESET."

pylint -d C0103,C0114,C0301,R0801,C0410,W0611,W1514,E0401,C0115,C0116,C0200,W1309,W0401,E0602,R1720,W0614,E1101 $(pwd)/benchmarks/*/case.py

# Run toolchain unit tests unless --no-test is specified
if [ "$RUN_TESTS" = true ]; then
    log "(venv) Running$MAGENTA unit tests$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""toolchain$COLOR_RESET."

    # Run tests as modules from the toolchain directory to resolve relative imports
    cd "$(pwd)/toolchain"
    python3 -m unittest mfc.params_tests.test_registry mfc.params_tests.test_definitions mfc.params_tests.test_validate mfc.params_tests.test_integration -v
    python3 -m unittest mfc.cli.test_cli -v
    cd - > /dev/null
fi

exit 0
