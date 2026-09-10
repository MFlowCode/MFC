#!/bin/bash

# Function to display help message
show_help() {
  echo "Usage: $(basename "$0") [OPTIONS] [PATHS...]"
  echo "Format code in the MFC repository."
  echo ""
  echo "Options:"
  echo "  -h, --help            Display this help message and exit."
  echo "  -j, --jobs JOBS       Runs JOBS number of parallel jobs."
  echo ""
  echo "Arguments:"
  echo "  PATHS                 Paths to format. Defaults to src/, examples/, benchmarks/"
  echo ""
  echo "Examples:"
  echo "  ./mfc.sh format                          Format all code"
  echo "  ./mfc.sh format -j 8                     Format with 8 parallel jobs"
  echo "  ./mfc.sh format examples/1D_sodshocktube Format a specific directory"
  echo ""
  exit 0
}

PATHS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            JOBS="$2"
            shift
            ;;
        -h | --help)
            show_help
            ;;
        *)
            # Collect positional arguments as paths
            PATHS+=("$1")
            ;;
    esac

    shift
done

log "Formatting MFC:"

if [[ ${#PATHS[@]} -gt 0 ]]; then
    SEARCH_PATHS="${PATHS[@]}"
    FORTRAN_DIRS="$SEARCH_PATHS"
    PYTHON_DIRS="$SEARCH_PATHS"
else
    FORTRAN_DIRS="src"
    PYTHON_DIRS="toolchain/ examples/ benchmarks/"
fi

# Format Fortran files with ffmt (single-pass, idempotent).
# ffmt reports structure it cannot match (e.g. "unmatched FortranBlockClose") as a WARNING on stderr
# and still exits 0. Discarding stderr therefore hid malformed Fortran -- an unbalanced end do /
# end block formatted and linted clean, and the failure only surfaced at the compiler. Capture
# stderr, always show it, and treat the structural warnings as fatal: they mean the source will not
# compile. A clean tree emits none of these, so this cannot fire on well-formed code.
FFMT_ERR="$(mktemp)"
ffmt -j ${JOBS:-1} $FORTRAN_DIRS 2>"$FFMT_ERR"
FFMT_RC=$?
[ -s "$FFMT_ERR" ] && cat "$FFMT_ERR" >&2
if [ $FFMT_RC -ne 0 ]; then
    rm -f "$FFMT_ERR"
    error "Formatting Fortran files failed: ffmt."
    exit 1
fi
if grep -q "unmatched" "$FFMT_ERR"; then
    rm -f "$FFMT_ERR"
    error "ffmt reported unmatched Fortran block structure: the source will not compile."
    exit 1
fi
rm -f "$FFMT_ERR"

# Apply safe auto-fixes (import sorting, etc.) before formatting.
# --fix-only exits 0 even when unfixable violations remain — those are
# caught later by `ruff check` in lint.sh.  This only errors if ruff
# itself fails to run.
if ! ruff check --fix-only $PYTHON_DIRS; then
    error "ruff failed to run. Check your ruff installation."
    exit 1
fi
if ! ruff format $PYTHON_DIRS; then
    error "Formatting Python files failed."
    exit 1
fi

ok "Done. MFC has been formatted."
