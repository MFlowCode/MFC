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
# ffmt prints unmatched-structure warnings to stderr yet exits 0: show stderr and treat them as fatal.
FFMT_ERR="$(mktemp)"
trap 'rm -f "$FFMT_ERR"' EXIT
ffmt -j ${JOBS:-1} $FORTRAN_DIRS 2>"$FFMT_ERR"
FFMT_RC=$?
[ -s "$FFMT_ERR" ] && cat "$FFMT_ERR" >&2
if [ $FFMT_RC -ne 0 ]; then
    error "Formatting Fortran files failed: ffmt."
    exit 1
fi
# match the construct name, not the bare word, which a diagnostic or path could carry
if grep -q "unmatched Fortran" "$FFMT_ERR"; then
    error "ffmt reported unmatched Fortran block structure: the source will not compile."
    exit 1
fi

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
