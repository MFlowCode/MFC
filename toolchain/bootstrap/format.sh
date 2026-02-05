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
    # Custom paths provided - format all file types in those paths
    SEARCH_PATHS="${PATHS[@]}"

    # Format Fortran files (.f90, .fpp)
    if ! find $SEARCH_PATHS -type f 2>/dev/null | grep -Ev 'autogen' | grep -E '\.(f90|fpp)$' \
            | xargs --no-run-if-empty -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_file.sh; then
        error "Formatting Fortran files failed."
        exit 1
    fi

    # Format Python files
    if ! find $SEARCH_PATHS -type f 2>/dev/null | grep -E '\.(py)$' \
            | xargs --no-run-if-empty -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_python.sh; then
        error "Formatting Python files failed."
        exit 1
    fi
else
    # Default: format src/, examples/, and benchmarks/

    # Format Fortran files (.f90, .fpp) in src/
    if ! find src -type f 2>/dev/null | grep -Ev 'autogen' | grep -E '\.(f90|fpp)$' \
            | xargs --no-run-if-empty -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_file.sh; then
        error "Formatting MFC source failed."
        exit 1
    fi

    # Format Python files in examples/
    if ! find examples -type f 2>/dev/null | grep -E '\.(py)$' \
            | xargs --no-run-if-empty -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_python.sh; then
        error "Formatting MFC examples failed."
        exit 1
    fi

    # Format Python files in benchmarks/
    if ! find benchmarks -type f 2>/dev/null | grep -E '\.(py)$' \
            | xargs --no-run-if-empty -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_python.sh; then
        error "Formatting MFC benchmarks failed."
        exit 1
    fi
fi

ok "Done. MFC has been formatted."
