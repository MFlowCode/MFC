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

# Format Fortran files (.f90, .fpp)
FORTRAN_FILES=$(find $FORTRAN_DIRS -type f 2>/dev/null | grep -Ev 'autogen' | grep -E '\.(f90|fpp)$' || true)
if [[ -n "$FORTRAN_FILES" ]]; then
    FPRETTIFY_OPTS="--silent --indent 4 --c-relations --enable-replacements --enable-decl --whitespace-comma 1 --whitespace-multdiv 0 --whitespace-plusminus 1 --case 1 1 1 1 --strict-indent --line-length 1000"

    # Skip files unchanged since last format (hash cache in build/.cache/format/)
    CACHE_DIR="build/.cache/format"
    mkdir -p "$CACHE_DIR"
    DIRTY_FILES=""
    for f in $FORTRAN_FILES; do
        cache_key=$(echo "$f" | tr '/' '_')
        current_hash=$(md5sum "$f" | cut -d' ' -f1)
        cached_hash=$(cat "$CACHE_DIR/$cache_key" 2>/dev/null || true)
        if [[ "$current_hash" != "$cached_hash" ]]; then
            DIRTY_FILES+="$f"$'\n'
        fi
    done
    DIRTY_FILES=$(echo "$DIRTY_FILES" | sed '/^$/d')

    if [[ -n "$DIRTY_FILES" ]]; then
        for niter in 1 2 3 4; do
            old_hash=$(echo "$DIRTY_FILES" | xargs cat | md5sum)

            # Run indenter on dirty files in one process
            if ! echo "$DIRTY_FILES" | xargs python3 toolchain/indenter.py; then
                error "Formatting Fortran files failed: indenter.py."
                exit 1
            fi

            # Run fprettify in parallel (one process per file)
            if ! echo "$DIRTY_FILES" | xargs -P ${JOBS:-1} -L 1 fprettify $FPRETTIFY_OPTS; then
                error "Formatting Fortran files failed: fprettify."
                exit 1
            fi

            new_hash=$(echo "$DIRTY_FILES" | xargs cat | md5sum)
            if [[ "$old_hash" == "$new_hash" ]]; then
                break
            fi
            if [[ "$niter" -eq 4 ]]; then
                error "Formatting Fortran files failed: no steady-state after $niter iterations."
                exit 1
            fi
        done

        # Update hash cache for formatted files
        for f in $DIRTY_FILES; do
            cache_key=$(echo "$f" | tr '/' '_')
            md5sum "$f" | cut -d' ' -f1 > "$CACHE_DIR/$cache_key"
        done

        echo "$DIRTY_FILES" | while read -r f; do echo "> $f"; done
    fi
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
