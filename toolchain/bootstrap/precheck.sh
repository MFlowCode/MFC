#!/bin/bash
set -e
set -o pipefail

# Function to display help message
show_help() {
  echo "Usage: ./mfc.sh precheck [OPTIONS]"
  echo "Run the same fast checks that CI runs before expensive tests start."
  echo "Use this locally before pushing to catch issues early."
  echo ""
  echo "Options:"
  echo "  -h, --help            Display this help message and exit."
  echo "  -j, --jobs JOBS       Runs JOBS number of parallel jobs for formatting."
  echo ""
  exit 0
}

# Cross-platform hash function (macOS uses md5, Linux uses md5sum)
compute_hash() {
    if command -v md5sum > /dev/null 2>&1; then
        md5sum | cut -d' ' -f1
    elif command -v md5 > /dev/null 2>&1; then
        md5 -q
    else
        # Fallback: use cksum if neither available
        cksum | cut -d' ' -f1
    fi
}

JOBS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            if [[ -z "$2" || "$2" == -* ]]; then
                echo "Precheck: -j/--jobs requires a value."
                exit 1
            fi
            if ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "Precheck: jobs value '$2' is not a valid number."
                exit 1
            fi
            JOBS="$2"
            shift
            ;;
        -h | --help)
            show_help
            ;;
        *)
            echo "Precheck: unknown argument: $1."
            exit 1
            ;;
    esac

    shift
done

# Skip slow rendering tests (matplotlib/imageio) during local precheck.
# CI runs the full suite via ./mfc.sh lint without this variable.
export MFC_SKIP_RENDER_TESTS=1

# Detect whether any Fortran sources or CMakeLists changed — if so, run a build.
# This catches compilation errors (duplicate declarations, missing symbols, etc.)
# cheaply: Python-only changes pay zero build cost.
if git diff HEAD --name-only 2>/dev/null | grep -qE '\.(fpp|f90)$|CMakeLists\.txt'; then
    BUILD_FORTRAN=1
else
    BUILD_FORTRAN=0
fi

if [ "$BUILD_FORTRAN" = "1" ]; then
    NCHECK=9
else
    NCHECK=8
fi

log "Running$MAGENTA precheck$COLOR_RESET (same checks as CI lint-gate)..."
echo ""

# Temp files for collecting results from parallel jobs
TMPDIR_PC=$(mktemp -d)
trap "rm -rf $TMPDIR_PC" EXIT

# --- Phase 1: Format (modifies files, must run alone) ---
BEFORE_HASH=$(git diff -- '*.f90' '*.fpp' '*.py' 2>/dev/null | compute_hash)
if ! ./mfc.sh format -j "$JOBS" > /dev/null 2>&1; then
    FORMAT_OK=1
else
    AFTER_HASH=$(git diff -- '*.f90' '*.fpp' '*.py' 2>/dev/null | compute_hash)
    if [ "$BEFORE_HASH" != "$AFTER_HASH" ]; then
        FORMAT_OK=2
    else
        FORMAT_OK=0
    fi
fi

# --- Phase 2: All fast checks in parallel (read-only) ---

# Spell check
(
    if ./mfc.sh spelling > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/spell_exit"
    else
        echo "1" > "$TMPDIR_PC/spell_exit"
    fi
) &
PID_SPELL=$!

# Lint (ruff + unit tests) — safe after format since files are stable
(
    if ./mfc.sh lint > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/lint_exit"
    else
        echo "1" > "$TMPDIR_PC/lint_exit"
    fi
) &
PID_LINT=$!

# Source lint (all checks consolidated in lint_source.py)
(
    if python3 toolchain/mfc/lint_source.py > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/source_exit"
    else
        echo "1" > "$TMPDIR_PC/source_exit"
    fi
) &
PID_SOURCE=$!

# Doc reference check
DOC_FAILED=0
if ! python3 toolchain/mfc/lint_docs.py > /dev/null 2>&1; then
    DOC_FAILED=1
fi

# Parameter documentation check
(
    if python3 toolchain/mfc/lint_param_docs.py > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/param_docs_exit"
    else
        echo "1" > "$TMPDIR_PC/param_docs_exit"
    fi
) &
PID_PARAM_DOCS=$!

# Generated files check (JSON schema, completions, docs)
(
    if ./mfc.sh generate --check > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/generate_exit"
    else
        echo "1" > "$TMPDIR_PC/generate_exit"
    fi
) &
PID_GENERATE=$!

# Example case validation
(
    failed=0
    for case in examples/*/case.py; do
        [ -f "$case" ] || continue
        if ! ./mfc.sh validate "$case" > /dev/null 2>&1; then
            failed=$((failed + 1))
        fi
    done
    echo "$failed" > "$TMPDIR_PC/examples_exit"
) &
PID_EXAMPLES=$!

# --- Collect results (fast checks) ---

FAILED=0

log "[$CYAN 1/$NCHECK$COLOR_RESET] Checking$MAGENTA formatting$COLOR_RESET..."
if [ "$FORMAT_OK" = "1" ]; then
    error "Formatting check failed to run."
    FAILED=1
elif [ "$FORMAT_OK" = "2" ]; then
    error "Code was not formatted. Files have been auto-formatted; review and stage the changes."
    echo ""
    git diff --stat -- '*.f90' '*.fpp' '*.py' 2>/dev/null || true
    echo ""
    FAILED=1
else
    ok "Formatting check passed."
fi

wait $PID_SPELL
log "[$CYAN 2/$NCHECK$COLOR_RESET] Running$MAGENTA spell check$COLOR_RESET..."
SPELL_RC=$(cat "$TMPDIR_PC/spell_exit" 2>/dev/null || echo "1")
if [ "$SPELL_RC" = "0" ]; then
    ok "Spell check passed."
else
    error "Spell check failed. Run$MAGENTA ./mfc.sh spelling$COLOR_RESET for details."
    FAILED=1
fi

wait $PID_LINT
log "[$CYAN 3/$NCHECK$COLOR_RESET] Running$MAGENTA toolchain lint$COLOR_RESET..."
LINT_RC=$(cat "$TMPDIR_PC/lint_exit" 2>/dev/null || echo "1")
if [ "$LINT_RC" = "0" ]; then
    ok "Toolchain lint passed."
else
    error "Toolchain lint failed. Run$MAGENTA ./mfc.sh lint$COLOR_RESET for details."
    FAILED=1
fi

wait $PID_SOURCE
log "[$CYAN 4/$NCHECK$COLOR_RESET] Running$MAGENTA source lint$COLOR_RESET..."
SOURCE_RC=$(cat "$TMPDIR_PC/source_exit" 2>/dev/null || echo "1")
if [ "$SOURCE_RC" = "0" ]; then
    ok "Source lint passed."
else
    error "Source lint failed. Run$MAGENTA python3 toolchain/mfc/lint_source.py$COLOR_RESET for details."
    FAILED=1
fi

log "[$CYAN 5/$NCHECK$COLOR_RESET] Checking$MAGENTA doc references$COLOR_RESET..."
if [ $DOC_FAILED -eq 0 ]; then
    ok "Doc references are valid."
else
    error "Doc reference check failed. Run$MAGENTA python3 toolchain/mfc/lint_docs.py$COLOR_RESET for details."
    FAILED=1
fi

wait $PID_PARAM_DOCS
log "[$CYAN 6/$NCHECK$COLOR_RESET] Checking$MAGENTA parameter docs$COLOR_RESET..."
PARAM_DOCS_RC=$(cat "$TMPDIR_PC/param_docs_exit" 2>/dev/null || echo "1")
if [ "$PARAM_DOCS_RC" = "0" ]; then
    ok "Parameter documentation check passed."
else
    error "Parameter documentation check failed. Run$MAGENTA python3 toolchain/mfc/lint_param_docs.py$COLOR_RESET for details."
    FAILED=1
fi

wait $PID_GENERATE
log "[$CYAN 7/$NCHECK$COLOR_RESET] Checking$MAGENTA generated files$COLOR_RESET..."
GENERATE_RC=$(cat "$TMPDIR_PC/generate_exit" 2>/dev/null || echo "1")
if [ "$GENERATE_RC" = "0" ]; then
    ok "Generated files are up to date."
else
    error "Generated files are out of date. Run$MAGENTA ./mfc.sh generate$COLOR_RESET to update."
    FAILED=1
fi

wait $PID_EXAMPLES
log "[$CYAN 8/$NCHECK$COLOR_RESET] Validating$MAGENTA example cases$COLOR_RESET..."
EXAMPLES_FAILED=$(cat "$TMPDIR_PC/examples_exit" 2>/dev/null || echo "1")
if [ "$EXAMPLES_FAILED" = "0" ]; then
    ok "All example cases are valid."
else
    error "$EXAMPLES_FAILED example case(s) failed validation. Run$MAGENTA ./mfc.sh validate examples/\*/case.py$COLOR_RESET for details."
    FAILED=1
fi

# --- Phase 3: Build (only when Fortran sources changed) ---

if [ "$BUILD_FORTRAN" = "1" ]; then
    log "[$CYAN 9/9$COLOR_RESET] Building$MAGENTA (Fortran changes detected)$COLOR_RESET..."
    if ./mfc.sh build -j "$JOBS" > "$TMPDIR_PC/build_out" 2>&1; then
        ok "Build passed."
    else
        error "Build failed. Output:"
        cat "$TMPDIR_PC/build_out" | tail -20
        FAILED=1
    fi
fi

echo ""

if [ $FAILED -eq 0 ]; then
    ok "All precheck tests passed! Safe to push."
    exit 0
else
    error "Some precheck tests failed. Fix issues before pushing."
    exit 1
fi
