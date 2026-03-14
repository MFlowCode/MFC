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

# --- Phase 2: All remaining checks in parallel (read-only) ---

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

# Source lint (fast grep checks — run inline)
SOURCE_FAILED=0
SOURCE_MSGS=""
if grep -qiR '!\$acc\|!\$omp' --exclude="parallel_macros.fpp" --exclude="acc_macros.fpp" --exclude="omp_macros.fpp" --exclude="shared_parallel_macros.fpp" --exclude="syscheck.fpp" ./src/* 2>/dev/null; then
    SOURCE_MSGS+="Found raw OpenACC/OpenMP directives. Use macros instead.\n"
    SOURCE_FAILED=1
fi
if grep -qiR 'double_precision\|dsqrt\|dexp\|dlog\|dble\|dabs\|double\ precision\|real(8)\|real(4)\|dprod\|dmin\|dmax\|dfloat\|dreal\|dcos\|dsin\|dtan\|dsign\|dtanh\|dsinh\|dcosh\|d0' --exclude-dir=syscheck --exclude="*nvtx*" --exclude="*precision_select*" ./src/* 2>/dev/null; then
    SOURCE_MSGS+="Found double precision intrinsics. Use generic intrinsics.\n"
    SOURCE_FAILED=1
fi
if grep -qiR -e '\.\.\.' -e '\-\-\-' -e '===' ./src/* 2>/dev/null; then
    SOURCE_MSGS+="Found junk code patterns (..., ---, ===) in source.\n"
    SOURCE_FAILED=1
fi

# Doc reference check
DOC_FAILED=0
if ! python3 toolchain/mfc/lint_docs.py > /dev/null 2>&1; then
    DOC_FAILED=1
fi

# --- Collect results ---

FAILED=0

log "[$CYAN 1/5$COLOR_RESET] Checking$MAGENTA formatting$COLOR_RESET..."
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
log "[$CYAN 2/5$COLOR_RESET] Running$MAGENTA spell check$COLOR_RESET..."
SPELL_RC=$(cat "$TMPDIR_PC/spell_exit" 2>/dev/null || echo "1")
if [ "$SPELL_RC" = "0" ]; then
    ok "Spell check passed."
else
    error "Spell check failed. Run$MAGENTA ./mfc.sh spelling$COLOR_RESET for details."
    FAILED=1
fi

wait $PID_LINT
log "[$CYAN 3/5$COLOR_RESET] Running$MAGENTA toolchain lint$COLOR_RESET..."
LINT_RC=$(cat "$TMPDIR_PC/lint_exit" 2>/dev/null || echo "1")
if [ "$LINT_RC" = "0" ]; then
    ok "Toolchain lint passed."
else
    error "Toolchain lint failed. Run$MAGENTA ./mfc.sh lint$COLOR_RESET for details."
    FAILED=1
fi

log "[$CYAN 4/5$COLOR_RESET] Running$MAGENTA source lint$COLOR_RESET checks..."
if [ $SOURCE_FAILED -eq 0 ]; then
    ok "Source lint passed."
else
    echo -e "$SOURCE_MSGS" | while read -r msg; do
        [ -n "$msg" ] && error "$msg"
    done
    FAILED=1
fi

log "[$CYAN 5/5$COLOR_RESET] Checking$MAGENTA doc references$COLOR_RESET..."
if [ $DOC_FAILED -eq 0 ]; then
    ok "Doc references are valid."
else
    error "Doc reference check failed. Run$MAGENTA python3 toolchain/mfc/lint_docs.py$COLOR_RESET for details."
    FAILED=1
fi

echo ""

if [ $FAILED -eq 0 ]; then
    ok "All precheck tests passed! Safe to push."
    exit 0
else
    error "Some precheck tests failed. Fix issues before pushing."
    exit 1
fi
