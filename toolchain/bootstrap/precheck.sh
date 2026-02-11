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

FAILED=0

log "Running$MAGENTA precheck$COLOR_RESET (same checks as CI lint-gate)..."
echo ""

# 1. Check formatting
log "[$CYAN 1/4$COLOR_RESET] Checking$MAGENTA formatting$COLOR_RESET..."
# Capture state before formatting
BEFORE_HASH=$(git diff -- '*.f90' '*.fpp' '*.py' 2>/dev/null | compute_hash)
if ! ./mfc.sh format -j "$JOBS" > /dev/null 2>&1; then
    error "Formatting check failed to run."
    FAILED=1
else
    # Check if formatting changed any Fortran/Python files
    AFTER_HASH=$(git diff -- '*.f90' '*.fpp' '*.py' 2>/dev/null | compute_hash)
    if [ "$BEFORE_HASH" != "$AFTER_HASH" ]; then
        error "Code was not formatted. Files have been auto-formatted; review and stage the changes."
        echo ""
        git diff --stat -- '*.f90' '*.fpp' '*.py' 2>/dev/null || true
        echo ""
        FAILED=1
    else
        ok "Formatting check passed."
    fi
fi

# 2. Spell check
log "[$CYAN 2/4$COLOR_RESET] Running$MAGENTA spell check$COLOR_RESET..."
if ./mfc.sh spelling > /dev/null 2>&1; then
    ok "Spell check passed."
else
    error "Spell check failed. Run$MAGENTA ./mfc.sh spelling$COLOR_RESET for details."
    FAILED=1
fi

# 3. Lint toolchain (Python)
log "[$CYAN 3/4$COLOR_RESET] Running$MAGENTA toolchain lint$COLOR_RESET..."
if ./mfc.sh lint > /dev/null 2>&1; then
    ok "Toolchain lint passed."
else
    error "Toolchain lint failed. Run$MAGENTA ./mfc.sh lint$COLOR_RESET for details."
    FAILED=1
fi

# 4. Source code lint checks
log "[$CYAN 4/4$COLOR_RESET] Running$MAGENTA source lint$COLOR_RESET checks..."
SOURCE_FAILED=0

# Check for raw OpenACC/OpenMP directives
if grep -qiR '!\$acc\|!\$omp' --exclude="parallel_macros.fpp" --exclude="acc_macros.fpp" --exclude="omp_macros.fpp" --exclude="shared_parallel_macros.fpp" --exclude="syscheck.fpp" ./src/* 2>/dev/null; then
    error "Found raw OpenACC/OpenMP directives. Use macros instead."
    SOURCE_FAILED=1
fi

# Check for double precision intrinsics
if grep -qiR 'double_precision\|dsqrt\|dexp\|dlog\|dble\|dabs\|double\ precision\|real(8)\|real(4)\|dprod\|dmin\|dmax\|dfloat\|dreal\|dcos\|dsin\|dtan\|dsign\|dtanh\|dsinh\|dcosh\|d0' --exclude-dir=syscheck --exclude="*nvtx*" --exclude="*precision_select*" ./src/* 2>/dev/null; then
    error "Found double precision intrinsics. Use generic intrinsics."
    SOURCE_FAILED=1
fi

# Check for junk code patterns
if grep -qiR -e '\.\.\.' -e '\-\-\-' -e '===' ./src/* 2>/dev/null; then
    error "Found junk code patterns (..., ---, ===) in source."
    SOURCE_FAILED=1
fi

if [ $SOURCE_FAILED -eq 0 ]; then
    ok "Source lint passed."
else
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
