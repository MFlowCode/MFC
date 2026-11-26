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
  echo "      --staged          Skip the checks no staged file can affect (the pre-commit hook)."
  echo ""
  exit 0
}


JOBS=1
STAGED=0

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
        --staged)
            STAGED=1
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
# Skip tests that compile Fortran: they depend on the modules loaded in the committer's
# shell (./mfc.sh load), not on the commit. CI still runs them via ./mfc.sh lint.
export MFC_SKIP_COMPILER_TESTS=1
# Skip the test-suite linters: they materialize every test case, which takes about five of the
# seven minutes this check otherwise needs. CI still runs them via ./mfc.sh lint.
export MFC_SKIP_SUITE_LINT=1

NCHECK=7

# --staged: the unit tests read src/ and examples/ as well, so only a docs-only commit skips them;
# the example decks depend only on examples/ and toolchain/. CI runs everything.
LINT_ARGS=()
CHECK_EXAMPLES=1
if [ "$STAGED" = 1 ]; then
    staged=$(git diff --cached --name-only)
    if [ -n "$staged" ] && ! grep -qv '^docs/' <<< "$staged"; then LINT_ARGS=(--no-test); fi
    grep -qE '^(examples|toolchain)/' <<< "$staged" || CHECK_EXAMPLES=0
fi

log "Running$MAGENTA precheck$COLOR_RESET (same checks as CI lint-gate)..."
echo ""

# Temp files for collecting results from parallel jobs
TMPDIR_PC=$(mktemp -d)
export TMPDIR_PC
trap "rm -rf $TMPDIR_PC" EXIT

# --- Phase 1: Format check (non-mutating; mirrors CI's format+diff check) ---
# Use --check mode so staged-but-unformatted code is caught even if the
# working tree was already reformatted by a prior ./mfc.sh format run.
FORMAT_OK=0
if ! ruff format --check toolchain/ examples/ benchmarks/ > /dev/null 2>&1; then
    FORMAT_OK=2
elif ! ffmt --check -j "$JOBS" src/ > /dev/null 2>&1; then
    FORMAT_OK=2
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
    if ./mfc.sh lint "${LINT_ARGS[@]}" > /dev/null 2>&1; then
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
(
    if python3 toolchain/mfc/lint_docs.py > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/doc_exit"
    else
        echo "1" > "$TMPDIR_PC/doc_exit"
    fi
) &
PID_DOC=$!

# Parameter documentation check
(
    if python3 toolchain/mfc/lint_param_docs.py > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/param_docs_exit"
    else
        echo "1" > "$TMPDIR_PC/param_docs_exit"
    fi
) &
PID_PARAM_DOCS=$!

# Example case validation - validate ALL cases in ONE process (validate accepts multiple files) instead of a
# per-case './mfc.sh validate' loop. Each invocation re-paid the mfc.sh/venv bootstrap (~2.4s), which dominated
# (~90%) the ~6 min serial loop over 155 cases; batching drops it to well under a minute.
# Every example directory holds its deck as case.py, so the glob below sees it (scaling/ is a benchmark harness).
(
    [ "$CHECK_EXAMPLES" = 1 ] || { echo "skip" > "$TMPDIR_PC/examples_exit"; exit 0; }
    missing=$(git ls-files 'examples/*/*' | cut -d/ -f2 | sort -u | while read -r d; do
        [ -f "examples/$d/case.py" ] || [ "$d" = scaling ] || echo "examples/$d"
    done)
    if [ -n "$missing" ]; then echo "No case.py in: $missing" > "$TMPDIR_PC/examples_missing"; fi
    if [ -z "$missing" ] && ./mfc.sh validate examples/*/case.py > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/examples_exit"
    else
        echo "1" > "$TMPDIR_PC/examples_exit"
    fi
) &
PID_EXAMPLES=$!

# --- Collect results (fast checks) ---

FAILED=0

log "[$CYAN 1/$NCHECK$COLOR_RESET] Checking$MAGENTA formatting$COLOR_RESET..."
if [ "$FORMAT_OK" = "2" ]; then
    error "Code is not formatted. Run$MAGENTA ./mfc.sh format$COLOR_RESET and re-stage the changes."
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
    ok "Toolchain lint passed${LINT_ARGS:+ (unit tests skipped: only docs/ staged)}."
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

wait $PID_DOC
log "[$CYAN 5/$NCHECK$COLOR_RESET] Checking$MAGENTA doc references$COLOR_RESET..."
if [ "$(cat "$TMPDIR_PC/doc_exit" 2>/dev/null || echo 1)" = "0" ]; then
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

wait $PID_EXAMPLES
log "[$CYAN 7/$NCHECK$COLOR_RESET] Validating$MAGENTA example cases$COLOR_RESET..."
EXAMPLES_FAILED=$(cat "$TMPDIR_PC/examples_exit" 2>/dev/null || echo "1")
if [ "$EXAMPLES_FAILED" = "skip" ]; then
    ok "Example cases skipped: nothing staged under examples/ or toolchain/."
elif [ "$EXAMPLES_FAILED" = "0" ]; then
    ok "All example cases are valid."
else
    [ -f "$TMPDIR_PC/examples_missing" ] && error "$(cat "$TMPDIR_PC/examples_missing") (name the deck case.py)."
    error "Example case validation failed. Run$MAGENTA ./mfc.sh validate examples/\*/case.py$COLOR_RESET for the failing case(s)."
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
