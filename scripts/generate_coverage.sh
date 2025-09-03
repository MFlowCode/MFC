#!/bin/bash

# MFC Coverage Generation Script
# This script generates coverage reports locally using gcovr with the same
# configuration as the CI workflow.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MFC_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "MFC Coverage Generation"
echo "======================"
echo "Root directory: $MFC_ROOT"
echo

# Check if gcovr is available
if ! command -v gcovr &> /dev/null; then
    echo "Error: gcovr is not installed. Install it with:"
    echo "  pip install gcovr"
    exit 1
fi

cd "$MFC_ROOT"

# Check if MFC was built with coverage
if [ ! -d "build" ]; then
    echo "Error: No build directory found. Build MFC with coverage first:"
    echo "  ./mfc.sh build --gcov -j \$(nproc)"
    exit 1
fi

# Look for gcov files
GCOV_FILES=$(find build -name "*.gcda" 2>/dev/null | wc -l)
if [ "$GCOV_FILES" -eq 0 ]; then
    echo "Warning: No coverage data found. Run tests first:"
    echo "  ./mfc.sh test -j \$(nproc) -f EA8FA07E -t 9E2CA336"
    echo
fi

echo "Generating coverage reports..."

# Generate XML coverage report (same as CI)
gcovr \
    --root "$MFC_ROOT" \
    --object-directory "$MFC_ROOT/build" \
    --filter '.*/fypp/.+\.f90$' \
    --exclude '.*(tests|examples|toolchain|docs|benchmarks)/.*' \
    --gcov-ignore-parse-errors \
    --xml-pretty -o coverage.xml

echo "Generated coverage.xml"

# Also generate a broader report for debugging
gcovr \
    --root "$MFC_ROOT" \
    --object-directory "$MFC_ROOT/build" \
    --exclude '.*(tests|examples|toolchain|docs|benchmarks)/.*' \
    --gcov-ignore-parse-errors \
    --xml-pretty -o coverage_all.xml

echo "Generated coverage_all.xml (broader scope for debugging)"

# Generate HTML coverage report for local viewing
gcovr \
    --root "$MFC_ROOT" \
    --object-directory "$MFC_ROOT/build" \
    --filter '.*/fypp/.+\.f90$' \
    --exclude '.*(tests|examples|toolchain|docs|benchmarks)/.*' \
    --gcov-ignore-parse-errors \
    --html-details -o coverage.html

echo "Generated coverage.html"

# Rewrite paths from generated .fpp.f90 -> repo .fpp (same as CI)
python3 - <<'PY'
import os, xml.etree.ElementTree as ET

print("Rewriting XML paths from generated .fpp.f90 to source .fpp files...")

fn = "coverage.xml"
if not os.path.exists(fn):
    print(f"Error: {fn} not found")
    exit(1)

t = ET.parse(fn)
r = t.getroot()
remapped = 0

for c in r.iter('class'):
    p = c.get('filename') or ''
    if p.endswith('.fpp.f90'):
        base = os.path.basename(p).rsplit('.fpp.f90', 1)[0] + '.fpp'
        for sub in ('pre_process','simulation','post_process','common'):
            cand = os.path.join('src', sub, base)
            if os.path.exists(cand):
                c.set('filename', cand)
                remapped += 1
                break

t.write(fn, encoding='utf-8', xml_declaration=True)
print(f"Remapped {remapped} file paths in coverage.xml")
PY

echo
echo "Coverage reports generated:"
echo "  - coverage.xml (for Codecov upload)"
echo "  - coverage.html (for local viewing)"
echo
echo "To view the HTML report, open coverage.html in your browser."
