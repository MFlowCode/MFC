#!/bin/bash
# Generate case constraints documentation from case_validator.py

set -e

REPO_ROOT="$1"

if [ -z "$REPO_ROOT" ]; then
    echo "Usage: $0 <repo_root>"
    exit 1
fi

echo "Generating case constraints documentation..."
python3 "$REPO_ROOT/toolchain/mfc/gen_case_constraints_docs.py" > "$REPO_ROOT/docs/documentation/case_constraints.md"
echo "✓ Generated docs/documentation/case_constraints.md"


