#!/bin/bash
# Generate parameters documentation from parameter registry

set -e

REPO_ROOT="$1"

if [ -z "$REPO_ROOT" ]; then
    echo "Usage: $0 <repo_root>"
    exit 1
fi

echo "Generating parameters documentation..."
cd "$REPO_ROOT/toolchain"
python3 -c "
import sys; sys.path.insert(0, '.')
from mfc.params.generators.docs_gen import generate_parameter_docs
print(generate_parameter_docs())
" > "$REPO_ROOT/docs/documentation/parameters.md"
echo "Generated docs/documentation/parameters.md"
