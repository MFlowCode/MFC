#!/bin/bash
# Generate CLI reference documentation from cli/commands.py

set -e

REPO_ROOT="$1"

if [ -z "$REPO_ROOT" ]; then
    echo "Usage: $0 <repo_root>"
    exit 1
fi

echo "Generating CLI reference documentation..."
cd "$REPO_ROOT/toolchain"
python3 -c "
import sys; sys.path.insert(0, '.')
from mfc.cli.commands import MFC_CLI_SCHEMA
from mfc.cli.docs_gen import generate_cli_reference
print(generate_cli_reference(MFC_CLI_SCHEMA))
" > "$REPO_ROOT/docs/documentation/cli-reference.md"
echo "Generated docs/documentation/cli-reference.md"
