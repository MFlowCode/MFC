#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CASES=(
  "examples/1D_sodshocktube_muscl/case.py"
  "examples/1D_laxshocktube/case.py"
  "examples/1D_jwl_single_material_shocktube/case.py"
)

if [[ "${1:-}" == "--build" ]]; then
  ./mfc.sh build -j "${MFC_BUILD_JOBS:-2}"
fi

for case_cmd in "${CASES[@]}"; do
  echo "============================================================"
  echo "Running: $case_cmd"
  ./mfc.sh run $case_cmd --no-build
  echo "Completed: $case_cmd"
  echo "============================================================"
  echo
done

echo "All 1D verification cases completed."
