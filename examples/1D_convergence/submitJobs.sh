#!/usr/bin/env bash
set -e  # Exit on error
set -u  # Exit on undefined variable

Nx=(32 64 128 256 512 1024)
Order=(1 3 5)

ME=2 # Model equations = 2 for five-equation model
RS=2 # Riemann solver = 2 for HLLC

# Get the directory of the script itself
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Assume the script is in examples/1D_convergence, so MFC_DIR is two levels up
MFC_DIR="$(dirname "$(dirname "$ROOT_DIR")")"

for i in "${Nx[@]}"; do
    for j in "${Order[@]}"; do
        rm -rf "N${i}_O${j}"
        mkdir -p "N${i}_O${j}"
        cp case.py "N${i}_O${j}/"
    done
done

cd "$MFC_DIR" || exit 1

for i in "${Nx[@]}"; do
    for j in "${Order[@]}"; do
        ./mfc.sh run "$ROOT_DIR/N${i}_O${j}/case.py" --case-optimization --no-debug -- --order "$j" -N "$i" --meqns "$ME" --rs "$RS" || {
            echo "Error: mfc.sh failed for N=${i}, Order=${j}" >&2
            exit 1
        }
    done
done
