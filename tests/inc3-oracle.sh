#!/usr/bin/env bash
#
# tests/inc3-oracle.sh -- dev oracle for L0/AMR Increment-3 coexist
#
# Runs a "coexist" case (amr=T .and. l0_ntile>0, guard opens in a later Increment-3 task)
# against a "monolithic" reference case (the single-mode run it should reproduce) at
# np=1, then compares their final-step D/*.dat output. Modeled on the tolerance compare
# in toolchain/mfc/packer/{pack,tol}.py: each D/*.dat line is "<x> [<y> <z>] <value>",
# so only the trailing <value> column is compared, at an absolute tolerance of 1e-13.
#
# This is a developer convenience script, not wired into ./mfc.sh test: it produces no
# golden files and is not run by CI.
#
# Usage: tests/inc3-oracle.sh <coexist_case> <monolithic_case>
#   <case> may be a case directory or a path to its case.py.

set -euo pipefail

TOL=1e-13

usage() {
    echo "usage: $0 <coexist_case> <monolithic_case>" >&2
    exit 1
}

[[ $# -eq 2 ]] || usage

# Resolve a case argument (directory or case.py path) to both its case.py and its
# containing directory (where D/ will be written).
resolve_case_py() {
    local c="$1"
    if [[ -f "$c" ]]; then
        echo "$c"
    elif [[ -f "$c/case.py" ]]; then
        echo "$c/case.py"
    else
        echo "error: no case.py found at or under '$c'" >&2
        exit 1
    fi
}

resolve_case_dir() {
    local c="$1"
    if [[ -d "$c" ]]; then echo "$c"; else dirname "$c"; fi
}

COEXIST_CASE=$(resolve_case_py "$1")
MONO_CASE=$(resolve_case_py "$2")
COEXIST_DIR=$(resolve_case_dir "$1")
MONO_DIR=$(resolve_case_dir "$2")

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT"

echo "==> running coexist case:    $COEXIST_CASE"
./mfc.sh run "$COEXIST_CASE" -n 1 --gpu mp -g 0 -t pre_process simulation

echo "==> running monolithic case: $MONO_CASE"
./mfc.sh run "$MONO_CASE" -n 1 --gpu mp -g 0 -t pre_process simulation

COEXIST_D="$COEXIST_DIR/D"
MONO_D="$MONO_DIR/D"

[[ -d "$COEXIST_D" ]] || { echo "error: $COEXIST_D not found (did the run fail?)" >&2; exit 1; }
[[ -d "$MONO_D" ]] || { echo "error: $MONO_D not found (did the run fail?)" >&2; exit 1; }

# D/*.dat filenames are <field>.<var>.<rank>.<step>.dat (see src/simulation/m_data_output.fpp,
# e.g. cons.1.00.000000.dat). For each (field, var, rank) key, keep only the file at the
# largest step -- i.e. the final write.
final_step_files() {
    local d="$1"
    find "$d" -maxdepth 1 -regextype posix-extended \
        -regex '.*/(cons|prim)\.[0-9]+\.[0-9]{2}\.[0-9]{6}\.dat' -printf '%f\n' |
        awk -F. '
            { key = $1"."$2"."$3; step = $4 + 0
              if (!(key in best) || step > best[key]) { best[key] = step; file[key] = $0 }
            }
            END { for (k in file) print file[k] }
        '
}

max_delta="0"
n_compared=0

while read -r fname; do
    [[ -z "$fname" ]] && continue
    mine="$COEXIST_D/$fname"
    other="$MONO_D/$fname"
    if [[ ! -f "$other" ]]; then
        echo "WARN: $fname present in coexist D/ but not in monolithic D/ (skipping)" >&2
        continue
    fi
    n_compared=$((n_compared + 1))
    d=$(paste "$mine" "$other" | awk '
        { nf = NF / 2; v1 = $nf; v2 = $NF
          diff = v1 - v2; if (diff < 0) diff = -diff
          if (diff > m) m = diff
        }
        END { printf "%.17e", m + 0 }
    ')
    max_delta=$(awk -v a="$max_delta" -v b="$d" 'BEGIN { print (a + 0 > b + 0) ? a : b }')
done < <(final_step_files "$COEXIST_D")

if [[ "$n_compared" -eq 0 ]]; then
    echo "error: no matching final-step D/*.dat files found to compare" >&2
    exit 1
fi

echo "compared $n_compared file(s); max abs delta = $max_delta (tol = $TOL)"
if awk -v m="$max_delta" -v t="$TOL" 'BEGIN { exit !(m + 0 <= t + 0) }'; then
    echo "PASS"
    exit 0
else
    echo "FAIL"
    exit 1
fi
