#!/bin/bash

. toolchain/util.sh

echo "> $1"

niter=0
old_file=""

while :
do
    niter=$((niter+1))
    new_file=`cat "$1"`
    if [[ "$new_file" == "$old_file" ]]; then
        break
    fi
    old_file="$new_file"

    if [[ "$niter" -gt 4 ]]; then
        error "Failed to format $1: No steady-state (after $niter iterations)."
        exit 1
    fi

    if ! python3 toolchain/indenter.py "$1"; then
        error "Failed to format $1: indenter.py."
        exit 1
    fi

    if ! fprettify "$1" --silent --indent 4 --c-relations --enable-replacements \
            --enable-decl --whitespace-comma 1 --whitespace-multdiv 0           \
            --whitespace-plusminus 1 --case 1 1 1 1 --strict-indent             \
            --line-length 1000; then
        error "Failed to format $1: fprettify."
        exit 1
    fi
done

