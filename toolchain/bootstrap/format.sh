#!/bin/bash

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            JOBS="$2"
            shift
            ;;
        *)
            echo "Format, unknown argument: $1."
            exit 1
            ;;
    esac

    shift
done

log "Formatting MFC:"

if ! find ${@:-src} -type f | grep -Ev 'autogen' | grep -E '\.(f90|fpp)$' | \
         parallel --jobs ${JOBS:-1} -- \
             echo      "\> {}" \&\& \
             python3   toolchain/indenter.py "{}" \&\& \
             fprettify "{}" --silent --indent 4 --c-relations --enable-replacements \
                       --enable-decl --whitespace-comma 1 --whitespace-multdiv 0    \
                       --whitespace-plusminus 1 --case 1 1 1 1 --strict-indent      \
                       --line-length 1000\;; then
    error "Formatting MFC failed."
    exit 1
fi

ok "Done. MFC has been formatted."

