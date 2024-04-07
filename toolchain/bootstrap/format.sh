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

if ! find ${@:-src} -type f | grep -Ev 'autogen' | grep -E '\.(f90|fpp)$' \
        | xargs -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_file.sh; then
    error "Formatting MFC failed."
    exit 1
fi

ok "Done. MFC has been formatted."

