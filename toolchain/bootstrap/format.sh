#!/bin/bash

# Function to display help message
show_help() {
  echo "Usage: $(basename "$0") [OPTIONS]"
  echo "This function formats all code in and below the current execution directory."
  echo ""
  echo "Options:"
  echo "  -h, --help            Display this help message and exit."
  echo "  -j, --jobs JOBS       Runs JOBS number of jobs."
  echo ""
  exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            JOBS="$2"
            shift
            ;;
        -h | --help)
            show_help
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
    error "Formatting MFC source failed."
    exit 1
fi

if ! find ${@:-examples} -type f | grep -E '\.(py)$' \
        | xargs -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_python.sh; then
    error "Formatting MFC examples failed."
    exit 1
fi

if ! find ${@:-benchmarks} -type f | grep -E '\.(py)$' \
        | xargs -L 1 -P ${JOBS:-1} $SHELL toolchain/bootstrap/format_python.sh; then
    error "Formatting MFC examples failed."
    exit 1
fi

ok "Done. MFC has been formatted."
