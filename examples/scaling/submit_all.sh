#!/usr/bin/env sh

# Mandatory argument
ACCOUNT=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --account)
            ACCOUNT="$2"
            shift
            shift
            ;;
    esac
done

# Check mandatory argument
if [[ -z "$ACCOUNT" ]]; then
    echo "Error: --account is required"
    exit 1
fi

./examples/scaling/submit_weak.sh --account $ACCOUNT
./examples/scaling/submit_strong.sh --account $ACCOUNT
./examples/scaling/submit_grind.sh --account $ACCOUNT

