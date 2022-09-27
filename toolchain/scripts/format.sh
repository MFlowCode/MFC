#!/usr/bin/env bash

if [ ! -f "$(pwd)/toolchain/util.sh" ]; then
    echo ""
    echo "[toolchain/scripts/prettify.sh] Please execute this script from the top-level MFC directory (currently in: $(pwd))."

    exit 1
fi

source "$(pwd)/toolchain/util.sh"

if ! command -v fprettify; then
    pip3 install --upgrade fprettify
fi

fprettify src --exclude "src/*/autogen" --recursive \
    --indent 4 --c-relations --enable-replacements --enable-decl \
    --whitespace-comma 1 --whitespace-multdiv 0 --whitespace-plusminus 1 \
    --case 1 1 1 1 --strict-indent
ret="$?"

if [ "$ret" != "0" ]; then
    echo -en "$RED"
    echo "[mfc.sh] Error: failed to execute fprettify."
    echo "[mfc.sh] Error: MFC has not been fprettify'ied."
    echo -en "$COLOR_RESET"

    exit 1
fi

echo -en "$GREEN"
echo "[mfc.sh] MFC has been fprettify'ied."
echo -en "$COLOR_RESET"
