#!/bin/bash

log "Formatting MFC with fprettify..."

fprettify ${@:-src} --exclude "src/*/autogen" --recursive --silent \
    --indent 4 --c-relations --enable-replacements --enable-decl \
    --whitespace-comma 1 --whitespace-multdiv 1 --whitespace-plusminus 1 \
    --case 1 1 1 1 --strict-indent
ret="$?"

if [ "$ret" != '0' ]; then
    error "failed to execute fprettify."
    error "MFC has not been fprettify'ied."

    exit 1
fi

ok "MFC has been fprettify'ied."