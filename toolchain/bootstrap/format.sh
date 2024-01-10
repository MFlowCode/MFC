#!/bin/bash

log "Formatting MFC:"

fortran_files=$(find ${@:-src} -type f | grep -Ev 'src/.+/autogen/')

longest=0
for filepath in $fortran_files; do
    if [ "${#filepath}" -gt "$longest" ]; then
        longest="${#filepath}"
    fi
done

for filepath in $fortran_files; do
    echo -n " > $filepath $(printf '%*s' "$((longest - ${#filepath}))" '')"

    before=$(sha256sum "$filepath" | cut -d' ' -f1)

    python3 toolchain/indenter.py "$filepath"

    if ! fprettify "$filepath" \
        --silent --indent 4 --c-relations --enable-replacements --enable-decl \
        --whitespace-comma 1 --whitespace-multdiv 0 --whitespace-plusminus 1 \
        --case 1 1 1 1 --strict-indent --line-length 1000; then
        error "failed to execute fprettify."
        error "MFC has not been fprettify'ied."
        exit 1
    fi

    after=$(sha256sum "$filepath" | cut -d' ' -f1)

    if [ "$before" != "$after" ]; then
        echo -e "$YELLOW[formatted]$COLOR_RESET"
    else
        echo -e "$GREEN[unchanged]$COLOR_RESET"
    fi
done

ok "Done. MFC has been formatted."
