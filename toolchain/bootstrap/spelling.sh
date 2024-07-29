#!/bin/bash

log "Running$MAGENTA typos$COLOR_RESET on$MAGENTA MFC$COLOR_RESET's $MAGENTA""toolchain$COLOR_RESET..."

if typos --config .typos.toml; then
    ok "Done. No typos found."
else
    error "Typos found."
    exit 1
fi