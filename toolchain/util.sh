#!/bin/bash

if [ -t 1 ]; then
    RED="\x1B[31m";    CYAN="\x1B[36m";    GREEN="\x1B[32m"
    YELLOW="\x1B[33m"; MAGENTA="\x1B[35m"; COLOR_RESET="\033[m"
fi

log()     { echo -e "$CYAN"mfc"$COLOR_RESET: $1$COLOR_RESET"; }
ok()      { log "$GREEN""OK$COLOR_RESET > $1";                }
warning() { log "$YELLOW""WARNING$COLOR_RESET > $1";          }
error()   { log "$RED""ERROR$COLOR_RESET > $1";               }
