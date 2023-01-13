#!/bin/bash

if [ -t 1 ]; then
    RED="\x1B[31m";    CYAN="\x1B[36m";    GREEN="\x1B[32m"
    YELLOW="\x1B[33m"; MAGENTA="\x1B[35m"; COLOR_RESET="\033[m"

    R=$RED;    C=$CYAN;    G=$GREEN
    Y=$YELLOW; M=$MAGENTA; CR=$COLOR_RESET; W=$CR
fi

log()   { echo -e  "$CYAN"mfc"$COLOR_RESET: $1$COLOR_RESET"; }
log_n() { echo -en "$CYAN"mfc"$COLOR_RESET: $1$COLOR_RESET"; }
ok()    { log "$GREEN""OK$COLOR_RESET > $1";                 }
warn()  { log "$YELLOW""WARNING$COLOR_RESET > $1";           }
error() { log "$RED""ERROR$COLOR_RESET > $1";                }
