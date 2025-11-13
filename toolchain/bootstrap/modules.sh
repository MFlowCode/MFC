#!/bin/bash

# Function to display help message
show_help() {
  echo "Usage: $(basename "$0") [OPTIONS]"
  echo "This function loads in modules for a specific computer and mode"
  echo ""
  echo "Options:"
  echo "  -h, --help                  Display this help message and exit."
  echo "  -c, --computer COMPUTER     Configures for COMPUTER environment."
  echo "                 Options:     Ascent (a) | Frontier (f) | Summit (s) | Wombat (w)"
  echo "                              Bridges2 (b) | Expanse (e) | Delta (d) | DeltaAI (dai)"
  echo "                              Phoenix (p) | Richardson (r) | Oscar (o)"
  echo "                              Carpenter Cray (cc) | Carpenter GNU (c) |  Nautilus (n)"
  echo "  -m, --mode MODE             Configures into MODE."
  echo "                 Options:     gpu (g) | cpu (c)"
  echo ""
  exit 0
}

unset u_c
unset u_cg

# Reset u_computer & u_cg to known values since this script is run with "source"
# Therefore, values of variables defined here are kept when the script runs again.
u_computer=""; u_cg=""

# If there are command-line arguments, parse them:
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--computer) u_c="$2";  shift; shift; ;;
        -m|--mode)     u_cg="$2"; shift; shift; ;;
        -h|--help)     show_help; shift; shift; ;;
        -*|--*)        echo "Unknown option $1"; return; ;;
    esac
done

# Get computer (if not supplied in command line)
if [ -v $u_c ]; then
    log   "Select a system:"
    log   "$G""ORNL$W:    Ascent     (a) | Frontier (f) | Summit (s) | Wombat (w)"
    log   "$C""ACCESS$W:  Bridges2   (b) | Expanse (e) | Delta  (d) | DeltaAI (dai)"
    log   "$Y""Gatech$W:  Phoenix    (p)"
    log   "$R""Caltech$W: Richardson (r)"
    log   "$BR""Brown$W: Oscar (o)"
    log   "$B""DoD$W:     Carpenter Cray (cc) | Carpenter GNU (c) |  Nautilus (n)"
    log   "$OR""Florida$W: HiPerGator (h)"
    log_n "($G""a$W/$G""f$W/$G""s$W/$G""w$W/$C""b$W/$C""e$CR/$C""d/$C""dai$CR/$Y""p$CR/$R""r$CR/$B""cc$CR/$B""c$CR/$B""n$CR/$BR""o"$CR"/$OR""h"$CR"): "
    read u_c
    log
fi

# Get CPU/GPU (if not supplied in command-line)
if [ -v $u_cg ]; then
    log   "Select configuration:"
    log   " - CPU (c | cpu)"
    log   " - GPU (g | gpu)"
    log_n "(c/cpu/g/gpu): "
    read u_cg
    log
fi

# User input to lowercase
u_c=$(echo "$u_c" | tr '[:upper:]' '[:lower:]')
u_cg=$(echo "$u_cg" | tr '[:upper:]' '[:lower:]')

if [ "$u_cg" '==' 'c' ] || [ "$u_cg" '==' 'cpu' ]; then
    CG='CPU'; cg='cpu'
elif [ "$u_cg" '==' "g" ] || [ "$u_cg" '==' 'gpu' ]; then
    CG='GPU'; cg='gpu'
fi

__combine() {
    echo -n $@ | sed 's/\[\n\r\s\]\+/\ /'
}

__extract() {
    __combine "$(grep -E "^$1\s+" toolchain/modules | sed "s/^$1\s\+//")"
}

COMPUTER="$(__extract "$u_c")"

if [[ -z "$COMPUTER" ]]; then
    error "Computer $M$u_cg$CR not recognized."
    return
fi

log "Loading modules (& env variables) for $M$COMPUTER$CR on $M$CG$CR"'s:'

# Reset modules to default system configuration (unless Carpenter)
if [ "$u_c" '!=' 'cc' ] && [ "$u_c" '!=' 'c' ]; then
    module reset > /dev/null 2>&1
    code="$?"

    # Purge if reset is not available
    if [ "$code" -ne '0' ]; then
        module purge > /dev/null 2>&1
    fi
else
    module purge > /dev/null 2>&1
fi

ELEMENTS="$(__extract "$u_c-all") $(__extract "$u_c-$cg")"
MODULES=`echo "$ELEMENTS" | tr ' ' '\n' | grep -v = | xargs`
VARIABLES=`echo "$ELEMENTS" | tr ' ' '\n' | grep = | xargs`

log " $ module load $MODULES"
if ! module load $MODULES; then
    error "Failed to load modules."

    return
fi

if [ $(echo "$VARIABLES" | grep = | wc -c) -gt 0 ]; then
    log " $ export $(eval "echo $VARIABLES")"
    export $(eval "echo $VARIABLES") > /dev/null
fi

# Don't check for Cray paths on Carpenter, otherwise do check if they exist
if [ ! -z ${CRAY_LD_LIBRARY_PATH+x} ] && [ "$u_c" '!=' 'c' ]; then
    ok "Found $M\$CRAY_LD_LIBRARY_PATH$CR. Prepending to $M\$LD_LIBRARY_PATH$CR."
    export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"
fi

ok 'All modules and environment variables have been loaded.'

return
