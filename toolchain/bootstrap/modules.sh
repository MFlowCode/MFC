#!/bin/bash

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
        -*|--*)        echo "Unknown option $1"; return; ;;
    esac
done

# Get computer (if not supplied in command-line)
if [ -v $u_c ]; then
    log   "Select a system:"
    log   "$G""ORNL$W:    Ascent     (a) | Frontier (f) | Summit (s) | Wombat (w)"
    log   "$C""ACCESS$W:  Bridges2   (b) | Expanse (e) | Delta  (d)"
    log   "$Y""Gatech$W:  Phoenix    (p)"
    log   "$R""Caltech$W: Richardson (r)"
    log_n "($G""a$W/$G""f$W/$G""s$W/$G""w$W/$C""b$W/$C""e$CR/$C""d$CR/$Y""p$CR/$R""r$CR): "
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

if [ "$u_cg" == 'c' ] || [ "$u_cg" == 'cpu' ]; then
    CG='CPU'; cg='cpu'
elif [ "$u_cg" == "g" ] || [ "$u_cg" == 'gpu' ]; then
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

# Reset modules to default system configuration
if [ "$u_c" != 'p' ]; then
    module reset > /dev/null 2>&1
    code="$?"

    # Purge if reset is not available
    if [ "$code" -ne '0' ]; then
        module purge > /dev/null 2>&1
    fi
else
    module purge > /dev/null 2>&1
fi

ELEMENTS=($(__extract "$u_c-all") $(__extract "$u_c-$cg"))

for element in ${ELEMENTS[@]}; do
    if [[ "$element" != *'='* ]]; then
        log " $ module load $M$element$CR"
        module load "$element" > /dev/null 2>&1

        # Handle Success / Failure
        code=$?
        if [ "$code" != '0' ]; then
            error "Failed to load module $M$element$CR:"

            # Run load again to show error message
            module load "$element"

            return
        fi
    else
        log " $ export $M$element$CR"
        export $element
    fi
done

if [ ! -z ${CRAY_LD_LIBRARY_PATH+x} ]; then
    ok "Found $M\$CRAY_LD_LIBRARY_PATH$CR. Prepending to $M\$LD_LIBRARY_PATH$CR."
    export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"
fi

ok 'All modules and environment variables have been loaded.'

return
