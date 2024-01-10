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
    log   "$G""ORNL$W:    Ascent     (a) | Crusher (c) | Summit (s) | Wombat (w)"
    log   "$C""ACCESS$W:  Bridges2   (b) | Expanse (e) | Delta  (d)"
    log   "$Y""GaTech$W:  Phoenix    (p)"
    log   "$R""CALTECH$W: Richardson (r)"
    log_n "($G""a$W/$G""c$W/$G""s$W/$G""w$W/$C""b$W/$C""e$CR/$C""d$CR/$Y""p$CR/$R""r$CR): "
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

MODULES=($(__extract "$u_c-$cg") $(__combine $(__extract "$u_c-all")))

log "Loading modules for $M$COMPUTER$CR on $M$CG$CR"'s:'

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

# Find length of longest module_name in $MODULES for $COMPUTER
max_module_length="0"
for module_name in ${MODULES[@]}; do
    module_length="${#module_name}"

    if [ "$module_length" -gt "$max_module_length" ]; then
        max_module_length="$module_length"
    fi
done

# Load modules ($MODULES)
for module_name in ${MODULES[@]}; do
    log_n " - $CYAN$module_name$COLOR_RESET "

    # Add padding spaces
    module_length="${#module_name}"
    delta="$((max_module_length-module_length-1))"
    if [ "$delta" -ge "0" ]; then
        printf "%0.s " $(seq 0 $delta)
    fi

    # Load the module
    module load "$module_name" > /dev/null 2>&1

    # Handle Success / Failure
    code=$?
    if [ "$code" == '0' ]; then
        echo -e "[$G""SUCCESS$W]"
    else
        echo -e "[$R""FAILURE$W]"

        # Run load again to show error message
        module load "$module_name"

        return
    fi
done

if [ "$cg" == 'gpu' ]; then
    isnv=$($FC --version | grep NVIDIA | wc -l)
    if [ $isnv -eq 0 ]; then
        export CC=nvc CXX=nvc++ FC=nvfortran
    fi
fi
ok 'All modules and environment variables have been loaded.'


return
