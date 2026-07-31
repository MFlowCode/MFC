#!/bin/bash

# Function to display help message
show_help() {
  echo "Usage: $(basename "$0") [OPTIONS]"
  echo "This function loads in modules for a specific computer and mode"
  echo ""
  echo "Options:"
  echo "  -h, --help                  Display this help message and exit."
  echo "  -c, --computer COMPUTER     Configures for COMPUTER environment."
  echo "                 Options:     Ascent (a) | Frontier (f) | Frontier_amd (famd) | Summit (s) | Wombat (w)"
  echo "                              AMD HPCFund (amdfund)"
  echo "                              Bridges2 (b) | Expanse (e) | Delta (d) | DeltaAI (dai)"
  echo "                              Phoenix (p) | Richardson (r) | Oscar (o)"
  echo "                              Carpenter Cray (cc) | Carpenter GNU (c) |  Nautilus (n)"
  echo "                              Isaac (i)"
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
    log   "$G""ORNL$W:    Ascent     (a) | Frontier (f) | Frontier_amd (famd) | Summit (s) | Wombat (w)"
    log   "$B""LLNL $W:   Tuolumne   (tuo)"
    log   "$C""ACCESS$W:  Bridges2   (b) | Expanse (e) | Delta  (d) | DeltaAI (dai)"
    log   "$Y""Gatech$W:  Phoenix    (p)"
    log   "$R""Caltech$W: Richardson (r)"
    log   "$BR""Brown$W: Oscar (o)"
    log   "$BR""Purdue$W: Anvil (pa)"
    log   "$OR""Tennessee$W: Isaac (i)"
    log   "$B""DoD$W:     Carpenter Cray (cc) | Carpenter GNU (c) |  Nautilus (n)"
    log   "$OR""Florida$W: HiPerGator (h)"
    log   "$C""WPI $W:   Turing   (t)"
    log   "$R""AMD$W:     HPCFund  (amdfund)"
    log_n "(${G}a${W}/${G}f${W}/${G}s${W}/${G}w${W}/${B}tuo${W}/${C}b${W}/${C}e${CR}/${C}d/${C}dai${CR}/${Y}p${CR}/${R}r${CR}/${B}cc${CR}/${B}c${CR}/${B}n${CR}/${BR}o${CR}/${BR}pa${CR}/${OR}i${CR}/${OR}h${CR}/${C}t${CR}/${R}amdfund${CR}): "
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

log " $ module load $MODULES"
if ! module load $MODULES; then
    error "Failed to load modules."

    return
fi

# Export variables one line at a time so each can reference previously exported vars
# (e.g. PATH="${OLCF_AFAR_ROOT}/..." requires OLCF_AFAR_ROOT to already be set)
for _suffix in "all" "$cg"; do
    while IFS= read -r _entry; do
        if echo "$_entry" | grep -q '='; then
            log " $ export $(eval "echo \"$_entry\"")"
            eval "export $_entry"
        fi
    done < <(grep -E "^$u_c-$_suffix\s+" toolchain/modules | sed "s/^$u_c-$_suffix\s\+//")
done
unset _suffix _entry

UNLOAD_MODULES="$(__extract "$u_c-all-unload") $(__extract "$u_c-$cg-unload")"
UNLOAD_MODULES=$(echo "$UNLOAD_MODULES" | xargs)
if [ -n "$UNLOAD_MODULES" ]; then
    log " $ module unload $UNLOAD_MODULES"
    module unload $UNLOAD_MODULES
fi

# Don't check for Cray paths on Carpenter, otherwise do check if they exist
if [ ! -z ${CRAY_LD_LIBRARY_PATH+x} ] && [ "$u_c" '!=' 'c' ] &&  [ "$u_c" '!=' 'famd' ] ; then
    ok "Found $M\$CRAY_LD_LIBRARY_PATH$CR. Prepending to $M\$LD_LIBRARY_PATH$CR."
    export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"
fi

if [ "$u_c" '==' 'famd' ]; then 
    export OLCF_AFAR_ROOT="/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6"

    export PATH=${OLCF_AFAR_ROOT}/lib/llvm/bin:${PATH}
    export LD_LIBRARY_PATH=${OLCF_AFAR_ROOT}/lib:${OLCF_AFAR_ROOT}/lib/llvm/lib:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/rocm-7.0.2/lib/llvm/lib:/opt/rocm-7.0.2/lib/

    export CRAY_MPICH_INC="-I${OLCF_AFAR_ROOT}/include/mpich3.4a2"
    export CRAY_HIPFORT_INC="-I${OLCF_AFAR_ROOT}/include/hipfort/amdgcn"
    export CRAY_HIPFORT_LIB="-L${OLCF_AFAR_ROOT}/lib -lhipfort-amdgcn -lhipfft -lamdhip64"
    export CRAY_HIP_INC="-I${OLCF_AFAR_ROOT}/include/hip"
    export CRAY_MPICH_LIB="-L${CRAY_MPICH_PREFIX}/lib \
                        ${CRAY_PMI_POST_LINK_OPTS} \
                        -lmpifort_amd -lmpi_amd -lmpi -lpmi -lpmi2"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${CRAY_LD_LIBRARY_PATH}"
    export CMAKE_PREFIX_PATH="${OLCF_AFAR_ROOT}:${CMAKE_PREFIX_PATH}"
    export FC="${OLCF_AFAR_ROOT}/bin/amdflang"

    unset MPICH_GPU_SUPPORT_ENABLED
fi

# AMD HPCFund: AFAR LLVM Flang (amdflang) + system GNU-built OpenMPI, gfx90a.
# The stock mpi.mod is GNU-built and binary-incompatible with amdflang, so we
# point at a pre-generated flang-compatible one under $OLCF_AFAR_ROOT/include/mpi.
# amdflang is only supported for OpenMP-offload GPU builds; CPU-mode (-m c) builds
# use the gfortran + OpenMPI provided by the amdfund-all modules, so gate on GPU.
if [ "$u_c" '==' 'amdfund' ] && [ "$cg" '==' 'gpu' ]; then
    # Direct WORK-filesystem path (mounted at the same /work1 point on login and
    # compute nodes); override OLCF_AFAR_ROOT before loading to use another drop.
    export OLCF_AFAR_ROOT="${OLCF_AFAR_ROOT:-/work1/spencerbryngelson/sbryngelson/software/therock-afar-23.2.1-gfx90a-7.13.0-7357b5084b}"
    # Track the loaded openmpi4 module's prefix (OpenHPC sets MPI_DIR) instead of
    # pinning a path that can drift from the module.
    _mfc_mpi_lib="${MPI_DIR:-/opt/ohpc/pub/mpi/openmpi4-gnu12/4.1.8}/lib"

    export PATH="${OLCF_AFAR_ROOT}/lib/llvm/bin:${OLCF_AFAR_ROOT}/bin:${PATH}"
    export LD_LIBRARY_PATH="${OLCF_AFAR_ROOT}/lib:${OLCF_AFAR_ROOT}/lib/llvm/lib:${_mfc_mpi_lib}:${LD_LIBRARY_PATH}"
    export CMAKE_PREFIX_PATH="${OLCF_AFAR_ROOT}:${CMAKE_PREFIX_PATH}"

    export MFC_FLANG_MPI_INC="-I${OLCF_AFAR_ROOT}/include/mpi"
    export MFC_FLANG_MPI_LIB="${_mfc_mpi_lib}/libmpi.so;${_mfc_mpi_lib}/libmpi_mpifh.so"
    export MFC_FLANG_HIPFORT_LIB="${OLCF_AFAR_ROOT}/lib/llvm/lib/libhipfort-amdgcn.a;${OLCF_AFAR_ROOT}/lib/libhipfft.so"

    export FC="${OLCF_AFAR_ROOT}/bin/amdflang"

    unset _mfc_mpi_lib

    if [ ! -x "$FC" ]; then
        error "amdfund: amdflang not found at $M$FC$CR."
        error "Set $M\$OLCF_AFAR_ROOT$CR to your AFAR drop, then reload, e.g.:"
        error "  export OLCF_AFAR_ROOT=/path/to/therock-afar-<ver>-gfx90a-...; source ./mfc.sh load -c amdfund -m g"
        return
    fi
fi

ok 'All modules and environment variables have been loaded.'

return
