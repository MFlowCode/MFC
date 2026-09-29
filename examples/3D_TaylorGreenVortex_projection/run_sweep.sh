#!/usr/bin/env bash
# Explicit vs. all-Mach projection sweep on the 3D Taylor-Green vortex, on 1, 2 and 4 GPUs.
#
#   timing:   a fixed number of steps per run for the cost per step. Explicit cost per step does not
#             depend on Mach, so it runs at EXPLICIT_MACHS only; its steps per tC scale as 1 + 1/M.
#   accuracy: full runs to ACC_TEND convective times on ACC_NGPU GPUs, compared with analyze.py ke.
#   weak:     WEAK_NS^3 cells per GPU on each of WEAK_NGPUS, the vortex repeated along x once per GPU;
#             analyze.py weak gives the efficiency against one GPU.
#
# usage: run_sweep.sh [timing|accuracy|weak|all] [-g GPU_IDS] [-o OUT_DIR] [--dry-run]
#
# Every run gets a new directory under OUT_DIR (default runs/<suite>_<date>); existing ones are never
# overwritten. Matrix entries are overridable through the environment, e.g.
#   NS="64" NGPUS="1 2" MACHS="0.01" ./run_sweep.sh timing
# MFC_ENV, if set, is sourced first (compiler/MPI environment). Runs wait for the chosen GPUs to be idle. Runs are
# case-optimized, as MFC is for performance (CASE_OPT=0 uses one generic build); each builds its binary on first use, so
# wall_s includes that build while s_per_step does not. RDMA=1 runs every case with GPU-aware MPI (case.py --rdma).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)

NS=${NS:-"128 256"}
NGPUS=${NGPUS:-"1 2 4"}
MACHS=${MACHS:-"0.1 0.03 0.01 0.003 0.001"}
EXPLICIT_MACHS=${EXPLICIT_MACHS:-"0.1 0.01"}
STEPS=${STEPS:-200}
ACC_N=${ACC_N:-128}
ACC_NGPU=${ACC_NGPU:-4}
ACC_TEND=${ACC_TEND:-10}
ACC_SAVES=${ACC_SAVES:-20}
ACC_MACHS=${ACC_MACHS:-"0.1 0.01 0.001"}
ACC_EXPLICIT_MACHS=${ACC_EXPLICIT_MACHS:-"0.1 0.01"}
WEAK_NS=${WEAK_NS:-"128 256"}
WEAK_NGPUS=${WEAK_NGPUS:-"1 2 4"}
WEAK_MACHS=${WEAK_MACHS:-"0.1 0.01"}
WEAK_EXPLICIT_MACHS=${WEAK_EXPLICIT_MACHS:-"0.01"}
CASE_OPT=${CASE_OPT:-1}
RDMA=${RDMA:-0}

suite=${1:-timing}
[[ $suite =~ ^(timing|accuracy|weak|all)$ ]] || { sed -n '2,18p' "$0"; exit 1; }
shift || true
gpus=0,1,2,3 out="" dry=0
while (($#)); do
    case $1 in
        -g) gpus=$2; shift 2 ;;
        -o) out=$2; shift 2 ;;
        --dry-run) dry=1; shift ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done
out=${out:-$HERE/runs/${suite}_$(date +%Y%m%d_%H%M%S)}
csv=$out/results.csv

[[ -n ${MFC_ENV:-} ]] && source "$MFC_ENV"

wait_idle() {  # $1: comma-separated GPU ids
    while [[ -n $(nvidia-smi -i "$1" --query-compute-apps=pid --format=csv,noheader) ]]; do
        echo "GPUs $1 busy, waiting"; sleep 60
    done
}

# run SOLVER N NGPU MACH [case.py args...]
run() {
    local solver=$1 N=$2 ng=$3 mach=$4; shift 4
    local d=$out/${solver}_N${N}_M${mach}_g${ng}
    local ids; ids=$(cut -d, -f1-"$ng" <<< "$gpus")
    local cargs=(--N "$N" --mach "$mach" "$@")
    [[ $solver == explicit ]] && cargs+=(--explicit)
    ((RDMA)) && cargs+=(--rdma)
    local spt; spt=$(python3 "$HERE/case.py" "${cargs[@]}" 2>&1 >/dev/null | sed -E 's/.*steps per tC = ([0-9.]+).*/\1/')
    echo "$solver N=$N GPUs=$ids Mach=$mach: ${cargs[*]}"
    ((dry)) && return
    mkdir "$d"  # fails rather than reuse a directory
    cp "$HERE/case.py" "$d/"
    wait_idle "$ids"
    local t0=$SECONDS rc=0
    local build=(--no-build)
    ((CASE_OPT)) && build=(--case-optimization -j 32)
    (cd "$ROOT" && CUDA_VISIBLE_DEVICES=$ids ./mfc.sh run "$d/case.py" --gpu acc --no-debug "${build[@]}" -n "$ng" \
        -t pre_process simulation -- "${cargs[@]}") > "$d/log" 2>&1 || rc=$?
    local avg; avg=$(sed 's/\x1b\[[0-9;]*m//g' "$d/log" | grep -oE 'avg +[0-9.]+E[-+][0-9]+' | tail -1 | awk '{print $2}')
    echo "$solver,$N,$ng,$mach,$spt,${avg:-},$((SECONDS - t0)),$rc" >> "$csv"
    echo "  rc=$rc  s/step=${avg:-?}  wall=$((SECONDS - t0))s"
}

if ((!dry)); then
    mkdir -p "$out"
    echo "solver,N,ngpu,mach,steps_per_tC,s_per_step,wall_s,rc" > "$csv"
    ((CASE_OPT)) || (cd "$ROOT" && ./mfc.sh build --gpu acc --no-debug -t pre_process simulation -j 32) > "$out/build.log" 2>&1
fi

if [[ $suite =~ ^(timing|all)$ ]]; then
    for N in $NS; do
        for ng in $NGPUS; do
            for M in $MACHS; do run projection "$N" "$ng" "$M" --steps "$STEPS" --saves 1; done
            for M in $EXPLICIT_MACHS; do run explicit "$N" "$ng" "$M" --steps "$STEPS" --saves 1; done
        done
    done
fi

if [[ $suite =~ ^(accuracy|all)$ ]]; then
    for M in $ACC_MACHS; do run projection "$ACC_N" "$ACC_NGPU" "$M" --tend "$ACC_TEND" --saves "$ACC_SAVES"; done
    for M in $ACC_EXPLICIT_MACHS; do run explicit "$ACC_N" "$ACC_NGPU" "$M" --tend "$ACC_TEND" --saves "$ACC_SAVES"; done
fi

if [[ $suite == weak ]]; then
    for N in $WEAK_NS; do
        for ng in $WEAK_NGPUS; do
            for M in $WEAK_MACHS; do run projection "$N" "$ng" "$M" --copies "$ng" --steps "$STEPS" --saves 1; done
            for M in $WEAK_EXPLICIT_MACHS; do run explicit "$N" "$ng" "$M" --copies "$ng" --steps "$STEPS" --saves 1; done
        done
    done
fi

((dry)) || echo "results: $csv  (python3 $HERE/analyze.py timing $csv; analyze.py ke $out/explicit_N${ACC_N}_M0.01_g${ACC_NGPU} $out/projection_N${ACC_N}_M*)"
