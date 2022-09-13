#!/usr/bin/env bash

RED="\u001b[31m"
CYAN="\u001b[36m"
GREEN="\u001b[32m"
YELLOW="\u001b[33m"
MAGENTA="\u001b[35m"
COLOR_RESET="\033[m"

TIMEOUT="60s"

./mfc.sh build -j $(nproc) -t pre_process simulation

index=1
failed=()

files=samples/*/*.py

nFiles=$(find samples/*/*.py | wc -l) 

for f in samples/*/*.py; do

    echo -e "($index/$nFiles) Running $MAGENTA$f$COLOR_RESET:" 
    echo -en " - pre_process... "
    ./mfc.sh run "$f" -j $(nproc) -n 4 -b mpirun -t pre_process --no-build > /dev/null 2>&1
    
    code="$?"
    if [[ "$code" != "0" ]]; then
        echo -e "$RED""FAILED$COLOR_RESET"
        failed=(${failed[@]} "$f")
        index=$((index+1))
        continue
    else
        echo -e "$GREEN""PASSED$COLOR_RESET"
    fi

    echo -en " - simulation ($TIMEOUT timeout)... "
    timeout "$TIMEOUT" ./mfc.sh run "$f" -j $(nproc) -n 4 -b mpirun -t simulation --no-build > /dev/null 2>&1
    
    code="$?"
    if [[ "$code" != "0" ]] && [[ "$code" != "124" ]]; then
        echo -e "$RED""FAILED$COLOR_RESET"
        failed=(${failed[@]} "$f")
    else
        echo -e "$GREEN""PASSED$COLOR_RESET"
    fi

    index=$((index+1))

done

echo -e "\n${#failed[@]} samples failed:"
for value in "${failed[@]}"; do
    echo " - $value"
done

exit "${#failed[@]}"
