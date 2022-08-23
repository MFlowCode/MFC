#!/usr/bin/env bash

RED="\u001b[31m"
CYAN="\u001b[36m"
GREEN="\u001b[32m"
YELLOW="\u001b[33m"
MAGENTA="\u001b[35m"
COLOR_RESET="\033[m"

TIMEOUT="10s"

./mfc.sh build -j $(nproc) -t pre_process simulation

index=1
failed=()

files=samples/*/*.py

echo "$files"

for f in samples/*/*.py; do

    echo -en "($index) Running $MAGENTA$f$COLOR_RESET - ($TIMEOUT timeout): "
    timeout "$TIMEOUT" ./mfc.sh run "$f" -j $(nproc) -n 4 -b mpirun -t pre_process simulation > /dev/null 2>&1
    code="$?"
    if [[ "$code" != "0" ]] && [[ "$code" != "124" ]]; then
        echo -e "[$RED""Failed$COLOR_RESET]"
        failed=(${failed[@]} "$f")
    else
        echo -e "[$GREEN""Passed$COLOR_RESET]"
    fi

    index=$((index+1))

done

echo -e "\n${#failed[@]} samples failed:"
for value in "${failed[@]}"; do
    echo " - $value"
done

exit "${#failed[@]}"
