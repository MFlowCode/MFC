#!/usr/bin/env bash

if [ ! -f "$(pwd)/toolchain/util.sh" ]; then
    echo ""
    echo "[samples/check_samples.sh] Please execute this script from the top-level MFC directory (currently in: $(pwd))."

    exit 1
fi

source $(pwd)/toolchain/util.sh

TIMEOUT="60s"

if ! ./mfc.sh build -j $(nproc) -t pre_process simulation; then
    exit 1
fi

echo ""
echo -e "Simulation's successful timeout is $TIMEOUT."
echo ""

index=1
failed=()

files=examples/*/*.py

nFiles=$(find examples/*/*.py | wc -l)

for f in examples/*/*.py; do

    echo -e "($index/$nFiles) Running $MAGENTA$f$COLOR_RESET:"
    echo -en " - pre_process "
    if ! ./mfc.sh run "$f" -j $(nproc) -n 4 -t pre_process --no-build > /dev/null 2>&1; then
        echo -e "$RED""CRASHED$COLOR_RESET"

        failed=(${failed[@]} "$f")
        index=$((index+1))

        continue
    else
        echo -e "$GREEN""RAN$COLOR_RESET"
    fi

    echo -en " - simulation  "
    timeout "$TIMEOUT" ./mfc.sh run "$f" -j $(nproc) -n 4 -t simulation --no-build > /dev/null 2>&1

    code="$?"
    if [[ "$code" != "0" ]] && [[ "$code" != "124" ]]; then
        echo -e "$RED""CRASHED$COLOR_RESET"
        failed=(${failed[@]} "$f")
    else
        echo -e "$GREEN""RAN$COLOR_RESET"
    fi

    index=$((index+1))

done

echo -e "\n${#failed[@]} samples failed:"
for value in "${failed[@]}"; do
    echo " - $value"
done

exit "${#failed[@]}"

