#!/usr/bin/env bash

FILENAME="output_185.txt"

COUNT=`cat $FILENAME | grep 'Present table' | wc -l`

declare -a FUNCTION_NAMES=(             \
"s_initialize_global_parameters_module" \
"s_initialize_riemann_solvers_module"   \
"s_initialize_rhs_module"               \
"s_initialize_time_steppers_module"     \
"s_initialize_weno_module"              \
)

declare -a TOTAL_MEMORY
TO_COLUMN_STR=""

for i in $(seq 1 $COUNT); do
    line_no=`cat $FILENAME | grep -n -m $i 'Present table' | tail -n 1 | tr ':' '\n' | head -n 1`
    line_count=`cat $FILENAME | wc -l`

    total_size=0

    j=$((line_no+1)); while [ $j -le $line_count ]; do
        line=`cat $FILENAME | head -n $j | tail -n 1`

        if [[ "$line" == *"Present table"* ]]; then
            break
        fi

        if [[ "$(echo "$line" | grep '^allocated' | wc -l)" == "1" ]]; then
            TOTAL_MEMORY[$((i-1))]=`expr ${TOTAL_MEMORY[$((i-1))]} + $(echo "$line" | tr " " "\n" | head -n 4 | tail -n 1 | sed s/size://)`
        fi

        j=$((j+1))
    done

    func="${FUNCTION_NAMES[$i-1]}"
    mem="${TOTAL_MEMORY[$i-1]}"
    delta="N/A"
    if [ "$i" -ne "1" ]; then
        mem_prev="${TOTAL_MEMORY[$i-2]}"
        delta=`expr $mem - $mem_prev`
    fi

    TO_COLUMN_STR="$TO_COLUMN_STR\n$func $(bc <<< "scale=5; $mem/1000000000") $(bc <<< "scale=5; $delta/1000000000")"
done

echo -e "$TO_COLUMN_STR" | column -t -N "Last Function Called,ACC Program Memory (GB),Delta (GB)" -o " | "

