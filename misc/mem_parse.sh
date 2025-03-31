#!/usr/bin/env bash

FILENAME="$1"

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

    CURRENT_VAR_IDX=1
    CURRENT_VARIABLE_COLS=""

    j=$((line_no+1)); while [ $j -le $line_count ]; do
        line=`cat $FILENAME | head -n $j | tail -n 1`

        if [[ "$line" == *"Present table"* ]]; then
            break
        fi

        if [[ "$(echo "$line" | grep '^allocated' | wc -l)" == "1" ]]; then
            TOTAL_MEMORY[$((i-1))]=`expr ${TOTAL_MEMORY[$((i-1))]} + $(echo "$line" | tr " " "\n" | head -n 4 | tail -n 1 | sed s/size://)`
        elif [[ "$(echo "$line" | grep '^host' | wc -l)" == "1" ]]; then
            var_name=`echo "$line" | tr " " "\n" | tail -n 1 | sed s/name://`
            var_size=`echo "$line" | tr " " "\n" | tail -n 4 | head -n 1 | sed s/size://`
            var_line=`echo "$line" | tr " " "\n" | tail -n 2 | head -n 1 | sed s/line://`

            CURRENT_VARIABLE_COLS="$CURRENT_VARIABLE_COLS\n$CURRENT_VAR_IDX $var_line $var_name $(bc <<< "scale=5; $var_size/1000000000") $(bc <<< "scale=5; $var_size/1000000")"

            CURRENT_VAR_IDX=$((CURRENT_VAR_IDX+1))
        fi

        j=$((j+1))
    done

    func="${FUNCTION_NAMES[$i-1]}"

    echo -e ""
    echo -e "---- Last Function Call: $func ----"
    echo -e ""
    echo -e "$CURRENT_VARIABLE_COLS" | column -t -N "#,Line #,Variable,Size (GB),Size (MB)" -o " | "
    echo -e ""

    mem="${TOTAL_MEMORY[$i-1]}"
    delta="N/A"
    if [ "$i" -ne "1" ]; then
        mem_prev="${TOTAL_MEMORY[$i-2]}"
        delta=`expr $mem - $mem_prev`
    fi

    TO_COLUMN_STR="$TO_COLUMN_STR\n$i $func $(bc <<< "scale=5; $mem/1000000000") $(bc <<< "scale=5; $delta/1000000000")"
done

echo -e "\n\n"
echo -e "$TO_COLUMN_STR" | column -t -N "#,Last Function Called,ACC Program Memory (GB),Delta (GB)" -o " | "
echo -e "\n"

