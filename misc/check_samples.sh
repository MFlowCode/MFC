#!/usr/bin/env bash

failed=()

for d in samples/*/; do
    echo "Running $d.."
    # 2D_shearlayer a large amount of timesteps so one might want to skip it
    #if [[ ! "$d" == *"2D_shearlayer"* ]]; then
        mkdir -p "$d""D"
        ./mfc.sh run "$d"case.py -j $(nproc) -c 4 -b mpirun -t pre_process simulation
        if [[ "$?" != "0" ]]; then
            failed=(${failed[@]} "$d""case.py")
        fi
    #fi
done

echo "Samples that failed:"
for value in "${failed[@]}"; do
    echo " - $value"
done

exit "${#failed[@]}"
