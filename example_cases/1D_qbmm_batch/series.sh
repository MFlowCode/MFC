#!/bin/bash

# Get current directory
curdir=$(pwd)   

# Volume fractions to go through
vfcs=( 0.0001 0.001 )

for vfrac in "${vfcs[@]}"; do
    workdir='bubscr1d-qbmm-'$vfrac'vf0'
    echo Working directory: $workdir

    # Copy dummy directory to $workdir
    cp -r dummy/ $workdir
    
    # Replace input_qbmm.dummy with input.py, substituting XXVF for $vfrac
    cat ./dummy/input_qbmm.dummy | \
        sed -e s/XXNAME/$workdir/g -e s/XXVF/$vfrac/g > ./$workdir/input.py
    cd $workdir
        rm -f *.dummy
        chmod +x input.py
        ./input.py pre_process  > out/pre_process.out
        sleep 1
        ./input.py simulation  > out/simulation.out
        sleep 1
    cd $curdir
done
