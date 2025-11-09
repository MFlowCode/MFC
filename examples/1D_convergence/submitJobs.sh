#!/usr/bin/env sh
Nx=(32 64 128 256 512 1024)
Order=(1 3 5)
ME=2 # Model equations = 2 for five-equation model
RS=2 # Riemann solver = 2 for HLLC

#ROOT_DIR=<WORKING DIRECTORY>
#MFC_DIR=<MFC ROOT DIR>
ROOT_DIR="/Users/benwilfong/Documents/software/MFC-Wilfong/examples/1D_convergence"
MFC_DIR="/Users/benwilfong/Documents/software/MFC-Wilfong"

for i in "${Nx[@]}"; do
    for j in "${Order[@]}"; do
        rm -rf N${i}_O${j}
        mkdir N${i}_O${j}
        cp case.py N${i}_O${j}/
    done
done

cd $MFC_DIR

for i in "${Nx[@]}"; do
    for j in "${Order[@]}"; do
        ./mfc.sh run $ROOT_DIR/N${i}_O${j}/case.py --case-optimization --no-debug -- --order $j -N $i --meqns $ME --rs $RS
    done
done

