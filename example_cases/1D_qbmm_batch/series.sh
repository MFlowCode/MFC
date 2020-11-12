#!/bin/bash

# Get current directory
curdir=$(pwd)   

cfl=0.01
nx=400
npp=24
nfiles=50

Re=10
#R0-quad: 1 = Simpson; 2 = Wheeler
r0type=1

sigrs=( 0.2 )
sigvs=( 0.2 )
sig0s=( 0.2 )
# sigrs=( 0.1 0.2 0.3 )
# sigvs=( 0.1 0.2 0.3 )
# sig0s=( 0.1 0.2 0.3 )
nbs=( 3 )
# nbs=( 1 5 9 21 41 71 )
vfcs=( 0.0001 )
pfacs=( 5 )

for sigr  in "${sigrs[@]}"; do
for sigv  in "${sigvs[@]}"; do
for sig0  in "${sig0s[@]}"; do
for vfrac in "${vfcs[@]}";  do
for pfac  in "${pfacs[@]}"; do
for nb    in "${nbs[@]}";   do
    workdir='scr1d-'$vfrac'vf0-'$pfac'pf-'$nb'nr0-'$r0type'r0typ-'$sigr'sigr-'$sigv'sigv-'$sig0'sig0-'$Re'Re'
    # workdir='bubscr1d-qbmm-'$vfrac'vf0-'$pfac'pf-'$nb'nr0-'$r0type'r0typ-dtsm'
    echo Working directory: $workdir

    if [ -d $workdir ]; then
        echo $workdir already exists!
        read -p "Delete and continue? " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Nn]$ ]]; then
            echo Exit...
            exit
        fi
        rm -rf $workdir
    fi


    cp -r dummy/ $workdir
    
    poly='T'
    if [ $nb -eq 1 ]; then
        poly='F'
    fi

    cat ./dummy/input_qbmm.dummy | \
        sed \
        -e s/XXNAME/$workdir/g \
        -e s/XXNP/$npp/g \
        -e s/XXNX/$nx/g \
        -e s/XXNFILES/$nfiles/g \
        -e s/XXCFL/$cfl/g \
        -e s/XXPFAC/$pfac/g \
        -e s/XXRE/$Re/g \
        -e s/XXSIGR/$sigr/g \
        -e s/XXSIGV/$sigv/g \
        -e s/XXSIG0/$sig0/g \
        -e s/XXPOLY/$poly/g \
        -e s/XXR0TYPE/$r0type/g \
        -e s/XXNB/$nb/g \
        -e s/XXVF/$vfrac/g \
        > ./$workdir/input.py

    cd $workdir
        rm -f *.dummy
        mkdir out D
        chmod +x input.py
        ./input.py pre_process  > out/pre_process.out
        sleep 1
        ./input.py simulation  > out/simulation.out
        sleep 1
    cd $curdir
done
done
done
done
done
done
