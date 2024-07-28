#!/bin/bash

./mfc.sh run /users/mrodri97/scratch/lungwave/case0/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pre -t pre_process -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case0/case.py -e batch -N 1 -n 4 -w 00:30:00 -# sim -t simulation -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case0/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pos -t post_process -c oscar

./mfc.sh run /users/mrodri97/scratch/lungwave/case1/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pre -t pre_process -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case1/case.py -e batch -N 1 -n 4 -w 00:30:00 -# sim -t simulation -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case1/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pos -t post_process -c oscar

./mfc.sh run /users/mrodri97/scratch/lungwave/case2/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pre -t pre_process -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case2/case.py -e batch -N 1 -n 4 -w 00:30:00 -# sim -t simulation -c oscar
./mfc.sh run /users/mrodri97/scratch/lungwave/case2/case.py -e batch -N 1 -n 4 -w 00:30:00 -# pos -t post_process -c oscar

