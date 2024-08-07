#!/bin/bash

./mfc.sh run ./examples/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run ./examples/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run ./examples/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

