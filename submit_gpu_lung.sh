#!/bin/bash

./mfc.sh run ./examples/2D_lungwave_horizontal/case.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run ./examples/2D_lungwave_horizontal/case.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

