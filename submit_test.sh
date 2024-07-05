#!/bin/bash

./mfc.sh run ./tests/AED93D34/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
./mfc.sh run ./tests/AED93D34/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta

#./mfc.sh run ./tests/DA8AF07E/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./tests/DA8AF07E/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./tests/6F296065/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./tests/6F296065/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./tests/D3C860B9/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./tests/D3C860B9/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta

#./mfc.sh run ./tests/DA8AF07E/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./tests/DA8AF07E/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./tests/18431ACB/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./tests/18431ACB/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t simulation -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta


