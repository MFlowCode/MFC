#!/bin/bash

./mfc.sh build -t pre_process simulation --case-optimization -i examples/scaling/case.py \
               -j 8 --gpu --mpi --no-debug -- -s strong -m 512