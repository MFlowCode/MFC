#!/bin/bash

gpus=`rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' '`
ngpus=`echo "$gpus" | tr -d '[:space:]' | wc -c`

./mfc.sh test --max-attempts 3 -j $ngpus --sys-hdf5 --sys-fftw -- -c frontier

