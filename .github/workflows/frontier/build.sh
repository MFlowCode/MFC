#!/bin/bash

. ./mfc.sh load -c f -m g
./mfc.sh build -j 8 --gpu --sys-hdf5 --sys-fftw
