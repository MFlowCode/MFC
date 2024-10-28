#!/bin/bash

. ./mfc.sh load -c f -m g
./mfc.sh build strict -j 8 --gpu
