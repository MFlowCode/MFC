#!/bin/bash

. ./mfc.sh load -c f -m g
./mfc.sh test --dry-run -j 8 --gpu
