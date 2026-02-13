#!/bin/bash

# Fetch dependencies on login node (internet access).
# Source code is built on compute nodes via test.sh / bench.sh.

# Ignore SIGHUP to survive login node session drops
trap '' HUP

job_device=$1
job_interface=$2
build_opts=""
if [ "$job_device" = "gpu" ]; then
  build_opts+="--gpu"
  if [ "$job_interface" = "acc" ]; then
      build_opts+=" acc"
  elif [ "$job_interface" = "omp" ]; then
      build_opts+=" mp"
  fi
fi

. ./mfc.sh load -c f -m g

./mfc.sh build --deps-only -j 8 $build_opts
