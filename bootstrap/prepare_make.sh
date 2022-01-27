#!/usr/bin/env bash

BUILD_DIR="./.build"

mkdir -p "$BUILD_DIR"
cp -u *.f90 "$BUILD_DIR"
cp -u Makefile "$BUILD_DIR"

for fpp_f in *.fpp; do
	f90_f=$(echo "$fpp_f" | sed s/.fpp/.f90/g)
	fypp "$fpp_f" "$BUILD_DIR/$f90_f"
done

make -C "$BUILD_DIR" "$@"
