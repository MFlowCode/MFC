#!/usr/bin/env bash

BUILD_DIR="./.build"

[ -d "$BUILD_DIR" ] && rm -rf "$BUILD_DIR"

mkdir -p "$BUILD_DIR"
cp *.f90 "$BUILD_DIR"
cp Makefile "$BUILD_DIR"

for fpp_f in *.fpp; do
	f90_f=$(echo "$fpp_f" | sed s/.fpp/.f90/g)
	fypp "$fpp_f" "$BUILD_DIR/$f90_f"
done

make -C "$BUILD_DIR" "$@"
