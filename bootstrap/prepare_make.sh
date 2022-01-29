#!/usr/bin/env bash

BUILD_DIR="./.build"

mkdir -p "$BUILD_DIR"
cp -u *.f90    "$BUILD_DIR"
cp -u Makefile "$BUILD_DIR"

for fpp_f in "$(find . -name "*.fpp" -not -path "$BUILD_DIR/*")"; do
	if [ -f "$fpp_f" ]; then
		f90_f=$(echo "$fpp_f" | sed s/.fpp/.f90/g)
		fypp "$fpp_f" "$BUILD_DIR/$f90_f"
		if (($?)); then
			echo "[bootstrap/prepare_make.sh] Failed to execute fypp with:"
			echo " - Base:        '$(pwd)'"
			echo " - Source:      '$fpp_f'"
			echo " - Destination: '$BUILD_DIR/$f90_f'"
			exit 1
		fi
	fi
done

make -C "$BUILD_DIR" "$@"
exit $?
