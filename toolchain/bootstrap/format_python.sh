#!/bin/bash

. toolchain/util.sh

echo "> $1"

black -l 200 -x "$1" &>/dev/null 2>&1

