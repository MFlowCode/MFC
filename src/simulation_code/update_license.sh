#!/bin/bash

FILES="*.f90"

for f in $FILES; do
  echo "Processing $f file..."
  # cat "$f"
    tail -n +25 "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
