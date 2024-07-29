#!/bin/bash

for file in $(find src -type f | grep -Ev 'autogen' | grep -E '\.fpp$'); do
    echo "$file"
    mv "$file" "$(echo "$file" | sed s/\.fpp/\.fypp/)"
done
