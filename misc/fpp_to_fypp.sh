#!/bin/bash

for file in $(find src -type f | grep -Ev 'autogen' | grep -E '\.fypp$'); do
    echo "$file"
    mv "$file" "$(echo "$file" | sed s/\.fypp/\.fypp/)"
done
