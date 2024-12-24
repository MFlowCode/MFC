#!/bin/bash

MAX_HEIGHT=400

set -e -x

examples_md="$1/docs/documentation/examples.md" 
rm "$examples_md" || true

echo -e "# Example Cases\n" > "$examples_md"

for casedir in $(find "$1/examples/" -mindepth 1 -maxdepth 1 -type d); do
    casename="$(basename "$casedir")"

    if [ -f "$casedir/README.md" ]; then
        sed -e "s/\.png/-$casename-example\.png/g" "$casedir/README.md" | sed 's/^#/##/g' | sed "s/MAX_HEIGHT/$MAX_HEIGHT/g" >> "$examples_md"
        echo '' >> "$examples_md"

        for png in $(find "$casedir" -maxdepth 1 -name '*.png'); do
            cp "$png" "$1/docs/documentation/$(basename "$png" | sed s/\.png//g)-$casename-example.png"
        done
    fi
done
