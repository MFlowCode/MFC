#!/bin/bash

# Use gawk if available, otherwise fall back to awk
if command -v gawk > /dev/null; then
    AWK_CMD="gawk"
    IGNORECASE_BLOCK='BEGIN { IGNORECASE = 1 }'
else
    AWK_CMD="awk"
    IGNORECASE_BLOCK=''
    echo "Warning: gawk not found. Case-insensitive matching may not work as expected." >&2
fi

find . -type f \( -name "*.f90" -o -name "*.fpp" \) | while read file; do
    "$AWK_CMD" "
    $IGNORECASE_BLOCK
    /^[ \t]*((pure|elemental|impure)[ \t]+)*subroutine[ \t]+[a-zA-Z_][a-zA-Z0-9_]*[ \t]*\\(/ {
        in_sub = 1
        match(\$0, /subroutine[ \t]+([a-zA-Z_][a-zA-Z0-9_]*)/, arr)
        sub_name = arr[1]
        start_line = NR
        next
    }
    /^[ \t]*end[ \t]+subroutine[ \t]*([a-zA-Z_][a-zA-Z0-9_]*)?[ \t]*\$/ && in_sub {
        end_line = NR
        print (end_line - start_line + 1) \"\t\" FILENAME \": \" sub_name
        in_sub = 0
    }
    " "$file"
done | sort -nr | awk -F'\t' '{print $2 " : " $1 " lines"}' | head -20
