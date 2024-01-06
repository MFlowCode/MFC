#>
#> The MFC epilogue stops the timer and prints the execution summary. It also
#> performs some cleanup and housekeeping tasks before exiting.
#>

code=$?

t_stop="$(date +%s)"

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished"    "$MAGENTA${name}$COLOR_RESET:"
printf "$TABLE_FORMAT_LINE"  "Total-time:" "$(expr $t_stop - $t_start)s" "Exit Code:" "$code"
printf "$TABLE_FORMAT_LINE"  "End-time:"   "$(date +%T)"                 "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

exit $code