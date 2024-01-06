#>
#> The MFC prologue prints a summary of the running job and starts a timer.
#>

<%!
import os
%>\

. "${rootdir}/toolchain/util.sh"

TABLE_FORMAT_LINE="| * %-14s $MAGENTA%-35s$COLOR_RESET * %-14s $MAGENTA%-35s$COLOR_RESET |\\n"
TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_TITLE_FORMAT="| %-105s |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time"   "$(date +%T)"    "Start-date" "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition"    "${partition}"   "Walltime"   "${walltime}")
$(printf "$TABLE_FORMAT_LINE" "Account"      "${account}"     "Nodes"      "${nodes}")
$(printf "$TABLE_FORMAT_LINE" "Job Name"     "${name}"        "Engine"      "${engine}")
$(printf "$TABLE_FORMAT_LINE" "Queue System" "{qsystem.name}" "Email"      "${email}")
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "MFC case # ${name} @ ${input}:"
printf "$TABLE_HEADER"
printf "$TABLE_CONTENT\\n"
printf "$TABLE_FOOTER\\n"

% for key, value in env.items():
    export ${key}='${value}'
% endfor

cd "${os.path.dirname(input)}"

t_start=$(date +%s)