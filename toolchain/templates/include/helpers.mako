<%def name="template_prologue()">
#>
#> The MFC prologue prints a summary of the running job and starts a timer.
#>

<%!
import os
%>\

. "${MFC_ROOTDIR}/toolchain/util.sh"

TABLE_FORMAT_LINE="| * %-14s $MAGENTA%-35s$COLOR_RESET * %-14s $MAGENTA%-35s$COLOR_RESET |\\n"
TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_TITLE_FORMAT="| %-105s |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time"   "$(date +%T)"                    "Start-date" "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition"    "${partition or 'N/A'}"          "Walltime"   "${walltime}")
$(printf "$TABLE_FORMAT_LINE" "Account"      "${account   or 'N/A'}"          "Nodes"      "${nodes}")
$(printf "$TABLE_FORMAT_LINE" "Job Name"     "${name}"                        "Engine"     "${engine}")
$(printf "$TABLE_FORMAT_LINE" "QoS"          "${quality_of_service or 'N/A'}" "Binary"     "${binary or 'N/A'}")
$(printf "$TABLE_FORMAT_LINE" "Queue System" "${qsystem.name}"                "Email"      "${email or 'N/A'}")
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

t_start=$(date +%s)
</%def>

<%def name="template_epilogue()">
#>
#> The MFC epilogue stops the timer and prints the execution summary. It also
#> performs some cleanup and housekeeping tasks before exiting.
#>

code=$?

t_stop="$(date +%s)"

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished ${name}:"
printf "$TABLE_FORMAT_LINE"  "Total-time:" "$(expr $t_stop - $t_start)s" "Exit Code:" "$code"
printf "$TABLE_FORMAT_LINE"  "End-time:"   "$(date +%T)"                 "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

exit $code
</%def>

<%def name="run_prologue(target)">
ok ":) Running$MAGENTA ${target.name}$COLOR_RESET:\n"

if [ '${target.name}' == 'simulation' ]; then
    export CRAY_ACC_MODULE='${target.get_staging_dirpath()}/simulation-wg256.lld.exe'
fi

cd "${os.path.dirname(input)}"

t_${target.name}_start=$(date +%s)
</%def>

<%def name="run_epilogue(target)">
code=$?

t_${target.name}_stop=$(date +%s)

if [ $code -ne 0 ]; then
    echo
    error ":( $MAGENTA${target.get_install_binpath()}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
    echo
    exit 1
fi

unset CRAY_ACC_MODULE

% if output_summary:

cd "${MFC_ROOTDIR}"

cat >>"${output_summary}" <<EOL
${target.name}: $(expr $t_${target.name}_stop - $t_${target.name}_start)
EOL

cd - > /dev/null

% endif
</%def>
