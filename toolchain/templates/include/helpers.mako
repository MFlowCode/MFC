<%! import os %>

<%def name="template_prologue()">
% if os.name != 'nt':
    #>
    #> The MFC prologue prints a summary of the running job and starts a timer.
    #>

    . "${MFC_ROOT_DIR}/toolchain/util.sh"

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
% else:
    echo MFC case # ${name} @ ${input}:
    echo.
    echo Start Date: %DATE%
    echo Start Time: %TIME%
    echo Job Name:   ${name}
    echo Wall Time:  ${walltime}
    echo Nodes:      ${nodes}
    echo.
% endif
</%def>

<%def name="template_epilogue()">
% if os.name != 'nt':
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
% endif
</%def>

<%def name="run_prologue(target)">
% if os.name != 'nt':
    ok ":) Running$MAGENTA ${target.name}$COLOR_RESET:\n"

    cd '${os.path.dirname(input)}'

    t_${target.name}_start=$(python3 -c 'import time; print(time.time())')
% else:
    echo ^:) Running ${target.name}.
    echo.

    cd "${os.path.dirname(input)}"
% endif
</%def>

<%def name="run_epilogue(target)">
% if os.name != 'nt':
    code=$?

    t_${target.name}_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA${target.get_install_binpath(case)}$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    % if output_summary:

        cd '${MFC_ROOT_DIR}'

        cat >>'${output_summary}' <<EOL
${target.name}:
    exec:  $(echo "$t_${target.name}_stop - $t_${target.name}_start" | bc -l)
% if target == SIMULATION:
    grind: $(cat '${os.path.join(os.path.dirname(input), 'time_data.dat')}' | tail -n 1 | awk '{print $NF}')
% endif
EOL

        cd - > /dev/null

    % endif
% else:
    set code=%errorlevel%

    if %code% neq 0 (
        echo.
        echo ^:( ${target.get_install_binpath(case)} failed with exit code %code%.
        echo.

        exit /b %code%
    )
% endif
</%def>
