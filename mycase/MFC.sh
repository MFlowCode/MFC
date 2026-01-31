#!/usr/bin/env bash





    #>
    #> The MFC prologue prints a summary of the running job and starts a timer.
    #>

    . "/storage/scratch1/6/sbryngelson3/MFC-spencer/toolchain/util.sh"

    TABLE_FORMAT_LINE="| * %-14s $MAGENTA%-35s$COLOR_RESET * %-14s $MAGENTA%-35s$COLOR_RESET |\\n"
    TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
    TABLE_TITLE_FORMAT="| %-105s |\\n"
    TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time"   "$(date +%T)"                    "Start-date" "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition"    "N/A"          "Walltime"   "01:00:00")
$(printf "$TABLE_FORMAT_LINE" "Account"      "N/A"          "Nodes"      "1")
$(printf "$TABLE_FORMAT_LINE" "Job Name"     "MFC"                        "Engine"     "interactive")
$(printf "$TABLE_FORMAT_LINE" "QoS"          "N/A" "Binary"     "N/A")
$(printf "$TABLE_FORMAT_LINE" "Queue System" "Interactive"                "Email"      "N/A")
END
)

    printf "$TABLE_HEADER"
    printf "$TABLE_TITLE_FORMAT" "MFC case # MFC @ /storage/scratch1/6/sbryngelson3/MFC-spencer/mycase/case.py:"
    printf "$TABLE_HEADER"
    printf "$TABLE_CONTENT\\n"
    printf "$TABLE_FOOTER\\n"


    t_start=$(date +%s)


ok ":) Loading modules:\n"
cd "/storage/scratch1/6/sbryngelson3/MFC-spencer"
. ./mfc.sh load -c p -m c
cd - > /dev/null
echo

    
    ok ":) Running$MAGENTA syscheck$COLOR_RESET:\n"

    cd '/storage/scratch1/6/sbryngelson3/MFC-spencer/mycase'

    t_syscheck_start=$(python3 -c 'import time; print(time.time())')


        (set -x;                 mpirun -np 2                               --bind-to none                                            "/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/36dd4af3e9/bin/syscheck")

    
    code=$?

    t_syscheck_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/36dd4af3e9/bin/syscheck$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi



    echo
    
    ok ":) Running$MAGENTA pre_process$COLOR_RESET:\n"

    cd '/storage/scratch1/6/sbryngelson3/MFC-spencer/mycase'

    t_pre_process_start=$(python3 -c 'import time; print(time.time())')


        (set -x;                 mpirun -np 2                               --bind-to none                                            "/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/52b60342be/bin/pre_process")

    
    code=$?

    t_pre_process_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/52b60342be/bin/pre_process$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi



    echo
    
    ok ":) Running$MAGENTA simulation$COLOR_RESET:\n"

    cd '/storage/scratch1/6/sbryngelson3/MFC-spencer/mycase'

    t_simulation_start=$(python3 -c 'import time; print(time.time())')


        (set -x;                 mpirun -np 2                               --bind-to none                                            "/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/51b5b52d6c/bin/simulation")

    
    code=$?

    t_simulation_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/storage/scratch1/6/sbryngelson3/MFC-spencer/build/install/51b5b52d6c/bin/simulation$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi



    echo


    #>
    #> The MFC epilogue stops the timer and prints the execution summary. It also
    #> performs some cleanup and housekeeping tasks before exiting.
    #>

    code=$?

    t_stop="$(date +%s)"

    printf "$TABLE_HEADER"
    printf "$TABLE_TITLE_FORMAT" "Finished MFC:"
    printf "$TABLE_FORMAT_LINE"  "Total-time:" "$(expr $t_stop - $t_start)s" "Exit Code:" "$code"
    printf "$TABLE_FORMAT_LINE"  "End-time:"   "$(date +%T)"                 "End-date:"  "$(date +%T)"
    printf "$TABLE_FOOTER"

    exit $code

