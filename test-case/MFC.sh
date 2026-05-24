
#!/usr/bin/env bash




    #>
    #> The MFC prologue prints a summary of the running job and starts a timer.
    #>

    . "/Users/fahadnabid/Documents/default/MFC/toolchain/util.sh"

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
    printf "$TABLE_TITLE_FORMAT" "MFC case # MFC @ /Users/fahadnabid/Documents/default/MFC/test-case/case2.py:"
    printf "$TABLE_HEADER"
    printf "$TABLE_CONTENT\\n"
    printf "$TABLE_FOOTER\\n"


    t_start=$(date +%s)




    warn "This is the$MAGENTA default$COLOR_RESET template."
    warn "It is not intended to support all systems and execution engines."
    warn "Consider using a different template via the $MAGENTA--computer$COLOR_RESET option if you encounter problems."

        # Find a suitable MPI launcher and store it in the variable "binary".
        for binary in  jsrun srun mpirun mpiexec; do
            if command -v $binary > /dev/null; then
                break
            fi
        done

        if ! command -v $binary > /dev/null; then
            error ":( Could not find a suitable MPI launcher.\n"
            exit 1
        else
            ok ":) Selected MPI launcher $MAGENTA$binary$COLOR_RESET. Use$MAGENTA --binary$COLOR_RESET to override."
        fi

        
    ok ":) Running$MAGENTA syscheck$COLOR_RESET:\n"

    cd '/Users/fahadnabid/Documents/default/MFC/test-case'

    t_syscheck_start=$(python3 -c 'import time; print(time.time())')


            if [ "$binary" == "jsrun" ]; then
                (set -x;                        jsrun --nrs          6                         --cpu_per_rs   1                                               --gpu_per_rs   0                              --tasks_per_rs 1                                               "/Users/fahadnabid/Documents/default/MFC/build/install/ccca43eb0a/bin/syscheck")
            elif [ "$binary" == "srun" ]; then
                (set -x;                                                 srun --ntasks 6                               "/Users/fahadnabid/Documents/default/MFC/build/install/ccca43eb0a/bin/syscheck")
            elif [ "$binary" == "mpirun" ]; then
                (set -x;                          $binary -np 6                                        "/Users/fahadnabid/Documents/default/MFC/build/install/ccca43eb0a/bin/syscheck")
            elif [ "$binary" == "mpiexec" ]; then
                (set -x;                                                    $binary --ntasks 6                                   "/Users/fahadnabid/Documents/default/MFC/build/install/ccca43eb0a/bin/syscheck")
            fi

        
    code=$?

    t_syscheck_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/Users/fahadnabid/Documents/default/MFC/build/install/ccca43eb0a/bin/syscheck$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi




        echo
        
    ok ":) Running$MAGENTA pre_process$COLOR_RESET:\n"

    cd '/Users/fahadnabid/Documents/default/MFC/test-case'

    t_pre_process_start=$(python3 -c 'import time; print(time.time())')


            if [ "$binary" == "jsrun" ]; then
                (set -x;                        jsrun --nrs          6                         --cpu_per_rs   1                                               --gpu_per_rs   0                              --tasks_per_rs 1                                               "/Users/fahadnabid/Documents/default/MFC/build/install/77ad0e446d/bin/pre_process")
            elif [ "$binary" == "srun" ]; then
                (set -x;                                                 srun --ntasks 6                               "/Users/fahadnabid/Documents/default/MFC/build/install/77ad0e446d/bin/pre_process")
            elif [ "$binary" == "mpirun" ]; then
                (set -x;                          $binary -np 6                                        "/Users/fahadnabid/Documents/default/MFC/build/install/77ad0e446d/bin/pre_process")
            elif [ "$binary" == "mpiexec" ]; then
                (set -x;                                                    $binary --ntasks 6                                   "/Users/fahadnabid/Documents/default/MFC/build/install/77ad0e446d/bin/pre_process")
            fi

        
    code=$?

    t_pre_process_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/Users/fahadnabid/Documents/default/MFC/build/install/77ad0e446d/bin/pre_process$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi




        echo
        
    ok ":) Running$MAGENTA simulation$COLOR_RESET:\n"

    cd '/Users/fahadnabid/Documents/default/MFC/test-case'

    t_simulation_start=$(python3 -c 'import time; print(time.time())')


            if [ "$binary" == "jsrun" ]; then
                (set -x;                        jsrun --nrs          6                         --cpu_per_rs   1                                               --gpu_per_rs   0                              --tasks_per_rs 1                                               "/Users/fahadnabid/Documents/default/MFC/build/install/503a859d02/bin/simulation")
            elif [ "$binary" == "srun" ]; then
                (set -x;                                                 srun --ntasks 6                               "/Users/fahadnabid/Documents/default/MFC/build/install/503a859d02/bin/simulation")
            elif [ "$binary" == "mpirun" ]; then
                (set -x;                          $binary -np 6                                        "/Users/fahadnabid/Documents/default/MFC/build/install/503a859d02/bin/simulation")
            elif [ "$binary" == "mpiexec" ]; then
                (set -x;                                                    $binary --ntasks 6                                   "/Users/fahadnabid/Documents/default/MFC/build/install/503a859d02/bin/simulation")
            fi

        
    code=$?

    t_simulation_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/Users/fahadnabid/Documents/default/MFC/build/install/503a859d02/bin/simulation$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi




        echo
        
    ok ":) Running$MAGENTA post_process$COLOR_RESET:\n"

    cd '/Users/fahadnabid/Documents/default/MFC/test-case'

    t_post_process_start=$(python3 -c 'import time; print(time.time())')


            if [ "$binary" == "jsrun" ]; then
                (set -x;                        jsrun --nrs          6                         --cpu_per_rs   1                                               --gpu_per_rs   0                              --tasks_per_rs 1                                               "/Users/fahadnabid/Documents/default/MFC/build/install/c084270971/bin/post_process")
            elif [ "$binary" == "srun" ]; then
                (set -x;                                                 srun --ntasks 6                               "/Users/fahadnabid/Documents/default/MFC/build/install/c084270971/bin/post_process")
            elif [ "$binary" == "mpirun" ]; then
                (set -x;                          $binary -np 6                                        "/Users/fahadnabid/Documents/default/MFC/build/install/c084270971/bin/post_process")
            elif [ "$binary" == "mpiexec" ]; then
                (set -x;                                                    $binary --ntasks 6                                   "/Users/fahadnabid/Documents/default/MFC/build/install/c084270971/bin/post_process")
            fi

        
    code=$?

    t_post_process_stop=$(python3 -c 'import time; print(time.time())')

    if [ $code -eq 22 ]; then
        echo
        error "$YELLOW CASE FILE ERROR$COLOR_RESET > $YELLOW Case file has prohibited conditions as stated above.$COLOR_RESET"
    fi

    if [ $code -ne 0 ]; then
        echo
        error ":( $MAGENTA/Users/fahadnabid/Documents/default/MFC/build/install/c084270971/bin/post_process$COLOR_RESET failed with exit code $MAGENTA$code$COLOR_RESET."
        echo
        exit 1
    fi

    python3 "/Users/fahadnabid/Documents/default/MFC/toolchain/templates/include/generate_silo_series.py" '/Users/fahadnabid/Documents/default/MFC/test-case'



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

