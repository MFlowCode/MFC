#!/bin/bash

### Note: Designed to be run from MFC root directory, not tests/

rm -rf ./tests/*/D ./tests/*/*.inp ./tests/*/p_all ./tests/*/*.out

mytests=( 1d_bubscreen kapila_shocktube sod_shocktube vacuum_generation )

ntest=${#mytests[@]}

npass=0
nfail=0
i=1

echo -----------------------------------------------

cd tests
for mytest in "${mytests[@]}"; do
    cd $mytest
        
        #Run test case
        ./input.py pre_process > pre_process.out
        ./input.py  simulation > simulation.out
        
        cd check
            check_file=$(echo *)
        cd ..
        
        #Check that the files are the same
        rm -f diff.out
        diff check/$check_file D/$check_file > diff.out

        mytest="Test $i of $ntest: $mytest"
        #Print if not
        if [ -s diff.txt ]; then
            echo $mytest: Test failed! Output files are different.
            ((++nfail))
        elif [ ! -f "D/$check_file" ]; then
            echo $mytest: Test failed! Output file was not found.
            ((++nfail))
        else
            echo $mytest: Test passed!
            ((++npass))
        fi

        ((++i))

    cd ..
done

echo -----------------------------------------------
echo ---- $nfail Tests failed
echo ---- $npass Tests passed
