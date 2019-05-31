rm -rf ./tests/*/D ./tests/*/*.inp ./tests/*/p_all ./tests/*/*.out

mytests=( 1d_bubscreen kapila_shocktube sod_shocktube vacuum_generation )

cd tests
for mytest in "${mytests[@]}"; do
    cd $mytest
        python ./input.py pre_process > pre_process.out
        python ./input.py  simulation > simulation.out
        cd check
            check_file=$(echo *)
        cd ..
        rm -f diff.out
        diff check/$check_file D/$check_file > diff.out
        if [ -s diff.txt ]; then
            echo $mytest: Check failed!
        else
            echo $mytest: Check passed!
        fi
    cd ..
done
