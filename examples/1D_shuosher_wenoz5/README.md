# Shu-Osher Problem (1D, WENO5-Z)

A Mach 3 shock runs into the density wave $\rho = 1 + 0.2\sin 5x$ (Shu & Osher 1989), here on $0 \le x \le 10$ with the shock at
$x = 1$, to $t = 1.8$. The other `1D_shuosher_*` examples run the same problem with WENO5-JS, WENO5-M, TENO5 and TENO7. `-N` sets
the cell count (default 1001; the step count scales with it).

`check_shuosher.py` compares one or more runs of any of these examples with a fine reference run (`-N 16016`), averaged onto each
run's cells, and prints the L1 density error over the domain and over the shocked entropy waves:

    ./mfc.sh run examples/1D_shuosher_wenoz5/case.py -t pre_process simulation -n 1 -- -N 16016
    python3 examples/1D_shuosher_wenoz5/check_shuosher.py <reference run> <run> [<run> ...]

Results (all five schemes at 1001 cells, WENO5-Z from 143 to 8008 cells, CPU and A100) and limits: @ref validation_shock_tubes.
