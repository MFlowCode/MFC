# 1D Sod Shock Tube with AMR Against the Exact Riemann Solution

The Sod problem (rho, u, p = 1, 0, 1 left and 0.125, 0, 0.1 right of x = 0.5, gamma = 1.4) run to t = 0.2 with one or two refined 2:1 levels that regrid onto the shock and the contact (`amr_tag_eps = 0.05` on the normalized density jump; the rarefaction stays at the base resolution). `case.py` uses mapped WENO5, HLLC and SSP-RK3 at CFL 0.4 on the cells that set `dt` (the finest level lock-step, the base grid subcycled), outflow ends.

    ./mfc.sh run examples/1D_amr_sod/case.py -t pre_process simulation -n 1 -- -N 200 --levels 2
    python3 examples/1D_amr_sod/check_sod.py examples/1D_amr_sod

`case.py` takes `-N` (base cells), `--levels` (refined levels; 0 is a uniform grid), `--fine` (a uniform grid at the finest AMR spacing, `N * 2**levels` cells), `--subcycle`, `--buf` (`amr_buf`, default 8), `--regrid` (`amr_regrid_int`), `--tag-eps`, `--cfl` and `--t-end`. `check_sod.py` takes one or more run directories and prints, from the composite solution (the finest level covering each point, read from the restart files), the L1 density error against the exact solution's cell averages and the signed contact and shock position errors (where rho crosses the middle of its exact jump). The exact solution is computed in the checker (the pressure function of Toro, ch. 4).

Measured on `master` at 5838b641c, GNU CPU, double precision. L1(rho):

| Base N | uniform N | 1 level | 1 level, subcycled | uniform 2N | 2 levels | 2 levels, subcycled | uniform 4N |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 3.42e-3 | 1.90e-3 | 1.90e-3 | 1.89e-3 | 9.26e-4 | 9.73e-4 | 9.37e-4 |
| 200 | 1.89e-3 | 9.38e-4 | 9.39e-4 | 9.37e-4 | 5.01e-4 | 4.99e-4 | 4.50e-4 |
| 400 | 9.37e-4 | 4.45e-4 | 4.45e-4 | 4.50e-4 | 2.75e-4 | 2.73e-4 | 2.62e-4 |

With one level the error is the uniform 2N error within 1.1 %; with two levels it is within 11 % of the uniform 4N error, on 40 to 48 % of its cells. With two levels and `--buf 4` the run keeps the toolchain's `amr_snap = 2`, and the shock runs ahead by 3.3e-3 at every N (see the validation page).

The suite registers `AMR -> 1D -> Example -> Sod two levels` (`-N 32`, two levels, lock-step). Results and limits: @ref validation.
