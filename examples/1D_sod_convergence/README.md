# 1D Riemann Problems Against the Exact Solution

Three ideal-gas ($\gamma = 1.4$) shock tubes on $0 \le x \le 1$ with the jump at $x = 0.5$ and outflow ends: Sod (default; $t = 0.2$),
Lax (`--problem lax`; $t = 0.14$) and Toro's test 2, the "123" double rarefaction (`--problem 123`; $t = 0.15$). `case.py` also
takes the grid (`-N`), the reconstruction (`--order`, `--muscl`, `--muscl-lim`, `--teno`), `--riemann hll|hllc|lax_friedrichs`,
`--cfl`, and `--igr 3|5`, which runs information geometric regularization on a 2D strip 8 cells high, periodic in $y$
(IGR has no 1D path), with elliptic smoothing of the initial jump.

`check_riemann.py` reads one or more finished runs and prints the L1 errors of density, velocity and pressure against the exact
cell averages (`riemann_exact.py`, Toro chapter 4), and the observed order between successive grids. It reads the states, the
interface, $\gamma$ and the time from the run's input files, so it also checks `examples/1D_sodshocktube` and
`examples/1D_laxshocktube`.

    ./mfc.sh run examples/1D_sod_convergence/case.py -t pre_process simulation -n 1 -- -N 512 --problem lax
    python3 examples/1D_sod_convergence/check_riemann.py <run 128> <run 256> <run 512> [--tol 1e-3]

The output is written next to `case.py`, so copy the directory per run. Results (WENO1-7, TENO5, MUSCL, HLL/HLLC/Lax-Friedrichs,
IGR, 128 to 2048 cells) and limits: @ref validation_shock_tubes. The suite runs the L1 self-convergence of this case for ten
reconstructions (`Convergence -> Sod`).
