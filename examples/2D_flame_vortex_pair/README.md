# Premixed H2/Air Flame and a Vortex Pair (2D)

A stoichiometric H2/air flame (1 atm, 300 K, `h2o2.yaml`, mixture-averaged diffusion and conduction, viscous) in a 3.2 mm square is hit by a counter-rotating Lamb-Oseen vortex pair. The flame starts at x = 1.8 mm as the Cantera `FreeFlame` profile (`flame_profile.txt`, written by `make_profile.py`, flame speed S_L = 2.3315 m/s, thermal thickness 0.33 mm); the pair (circulation +-0.05 m^2/s, core radius 0.25 mm, centres 0.8 mm apart, y-periodic image sums, cyclostrophic pressure) starts at x = 0.9 mm and moves toward it. The boundaries are a slip wall at x = 0, the characteristic subsonic outflow at x = 3.2 mm (`bc_x%end = -8`, GRCBC relaxation to 1 atm) and periodic in y. The time step is fixed (2.632e-9 s) on every grid, so the first-order error of the split chemistry is the same on all grids. `case.py` writes the initial condition for the chosen grid to `IC/` at the first run (read with `hcid = 274`).

    ./mfc.sh run examples/2D_flame_vortex_pair/case.py -n 1 -t pre_process simulation post_process -- --N 128 --tend 1e-4
    python3 examples/2D_flame_vortex_pair/analyze_fv.py RUN64 RUN128 RUN256                # convergence table: copies of this directory, one per grid
    ./mfc.sh run examples/2D_flame_vortex_pair/case.py -n 1 -t pre_process simulation post_process -- --N 256 --ny 32 --gamma 0 --tend 1e-4
    python3 examples/2D_flame_vortex_pair/analyze_fv.py RUN --flat                          # flat-flame speed against S_L

The case writes its output next to `case.py`: copy the directory (keeping `fv_common.py` and `flame_profile.txt`) for each grid. `analyze_fv.py` reads the Silo output and prints, at 40, 70 and 100 us, the heat release Q, the H2O mass and the length of the T = 1200 K isoline per grid with the observed order and Richardson value; with `--flat` it prints the displacement speed of the vortex-free strip against S_L. `case.py` takes `--N`, `--ny`, `--gamma`, `--tend`, `--frames`, `--dt` and `--set KEY=VAL`.

Measured on `master` at 0cdfdf02a (fixed dt = 2.632e-9 s, to 70 us):

| Quantity | t [us] | 64 x 64 | 128 x 128 | 256 x 256 | earlier (64, 128, 256, 512) |
| :--- | ---: | ---: | ---: | ---: | :--- |
| Q [W/m] | 40 | 18674.5 | 18660.7 | 18632.6 | 18675.1, 18662.2, 18632.8, 18618.6 |
| | 70 | 23317.2 | 23367.4 | 23401.6 | 23332.1, 23367.6, 23403.3, 23411.9 |
| H2O mass [1e-7 kg/m] | 40 | 1.52564 | 1.52500 | 1.52421 | 1.525664, 1.525013, 1.524224, 1.523798 |
| | 70 | 1.61155 | 1.61153 | 1.61110 | 1.611683, 1.611551, 1.611127, 1.610895 |
| flame length [mm] | 40 | 3.5284 | 3.5239 | 3.5243 | 3.52840, 3.52388, 3.52434, 3.52454 |
| | 70 | 5.0647 | 5.0801 | 5.0845 | 5.06566, 5.08060, 5.08497, 5.08562 |

The 64 x 64 run to 100 us gives Q = 32008.3 W/m, H2O mass 1.73517e-7 kg/m and flame length 7.14806 mm (earlier 32013.4, 1.735351e-7, 7.14902). Every value agrees with the earlier study to 6e-4 or better. The 64, 128 and 256 runs are the 64 x 64 and 128 x 128 runs on the CPU (NVHPC 25.11) and the 256 x 256 run on an NVIDIA A100 (NVHPC 25.11), all in double precision.

Observed orders from 64, 128 and 256 are not in the asymptotic range (the heat release at 40 us is non-monotone in the grid: -1.0 at 40 us, +0.6 at 70 us; the flame length at 70 us gives 1.8), as in the earlier study, where the orders 1.0 (40 us) and 2.1 (70 us) came from the three finest grids (128, 256, 512). The 512 x 512 grid was not run here.

Flat flame (`--N 256 --ny 32 --gamma 0`, a 256 x 32 strip at the 256 x 256 spacing, 26.6 thousand steps to 70 us, NVIDIA A100, NVHPC 25.11): the displacement speed over 30 to 70 us (41 saves one microsecond apart, x_f fitted linearly, u_a and rho_a at T_up + 10 K) is S_d = 2.3393 m/s against the Cantera S_L = 2.3315 m/s, +0.33 % (earlier: +0.07 % over 30 to 90 us, +0.05 % in 3D).

The suite registers `2D -> Chemistry -> Flame Vortex Pair` (64 x 64, 100 steps). Results and limits: @ref validation_flame_vortex.
