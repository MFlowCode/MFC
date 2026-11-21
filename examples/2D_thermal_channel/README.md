# 2D Thermally Developed Plane Channel

Laminar flow between two isothermal no-slip walls at $T_s$, periodic in $x$ and driven at a fixed bulk mass flux
(`mean_mflux`, `mean_dpdx`). The gas is air (`h2o2.yaml`, reactions off, mixture-averaged transport) at 1 atm, bulk Mach number 0.1,
$\mathrm{Re}_{D_h} = 150$ ($D_h = 2H$). The mean-temperature source with `mean_temp_wall = T_s` replaces the
heat that the walls remove or add by the thermally developed form $T = T_s + \psi(y)\,e^{-\beta x}$ with $\psi$ periodic, and the controller
holds the mixing-cup temperature $T_m$ at `mean_temp`. The developed Nusselt number of constant-property laminar flow between
parallel plates at constant wall temperature is exact (7.541, Shah & London 1978), so the case tests the source, its decay rate
$\beta$, and the wall heat transfer together.

`check_nu.py` reads the finished run(s) and prints, for each: the hold error of $T_m$, $\beta H$ against the exact eigenvalue,
and the Nusselt number two ways, $\mathrm{Nu}_\mathrm{src}$ from the source total (equal to the wall heat at steady state) and
$\mathrm{Nu}_\mathrm{grad}$ from the wall temperature gradient of the last save, which does not use the source. With a hot and a
cold run it also prints the pair value from half-differences, which cancels viscous heating and the first-order property asymmetry.
The exact value includes axial conduction at the run's Peclet number ($\mathrm{Pe}_H = 53.3$ gives 7.5417).

    python3 examples/2D_thermal_channel/check_nu.py <hot run> [<cold run>] [--tol 0.01]

Measured with GNU 13.3 (CPU, one rank), 32 cells across the gap, 25 along $x$, $T_m = 300$ K, run time 0.3 $H^2/\alpha$
(14070 steps; the state is developed by about 0.1):

| run | $T_m$ hold | $\beta H$ (exact 0.1412) | $\mathrm{Nu}_\mathrm{src}$ | $\mathrm{Nu}_\mathrm{grad}$ |
|---|---:|---:|---:|---:|
| $T_s$ = 330 K | 2e-8 | 0.1457 | 7.507 (-0.46 %) | 7.141 (-5.3 %) |
| $T_s$ = 270 K | 4e-9 | 0.1371 | 7.600 (+0.77 %) | 7.998 (+6.1 %) |
| pair | | 0.1414 (mean) | 7.552 (+0.13 %) | 7.554 (+0.16 %) |

A single run is biased by viscous heating and by the property change across the gap, with opposite signs for the hot and the cold
wall; the pair cancels both. With 64 cells across the gap (run time 0.15 $H^2/\alpha$) both pair values are 7.546 (+0.05 %). The source and the gradient estimates agree within 0.1 % after pairing, so the
source removes the heat that crosses the walls.

## Run

    ./mfc.sh run examples/2D_thermal_channel/case.py -t pre_process simulation -n 1 -- --ny 32 --nx 25 --tw 330 --tend 0.3

(copy the directory for the second run with `--tw 270`, since the output is written next to `case.py`). `case.py` takes `--ny`, `--nx`, `--tw`, `--t0`,
`--mach`, `--re`, `--tend` (in $H^2/\alpha$) and `--saves`. The suite registers the first 100 steps of the 330 K case, with non-mapped WENO3 (the golden of the default WENO5 differs by a few percent
between compilers), as `2D -> Chemistry -> Thermal Channel -> Wall Temperature`.

## Limits

Only the plane channel is verified. In periodic arrays of particles $\beta$ is not converged and the hot and cold runs disagree
(see the note under `mean_temp_wall` in @ref case). The streamwise velocity in the first four cells next to each wall alternates about the
Poiseuille profile by up to 12 % of the bulk velocity (cause not verified; the Nusselt numbers are not affected at the 0.2 % level). See @ref validation.
