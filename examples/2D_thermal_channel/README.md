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
cold run it also prints the pair value from half-differences, which cancels the first-order property asymmetry.
The exact value includes axial conduction at the run's Peclet number ($\mathrm{Pe}_H = 53.3$ gives 7.5417).

    python3 examples/2D_thermal_channel/check_nu.py <hot run> [<cold run>] [--tol 0.01]

Measured with NVHPC 25.11 (CPU, one rank), 32 cells across the gap, 25 along $x$, $T_m = 300$ K, run time 0.3 $H^2/\alpha$
(14070 steps; the state is developed by about 0.1). The source removes the viscous dissipation from the energy equation, which
the developed state excludes:

| run | $T_m$ hold | $\beta H$ (exact 0.1412) | $\mathrm{Nu}_\mathrm{src}$ | $\mathrm{Nu}_\mathrm{grad}$ |
|---|---:|---:|---:|---:|
| $T_s$ = 330 K | 2e-8 | 0.1432 | 7.375 (-2.2 %) | 7.342 (-2.6 %) |
| $T_s$ = 270 K | 6e-9 | 0.1397 | 7.746 (+2.7 %) | 7.771 (+3.0 %) |
| pair | | | 7.553 (+0.16 %) | 7.549 (+0.09 %) |

A single run is biased by the property change across the gap, with opposite signs for the hot and the cold wall; the pair cancels
it. The source and the gradient estimates of each run agree within 0.4 %, so the source removes the heat that crosses the walls.
At Mach 0.3 the single runs stay within 3.2 % (pair 7.556 and 7.527); before the dissipation was removed they were off by 14 % and
29 %, and at Mach 0.1 the two estimates of a run differed by 5 %.

## Run

    ./mfc.sh run examples/2D_thermal_channel/case.py -t pre_process simulation -n 1 -- --ny 32 --nx 25 --tw 330 --tend 0.3

(copy the directory for the second run with `--tw 270`, since the output is written next to `case.py`). `case.py` takes `--ny`, `--nx`, `--tw`, `--t0`,
`--mach`, `--re`, `--tend` (in $H^2/\alpha$) and `--saves`. The suite registers the first 100 steps of the 330 K case, with non-mapped WENO3 (the golden of the default WENO5 differs by a few percent
between compilers), as `2D -> Chemistry -> Thermal Channel -> Wall Temperature`.

## Limits

Only the plane channel is verified. In periodic arrays of particles $\beta$ is not converged and the hot and cold runs disagree
(see the note under `mean_temp_wall` in @ref case). The streamwise velocity in the first four cells next to each wall alternates about the
Poiseuille profile by up to 12 % of the bulk velocity (cause not verified; the Nusselt numbers are not affected at the 0.2 % level). See @ref validation.
