# Water Disk Advected in Air

A water disk (stiffened gas, pi_inf 1e9 Pa, 1000 kg/m^3) in air (1 kg/m^3), both at 1e5 Pa and u = v = 100 m/s, crosses a
periodic unit box once along the diagonal with the 6-equation model. The exact solution keeps p and u uniform and returns the disk
to its initial field. `check_advection.py` prints the largest pressure and velocity disturbance over all saves, the L1 error of
the volume fraction per unit interface length and the width of the mixed zone, for one or several grids.

    ./mfc.sh run examples/2D_advection/case.py -t pre_process simulation post_process -n 1 -- -N 100
    python3 examples/2D_advection/check_advection.py examples/2D_advection --case-args "-N 100"

`case.py` takes `-N` (cells per direction, default 100), `--periods`, `--cfl` (0.4), `--saves` (20) and `--eps` (trace volume
fraction, 1e-8). Measured with GNU 13.3 (CPU, one rank):

| N | max dp/p0 | max du/U | L1 [cells] | mixed-zone width, initial / final [cells] |
|---:|---:|---:|---:|---:|
| 50 | 2.5e-3 | 1.7e-4 | 0.066 | 9.17 / 9.25 |
| 100 | 3.1e-6 | 1.9e-9 | 0.095 | 9.18 / 9.26 |
| 200 | 3.2e-6 | 2.2e-9 | 0.142 | 9.16 / 9.54 |

At 100 cells and above the interface equilibrium holds to 3e-6. At 50 and 64 cells a localized pressure error of 2.5e-3 p0
appears after 0.4 flow-throughs in water cells with traces of air (open, cause not identified). An NVIDIA A100 run (NVHPC 25.11)
at N = 100 agrees with the CPU values. The example is not in the test suite. See @ref validation_multiphase_riemann.
