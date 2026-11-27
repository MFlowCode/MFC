# Water Disk Advected in Air

A water disk (stiffened gas, pi_inf 1e9 Pa, 1000 kg/m^3) in air (1 kg/m^3), both at 1e5 Pa and u = v = 100 m/s, crosses a
periodic unit box once along the diagonal with the 6-equation model. The exact solution keeps p and u uniform and returns the disk
to its initial field. `check_advection.py` prints the largest pressure and velocity disturbance over all saves, the L1 error of
the volume fraction per unit interface length and the width of the mixed zone, for one or several grids.

    ./mfc.sh run examples/2D_advection/case.py -t pre_process simulation post_process -n 1 -- -N 100
    python3 examples/2D_advection/check_advection.py examples/2D_advection --case-args "-N 100"

`case.py` takes `-N` (cells per direction, default 100), `--periods`, `--cfl` (0.4), `--saves` (20) and `--eps` (trace volume
fraction, 1e-8):

| N | max dp/p0 | max du/U | L1 [cells] | mixed-zone width, initial / final [cells] |
|---:|---:|---:|---:|---:|
| 50 | 3.7e-11 | 1.1e-11 | 0.066 | 9.17 / 9.25 |
| 100 | 2.5e-11 | 5.6e-12 | 0.095 | 9.18 / 9.26 |
| 200 | 1.5e-10 | 1.9e-10 | 0.142 | 9.16 / 9.54 |

The interface equilibrium holds to round-off at every grid. The deck runs with `mpp_lim = F`: the volume-fraction clip of
`mpp_lim` is not matched by the energy, so where a trace fraction undershoots it shifts the water pressure (to 2.5e-3 p0 at 50
cells, 3e-6 at 100; #2476). Measured on the CPU with NVHPC 25.11 (four ranks); an NVIDIA A100 (NVHPC 25.5) gives 3.8e-11 and 3.4e-11 at N = 50
and 100. The example is not in the test suite. See @ref validation_multiphase_riemann.
