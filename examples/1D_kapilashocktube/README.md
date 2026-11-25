# Water-Air Shock Tube

Water at 1e9 Pa (stiffened gas, gamma 4.4, pi_inf 6e8 Pa, 1000 kg/m^3) against air at 1e5 Pa (ideal gas, 50 kg/m^3), with the
interface at x = 0.7 m, run to t = 241 us with the 5-equation model (`alt_soundspeed`, HLLC, WENO5-M, SSP-RK3). A shock runs into
the air and a rarefaction into the water. `check_riemann.py` samples the exact two-material Riemann solution (Toro's solver with
p + pi_inf in each material) at the cell centres and prints the L1 errors of rho, u and p, the observed order between grids,
the material interface position from the volume fraction, and the star-plateau errors.

    ./mfc.sh run examples/1D_kapilashocktube/case.py -t pre_process simulation post_process -n 1 -- -N 1000
    python3 examples/1D_kapilashocktube/check_riemann.py <run 250> <run 500> <run 1000> <run 2000> --case-args "-N 250" "-N 500" "-N 1000" "-N 2000"

`case.py` takes `-N` (cells, default 1000); dt scales with 1/N at the fixed end time. Exact star state: p* = 1.4190e7 Pa,
u* = 482.61 m/s, rho* = 804.44 (water) and 288.17 kg/m^3 (air). Measured with GNU 13.3 (CPU, one rank); L1 normalized by the
jump across the tube:

| N | L1 rho | L1 u | L1 p | interface error [cells] |
|---:|---:|---:|---:|---:|
| 250 | 3.69e-3 | 7.98e-3 | 2.08e-3 | 0.02 |
| 500 | 2.32e-3 | 4.52e-3 | 1.06e-3 | 0.01 |
| 1000 | 1.11e-3 | 1.61e-3 | 5.28e-4 | 0.01 |
| 2000 | 7.18e-4 | 1.16e-3 | 2.69e-4 | 0.01 |

The pressure error is first order; density and velocity converge at order 0.5 to 1.5 between grids (the contact and the shock
dominate). An NVIDIA A100 run (NVHPC 25.11, OpenACC) at N = 1000 gives the same L1 errors to the digits shown. See @ref validation_multiphase_riemann.
