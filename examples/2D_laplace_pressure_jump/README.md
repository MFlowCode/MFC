# Laplace Pressure Jump of a 2D Drop

A quarter of a 2D water drop (r0 = 0.15 m, sigma = 8 N/m) in air at rest, with symmetry planes at x = 0 and y = 0, relaxes to
equilibrium under surface tension (6-equation model with `surface_tension`). `check_laplace.py` averages the pressure inside and
outside the drop over the last 20 % of the saves and compares the jump with sigma/R, R the equivalent radius of the water area,
and prints the largest residual velocity (spurious currents).

    ./mfc.sh run examples/2D_laplace_pressure_jump/case.py -t pre_process simulation post_process -n 1 -- -N 50
    python3 examples/2D_laplace_pressure_jump/check_laplace.py examples/2D_laplace_pressure_jump --case-args "-N 50"

`case.py` takes `-N` (cells per direction, default 50), `--tend` (0.5 s) and `--saves` (100). Measured with GNU 13.3 (CPU, one rank)
and an NVIDIA A100 (NVHPC 25.11, `--case-optimization`):

| N | platform | dp [Pa] | vs sigma/R | max abs(u) [m/s] |
|---:|---|---:|---:|---:|
| 25 | CPU | 53.633 | +1.02 % | 0.053 |
| 50 | CPU | 53.279 | +0.010 % | 0.038 |
| 50 | A100 | 53.280 | +0.013 % | 0.037 |
| 100 | A100 | 53.185 | -0.25 % | 0.042 |

At 100 cells the jump still varies by 0.9 % over the averaging window. The residual velocity (spurious currents) is about 20 % of
the capillary speed sqrt(sigma/(rho_l R)) = 0.23 m/s on every grid.

See @ref validation_surface_tension.
