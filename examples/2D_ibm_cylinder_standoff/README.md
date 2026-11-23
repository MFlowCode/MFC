# Bow-Shock Stand-Off of a Cylinder and a Sphere, Mach 1.5-3

A fixed no-slip body (immersed boundary) in supersonic flow, an impulsive start with extrapolation on every face, inviscid. The
detached bow shock settles at a stand-off distance from the body that depends on the Mach number and the body shape, and the
resolution needed to place it is high. `standoff.py` measures the distance on the stagnation line (where the pressure first reaches
the mean of the free-stream and post-shock pressures, scanning inward from the inflow) and compares it with Billig's correlation.
`bodyfitted_sphere.py` is an independent reference: an axisymmetric Euler solver on a body-fitted polar mesh around the sphere
(NumPy, no immersed boundary). The full results and limits are on the @ref validation_standoff page.

## Run

    ./mfc.sh run examples/2D_ibm_cylinder_standoff/case.py -t pre_process simulation -n 1 --gpu -- \
        --dim 3 --cells_per_D 40 --steps 0 --tstar 40 --frames 8
    python3 examples/2D_ibm_cylinder_standoff/standoff.py examples/2D_ibm_cylinder_standoff
    python3 examples/2D_ibm_cylinder_standoff/bodyfitted_sphere.py --mach 2 --nr 160 --nphi 160

`--dim 2` is the cylinder (the default), `--dim 3` the sphere; `--mach`, `--xu`, `--xd` and `--hy` set the Mach number and the domain
(a cylinder at Mach 1.5 needs `--xu 3.5`, because its shock stands about 1.5 D from the body). The output goes next to `case.py`: copy
the directory for each run. A build without MPI adds `--serial_io`. Without arguments the case is the CI-sized test: the cylinder at
Mach 2 with 16 cells per diameter in a 4 D box and 150 fixed steps, before the shock has settled; it runs in the suite as
`2D -> IBM -> Bow Shock Standoff`.

## Results

NVIDIA A100, NVHPC 25.11, double precision, `master` at 0cdfdf02a. Sphere, inviscid, no-slip, `ib_ip_interp = 1`, domain
[-1.5 D, 2.5 D] x [-2 D, 2 D]^2, run to t u/R = 40, the stand-off \f$\Delta/R\f$ averaged over the last three of eight saves. The shock position
oscillates by 0.003-0.008 R between saves at every grid. "Body-fitted" is `bodyfitted_sphere.py` at 160 x 160 cells (M 1.5: 0.6131, M 2: 0.3498, M 3: 0.2179).

| M | cells/D | Delta/R | Billig | error vs Billig | error vs body-fitted |
|---:|---:|---:|---:|---:|---:|
| 1.5 | 20 | 0.6988 | 0.6036 | +15.8 % | +14.0 % |
| 1.5 | 40 | 0.6650 | 0.6036 | +10.2 % | +8.5 % |
| 2 | 20 | 0.4167 | 0.3215 | +29.6 % | +19.1 % |
| 2 | 30 | 0.4025 | 0.3215 | +25.2 % | +15.1 % |
| 2 | 40 | 0.3939 | 0.3215 | +22.5 % | +12.6 % |
| 2 | 60 | 0.3776 | 0.3215 | +17.5 % | +7.9 % |
| 3 | 20 | 0.2482 | 0.2050 | +21.1 % | +13.9 % |
| 3 | 40 | 0.2469 | 0.2050 | +20.4 % | +13.3 % |

The body-fitted solver gives 0.3489 and 0.3498 R at M 2 on 80 and 160 cells (Billig is 8.8 % lower), 0.6130 and 0.6131 at M 1.5, and
0.2153 and 0.2179 at M 3 (not converged to better than about 1 %).

The cylinder (`--dim 2`, GNU 13 CPU, run to t u/R = 60) does not settle in that time: at Mach 2 the stand-off is 1.414 R at 20 cells/D
(Billig 1.241 R, +14.0 %), 1.464 R at 40 cells/D, still rising (1.425, 1.458, 1.508 R over the last three saves), and 1.377 R at 80 cells/D (+11.0 %, also still rising). The sphere is
the validated quantity; the cylinder is the CI-sized regression case and a slower, larger-subsonic-region problem.

## Limits

Resolution-limited, not a measurement problem: the shock stands 4-12 cells from the sphere at these grids. Mach 3 is unsteady (the
stand-off oscillates by 0.007-0.036 R) and was not examined further. Inviscid only; no reacting variant (see `2D_ibm_reacting_surface`).
