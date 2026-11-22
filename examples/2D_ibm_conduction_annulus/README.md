# 2D Isothermal Conduction Annulus: Image-Point Interpolation

Steady conduction between two isothermal immersed boundaries: a circle (radius `R`, wall temperature `T1 = 2`) inside a flat ring read
from `ring.stl` whose inner edge is the circle `r = 2R` (wall temperature `T2 = 3`). The gas is quiescent with constant conductivity
(`k_therm = 0.1`, `cv = 1`, `gamma = 1.4`, nondimensional), so the steady temperature is exactly

    T(r) = T1 + (T2 - T1) ln(r / R) / ln(r2 / R).

The case compares the two image-point weightings of the ghost-cell method: inverse distance (the solver default) and triangular-shaped
cloud (`ib_ip_interp = 2`, TSC). The gas starts from the exact profile on 10 concentric constant shells (uniform pressure) and runs 8 slowest
conduction times, so the error that remains is the steady discretization error of the wall condition.
`check_annulus.py` reads the last `restart_data/lustre_<k>.dat` and prints, relative to `T2 - T1`:

- the wall temperature error: a quadratic least-squares fit through the gas cells `0.5 dx < s < 3.5 dx` from the wall, extrapolated to the wall, minus the wall value;
- the wall-normal gradient error of the same fit, relative to the exact gradient;
- the relative L2 and the maximum error of `T` over the gas cells at least half a cell from either wall.

## Results

Signed wall error (MFC - exact) / (T2 - T1), inverse distance (`idw`) and TSC, 1 rank, 8 conduction times:

| cells / R | idw inner | idw outer | TSC inner | TSC outer |
|---:|---:|---:|---:|---:|
| 10 | -1.83e-2 | +1.18e-2 | +3.68e-3 | +1.31e-4 |
| 20 | -1.18e-2 | +7.71e-3 | +7.67e-4 | +6.05e-5 |
| 40 | -5.00e-3 | +3.42e-3 | +1.87e-4 | +2.96e-5 |
| observed order 10-20, 20-40 | 0.6, 1.2 | 0.6, 1.2 | 2.3, 2.0 | 1.1, 1.0 |

Whole-gap error of `T` over the gas cells:

| cells / R | idw L2 | TSC L2 | idw max | TSC max |
|---:|---:|---:|---:|---:|
| 10 | 9.36e-3 | 1.45e-3 | 4.57e-2 | 2.66e-3 |
| 20 | 5.32e-3 | 3.90e-4 | 1.98e-2 | 8.68e-4 |
| 40 | 2.02e-3 | 1.86e-4 | 8.94e-3 | 3.39e-4 |

Wall-normal gradient error (fit estimate, relative), inner wall: idw -5.0e-3, +2.7e-3, +1.4e-2; TSC -2.3e-2, -5.8e-3, -4.3e-4 at 10, 20, 40 cells/R.
The inverse-distance gradient does not converge and TSC converges at about second order on the inner wall; the estimator is a fit through the
first cells and is not a reference for the true wall flux.

TSC lowers the imposed wall temperature error by a factor of 5 to 130 and the whole-gap L2 error by a factor of 6 to 14 at these resolutions.
The maximum error converges at order 1.2 (inverse distance) and 1.4 to 1.6 (TSC). With inverse distance the wall error order is erratic between 0.6 and 1.2.

## Run

    ./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 --gpu                                  # inverse distance, 12 cells/R
    ./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 --gpu -- --interp 2 --cells_per_R 20   # TSC
    python3 examples/2D_ibm_conduction_annulus/check_annulus.py examples/2D_ibm_conduction_annulus

`case.py` takes `--cells_per_R` (default 12), `--interp` (1 or 2) and `--steps` (default: 8 slowest conduction times). The table above uses
10, 20 and 40 cells/R. Single-precision builds write single-precision restart files: pass `--dtype float32` to the checker.
The suite registers `2D -> IBM -> Conduction Annulus -> Inverse Distance` and `-> TSC` at 12 cells/R and 300 steps (about 3 conduction times).
`ring.stl` is `ring_stl.py ring.stl 2 5 1024`.

## Limits

- Constant conductivity and a quiescent gas. The temperature-dependent conductivity of air, a moving gas and the isothermal wall heat flux
  (`ib_heat_wrt` needs chemistry) are not covered here.
- The run is inviscid. With viscosity (`fluid_pp(1)%Re(1) = 10`, `viscous = T`) the 20 cells/R runs of both weightings returned NaN in the last step
  (cause not investigated); the 12 cells/R viscous runs completed.
- Starting from a uniform temperature instead of the exact profile leaves acoustic ringing in the closed gap (Mach number about 1e-2) that an
  inviscid run does not damp; the whole-gap error then stays between 2.6e-3 and 6.0e-3 for TSC at 10, 20 and 40 cells/R and hides the discretization error.

See @ref validation_tsc_annulus.
