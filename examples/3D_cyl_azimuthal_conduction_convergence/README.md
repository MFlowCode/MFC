# Azimuthal Fourier conduction convergence (3D cylindrical)

Third member of the conduction verification set, after `examples/1D_conduction_convergence`
(Cartesian) and `examples/2D_axisym_conduction_convergence` (radial, plus the axis cell).
Both of those run at `p = 0` and therefore have no azimuthal direction at all. This case
covers it.

In 3D cylindrical mode (`grid_geometry == 3`) the third coordinate is the azimuth and `dz` is
in **radians**, so the azimuthal term carries two metric factors, `(1/r^2) d2T/dtheta2`:
one for the gradient and one for the divergence. `m_rhs.fpp` divides the conduction flux
difference by `dz(l)` alone, so both factors have to come out of `grid_spacing` in
`m_conduction.fpp`, which is why that direction uses `y_cc(y)**2 * (z_cc(z+1) - z_cc(z))`.

## Exact solution

Ideal gas, so `T = p / ((Gamma - 1) * rho * cv)`. Setting

    rho(r, theta) = RHO0 / (1 + A (r/R)^2 cos(theta))

at uniform pressure gives exactly

    T(r, theta) = T0 (1 + A (r/R)^2 cos(theta)),   T0 = P0 / ((Gamma - 1) RHO0 cv),

which is single-valued and smooth on the axis. Its cylindrical Laplacian splits into

    (1/r) d/dr (r dT/dr) = +4 A T0 cos(theta) / R^2      (radial half)
    (1/r^2) d2T/dtheta2  = -1 A T0 cos(theta) / R^2      (azimuthal half)

and the sum, `3 A T0 cos(theta) / R^2`, is **independent of `r`**. At `t = 0` the velocity is
zero and the pressure uniform, so every Euler flux and every cylindrical geometric source
vanishes and the whole right-hand side is conduction:

    (rho E(dt) - rho E(0)) / dt = k grad^2 T + O(dt).

The radial half is reproduced *exactly* by the discretization: the 3D cylindrical `r`-grid is
uniform (`m_grid.f90` only half-cells the axis for `grid_geometry == 2`), `T` is quadratic in
`r`, and both the two-point face gradient and the two-point average of `F` are exact for that.
Subtracting it therefore isolates the azimuthal half, and `compare_analytic.py` reports it ring
by ring. Dropping either metric factor multiplies the azimuthal half by `r^2`, which no
refinement removes.

Boundaries: periodic in `x` and in `theta`, so the azimuthal direction has no boundary error at
all. Excluded from the bulk error are cell `k = 0`, which reads the ghost across the axis (the
documented non-converging axis error), and the two outermost cells: `k = NR-1` reads the
mirrored wall ghost where `dT/dr` is not actually zero, and `k = NR-2` picks up an `O(dt)`
share of that through the later Runge-Kutta stages.

## Running

    for np_ in 32 64 128 256; do
        NP=$np_ ./mfc.sh run examples/3D_cyl_azimuthal_conduction_convergence/case.py -n 1
        NP=$np_ ./build/venv/bin/python3 examples/3D_cyl_azimuthal_conduction_convergence/compare_analytic.py
    done

## Observed results

`NR = 32`, `NX = 32`, exact `= 7.5e-01 cos(theta)`. The azimuthal truncation error of a
centered second difference is `-dtheta^2/12` of the azimuthal half, i.e. `dtheta^2/36` of the
total, and nothing else contributes:

| NP  | bulk L2 rel err | order | predicted `dtheta^2/36` |
|-----|-----------------|-------|-------------------------|
| 32  | 1.070e-03       | --    | 1.071e-03               |
| 64  | 2.676e-04       | 2.00  | 2.677e-04               |
| 128 | 6.689e-05       | 2.00  | 6.693e-05               |
| 256 | 1.670e-05       | 2.00  | 1.673e-05               |

Azimuthal half recovered ring by ring, as a fraction of its exact value (`NP = 64`, so the
expected value is `1 - dtheta^2/12 = 0.999197` at every radius):

| r      | 0.094    | 0.531    | 1.031    | 1.844    |
|--------|----------|----------|----------|----------|
| got/exact | 0.999199 | 0.999197 | 0.999197 | 0.999197 |

Flat across a 20x span in `r`, i.e. a 380x span in `r^2`. Refining `NR` at fixed `NP` leaves
the bulk error unchanged (2.676e-04 at `NR = 32`, 2.656e-04 at `NR = 64`), confirming that the
radial half contributes no error and that the table above is a pure azimuthal measurement.

Without the `r^2` in `grid_spacing` the same `NP = 64` run gives a bulk L2 of 3.542e-01 and a
ring-by-ring azimuthal ratio of 0.0088, 0.282, 1.063, 3.397 at those four radii -- exactly
`r^2`, and independent of resolution.
