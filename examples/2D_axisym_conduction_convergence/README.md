# Axisymmetric Fourier conduction convergence

Axisymmetric twin of `examples/1D_conduction_convergence`. It verifies the radial part of
`div(k grad T)` on a cylindrical grid, and in particular the axis cell `k = 0`, which has no
face pair and is therefore the one cell fed by `s_compute_conduction_axis_source`.

## Exact solution

Ideal gas, so `T = p / ((Gamma - 1) * rho * cv)`. Setting

    rho(r) = RHO0 / (1 + A (1 - r^2/R^2) + B (1 - r^4/R^4))

at uniform pressure gives exactly

    T(r) = T0 (1 + A (1 - r^2/R^2) + B (1 - r^4/R^4)),   T0 = P0 / ((Gamma - 1) RHO0 cv),

whose cylindrical Laplacian is finite on the axis:

    (1/r) d/dr (r dT/dr) = -4 A T0 / R^2 - 16 B T0 r^2 / R^4.

At `t = 0` the velocity is zero and the pressure is uniform, so every Euler flux and every
cylindrical geometric source vanishes and the entire right-hand side is conduction:

    d(rho E)/dt = k (1/r) d/dr (r dT/dr).

`compare_analytic.py` forms `(rho E(dt) - rho E(0)) / dt` from `restart_data/lustre_{0,1}.dat`
and compares it to that expression.

`B = 0` (the default) makes the exact answer a **constant**: every radial cell, axis included,
must return the same number. That is the point of the profile — a sign error or a stray factor
in the axis cell is unmissable against a flat field, where a sinusoid would hide it. `B > 0`
makes the answer vary with `r` and turns the same case into an ordinary convergence test.

## The grid at the axis

`src/pre_process/m_grid.f90` gives the axisymmetric grid a **half-width cell at `r = 0`**:
cell 0 spans `[0, h/2]` with center `h/4`, and every cell above it has width `h`. Two cells
therefore behave differently from the bulk:

- **Cell 0**, the axis cell. `s_compute_conduction_axis_source` uses a distance-weighted
  central difference rather than `(T(k+1) - T(k-1)) / (y_cc(k+1) - y_cc(k-1))`; on a uniform
  grid the two are identical, but over this stencil the plain form loses an order and lands
  35% high. The weighted form is exact for a quadratic and converges at second order below.
- **Cell 1**, whose inner face lies on that half-width cell. The two-point face gradient
  `(T(k+1) - T(k)) / (y_cc(k+1) - y_cc(k))` evaluates the gradient at the midpoint of the two
  cell centers, which is the face only on a uniform grid. This face flux is shared by every
  MFC diffusive term (`m_chemistry.fpp` uses the same expression), so the resulting 1.2%
  error at this one cell is a property of that shared discretization, not of the axis source.
  It does not shrink with `NR`, and it is reported separately below.

## Running

    # constant exact answer: the axis-cell check
    for nr in 25 50 100 200 400; do
        NR=$nr ./mfc.sh run examples/2D_axisym_conduction_convergence/case.py -n 1
        NR=$nr ./build/venv/bin/python3 examples/2D_axisym_conduction_convergence/compare_analytic.py
    done

    # r-dependent exact answer: the second-order convergence table
    for nr in 25 50 100 200; do
        NR=$nr BAMP=0.4 ./mfc.sh run examples/2D_axisym_conduction_convergence/case.py -n 1
        NR=$nr BAMP=0.4 ./build/venv/bin/python3 examples/2D_axisym_conduction_convergence/compare_analytic.py
    done

## Observed results

Constant exact answer (`B = 0`, exact `= -4.000000e-02` in every cell):

| NR  | axis cell     | axis rel err | cell 1 rel err | max bulk rel err |
|-----|---------------|--------------|----------------|------------------|
| 25  | -4.00000033e-02 | -8.274e-08 | +3.125e-02     | 8.274e-08        |
| 50  | -4.00000033e-02 | -8.274e-08 | +3.125e-02     | 1.193e-06        |
| 100 | -4.00000033e-02 | -8.274e-08 | +3.125e-02     | 3.413e-06        |
| 200 | -4.00000033e-02 | -8.274e-08 | +3.125e-02     | 7.854e-06        |
| 400 | -3.99999589e-02 | +1.027e-06 | +3.125e-02     | 1.007e-05        |

The axis cell matches the interior to the `O(dt)` time-splitting floor at every resolution:
for a quadratic profile the discretization is exact in space, so nothing else is left. Before
the distance weighting was added the axis cell read `-5.400e-02`, 35% high and independent of
`NR`, which is what this table is built to expose.

Second-order convergence (`B = 0.4`, `r`-dependent exact answer):

| NR  | axis rel err | order | bulk L2 rel err | order |
|-----|--------------|-------|-----------------|-------|
| 25  | 1.171e-03    | --    | 9.623e-04       | --    |
| 50  | 2.856e-04    | 2.04  | 2.348e-04       | 2.03  |
| 100 | 7.062e-05    | 2.02  | 5.799e-05       | 2.02  |
| 200 | 1.757e-05    | 2.01  | 1.448e-05       | 2.00  |

The axis cell converges at the same second order as the bulk. Cell 1 sits at 1.19e-02 at every
resolution, as expected from the shared face-gradient form described above.
