# 1D Fourier conduction convergence

Verifies the Fourier heat conduction term `div(k grad T)` on the energy equation against
the exact Laplacian of a sinusoidal temperature field.

## Exact solution

Ideal gas (`pi_inf = 0`), so `T = p / ((Gamma - 1) * rho * cv)`. Setting

    rho(x) = RHO0 / (1 + A sin(2 pi x / L))

at uniform pressure gives exactly

    T(x) = T0 (1 + A sin(2 pi x / L)),   T0 = P0 / ((Gamma - 1) RHO0 cv).

At `t = 0` the velocity is zero and the pressure is uniform, so every Euler flux vanishes
and the entire right-hand side is conduction:

    d(rho E)/dt = k d2T/dx2 = -k T0 A (2 pi / L)^2 sin(2 pi x / L).

One time step therefore measures the conduction term in isolation, up to `O(dt)` time
error and `O(dx^2)` space error. `compare_analytic.py` forms `(rho E(dt) - rho E(0)) / dt`
from `restart_data/lustre_{0,1}.dat` and compares it to that expression.

## Running

    for nx in 50 100 200 400; do
        NX=$nx ./mfc.sh run examples/1D_conduction_convergence/case.py -n 1
        NX=$nx ./build/venv/bin/python3 examples/1D_conduction_convergence/compare_analytic.py
    done

## Observed convergence

| NX  | L2 error     | relative     | order |
|-----|--------------|--------------|-------|
| 50  | 9.174206e-06 | 1.314570e-03 | --    |
| 100 | 2.293085e-06 | 3.285757e-04 | 2.000 |
| 200 | 5.775688e-07 | 8.275971e-05 | 1.989 |
| 400 | 1.442740e-07 | 2.067299e-05 | 2.001 |
| 800 | 4.396086e-08 | 6.299142e-06 | 1.715 |

Second order through NX=400. The drop at NX=800 is the `O(dt)` time-splitting floor, which
by then is comparable to the spatial error; refining `dt` restores second order.

The result is unchanged on 2 ranks (`-n 2`), which exercises the MPI temperature halo
exchange.
