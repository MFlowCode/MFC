# 2D Constant Acceleration of a Two-Way Immersed Circle

A circle (`moving_ibm = 2`, 1000 kg/m) in air at rest is accelerated by the body force `g_x = 1e5 m/s^2`. The gas force on it is
negligible against the body force, so the centroid must follow `x = g t^2 / 2` and `v = g t`. The check is sensitive to the
Runge-Kutta stage weights of the position update and to nothing else: a wrong stage velocity gives an exact velocity and a wrong
position.

`check_kinematics.py` reads the last `restart_data/ib_state_*.dat` and prints the relative errors of the centroid and the velocity:

    python3 examples/2D_mibm_constant_acceleration/check_kinematics.py <case directory> [--tol 1e-10]

| cells per radius | steps | t [s] | position error | velocity error |
|---:|---:|---:|---:|---:|
| 5 | 200 | 2.34e-6 | -3.6e-13 | -2.3e-13 |
| 10 | 200 | 1.17e-6 | -2.9e-13 | -1.9e-13 |

Both are round-off. A build with the stage weights of the position update before the fix of #2040 gave a position error of 1.2e-3
with the same exact velocity (measured earlier, not repeated here).

## Run

    ./mfc.sh run examples/2D_mibm_constant_acceleration/case.py -t pre_process simulation -n 1 --gpu mp
    python3 examples/2D_mibm_constant_acceleration/check_kinematics.py examples/2D_mibm_constant_acceleration

`case.py` takes `--cells_per_R` (default 5), `--steps` (200), `--g` (1e5) and `--weno_order` (default 1: first-order reconstruction keeps the check free of WENO weights, which amplify round-off; 5 gives WENO5 with the same kinematics). The suite registers it as
`2D -> IBM -> Moving Circle -> Constant Acceleration` with the golden `ib_forces.dat` history (force, velocity and centroid); it is skipped in the single-precision lane, where pressure round-off at 1e5 Pa moves the quiescent gas by about 1e-3 of the signal.

## Limits

Only a constant acceleration is checked, not a force that varies in time. With `moving_ibm = 1` and an analytic velocity such as
`U0 sin(w t)`, the velocity is evaluated at the step time rather than the stage time, and the centroid is off by about 1e-3
relative (cause not verified). See @ref validation.
