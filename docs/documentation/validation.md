@page validation Validation

# Validation

This page records what MFC has been checked against: exact solutions, reference data and conservation laws.
Each section is one family of checks. It gives the status, the code version, the setup, the reference, the measured values with signed errors, and the limits.
Every check has an example in `examples/` with a checker script and a README, and a small version of it runs in the test suite (`./mfc.sh test`).

Status words: **validated** (agrees with the reference within the stated error), **partially validated** (agrees where measured; a regime or quantity is open), **open discrepancy** (disagrees, cause known or not), **comparison only** (reference data are qualitative).
Errors are (MFC − reference)/reference.

## Moving-body kinematics {#sec-validation-moving-kinematics}

**Status:** validated. **Code:** `master` at [e29d3c25](https://github.com/MFlowCode/MFC/commit/e29d3c25d), with the stage-consistent position update of [#2040](https://github.com/MFlowCode/MFC/pull/2040). **Examples:** `examples/2D_mibm_constant_acceleration`, `examples/3D_mibm_wall_collision`.

Two-way coupled bodies (`moving_ibm = 2`) are advanced from the force on them. Particle-resolved runs depend on the integration being exact where an exact answer exists, and on collisions returning the prescribed restitution.

### Constant acceleration {#sec-validation-moving-kinematics-accel}

| Item | Value |
| :--- | :--- |
| Body | circle, R = 50 µm, 1000 kg/m, `moving_ibm = 2` |
| Gas | air at rest, 1e5 Pa, 1.2 kg/m³, inviscid; the gas force is negligible against the body force |
| Drive | `bf_x = T`, `g_x = 1e5` m/s², `k_x = w_x = p_x = 0` |
| Numerics | WENO5, HLLC, SSP-RK3, fixed \f$\Delta t\f$ at acoustic CFL 0.4, 200 steps, 6R half-width, extrapolation boundaries |
| Reference | exact kinematics: \f$x = g t^2/2\f$, \f$v = g t\f$ |
| Measured | final centroid and velocity from `restart_data/ib_state_*.dat`, checked by `check_kinematics.py` |

| Cells per radius | Position error | Velocity error |
| ---: | ---: | ---: |
| 5 | −1.1e-13 | −6.1e-14 |
| 10 | −9.2e-14 | −5.2e-14 |

Both are round-off. A build with the stage weights of the position update from before #2040 gave a position error of 1.2e-3 with an exact velocity (measured earlier on a build without #2040; not repeated). The check cannot see a wrong stage velocity in the velocity, only in the position, which is why it uses a position.

Limits: only constant acceleration is checked. For a prescribed body (`moving_ibm = 1`) with an analytic velocity such as \f$U_0\sin\omega t\f$, the velocity is evaluated at the step time rather than the stage time, and the centroid is off by about 1e-3 relative (measured earlier; the cause is not verified and is not addressed here).

### Sphere-wall collision {#sec-validation-moving-kinematics-collision}

| Item | Value |
| :--- | :--- |
| Body | sphere, R = 5 m, mass 1e6 (gas force negligible), impact speed 3.9 m/s |
| Wall | `y = 0`, `bc_y%%beg = -15`; other boundaries extrapolation |
| Model | `collision_model = 1`, `collision_time` = N\f$\Delta t\f$, \f$\Delta t\f$ = 5e-6 s, `coefficient_of_restitution` e, `ib_coefficient_of_friction` \f$\mu\f$ |
| Numerics | WENO5, HLLC, SSP-RK3, fixed \f$\Delta t\f$, 24³ cells |
| Reference | rigid-body impact of a sphere on a fixed wall with Coulomb friction and no tangential elasticity (\cite Maw1976): \f$v_n' = -e v_n\f$; sliding throughout (\f$\tan\theta \ge 3.5\mu(1+e)\f$): \f$v_t' = v_t - \mu(1+e)|v_n|\f$; sticking (\f$\tan\theta < 3.5\mu(1+e)\f$): \f$v_t' = \frac{5}{7}v_t\f$. Setup after \cite DiMaio2004 |
| Measured | \f$e_n = -v_y^\text{out}/v_y^\text{in}\f$ and \f$v_x^\text{out}\f$ from the last `ib_state` record, by `check_collision.py` |

| Case | Prescribed | Measured | Error |
| :--- | ---: | ---: | ---: |
| normal, e = 0.5 | 0.5000 | 0.5108 | +2.17 % |
| normal, e = 0.8 | 0.8000 | 0.8055 | +0.69 % |
| normal, e = 0.98 | 0.9800 | 0.97984 | −0.02 % |
| 30°, \f$\mu\f$ = 0, e = 0.98: \f$v_x^\text{out}\f$ | 1.9500 | 1.95000 | 0 |
| 30°, \f$\mu\f$ = 0.092, e = 0.98: \f$v_x^\text{out}\f$ (sticking) | 1.3929 | 1.3906 | −0.16 % |
| 50°, \f$\mu\f$ = 0.092, e = 0.98: \f$v_x^\text{out}\f$ (sliding) | 2.5309 | 2.5308 | −0.00 % |

Error in \f$e_n\f$ against the steps per collision:

| e | 20 | 40 | 80 |
| ---: | ---: | ---: | ---: |
| 0.5 | +2.17 % | +1.12 % | +0.57 % |
| 0.8 | +0.69 % | not run | +0.18 % |

The restitution error is first order in the steps per collision. At the default of 20 steps, expect up to about 2 % error in e at low e. Without friction the tangential velocity is conserved exactly. The sticking and sliding regimes both match the rigid-body result.

Limits: one sphere against a wall, three angles. Body-body collisions, spin and other speeds are not covered. The first-order behavior is in the collision force; whether it is the force evaluation or the contact detection is not separated.

### Reproduce {#sec-validation-moving-kinematics-run}

```
./mfc.sh run examples/2D_mibm_constant_acceleration/case.py -t pre_process simulation -n 1
python3 examples/2D_mibm_constant_acceleration/check_kinematics.py examples/2D_mibm_constant_acceleration
./mfc.sh run examples/3D_mibm_wall_collision/case.py -t pre_process simulation -n 1 -- --angle 30 --e 0.98 --mu 0.092
python3 examples/3D_mibm_wall_collision/check_collision.py examples/3D_mibm_wall_collision --angle 30 --e 0.98 --mu 0.092
```

The two small cases run in the test suite as `2D -> IBM -> Moving Circle -> Constant Acceleration` and `3D -> IBM -> Wall Collision`.
The tables above were measured with GNU 13.3 (CPU, one rank). The golden tolerances are 1e-10 (constant acceleration) and 1e-6 (collision).

## Adding a section

One `##` section per family, with an anchor `{#sec-validation-<family>}`: status and code version, one setup table per check, the reference and how it was obtained, a results table with signed errors, and the limits.
Put each check's example, checker script and README under `examples/`, and register a small version in `toolchain/mfc/test/cases.py`.
Cite references with `\cite` keys from `docs/references.bib`.
