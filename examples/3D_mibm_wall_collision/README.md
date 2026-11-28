# 3D Soft-Sphere Collision of a Heavy Sphere with a Wall

A sphere (R = 5 m, mass 1e6, heavy enough that the gas force is negligible) moves at 3.9 m/s toward the wall `y = 0`
(`bc_y%beg = -15`) and rebounds under `collision_model = 1`. For a homogeneous sphere on a fixed wall with Coulomb friction
coefficient mu and normal restitution e, rigid-body impact theory (\cite Maw1976) gives, for incoming normal and tangential
speeds v_n and v_t:

- normal: `v_n' = -e v_n`;
- sliding throughout the contact (`tan(angle) >= 3.5 mu (1 + e)`): `v_t' = v_t - mu (1 + e) |v_n|`;
- sticking before the end of the contact (`tan(angle) < 3.5 mu (1 + e)`): `v_t' = (5/7) v_t`.

The threshold and the two tangential results follow from the impulse balance of a sphere with moment of inertia 2/5 m R^2 and
no initial spin. The setup follows the verification cases of \cite DiMaio2004.

`check_collision.py` reads the last `restart_data/ib_state_*.dat` and prints `e_n = -vy_out / vy_in` and the tangential speed
against the applicable rigid-body result:

    python3 examples/3D_mibm_wall_collision/check_collision.py <case directory> --angle 30 --e 0.98 --mu 0.092

Measured with 16^3 cells, first-order reconstruction, air at 1e5 Pa, a contact onset half a step into a step and `collision_time` = 20 steps
(the gas plays no role, so the grid and the reconstruction do not matter):

| Case | Prescribed | Measured | Error |
|---|---|---|---|
| normal, e = 0.5 | 0.5000 | 0.4765 | -4.70 % |
| normal, e = 0.8 | 0.8000 | 0.7877 | -1.54 % |
| normal, e = 0.98 | 0.9800 | 0.97948 | -0.05 % |
| 30 deg, mu = 0, e = 0.98: `vx_out` | 1.9500 | 1.94998 | 0.00 % |
| 30 deg, mu = 0.092, e = 0.98: `vx_out` (sticking, tan 30 deg = 0.577 < 0.638) | 1.3929 | 1.3958 | +0.21 % |
| 50 deg, mu = 0.092, e = 0.98: `vx_out` (sliding, tan 50 deg = 1.19 > 0.638) | 2.5309 | 2.5310 | +0.00 % |

Error in e_n against the steps per collision (`--steps-per-collision`):

| e | 20 | 40 | 80 |
|---|---|---|---|
| 0.5 | -4.70 % | -2.33 % | -1.16 % |
| 0.8 | -1.54 % | not run | -0.37 % |

With `--wall-bc -2` the face y = 0 is reflective (a symmetry plane). The sphere then meets its mirror image: the image contact has
half the effective mass, twice the overlap and twice the closing speed, so the normal force equals the wall's, and the image's
surface moves with the body's along the face, so there is no friction. Measured on an MI300A (CCE 19, OpenMP offload), same setup:

| Case | slip wall (`-15`) | reflective (`-2`) | expected |
|---|---|---|---|
| normal, e = 0.5, 20 steps: e_n | 0.47651 | 0.47651 | 0.5 (same discretisation error) |
| normal, e = 0.5, 80 steps: e_n | 0.49420 | 0.49420 | 0.5 |
| normal, e = 0.8, 20 steps: e_n | 0.78767 | 0.78767 | 0.8 |
| 50 deg, e = 0.98: `vx_out` (slip with mu = 0, reflective with mu = 0.092) | 2.987518 | 2.987518 | 2.98757 (v_t unchanged) |

A sphere of mass 3000 (the gas force now matters) over the reflective face and the same sphere paired with its mirror image in a
domain twice as tall, with slip walls far away, follow the same path: velocities agree to 1.7e-5 of the impact speed and centroids
to 2e-9 R.

The restitution error is first order in the steps per collision (it halves at each doubling). At the default 20 steps expect a few
percent at low e. Its size and sign depend on where in a step the contact begins: an earlier run with the contact beginning on a step
boundary gave +2.17 % and +0.69 % for e = 0.5 and 0.8 at 20 steps. The test registration therefore starts the contact half a step into
a step, which keeps the contact timing the same on every compiler.

## Run

    ./mfc.sh run examples/3D_mibm_wall_collision/case.py -t pre_process simulation -n 1 --gpu mp
    python3 examples/3D_mibm_wall_collision/check_collision.py examples/3D_mibm_wall_collision

`case.py` takes `--angle`, `--e`, `--mu`, `--wall-bc`, `--speed`, `--cells`, `--steps-per-collision`, `--lead-steps` and `--weno_order`. The suite registers
the normal impact, the 50 degree sliding impact and the 50 degree impact on a reflective face as `3D -> IBM -> Wall Collision -> ...`. The sticking case is not a golden: where
the contact point stops sliding, the friction direction is the sign of a velocity that goes to zero, so the rebound depends on round-off.
The tests are skipped in the single-precision lane, where one-step shifts of the contact onset flip the rebound.

## Limits

One sphere against a wall, at 3.9 m/s and three angles. Body-body collisions, spin and other speeds are not covered. The
collision force is first order in the time step; this check does not separate the force evaluation from the detection. The
`3D_mibm_sphere_head_on_collision` example places the start position with `sin(angle)` in both offsets, where a straight approach
needs `cos(angle)` for the y offset; this example uses `cos`. See @ref validation.
