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

Measured with 24^3 cells, 30 lead steps and `collision_time` = 20 steps (the gas plays no role, so the grid does not matter):

| Case | Prescribed | Measured | Error |
|---|---|---|---|
| normal, e = 0.5 | 0.5000 | 0.5108 | +2.17 % |
| normal, e = 0.8 | 0.8000 | 0.8055 | +0.69 % |
| normal, e = 0.98 | 0.9800 | 0.97984 | -0.02 % |
| 30 deg, mu = 0, e = 0.98: `vx_out` | 1.9500 | 1.95000 | exact; e_n equals the normal impact |
| 30 deg, mu = 0.092, e = 0.98: `vx_out` (sticking, tan 30 deg = 0.577 < 0.638) | 1.3929 | 1.3906 | -0.16 % |
| 50 deg, mu = 0.092, e = 0.98: `vx_out` (sliding, tan 50 deg = 1.19 > 0.638) | 2.5309 | 2.5308 | -0.00 % |

Error in e_n against the steps per collision (`--steps-per-collision`):

| e | 20 | 40 | 80 |
|---|---|---|---|
| 0.5 | +2.17 % | +1.12 % | +0.57 % |
| 0.8 | +0.69 % | not run | +0.18 % |

The restitution error is first order in the steps per collision (it halves at each doubling). At the default 20 steps expect up
to about 2 % error in e at low e.

## Run

    ./mfc.sh run examples/3D_mibm_wall_collision/case.py -t pre_process simulation -n 1 --gpu mp
    python3 examples/3D_mibm_wall_collision/check_collision.py examples/3D_mibm_wall_collision

`case.py` takes `--angle`, `--e`, `--mu`, `--speed`, `--cells`, `--steps-per-collision` and `--lead-steps`. The suite registers
the normal impact and the 30 degree sticking-friction impact as `3D -> IBM -> Wall Collision -> ...`.

## Limits

One sphere against a wall, at 3.9 m/s and three angles. Body-body collisions, spin and other speeds are not covered. The
collision force is first order in the time step; this check does not separate the force evaluation from the detection. The
`3D_mibm_sphere_head_on_collision` example places the start position with `sin(angle)` in both offsets, where a straight approach
needs `cos(angle)` for the y offset; this example uses `cos`. See @ref validation.
