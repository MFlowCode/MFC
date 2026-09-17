# 3D Flapping Flat Plate (prescribed immersed-boundary kinematics)

## Case Set Up

A rigid rectangular plate of aspect ratio 4 glides at a fixed incidence and then begins to flap. Roll is about
the streamwise axis through a hinge at the wing root; pitch is about the body spanwise axis through the same
hinge, leading roll by a quarter cycle so the wing feathers through stroke reversal. Half the span is simulated
with a symmetry plane at the root and a free tip.

| | |
|---|---|
| Chord Reynolds number | 1 000 |
| Mach number | 0.2 |
| Roll amplitude | 30 degrees |
| Pitch amplitude | 20 degrees, mean 5 degrees (the glide incidence) |
| Strouhal number | 0.3, on the peak-to-peak tip excursion |
| Onset | flapping ramps up over half a cycle after 2 convective times of glide |

The motion is set by `patch_ib(1)%kin_model = 1` — the body state is evaluated from the current time at each
Runge-Kutta stage, rather than integrated from analytic velocity expressions. Two consequences worth knowing:
one case-optimized binary serves every Strouhal number and every ensemble member, because nothing about the
motion is compiled in; and restarts reproduce the trajectory exactly, because nothing is integrated.

Note the plate is defined to reach 0.5 c **past** the symmetry plane. An immersed boundary marks only cells
inside the body, so a plate ending exactly on the plane would have fluid on the far side of its root and the root
would be reconstructed as a no-slip face rather than the middle of a continuous wing.

## Numerics

WENO5 with the HLLC Riemann solver and the Thornber low-Mach velocity correction; `mp_weno` off, since the flow
is smooth and the monotonicity clipping only costs accuracy here. Characteristic inflow and outflow with the
generalized relaxation treatment. Fourth-order finite differences for the immersed boundary.

The section is set to four cells thick deliberately. An immersed body reconstructs its interior cells from the
fluid outside, and with fewer than about four cells across there is no real interior left: on a 2D pitching plate
the lift through a ramp was 30 to 40 percent high at two cells and converged by four. Keep that ratio if you
change `dx`.

## Verifying it

The prescribed motion has a closed form, so this case can be checked without any reference data. Read
`D/ib1_forces.dat`, which carries one line per time step, and compare columns `ax, ay` (roll and pitch angles),
`xc, yc, zc` (centroid), `vx, vy, vz` (centroid velocity) and `wx, wy, wz` (lab-frame angular velocity) against

```
phi(t)   = A(tau) phi_0 sin(2 pi f tau),        tau = t - t_0
theta(t) = theta_m + A(tau) theta_0 sin(2 pi f tau + pi/2)
centroid = hinge + Rx(phi) Ry(theta) offset
omega    = phi' e_x + theta' Rx(phi) e_y
```

with `A` the raised-cosine onset envelope. All of these agree to round-off (1e-15) at every step.

The angular velocity is the interesting one. It is **not** the vector of Euler-angle rates: with the rotation
composed as `R = Rx(phi) Ry(theta)`, the lab-frame angular velocity is `phi' e_x + theta' Rx(phi) e_y`. Using
`(phi', theta', 0)` instead — which is what the analytic-expression path effectively does, since it advances the
angles componentwise and then uses the same array as a vector — is wrong once both angles are moving. Sampling
the velocity carried by the body cells adjacent to the fluid and comparing with the rigid-body velocity
distinguishes them: 8e-13 with the correct vector, and 4e-3 to 1e-1 with the Euler rates, at roll angles near 30
degrees.

Two physical checks are worth running as well, because they need no reference either:

- During the glide, before flapping starts, the wing should carry the lift of a finite wing at its incidence.
  Lifting-line theory gives `C_L = 2 pi alpha / (1 + 2/AR)` = 0.37 at 5 degrees and AR 4; the case measures 0.33.
- Drag during the glide must be positive.

## A caution on the forces

The lift from `D/ib1_forces.dat` behaves sensibly here, but **the streamwise force on a thin inclined plate is
currently unreliable** — see issue #1849. The transverse component is unaffected. If you need drag or thrust
from a case like this one, cross-check it against a control-volume momentum balance rather than taking the
volume integral at face value.

## Running

```bash
./mfc.sh run examples/3D_ibm_flapping_plate/case.py -n 8 --gpu
```

0.48 M cells and 1 750 steps, which is three quarters of a flapping cycle: enough to see the leading-edge and
tip vortices form after onset, and to check the kinematics. Raise `t_end` for the approach to a periodic state.

## References

Prescribed kinematics of this form, and the Strouhal range, follow the animal-cruise literature; the canonical
pitch-ramp companion case is `examples/3D_ibm_pitchup_plate`.
