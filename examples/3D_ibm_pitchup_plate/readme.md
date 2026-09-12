# 3D Canonical Pitch-Up Plate (validation against published DNS)

## Case Set Up

A rectangular flat plate of aspect ratio 4 pitches about its **leading edge** from 0 to 45 degrees, following the
smoothed linear ramp adopted by the AIAA Fluid Dynamics Technical Committee's Low Reynolds Number Discussion
Group as a canonical maneuver. Half the span is simulated with a symmetry plane at the root and a free tip. The
plate holds at 45 degrees after the ramp.

| | |
|---|---|
| Chord Reynolds number | 300 |
| Mach number | 0.2 |
| Reduced pitch rate `K = Omega c / 2 U` | pi/8 — the plate pitches over one chord of travel (case *C1*) |
| Eldredge smoothing `a` | 21 |
| Pivot | leading edge |

The pitch angle follows

```
alpha(t) = (Omega / 2a) log[ cosh(a(t - t1)) / cosh(a(t - t2)) ] + alpha_max / 2
```

which is Eq. (1) of Jantzen et al. (2014). It is supplied by `patch_ib(1)%kin_model = 2`, so the ramp is set at
run time from `kin_pitch_rate`, `kin_smooth`, `kin_theta0` and `kin_t0` rather than compiled into the case. The
same binary therefore runs the whole `C1`/`C2`/`C4`/`C6` family, which differ only in pitch rate and smoothing.

The plate is defined to reach 0.5 c past the symmetry plane, so its root is a continuation of the wing rather
than a no-slip face — an immersed boundary marks only cells inside the body, and a plate ending exactly on the
plane would see fluid on the far side of its root.

## Numerics

WENO5 with HLLC and the Thornber low-Mach velocity correction, `mp_weno` off, characteristic inflow and outflow
with the generalized relaxation treatment, fourth-order finite differences for the immersed boundary. 80 cells
per chord with the section 0.05 c thick, so it holds four cells across.

**The four-cell rule matters and is the reason for those numbers.** An immersed body reconstructs its interior
from the fluid outside; with too few cells across there is no interior left. A companion 2D study of this same
maneuver, holding the physical thickness fixed at 2.5 percent of chord and varying only the mesh, found the lift
through the ramp 30 to 40 percent high at two cells across and converged by four. Two cells *appeared* to agree
better with the reference than four did, because two errors cancelled — so do not tune this by the agreement.

## Validation

Reference: Jantzen, Taira, Granlund & Ol (2014), Fig. 10, AR 4 panel, curve C1 — an immersed-boundary projection
DNS at the same Reynolds number, with a plate of zero thickness.

**What agrees.** The pitch history matches the closed-form ramp to 1e-11. The lift reaches 76 percent of the
reference after the ramp. The residual is consistent with the 2D study, which found a converged offset of about
25 percent and attributes it to the finite, square-edged section modelled here against the reference's
infinitely thin plate.

**What does not.** The streamwise force is wrong. A plate held at 45 degrees carries a force essentially normal
to its surface, so `C_D / C_L` should be about 1. It is at the instant the rotation stops, and then collapses:

| time after ramp start | C_D / C_L |
|---|---|
| 1.0 (ramp end) | 1.04 |
| 2.0 | 0.11 |
| 4.0 | 0.10 |

A control-volume momentum balance on the saved fields — evaluated on grid planes away from the body, so
independent of the immersed-boundary reconstruction — gives `C_D = 1.04` where the volume integral reports 0.12,
while the two agree on lift. At half this resolution the reported drag goes negative, which is impossible. This
is a defect in the force diagnostic for thin inclined plates, tracked as issue #1849; the `C_D / C_L` ratio after
the ramp is a reproducer that needs no reference data at all.

So: **use this case to validate lift and kinematics; do not take its drag at face value** until #1849 is
resolved.

## Running

```bash
./mfc.sh run examples/3D_ibm_pitchup_plate/case.py -n 8 --gpu
```

34 M cells, 9 600 steps — roughly 40 minutes on eight MI250X GCDs. `ib_neighborhood_radius` is set explicitly
because the plate spans more than one rank's subdomain at larger decompositions, where leaving it automatic
faults.

To run the rest of the family, change `K` and the matching smoothing `a` (21, 16, 11, 4 for C1, C2, C4, C6);
no rebuild is needed.

## References

Jantzen, R. T., Taira, K., Granlund, K. and Ol, M. V., "Vortex dynamics around pitching plates,"
*Physics of Fluids* **26**, 053606 (2014). doi:10.1063/1.4879035

Eldredge, J. D., Wang, C. and Ol, M. V., "A computational study of a canonical pitch-up, pitch-down wing
maneuver," AIAA Paper 2009-3687 (2009). doi:10.2514/6.2009-3687

Ol, M. V., Altman, A., Eldredge, J. D., Garmann, D. J. and Lian, Y., "Resume of the AIAA FDTC Low Reynolds
Number Discussion Group's canonical cases," AIAA Paper 2010-1085 (2010). doi:10.2514/6.2010-1085
