# PETN (Kuhl) JWL validation summary

Methodology follows the Blast-Afterburn validation checklist
(`numerical_checks_reactive_jwl.tex`, sections "Initiation and Propagation"
and "Breakout"): plot p(r,t), confirm a coupled shock-reaction structure
(not a decaying blast), confirm the front approaches the intended
detonation velocity, and characterize breakout.

## Parameters (Kuhl et al., UCRL-PROC-225822, 2006, Appendix)

    A = 5.8e11 Pa, B = 9.3e9 Pa, R1 = 7.0, R2 = 1.7, omega = 0.246
    rho0 = 1000 kg/m^3  (Kuhl's own pressed-charge reference density,
                          NOT full-density PETN TMD 1770 kg/m^3)
    Q = 5.9538e6 J/kg   (from the reported heat of detonation, 1423 cal/g)
    cv = 1174.9 J/(kg K) (from R_DP = 28.76 g/mol via cv = R/omega)
    jwl_delta_e = -7.98e6 J/kg (Garno Eq. 17 offset; closed-form, not fit)

## 1D calibration: `petn_1d_calibration.png`

Single-fluid planar detonation, dx = 10 um, direct initiation via a
rxn_val = 1 booster (a pressure-only hot spot was found to fizzle for this
material -- see below). Left panel: pressure profiles at successive times,
showing a von Neumann spike relaxing to a steady post-reaction plateau, not
a monotonically decaying pulse. Right panel: front trajectory (r vs t) is
linear (constant speed) once past the initial transient, fit to
D = 6085 m/s against the independently measured D_CJ = 6245 m/s (2.6%
difference, consistent with the fit window still containing some
relaxation).

**CJ eigenvalue check**: increasing `jwl_G` by > 3x (3.0e-14 to 1.0e-13)
left the fitted front speed and the post-reaction plateau pressure
unchanged to < 0.02% -- the operational definition of a converged CJ state
(the reaction rate sets zone width only, never propagation speed, once the
wave is truly self-sustaining). This measured pair (D_CJ = 6245 m/s,
P_CJ = 16.7 GPa) is used as the reference in the 3D octant case rather than
a hand-derived Hugoniot/Rayleigh tangency construction, which gave an
inconsistent result for this material's unusually steep coefficients
(R1 = 7.0, R2 = 1.7) and needs its own fix (flagged, not blocking).

## 3D octant spherical case: `probe_histories_dx1p0mm.png`

Radial probe pressure histories at r/Rc = 0.3, 0.6, 1.0 (charge edge), 1.5,
2.0. Confirms: (1) a coupled compression/reaction wave sweeps outward and
breaks out into the surrounding air -- not simply a decaying blast; (2) peak
pressures are of the right order relative to the measured P_CJ = 16.7 GPa.
Front-speed extraction from these probes is NOT yet trustworthy (see
`README.md`'s Status section for the two diagnostic issues -- coarse
ASCII probe time resolution and a threshold-based arrival detector that
mixes the true front with a resolution-dependent numerical precursor) and
should not be cited as a validated grid-convergence result yet.

## 2D reactive products-air mixing case

`examples/2D_jwl_reactive_air_mixing` was swapped from TNT to PETN_KUHL and
re-run end to end (400x400 mesh, full 4000-step run) with no NaN/abort.
Direct initiation again required a rxn_val = 1 booster rather than a
pressure-only hot spot. See `examples/2D_jwl_reactive_air_mixing/README.md`'s
"What the run shows" section for measured front/pressure numbers.

This case's `visualize.py` had a pre-existing bug (not introduced by this
PETN work): `load()` globbed only rank 0's `.dat` chunk under
`parallel_io=F`, silently plotting a small corner sliver of the domain
instead of the full field. Fixed by concatenating every rank's chunk before
reshaping. Also fixed: the products mass fraction panel showing `Y ~ 1`
across nearly the whole domain is real, not a residual bug -- PETN products
at Kuhl's reference density (1000 kg/m^3) outmass the swept air
(~1.2 kg/m^3) by ~3 orders of magnitude, so `Y` stays near 1 almost
everywhere the front has passed even after significant dilution.

The full Blast-Afterburn checklist (`numerical_checks_reactive_jwl.tex`,
all four sections) is run on this case by
`examples/2D_jwl_reactive_air_mixing/analyze_checklist.py`, which writes
`checklist_1_propagation.png` ... `checklist_4_closure.png`. Results:
(1) a coupled shock-reaction front with `λ` stepping 0→1 across it, not a
decaying blast; (2) explosive mass conserved to machine precision and total
energy rising x2.01 from the EOS-embedded chemical release (conservative,
no external source), both until the blast vents through the open boundary;
(3) charge breakout at ~7.5 us with the diverging front running below the
planar D_CJ ray (curvature velocity deficit); (4) the weighted-composition
closure active in a thin `0 < Y < 1` shell at the products/air contact.

## Open items

1. Fix or replace the hand-rolled analytic Hugoniot/Rayleigh-tangency CJ
   construction for this material (works for LX-10 in
   `examples/1D_jwl_znd_detonation/analyze_znd.py`, disagrees with the
   numerically measured CJ point here).
2. Build a lambda-crossing (reaction-completion) front tracker for the 3D
   octant case instead of a raw pressure threshold, then complete the
   dx = 0.25 mm grid level and report an actual observed order of
   convergence.
3. Increase ASCII probe output time precision (or write more frequently)
   so front-arrival fits are not quantized to ~1 us.
