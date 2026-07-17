# PETN (Kuhl) JWL validation summary

Methodology: plot p(r,t), confirm a coupled shock-reaction structure
(not a decaying blast), compare the front against the material's true
Chapman-Jouguet (CJ) state, and characterize breakout.

## Parameters and provenance

    A = 5.8e11 Pa, B = 9.3e9 Pa, R1 = 7.0, R2 = 1.7, omega = 0.246
    rho0 = 1000 kg/m^3   Q = 5.9538e6 J/kg   cv = 1174.9 J/(kg K)
    jwl_delta_e = -7.98e6 J/kg (Garno Eq. 17 offset; closed-form, not fit)

A, B, R1, R2, omega, rho0 are quoted verbatim by Kuhl et al.,
UCRL-PROC-225822 (2006), p. 6. Q = 1423 cal/g (theoretical heat of
detonation) and T_CJ ~ 4600 K are from Kuhl & Khasainov, 38th ICT (2007).
rho0 = 1 g/cc is the actual loading density of Kuhl's pressed PETN booster
spheres, not full-density PETN TMD (1770 kg/m^3); the constants are only
self-consistent at 1 g/cc.

Two honesty notes on what these constants are:

- They are a low-density CHEETAH fit to the **detonation-products expansion
  isentrope** (the CJ adiabat), which Kuhl uses to model products that are
  **already detonated** and are expanding and afterburning with air. Kuhl
  does not self-propagate a detonation; he imposes a Taylor CJ initial
  condition. These MFC examples instead repurpose the same constants to
  drive a self-propagating reactive burn (`jwl_reactive`) -- a legitimate
  use of the products EOS, but not what Kuhl did.
- Kuhl's own EOS carries an ideal-gas thermal term `rho R_DP T`
  (R_DP = 28.76 g/mol, gamma = 1 + omega = 1.246) in place of the standard
  JWL `omega E / v` energy term, and his listed C constant is unused. MFC
  uses the standard energy-form JWL. So MFC borrows Kuhl's A/B/R1/R2/omega
  and 1 g/cc reference, but is not a bit-for-bit reproduction of his
  thermodynamics.

## The material's true CJ state (computed, not reported by Kuhl)

Kuhl reports only T_CJ ~ 4600 K -- **no P_CJ, no D_CJ, no CJ density**. The
self-sustained CJ state of the standard-form JWL with these constants, from
the Hugoniot/Rayleigh tangency (rate-independent, EOS only):

    D_CJ ~ 5384 m/s     P_CJ ~ 7.35 GPa     P_vN ~ 13.6 GPa (offset model)

This is trustworthy: the same tangency construction reproduces literature CJ
for full-density PETN (33.5 GPa), TNT (23 GPa) and Comp B (29.5 GPa) to a
few percent, and empirical low-density PETN at 1 g/cc gives D ~ 5.5 km/s,
P_CJ ~ 7.8 GPa -- both bracket the computed pair.

## 1D calibration: `petn_1d_calibration.png` -- an OVERDRIVEN front

Single-fluid planar detonation, dx = 10 um, direct initiation via a 30 GPa
`rxn_val = 1` booster (a pressure-only hot spot fizzles for this material).
The captured front is steady but **overdriven**, not at CJ: peak ~19.8 GPa,
a ~16.8 GPa state just behind, propagating at D ~ 6245 m/s. That pair is
exactly the overdriven lambda = 1 Hugoniot/Rayleigh strong-branch state at
D = 6245 m/s (16.7 GPa computed), i.e. arithmetic proof it is overdriven,
not CJ. The overdrive is driven by the oversized booster plus the closed
rear boundary acting as a sustaining piston; a planar wave with no rear
rarefaction cannot relax to CJ.

An earlier note here claimed a "G-invariance -> CJ" result (front speed and
plateau pressure unchanged under a >3x change in `jwl_G`). That argument was
wrong: a piston-supported overdriven front is **also** G-independent (its
speed is set by the piston, not the reaction rate), so G-invariance does not
distinguish CJ from overdriven. The correct reference is the EOS tangency
above (D_CJ ~ 5384 m/s, P_CJ ~ 7.35 GPa), which agrees with the earlier
hand-rolled analytic construction (~5450 m/s, ~6.3 GPa) to within its
bracket -- that construction was right; the 16.7 GPa "measured CJ" was the
outlier.

To actually calibrate to CJ, soften the booster and give the rear a
rarefiable (outflow) boundary so a Taylor wave forms and the front relaxes
to the self-sustained state. This is the recommended follow-up; the shipped
1D run demonstrates a coupled ZND structure but at an overdriven speed.

## 2D and 3D cases -- the diverging front DOES relax

In 2D and 3D the geometry diverges and vents, so the overdrive relaxes.
`examples/2D_jwl_reactive_air_mixing/analyze_checklist.py` runs a four-part
check (writes `checklist_1_propagation.png` ... `checklist_4_closure.png`):

1. **Propagation**: a coupled shock-reaction front with lambda stepping 0->1
   across it, not a decaying blast. Near breakout the front is mildly
   overdriven from the booster, then falls **below** the self-sustained
   planar CJ ray (D ~ 5.4 km/s) as curvature/divergence take over -- the
   physical velocity deficit of a curved front.
2. **Energy**: explosive mass conserved to machine precision; total energy
   rises x2.01 from the EOS-embedded chemical release (conservative, no
   external source); both hold until the blast vents through the open
   boundary.
3. **Breakout**: charge breakout at ~7.5 us; near-field peak-pressure decay
   fitted over the uncontaminated window.
4. **Closure**: the weighted-composition closure active in a thin 0 < Y < 1
   shell at the products/air contact.

The 3D octant probe histories (`probe_histories_dx1p0mm.png`, r/Rc = 0.3,
0.6, 1.0, 1.5, 2.0) confirm a coupled compression/reaction wave that breaks
out into air. Front-speed extraction from those probes is NOT yet
trustworthy (coarse ASCII probe time resolution; a threshold arrival
detector that mixes the true front with a resolution-dependent numerical
precursor) and is not cited as a converged grid-convergence result.

`visualize.py` in the 2D case had a pre-existing bug: `load()` read only
rank 0's `.dat` chunk under `parallel_io=F`, plotting a corner sliver of the
domain. Fixed by concatenating all ranks. The `Y ~ 1` products fraction over
most of the domain is real, not a bug -- PETN products at 1 g/cc outmass the
swept air by ~3 orders of magnitude, so Y stays near 1 wherever the front
has passed even after dilution.

## Verdict

The **solver physics is correct**: a genuine coupled shock-reaction
detonation with the right energy release and closure behavior. The earlier
**documentation was not**: it labeled an overdriven 1D state (16.7 GPa,
6245 m/s) as the "measured CJ" and attributed it as if it were a material
property. Corrected here -- the true CJ is ~7.35 GPa / ~5384 m/s (computed),
Kuhl reports no CJ pressure, and the 1D run is overdriven by its booster.

## Open items

1. Re-calibrate the 1D case to the self-sustained CJ (softer booster,
   rarefiable rear boundary) so it demonstrates CJ, not an overdriven front.
2. Build a lambda-crossing front tracker for the 3D octant case, complete
   the dx = 0.25 mm grid level, and report an observed order of convergence.
3. Increase ASCII probe output time precision so front-arrival fits are not
   quantized to ~1 us.
