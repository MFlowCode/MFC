# 2D JWL Reactive Detonation with Full Products–Air Mixing (PETN)

A real-world 2D free-air burst: a cylindrical PETN charge (`rho0 = 1000 kg/m^3`,
5 cm radius, Kuhl et al.'s own pressed-charge reference density — see
`examples/3D_jwl_znd_spherical_octant/VALIDATION.md` for the full parameter
provenance) at rest at ambient pressure in air at STP. An 8 mm booster patch,
seeded already-reacted (`rxn_val = 1`, not just pressurized — see below),
ignites a self-propagating JWL++ (`jwl_reactive`, Souers 2000) detonation,
`dλ/dt = G p^b (1-λ)`, that sweeps outward through the charge, then the
products expand violently into the surrounding air.

Direct initiation with a pressure-only hot spot (no `rxn_val`) was tried first
and fizzles for this material: PETN_KUHL's CJ pressure (≈ 7.35 GPa, see
below) is low enough that a pressure spike alone decays without reaching a
self-sustained reaction wave. Seeding the booster patch as already-reacted
gives robust initiation instead.

This case exercises two things the shipped program-burn example
(`examples/2D_jwl_detonation`) does not:

1. **Genuine reactive burn with a ZND offset.** `fluid_pp(1)%jwl_delta_e`
   (Garno et al. 2020, Eq. 17) puts the unreacted explosive on a stiffer
   Hugoniot via `e_eff = e + Y(1-λ)Δe`, so the front carries a von Neumann
   pressure spike rather than a monotonic energy source. The reaction is gated
   by `Y(1-λ)` — only unreacted **explosive** burns; the air never reacts, and
   the `Y` factor in the offset keeps pure air cells on their own energy. (The
   sub-mm reaction zone is only partially resolved at this 1 mm engineering
   grid; the fully resolved spike-vs-analytic validation is the 1D companion
   `examples/1D_jwl_znd_detonation`.)

2. **Full products–air mixing.** Once the detonation reaches the charge edge,
   the products (`Y=1`) drive into the air (`Y=0`) and the whole `0 < Y < 1`
   range of the weighted-composition mixture closure is exercised across a
   Richtmyer–Meshkov-unstable contact — the mixture EOS's real workload, not
   just the pure endpoints.

The PETN JWL fit (`A`, `B`, `R1`, `R2`, `ω`, `E0`) is `PETN_KUHL` from
`toolchain/mfc/jwl_products.py`, extracted from Kuhl et al. (UCRL-PROC-225822,
2006, Appendix); `jwl_delta_e ≈ -7.98e6 J/kg` is computed in `case.py` (Garno
et al. 2020, Eq. 17, closed-form, not fit) so unreacted PETN at `(rho0, e=0)`
sits exactly at ambient pressure. `jwl_G = 1.0e-13` was calibrated on the 1D
companion case via a CJ-eigenvalue (G-invariance) test — see
`examples/3D_jwl_znd_spherical_octant/VALIDATION.md`.

## Run

```bash
./mfc.sh run examples/2D_jwl_reactive_air_mixing/case.py -n 4
./examples/2D_jwl_reactive_air_mixing/visualize.py examples/2D_jwl_reactive_air_mixing
```

`visualize.py` writes two figures from the `D/` output:

- `contours_final.png` — pressure (log), products mass fraction `Y`, reaction
  progress `λ`, and a synthetic schlieren at the final save. The schlieren shows
  the fingered, unstable products/air interface.
- `contours_evolution.png` — pressure at four times: initiation → outward
  detonation → products-air expansion.

`analyze_checklist.py` runs four quantitative validation checks by binning the
field onto radius `r = hypot(x, y)`:

- `checklist_1_propagation.png` — radial profiles `p, ρ, u_r, λ, Y` at six times
  (initiation and propagation: a coupled shock-reaction front, not a decaying
  blast; `λ` steps 0→1 across the front).
- `checklist_2_energy.png` — explosive mass (flat to machine precision until the
  blast reaches the open boundary), unreacted mass `∫ρY(1-λ)dV → 0`, and total
  energy (rises ×2.01 from the EOS-embedded chemical release, then flat until
  efflux — a conservative release, no external source).
- `checklist_3_breakout.png` — shock trajectory `R_s(t)` with charge breakout,
  the peak-pressure envelope `p_shock(r)` fitted over the uncontaminated window,
  and the radial specific impulse `I(r)`.
- `checklist_4_closure.png` — the mixing cells `1e-6 < Y < 1-1e-6` and the
  `(ρ, p, e)` state there, i.e. where the weighted-composition closure is active.

## What the run shows

Measured directly from the `D/` output of the shipped run (400×400, 1 mm cells):

- **t ≈ 4 µs** — detonation propagating outward, `Pmax ≈ 9.1 GPa` near the
  booster/charge center, front (`p > 2 atm`) at ~33 mm; the diverging
  cylindrical front runs below the self-sustained planar CJ state
  (`D_CJ ≈ 5384 m/s`, `P_CJ ≈ 7.35 GPa`, computed from the JWL EOS; Kuhl
  reports only `T_CJ ≈ 4600 K`) as expected from front curvature.
- **t ≈ 22–37 µs** — the charge has burnt through, the detonation has broken
  out of the charge, and the front (`p > 2 atm`) has advanced to ~116–185 mm
  while `Pmax` decays through ~1.8 GPa → 0.4 GPa as the blast disperses into
  air.
- **t ≈ 60 µs** — the front (`p > 2 atm`) reaches ~283 mm, i.e. the domain
  corner — the blast has essentially exited through the non-reflecting
  boundaries by the final save, `Pmax ≈ 0.075 GPa`.
- **Products mass fraction `Y`** is close to 1 over nearly the whole swept
  region by the final save. This is a genuine consequence of the density
  contrast, not a diffusion artifact: PETN products at Kuhl's reference
  density (`rho0 = 1000 kg/m^3`) outmass the swept air (`~1.2 kg/m^3`) by
  ~3 orders of magnitude, so even after the products expand and dilute by
  1-2 orders of magnitude within the domain, `Y = c1/(c1+c2)` stays near 1
  almost everywhere the front has passed. A visibly partial `0 < Y < 1`
  mixing layer is confined to a thin band at the products/air contact
  (see the schlieren panel).

The domain (40 cm square) is sized so the blast reaches the non-reflecting
outer boundaries near the final time; shorten `t_step_stop` or enlarge the
domain to keep the front interior for longer. A uniform 400×400 mesh is used
deliberately — grid stretching is known to go NaN with strong JWL reactive
sources (see the note in `examples/2D_jwl_detonation`).

`visualize.py`'s `load()` concatenates every rank's `.dat` chunk
(`parallel_io=F` writes one file per MPI rank) before reshaping onto the
grid — reading only rank 0 silently returns a small corner sliver of the
domain rather than an error, so always run with the fixed loader.
