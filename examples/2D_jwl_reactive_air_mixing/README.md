# 2D JWL Reactive Detonation with Full Products–Air Mixing (TNT)

A real-world 2D free-air burst: a full-density cylindrical TNT charge
(`rho0 = 1630 kg/m^3`, 5 cm radius) at rest at ambient pressure in air at STP.
A small central hot spot (8 mm, 25 GPa) ignites a self-propagating JWL++
(`jwl_reactive`, Souers 2000) detonation, `dλ/dt = G p^b (1-λ)`, that sweeps
outward through the charge, then the products expand violently into the
surrounding air.

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

The TNT JWL fit (`A`, `B`, `R1`, `R2`, `ω`, `E0`) matches
`examples/2D_jwl_detonation`; `jwl_delta_e ≈ -1.29e7 J/kg` is computed in
`case.py` so unreacted TNT at `(rho0, e=0)` sits exactly at ambient pressure.

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

## What the run shows

- **t ≈ 4 µs** — detonation propagating outward, `Pmax ≈ 18 GPa` (near TNT's
  ~21 GPa CJ pressure), front at ~29 mm; the diverging cylindrical front runs
  slightly below the planar CJ speed as expected from front curvature.
- **t ≈ 34 µs** — the whole charge has burnt (`λ→1`), the detonation has broken
  out of the charge, and the blast has decayed to ~0.7 GPa while a partially
  mixed `0.05 < Y < 0.95` layer develops at the products/air contact.
- **t ≈ 60 µs** — far-field air blast (~0.1 GPa) with a developed,
  Richtmyer–Meshkov-fingered mixing layer.

The domain (40 cm square) is sized so the blast reaches the non-reflecting
outer boundaries near the final time; shorten `t_step_stop` or enlarge the
domain to keep the front interior for longer. A uniform 400×400 mesh is used
deliberately — grid stretching is known to go NaN with strong JWL reactive
sources (see the note in `examples/2D_jwl_detonation`).
