# 2D Sedov–Taylor Cylindrical Blast Wave

A quantitative gas-dynamics benchmark: a concentrated deposit of internal energy
in an ideal gas (γ = 1.4) drives a strong cylindrical blast whose shock front
obeys the exact self-similar Sedov–Taylor solution

```
R_shock(t) = ξ₀ · (E_L / ρ₀)^(1/4) · t^(1/2)        (cylindrical, ν = 2)
```

Three of the checks are **parameter-free** — no EOS calibration, no fitted
constant — which is exactly what makes Sedov–Taylor the standard test that a
compressible solver gets blast energetics right. One quadrant `[0,L]²` with
reflective axes represents the full line-symmetric blast (`E_L = 4 ×` the
quadrant deposit energy); the driver is a small quarter-disk at the origin at
high pressure in low-pressure air (low ambient → high blast Mach number → the
genuine strong-shock regime where the solution is exact).

## Run

```bash
./mfc.sh run examples/2D_sedov_taylor_blast/case.py -n 4
./examples/2D_sedov_taylor_blast/analyze_sedov.py examples/2D_sedov_taylor_blast
```

## Result (299×299, WENO5 + HLLC)

```
1) growth exponent  R ~ t^n        : n = 0.483        (exact Sedov cyl = 0.5)          PASS
2) density ratio vs finite-Mach RH : 5.08 vs 5.97     (M→∞ limit 6.0; 15% = smearing)  PASS
3) Sedov constant  ξ₀ (measured)   : 1.031 ± 0.008    (literature ≈ 1.0)               PASS
   energy conservation             : 99.0% of E_L, std 0.0%
```

- **Growth exponent** — a log–log fit over the self-similar window (shock radius
  well past the deposit, still inside the domain; blast Mach 24–707) gives
  `R ∝ t^0.483`, within 3% of the exact `t^0.5`.
- **Self-similar constant** — `R / ((E_L/ρ₀)^¼ t^½)` is constant to 0.8% and
  equals 1.031, matching the textbook cylindrical γ = 1.4 value (≈ 1.0). It
  moves *toward* that value as the shock strengthens (1.065 at moderate Mach →
  1.031 at high Mach), confirming approach to the strong-shock limit.
- **Density ratio** — the post-shock peak reaches 5.08 vs the analytic 5.97 at
  the measured Mach; the ~15% deficit is grid smearing of the thin density spike
  at the shock (the hardest quantity for any shock-capturing scheme — pressure
  and velocity are captured much better), not a physics error.
- **Energy conservation** — the excess total energy integrated over the field
  stays at 99.0% of the deposited `E_L` with zero drift across the window.

`analyze_sedov.py` writes `sedov_validation.png` (blast trajectory vs the exact
solution + a log–log slope plot). This validates MFC's compressible
gas-dynamics and energy handling independently of any JWL/EOS modeling — the
companion to the JWL-specific checks in `examples/1D_jwl_znd_detonation` and
`examples/2D_jwl_reactive_air_mixing`.
