# 1D JWL ZND Detonation (LX-10-0)

A self-propagating JWL++ (`jwl_reactive`) detonation with the reactant/product
energy offset `fluid_pp(i)%jwl_delta_e` enabled, so the wave develops genuine
ZND structure: a von Neumann pressure spike decaying through a finite reaction
zone to the Chapman-Jouguet state, followed by a Taylor rarefaction.

The explosive is LX-10-0 with the JWL products parameters, detonation energy,
and CJ reference values of Garno et al., *J. Appl. Phys.* **128**, 195903
(2020), Table I (D_CJ = 8.821 km/s, P_CJ = 37.5 GPa). `jwl_delta_e` is
computed in `case.py` so that unreacted material (lambda = 0) at
(rho0, e = 0) sits at ambient pressure -- Garno's Eq. (17) offset. A 1 mm,
45 GPa hot spot at the left wall ignites the charge; `jwl_reactive`
initializes lambda = 0 everywhere else.

Run and check against the exact analytic Hugoniot / Rayleigh-line / CJ
construction (computed independently of the solver):

```bash
./mfc.sh run examples/1D_jwl_znd_detonation/case.py -n 4
./examples/1D_jwl_znd_detonation/analyze_znd.py examples/1D_jwl_znd_detonation 2.0e-10
```

`analyze_znd.py` verifies that (1) the measured front speed matches the
analytic D_CJ, (2) the pressure at reaction completion matches the analytic
P_CJ, and (3) the peak pressure is a real spike above P_CJ, approaching this
model's analytic von Neumann pressure from below as the grid is refined (the
spike is the least-resolved feature; at the default 10 micron grid it captures
part of the vN excursion). It also writes `znd_profile.png` with the p(x)
history and the Hugoniot/Rayleigh construction.

Note the distinction between this model's exact vN pressure (~75 GPa, from the
offset-JWL unreacted Hugoniot) and measured LX-10 spike values (~55 GPa): the
single-fluid offset makes the unreacted Hugoniot a stiffened products curve,
not a calibrated reactant EOS, so the spike height is model-specific. The CJ
speed and pressure, by contrast, are properties of the products EOS + Q alone
and match the calibrated values.
