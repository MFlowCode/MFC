# H2/air Constant-Volume Ignition Against Cantera

A uniform periodic box (26 cells, no fluxes) of stoichiometric H2/air at 1 atm ignites at constant volume. `compare_cantera.py` advances a
Cantera `IdealGasReactor` on the same mechanism (`h2o2.yaml`) and initial state and compares the ignition delay (time of max dT/dt), T(t),
the end state against the UV equilibrium, and the sum of the mass fractions.

    ./mfc.sh run examples/1D_h2_air_ignition/case.py -t pre_process simulation post_process -n 1 -- --T0 1100
    python3 examples/1D_h2_air_ignition/compare_cantera.py examples/1D_h2_air_ignition --case-args "--T0 1100"

`case.py` takes `--T0` (1100 K), `--dt` (1e-8 s), `--tend` (0.8 ms), `--nsaves` (200) and `--substeps` (`chem_params%reaction_substeps`, default 2).
The checker needs Cantera and the Silo output (`post_process`).

Measured with GNU 13 (CPU, one rank, double), `--nsaves 2000`, dt = 1e-8 s, `--substeps 2` (error = (MFC - Cantera)/Cantera for tau):

| T0 [K] | tau, MFC [us] | tau error | max abs(dT) [K] | T_end - T_eq [K] |
|---:|---:|---:|---:|---:|
| 1000 | 305.37 | +2.8e-6 | 0.29 | +0.060 |
| 1100 | 86.167 | +1.1e-5 | 0.29 | +0.051 |
| 1200 | 44.228 | +2.0e-5 | 0.28 | +0.043 |

The ignition delay is estimated from the saved samples (parabola through the maximum of dT/dt), so use many saves (`--nsaves 2000`, one save per 0.4 us) for a delay accurate to 1e-5; with 200 saves the absolute tau shifts by about 0.5 % (the MFC-Cantera difference is unaffected). The error column compares MFC and Cantera at the same sample times. The checker's default tolerances (tau 5e-5, max abs(dT) 1 K) are for dt = 1e-8 s.

Time-step study at T0 = 1100 K, `--substeps 2`: dt = 2e-8 s gives tau error +6.9e-6 (200 saves) and max abs(dT) = 1.13 K, and dt = 5e-8 s gives +4.5e-5 and 6.46 K, against 0.29 K at 1e-8 s: second order in dt. With `--substeps 0` (reaction source inside the Runge-Kutta stages) and dt = 5e-10 s at T0 = 1200 K the ignition delay agrees to 2e-11 and max abs(dT) < 1e-4 K. That integrator is stable only for dt below 2.51/(largest Jacobian eigenvalue), about 6e-10 s for burnt H2/air.

The suite registers `1D -> Chemistry -> H2-Air Ignition -> Split Reactions` (T0 = 1200 K, dt = 2e-8 s, 100 us). See @ref validation_gas_chemistry.
