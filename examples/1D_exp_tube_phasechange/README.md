# Expansion Tube with Phase Change

Liquid water, water vapor and air (stiffened-gas parameters of Le Metayer et al.) at 1.175e5 Pa and the saturation temperature,
pulled apart at 2 m/s from the centre, with infinitely fast pressure-temperature-Gibbs relaxation (`relax_model = 6`) after each
step. `check_saturation.py` reads the relaxed state and prints, in every cell, the temperature of each phase from its own EOS
(thermal equilibrium) and the distance of T from T_sat(p), the root of g_liquid = g_vapor (chemical equilibrium).

    ./mfc.sh run examples/1D_exp_tube_phasechange/case.py -t pre_process simulation post_process -n 1 -- -me 3
    python3 examples/1D_exp_tube_phasechange/check_saturation.py examples/1D_exp_tube_phasechange --case-args "-me 3"

`-me 3` (6-equation model, default) or `-me 2` (5-equation). Measured with GNU 13.3 (CPU, one rank), last save (t = 511 us):

| model | max abs(T_k/T_v - 1) | max abs(T - T_sat(p)) [K] | vapor mass change |
|---|---:|---:|---:|
| 6-equation | 9.4e-7 | 7.0e-4 | +6.76 % |
| 5-equation | 9.4e-7 | 6.2e-4 | +6.76 % |

The initial state lies on the saturation curve to 2.5e-7 K. The fluid parameters give T_sat(1 atm) = 373.79 K (373.12 K
measured). Only the equilibrium of the relaxed state is checked, not the wave structure. See @ref validation_multiphase_riemann.
