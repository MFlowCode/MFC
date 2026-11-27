# Freely Propagating H2/air Flame Speed Against Cantera

Stoichiometric H2/air at 300 K and 1 atm, `h2o2.yaml`, mixture-averaged transport, no Soret. `case.py` solves a Cantera `FreeFlame` (written to `IC/`,
read through `hcid = 170`), shifts it to x = 6.5 mm of an 8 mm domain and starts it in the frame of the unburnt gas. The flame moves toward x = 0 (a wall);
the outflow is a GRCBC subsonic boundary (`bc_x%end = -8`) with the outlet pressure held at 1 atm.

    ./mfc.sh run examples/1D_h2_air_flame_speed/case.py -t pre_process simulation post_process -n 1 -- --dx_um 20
    python3 examples/1D_h2_air_flame_speed/compare_flame_speed.py examples/1D_h2_air_flame_speed --case-args "--dx_um 20"

The case is not registered in the test suite (165 thousand steps at 40 um). Measured on master: dx = 40 um gives S_c = 2.3140 m/s (-0.76 % against the Cantera 2.3317), dx = 20 um gives 2.3268 m/s (-0.21 %). Results and limits: @ref validation_gas_chemistry.
