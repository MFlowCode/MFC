# Bubble Collapse in a Uniform Bubbly Liquid (Keller-Miksis)

A periodic, uniform bubbly mixture (`bubbles_euler`, one bubble size R0 = 50 um, void fraction 1e-5, water) is held at 1000 times
the bubble's equilibrium pressure, so every cell holds the same collapsing bubble. The bubble ODE (Keller-Miksis, polytropic gas,
surface tension, liquid viscosity) is sub-stepped adaptively (`adap_dt`). `check_collapse.py` integrates the same Keller-Miksis
equation in Python at a constant far-field pressure, and again driven by MFC's own cell pressure, and prints the collapse time,
the minimum radius and max |R_MFC - R_KM|/R0, with Rayleigh's collapse time 0.915 R0 sqrt(rho_l/(p_inf - p_v)) for scale.

    ./mfc.sh run examples/0D_bubblecollapse_adap/case.py -t pre_process simulation post_process -n 1
    python3 examples/0D_bubblecollapse_adap/check_collapse.py examples/0D_bubblecollapse_adap

`case.py` takes `--dt` (time step and save interval, default 1e-4 in R0/sqrt(p0/rho0)) and `--adap-tol` (default 1e-4).
Measured with GNU 13.3 (CPU, one rank):

| run | collapse time / t_Rayleigh | vs KM, constant p | max abs(R - R_KM)/R0, constant p | same, KM driven by MFC's p |
|---|---:|---:|---:|---:|
| default | 1.0425 | +7.2e-4 | 4.0e-3 | 8.9e-5 |
| `--adap-tol 1e-6` | 1.0425 | +7.2e-4 | 4.0e-3 | not run |
| `--dt 2.5e-5` | 1.0414 | +1.5e-4 | 4.2e-3 | 8.0e-5 |

The collapse time is the parabola through the three saves around the minimum, so it moves with the save interval for both MFC
and the ODE. The 4e-3 difference from the constant-pressure ODE is the two-way coupling: as the void collapses the mixture
pressure falls by 0.27 %, and the ODE driven by that pressure agrees to 9e-5 R0. An NVIDIA A100 run (NVHPC 25.11, OpenACC) gives
the same numbers to the digits shown. See @ref validation_bubbles.
