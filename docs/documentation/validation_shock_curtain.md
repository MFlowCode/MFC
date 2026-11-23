@page validation_shock_curtain Validation: Shock-Particle Interaction and Dense Curtains

# Shock-particle interaction: unsteady drag of a sphere and a shock through a dense curtain

**Status:** partially validated: the mean gas velocity through a dense curtain agrees with a published particle-resolved DNS to 0.04 to 0.05 \f$u_{ps}\f$ (3D, a reduced cross-section, two realizations); the single-sphere drag peak is within 3 % of the inviscid reference. **Code:** `master` at 0cdfdf02a. **Example:** `examples/2D_ibm_shock_curtain`.
**Last run:** 2026-11-23.

## Purpose

A shock crossing a cloud of particles loads each particle with an unsteady force and drives a mean gas flow through the cloud. This checks
both with fixed immersed bodies: the drag history of a single sphere after shock passage against the inviscid result of Mehta et al.
\cite Mehta2016, and the Favre-averaged gas velocity through a random curtain at a solid fraction of 0.2 against the particle-resolved DNS
of Sridhar, Fox and Capecelatro \cite Sridhar2025, including which wall model is needed. A short note records a qualitative comparison with
a shock-tube experiment.

## Setup

| Item | Value |
| :--- | :--- |
| Gas | air, \f$\gamma\f$ = 1.4, \f$p_1\f$ = 82.7 kPa, \f$T_1\f$ = 291.9 K; single-fluid ideal gas, inviscid; the post-shock state is the Rankine-Hugoniot jump, held at the inflow face (`bc_x%%beg = -17`) |
| Bodies | fixed, D = 115 um; sphere: one body at the shock's arrival point; curtain: random sequential addition, minimum gap 0.2 D, periodic across the flow, bodies kept off the curtain faces by half a diameter |
| Walls | inviscid no slip (`slip = F`, `viscous = F`: the ghost velocity enforces no slip) or inviscid slip (`--slip`) |
| Grid | sphere: 12 and 14 cells/D in [-3 D, 7 D] x [-3 D, 3 D]^2; curtain: 10 cells/D, 30 D x 6 D x 6 D (1.08 M cells, 239 spheres, curtain 17.4 D thick, \f$\phi\f$ = 0.2) |
| Numerics | WENO5-M, HLLC, SSP-RK3, `fd_order = 4`, adaptive step at CFL 0.4, `ib_ip_interp = 1` (the default) |
| Hardware | sphere runs: GNU 13 on an x86 CPU, one rank; curtain: NVIDIA A100 80GB, NVHPC 25.11; double precision |

## Reference data

- Mehta, Neal, Jackson, Balachandar and Thakur \cite Mehta2016, Fig. 6(a): the inviscid single-sphere drag at Ms 1.22, peak 7.49 at \f$t/\tau\f$ = 0.36 (\f$\tau = D/W_s\f$, \f$t = 0\f$ at the shock's arrival at the nose), read from the figure. The other peaks in the figure (3.72, 2.34, 1.71 and 1.44 at Ms 1.5, 2, 3 and 6) are tabulated in `check_curtain.py`.
- Sridhar, Fox and Capecelatro \cite Sridhar2025, Case 5 (Ms 1.66, \f$\phi\f$ = 0.2), Fig. 3(d)-(f): the Favre-averaged gas velocity over the post-shock velocity against \f$(x-x_0)/L\f$ at \f$t/\tau_L\f$ = 0.5, 1 and 2, \f$\tau_L = L/W_s\f$ with \f$L\f$ the curtain thickness. Their DNS is viscous with a 6th-order central scheme at 40 cells/D. The values were read from the published figure (raster panels, curve thickness about 0.01 in \f$u/u_{ps}\f$) and ship as binned medians in `examples/2D_ibm_shock_curtain/sridhar2025_ms166_phi02.txt`.

## Quantities

- **Cd(t)** = \f$F_x/(\tfrac12\rho_2 u_2^2 \pi D^2/4)\f$ with the post-shock state, from `D/ib_forces.dat` (`check_curtain.py --drag`).
- **ũ(x)** = \f$\langle I\rho u\rangle/\langle I\rho\rangle\f$ over the cross-stream slabs (I = 1 for cell centres outside every body), over \f$u_{ps}\f$. The comparison error is the mean \f$|\tilde u - \tilde u_\mathrm{ref}|/u_{ps}\f$ over \f$-0.35 < (x-x_0)/L < 1\f$, excluding \f$\pm 0.05 L\f$ around the reflected and transmitted shocks. Saves fall at \f$t/\tau_L\f$ within 0.25 of a reference panel (the nearest panel is used).

## Results

Single sphere, Ms 1.22, peak Cd (reference 7.49):

| Wall | cells/D | peak Cd | error | \f$t/\tau\f$ of the peak | mean over the last \f$\tau\f$ |
| :--- | ---: | ---: | ---: | ---: | ---: |
| inviscid slip | 12 | 7.636 | +1.9 % | 0.405 | 0.098 |
| inviscid slip | 14 | 7.724 | +3.1 % | 0.348 | 0.033 |
| inviscid no slip | 12 | 8.917 | +19.1 % | 0.404 | 0.848 |
| inviscid no slip | 14 | 8.955 | +19.6 % | 0.396 | 0.848 |

Curtain, Ms 1.66, 239 spheres at 10 cells/D, inviscid no slip, saves nearest each reference panel; mean absolute difference from the DNS profile in units of \f$u_{ps}\f$, and the mean velocity inside the curtain (0.1 < \f$(x-x_0)/L\f$ < 0.9):

| realization | \f$t/\tau_L\f$ (panel) | difference | mean \f$\tilde u/u_{ps}\f$ in the curtain |
| :--- | ---: | ---: | ---: |
| seed 1 | 0.479 (0.5) | 0.047 | 0.307 |
| seed 1 | 0.987 (1) | 0.050 | 0.575 |
| seed 1 | 2.002 (2) | 0.047 | 0.514 |
| seed 2 | 0.479 (0.5) | 0.038 | 0.312 |
| seed 2 | 0.987 (1) | 0.050 | 0.569 |
| seed 2 | 2.002 (2) | 0.052 | 0.505 |

The profiles of seed 1 are in `examples/2D_ibm_shock_curtain/curtain_profiles.png` (made by `check_curtain.py --png`): the transmitted and reflected shocks are where the DNS puts them to within a few percent of the curtain thickness, and the plateau inside the curtain is about 0.03 to 0.06 \f$u_{ps}\f$ below the DNS at \f$t/\tau_L\f$ = 1 and 2.

## Assessment

- The inviscid slip sphere lies within 2 to 3 % of Mehta's peak at 12 and 14 cells/D, with the late drag near zero as in the inviscid reference. The inviscid no-slip ghost velocity is not a substitute for it for a single body: its peak is 19 % high and its late drag 0.85, because no-slip with inviscid fluxes imposes a wall layer the reference does not have.
- Through the curtain the inviscid no-slip spheres reproduce the shape and the level of the DNS mean velocity profile to 0.04 to 0.05 \f$u_{ps}\f$ at the three panel times, for both realizations. This is a little larger than the 0.02 to 0.03 of earlier runs at \f$t/\tau_L\f$ = 2 (below); the plateau of these runs is slightly low, i.e. the curtain is slightly more opaque, which the reduced 6 D cross-section and the 10 cells/D resolution may explain (not verified).
- Qualitative comparison with a shock-tube experiment (Wagner et al. \cite Wagner2012, Beresh et al. \cite Beresh2011): earlier runs of a moving, colliding curtain with no-slip bodies followed the measured upstream and downstream curtain edges to about 150 microseconds after the shock's arrival (the downstream edge within the spread of the repeated experiments at 100 microseconds); later times, the transmitted and reflected wave strengths were not reconciled with the experiment. The published data are hard to reconcile with each other and with the stated conditions, so this is a comparison of trends, not a validation, and no benchmark is built on it. The example's `--moving` option is that configuration.

## Known limitations

- The curtain comparison is one realization on a 6 D x 6 D cross-section at 10 cells/D, not the 12 D x 12 D cross-section of the reference DNS and not at its 40 cells/D. Earlier runs on a pre-merge build with the full 12 D cross-section at 12, 20 and 40 cells/D gave errors of 0.094 to 0.096, 0.062 to 0.068 and 0.021 to 0.030 at \f$t/\tau_L\f$ = 0.5, 1 and 2, independent of resolution; an inviscid no-slip run at 12 cells/D gave 0.095, 0.066 and 0.023, and an inviscid slip run 0.124, 0.139 and 0.054 (gas through the curtain about 30 % too fast at \f$t/\tau_L\f$ = 1). The slip result is why the curtain uses no slip. These earlier numbers are not repeated on `master`.
- Shock positions: the reflected and transmitted shocks of the earlier runs trailed the reference by about 1 D; the definition of \f$x_0\f$ and \f$t = 0\f$ against the reference was not verified.
- Sridhar et al. use a power-law viscosity and Pr 0.7; the runs here are inviscid.
- Slip walls on moving bodies abort in the curtain (a NaN within the first steps): use no slip for moving bodies. This is a known solver limit documented with the slip-wall limits and is not fixed by this example; the registered test uses fixed bodies.
- The 3D sphere at Ms 3: the inviscid peak converged to 5.5 % above Mehta's and the late level was low in earlier runs (cause not verified); only Ms 1.22 is repeated here.

## Reproduce

```
./mfc.sh run examples/2D_ibm_shock_curtain/case.py -n 1 --gpu -- --dim 3 --single --slip --mach 1.22 --xu 3 --xd 7 \
    --width 6 --cells_per_D 12 --steps 0 --tau 8
python3 examples/2D_ibm_shock_curtain/check_curtain.py examples/2D_ibm_shock_curtain --drag --mach 1.22

./mfc.sh run examples/2D_ibm_shock_curtain/case.py -t pre_process simulation post_process -n 1 --gpu -- --dim 3 --h 17.4 \
    --width 6 --xu 7 --xd 23 --cells_per_D 10 --steps 0 --tau 2 --frames 20 --silo
python3 examples/2D_ibm_shock_curtain/check_curtain.py examples/2D_ibm_shock_curtain
```

The suite runs the CI-sized 2D case (12 cylinders, 240 x 96 cells, 300 steps, one rank) as `2D -> IBM -> Shock Curtain`.
