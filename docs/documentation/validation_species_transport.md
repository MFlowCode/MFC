@page validation_species_transport Validation: Species Diffusion, Heat Conduction and the Transport Time Step

# 1D species diffusion, conduction and `chem_params%%transport_dt`

**Status:** partially validated. **Code:** `master` at 0cdfdf02a. **Example:** `examples/1D_diffusion_verification`.
**Last run:** 2026-11-23.

## Purpose

Fuel, oxidizer and heat reach a reacting surface by mixture-averaged species diffusion (with the correction velocity), Fourier
conduction with the species-enthalpy flux, and the explicit transport time-step limit. This checks each against an exact or an
independent solution: the erfc solution of a step for trace-species diffusion and for conduction (order of accuracy), an
independent 1D solver for a multicomponent mixing layer (the face composition of the correction flux was made second order in
[#2205](https://github.com/MFlowCode/MFC/pull/2205)), and the stability threshold of `chem_params%%transport_dt` for both transport models.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | 1D, uniform grid, 1 atm, gas at rest, reactions and viscosity off, `h2o2.yaml` (10 species) |
| Trace species | H2 at mass fraction 1e-6 in N2 (left half) against pure N2, 300 K, L = 1 mm; the step starts on a cell face; ends when \f$\sqrt{4Dt}\f$ = 0.05 L |
| Conduction | pure N2, 300.01 K against 300 K, same domain and end time rule |
| Mixing layer | x < 1 mm: X(H2, N2) = (0.3, 0.7) at 300 K; x > 1 mm: X(O2, H2O, N2) = (0.2, 0.2, 0.6) at 1000 K; L = 3 mm; wall at x = 0, extrapolation at x = 3 mm; outputs at 10 and 20 us; 150-1200 cells |
| Stability | periodic, 64 cells of \f$\Delta x = D_\mathrm{max}/(10c)\f$ (the transport limit is ten times tighter than the acoustic one), 0.1 K spike in one cell, 400 steps at the transport-limited step |
| Numerics | five-equation model, WENO5 (mapped, `mp_weno`), HLLC, SSP-RK3; fixed \f$\Delta t = 0.4\,\Delta x/c\f$ (erf, mixing layer), `cfl_adap_dt` with `chem_params%%transport_dt` (stability); `chem_params%%transport_model` 1 (mixture-averaged) unless stated |
| Hardware | GNU 13, CPU, one rank, double precision |

## Reference data

- **Erfc solution.** A step of height A diffusing with a constant coefficient: \f$(A/2)\,\mathrm{erfc}((x-x_0)/\sqrt{4Dt})\f$ \cite Crank1975,
  averaged over each cell with the antiderivative of erfc. D is Cantera's mixture-averaged value for H2 in the N2 state (7.790e-5 m^2/s;
  the mixture-averaged and binary values agree in the trace limit) or the thermal diffusivity \f$\lambda/(\rho c_p)\f$ of N2 (2.241e-5 m^2/s).
- **Mixing-layer reference.** `examples/1D_diffusion_verification/ref_mixing_layer.py` is an independent solver for isobaric multicomponent
  diffusion in the Lagrangian mass coordinate \f$\psi = \int \rho\,dx\f$: \f$dY_k/dt = -\partial j_k/\partial\psi\f$ and
  \f$c_p\,dT/dt = -\partial q/\partial\psi + \sum_k h_k\,\partial j_k/\partial\psi\f$, with \f$j_k = -\rho D_k (W_k/W)\,\partial X_k/\partial x + Y_k \sum_j \rho D_j (W_j/W)\,\partial X_j/\partial x\f$ and
  \f$q = -\lambda\,\partial T/\partial x + \sum_k h_k j_k\f$, Cantera mixture-averaged properties (Soret off), fixed-mass finite volumes with central
  face averages and SSP-RK3. It shares no code with MFC.
- **Time-step limit.** The explicit SSP-RK3 stability region on the negative real axis is \f$|z| \le 2.513\f$, and the two-point diffusion operator has the largest
  eigenvalue \f$4D/\Delta x^2\f$, so the diffusion number \f$D\,\Delta t/\Delta x^2\f$ must stay below 0.628. With `chem_params%%transport_dt` the `cfl_target` is that number on
  \f$D_\mathrm{max}\f$, so the threshold is expected at 0.628. Derived, not measured.

## Quantities

- **L1, Linf (erf).** \f$L_1 = \sum |e_i|\,\Delta x/(A s)\f$ and \f$L_\infty = \max|e_i|/A\f$ at the last save, with e the MFC cell value minus the exact cell average, A the step (1e-6 in Y or 0.01 K) and \f$s = \sqrt{4Dt}\f$ (`check_diffusion.py`).
- **Mixing layer.** Errors at 10 and 20 us in the mass coordinate (the MFC \f$\psi_i = \sum_{j<i}\rho_j\Delta x + \rho_i\Delta x/2\f$), reference interpolated linearly in \f$\psi\f$:
  \f$L_1 = \overline{|e|}\,(3\ \mathrm{mm})/(|\Delta Y_k|\,(1\ \mathrm{mm}))\f$ for each species and for T (jump 700 K). MFC's T is first corrected isentropically to the reference
  pressure, because the extrapolation outflow leaves its pressure up to 3.4e-4 above 1 atm (`check_mixing_layer.py`).
- **Growth.** The amplitude of the grid-scale (Nyquist) mode of T, \f$|\sum_i (-1)^i (T_i - \bar T)|/N\f$, at the end over its first value; above 1 is unstable.

## Results

**Trace H2 and conduction against the erfc solution** (front 1.6 to 25.6 cells wide):

| N | species L1 | order | Linf | heat L1 | order | Linf |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 7.07e-2 | | 3.64e-2 | 7.07e-2 | | 3.64e-2 |
| 64 | 1.58e-2 | 2.16 | 8.43e-3 | 1.58e-2 | 2.16 | 8.44e-3 |
| 128 | 3.85e-3 | 2.04 | 2.08e-3 | 3.84e-3 | 2.04 | 2.08e-3 |
| 256 | 9.57e-4 | 2.01 | 5.16e-4 | 9.51e-4 | 2.01 | 5.18e-4 |
| 512 | 2.39e-4 | 2.00 | 1.29e-4 | 2.33e-4 | 2.03 | 1.31e-4 |

Species diffusion and conduction use the same two-point operator and agree to three digits at equal cells per front width.

**Mixing layer against the reference** (20 us, mass coordinate; reference of 3000 cells, which a 1500-cell run of the same solver differs from by an L1 of 1.8e-5 for H2, 1.1e-5 for O2, 1.5e-5 for H2O, 1.2e-5 for N2 and 2.7e-6 for T, so the 3000-cell error is about a quarter of that):

| N | dx [um] | L1 H2 | L1 O2 | L1 H2O | L1 N2 | L1 T | Linf T |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 150 | 20.0 | 8.67e-04 | 6.71e-04 | 8.04e-04 | 7.15e-04 | 3.33e-04 | 1.83e-03 |
| 300 | 10.0 | 2.84e-04 | 2.14e-04 | 2.57e-04 | 2.28e-04 | 8.05e-05 | 3.96e-04 |
| 600 | 5.0 | 9.68e-05 | 6.82e-05 | 8.14e-05 | 7.20e-05 | 1.95e-05 | 7.53e-05 |
| 1200 | 2.5 | 3.85e-05 | 2.35e-05 | 2.79e-05 | 2.44e-05 | 1.22e-05 | 6.97e-05 |
| order 150 to 300 | | 1.61 | 1.65 | 1.65 | 1.65 | 2.05 | 2.20 |
| order 300 to 600 | | 1.56 | 1.65 | 1.66 | 1.66 | 2.04 | 2.40 |
| order 600 to 1200 | | 1.33 | 1.54 | 1.54 | 1.56 | 0.68 | 0.11 |

**Stability of the transport-limited step** (growth after 400 steps):

| mixture | model | cfl_target 0.5 | 0.6 | 0.65 | 0.7 |
| :--- | :--- | ---: | ---: | ---: | ---: |
| N2, 300 K | 1 mixture-averaged | 3.5e-5 | 1.1e-4 | 7.4e3 | 1.7e4 |
| N2, 300 K | 2 unity Lewis | 3.5e-5 | 1.1e-4 | 9.9e3 | 2.8e4 |
| H2/O2/H2O/N2, 1050 K | 1 | 1.5e-4 | 8.9e-5 | 1.4e-4 | 1.1e-4 |
| H2/O2/H2O/N2, 1050 K | 2 | 7.8e-5 | 3.5e-5 | 1.1e4 | 7.1e4 |

An A100 (NVHPC 25.11, OpenACC) gives the same growth to four digits in the seven cases run there (N2 under unity Lewis at 0.5, 0.6 and 0.65, N2 mixture-averaged at 0.65, the mixture under unity Lewis at 0.6 and 0.65 and mixture-averaged at 0.7).

The mixture was run at 1050 K: at 1000 K, the switch temperature of the NASA-7 fits, a grid-scale temperature sawtooth is not damped and lowers the threshold.

## Assessment

- Trace-species diffusion and conduction converge to the erfc solutions at second order (2.0 to 2.2 between 32 and 512 cells), with an L1 of 9.6e-4 of the step at 12.8 cells per front width.
- The multicomponent mixing layer agrees with the independent solver to 3.9e-5 (H2) and 1.2e-5 (T) of the jump in L1 at 2.5 um cells (1200 cells). Since the face composition of the correction flux became second order ([#2205](https://github.com/MFlowCode/MFC/pull/2205)) the species errors fall at 1.5 to 1.7 between 150 and 600 cells and at 1.3 to 1.6 from 600 to 1200 (they fell at about 1.0 beyond 300 cells before it); the temperature error falls at 2.0 to 600 cells and then stalls at 1.2e-5, 0.7 from 600 to 1200. Neither is the second order of the trace case.
- `chem_params%%transport_dt` is stable at its default `cfl_target` 0.5 for both transport models and goes unstable between 0.6 and 0.65, the SSP-RK3 limit 0.628, whenever conduction sets the step: N2 under both models (a single species has no diffusion term, so the two models are the same) and the mixture under unity Lewis. For the mixture-averaged mixture, \f$D_{H_2}\f$ (7.6e-4 m^2/s) is twice \f$\lambda/(\rho c_v)\f$ (3.7e-4), so the step is half the conduction limit and the temperature-spike test stays stable to 0.7; it perturbs T only, so it does not excite the species mode that \f$D_{H_2}\f$ limits and does not test it.
- Partially validated: the transport limit was confirmed on the CPU and on an NVIDIA A100 (NVHPC); the order study of the mixing layer is on one grid family and one reference solver.

## Known limitations

- The step initial condition makes the first steps of the mixing layer an under-resolved start; the errors above include it and converge as the grid is refined.
- The extrapolation outflow (`bc_x%%end = -3`) lets the pressure of the column drift by about 3.4e-4 in the run, a boundary-condition effect that shifts the profiles in x. The comparison is made in the mass coordinate with an isentropic temperature correction for it. A characteristic outflow holding the pressure (`bc_x%%end = -8` with `bc_x%%grcbc_out`) was tried and drifts further (up to 6e-3 at 300 cells), so it is not used.
- A species with \f$Y = 0\f$ does not limit the step (it has no flux), so the limit jumps when a trace of it appears: in air a trace of H sets \f$D_\mathrm{max}\f$ to \f$D_H\f$ and \f$\Delta t\f$ shrinks by the ratio (about 4 for air at 300 K). Reacting runs, in which every radical is present, do not gain from the rule.
- Viscous transport (\f$\mu/\rho\f$ in the limit) was not tested here; the trace and small-\f$\Delta T\f$ cases use amplitudes of 1e-6 and 0.01 K so the variable-property terms are below the discretization error, and the variable-coefficient part is covered by the mixing layer only.
- The mixing layer does not reach second order for the species, and its temperature error stalls beyond 600 cells (L1 1.2e-5, Linf 7e-5). Three causes are candidates and were not separated: the first-order start of the step, the pressure drift of the outflow (the isentropic temperature correction removes the uniform part), and the species-enthalpy flux, which uses the same limited face composition. Earlier, with the upwind face composition, the species stalled at first order and the stall of T was attributed to it; here both are better but the stall of T remains.

## Reproduce

```
./mfc.sh run examples/1D_diffusion_verification/case.py -t pre_process simulation post_process -n 1 -- --mode erf_species --n 256
python3 examples/1D_diffusion_verification/check_diffusion.py examples/1D_diffusion_verification --case-args '--mode erf_species --n 256'
python3 examples/1D_diffusion_verification/ref_mixing_layer.py 3000 ref3000.npz
./mfc.sh run examples/1D_diffusion_verification/case.py -t pre_process simulation post_process -n 1 -- --mode mixing --n 300
python3 examples/1D_diffusion_verification/check_mixing_layer.py examples/1D_diffusion_verification ref3000.npz --case-args '--mode mixing --n 300'
./mfc.sh run examples/1D_diffusion_verification/case.py -t pre_process simulation post_process -n 1 -- --mode stab --tm 2 --cfl 0.65
```

`--mode erf_heat`, `--tm`, `--mixture mix` and `--cfl` select the other cases. The reference run of 3000 cells takes about an hour on one core. The suite runs the CI-sized versions
(`1D -> Chemistry -> Diffusion Verification -> Species Erf`, `Heat Erf` and `Unity-Lewis Transport-Limited dt`).
