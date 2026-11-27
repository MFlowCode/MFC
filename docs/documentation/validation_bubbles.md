@page validation_bubbles Validation: Bubble Dynamics

# Bubble dynamics: Keller-Miksis collapse in a uniform bubbly liquid

**Status:** validated (Euler-Euler bubbles, one bubble size, polytropic gas). **Code:** `master` at 5838b641c. **Example:** `examples/0D_bubblecollapse_adap`. **Last run:** 2026-11-24. See also @ref validation.

## Purpose

The Euler-Euler bubble model (`bubbles_euler`) carries a bubble radius and wall velocity per cell and advances them with a bubble
ODE, here Keller-Miksis \cite Keller80, sub-stepped adaptively (`adap_dt`). In a uniform periodic mixture every cell holds the same
bubble, so the radius history can be compared with the ODE integrated on its own. This checks the ODE, the adaptive sub-steps and
the coupling to the mixture pressure, and measures the collapse time against Rayleigh's \cite Rayleigh17.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | 1D periodic, 31 cells (uniform state: no fluxes) |
| Liquid | water, stiffened gas (gamma 7.1, pi_inf 306 MPa), 1000 kg/m^3, mu 1.002e-3 Pa s, sigma 0.07275 N/m, p_v 2338.8 Pa |
| Bubbles | R0 = 50 um, at rest, equilibrium at p0 = 8236 Pa, polytropic air (gamma 1.4), one size (`nb = 1`), void fraction 1e-5 |
| Drive | mixture pressure 1000 p0 at t = 0 |
| Numerics | 5-equation model, WENO5-M, HLLC, SSP-RK3, fixed dt = 1e-4 t0 (t0 = R0/sqrt(p0/rho_l) = 17.4 us), `adap_dt_tol` 1e-4; 501 steps, a save every step |
| Reference | Keller-Miksis in the units of the case (the same equation as `m_bubbles.fpp`: polytropic gas with surface tension and viscosity, liquid sound speed of the initial state), SciPy LSODA, rtol 1e-11; driven by (a) the constant initial pressure, (b) MFC's cell pressure interpolated between saves |
| Measured | R = nR001/n from the saves, by `check_collapse.py` |

## Results

Rayleigh time \f$t_c = 0.915 R_0\sqrt{\rho_l/(p_\infty - p_v)}\f$ = 0.02894 t0 (0.504 us). GNU 13.3 (CPU, one rank):

| Run | First minimum, MFC [t_c] | vs KM (a) | max abs(R - R_KM)/R0 (a) | same, (b) |
| :--- | ---: | ---: | ---: | ---: |
| default | 1.0425 | +7.2e-4 | 4.0e-3 | 8.9e-5 |
| `--adap-tol 1e-6` | 1.0425 | +7.2e-4 | 4.0e-3 | not run |
| `--dt 2.5e-5` | 1.0414 | +1.5e-4 | 4.2e-3 | 8.0e-5 |

An NVIDIA A100 (NVHPC 25.11, OpenACC, not case-optimized) gives the same values to the digits shown.

## Assessment

- Driven by its own cell pressure, the ODE reproduces MFC's radius to 9e-5 R0 through the collapse and rebound (1.7 t_c), and the
  result does not move with the sub-step tolerance or the time step. The bubble ODE and its adaptive integration are verified.
- Against the constant-pressure ODE the radius differs by 4e-3 R0. This is the two-way coupling, not an error: the mixture
  pressure falls by 0.27 % as the void fraction collapses from 1e-5 to 3e-8 (stiff liquid, fixed mixture density), and the bubble
  sees the lower pressure.
- The collapse takes 1.041 to 1.043 Rayleigh times: the gas, surface tension, viscosity and liquid compressibility of
  Keller-Miksis lengthen it by 4 % over the empty-cavity value. The time of the minimum is a parabola through the three saves
  around it and moves by 0.1 % with the save interval; the sampled minimum radius (0.04 to 0.05 R0) is not resolved by the saves.

## Known limitations

- One bubble size, polytropic gas, no wave propagation. Polydisperse and QBMM closures (`examples/1D_qbmm`), the non-polytropic
  thermal model, bubble screens with acoustic coupling and the Lagrangian sub-grid bubbles (`bubbles_lagrange`) are not checked here.
  The Lagrangian model requires `thermal = 3` with the non-polytropic gas model, so the polytropic ODE above does not apply to it.
- Only the Keller-Miksis model is run; Rayleigh-Plesset and Gilmore are not.

## Reproduce

```
./mfc.sh run examples/0D_bubblecollapse_adap/case.py -t pre_process simulation post_process -n 1
python3 examples/0D_bubblecollapse_adap/check_collapse.py examples/0D_bubblecollapse_adap
```

Options of `case.py`: `--dt`, `--adap-tol`. The suite runs the example at a reduced size (every example not on the skip list in `toolchain/mfc/test/cases.py` is registered with at most 25 cells and 50 steps), which guards the code path; the checker validates the physics. No test was added.
