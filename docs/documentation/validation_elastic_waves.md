@page validation_elastic_waves Validation: Hypoelastic Wave Speeds

# Hypoelastic solids: longitudinal and shear Riemann problems

**Status:** validated (first order in the jump, at small jumps). **Code:** `master` at 5838b641c. **Examples:** `examples/1D_sodHypo`, `examples/2D_hypo_shear_contact`. See also @ref validation.

With `hypoelasticity = T` the deviatoric stress obeys \f$\dot\tau = 2G\,\mathrm{dev}(D)\f$ plus the objective rotation terms, so a 1D compression changes \f$\tau_{xx}\f$ at the rate \f$\frac{4}{3}G\,\partial_x u\f$ and a shear changes \f$\tau_{xy}\f$ at \f$G\,\partial_x v\f$.
Linear elastodynamics then gives longitudinal waves at \f$c_L = \sqrt{(K + 4G/3)/\rho}\f$, with \f$K = \gamma(p + \pi_\infty)\f$ the bulk modulus of the stiffened gas, and shear waves at \f$c_S = \sqrt{G/\rho}\f$.
For a jump at rest, the Riemann solution of the linear system is

- normal stress \f$\sigma = -p + \tau_{xx}\f$: \f$u^* = (\sigma_R - \sigma_L)/(Z_l + Z_r)\f$ with \f$Z = \rho c_L\f$; across each longitudinal wave \f$dp = -(K/M)\,d\sigma\f$ and \f$d\tau_{xx} = (4G/3)/M\,d\sigma\f$, \f$M = K + 4G/3\f$;
- shear stress: \f$v^* - v_0 = (\tau_{xy,R} - \tau_{xy,L})/(Z_l + Z_r)\f$ with \f$Z = \rho c_S\f$, and \f$\tau_{xy}^* = \tau_{xy,L} + Z_l (v^* - v_0)\f$.

`examples/1D_sodHypo/check_elastic.py` measures the plateau values (median over the middle 60 % of each plateau) and the front positions (where u, or \f$v - v_0\f$, crosses half its plateau value) and prints their errors against these formulas.
All tables were measured with GNU 13.3 (CPU, one rank, double precision) on `master` at 5838b641c; one resolution of each was repeated on an NVIDIA A100 (NVHPC 25.11, OpenACC) with the same digits shown, except where noted.

## Longitudinal waves {#sec-validation-elastic-longitudinal}

| Item | Value |
| :--- | :--- |
| Material | water, stiffened gas \f$\gamma\f$ = 4.4, \f$\pi_\infty\f$ = 6e8 Pa, 1000 kg/m³, G = 10 GPa: \f$c_L\f$ = 3997 m/s against a sound speed of 1625 m/s |
| Riemann problem | \f$x \in [0, 1]\f$ m, jump at 0.5 m: left pressure `--pL` (1e8 Pa in the example, 1e6 Pa for the small-jump sweep), right 1e5 Pa, \f$\tau = 0\f$, at rest |
| Numerics | WENO5, HLL, SSP-RK3, `fd_order` = 4, 200 steps to 64 µs (202 cells; dt scales with the cells); extrapolation ends |
| Measured | \f$u^*\f$, \f$p^*\f$ and \f$\tau_{xx}^*\f$ on both sides of the contact (errors relative to the jump: \f$|u^*|\f$ for \f$u^*\f$, \f$|\Delta\sigma|\f$ on that side otherwise), front speeds |

Small jump (`--pL 1e6`, \f$\Delta\sigma/M\f$ = 6e-5):

| cells | \f$u^*\f$ error | worst \f$p^*\f$, \f$\tau_{xx}^*\f$ error | left front speed − 1 | right front speed − 1 |
| ---: | ---: | ---: | ---: | ---: |
| 202 | +6.6e-5 | +1.2e-4 | +5.9e-4 | +6.2e-4 |
| 404 | +1.0e-5 | +3.9e-5 | +1.7e-4 | +2.0e-4 |
| 808 | +7.6e-7 | +2.6e-5 | +8.6e-5 | +1.1e-4 |
| 1616 | +1.2e-7 | +2.4e-5 | +3.5e-5 | +6.1e-5 |

The front speeds converge to \f$c_L\f$ (to 4e-5 at 1616 cells); the plateau stresses level off at 2.4e-5 of the jump, which is the size of the neglected second-order terms (\f$\Delta\sigma/M\f$ = 6e-5). The 404-cell run on the A100 gives the same values to the digits printed.

Example jump (1e8 Pa, \f$\Delta\sigma/M\f$ = 6e-3), 202 cells: \f$u^*\f$ +6.1e-5, worst plateau stress +2.8e-3, front speeds −1.4e-3 and +2.3e-3. The errors are O(\f$\Delta\sigma/M\f$), as expected of a first-order prediction for a jump this size; the speed of the compressive (right-moving) wave is above \f$c_L\f$ and that of the expansion below it.

## Shear waves {#sec-validation-elastic-shear}

| Item | Value |
| :--- | :--- |
| Material | ideal gas \f$\gamma\f$ = 1.4, G = 1, p = 1, density 1 \| 2 at x = 0.5 (two identical components), \f$c_S\f$ = 1 \| 0.7071 |
| Riemann problem | \f$\tau_{xy}\f$ = 0 \| 1e-3, uniform tangential velocity \f$v_0\f$ = 1, periodic in y (8 cells) |
| Numerics | first-order MUSCL, HLLD, forward Euler, CFL 0.5 on \f$c_L\f$, t = 0.25 (`--steps` = ceil(0.8265 N)) |
| Measured | \f$v^* - v_0\f$, \f$\tau_{xy}^*\f$ (relative errors), shear front speeds |

| cells | \f$v^* - v_0\f$ error | \f$\tau_{xy}^*\f$ error | left front speed − 1 | right front speed − 1 |
| ---: | ---: | ---: | ---: | ---: |
| 64 | −3.6e-5 | +4.9e-5 | +1.37e-2 | −1.50e-2 |
| 128 | +8.9e-11 | +1.4e-8 | +7.2e-3 | −7.4e-3 |
| 256 | +1.3e-10 | −3.6e-11 | +3.7e-3 | −3.7e-3 |
| 512 | +7.8e-11 | −3.2e-11 | +1.9e-3 | −1.8e-3 |

The plateau is exact to round-off once the waves have left the contact region (at 64 cells the plateau is only a few cells wide).
The front speeds converge at first order, the order of the scheme: both measured fronts sit about 0.2 cells to the left of the linear positions at every resolution, which makes the left wave look fast and the right one slow by the same amount. The 128-cell run on the A100 gives the same values.

Limits: linear waves only (the jumps are small or their nonlinear error is reported); one material per test; no plastic or damage models; the hyperelastic model (`hyperelasticity`) is not covered.

## Reproduce {#sec-validation-elastic-run}

```
./mfc.sh run examples/1D_sodHypo/case.py -t pre_process simulation post_process -n 1 -- -N 404 --pL 1e6
python3 examples/1D_sodHypo/check_elastic.py examples/1D_sodHypo --case-args "-N 404 --pL 1e6"
./mfc.sh run examples/2D_hypo_shear_contact/case.py -t pre_process simulation post_process -n 1 -- -N 128 --steps 106
python3 examples/1D_sodHypo/check_elastic.py examples/2D_hypo_shear_contact --case-args "-N 128 --steps 106"
```
