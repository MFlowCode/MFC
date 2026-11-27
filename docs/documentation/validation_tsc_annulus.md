@page validation_tsc_annulus Isothermal Conduction Annulus: Image-Point Interpolation

# Isothermal Conduction Annulus: Image-Point Interpolation

**Status:** partially validated (wall temperature and whole-gap error validated at three resolutions; wall gradient measured only as a fit estimate).
**Code:** `master` at d1e6852c1, with `ib_ip_interp` ([#2165](https://github.com/MFlowCode/MFC/pull/2165)). TSC was opt-in when this page was written; it is now the default (`ib_ip_interp = 1` gives the inverse-distance weights).
**Example:** `examples/2D_ibm_conduction_annulus`. **Last run:** 2026-11-21.

## Purpose

The ghost-cell immersed boundary method imposes an isothermal wall temperature through the state interpolated at an image point inside the gas.
The weights of that interpolation set the accuracy of the imposed wall temperature. This page compares the inverse-distance weights (`ib_ip_interp = 1`, the default when this page was written)
with the triangular-shaped-cloud weights (`ib_ip_interp = 2`; the kernel is described by Ezra, Peles and Kozak, Int. J. Heat Mass Transf. 252:127314, 2025)
against the exact conduction solution of an annulus.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | inner circle \f$R = 1\f$ (`geometry = 2`), outer flat ring (`geometry = 5`, `ring.stl`, 1024 segments) with inner edge \f$r_2 = 2R\f$; both static |
| Gas | quiescent ideal gas, \f$\gamma = 1.4\f$, \f$c_v = 1\f$, \f$\rho = 1\f$ at the mean wall temperature, constant `k_therm` = 0.1 (diffusivity 0.1); inviscid |
| Walls | `thermal_bc = 1`, no slip; \f$T_1 = 2\f$ (inner), \f$T_2 = 3\f$ (outer) |
| Grid | 10, 20, 40 cells per \f$R\f$, square domain just larger than the ring's inner edge, extrapolation boundaries |
| Numerics | WENO5, HLLC, SSP-RK3, `fd_order = 2`, fixed \f$\Delta t = \min(0.25\Delta x/c, 0.15\Delta x^2/\alpha)\f$, 8 slowest conduction times \f$\tau = ((r_2 - R)/\pi)^2/\alpha\f$ |
| Initial state | the exact profile on 10 concentric constant shells at uniform pressure |
| Settings | `ib_ip_interp` 1 or 2; every other `ib_*` option at its default |

## Reference

Exact steady conduction at constant conductivity: \f$T(r) = T_1 + (T_2 - T_1)\ln(r/R)/\ln(r_2/R)\f$. It holds for a gas at rest at any uniform pressure.

## Quantities

Errors are (MFC - exact) divided by \f$T_2 - T_1\f$ (`check_annulus.py`).

- **Wall error:** a quadratic least-squares fit of \f$T(s)\f$ through the gas cells \f$0.5\Delta x < s < 3.5\Delta x\f$ (s the distance into the gas), extrapolated to \f$s = 0\f$, minus the wall value.
- **Gradient error:** \f$dT/dr\f$ of the same fit at the wall, relative to the exact gradient.
- **Gap error:** relative L2 and maximum error of \f$T\f$ over the gas cells more than \f$0.5\Delta x\f$ from either wall.

## Results

Measured on AMD MI300A (amdflang), one rank, the code version above.

| Cells per R | Inner wall, idw | Inner wall, TSC | Outer wall, idw | Outer wall, TSC |
| ---: | ---: | ---: | ---: | ---: |
| 10 | −1.83e-2 | +3.68e-3 | +1.18e-2 | +1.31e-4 |
| 20 | −1.18e-2 | +7.67e-4 | +7.71e-3 | +6.05e-5 |
| 40 | −5.00e-3 | +1.87e-4 | +3.42e-3 | +2.96e-5 |
| Order (10→20, 20→40) | 0.6, 1.2 | 2.3, 2.0 | 0.6, 1.2 | 1.1, 1.0 |

| Cells per R | L2, idw | L2, TSC | Max, idw | Max, TSC | Inner gradient, idw | Inner gradient, TSC |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 9.36e-3 | 1.45e-3 | 4.57e-2 | 2.66e-3 | −5.0e-3 | −2.3e-2 |
| 20 | 5.32e-3 | 3.90e-4 | 1.98e-2 | 8.68e-4 | +2.7e-3 | −5.8e-3 |
| 40 | 2.02e-3 | 1.86e-4 | 8.94e-3 | 3.39e-4 | +1.4e-2 | −4.3e-4 |

Comparison with an earlier measurement of the same geometry (a start from a uniform temperature, a build before the current master; not repeated here): the wall errors agree within a factor of 1.5
(inner idw 2.6e-2, 1.7e-2, 5.0e-3; inner TSC 5.8e-3, 1.0e-3, 3.0e-4; outer idw 1.1e-2, 7.6e-3, 3.5e-3; outer TSC 3.1e-4, 4.5e-5, 2.0e-5). The gap errors do not: that measurement
gave a whole-gap L2 of 5.1e-2, 2.5e-2, 1.3e-2 for both weightings. The present table starts from the exact profile; a uniform start on the current master gives a TSC gap L2 of 2.6e-3, 5.5e-3, 6.0e-3
and an inverse-distance one of 7.8e-3, 3.3e-3, 4.5e-3.

## Assessment

- TSC lowers the imposed wall temperature error by a factor of 5 to 130 at every resolution; on the inner wall it converges at about second order, on the outer wall at first order from a much lower level. Inverse distance is first order or worse.
- The whole-gap L2 error is 6 to 14 times lower with TSC and converges at order 1.9 then 1.1; inverse distance gives 0.8 then 1.4.
- The wall gradient from a fit through the first cells converges for TSC on the inner wall and not for inverse distance. This is an estimator, not the heat flux the solver applies to the wall.

## Known limitations

- Constant conductivity and a quiescent gas only. A temperature-dependent conductivity and the wall heat flux need chemistry (`ib_heat_wrt`) and are not covered.
- A start from a uniform temperature leaves acoustic ringing in the closed gap (Mach number near 1e-2, no damping in an inviscid run) that sets a floor of a few 1e-3 on the gap error. The example starts from the exact profile for that reason.
- Viscous runs (`viscous = T`, `Re(1) = 10`) completed at 12 cells/R and returned NaN in the last step at 20 cells/R for both weightings (cause not investigated).
- Rotating (Taylor-Couette) walls, the minimum image-point distance (`ib_ip_min_dist`) and the quadratic wall fill are not covered.
- Isothermal ghost conduction through the unbounded mirror temperature is under review in [#2248](https://github.com/MFlowCode/MFC/pull/2248); the numbers above are without it.

## Reproduce

```
./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 -- --cells_per_R 20 --interp 2
python3 examples/2D_ibm_conduction_annulus/check_annulus.py examples/2D_ibm_conduction_annulus
```

The suite runs `2D -> IBM -> Conduction Annulus -> Inverse Distance` and `-> TSC` at 12 cells/R and 300 steps; the goldens are regression guards and the checker validates the physics.
