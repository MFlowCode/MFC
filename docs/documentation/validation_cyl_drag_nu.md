@page validation_cyl_drag_nu Validation: Circular Cylinder Drag and Nusselt Number

# Circular cylinder, Re 20-100: drag, wake, shedding and heat transfer

**Status:** partially validated. **Code:** `master` at 615348a25. **Example:** `examples/2D_ibm_cylinder_drag_nu`.
**Last run:** 2026-11-21.

## Purpose

A ghost-cell immersed boundary with chemistry-based thermodynamics and transport is the code path of particle-resolved reacting
runs. This checks the drag from the face-flux force of `ib_force_wrt` (it replaced a volume integral that missed part of the
friction drag, [#2013](https://github.com/MFlowCode/MFC/issues/2013)), the wake and shedding frequency, and the Nusselt number
from the face heat of `ib_heat_wrt`, against references and against a force that does not use the immersed-boundary code.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | 2D cylinder (`patch_ib` geometry 2), centred off the symmetry planes by (dx/3, dx/5) to seed shedding |
| Flow | M = 0.2, Re_D = 20, 40, 100, 300 K, 1 atm, dry air (`h2o2.yaml`: O2 0.21, N2 0.79 by mole), chemistry on for thermodynamics and mixture-averaged transport, reactions off; constant \f$\mu\f$ from the free stream |
| Grid | 10 and 20 cells/D, uniform; domain [-15 D, 25 D] x [-16 D, 16 D] (blockage 3.1 %); Re 100: [-8 D, 16 D] x [-8 D, 8 D] (6.25 %) |
| Boundaries | subsonic characteristic inflow, characteristic outflow holding the free-stream pressure (extrapolation at the outflow lets the far-field pressure creep up), extrapolated lateral faces |
| Wall | fixed, no slip, `thermal_bc = 1`, Twall = T_inf (drag runs) or 1.1 T_inf (heat runs) |
| Numerics | WENO5-M, HLLC, SSP-RK3, `weno_Re_flux`, `fd_order = 4`, adaptive step with CFL min(0.5, 0.4 Re/N) and `chem_params%%transport_dt`; `ib_ip_interp = 1`; these runs used `ib_wall_order = 1` and `ib_visc_compact = F`, the defaults then (the Reproduce commands pass `--wall_order 1 --visc_compact 0`; the example now follows the solver defaults 2 and T) |
| Length | \f$t U/D\f$ = 50 (Re 20, 40), 120 (Re 100); averages over the last 20 % |
| Hardware | AMD MI300A, OpenMP offload, double precision |

## Reference data

- Dennis and Chang \cite DennisChang1970: Cd 2.05 (Re 20) and 1.52 (Re 40); L/D 0.94 and 2.35. Fornberg \cite Fornberg1980: Cd 2.00 and 1.50 (values as compiled in the literature; the original tables were not re-read).
- Henderson \cite Henderson1995 (Cd 1.35) and Park et al. \cite Park1998 (Cd 1.33, Cl amplitude 0.332, St 0.165) at Re 100.
- Williamson \cite Williamson1989: St = -3.3265/Re + 0.1816 + 1.6e-4 Re (0.1643 at Re 100).
- Nu: Lange et al. \cite Lange1998 (2D DNS fit, Pr = 0.7) and Churchill and Bernstein \cite ChurchillBernstein1977 (evaluated at the film Pr = 0.712).

## Quantities

- **Cd_IB**: \f$\langle F_x\rangle/(\tfrac12\rho U^2 D)\f$ from `D/ib_forces.dat`.
- **Cd_CV**: the same from a steady momentum balance on the square of half-width 1, 2 or 4 D, from the saved fields (`check_cylinder.py`); it drops the unsteady and acoustic terms.
- **L/D**: first u = 0 crossing on the wake centreline behind the rear stagnation point, last save.
- **St** from the dominant frequency of the lift history (second half of the run).
- **Nu** = \f$-(Q_\mathrm{hot}-Q_\mathrm{cold})\,D/(\pi D\,k_f\,(T_w-T_\infty))\f$ with Q the heat the body takes (`D/ib_heat.dat`, face fluxes), \f$k_f\f$ at the film temperature, and the cold (Tw = T_inf) twin removing the viscous heating of the wall.

## Results

| Re | cells/D | Cd_IB | err | Cd_CV 1 D | 2 D | 4 D | L/D | ref L/D |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 10 | 1.982 | -3.3 % | 2.023 | 1.983 | 1.975 | 0.721 | 0.94 |
| 20 | 20 | 2.029 | -1.0 % | 2.051 | 2.029 | 2.022 | 0.818 | 0.94 |
| 40 | 10 | 1.500 | -1.3 % | 1.525 | 1.498 | 1.490 | 1.932 | 2.35 |
| 40 | 20 | 1.521 | +0.1 % | 1.534 | 1.519 | 1.511 | 2.064 | 2.35 |

Errors of Cd_IB are against Dennis and Chang (Re 20: 2.05, Re 40: 1.52).

Re 100, 20 cells/D, \f$t U/D\f$ = 120 (smaller domain, Cl amplitude saturated by about 90): mean Cd_IB 1.354 (+0.3 % against 1.35,
+1.8 % against 1.33); Cd_CV 1.344, 1.332, 1.302 for 1, 2, 4 D; Cl amplitude 0.287 (references 0.332 and 0.339); St 0.1668
(+1.5 % against Williamson).

Heat transfer, Tw/T_inf = 1.1, 20 cells/D:

| Re | Nu | Lange | err | Churchill-Bernstein | err |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 2.442 | 2.409 | +1.4 % | 2.479 | -1.5 % |
| 40 | 3.243 | 3.281 | -1.2 % | 3.386 | -4.2 % |
| 100 | 4.962 | 5.128 | -3.2 % | 5.190 | -4.4 % |

## Assessment

- The face-flux force agrees with the momentum balance (box 2 D) to 0.1 % at Re 20 and 0.1-0.3 % at Re 40 at 20 cells/D, and with the references to within 1.0 % at Re 20 and 40 and 0.3-1.8 % at Re 100. It moves toward the reference from 10 to 20 cells/D (Cd changes by 2.4 % at Re 20 and 1.4 % at Re 40).
- The Strouhal number is +1.5 % against Williamson at 20 cells/D.
- Nu is within 1.4 % of Lange at Re 20 and 40 and 3.2 % low at Re 100.
- Status is partially validated: 40 cells per diameter, the grid-convergence order and Re 100 on the larger domain were not measured here.

## Known limitations

- **L/D is first order and low** (0.818 against 0.94 at Re 20, 2.064 against 2.35 at Re 40: 13 % and 12 % low at 20 cells/D, 23 % and 18 % low at 10). These runs used the wall-value ghost (`ib_wall_order = 1`), not a reflection, so the effective wall sits up to a cell inside the surface ([#2014](https://github.com/MFlowCode/MFC/issues/2014)). The quadratic fill, `ib_wall_order = 2`, is now the default; these tables were not repeated with it.
- Not re-measured (cost, about 5 GPU-hours per run on the 40 D x 32 D domain): 40 cells per diameter, and Re 100 on the domain of Re 20 and 40.
- Re 100 uses a 6.25 % blockage; Cd_CV depends on the box there (1.302 for 4 D, which comes within 4 D of the inflow and lateral faces).
- Nu varies by 2-4 % between correlations. The face-heat Nu carries no wall-gradient fit: an earlier drift of Nu with refinement at Re 40 came from fitting the wall gradient of the temperature field, not from the solver.
- The film Prandtl number of dry air in `h2o2.yaml` is 0.712; Lange's fit is for 0.7.

## Reproduce

```
./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
    --silo --Re 40 --cells_per_D 20 --xu 15 --xd 25 --hy 16 --steps 0 --tstar 50 --wall_order 1 --visc_compact 0
./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
    --silo --Re 40 --cells_per_D 20 --xu 15 --xd 25 --hy 16 --steps 0 --tstar 50 --tw_ratio 1.1 --wall_order 1 --visc_compact 0
python3 examples/2D_ibm_cylinder_drag_nu/check_cylinder.py <cold run> --hot <hot run>
```

Re 100: `--Re 100 --tstar 120 --xu 8 --xd 16 --hy 8`. A 20 cells/D run at Re 40 takes about 40 minutes on one MI300A
(four runs ran concurrently on four GPUs).
The test suite runs the CI-sized case (`--tw_ratio 1.1`, 16 cells/D, 200 steps) as `2D -> Chemistry -> IBM -> Cylinder Drag Nu`.
