@page validation_sphere_drag_nu Validation: Sphere Drag and Nusselt Number

# Sphere at Mach 0.2, Re 20 and 100: drag, wake and heat transfer

**Status:** partially validated. **Code:** `master` at 0cdfdf02a. **Example:** `examples/3D_ibm_sphere_drag_nu`.
**Last run:** 2026-11-23.

## Purpose

The 3D ghost-cell immersed boundary with chemistry-based thermodynamics and transport is the code path of particle-resolved reacting
runs. This checks the drag from the face-flux force of `ib_force_wrt` and the Nusselt number from the face heat of `ib_heat_wrt` of a
sphere against the standard correlations, and the drag against a momentum balance that does not use the immersed-boundary code. It is the
3D companion of the cylinder check (@ref validation_cyl_drag_nu).

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | sphere (`patch_ib` geometry 8) off the symmetry planes by (dx/3, dx/5, dx/7) |
| Flow | M = 0.2, Re_D = 20 and 100, 300 K, 1 atm, dry air (`h2o2.yaml`: O2 0.21, N2 0.79 by mole), chemistry on for thermodynamics and mixture-averaged transport, reactions off; constant \f$\mu\f$ from the free stream |
| Grid | 8 and 12 cells/D, uniform (12 cells/D: 108 x 72 x 72 cells); domain [-3 D, 6 D] x [-3 D, 3 D]^2 (area blockage 2.2 %) |
| Boundaries | subsonic characteristic inflow, characteristic outflow holding the free-stream pressure, extrapolated lateral faces |
| Wall | fixed, no slip, `thermal_bc = 1`, Twall = T_inf (drag runs) or 1.1 T_inf (heat runs); `ib_ip_interp = 1` |
| Numerics | WENO5-M, HLLC, SSP-RK3, `weno_Re_flux`, `fd_order = 4`, adaptive step with CFL 0.5 and `chem_params%%transport_dt` |
| Length | \f$t U/D\f$ = 15 (Re 20) and 20 (Re 100) at 12 cells/D, 12 at 8 cells/D; the force and the heat are constant to 1e-4 over the last saves at 12 cells/D (the 8 cells/D runs stop at 12 and 12.3, where Cd of Re 100 still falls by 0.2 % per two units) |
| Hardware | NVIDIA A100 80GB (Georgia Tech), NVHPC 25.11, OpenACC, double precision |

## Reference data

- Schiller and Naumann \cite SchillerNaumann1933: \f$C_d = 24/Re\,(1 + 0.15\,Re^{0.687})\f$ (quoted accuracy 5 % for Re < 800); incompressible, so the Prandtl-Glauert effect of about 2 % at M = 0.2 is not removed.
- Ranz and Marshall \cite RanzMarshall1952: \f$Nu = 2 + 0.6\,Re^{1/2}Pr_f^{1/3}\f$; Whitaker \cite Whitaker1972: \f$Nu = 2 + (0.4\,Re^{1/2} + 0.06\,Re^{2/3})\,Pr_f^{0.4}\f$ (viscosity ratio 1); at the free-stream Re and the film Prandtl number \f$Pr_f\f$ = 0.712.
- Johnson and Patel \cite JohnsonPatel1999: no recirculation below Re of about 24; recirculation length L/D = 0.88 at Re 100 (value as quoted, the paper's table was not re-read).

## Quantities

- **Cd_IB**: \f$\langle F_x\rangle/(\tfrac18\rho U^2\pi D^2)\f$ from `D/ib_forces.dat`, mean over the last 20 % of the run.
- **Cd_CV**: the same drag from a steady momentum balance on cubes of half-width 1 and 2 D (the fields of the last saves; unsteady and acoustic terms dropped).
- **L/D**: first u = 0 crossing on the wake axis behind the rear stagnation point, last save.
- **Nu** = \f$-(Q_\mathrm{hot}-Q_\mathrm{cold})\,D/(\pi D^2 k_f (T_w - T_\infty))\f$ with Q the heat the body takes (`D/ib_heat.dat`, face fluxes), \f$k_f\f$ at the film temperature, and the cold (Tw = T_inf) twin removing the viscous heating of the wall.
- Script: `examples/3D_ibm_sphere_drag_nu/check_sphere.py`.

## Results

| Re | cells/D | Cd_IB | error vs SN | Cd_CV 1 D | 2 D | L/D | Nu | error vs RM | error vs W |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 8 | 2.532 | -3.0 % | 2.646 | 2.534 | 0 | 4.083 | -7.1 % | +3.4 % |
| 20 | 12 | 2.617 | +0.3 % | 2.694 | 2.618 | 0 | 4.059 | -7.7 % | +2.8 % |
| 100 | 8 | 1.114 | +2.1 % | 1.121 | 1.102 | 0.755 | 6.939 | -5.7 % | +4.8 % |
| 100 | 12 | 1.101 | +0.9 % | 1.105 | 1.100 | 0.795 | 6.905 | -6.2 % | +4.3 % |

Schiller-Naumann: 2.610 (Re 20) and 1.092 (Re 100); Ranz-Marshall 4.396 and 7.359, Whitaker 3.948 and 6.621.

## Assessment

- The face-flux drag is within 1 % of Schiller-Naumann at 12 cells/D (+0.3 % at Re 20, +0.9 % at Re 100) and agrees with the momentum balance on the 2 D cube to 0.03 % at Re 20 and 0.15 % at Re 100; the 1 D cube is 2.9 % and 0.5 % high.
- The Nusselt number lies between the Ranz-Marshall and Whitaker correlations, which differ by 11 % from each other: 6-8 % below Ranz-Marshall and 3-4 % above Whitaker.
- The Re 20 sphere has no recirculation, as expected, and the Re 100 wake length (0.795 D) is 10 % below 0.88.
- Status is partially validated: Re 300 and the finer grid (20 cells/D) were not re-run; the 8-to-12 cells/D change (Cd +3.4 % at Re 20, -1.2 % at Re 100; Nu -0.6 % and -0.5 %) is the only convergence evidence, and at Re 100 the 8 cells/D run is not fully settled.

## Known limitations

- **First-order wall.** The no-slip ghost value is the wall value, so the effective wall sits up to a cell inside the surface ([#2389](https://github.com/MFlowCode/MFC/issues/2389)); the wake length is low and the Nusselt number may drift with refinement (not verified).
- **Domain.** The 2.2 % blockage of the 6 D x 6 D cross-section raises Cd; not measured here.
- Not re-run on master: Re 300 (the sphere sheds late and needs t U/D of about 100), the Mach 0.3-2 rows, and 20 and 40 cells/D. Earlier measurements of this configuration at 20 cells/D on a pre-merge build gave Cd within -1.6 to +3.9 % of Schiller-Naumann and Nu between Ranz-Marshall and Whitaker, the same bands; they are not repeated here.
- Schiller-Naumann is itself a fit of about 5 %, and the references are incompressible; Whitaker's and Ranz-Marshall's fits differ by 11 %, so Nu is a consistency check, not a validation to better than that.
- The CPU and GPU builds of this case differ by 3e-4 in a CI-sized run (see the test below), which is the size of the cross-architecture noise for a body in no-slip flow; the table is from one GPU.

## Reproduce

```
./mfc.sh run examples/3D_ibm_sphere_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
    --silo --Re 100 --cells_per_D 12 --xu 3 --xd 6 --hy 3 --steps 0 --tstar 20 --frames 5
./mfc.sh run examples/3D_ibm_sphere_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
    --silo --Re 100 --cells_per_D 12 --xu 3 --xd 6 --hy 3 --steps 0 --tstar 20 --frames 5 --tw_ratio 1.1
python3 examples/3D_ibm_sphere_drag_nu/check_sphere.py <cold run> --hot <hot run>
```

Re 20: `--Re 20 --tstar 15 --frames 5`. The runs were stopped after their last save of interest, so give `--tstar` and `--frames` of the saves you keep. Each run writes next to `case.py`, so copy the directory for each. A 12 cells/D run takes about 2-3 hours on one A100 shared with other jobs.
The test suite runs the CI-sized case (`--tw_ratio 1.1`, 6 cells/D, 30 steps) as `3D -> Chemistry -> IBM -> Sphere Drag Nu`.
