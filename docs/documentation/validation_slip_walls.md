@page validation_slip_walls Validation: Slip Immersed Walls

# Slip immersed walls: mass conservation, curvature pressure and the CCST ghost state

**Status:** partially validated (mass budget and exact vortex); open discrepancy for slip walls in unsteady shock flow and close approach.
**Code:** `master` at 0cdfdf02a. **Example:** `examples/2D_ibm_slip_cylinder`.
**Last run:** 2026-11-23.

## Purpose

Gas at tens to hundreds of m/s slides along 50-micron bodies in particle-resolved shock flows, and the slip wall has to conserve the gas
mass and give the right pressure on a curved surface. This checks three things of the slip-wall ghost state of an immersed boundary
against exact results: the gas mass in a closed box (with and without `ib_mass_relax`), the near-wall pressure and velocity of an exact
steady vortex (the curvature term of the ghost pressure, the image-point interpolation `ib_ip_interp` and `ib_slip_ccst`, the
curvature-corrected symmetry technique of \cite DadoneGrossman2004), and records two cases of slip walls that do not work.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | fixed circle R = 50 um (`patch_ib(1)%%geometry = 2`, `patch_ib(1)%%slip = T`, `patch_ib(1)%%moving_ibm = 0`), centred in a square box |
| Gas | inviscid, gamma = 1.4, stagnation state 1e5 Pa and 1.2 kg/m^3 |
| Mass budget | uniform flow at Mach 0.5; box of side 16 R, fully periodic; one flow-through; 16 and 32 cells per R |
| Vortex | \f$u_\theta = M c_0 R/r\f$ with M = 0.3 and 0.6 at the wall, isentropic density and pressure; box of side 6 R, extrapolation boundaries; run time \f$0.4\,(3R)/c_0\f$ before waves from the boundaries return; 10, 20 and 40 cells per R |
| Numerics | WENO5 (mapped, `mp_weno`), HLLC, SSP-RK3, adaptive CFL 0.5, `fd_order = 2` |
| Options set | `ib_mass_relax` 0, 1 or 30 (mass budget); `ib_ip_interp` 1 or 2 and `ib_slip_ccst` F or T (vortex); the curvature term of the ghost pressure is always on |
| Hardware | GNU 13, CPU, one rank, double precision |

## Reference data

- **Mass.** In a periodic box with no source the total gas mass is constant.
- **Exact vortex.** \f$u_\theta = M c_0 R/r\f$, \f$T/T_0 = 1 - \frac{\gamma-1}{2}(u_\theta/c_0)^2\f$, \f$\rho = \rho_0 (T/T_0)^{1/(\gamma-1)}\f$, \f$p = p_0 (T/T_0)^{\gamma/(\gamma-1)}\f$.
  An exact steady Euler solution tangent to the wall, with \f$dp/dr = \rho u_\theta^2/r\f$ (the balance the curvature term of the ghost pressure imposes).

## Quantities

- **Mass change.** The sum of the density over the fluid cells (centres outside the body) at the saves, relative to the first save (`check_slip.py`).
- **Vortex error.** Root-mean-square over the cells of a wall-distance band (in cells) of \f$(p-p_\mathrm{ex})/(\frac12\rho_0 (M c_0)^2)\f$, \f$(u_\theta-u_{\theta,\mathrm{ex}})/(M c_0)\f$ and \f$u_r/(M c_0)\f$ (zero in the exact solution), at the end of the run. Only the first band (0 to 1 cells) is tabulated.

## Results

**Gas mass in a closed box, one flow-through** (changes at half and at the end of the flow-through):

| cells per R | `ib_mass_relax` | half a flow-through | one flow-through |
| ---: | :--- | ---: | ---: |
| 16 | 0 (what -1 gives for a fixed body) | +1.69e-3 | +3.61e-3 |
| 16 | 30 | +4.3e-5 | +5.8e-5 |
| 16 | 1 | -5.2e-14 | -1.2e-13 |
| 32 | 0 | +1.00e-3 | +2.08e-3 |
| 32 | 1 | -1.0e-13 | -2.3e-13 |

**Exact vortex, M = 0.3 at the wall, first-band rms error at 10 / 20 / 40 cells per R** (observed orders 10 to 20, 20 to 40):

| Scheme | Pressure | Order | Tangential velocity | Order | Radial velocity | Order |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| default image point | 3.28e-02 / 1.58e-02 / 6.38e-03 | 1.1, 1.3 | 3.19e-02 / 2.00e-02 / 1.07e-02 | 0.7, 0.9 | 1.31e-02 / 6.73e-03 / 2.61e-03 | 1.0, 1.4 |
| TSC (`ib_ip_interp = 2`) | 1.23e-02 / 4.53e-03 / 1.97e-03 | 1.4, 1.2 | 3.26e-02 / 2.25e-02 / 1.39e-02 | 0.5, 0.7 | 2.92e-03 / 1.29e-03 / 5.37e-04 | 1.2, 1.3 |
| `ib_slip_ccst`, default image point | 5.96e-02 / 3.56e-02 / 1.23e-02 | 0.7, 1.5 | 2.68e-02 / 1.84e-02 / 1.04e-02 | 0.5, 0.8 | 1.94e-02 / 1.23e-02 / 4.39e-03 | 0.7, 1.5 |
| `ib_slip_ccst`, TSC | 4.51e-03 / 1.22e-03 / 2.40e-04 | 1.9, 2.3 | 1.13e-03 / 7.27e-04 / 2.49e-04 | 0.6, 1.5 | 1.64e-03 / 5.38e-04 / 1.29e-04 | 1.6, 2.1 |

**Exact vortex, M = 0.6:**

| Scheme | Pressure | Order | Tangential velocity | Order | Radial velocity | Order |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| default image point | 2.57e-02 / 1.43e-02 / 5.11e-03 | 0.8, 1.5 | 3.37e-02 / 1.97e-02 / 9.95e-03 | 0.8, 1.0 | 1.92e-02 / 9.45e-03 / 3.74e-03 | 1.0, 1.3 |
| TSC (`ib_ip_interp = 2`) | 6.49e-03 / 2.74e-03 / 1.20e-03 | 1.2, 1.2 | 3.59e-02 / 2.50e-02 / 1.54e-02 | 0.5, 0.7 | 2.21e-03 / 1.11e-03 / 3.62e-04 | 1.0, 1.6 |
| `ib_slip_ccst`, default image point | 4.16e-02 / 2.09e-02 / 8.09e-03 | 1.0, 1.4 | 2.23e-02 / 1.44e-02 / 8.80e-03 | 0.6, 0.7 | 2.46e-02 / 1.44e-02 / 5.41e-03 | 0.8, 1.4 |
| `ib_slip_ccst`, TSC | 5.04e-03 / 1.37e-03 / 3.17e-04 | 1.9, 2.1 | 1.70e-03 / 7.58e-04 / 2.74e-04 | 1.2, 1.5 | 2.91e-03 / 8.46e-04 / 2.21e-04 | 1.8, 1.9 |

## Assessment

- With the default `ib_mass_relax` a fixed slip wall gains gas: 0.36 % of the box mass per flow-through at 16 cells per R. `ib_mass_relax = 1` holds the mass to round-off (1e-13), and 30 to 6e-5. The drift falls slowly with resolution (3.6e-3 at 16 cells per R, 2.1e-3 at 32: a factor 1.7 for doubling), so the option is the remedy.
- The near-wall pressure error of the exact vortex falls at first order (1.1 to 1.5) with the default image point and the curvature term, and 1.9 to 2.3 with `ib_slip_ccst` and TSC weights. TSC alone lowers the pressure error 3 to 4 times at 40 cells per R; with `ib_slip_ccst` and TSC it is 27 times lower at M = 0.3 and 16 times lower at M = 0.6 than the default.
- The tangential velocity error stays at first order or below (0.5 to 1.0) with the default ghost copy and with TSC alone, and TSC makes it larger (1.4e-2 against 1.1e-2 at 40 cells per R, M = 0.3); only `ib_slip_ccst` with TSC removes it (2.5e-4 at 40 cells per R, 43 times lower).
- `ib_slip_ccst` with the default image point is worse in pressure than the curvature term alone on all three grids; use it with `ib_ip_interp = 2`.
- The same vortex runs reproduce, to three digits, the first-band errors measured for these options before they were merged (pressure 3.28e-2 / 1.58e-2 / 6.38e-3 at M = 0.3 with the default options), so the merged code matches the study that motivated it.
- Status is partially validated: only a fixed circle, steady inviscid flow and 2D were run; the options were not measured on a moving or accelerating slip body, on a sphere, or in the shock-curtain flow.

## Known limitations

- **Slip curtains are unsupported.** A shock through a 2D curtain of slip cylinders (about 115 um diameter, 16 cells per diameter, area fraction 0.21, Mach 1.66) fails shortly after the shock enters the curtain (two-way coupled moving particles); the same curtain with no-slip walls runs to the end. The likely cause is the inviscid slip model in 2D: the gas stays attached to each cylinder and expands at the shoulders toward vacuum, deeper as the grid is refined, where published inviscid slip-wall bed simulations \cite Mehta2019 are 3D (the flow can relieve sideways) and body fitted. This was not verified in MFC (an independent prototype of the 2D Euler problem keeps a positive minimum pressure, so it does not reproduce the failure itself). Use no-slip walls for curtains.
- **Slip particles in close approach.** A gap of under about two cells between two slip bodies (or a body and a wall) drains to vacuum (mechanism not isolated).
- **Unsteady flow with `ib_slip_ccst`.** Conserving total enthalpy between the image point and the ghost holds only for steady flow. A version of the option failed in three of three runs of a 2D shock curtain of about 40 slip particles (NaN partway through the run) where the zero-gradient slip fill passed four of four; the merged option was not run in that curtain. A nearly stagnant image point is given a bounded tangential speed (0.5 to 2 times the image point's), and the option is refused with chemistry.
- **Local conservation.** `ib_mass_relax` returns each body's net ghost-layer intake, so the body's total gas mass is conserved; the face fluxes that carry the intake are not corrected, and the correction is local to two cells of the body. A conservative face flux (cut cells) is not available. Count the mass of all fluid cells: the relaxation puts the returned mass in the two cells next to the body, so a count that leaves them out sees the drift as smaller than it is and the correction as incomplete.
- **NVIDIA GPU.** On an A100 (NVHPC 25.11 and 25.5, OpenACC, with and without MPI, one rank) `ib_mass_relax = 1` had no effect on the CI-sized mass budget: the mass change equals the run without it to every printed digit (+1.6e-3 in 100 steps at 8 cells per R, where a CPU run gives -4e-15), also with extrapolation boundaries and at 16 cells per R. With `ib_state_wrt = T` added to the same case the A100 gives -3.7e-15, the CPU value, so the failure is specific to the relaxation-only path (the one taken when no body moves and neither `ib_state_wrt` nor `ib_heat_wrt` is set). The case documentation of `ib_mass_relax` records 2e-13 on an A100 and an MI300A; that was not reproduced. The cause was not traced. Reproducer: `./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 --gpu acc -- --mass_relax 1` compared with `--mass_relax 0` (identical on an A100), then `check_slip.py`. The suite does not register the `--mass_relax 1` case for this reason; `IBM -> Fixed Circle -> Slip -> ib_mass_relax` passes on the A100 with and without the option and does not discriminate.
- **Curvature.** The curvature term is exact for circles, spheres and the sides of cylinders; other shapes have none and keep the zero-gradient ghost pressure.
- Not re-measured here: the pressure of an accelerating slip cylinder (added mass), the uniform-flow pressure coefficient of a slip circle and the trajectory of a free slip cylinder behind a shock, which in the earlier study did not separate the options; and the comparison with the code before the curvature term, which cannot be made on `master`.

## Reproduce

```
./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 -- --steps 0 --cells_per_R 16 --half 8 --mass_relax 1
python3 examples/2D_ibm_slip_cylinder/check_slip.py examples/2D_ibm_slip_cylinder --case-args '--steps 0 --cells_per_R 16 --half 8 --mass_relax 1'
./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 -- --mode vortex --mach 0.3 --cells_per_R 20 --ip_interp 2 --ccst
python3 examples/2D_ibm_slip_cylinder/check_slip.py examples/2D_ibm_slip_cylinder --case-args '--mode vortex --mach 0.3 --cells_per_R 20 --ip_interp 2 --ccst'
```

`--mass_relax`, `--ip_interp`, `--ccst`, `--cells_per_R` and `--mach` select the cases of the tables. The vortex mode compiles its own
`pre_process` and `simulation` (the initial condition is an analytic patch). The suite runs the CI-sized mass budget (8 cells per R, a
box of side 12 R, 100 steps, default `ib_mass_relax`) as `2D -> IBM -> Slip Cylinder -> Mass Budget`.
