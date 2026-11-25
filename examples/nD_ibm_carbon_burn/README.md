# Burning Carbon Cylinder and Sphere Against a 1D Radial Reference

A fixed carbon body (D = 100 um) held at 2000 K gasifies in quiescent moist air at 10 atm and 2000 K, with a CO flame off the surface (the mechanisms of `2D_ibm_reacting_surface`: 11-species reduced GRI gas phase and the carbon surface mechanism). The body appears at t = 0 in a box of half-width 6 R with extrapolation faces; the surface chemistry runs with `ib_flux_wall` and `ib_heat_wrt` records the burn rate. `--ndim 2` is a cylinder, `--ndim 3` a sphere. The reference is a 1D transient solution of the same model (`ref1d.py`: finite volumes with geometric stretching, Cantera thermodynamics, kinetics and mixture-averaged transport, Stefan velocity from the ideal-gas constraint, the wall receiving the surface-solve fluxes); its burn rate and flame peak are tabulated in `ref_cylinder.txt` and `ref_sphere.txt` (N = 400 cells, outer radius 200 R).

    ./mfc.sh run examples/nD_ibm_carbon_burn/case.py -n 1 -t pre_process simulation post_process -- --cells_per_D 24 --tend 6e-5
    python3 examples/nD_ibm_carbon_burn/compare_burn_rate.py examples/nD_ibm_carbon_burn
    python3 examples/nD_ibm_carbon_burn/ref1d.py --geom 1 --tend 6e-5 --out ref_cylinder.txt   # regenerates the reference (43 min on one CPU core; it reproduces the shipped table to 5e-7)

`compare_burn_rate.py` prints the burn rate (`D/ib_heat.dat` column 4, the face mass the flux wall passes to the gas, and column 3, the surface-solve rate) against the reference at 10, 20, 30 and 60 us, and the flame peak at the last save. `case.py` takes `--ndim`, `--D`, `--p_atm`, `--T`, `--H_over_R`, `--cells_per_D`, `--tend`, `--frames`, `--cfl` and `--set KEY=VAL`.

NVIDIA A100 (NVHPC 25.11), error = (MFC - reference)/reference of the burn rate (`ib_heat.dat` column 4) at 10 and 20 us. "Shared" is `master` at 0cdfdf02a, which gave all species one ghost mirror weight; "per species" is the current rule (each species its own weight, the most abundant closing the sum). The recommended setting is `--set ib_ip_interp=2` (triangular-shaped-cloud image points); the suite test keeps the default.

| Body | cells per D | image points | shared | per species |
| :--- | ---: | :--- | ---: | ---: |
| cylinder | 12 | default | +1.87 %, +1.36 % | +1.87 %, +1.35 % |
| cylinder | 12 | `ib_ip_interp=2` | +1.39 %, +0.80 % | +0.90 %, +0.07 % |
| cylinder | 24 | default | +0.86 %, +0.64 % | not run |
| cylinder | 24 | `ib_ip_interp=2` | +0.32 %, +0.16 % | +0.15 %, -0.11 % |
| sphere | 12 | default | +7.35 %, +6.22 % | +7.34 %, +6.20 % |
| sphere | 12 | `ib_ip_interp=2` | +3.17 %, +2.49 % | +2.07 %, +1.15 % |

With the per-species weight and `ib_ip_interp=2` the runs reproduce the earlier pre-merge study (+0.92 and +0.08 % for the 12-cell cylinder, +0.16 and -0.11 % at 24, +2.07 and +1.14 % for the sphere). Compare at times at least 0.5 us before the end of the run: the burn rate falls steeply and the averaging window of `compare_burn_rate.py` is cut off at the last row. See the page below for the settings tried.

The suite registers `2D -> Chemistry -> IBM Carbon Cylinder -> Burn Rate` (12 cells per diameter, 5e-7 s). Results and limits: @ref validation_burning_body.
