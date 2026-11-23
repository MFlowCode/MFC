# 3D Sphere, Re 20-300, Mach 0.2: Drag and Nusselt Number

A fixed no-slip isothermal sphere (immersed boundary, `thermal_bc = 1`) in dry air at Mach 0.2: the 3D companion of
`2D_ibm_cylinder_drag_nu`. It measures the mean drag coefficient and the Nusselt number of a heated wall against the
Schiller-Naumann drag correlation and the Ranz-Marshall and Whitaker Nusselt correlations. The drag is taken two ways that share no
code: the face-flux force MFC writes with `ib_force_wrt` (`D/ib_forces.dat`), and a momentum balance on cubes of half-width 1 and 2
diameters around the body, from the saved fields. Chemistry is on for the thermodynamics and mixture transport only (`h2o2.yaml`, no
fuel, reactions off). The full results and limits are on the @ref validation_sphere_drag_nu page.

## Run

    ./mfc.sh run examples/3D_ibm_sphere_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
        --silo --Re 100 --cells_per_D 12 --xu 3 --xd 6 --hy 3 --steps 0 --tstar 40
    python3 examples/3D_ibm_sphere_drag_nu/check_sphere.py examples/3D_ibm_sphere_drag_nu [--hot <heated run>]

`--tw_ratio 1.1` gives the heated twin whose directory goes to `--hot` (the Nusselt number uses the difference of the two, which
removes the viscous heating of the wall). The output goes next to `case.py`, so copy the directory for each run. `--steps 0` selects the adaptive time step
and `--tstar` the run length in `t U / D`. Without arguments the case is the CI-sized test: Re 20, 6 cells per diameter, a 4.2 D cube
(25 x 25 x 25 cells), 30 fixed steps; it runs in the suite as `3D -> Chemistry -> IBM -> Sphere Drag Nu` (golden: the force and heat
histories and the conserved fields). The test is a smoke-level body 6 cells across, not a resolved sphere.

## Results

A100, NVHPC 25.11, double precision, `master` at 860663835. Domain [-3 D, 6 D] x [-3 D, 3 D]^2 (2.2 % area blockage),
`ib_ip_interp = 1`, adaptive step, 12 and 8 cells per diameter. Cd_IB is the face-flux force, Cd_CV the momentum balance on the cube of half-width 2 D (1 D in the
page), SN Schiller-Naumann (2.610 at Re 20, 1.092 at Re 100), RM Ranz-Marshall (4.396, 7.359), W Whitaker (3.948, 6.621).

| Re | cells/D | Cd_IB | vs SN | Cd_CV (2 D) | L/D | Nu | vs RM | vs W |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 8 | 2.532 | -3.0 % | 2.534 | 0 | 4.083 | -7.1 % | +3.4 % |
| 20 | 12 | 2.617 | +0.3 % | 2.618 | 0 | 4.059 | -7.7 % | +2.8 % |
| 100 | 8 | 1.114 | +2.1 % | 1.102 | 0.755 | 6.939 | -5.7 % | +4.8 % |
| 100 | 12 | 1.101 | +0.9 % | 1.100 | 0.795 | 6.905 | -6.2 % | +4.3 % |

L/D at Re 100: Johnson and Patel 0.88; no recirculation at Re 20 (below Re of about 24). The Nusselt number lies between the two
correlations, which differ by 11 % from each other. Not run here: Re 300 (the wake needs t U/D of about 100), 20 cells per diameter
(about 8 times the cost of 12), and the Mach 0.3-2 rows; the 8-cell runs stop at t U/D 12 (Re 100: not fully settled, Cd still falling by
0.2 % per two time units).

## Limits

The ghost value of a no-slip wall is the wall value (first order, [#2389](https://github.com/MFlowCode/MFC/issues/2389)), so the
wake is short and the Nusselt number may drift with refinement; the 2.2 % blockage raises Cd slightly; the correlations are
incompressible fits of 5 % (drag) and 11 % (between the two Nusselt correlations). A 12 cells/D run takes 2 to 3 hours on one shared A100.
