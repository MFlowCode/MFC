# 2D Circular Cylinder, Re 20-100: Drag, Wake, Shedding and Nusselt Number

A fixed no-slip isothermal cylinder (immersed boundary, `thermal_bc = 1`) in dry air at Mach 0.2. The case measures four
things against literature values: the mean drag coefficient, the steady wake length, the Strouhal number of the shedding wake,
and the Nusselt number of a heated wall. The drag is taken two ways that share no code: the face-flux force MFC writes with
`ib_force_wrt` (`D/ib_forces.dat`), and a momentum balance on squares of half-width 1, 2 and 4 diameters around the body, from
the saved fields. Chemistry is on for the thermodynamics and mixture transport only (`h2o2.yaml`, no fuel, reactions off).
The full results and limits are on the @ref validation_cyl_drag_nu page.

## Run

    ./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation post_process -n 1 --gpu -- \
        --silo --Re 40 --cells_per_D 20 --xu 15 --xd 25 --hy 16 --steps 0 --tstar 50
    python3 examples/2D_ibm_cylinder_drag_nu/check_cylinder.py examples/2D_ibm_cylinder_drag_nu [--hot <heated run>]

`--tw_ratio 1.1` gives the heated twin whose directory goes to `--hot` (the Nusselt number uses the difference of the two, which
removes the viscous heating of the wall). `--steps 0` selects the adaptive time step and `--tstar` the run length in
`t U / D`; Re 100 needs `--tstar 120` (periodic shedding is reached by about 90). Without arguments the case is the CI-sized
test: Re 20, 16 cells per diameter, a 12 D x 8 D domain, 200 fixed steps. It runs in the suite as
`2D -> Chemistry -> IBM -> Cylinder Drag Nu` (golden: the force and heat histories and the conserved fields).

## Results

AMD MI300A, OpenMP offload, double precision, `master` at 615348a25; 20 and 10 cells per diameter, domain
[-15 D, 25 D] x [-16 D, 16 D] (3 % blockage), `ib_ip_interp = 1`, adaptive step. Errors are against Dennis & Chang for
Re 20 and 40 and Henderson (mean Cd 1.35) for Re 100.

| Re | cells/D | Cd_IB | Cd_CV (2 D) | L/D | Nu |
|---:|---:|---:|---:|---:|---:|
| 20 | 10 | 1.982 | 1.983 | 0.721 | not run |
| 20 | 20 | 2.029 (-1.0 %) | 2.029 | 0.818 (ref 0.94) | 2.442 (Lange +1.4 %) |
| 40 | 10 | 1.500 | 1.498 | 1.932 | not run |
| 40 | 20 | 1.521 (+0.1 %) | 1.519 | 2.064 (ref 2.35) | 3.243 (Lange -1.2 %) |
| 100 | 20 | 1.354 (+0.3 %) | 1.332 | shedding | 4.962 (Lange -3.2 %) |

Re 100 used a smaller domain ([-8 D, 16 D] x [-8 D, 8 D], 6 % blockage) and `t U / D` = 120: Cl amplitude 0.287 (references
0.33-0.34) and St = 0.1668 (+1.5 % against Williamson's 0.1643). The Re 20 and 40 wakes are steady. Nu is the twin-corrected
face-heat value at Tw/T_inf = 1.1, with the film conductivity from Cantera at Pr = 0.712. Not run: 40 cells per diameter (about
5 GPU-hours per run on this domain).

Two limits are visible in the table. L/D is first order and low (the ghost-cell wall sits up to a cell inside the surface), and
Cd_CV depends on the box at Re 100 in the smaller domain (1.344, 1.332, 1.302 for 1, 2 and 4 D).
