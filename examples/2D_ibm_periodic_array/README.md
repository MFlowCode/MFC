# Periodic Arrays of Fixed Particles: Drag, Force Spread and Nusselt Number

A random array of fixed no-slip isothermal bodies (2D cylinders, or spheres with `--dim 3`) in a periodic box, driven at a held
bulk mass flux by the mean pressure gradient controller (`mean_mflux`). The case records the force on every body
(`D/ib_forces.dat`), the heat it takes (`D/ib_heat.dat`) and the controller state (`mean_dpdx.dat`). `--thermostat` adds the thermally
developed controller of Sun et al. (2015) (`mean_temp_wall`) for the Nusselt number. Gas: dry air with `h2o2.yaml` (thermodynamics and
transport only, no reactions). Results, references and limits are on the @ref validation_particle_array page.

## Run

    ./mfc.sh run examples/2D_ibm_periodic_array/case.py -n 2 --gpu
    python3 examples/2D_ibm_periodic_array/check_array.py examples/2D_ibm_periodic_array

Without arguments the case is the CI-sized run: 12 cylinders at a solid fraction of 0.2, Re_m 20, 14 cells per diameter
(96 x 96 cells), 200 fixed steps on 2 ranks, far from steady (`2D -> Chemistry -> IBM -> Periodic Array`; its golden holds the
conserved fields, the force and the heat histories). A drag run needs `--tstar 40` or more (about 12000 steps), and the 3D box
`--dim 3 --n_bodies 20 --mach 0.3 --cfl 0.2 --cells_per_D 16 --tstar 10`. The placement is a seeded random sequential
addition (`--seed`), identical on every Python version. `check_array.py` needs no post_process output.

## Results

NVIDIA A100 80GB, NVHPC 25.11, double precision, `master` at 0cdfdf02a; 20 spheres at a solid fraction of 0.2 in a cube of 3.74
diameters, `ib_ip_interp = 1`, last 40 % of `tW/D` = 10. F/F_St is the total force per sphere over 3 pi mu D (1 - phi) W
(Tenneti et al. 2011, Eqs. 33 and 38); the reference is their Eq. 43 at the realized Re_m.

| Re_m | seed | cells/D | F/F_St | Tenneti | error | sum(F_x) / (-dp/dx V_f) |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 1 | 16 | 7.075 | 7.366 | -4.0 % | 1.0001 |
| 20 | 2 | 16 | 7.354 | 7.355 | -0.0 % | 1.0001 |
| 100 | 1 | 16 | 13.973 | 12.621 | +10.7 % | 1.0001 |
| 100 | 2 | 16 | 14.337 | 12.605 | +13.7 % | 1.0001 |
| 100 | 1 | 24 | 13.894 | 12.606 | +10.2 % | 1.0000 |

The 2D run of the CI case to `tW/D` = 43 holds Q to 1e-6, closes the momentum balance to 1.0000 and gives a per-body force spread of
0.44 of the mean.

Nusselt number (`--thermostat`, 20 spheres, Re_m 20, 16 cells/D, Mach 0.1, bodies at 290 K and 310 K, `tW/D` = 8): 6.11 (cold) and
5.71 (hot) against Gunn 5.995; the pair mean 5.91 is 1.4 % below Gunn and 1.7 times the Sun et al. (2015) value on the same basis. The hot/cold
asymmetry (7 %, 28 % before the thermostat removed the viscous dissipation) is an open item, so treat the values as plus or minus 15 %; see the @ref validation_particle_array page.
