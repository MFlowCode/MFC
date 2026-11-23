# 2D Moving Immersed Boundary: Gas-Mass Budget and Galilean Twin

A circle translates through gas at rest in a doubly periodic box (non-dimensional ideal gas, no chemistry), and the same circle is
held while the gas streams past it at the same relative speed (the Galilean twin). In the closed box the gas mass is constant, so
whatever the box gains is created by the body: a ghost-cell immersed boundary does not conserve mass at a moving wall, because gas
flows through the faces between the fluid and the body into ghost cells whose values are then overwritten, and cells change side as
the body moves. `ib_mass_relax` (default 30 for a moving body) takes that mass back. The twin checks the conducted heat and the force
of the moving body against the held one. The full results and limits are on the @ref validation_moving_body_budget page.

`check_budget.py` reads the run directory:

    check_budget.py mass RUN ...                       gas mass (cells outside the body) at every save, and the steady
                                                       rate in % of the gas per diameter travelled after the first diameter
    check_budget.py twin --moving RUN --fixed RUN      moving/fixed ratio of the conducted heat (from the temperature field) and of
                                                       F_x (ib_forces.dat), and the RMS of F_moving(t) - F_fixed(t)

## Run

    ./mfc.sh run examples/2D_ibm_moving_mass_budget/case.py -t pre_process simulation -n 1 --gpu -- \
        --cells_per_D 24 --mach 0.3 --box 5 --travel 6 --nsave 24 --relax 0
    python3 examples/2D_ibm_moving_mass_budget/check_budget.py mass examples/2D_ibm_moving_mass_budget

The output goes next to `case.py`, so copy the directory for each run. `--relax 30` gives the correction; omitting `--relax` leaves
the solver default. The twin pair (cold wall, Re 40):

    ... -- --cells_per_D 16 --mach 0.3 --box 8 --Re 40 --Twall 0.2 --ip_interp 2 --travel 6 --nsave 12 --stride 1
    ... -- (the same) --frame fixed
    python3 examples/2D_ibm_moving_mass_budget/check_budget.py twin --moving <moving run> --fixed <fixed run>

Without arguments the case is the CI-sized test: 24 cells per diameter in a box of 5 D, 100 fixed steps (a quarter of a diameter of
travel). It runs in the suite as `2D -> IBM -> Moving Circle -> Mass Budget` (and `-> Cold Wall`, the viscous isothermal twin).
A build without MPI needs `--serial_io`; the checker then reads `p_all/`.

## Results

GNU 13 CPU, double precision, `master` at 0cdfdf02a. Inviscid, Mach 0.3, box 5 D, 6 D travelled, 24 saves, inverse-distance image
points (`ib_ip_interp = 1`). Steady rate of the gas mass, % of the gas per diameter travelled (least squares after the first D),
and the total change after 6 D:

| cells/D | `ib_mass_relax = 0`: rate | total | `ib_mass_relax = 30`: rate | total |
|---:|---:|---:|---:|---:|
| 16 | +0.0503 | +0.357 % | +0.0022 | +0.016 % |
| 24 | +0.0237 | +0.192 % | -0.0008 | +0.002 % |
| 32 | +0.0130 | +0.134 % | -0.0023 | +0.003 % |

With TSC image points (`ib_ip_interp = 2`): +0.0544 and +0.0166 %/D without the correction at 16 and 32 cells/D, +0.0024 and
-0.0022 with it. The uncorrected drift falls with refinement at roughly first order; the correction holds it to about 0.002 %/D.

Galilean twin, viscous (Re 40, Pr 0.7), Mach 0.3, box 8 D, TSC, 6 D travelled, means over the last 3 D. Moving over fixed:

| wall T / gas T | cells/D | heat, `ib_mass_relax = 0` | heat, 30 | F_x, 0 | F_x, 30 |
|---:|---:|---:|---:|---:|---:|
| 0.9 | 16 | 0.940 | 0.993 | 0.966 | 0.963 |
| 0.9 | 32 | 0.966 | 1.002 | 0.969 | 0.968 |
| 0.2 | 16 | 0.760 | 0.817 | 1.248 | 1.748 |
| 0.2 | 32 | 0.882 | 0.917 | 0.690 | 0.682 |

The heat of a wall close to the gas temperature agrees with its twin within 1 % with the correction; the force differs by 3 %. A
cold wall (0.2 of the gas temperature) is not like for like: the conducted heat is 18 % low at 16 cells/D and 8 % at 32, and the
force of the moving body differs from its twin by tens of percent with no consistent sign. See the page for the causes.

## Limits

The box is closed and the body is a prescribed (kinematic) circle: this checks conservation and Galilean invariance, not a drag
or heat-transfer value against a measurement. Not run here: reacting moving bodies and the 3D particle array (the examples
`2D_ibm_wall_mass_budget` and `2D_ibm_reacting_surface` hold the fixed-body budget).
