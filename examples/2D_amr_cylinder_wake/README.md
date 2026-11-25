# 2D Cylinder Wake with AMR: Drag and Strouhal Number, Re 100

A no-slip immersed cylinder in a low-Mach stream (U = 0.2, c = 1, M = 0.2, Re = 100 on the diameter D) sheds a von Karman
street. One static 2:1 AMR block covers the body and the near wake (-1.5 D <= x <= 6 D, |y| <= 2.5 D); the rest of the
36 D x 32 D domain (3.1 % blockage) stays at 16 cells per diameter. The case checks that the block reproduces the
uniform 32-cells-per-diameter grid in mean drag and shedding frequency, and both against published values. It is also the
body-wake row of "When AMR pays" in `docs/documentation/amr.md`.

## Run

    ./mfc.sh run examples/2D_amr_cylinder_wake/case.py -t pre_process simulation -n 1 --gpu
    python3 examples/2D_amr_cylinder_wake/check_wake.py examples/2D_amr_cylinder_wake [<uniform fine run>]

Variants are edits at the top of `case.py`: `amr = False` gives a uniform grid at `N` cells per diameter (`N = 32` for the
uniform fine grid), `subcycle = True` gives `amr_subcycle`. `dt` follows the finest cells and every variant runs to
t = 650 (130 D/U). `check_wake.py` reads `D/ib_forces.dat` and prints the mean Cd and the Strouhal number over the
whole shedding periods in the last 250 time units (`--window`), with signed errors against Henderson (1995), Park et al.
(1998) and Williamson's St(Re) fit, and against a second run directory when one is given. St is from the mean spacing of
the upward zero crossings of Cl; an FFT of the same window agrees to 0.2 %.

The case is not in the test suite (a run is 130 D/U long); the AMR immersed-boundary path is covered by `amr_golden_tests`.

## Results

One MI250X GCD (amdflang, OpenMP offload), master at 5838b641c. Window t = 400-650, 8 shedding periods.

| Run | Grid | Cd | St | Cd vs uniform fine | St vs uniform fine |
| :--- | :--- | ---: | ---: | ---: | ---: |
| uniform fine | 1152 x 1024 (32/D) | 1.3687 | 0.1674 | | |
| uniform coarse | 576 x 512 (16/D) | 1.3672 | 0.1647 | -0.11 % | -1.60 % |
| AMR, lock-step | 576 x 512 + block 240 x 160 | 1.3682 | 0.1673 | -0.03 % | -0.03 % |
| AMR, subcycled | same | (1.6970) | 0.1673 | | -0.03 % |

With `amr_subcycle` the force is summed on the coarse faces only: the subcycled Cd is 24 % above the lock-step one and is
not a drag measurement; its St, from the same flow, is. Against the references, the uniform fine grid is +1.4 % (Henderson, Cd 1.35) and +2.9 % (Park, 1.33) in
Cd and +1.9 % (Williamson, St 0.1643) and +1.5 % (Park, 0.165) in St. On a 24 D x 16 D domain (6.25 % blockage) the
same grids gave Cd and St 3 to 4 % higher; the AMR run still matched the uniform fine one (-0.31 %, -0.17 %).
