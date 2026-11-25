# 2D Conservation with AMR in a Periodic Box

A square of dense, high-pressure gas (rho = p = 2 in a background of 1, uniform flow (0.5, 0.25)) in a periodic unit box with two refined 2:1 levels that regrid every few steps. The density jump advects; the pressure jump sends shocks and rarefactions across block edges and the periodic boundaries. With periodic boundaries nothing enters or leaves, so the total mass, momentum and energy must stay at their initial values to round-off. `case.py` uses mapped WENO5, HLLC and SSP-RK3; the square's edges lie on base-grid faces, so the blocks start from exactly the level-0 state.

    ./mfc.sh run examples/2D_amr_conservation/case.py -t pre_process simulation -n 1 -- --subcycle
    python3 examples/2D_amr_conservation/check_conservation.py examples/2D_amr_conservation

`case.py` takes `-N` (base cells per side, default 64), `--levels` (default 2), `--subcycle`, `--steps` (default: to t = 0.4) and `--saves`. `check_conservation.py` reads the serial restart files of every save and prints the relative drift (S(t) - S(0))/S(0) of each total, summed exactly (`math.fsum`), and the gap between the level-0 totals and the totals over the composite grid (leaf cells of every level, from `amr_fine.dat`). It exits 1 if any drift exceeds `--tol` (1e-13).

Measured on `master` after #2482, #2474 and #2461, AMD flang on a Frontier CPU node, double precision, 64^2 base, t = 0.4 (up to 52544 refined cells):

| Mode | mass | x-momentum | y-momentum | energy | level 0 - composite |
| :--- | ---: | ---: | ---: | ---: | ---: |
| lock-step, 684 steps | -6.6e-16 | 0 | 0 | -1.6e-16 | 4.4e-16 |
| subcycled, 172 steps | 0 | 0 | +2.2e-16 | +1.6e-16 | 2.2e-16 |

The largest drift over the four saves is one to three units in the last place of the totals.

The suite registers `AMR -> 2D -> Example -> conservation subcycle` (32^2 base, 40 steps). Results and limits: @ref validation_amr.
