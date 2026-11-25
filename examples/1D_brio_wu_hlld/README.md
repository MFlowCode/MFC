# 1D Brio-Wu MHD shock tube (HLLD)

The coplanar MHD Riemann problem of Brio and Wu (1988): left $(\rho, p, B_y) = (1, 1, 1)$, right $(0.125, 0.1, -1)$,
$B_x = 0.75$, $\gamma = 2$, at rest, on $[0, 1]$ to t = 0.2, WENO5-M, HLLD (`1D_brio_wu` is the same deck with HLL). The
solution has a fast rarefaction, a slow compound wave (an intermediate shock followed by a slow rarefaction), a contact, a slow
shock and a fast rarefaction. With the compound wave the Riemann problem has no unique exact solution, so the check has two parts:

- grid convergence against a fine run of the same deck (L1 errors of $\rho$, $u$, $v$, $p$, $B_y$ and the observed order);
- with `--gamma 5/3` (the variant of Ryu and Jones 1995, section 5), the states on both sides of the intermediate shock against
  the values they read from an 8192-cell calculation.

## Run

    ./mfc.sh run examples/1D_brio_wu_hlld/case.py -t pre_process simulation post_process -n 1 -- -N 1602
    python3 examples/1D_brio_wu_hlld/check_brio_wu.py --ref <run with -N 25632> <run with -N 1602> [--rj]

`case.py` takes `-N` (cells, default 801; dt scales with 801/N) and `--gamma`. The reference's cell count must be an integer
multiple of each run's; it is block-averaged onto the coarser grids. Results: @ref validation_mhd.
