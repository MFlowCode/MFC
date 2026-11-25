# 1D circularly polarized Alfven wave

The smooth Alfven wave of Tóth (2000, section 6.3.1): $\rho = 1$, $p = 0.1$, $\gamma = 5/3$, $B_x = 1$,
$v_\perp = B_\perp = 0.1(\sin 2\pi x, \cos 2\pi x)$ in a periodic unit box (`hcid = 150`), HLLD. It is an exact nonlinear
solution of ideal MHD: the wave travels at the Alfven speed $B_x/\sqrt{\rho} = 1$ in $-x$ with $\rho$, $p$ and
$|B_\perp|$ constant, so at t = 1 the state equals the initial one.

`check_alfven.py` compares the last save of each run with the initial state shifted by $-t$ and prints the L1 errors of
$v_y$, $v_z$, $B_y$, $B_z$, the largest deviation of $\rho$ and $p$, and the observed order of the $B_y$ error over a sweep.

## Run

    ./mfc.sh run examples/1D_mhd_smooth_alfven_wave/case.py -t pre_process simulation post_process -n 1
    python3 examples/1D_mhd_smooth_alfven_wave/check_alfven.py examples/1D_mhd_smooth_alfven_wave

The default is the original deck (1024 cells, WENO3, dt = 1e-5, 100000 steps). `case.py` also takes `-N` (cells),
`--cfl` (fast-wave CFL number, which sets dt to end exactly at t = 1) and `--weno` (`weno_order`); the sweep on
@ref validation_mhd uses `-N 32 ... 512 --cfl 0.4` with `--weno 3` and `--weno 5`, one copy of the directory per run.
