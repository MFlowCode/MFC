# Smooth 1D waves: order of accuracy from exact cell averages

Two smooth periodic Euler problems that separate the reconstruction order from the primitive-variable conversion
(@ref validation_vortex_order):

- `--problem advect`: a density wave $\rho = 1 + 0.2\sin 2\pi x$ carried at u = 1 through p = 1 ($\gamma$ = 1.4) for one
  period. u and p are uniform, so the primitives are linear in the conserved variables.
- `--problem gamma3`: the isentropic $\gamma = 3$ wave $\rho_0 = 1 + 0.2\sin 2\pi x$, $u_0 = 0$, $p = \rho^3/3$, to t = 0.2.
  $u \pm c$ each obey Burgers' equation, so the exact solution follows from the characteristics.

`--ic avg` (default) starts from the exact conservative cell averages: `case.py` writes them to `IC/` next to itself and
`hcid = 170` reads them, so each run needs its own copy of this directory. `--ic point` starts from point values. Other
options: `-N`, `--order` (3, 5, 7), `--teno`, `--cfl` (0.02; use 0.005 for seventh order) and `--t-end`.

    for N in 32 64 128 256; do
      cp -r examples/1D_smooth_wave_order run$N
      ./mfc.sh run run$N/case.py -t pre_process simulation post_process -n 1 -- -N $N --problem gamma3
    done
    python3 examples/1D_smooth_wave_order/check_smooth.py run32 run64 run128 run256 \
        --case-args "-N 32 --problem gamma3;-N 64 --problem gamma3;-N 128 --problem gamma3;-N 256 --problem gamma3"

The linear wave reaches the design order (5.0 for WENO5 and TENO5, 7.0 for WENO7 and TENO7). The $\gamma = 3$ wave is second
order for every scheme, because MFC reconstructs primitive variables from cell averages. The tables are on
@ref validation_vortex_order. The example is not in the test suite.
