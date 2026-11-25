# 2D Orszag-Tang vortex with hyperbolic divergence cleaning

The Orszag-Tang vortex ($\gamma = 5/3$, `hcid = 250`) on a periodic unit square, 513 x 513 cells, WENO5-M, HLL, to t = 1,
with `hyper_cleaning` (speed 2.5, damping time 0.004).

`check_divb.py` prints, for each save, the normalized divergence error $\overline{|\nabla\cdot B|}\,\Delta x/\overline{|B|}$
(second-order central differences), its maximum and the change of the total mass; given a second run of the same deck with
`--no-clean`, it prints the two errors side by side. Without cleaning this deck stops with NaN at step 900 (t = 0.36), so
run the control with `--steps 800`.

## Run

    ./mfc.sh run examples/2D_orszag_tang_hyper_cleaning/case.py -t pre_process simulation post_process -n 1
    python3 examples/2D_orszag_tang_hyper_cleaning/check_divb.py examples/2D_orszag_tang_hyper_cleaning [<control run>]

Results: @ref validation_mhd.
