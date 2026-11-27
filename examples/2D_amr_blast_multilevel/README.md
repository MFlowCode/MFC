# 2D Circular Blast with Two AMR Levels

A disc of gas at $(\rho, p) = (1, 1)$, radius 0.1, expands into gas at $(0.125, 0.1)$ (the 2D Sod explosion). By default the
128² grid carries two nested 2:1 levels on the shock and the contact (see @ref amr); `--no-amr -N <N>` runs the same blast on a
uniform $N^2$ grid to the same time, $t = 0.15625$.

`check_blast.py` solves the radially symmetric problem on 16000 cells (second-order finite volume) and prints, per run, the L1
density error against it and the shock radius:

    ./mfc.sh run examples/2D_amr_blast_multilevel/case.py -t pre_process simulation -n 1
    python3 examples/2D_amr_blast_multilevel/check_blast.py <run> [<run> ...]

The AMR run is as accurate as a uniform 512² grid (L1 5.5e-4 against 5.9e-4) and the shock radius agrees with the reference within
a tenth of a cell. Results and limits (the AMR run aborts on an NVIDIA GPU, #2467): @ref validation_blast.
