# Isentropic vortex through AMR blocks

Shu's isentropic vortex (\f$\epsilon = 5\f$, \f$\gamma = 1.4\f$) carried by a uniform stream \f$u = 1\f$ across a periodic
\f$[-10, 10]^2\f$, from \f$x = -3\f$ to \f$x = 3\f$ at \f$t = 6\f$. The exact solution is the initial vortex translated, so the
density error at the end measures the accuracy of the coarse/fine coupling for a smooth flow. Every level starts from the exact
cell averages of the conserved variables (3x3-point Gauss in the patch expressions, which use `dx`, `dy`).

- `--mode static`: one 2:1 block over \f$[-1.25, 1.25] \times [-2.5, 2.5]\f$. The vortex enters it, crosses it and leaves it,
  so every seam is crossed by the swirl.
- `--mode regrid`: the blocks follow the vortex (density-jump tags down to \f$r \approx 2.5\f$, regrid every 10 steps).
- `--mode uniform`: no AMR at `--N` cells, or at 2N with `--fine`.
- `--prolong 1|2`: linear or quartic ghost fill and new-block fill; `--subcycle`: two fine steps per level-0 step.

All runs at one `--N` take the same fine-level step (CFL 0.4 on the 2N grid), so a uniform run at 2N, an AMR run and a uniform
run at N differ only in space.

`check_vortex.py` reads the last save (`restart_data/lustre_<step>.dat` and the fine blocks in `restart_data/lustre_amr_<step>.dat`;
`parallel_io = T`) and compares the density with the exact cell averages (4x4-point Gauss) on the composite grid: the fine cells where a block
covers a point, level 0 elsewhere. It prints L1 (area mean) and Linf, the largest error within 3 level-0 cells of a block edge
("seam") and elsewhere ("away"), the refined share of the domain, and the observed order between successive `--N` of one kind.

```
./mfc.sh run examples/2D_amr_isentropic_vortex/case.py -t pre_process simulation -n 1 -- --N 128 --mode static --prolong 2
python3 examples/2D_amr_isentropic_vortex/check_vortex.py <run directories in increasing N>
```

Each run writes next to `case.py`: copy the directory per run. Results and limits: the AMR section of the validation page.
