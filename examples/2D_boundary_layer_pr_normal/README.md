# Standalone normal pr boundary layer

All case parameters and thermodynamic calculations are in this folder's `case.py`.
It does not import another example or a shared case file. Edit `NX`, `NY`, `SAVES`,
`FLOW_TIMES`, `TEMPERATURE`, `PRESSURE`, and `WALL_TEMPERATURE` there.

Defaults: nitrogen at 300 K and 101325 Pa, isothermal wall at 270 K,
96 × 64 cells, Mach 0.2, Re_L = 10000, and four flow-through times. Viscosity and
conductivity are prescribed constants. Plate length is calculated to match Re_L;
therefore the three examples have different physical lengths. The wall-normal
grid is stretched. This is a coarse boundary-layer smoke test; quantitative wall
heat flux requires a grid-convergence study and appropriate transport properties.

## Run independently

From the MFChaocheyfork repository root (four MPI ranks):

```bash
./mfc.sh run examples/2D_boundary_layer_pr_normal/case.py -t pre_process simulation --mpi --no-gpu -n 4 -j 4
```

Your current `parallel_io = "T"` preference is preserved as `PARALLEL_IO = "T"`.
For the previously tested serial mode, set `PARALLEL_IO = "F"` in this file and
replace `--mpi -n 4` with `--no-mpi` in the command. PR remains CPU-only. The prior
physical runs were validated in serial mode; MPI execution has not been validated
as part of this folder reorganization.

## Save frames and make animations

Set `SAVES = 80` (or another desired count) **before running** for approximately
81 saved frames including the initial state. The timestep and stop/save steps
are calculated automatically. The default `SAVES = 8` is a small verification run.
Post-process using the same case settings and I/O mode as the simulation:

```bash
./mfc.sh run examples/2D_boundary_layer_pr_normal/case.py -t post_process --mpi --no-gpu -n 4 -j 4
./mfc.sh viz examples/2D_boundary_layer_pr_normal/ --list-vars --step 0
./mfc.sh viz examples/2D_boundary_layer_pr_normal/ --step all --mp4 --fps 12
```

The last command renders the available variables as a tiled animation. Add
`--var NAME` using a name from `--list-vars` to animate one field. MFC's MP4 encoder
uses `imageio` and `imageio-ffmpeg` (which bundles ffmpeg). For serial output, use `--no-mpi` instead of `--mpi -n 4` for
post-processing too. Keep case settings consistent with the saved data.

## Verify and plot

```bash
build/venv/bin/python examples/2D_boundary_layer_pr_normal/verify.py
```

The local verifier reads this folder only and writes profiles and metrics to
`viz/`. It supports double-precision serial output and combined MPI-I/O restart
files (`file_per_process = F`, the default). Its MPI file reader has a layout
check, but is not a substitute for validating an actual multi-rank PR run.
Re-running replaces matching outputs; preserve results elsewhere before changing
grid resolution. The previous shared example and its existing results are retained
in `examples/2D_peng_robinson_boundary_layer/`.
