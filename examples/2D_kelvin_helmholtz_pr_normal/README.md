# Kelvin–Helmholtz: pr normal

Standalone CPU case with xenon thermodynamics, periodic boundaries, a 2:1 density ratio,
and counterflow. All setup and thermodynamic functions are in this folder's `case.py`.
Edit `NX`, `NY`, `FLOW_TIMES`, and `SAVES` there. The default is 128 × 128 cells,
two flow times, and 101 frames. This is an inviscid EOS-sensitivity demonstration,
not a physical transport or grid-converged turbulence calculation.

From the repository root on a Mac with a Fortran compiler and MPI:

```bash
./mfc.sh run examples/2D_kelvin_helmholtz_pr_normal/case.py -t pre_process simulation --mpi --no-gpu -n 4 -j 4
build/venv/bin/python examples/2D_kelvin_helmholtz_pr_normal/verify.py --animate
./mfc.sh run examples/2D_kelvin_helmholtz_pr_normal/case.py -t post_process --mpi --no-gpu -n 4 -j 4
```

The verifier writes `viz/summary.png`, `metrics.json`, `history.csv`, `fields.npz`,
and optionally `density.gif`. Post-processing writes Silo files for VisIt animations.
Parallel I/O is enabled. For serial execution, set `PARALLEL_IO = "F"` and use
`--no-mpi` instead of `--mpi -n 4`. Do not mix serial and MPI results within a run.

A smaller serial trial can be run without editing the defaults:

```bash
MFC_KH_NX=64 MFC_KH_PARALLEL_IO=F ./mfc.sh run examples/2D_kelvin_helmholtz_pr_normal/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
build/venv/bin/python examples/2D_kelvin_helmholtz_pr_normal/verify.py
```

Verification uses the saved grid and timestep, so these grid/duration overrides do not
need to be repeated when verifying. Rerun after changing physical settings. Clear old
output before changing the number of saves; the verifier uses the current input schedule.
The first compilation can be much slower than subsequent runs.

See [comparison design](../2D_kelvin_helmholtz_comparison/README.md) for controls,
physical states, differences from the original example, and the comparison command.
