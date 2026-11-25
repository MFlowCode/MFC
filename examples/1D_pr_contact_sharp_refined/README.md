# Xenon contact advection — sharp refined

From the repository root on macOS, run

```bash
./mfc.sh run examples/1D_pr_contact_sharp_refined/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
build/venv/bin/python examples/1D_pr_contact_sharp_refined/verify.py
```

The profile translates through one periodic domain at 20 m/s along x. Its exact pressure is 10 MPa everywhere throughout. This folder uses 256 cells, interface width parameter 0.005 L, and dense-stream temperature 330 K. Edit these at the top of case.py. The two plateau densities have ratio 2.

The verifier reads 101 serial primitive output frames and writes pressure_history.csv. A nonzero exit after printing peak errors indicates loss of pressure equilibrium above the selected relative tolerance (default 1e-8). A solver abort or missing frame is a separate outcome. Compare peak errors at fixed width across grids and at fixed grid across widths. Compare the 305 K case with the sharp 330 K case. Use fresh output when changing the case so older frames cannot mask an interrupted run.
