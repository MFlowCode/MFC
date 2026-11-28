# 2D IB Accuracy Limits: Wall-Temperature Order Study

The order study behind the wall-temperature rows of @ref ib_accuracy_limits. `case.py` is the `2D_ibm_conduction_annulus` deck
(steady conduction between an isothermal circle, `T1 = 2`, and an isothermal STL ring, `T2 = 3`, with the exact solution
`T(r) = T1 + (T2 - T1) ln(r / R) / ln(r2 / R)`), run from this directory. `order_study.py` runs it at 10, 20 and 40 cells per inner
radius with inverse-distance (`ib_ip_interp = 1`) and TSC (`ib_ip_interp = 2`) image points, measures each run with
`../2D_ibm_conduction_annulus/check_annulus.py` (the wall value of a quadratic fit through the first gas cells, minus the wall
temperature, relative to `T2 - T1`), prints the observed orders and exits 1 if one differs from the page by more than 0.1.

## Results

Magnitude of the wall error relative to `T2 - T1`, 1 rank, 8 conduction times, on an MI300A (CCE 19.0.0, OpenMP offload), `master` `ef766cc1` plus the change that added this example:

| cells / R | idw inner | idw outer | TSC inner | TSC outer |
|---:|---:|---:|---:|---:|
| 10 | 1.83e-2 | 1.18e-2 | 3.68e-3 | 1.27e-4 |
| 20 | 1.18e-2 | 7.71e-3 | 7.67e-4 | 6.10e-5 |
| 40 | 5.00e-3 | 3.42e-3 | 1.87e-4 | 2.97e-5 |
| observed order 10-20, 20-40 | 0.64, 1.24 | 0.61, 1.17 | 2.26, 2.04 | 1.05, 1.04 |
| page | 0.6, 1.2 | 0.6, 1.2 | 2.3, 2.0 | 1.1, 1.0 |

All four match the page within 0.05, and the errors agree with the `2D_ibm_conduction_annulus` README table to within 1 %.
The finest grid (184 x 184 cells, 8648 steps) takes about 3 minutes on one MI300A; the whole study takes about 9.

## Run

    python3 examples/2D_ibm_accuracy_limits/order_study.py -- -n 1 --gpu mp -c tuo   # options for ./mfc.sh run after --

Without options after `--` each run uses `-n 1`. A single run: `./mfc.sh run examples/2D_ibm_accuracy_limits/case.py -n 1 -- --cells_per_R 20 --interp 2`.
The study is not a registered test (it is in `casesToSkip`): the suite registers the same deck as
`2D -> IBM -> Conduction Annulus -> Inverse Distance` and `-> TSC` at 12 cells/R.

## Limits

The ones of `2D_ibm_conduction_annulus`: constant conductivity, a quiescent inviscid gas, a prescribed wall temperature. The
wall-normal gradient and the rotating-wall shear rows of the page are not reproduced here (see the page and
`examples/2D_ibm_taylor_couette`).
