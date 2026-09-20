# Probe files across a re-run

`s_open_probe_files` appends whenever `D/probe*_prim.dat` already exists:

```fortran
if (file_exist) then
    open (..., STATUS='old', POSITION='append')
```

That is correct when a run is being **continued**. It is wrong when one is being **started over**, which is
what happens every time a case is re-run in place after a parameter change. The second run's rows land on top
of the first's, nothing in the file marks the join, and the time column simply resets partway down. A reader
sees one monotonic series and is silently wrong.

There is no warning, no header and no separator. The two runs need not even share a grid — a case whose
resolution changed between runs produces a file whose first half was recorded at different probe locations.

## Reproducing

```
./mfc.sh run examples/2D_probe_rerun/case.py -n 1      # 20 steps
./mfc.sh run examples/2D_probe_rerun/case.py -n 1      # same again, from scratch
wc -l D/probe1_prim.dat
```

| | rows in `D/probe1_prim.dat` |
| --- | --- |
| before | **40** — two runs of 20, spliced |
| after | **20** |

And the time column, read straight through, before the fix:

```
0.028977
0.030682
0.032386      <- end of run 1
0.000000      <- run 2 starts, time goes backwards
0.001705
0.003409
```

## The fix

Append only when continuing: `t_step_start > 0`, or `n_start > 0` under `cfl_dt`. A fresh start replaces the
file, which is what every other output MFC writes already does. `s_open_com_files` had the same pattern and
gets the same treatment.

## Why it matters beyond tidiness

This produced three separate wrong numbers in one project before it was noticed. The worst was a jet whose
probe files held a t = 15 run of 4,962 rows followed by a t = 40 run of 13,233 — **on different grids**.
Read together they manufactured a velocity drop of 0.99 U_j in a single sample, which was time running
backwards at the seam and was diagnosed as a physical instability first.

A related symptom is louder and easier to spot: if the probe output format changes between runs, the column
count changes partway down the file and `numpy.loadtxt` refuses it outright ("the number of columns changed
from 11 to 18"). That one at least announces itself. The time reset does not.
