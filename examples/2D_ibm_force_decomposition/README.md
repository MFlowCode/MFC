# Immersed-boundary force under a changing MPI decomposition

A static cylinder in uniform flow at Re 40, Ma 0.1, on a uniform 200 x 200 grid. The body sits at the
origin, which is also the center of the domain, so any decomposition with an even number of ranks in a
direction puts a subdomain edge straight through it.

The force on a rigid body is a property of the flow. Running the same case on a different number of ranks
must not change it. This case measures whether it does.

## Running it

```
./mfc.sh run examples/2D_ibm_force_decomposition/case.py -n 1
./mfc.sh run examples/2D_ibm_force_decomposition/case.py -n 4
```

Each run writes `restart_data/ib_state_200.dat`: 20 doubles, `[time, Fx, Fy, Fz, Tx, Ty, Tz, ...]`.

```
python3 -c "import numpy as np; a = np.fromfile('restart_data/ib_state_200.dat'); print(a[1], a[2])"
```

## What it shows

`fd_order = 4` here, so the force integral's stencil reaches two cells, and a body cell two cells from a
subdomain edge asks for finite-difference coefficients outside the interior. Before the coefficients were
defined over that range they were read off the end of the array:

| ranks | Fx before | Fx after |
| --- | --- | --- |
| 1 | 1.02352019 | 1.02352019 |
| 4 | 1.02933438 | 1.02351629 |

0.57 percent of drag appearing out of adjacent memory, against 0.0004 percent after. The single-rank
answer is unchanged, because a single rank never reaches outside its own interior.

The size of the discrepancy is set by whatever is adjacent in memory and is not bounded by anything: on a
3D sphere at Re 100 on 64 ranks the same read produced a transverse force of 1.08 times the drag on a body
that has none.

Lift is a second, independent check: the cylinder is symmetric about y = 0, so Fy must be zero. It is
about 8e-6 here in every configuration, which is the level at which the staircase representation of a
circle on this grid is symmetric, and it does not move.
