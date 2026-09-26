# Dam Break Against an Obstacle (2D, all-Mach pressure projection, immersed boundary)

A water column (a = 0.146 m wide, 2a high) collapses in a 4a x 4a air-filled tank and hits a
no-slip block 0.024 m wide and 0.048 m high centred on the floor at x = 2a, after the experiment of
Koshizuka, Tamako & Oka (Computational Fluid Dynamics Journal 4:29-46, 1995). The surge is thrown up over the
block into a jet that reaches the far side and the ceiling. The block is a stationary immersed
boundary: the ghost-cell method sets its boundary conditions, and the projection closes the faces
it touches. Water and air carry viscosity and surface tension, with MTHINC interface compression
(ic_beta = 0.6), and the run ends at t = 0.5 s with adaptive time stepping.

```shell
./mfc.sh run examples/2D_dam_break_obstacle/case.py -- --ppa 40
```

At a/40 (160 x 160) the run takes about 6,200 steps, about 1.5 minutes on one A100.
