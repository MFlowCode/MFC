# 2D Lumped Cylinder Cooling in a Gas Stream

A hot, inert cylinder (an immersed boundary with `thermal_bc = 3`) sits in a cold air stream. The body is one lumped
solid: a single temperature `Twall`, advanced every step by the heat the gas takes from it,

    m c_s dT/dt = Q_in + P - eps sigma A (T^4 - T_rad^4)

with `Q_in` the energy flux through the faces between the body's cells and the gas (see `thermal_bc` in the case
documentation). Here `P` (`heat_power`) and the radiation term are zero and the surface is inert, so the
body only cools. `ib_heat_wrt` writes the heat rate and `Twall` of every step to `D/ib_heat.dat`.

## What it shows

- `Twall` falls from 900 K toward the 300 K stream as a roughly exponential relaxation, with the heat rate largest at the
  start.
- The wall temperature follows the heat the faces report: `heat_balance.py` checks, step by step, that `m c_s dT` equals
  `Q dt` from the written file (with `ib_force_stride = 1`). This verifies the bookkeeping of the update, not that the gas
  gained the same energy and not the heat-transfer coefficient.
- `thermal_layer` starts the gas in the conduction solution for the hot wall, so the cylinder does not meet the cold gas as a
  temperature step at `t = 0`.

## Setup

Air at 1 atm and 300 K flows at 20 m/s (Re_D about 800 at the film temperature) past a cylinder of diameter 2 mm and
initial temperature 900 K. The gas is the reduced 11-species mechanism shared with `2D_ibm_reacting_surface`, with
transport and species diffusion on and reactions off; the flow numerics are the ones of that example. The 20-flow-through
run takes about 11000 steps on a 168 x 120 grid.

The heat capacity is deliberately tiny. A real metal rod of this size cools over seconds, about 10^4 flow-through times;
`case.py` instead sets `cp_solid` so the thermal time is a quarter of the simulated 2 ms, using a Churchill-Bernstein
Nusselt number for the film coefficient, and prints the estimate. The dashed curve in the plot is that exponential. It is an
order-of-magnitude estimate: the actual coefficient is not uniform around the cylinder and changes as the wake develops, so
the computed curve is not expected to lie on it.

## Run

    ./mfc.sh run examples/2D_ibm_lumped_wall/case.py -t pre_process simulation -n 1 -g 0 --gpu acc
    cd examples/2D_ibm_lumped_wall && python3 heat_balance.py

Use `--gpu mp` for an OpenMP-offload build, or drop the GPU flags on a CPU build. The plot is written to
`lumped_wall.png`.

## Limitations

One temperature per body is valid only while the Biot number is small; the geometry and mass are fixed (no recession);
the update is explicit; the body must be a non-moving circle, sphere or cylinder; and the chemistry module is required
even though no reaction occurs. The balance has been checked against first-order estimates, not an analytic solution.
