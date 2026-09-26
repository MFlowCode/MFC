# Interfacial Breakup (2D water-air multimode Rayleigh-Taylor / Faraday)

A single-valued water-air interface, seeded with a multimode random perturbation
(`perturbation2d.py`), is driven by a background plus oscillatory acceleration
`g_y + k_y sin(w_y t)` in a periodic-by-walled box, with viscosity and surface tension.
`case.py` writes `interface_profile.dat` next to itself; hardcoded IC 209 reads it and blends
the two fluids across `y = y_int(x)`. See the header of `case.py` for the scaling (We sets the
length scale; Re follows).

```shell
./mfc.sh run examples/2D_interface_breakup/case.py -n 4             # all-Mach pressure projection
./mfc.sh run examples/2D_interface_breakup/case.py -n 4 -- --explicit
```

The projection solves the acoustics implicitly, so water keeps its physical stiffness and the step
is set by the flow and capillarity. `--explicit` runs the HLLC solver at the acoustic limit with the
original case's settings (water softened to the gas sound speed, `mpp_lim`, `low_Mach = 2`). On the
1024 x 1024 grid the projection takes about 200 steps per 0.2 forcing periods, the explicit solver
about 20,000.
