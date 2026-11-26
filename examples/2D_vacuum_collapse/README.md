# Near-vacuum collapse

A disc of radius 0.3 centred at (0.48, 0.53) in the unit square holds gas at 1/1000 of the ambient density and pressure, so both
regions have the same temperature. The ambient gas (\f$\rho = p = 1\f$, \f$\gamma = 1.4\f$, at rest) rushes into the disc,
converges on its centre at t ≈ 0.07 and rebounds as a strong shock. This is the 2D analogue of the blast-centre collapse of
[#2469](https://github.com/MFlowCode/MFC/issues/2469): unlimited WENO5 sets a negative face pressure at the collapse and the
sound speed becomes NaN. `face_bound = 2` (the default here) scales each cell's face states toward the cell state only as far as
positivity needs, which carries the run through.

- `--fluids 2`: the disc holds a second gas (\f$\gamma = 5/3\f$); `--model 6eq` runs it with the 6-equation model.
- `--face-bound 0` reproduces the failure; `--face-bound 1` is the primitive bound (5-equation model only).

Mapped WENO5, HLLC, SSP-RK3, fixed \f$\Delta t\f$ at CFL 0.25 on \|u\| + c ≤ 8, extrapolation boundaries, `-N` = 128 cells per
direction to `--t-end` = 0.12.

`check_collapse.py` reads every save (`parallel_io = T`) and prints the minimum density and pressure and the drift of total mass
and energy. It passes when both minima are positive at every save and the drift stays below 1e-4. The limiter itself is
conservative; the drift is the weak outer rarefaction leaving through the boundaries.

```
./mfc.sh run examples/2D_vacuum_collapse/case.py -t pre_process simulation -n 8 -- --fluids 2 --model 6eq
python3 examples/2D_vacuum_collapse/check_collapse.py examples/2D_vacuum_collapse --fluids 2 --model 6eq
```

Results and limits: @ref validation_positivity.
