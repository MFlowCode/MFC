# Static Drop (2D, all-Mach pressure projection, surface tension)

A quarter water drop (R = 0.15, sigma = 8) at rest in air, as in
`examples/2D_laplace_pressure_jump` but with the five-equation model. The exact state is
rest with a pressure jump sigma/R = 53.33; any velocity is a spurious current.

```shell
./mfc.sh run examples/2D_laplace_projection/case.py
python3 examples/2D_laplace_projection/analyze.py examples/2D_laplace_projection
```

Surface tension enters through MFC's capillary stress tensor, explicitly. At N = 50 and
t = 0.5:

|                              | dt     | steps   | pressure jump | max spurious speed |
|------------------------------|--------|---------|---------------|--------------------|
| projection                   | 1e-3   | 500     | 52.9          | ~1 m/s             |
| explicit HLLC (`--explicit`) | 5e-6   | 100,000 | 53.3          | 0.02 m/s           |

The pressure jump is right, but the spurious currents are about 50x the explicit solver's
and do not shrink with dt: the capillary force reaches the cells as a stress divergence,
while the pressure that should cancel it acts on faces, so the two do not balance exactly.
