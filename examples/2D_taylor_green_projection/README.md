# Taylor-Green Vortex (2D, all-Mach pressure projection, viscous)

Decaying Taylor-Green vortex in water at Re = 100, Mach ~1e-3, on [-pi, pi]^2 with periodic
boundaries. The exact solution decays as exp(-2 nu t), its kinetic energy as exp(-4 nu t).
The projection steps at ~400x the water acoustic limit; viscosity is explicit.

```shell
./mfc.sh run examples/2D_taylor_green_projection/case.py -- --N 64
python3 examples/2D_taylor_green_projection/analyze.py examples/2D_taylor_green_projection
```

At t = 10 (kinetic energy decayed to 0.67):

| N  | KE relative error | velocity L2 error |
|----|-------------------|-------------------|
| 32 | 6.8e-3            | 3.6e-3            |
| 64 | 7.3e-4            | 3.9e-4            |
