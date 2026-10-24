# Taylor-Green Vortex at Low Mach Number (3D, all-Mach pressure projection)

The Re = 1600 Taylor-Green vortex of `examples/3D_TaylorGreenVortex` (hardcoded IC 380) at Mach 0.01, set through a stiffened-gas
equation of state with sound speed U0/M. Both solvers take a fixed step: the projection's is set by the flow (acoustic CFL about
50), the explicit solver's by the sound speed.

```shell
./mfc.sh run examples/3D_TaylorGreenVortex_projection/case.py --case-optimization                   # 204 steps to 10 tC
./mfc.sh run examples/3D_TaylorGreenVortex_projection/case.py --case-optimization -- --explicit     # 20,576 steps
```

`--explicit` runs the HLLC solver with the `low_Mach = 2` correction, without which its numerical dissipation grows as the Mach
number falls.
