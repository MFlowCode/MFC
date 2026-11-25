# Three isothermal laminar flat plates

Standalone versions with independent `case.py` files are now available in
[`2D_boundary_layer_ideal_gas`](../2D_boundary_layer_ideal_gas/README.md),
[`2D_boundary_layer_pr_normal`](../2D_boundary_layer_pr_normal/README.md), and
[`2D_boundary_layer_pr_supercritical`](../2D_boundary_layer_pr_supercritical/README.md).
Use those folders to edit, rerun, and animate each case independently. This shared
example and its existing results are retained for reference.

CPU/serial, single-component nitrogen with a fixed freestream inlet, an isothermal
no-slip lower wall, a fixed freestream upper boundary, and an extrapolated outlet. Uniform initial
freestream develops a boundary layer from the plate's leading edge. All units are SI.

| Directory | Equation of state | Freestream temperature | Pressure | Wall temperature |
| --- | --- | --- | --- | --- |
| `normal_ideal` | Ideal gas | 300 K | 101325 Pa | 270 K |
| `normal_pr` | Peng–Robinson | 300 K | 101325 Pa | 270 K |
| `supercritical_pr` | Peng–Robinson | 200 K | 7.5 MPa | 180 K |

The PR critical constants are 126.192 K and 3.3958 MPa. Both the freestream and the
cooled wall of the third case are above the critical temperature; verification
also checks that every saved final cell remains above the critical temperature
and pressure.

## Run on macOS

From the MFC repository root, using the existing MFC environment:

```bash
build/venv/bin/python examples/2D_peng_robinson_boundary_layer/run_all.py
```

Or run the three cases separately, followed by the comparison:

```bash
./mfc.sh run examples/2D_peng_robinson_boundary_layer/normal_ideal/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
./mfc.sh run examples/2D_peng_robinson_boundary_layer/normal_pr/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
./mfc.sh run examples/2D_peng_robinson_boundary_layer/supercritical_pr/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
build/venv/bin/python examples/2D_peng_robinson_boundary_layer/verify.py
```

A working Fortran compiler is needed; MPI and GPU hardware are not needed. The
first build can take substantially longer than subsequent runs. Each run uses
96 × 64 cells and approximately 21760 RK3 steps (four flow-through times), saving
every half flow-through time. The verifier reads double-precision serial restart
files directly; no Silo viewer or `post_process` run is required. Re-running a case
replaces that case's existing output.

## What is held fixed

All three use Mach 0.2, plate Reynolds number 10000, physical height/length ≈ 0.487, wall/free-stream
temperature ratio 0.9, and a freestream Prandtl number of 0.72. Dynamic viscosity is
prescribed as a constant 1.8e-5 Pa s. Conductivity is constant within each run and
set to `mu * cp_freestream / 0.72`. Ideal-gas cv is `2.5 R`; PR includes its residual
energy and heat capacities. These are controlled **constant-transport benchmarks**,
not a state-dependent nitrogen transport model. An EOS alone does not specify
viscosity or conductivity.

Plate length adjusts to match Reynolds number: approximately 2.240 mm for ideal
gas, 2.238 mm for normal PR, and 20.64 micrometers for supercritical PR. Thus the
comparison is at matched dimensionless flow conditions, not on a common physical
plate. The small dense-fluid plate keeps its boundary layer resolved at the same
cell count. Wall-normal stretching preserves fine near-wall spacing while moving
the upper boundary farther away. The input y-domain height is 0.2 L; MFC’s
stretching maps it to the physical height of approximately 0.487 L. Edit `setup_case.py` to change the common grid or physical choices.

## Verification and output

`verify.py` fails on nonfinite/nonphysical states, appreciable final profile drift,
missing velocity or thermal layers, lack of downstream growth of the 99%-velocity thickness, inconsistent
wall extrapolation, or excessive outer-flow disturbance. It checks that the normal
ideal/PR dimensionless profiles nearly coincide. Profiles at x/L ≈ 0.75 are plotted
in `viz/profiles.png`; `viz/fields.png` shows the two-dimensional velocity and
temperature layers near the wall. Numerical profiles and metrics are in the same directory as
CSV and JSON files.

The dilute cases also have a deliberately loose comparison against the classical
[Blasius flat-plate reference](https://www.grc.nasa.gov/www/wind/valid/fplam/fplam.html).
Blasius assumes incompressible constant-property flow: it is an approximate reference
for a cooled compressible gas, not an exact solution for the supercritical case.
The supercritical curve is displayed against it without an accuracy assertion.
These are reproducible boundary-layer smoke tests, not a grid-convergence study or
experimental validation of supercritical heat transfer. Refine NX and NY and repeat
before using wall shear or heat flux quantitatively.

The existing conduction timestep uses ideal cv. For the supported PR acentric
factor interval, the residual cv is positive, so this is a conservative diffusion
bound. The existing isothermal boundary reflects temperature about the specified
wall value for the conduction flux; density and pressure retain their standard
zero-gradient wall treatment.

## Reference CPU run

All three cases passed the serial macOS check. Maximum changes between the final
two saves, normalized by freestream velocity and by `T∞ − Tw`, respectively:

| Case | Velocity change | Temperature change |
| --- | --- | --- |
| Normal ideal gas | 0.363% | 0.0681% |
| Normal PR | 0.359% | 0.0674% |
| Supercritical PR | 0.257% | 0.0545% |

The largest normal ideal/PR profile difference was 2.36e-5 in normalized units
(0.00236%). The dilute-case velocity RMS difference from Blasius was about 6.5%;
this coarse finite-domain, cooled-wall check should not be read as an exact
Blasius validation. Serial runtimes after compilation were approximately 5–7
minutes per case on the test Mac. The first executable builds took longer.

To verify just one existing result, add, for example, `--case normal_pr` to the
verification command. Running verification without that option regenerates the
complete three-case comparison.
