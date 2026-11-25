# Peng–Robinson acoustic wave on a Mac

From the repository root, run

```sh
./mfc.sh run examples/1D_peng_robinson/case.py -t pre_process simulation --no-mpi --no-gpu -j 4
build/venv/bin/python examples/1D_peng_robinson/verify.py
```

The case uses 64 cells and 400 RK3 steps to propagate a periodic, small-amplitude acoustic wave through one wavelength. The base state represents nitrogen at 200 K and 150 kg/m³, approximately 7.52 MPa. Its PR sound speed is approximately 291.56 m/s. The initial density, velocity, and pressure perturbations form a right-going linear acoustic mode. The verifier compares the final wave to the initial wave and checks mass, momentum, and total-energy conservation. Small nonlinear corrections are expected at the relative density amplitude of 1e-5.

The initial implementation supports CPU, single-component, single-phase flow with the five-equation model, HLLC direct wave speeds, and periodic, reflective, extrapolation, no-slip wall, or fixed Dirichlet boundaries. The case uses fifth-order primitive WENO. GPU execution, phase equilibrium, reacting mixtures, and state-dependent real-fluid transport require additional implementation.

Select `fluid_pp(1)%eos = "peng_robinson"` and supply `pr_tc` in K, `pr_pc` in Pa, `pr_r` in J/(kg K), `pr_omega`, and `cv` in J/(kg K). Here `cv` is the constant **ideal-gas** heat capacity. The actual heat capacities include PR departure terms. This version uses the conventional quadratic acentric-factor correlation for values from 0 to 0.49. `qv`, when supplied, is an additive specific-energy reference. The ideal-gas energy convention is `e_ideal = cv*T + qv`.

The pressure, departure energy, heat capacities, and sound speed share one thermodynamic closure. Cell temperatures are recovered from density and internal energy. Face temperatures are recovered from reconstructed density and pressure before evaluating face energy and sound speed. Both inversions use the positive quadratic root in `sqrt(T)` for this constant-heat-capacity model. Invalid density, temperature, mechanical stability, pressure, or sound speed terminates the run with a PR diagnostic. Local mechanical stability checks do not perform a vapor–liquid equilibrium calculation.

The production Fortran thermodynamics have tests for temperature round trips, finite-difference derivative identities, invalid states, and comparison with Cantera's independent PR implementation. Cantera uses more precise critical constants, so that comparison allows their small coefficient differences. Run these checks with

```sh
build/venv/bin/python -m pytest toolchain/mfc/test_peng_robinson.py --basetemp=build/pr_validation/pytest -q
```

The Cantera model is documented in its [PengRobinson reference](https://www.cantera.org/stable/cxx/d3/ddc/classCantera_1_1PengRobinson.html).

Isothermal viscous cases with prescribed constant transport are provided in
[`../2D_peng_robinson_boundary_layer`](../2D_peng_robinson_boundary_layer/README.md).
