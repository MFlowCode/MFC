# JWL EOS in MFC

This branch adds Jones-Wilkins-Lee (JWL) equation-of-state support to MFC, including mixed-material JWL explosive products with ideal-gas air.

The implementation follows the ROCFLU-style JWL thermodynamic path for pressure, pressure-to-energy conversion, temperature, and sound speed, then wires that path into MFC primitive/conservative conversion, flux conversion, MPI parameter broadcast, GPU-visible parameter arrays, and the Python case schema.

## Current Capabilities

- Single JWL material in a case through `fluid_pp(i)%eos = 2`.
- Mixed JWL/air cases with `num_fluids = 2`, using the JWL material mass fraction.
- JWL pressure from density, specific internal energy, and JWL mass fraction.
- JWL specific internal energy from pressure, density, and JWL mass fraction.
- JWL temperature estimate from pressure, density, energy, mass fraction, and `cv`.
- JWL and mixed JWL/air sound speed.
- Primitive-to-conservative initialization from pressure.
- Conservative-to-primitive pressure recovery during the simulation.
- Flux-state energy conversion.
- MPI broadcast of all JWL parameters.
- GPU-visible storage/update for all JWL parameter arrays.
- Case-file schema and parameter descriptions for JWL inputs.
- Example 1D and 2D smoke tests.
- 2D high-pressure JWL blast test with a nearby light moving-IBM body initialized from rest.

## Limitations

The current implementation supports one JWL fluid per case. It explicitly rejects:

- more than one JWL fluid,
- MHD with JWL,
- Eulerian bubbles with JWL,
- chemistry with JWL,
- `model_eqns = 4` with JWL.

For mixed-material use, `fluid_pp(i)%jwl_E0` must be positive. JWL `R1`, `R2`, `omega`, `rho0`, `jwl_air_e0`, `jwl_air_rho0`, and `jwl_air_gamma` must also be positive.

## EOS Selector

Use:

```python
"fluid_pp(1)%eos": 2,
```

Supported values are:

```text
1 = stiffened-gas EOS
2 = JWL EOS
```

For the current mixed JWL/air examples:

```python
"model_eqns": 2,
"num_fluids": 2,

"fluid_pp(1)%eos": 2,  # JWL explosive products
"fluid_pp(2)%eos": 1,  # ideal-gas air
```

## JWL Parameters

The JWL material requires:

```python
"fluid_pp(1)%jwl_A": 3.712e11,
"fluid_pp(1)%jwl_B": 3.231e9,
"fluid_pp(1)%jwl_R1": 4.15,
"fluid_pp(1)%jwl_R2": 0.95,
"fluid_pp(1)%jwl_omega": 0.30,
"fluid_pp(1)%jwl_rho0": 1630.0,
```

Mixed JWL/air support also uses:

```python
"fluid_pp(1)%jwl_E0": 1.0089e10,
"fluid_pp(1)%jwl_air_e0": 2.5575e5,
"fluid_pp(1)%jwl_air_rho0": 1.225,
"fluid_pp(1)%jwl_air_gamma": 0.4,
```

Parameter meaning:

```text
jwl_A, jwl_B       JWL pressure constants
jwl_R1, jwl_R2     JWL exponential coefficients
jwl_omega          JWL Gruneisen coefficient
jwl_rho0           JWL reference density
jwl_E0             reference explosive energy for mixed JWL/air states
jwl_air_e0         reference air specific internal energy
jwl_air_rho0       reference air density
jwl_air_gamma      gamma - 1 for the reference air contribution
```

Keep units consistent. For example, if pressure is in Pa, use Pa for `A` and `B`; if density is in kg/m3, use kg/m3 for `rho0`.

## Equation Form

For a pure JWL state, the pressure form is:

```text
V = rho0 / rho

p = A * (1 - omega / (R1 * V)) * exp(-R1 * V)
  + B * (1 - omega / (R2 * V)) * exp(-R2 * V)
  + omega * rho * e
```

where:

```text
p     pressure
rho   mixture or material density
rho0  JWL reference density
V     relative specific volume
e     specific internal energy
A, B  pressure constants
R1,R2 exponential constants
omega JWL Gruneisen coefficient
```

In MFC, total energy and kinetic energy are stored per volume, so:

```text
rho * e = total_energy - kinetic_energy
```

For pressure initialization, the inverse is direct once the cold pressure term is known:

```text
e = (p - pcold) / (omega * rho)
```

For mixed JWL/air states, MFC passes the JWL mass fraction `Y` into the JWL helpers. The implementation uses ROCFLU-style mixture pressure, pressure-to-energy, temperature, and sound-speed helpers:

```text
s_jwl_pressure_er
s_jwl_energy_pr
s_jwl_temperature_pr
s_jwl_mixture_sound_speed_squared
```

## Code Map

Main implementation:

```text
src/common/m_variables_conversion.fpp
```

Key helpers:

```text
s_jwl_pcold
s_jwl_dpcold_drho
s_jwl_sound_speed_squared
s_jwl_mixture_sound_speed_squared
s_jwl_pressure_er
s_jwl_energy_pr
s_jwl_temperature_pr
```

Integration points in `m_variables_conversion.fpp`:

```text
s_compute_pressure
s_convert_conservative_to_primitive_variables
s_convert_primitive_to_conservative_variables
s_convert_primitive_to_flux_variables
s_compute_speed_of_sound
```

JWL material parameters are stored in:

```text
src/common/m_derived_types.fpp
```

as:

```text
fluid_pp(i)%eos
fluid_pp(i)%jwl_A
fluid_pp(i)%jwl_B
fluid_pp(i)%jwl_R1
fluid_pp(i)%jwl_R2
fluid_pp(i)%jwl_omega
fluid_pp(i)%jwl_rho0
fluid_pp(i)%jwl_E0
fluid_pp(i)%jwl_air_e0
fluid_pp(i)%jwl_air_rho0
fluid_pp(i)%jwl_air_gamma
```

Defaults are set in:

```text
src/pre_process/m_global_parameters.fpp
src/simulation/m_global_parameters.fpp
src/post_process/m_global_parameters.fpp
```

MPI broadcasts are handled in:

```text
src/pre_process/m_mpi_proxy.fpp
src/simulation/m_mpi_proxy.fpp
src/post_process/m_mpi_proxy.fpp
```

Case schema and descriptions are registered in:

```text
toolchain/mfc/params/definitions.py
toolchain/mfc/params/descriptions.py
```

## Example Cases

### 1D Mixed JWL/Air Smoke Test

```text
examples/1D_jwl_mixture_test/case.py
examples/1D_jwl_mixture_test/README.md
```

This is a small 1D JWL-rich pressure driver expanding into mostly air.

Current settings:

```text
grid          79x0x0
dt            1.0e-7
steps         20
save period   10
driver p      1.0e7 Pa
```

Run:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
```

### 2D Mixed JWL/Air Blast With Light IBM Body

```text
examples/2D_jwl_mixture_test/case.py
examples/2D_jwl_mixture_test/README.md
```

This is a 2D high-pressure JWL blast near a circular immersed body. The body uses the moving-IBM path, starts with zero translational and angular velocity, and has low mass so pressure forces can accelerate it.

Current settings:

```text
grid              63x31x0
dt                2.5e-8
steps             1000
save period       250
driver p          5.0e9 Pa
driver center     (0.20, 0.25)
driver radius     0.09
IB center         (0.34, 0.25)
IB radius         0.045
IB surface gap    about 0.005
IB initial vel    0
IB angular vel    0
IB mass           0.01
```

Run:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

Expected post-process output times:

```text
0, 250, 500, 750, 1000
```

Expected moving-IBM output includes:

```text
restart_data/ib_state_250.dat
restart_data/ib_state_500.dat
restart_data/ib_state_750.dat
restart_data/ib_state_1000.dat

p_all/p0/250/ib_data.dat
p_all/p0/500/ib_data.dat
p_all/p0/750/ib_data.dat
p_all/p0/1000/ib_data.dat
```

## Build and Verification

Build:

```bash
./mfc.sh build -j 2
```

Run the smoke tests:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

The current branch has been verified with:

```text
./mfc.sh build -j 2
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

The 2D case completed 1000 steps with the stronger `5.0e9` Pa driver and the light `0.01` mass IBM body.

On macOS, `mpirun` may occasionally fail before `pre_process` starts with a dynamic-loader `dyld` error. That failure is outside the case physics; rerunning the same command has completed cleanly.

## Files Changed

```text
README-JWL-EOS.md
src/common/m_derived_types.fpp
src/common/m_variables_conversion.fpp
src/pre_process/m_global_parameters.fpp
src/pre_process/m_mpi_proxy.fpp
src/simulation/m_global_parameters.fpp
src/simulation/m_mpi_proxy.fpp
src/post_process/m_global_parameters.fpp
src/post_process/m_mpi_proxy.fpp
toolchain/mfc/params/definitions.py
toolchain/mfc/params/descriptions.py
examples/1D_jwl_mixture_test/case.py
examples/1D_jwl_mixture_test/README.md
examples/2D_jwl_mixture_test/case.py
examples/2D_jwl_mixture_test/README.md
```
