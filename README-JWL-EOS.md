# JWL EOS Support in MFC

This branch adds Jones-Wilkins-Lee (JWL) equation-of-state support to MFC for explosive-product calculations, with the main tested use case being a mixed JWL-products/air flow.

The implementation was guided by ROCFLU's JWL thermodynamic path, but it is not yet a line-for-line port of every ROCFLU safety layer. The goal of the current work is to give MFC a usable JWL material model, wire it through the usual primitive/conservative conversion paths, and provide simple 1D and 2D tests that exercise the model.

## What Works Now

The current implementation supports one JWL fluid in a case:

```python
"fluid_pp(1)%eos": 2
```

The standard mixed setup is:

```python
"model_eqns": 2
"num_fluids": 2

"fluid_pp(1)%eos": 2  # JWL explosive products
"fluid_pp(2)%eos": 1  # ideal-gas air
```

The JWL path is used in:

- pressure recovery from density and internal energy,
- energy initialization from pressure and density,
- temperature estimation,
- sound-speed calculation,
- primitive-to-conservative conversion,
- conservative-to-primitive conversion,
- primitive-to-flux conversion,
- pre-process, simulation, and post-process parameter handling,
- MPI broadcast of JWL parameters,
- GPU-visible parameter storage,
- Python case parameter definitions and descriptions.

The branch also includes small example cases:

```text
examples/1D_jwl_mixture_test
examples/2D_jwl_mixture_test
```

The 2D example is currently a bottom-wall blast problem with a moving immersed-boundary particle bed. It is meant to test the coupled JWL + moving IBM workflow, not to serve as a calibrated blast experiment.

## JWL Parameters

A JWL fluid needs the usual pressure constants and reference density:

```python
"fluid_pp(1)%jwl_A": 3.712e11,
"fluid_pp(1)%jwl_B": 3.231e9,
"fluid_pp(1)%jwl_R1": 4.15,
"fluid_pp(1)%jwl_R2": 0.95,
"fluid_pp(1)%jwl_omega": 0.30,
"fluid_pp(1)%jwl_rho0": 1630.0,
```

The mixed JWL/air helper also uses:

```python
"fluid_pp(1)%jwl_E0": 1.0089e10,
"fluid_pp(1)%jwl_air_e0": 2.5575e5,
"fluid_pp(1)%jwl_air_rho0": 1.225,
"fluid_pp(1)%jwl_air_gamma": 0.4,
```

The meanings are:

```text
jwl_A, jwl_B       pressure constants
jwl_R1, jwl_R2     exponential constants
jwl_omega          Gruneisen coefficient
jwl_rho0           reference explosive density
jwl_E0             reference explosive specific energy
jwl_air_e0         reference air specific internal energy
jwl_air_rho0       reference air density
jwl_air_gamma      air gamma-minus-one coefficient used by the mixed helper
```

Use consistent units. The examples use SI units: Pa, kg, m, s, and J/kg.

## Exact Equation Being Solved

This implementation does not use one single formula for every mixed cell. It uses a piecewise JWL/air helper in `src/common/m_variables_conversion.fpp`, with the JWL mass fraction

```text
Y = alpha_rho_JWL / rho
rho_safe = max(rho, sgm_eps)
Y_safe = min(max(Y, 0), 1)
```

For cells with almost no JWL material, the code uses an ideal-gas air fallback:

```text
if Y_safe <= 1.0e-2:
    p = gamma_air * rho * e
```

Here `gamma_air` is stored as `jwl_air_gamma`; in the current examples it is `0.4`, meaning gamma minus one for air.

For JWL-rich states, the pressure is:

```text
p = A_eff * (1 - omega * rho / (R1 * rho0)) * exp(-R1 * rho0 / rho)
  + B_eff * (1 - omega * rho / (R2 * rho0)) * exp(-R2 * rho0 / rho)
  + omega * rho * e
```

For pure JWL behavior, the effective constants are just the input constants:

```text
A_eff = A
B_eff = B
omega = omega0
```

For mixed JWL/air cells with:

```text
1.0e-2 < Y_safe <= 0.99
E0 > 0
```

the code modifies the constants before evaluating the pressure. It first defines:

```text
eJ = E0 / rho0
ma = A / max(eJ - air_e0, sgm_eps)
mb = B / max(eJ - air_e0, sgm_eps)
```

It also builds an effective Gruneisen coefficient:

```text
if rho < air_rho0:
    mp = 0
else:
    mp = (omega0 - air_gamma) / max(rho0 - air_rho0, sgm_eps)

omega = max(mp * (rho - air_rho0) + air_gamma, omega0)

if rho > rho0:
    omega = omega0
```

Then it chooses effective JWL pressure constants from the specific internal energy:

```text
if e >= eJ:
    A_eff = A
    B_eff = B
elif e <= air_e0:
    A_eff = 0
    B_eff = 0
else:
    A_eff = A + ma * (e - eJ)
    B_eff = B + mb * (e - eJ)
```

Those `A_eff`, `B_eff`, and `omega` values are inserted into the pressure equation above.

The inverse pressure-to-energy path solves the same model. For the low-JWL air fallback:

```text
if Y_safe <= 1.0e-2:
    e = max(p / (air_gamma * rho), air_e0)
```

For mixed JWL/air cells, the code computes:

```text
C1 = (1 - omega * rho / (R1 * rho0)) * exp(-R1 * rho0 / rho)
C2 = (1 - omega * rho / (R2 * rho0)) * exp(-R2 * rho0 / rho)

e = (p + (ma * eJ - A) * C1 + (mb * eJ - B) * C2)
    / max(ma * C1 + mb * C2 + omega * rho, sgm_eps)
```

If that mixed inversion gives an energy outside the intended interval, or if the state is treated as pure JWL, the code falls back to the direct JWL inversion:

```text
e = p / (omega * rho)
  - A * (1 / (omega * rho) - 1 / (rho0 * R1)) * exp(-R1 * rho0 / rho)
  - B * (1 / (omega * rho) - 1 / (rho0 * R2)) * exp(-R2 * rho0 / rho)
```

Finally, if the recovered energy is non-positive:

```text
e = air_e0
```

The temperature helper uses the same effective constants and computes:

```text
T = (p - A_eff * exp(-R1 * rho0 / rho)
       - B_eff * exp(-R2 * rho0 / rho))
    / max(omega * cv * rho, sgm_eps)
```

If that gives a non-positive temperature in a JWL-containing cell, it currently resets `T` to `270 K`.

The sound-speed path is also specialized. For low-JWL cells, it reduces to:

```text
c^2 = max(air_gamma * (p/rho + p/(rho * air_gamma)), sgm_eps)
```

For JWL/mixed cells, it evaluates the implemented ROCFLU/Stanley-style derivative expression:

```text
c^2 =
  exp(-R1*rho0/rho)
  * [
      A_eff * (R1*rho0/rho^2 - omega/rho
               - omega/(R1*rho0) - rho*mp/(R1*rho0))
      + (ma*p/rho^2) * (1 - omega*rho/(R1*rho0))
    ]

  + exp(-R2*rho0/rho)
  * [
      B_eff * (R2*rho0/rho^2 - omega/rho
               - omega/(R2*rho0) - rho*mp/(R2*rho0))
      + (mb*p/rho^2) * (1 - omega*rho/(R2*rho0))
    ]

  + omega * (e + p/rho)
  + mp * rho * e
```

and then floors:

```text
c^2 = max(c^2, sgm_eps)
```

So the exact behavior is: ideal-gas air fallback at tiny JWL mass fraction, effective mixed JWL constants in partially mixed cells, and the standard JWL form in JWL-rich cells.

## Code Map

Most of the thermodynamic work lives in:

```text
src/common/m_variables_conversion.fpp
```

Important helper routines include:

```text
s_jwl_pcold
s_jwl_dpcold_drho
s_jwl_pressure_er
s_jwl_energy_pr
s_jwl_temperature_pr
s_jwl_sound_speed_squared
s_jwl_mixture_sound_speed_squared
```

Those helpers are called from the usual conversion and sound-speed paths:

```text
s_compute_pressure
s_convert_conservative_to_primitive_variables
s_convert_primitive_to_conservative_variables
s_convert_primitive_to_flux_variables
s_compute_speed_of_sound
```

JWL parameters are stored in the common physical-parameter type:

```text
src/common/m_derived_types.fpp
```

Defaults are initialized in:

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

Case-file schema support is in:

```text
toolchain/mfc/params/definitions.py
toolchain/mfc/params/descriptions.py
```

## Current Restrictions

The current implementation intentionally rejects combinations that have not been made safe yet:

- more than one JWL fluid in the same case,
- JWL with MHD,
- JWL with Eulerian bubbles,
- JWL with chemistry,
- JWL with `model_eqns = 4`.

Required JWL inputs must be positive where physically necessary. In particular:

```text
jwl_E0
jwl_R1
jwl_R2
jwl_omega
jwl_rho0
jwl_air_e0
jwl_air_rho0
jwl_air_gamma
```

## Relationship To ROCFLU

This branch brings over the main ideas from ROCFLU's JWL thermodynamic workflow:

- pressure from density and energy,
- energy from pressure and density,
- temperature support,
- sound-speed support,
- mixed JWL/air behavior through a JWL mass fraction,
- basic validation of unsupported model combinations.

There are still ROCFLU features that are not fully reproduced:

- a standalone `RFLU_JWL_Yfix()`-style mass-fraction correction pass,
- the full production `ComputePressureMixt()` and `ComputeEnergyMixt()` guard/recovery workflow,
- robust NaN recovery for every bad-state path,
- a dedicated JWL energy-floor stability layer,
- ROCFLU PICL particle-in-cell coupling,
- multiple JWL materials in the same case.

So the current answer is: MFC now has a useful JWL EOS path, but it does not yet do everything ROCFLU's mature production JWL implementation does.

## Examples

### 1D Smoke Test

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
```

This is the fastest test. It checks that a JWL-rich pressure driver can expand into air without breaking pressure, energy, or sound-speed conversion.

### 2D Blast With Moving Particle Bed

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

This is the more useful integrated test. It uses a compact bottom-wall JWL blast and a 40-particle moving IBM bed. It checks the coupled path: JWL EOS, shock propagation, IBM marker generation, particle motion, collision parameters, and IBM state output.

The 2D case is intentionally idealized. The particles are 2D circular bodies, which represent cylinder cross-sections. They are not true 3D spheres, and the blast is initialized as a high-pressure JWL products patch rather than a resolved reactive detonation.

## Practical Notes

Use a tiny volume-fraction floor in mixed JWL/air cases:

```python
eps = 1.0e-8
```

This is not a physical interface thickness. It keeps both materials present at a tiny level and avoids exactly zero volume fractions, which can otherwise cause singular mixture states.

For open boundaries such as:

```python
"bc_x%beg": -3
"bc_x%end": -3
"bc_y%end": -3
```

remember that these are fluid boundary conditions. They do not remove IBM particles and they do not stop the simulation when a particle leaves the domain. Particle removal or escape-based stopping would need separate IBM logic.
