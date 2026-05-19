# JWL EOS in MFC

This branch adds a single-fluid Jones-Wilkins-Lee (JWL) equation of state option.

## Select the EOS

Set the EOS selector for fluid 1 in a case file:

```python
"fluid_pp(1)%eos": 2,
```

Values are:

```text
1 = existing stiffened-gas EOS
2 = JWL EOS
```

The current JWL implementation supports single-fluid, non-chemistry, non-MHD, non-bubble cases. Use `model_eqns = 2` and `num_fluids = 1`.


## Equation Form

The JWL EOS implemented here uses the relative specific volume

```text
V = rho0 / rho
```

and computes pressure as

```text
p = A * (1 - omega / (R1 * V)) * exp(-R1 * V)
  + B * (1 - omega / (R2 * V)) * exp(-R2 * V)
  + omega * rho * e
```

where:

```text
p     = pressure
rho   = current density
rho0  = reference density
V     = relative specific volume
e     = specific internal energy
A, B  = pressure constants
R1,R2 = exponential constants
omega = JWL Gruneisen coefficient
```

In the MFC implementation, total energy and kinetic energy are stored per volume, so

```text
rho * e = total_energy - kinetic_energy
```

The code therefore evaluates the pressure as

```text
pcold = A * (1 - omega / (R1 * V)) * exp(-R1 * V)
      + B * (1 - omega / (R2 * V)) * exp(-R2 * V)

p = pcold + omega * (total_energy - kinetic_energy)
```

For primitive-to-conservative initialization from pressure, the inverse is direct:

```text
total_energy = kinetic_energy + (p - pcold) / omega
```

The JWL sound-speed branch uses

```text
c^2 = dpcold/drho + ((1 + omega) * p - pcold) / rho
```

## Required JWL Parameters

Add the six JWL constants to the case file:

```python
"fluid_pp(1)%jwl_A": 3.712e11,
"fluid_pp(1)%jwl_B": 3.231e9,
"fluid_pp(1)%jwl_R1": 4.15,
"fluid_pp(1)%jwl_R2": 0.95,
"fluid_pp(1)%jwl_omega": 0.30,
"fluid_pp(1)%jwl_rho0": 1630.0,
```

Keep units consistent: `A`, `B`, and pressure must use the same pressure unit; `rho0` and density must use the same density unit.

## Minimal Case Settings

A minimal JWL setup should include:

```python
"model_eqns": 2,
"num_fluids": 1,
"mhd": "F",
"bubbles_euler": "F",

"patch_icpp(1)%geometry": 1,
"patch_icpp(1)%x_centroid": 0.0005,
"patch_icpp(1)%length_x": 0.001,
"patch_icpp(1)%vel(1)": 0.0,
"patch_icpp(1)%pres": 2.0e10,
"patch_icpp(1)%alpha_rho(1)": 1630.0,
"patch_icpp(1)%alpha(1)": 1.0,

"fluid_pp(1)%eos": 2,
"fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
"fluid_pp(1)%pi_inf": 0.0,
"fluid_pp(1)%jwl_A": 3.712e11,
"fluid_pp(1)%jwl_B": 3.231e9,
"fluid_pp(1)%jwl_R1": 4.15,
"fluid_pp(1)%jwl_R2": 0.95,
"fluid_pp(1)%jwl_omega": 0.30,
"fluid_pp(1)%jwl_rho0": 1630.0,
```

`gamma` and `pi_inf` are still read by existing MFC setup paths, but the JWL pressure/energy branch uses the JWL constants.

Make sure the initial pressure is high enough that `(p - pcold) / omega` is positive. Otherwise the initialized internal energy is negative.

## Build and Run

From the MFC repository root:

```bash
./mfc.sh build --no-mpi -j 4
./mfc.sh run path/to/case.py --no-mpi --no-build
```

To write useful fields, enable output such as:

```python
"prim_vars_wrt": "T",
"rho_wrt": "T",
"E_wrt": "T",
"pres_wrt": "T",
```

## Where the EOS Is Implemented

The JWL pressure helper is in `src/common/m_variables_conversion.fpp`:

```fortran
subroutine s_jwl_pcold(rho, A, B, R1, R2, omega, rho0, pcold)
```

The full pressure uses:

```text
p = pcold + omega * (total_energy - kinetic_energy)
```

Since MFC stores total and kinetic energy per volume, `(total_energy - kinetic_energy)` is `rho * e`.


## Modified Code Map

The implementation touches these files:

```text
src/common/m_derived_types.fpp
```

Adds the EOS selector and JWL material constants to `type physical_parameters`:

```text
fluid_pp(i)%eos
fluid_pp(i)%jwl_A
fluid_pp(i)%jwl_B
fluid_pp(i)%jwl_R1
fluid_pp(i)%jwl_R2
fluid_pp(i)%jwl_omega
fluid_pp(i)%jwl_rho0
```

```text
src/common/m_variables_conversion.fpp
```

Adds the JWL thermodynamic implementation:

```text
s_jwl_pcold              computes the JWL cold pressure terms
s_jwl_dpcold_drho       computes d(pcold)/d(rho) for sound speed
s_compute_pressure      selects JWL when fluid_pp(1)%eos = 2
s_convert_primitive_to_conservative_variables
                         computes total energy from JWL pressure
s_convert_primitive_to_flux_variables
                         computes JWL total energy for fluxes
s_compute_speed_of_sound
                         uses the JWL sound-speed expression
```

It also allocates GPU-visible arrays for the EOS selector and JWL constants:

```text
eos_idxs, jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s
```

```text
src/pre_process/m_global_parameters.fpp
src/simulation/m_global_parameters.fpp
src/post_process/m_global_parameters.fpp
```

Set default values for the new fields. The default remains stiffened gas:

```text
fluid_pp(i)%eos = 1
```

```text
src/pre_process/m_mpi_proxy.fpp
src/simulation/m_mpi_proxy.fpp
src/post_process/m_mpi_proxy.fpp
```

Broadcast the new EOS selector and JWL constants for MPI runs.

```text
toolchain/mfc/params/definitions.py
toolchain/mfc/params/descriptions.py
```

Register the new case-file parameters so cases can use:

```python
"fluid_pp(1)%eos": 2,
"fluid_pp(1)%jwl_A": ...,
"fluid_pp(1)%jwl_B": ...,
"fluid_pp(1)%jwl_R1": ...,
"fluid_pp(1)%jwl_R2": ...,
"fluid_pp(1)%jwl_omega": ...,
"fluid_pp(1)%jwl_rho0": ...,
```

