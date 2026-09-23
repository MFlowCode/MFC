@page thermochemistry Thermochemistry implementation

# Mechanism generation

MFC owns the Fortran thermochemistry generator in `toolchain/mfc/thermochem/`.
Cantera loads a mechanism, and `generate_fortran` produces the Fypp source
`m_thermochem.fpp` in the target's build staging directory. The existing CMake build compiles that
module into MFC. Neither Cantera nor Python is called inside the flow solver.
Pyrometheus and JAX are not installation or runtime requirements.

The implementation was adapted from the MIT-licensed Pyrometheus 1.1.1 Fortran
generator. Its license is retained in `toolchain/mfc/thermochem/LICENSE`.
Changes are made directly in MFC; no external fork or upstream generator release
is required. Pymbolic constructs scalar expressions during generation and Mako
renders the Fortran template. Neither is a runtime solver dependency.

| File | Responsibility |
|---|---|
| `thermochem/expressions.py` | NASA7, reaction-rate, equilibrium and transport expressions |
| `thermochem/fortran.py` | Supported-feature checks, expression formatting and module generation |
| `thermochem/module.fpp.mako` | Fortran interface and numerical routines |
| `thermochem/fingerprint.py` | Mechanism and generator content identities for build reuse |
| `run/input.py` | Mechanism resolution and generation for each target |

The generator supports ideal-gas mixtures with NASA7 thermodynamics, elementary
and third-body Arrhenius reactions, and Troe or Lindemann falloff reactions.
It requires gas transport data and positive Arrhenius pre-exponential factors.
Other rate types, custom reaction orders and other thermodynamic models are
rejected before emission with a species or reaction identifier. Cantera's ability
to parse a mechanism does not imply that the generator supports every feature in it.

The generated module provides species metadata, caloric and ideal-gas properties,
temperature inversion, net production rates, fused creation/destruction rates,
mixture viscosity and thermal conductivity, and species diffusivities. MFC still
owns reaction time integration, including alpha-QSS, and spatial transport
discretization. The module follows MFC's source conventions: it uses `wp` from
`m_precision_select` and marks device routines with `$:GPU_ROUTINE`, so one generated
source serves every precision and offload configuration and compiler-specific
directive handling stays in MFC's macros. Nonchemistry builds retain the existing
dummy `h2o2.yaml` module to satisfy the shared Fortran interfaces.

The initial ownership change preserves the previous numerical formulas. In
particular, the pure-species diffusion limit retains the self-diffusion coefficient;
Cantera's mixture diffusion query may return zero for that degenerate state.
Temperature inversion also retains the existing Newton iteration and tolerance.
Real constants use the selected working precision: single for single-precision
builds and double for default and mixed-storage builds. Single-precision results
can therefore differ in roundoff from the former double-literal expressions.
Troe falloff guards the logarithm at zero reduced pressure, including when a
compiler evaluates both arguments of a Fortran `merge` expression.

# Mixing-layer initial conditions

The temporal 2D, spatial 2D and temporal 3D reacting mixing-layer examples share
`toolchain/mfc/flamelet.py`. Their local `flamelet_ic.py` files retain grid generation,
coordinate conventions, perturbations and MFC input-file layouts.

The default cold profile uses the prescribed hyperbolic-tangent mixture fraction,
linearly mixed stream mass fractions and specific enthalpy, and Cantera's HP
temperature recovery. Velocity is interpolated between the prescribed stream
velocities. Density follows the ideal-gas mixture equation of state.

With `--hot`, Cantera solves a counterflow diffusion flame using unity-Lewis transport.
The counterflow width is ten vorticity thicknesses. Each inlet velocity is
`flame_strain_rate * width / 2`, with mass flux determined by its stream density;
`flame_strain_rate` is a nominal inlet strain parameter in inverse seconds, not the
computed local strain. Its default is 100/s and is set in each example's `case.py`.
The solution is tabulated against Bilger mixture fraction normalized to the inlet
compositions, then mapped to the mixing layer's prescribed tanh profile. Temperature
is recovered at the linearly mixed stream enthalpy. Extinguished flames or invalid
profiles fail explicitly instead of silently producing a cold initialization.

This hot initialization is a counterflow-based seed for the evolving MFC flow. It
replaces the former JAX flamelet solver and its iterative scalar-dissipation matching;
it is not numerically identical to that initialization or a steady solution of the
mixing-layer equations. Changes to the initialization method and mechanism invalidate
the example caches. No precomputed flame data or optional Pyrometheus installation
is required.

# Validation

After bootstrapping the toolchain, run:

```sh
PYTHONPATH=toolchain build/venv/bin/pytest -q \
    toolchain/mfc/test_thermochem.py toolchain/mfc/test_flamelet.py
./mfc.sh test --no-mpi -j 8 --only Chemistry
```

The kernel tests compile generated Fortran with GNU Fortran and compare with
Cantera for `h2o2.yaml`, `gri30.yaml`, the bundled San Diego mechanism, and the
hydrogen/xenon mechanism. They cover thermodynamics, energy/enthalpy inversion,
reaction production and destruction, elemental conservation, transport, single
precision, mixed-storage working-precision compatibility, long species names,
zero-concentration falloff with floating-point exception traps, and compilation
with OpenACC and OpenMP directives. They run the generated source through Fypp with
MFC's macros, as the build does. The directive tests
execute on the host; they do not validate GPU offload on accelerator hardware.

Initialization tests check stream limits, normalization, elemental composition,
enthalpy, density and a burning hot profile for both mixing-layer mechanisms.
