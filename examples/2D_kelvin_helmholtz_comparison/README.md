# Matched Kelvin–Helmholtz EOS comparison

Run the three standalone folders, then from the repository root:

```bash
build/venv/bin/python examples/2D_kelvin_helmholtz_comparison/compare.py
```

The script runs each verifier, requires matching grids, normalized initial fields and
saved times, and writes `viz/density_comparison.png`, `growth.png`, and `comparison.json`.
It measures differences without imposing a desired supercritical result. Refine all
three grids together before interpreting small differences or late-time fine structures.

## Physical design

| Folder suffix | EOS | Pressure | Initial temperature |
|---|---|---|---|
| `ideal_gas` | gamma = 5/3 | 101325 Pa | approximately 300–600 K |
| `pr_normal` | Peng–Robinson | 101325 Pa | approximately 302–600 K |
| `pr_supercritical` | Peng–Robinson | 10 MPa | approximately 330–411 K |

All folders start with the same normalized density (1 to 2), streamwise velocity
(-0.5 to +0.5), transverse perturbation (amplitude 0.1), and pressure (2.5).
The physical domain is a 1 cm periodic square. Reference velocity is
`sqrt(p / (2.5*rho_ref))`; comparison time is `t* = t U_ref / L`.
Pressure, density, and velocity are matched as independent initial variables;
temperature and sound speed follow each EOS. Matching temperature as well would
violate one EOS or change the density contrast. Mach numbers are deliberately allowed
to respond to the EOS. The ideal-gas normalized solution also provides the ideal
control at the supercritical run's dimensional scales because the inviscid ideal-gas
Euler equations have this scaling symmetry. This does not imply equal dimensional
fields across pressures.

Xenon retains the original example's monatomic gamma = 5/3 in the dilute limit.
Its properties are R = 8.31446261815324 / 0.131293 J/(kg K), Tc = 289.732568 K,
Pc = 5.841914 MPa, and acentric factor 0.00363, from
[CoolProp's xenon fluid reference](https://coolprop.org/fluid_properties/fluids/Xenon.html).
The caloric closure uses constant ideal cv = 1.5 R plus the PR departure energy.
This is a cubic-EOS demonstration, not a claim of reference-quality near-critical xenon
properties. The verifier checks positive stable thermodynamics at every saved frame
and requires the supercritical case to remain above Tc and Pc at those frames.

## Differences from the original example

`../2D_kelvin_helmholtz` is preserved. It uses two labeled fluids, HLL, a discontinuous
shear layer, and a 513 × 513 grid. At setup time its directory contained reference
images but no raw simulation files. Its images cannot isolate an EOS effect from
these numerical/model differences.

The matched controls use one fluid (the current PR backend's supported model), HLLC,
and a smooth periodic shear profile of thickness parameter 0.025 L. The perturbation
uses the original two wavelengths and Gaussian width. All three matched runs use
identical numerical settings and exact matching nondimensional output times. The
finite shear thickness reduces the grid dependence of an unresolved vortex sheet;
128 × 128 is an accessible Mac starting point, not a convergence claim.

Periodic mass, both momenta, and total energy must be conserved within 1e-8 of their
specified characteristic scales. The comparison reports normalized density/velocity/
pressure RMS differences and transverse kinetic-energy growth. Reported density uses
`rho_ref`, velocity uses `U_ref`, and pressure uses the initial pressure as their scales;
kinetic energy uses `rho_ref * U_ref**2`. Energy and temperature
are not expected to match across EOS at the same initial pressure and density.
No physical viscosity, heat conduction, or phase-equilibrium model is included.

## Initial validation (2026-11-18)

All three examples completed on this Mac with 64 × 64 cells, serial CPU execution,
101 saves, and `t* = 2`. The saved case defaults remain 128 × 128 with MPI I/O enabled.
These are coarse demonstration results, not a grid-convergence study.

| Difference from matched ideal gas at t* = 2 | PR, 1 atm | PR, 10 MPa |
|---|---:|---:|
| Density RMS / rho_ref | 0.05015% | 4.5076% |
| Streamwise velocity RMS / U_ref | 0.03541% | 1.2367% |
| Transverse velocity RMS / U_ref | 0.04161% | 1.0311% |
| Pressure RMS / initial pressure | 0.05438% | 2.5958% |
| Relative change in transverse kinetic energy | -0.04674% | +4.1785% |

The dense case's density RMS difference is about 90 times the dilute case's.
All three runs conserved periodic mass, both momenta, and total energy to better than
6e-16 using the verifier's characteristic scales. In the supercritical run, the
minimum over saved frames was 320.697 K and 7.75726 MPa, above Tc and Pc.
The 30 existing PR/configuration checks and the added independent Cantera xenon check
passed. The local animations were also generated successfully.

The observed differences include the numerical discretization's response to the EOS.
Refine all three cases together (for example, 128 then 256 cells in each direction)
and check the growth histories and RMS differences before assigning physical significance
to fine vortex structures. This demonstration does not validate a real-fluid turbulence
or heat-transfer model.
