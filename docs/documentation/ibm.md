@page ibm Immersed Boundaries

# Immersed Boundaries

This page is a map of the immersed-boundary (IB) options.
The parameters and their constraints are in [Immersed Boundary Patches](@ref sec-immersed-boundary-patches) of the case-file documentation; this page only says which option does what and where it is checked.
Enable the model with `ib = T` and `num_ibs`; each body is a `patch_ib(i)` entry whose geometry is listed in [Immersed Boundary Patch Types](@ref immersed-boundary-patch-types).

## Motion

| Option | Effect |
| :--- | :--- |
| `moving_ibm = 0` | fixed body |
| `moving_ibm = 1` | one-way coupling: the body moves at its prescribed rate and acts on the fluid, which does not act back |
| `moving_ibm = 2` | two-way coupling: the body moves under pressure, viscous and body forces (`bf_x`, `bf_y`, `bf_z`) |
| `collision_model = 1` | soft-sphere collisions between circles or spheres and with walls; `coefficient_of_restitution`, `collision_time`, `collision_temporal_resolution`, `ib_coefficient_of_friction` |
| `ib_friction_slip_eps` | regularises the friction direction at small slip speed, so contacts split across MPI ranks agree with a single rank |

Moving-body kinematics and the restitution of a wall collision are checked in [Validation](@ref sec-validation-moving-kinematics).

## Wall treatment

| Option | Effect |
| :--- | :--- |
| `ib_second_order_vel`, `ib_ip_min_dist` | linear ghost-velocity extrapolation; minimum image-point distance |
| `ib_ip_interp = 2` | triangular-shaped-cloud image-point interpolation instead of inverse-distance weights |
| `ib_wall_order = 2` | two-probe quadratic wall-normal ghost velocity at no-slip walls |
| `ib_slip_ccst`, `ib_smooth_contact`, `ib_uncover_fill` | slip-wall ghost state, HLLC contact blending in the viscous work, fill of cells a moving body uncovers |
| `weno_Re_flux`, `ib_visc_compact` | compact velocity gradient on the viscous faces within two cells of a body |
| `ib_mass_relax` | returns the spurious gas mass of a body over about N steps; the default resolves at start-up (see the parameter entry) |

Each option has a default that leaves existing results unchanged, except `ib_mass_relax`, whose default depends on whether a body moves.
Drag and Nusselt number of a cylinder: @ref validation_cyl_drag_nu. Interpolation against an exact conduction solution: @ref validation_tsc_annulus. Slip walls (mass budget, exact vortex, known failures): @ref validation_slip_walls. Lumped cooling and the wall conduction column: @ref validation_wall_thermal.

## Heat transfer and chemistry

`thermal_bc` selects the wall condition: zero-gradient (0), prescribed `Twall` (1), reacting-surface energy balance (2), lumped solid with an evolving temperature (3).
`wall_column` adds a sub-grid conduction column to a fixed isothermal circle or sphere, and `ib_flux_wall` sets the face fluxes of a reacting wall from its surface solve.
Examples: `examples/2D_ibm_lumped_wall`, `examples/2D_ibm_wall_mass_budget`, `examples/2D_ibm_conduction_annulus`.

## Output and restart

- `ib_force_wrt` writes `ib_forces.dat` ([force history](@ref sec-ib-force-history)); `ib_heat_wrt` writes `ib_heat.dat` ([heat history](@ref sec-ib-heat-history)); `ib_surface_wrt` writes surface-point records for a chemistry run.
- `ib_state_wrt` writes the IB state as a point mesh in post-processing.
- The IB state is saved under `restart_data/` (`ib_state_*.dat`; the conduction temperature of ghost points in `ib_tcond_<step>.dat`).
