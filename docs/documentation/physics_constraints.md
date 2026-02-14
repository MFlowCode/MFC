@page physics_constraints Physics Constraints

# Physics Constraints Reference

This document catalogs the physics constraints enforced by MFC's case parameter validator.
Constraints are organized by physical category with mathematical justifications.
For parameter syntax and allowed values, see @ref case "Case Files" and the @ref parameters "Case Parameters" reference.

The validator lives in `toolchain/mfc/case_validator.py` and runs automatically before each MFC stage.
Hard constraint violations produce **errors** that abort the run.
Soft constraint violations produce **warnings** that flag likely mistakes without stopping execution.

---

## 1. Thermodynamic Constraints

### 1.1 Positive Pressure

\f[p > 0\f]

All initial patch pressures (`patch_icpp(i)%pres`) must be strictly positive.
Negative or zero pressure is unphysical for the stiffened gas and other equation-of-state models used by MFC.

**Stage:** pre_process | **Severity:** error

### 1.2 Non-negative Density

\f[\alpha_i \rho_i \geq 0\f]

Partial densities (`patch_icpp(i)%alpha_rho(j)`) must be non-negative.
Zero partial density is allowed for vacuum regions.

**Stage:** pre_process | **Severity:** error

### 1.3 EOS Parameter Sanity (Transformed Gamma)

MFC uses the **transformed** stiffened gas parameter:

\f[\Gamma = \frac{1}{\gamma - 1}\f]

where \f$\gamma\f$ is the physical specific heat ratio.
A common mistake is entering the physical \f$\gamma\f$ directly (e.g., 1.4 for air) instead of the transformed value \f$1/(1.4-1) = 2.5\f$.

The validator warns when:
- `fluid_pp(i)%gamma < 0.1` (implies physical \f$\gamma > 11\f$, unusually high)
- `fluid_pp(i)%gamma > 1000` (implies physical \f$\gamma \approx 1.001\f$, unusually close to 1)

Similarly, `pi_inf` is stored as \f$\gamma \pi_\infty / (\gamma - 1)\f$.

**Stage:** common (all stages) | **Severity:** warning

### 1.4 Stiffened EOS Positivity

\f[\Gamma > 0, \quad \Pi_\infty \geq 0, \quad c_v \geq 0\f]

The equation-of-state parameters `fluid_pp(i)%gamma`, `fluid_pp(i)%pi_inf`, and `fluid_pp(i)%cv` must satisfy basic positivity requirements for thermodynamic stability.

**Stage:** common | **Severity:** error

---

## 2. Mixture Constraints

### 2.1 Volume Fraction Sum

For multi-component models (`model_eqns` \f$\in \{2, 3, 4\}\f$), the volume fractions must satisfy the mixture constraint:

\f[\sum_{i=1}^{N_f} \alpha_i = 1\f]

The validator checks this per patch and warns if the deviation exceeds \f$10^{-6}\f$.

**Exceptions** (constraint does not apply):
- Single-fluid Euler-Euler bubble models (`bubbles_euler = T`, `num_fluids = 1`): \f$\alpha\f$ represents void fraction
- Lagrangian bubble models (`bubbles_lagrange = T`): Lagrangian phase is not tracked on the Euler grid
- IBM cases (`num_ibs > 0`): \f$\alpha\f$ acts as a level-set indicator
- Alter patches and hard-coded IC (hcid) patches: values are computed at runtime
- Analytical expressions (strings): cannot be validated statically

**Stage:** pre_process | **Severity:** warning

### 2.2 Alpha-Rho Consistency

The validator warns about physically inconsistent combinations:

- \f$\alpha_j = 0\f$ but \f$\alpha_j \rho_j \neq 0\f$: density assigned to an absent phase
- \f$\alpha_j > 10^{-10}\f$ but \f$\alpha_j \rho_j = 0\f$: present phase has zero density

These are not strictly errors (the solver can handle them) but usually indicate a configuration mistake.

**Stage:** pre_process | **Severity:** warning

### 2.3 Volume Fraction Bounds

\f[0 \leq \alpha_i \leq 1\f]

Individual volume fractions must be non-negative and (in non-IBM cases) at most 1.

**Stage:** pre_process | **Severity:** error

---

## 3. Domain and Geometry Constraints

### 3.1 Domain Bounds

For each active spatial dimension:

\f[x_{\mathrm{end}} > x_{\mathrm{beg}}, \quad y_{\mathrm{end}} > y_{\mathrm{beg}}, \quad z_{\mathrm{end}} > z_{\mathrm{beg}}\f]

The domain must have positive extent. A reversed or zero-width domain is always a configuration error.

**Stage:** common | **Severity:** error

### 3.2 Positive Patch Dimensions

Patch geometry parameters (`length_x`, `length_y`, `length_z`, `radius`) must be positive.
Exception: in cylindrical coordinates, `length_y` and `length_z` may use sentinel values.

**Stage:** pre_process | **Severity:** error

### 3.3 Patch Within Domain

For patches with centroid + length geometry (line segments, rectangles, cuboids), the validator checks that the patch bounding box is not entirely outside the computational domain.

Skipped when grid stretching is active (physical coordinates are transformed).

**Stage:** pre_process | **Severity:** error

### 3.4 Dimensionality

- \f$m > 0\f$ is required (x-direction must have cells)
- \f$n \geq 0\f$, \f$p \geq 0\f$
- If \f$n = 0\f$ then \f$p = 0\f$ (cannot have z without y)
- Cylindrical coordinates (\f$p > 0\f$): \f$p\f$ must be odd

**Stage:** common | **Severity:** error

---

## 4. Velocity and Dimensional Consistency

### 4.1 Velocity Components in Inactive Dimensions

\f[n = 0 \implies v_{2} = 0, \quad p = 0 \implies v_{3} = 0\f]

Setting velocity components in dimensions that do not exist is almost certainly a mistake.

**Exception:** MHD simulations legitimately use transverse velocity components in 1D because they carry transverse momentum coupled to the magnetic field.

**Stage:** pre_process | **Severity:** error

### 4.2 Momentum and Velocity Output Constraints

Post-process outputs `mom_wrt(2)`, `vel_wrt(2)` require \f$n > 0\f$; `mom_wrt(3)`, `vel_wrt(3)` require \f$p > 0\f$.

**Stage:** post_process | **Severity:** error

---

## 5. Model Equation Compatibility

### 5.1 Model Selection

| `model_eqns` | Name | \f$N_f\f$ | Key requirement |
|:---:|------|:---:|-------|
| 1 | \f$\gamma\f$-law (single-fluid) | not set | No `num_fluids`, no bubbles, no `fluid_pp` |
| 2 | Five-equation (\cite Allaire02) | \f$\geq 1\f$ | Primary workhorse model |
| 3 | Six-equation (\cite Saurel09) | \f$\geq 1\f$ | `riemann_solver = 2`, `avg_state = 2`, `wave_speeds = 1` |
| 4 | Four-equation | \f$= 1\f$ | Single-component with bubbles |

### 5.2 Key Incompatibilities

- `model_eqns = 1`: no `mpp_lim`, no viscosity Re parameters, no volume fraction output
- `model_eqns = 3`: no cylindrical 3D, no bubble models, requires HLLC solver
- `model_eqns = 4`: requires `num_fluids = 1`, bubble-specific

**Stage:** common | **Severity:** error

---

## 6. Boundary Conditions

### 6.1 Periodicity Matching

If one end of a dimension is periodic (`bc = -1`), the other end must also be periodic.

### 6.2 Boundary Condition Range

Valid BC values range from \f$-1\f$ to \f$-17\f$:

| Value | Boundary Type |
|:---:|------|
| -1 | Periodic |
| -2 | Reflective (slip wall) |
| -3 | Extrapolation (ghost cell) |
| -4 | Thompson (non-reflecting, sim only) |
| -5 to -12 | Characteristic BCs |
| -14 | Axis (cylindrical only) |
| -15 to -17 | Additional wall types |

### 6.3 Cylindrical Coordinate BCs

- 2D cylindrical (\f$p = 0\f$): `bc_y%beg = -2` (reflective at axis)
- 3D cylindrical (\f$p > 0\f$): `bc_y%beg = -14` (axis BC)
- 3D cylindrical z-direction: only periodic (-1) or reflective (-2)

**Stage:** common | **Severity:** error

---

## 7. Bubble Physics Constraints

### 7.1 Euler-Euler Bubbles (`bubbles_euler`)

- `nb >= 1` (number of bubble bins)
- Polydisperse requires odd `nb > 1` and `poly_sigma > 0`
- Not tested with `model_eqns` 1 or 3
- Reference quantities (`rhoref`, `pref`, `bub_pp%R0ref`, etc.) must be positive when set
- QBMM requires `nnode = 4`

### 7.2 Simulation-Specific Bubble Constraints

- Requires HLLC Riemann solver (`riemann_solver = 2`)
- Requires arithmetic average (`avg_state = 2`)
- Five-equation model does not support Gilmore (`bubble_model = 1`)
- Cannot use both `bubbles_euler` and `bubbles_lagrange` simultaneously

### 7.3 Euler-Lagrange Bubbles (`bubbles_lagrange`)

- 2D/3D only (`n > 0`)
- `file_per_process = F`
- Not compatible with `model_eqns = 3`
- Requires `polytropic = F` and `thermal = 3`

**Stage:** common + simulation | **Severity:** error

---

## 8. Feature Compatibility Matrix

Several physics models have mutual exclusion constraints. The key incompatibilities:

| Feature | Model | Riemann Solver | Other |
|---------|-------|---------------|-------|
| **MHD** | `= 2`, `num_fluids = 1` | HLL (1) or HLLD (4) | No relativity+HLLD |
| **Surface tension** | `= 2` or `= 3`, `num_fluids = 2` | — | — |
| **Hypoelasticity** | `= 2` | HLL (1) | — |
| **Hyperelasticity** | `= 2` or `= 3` | — | — |
| **Phase change** | `= 2` (relax 5,6) or `= 3` (relax 1,4,5,6) | — | — |
| **Alt sound speed** | `= 2`, `num_fluids` 2–3 | HLLC (2) | No bubbles |
| **IGR** | `= 2` | No characteristic BCs | No bubbles, MHD, elastic, etc. |

**Stage:** common + simulation | **Severity:** error

---

## 9. Numerical Scheme Constraints

### 9.1 WENO Reconstruction (`recon_type = 1`)

- `weno_order` \f$\in \{1, 3, 5, 7\}\f$
- Grid must have enough cells: \f$m + 1 \geq\f$ `num_stcls_min * weno_order`
- Schemes are mutually exclusive: only one of `mapped_weno`, `wenoz`, `teno`
- `teno` requires order 5 or 7; `mp_weno` requires order 5

### 9.2 MUSCL Reconstruction (`recon_type = 2`)

- `muscl_order` \f$\in \{1, 2\}\f$
- Second order requires `muscl_lim` \f$\in \{1, 2, 3, 4, 5\}\f$
- THINC interface compression (`int_comp`) requires MUSCL

### 9.3 Time Stepping

- `time_stepper` \f$\in \{1, 2, 3\}\f$
- `dt > 0` required for fixed time stepping
- CFL-based modes: `cfl_target` \f$\in (0, 1]\f$, `t_save \leq\f$ `t_stop`
- Adaptive dt (`adap_dt`): requires RK3, `polytropic = T` or `bubbles_lagrange = T`

### 9.4 Viscosity

- Reynolds numbers `Re(1)`, `Re(2)` must be positive
- Requires `viscous = T`
- Not supported with `model_eqns = 1`
- `weno_order = 1` without `weno_avg` does not support viscosity (unless IGR)

**Stage:** simulation | **Severity:** error

---

## 10. Acoustic Source Constraints

Acoustic sources have dimension-specific support types:

| Dimension | Allowed `support` values |
|:---------:|:----------------------:|
| 1D | 1 |
| 2D | 2, 5, 6, 9, 10 |
| 2D (cyl) | 2, 6, 10 |
| 3D | 3, 7, 11 |

Additional constraints:
- `pulse` \f$\in \{1, 2, 3, 4\}\f$
- Sinusoidal/square (1, 3): exactly one of `frequency` or `wavelength`
- Gaussian (2): exactly one of `gauss_sigma_time` or `gauss_sigma_dist`, plus `delay`
- Broadband (4): requires `bb_num_freq`, `bb_bandwidth`, `bb_lowest_freq`
- Non-planar sources (`support >= 5`): require `foc_length` and `aperture`

**Stage:** simulation | **Severity:** error

---

## 11. Post-Processing Constraints

### 11.1 Vorticity and Schlieren

- `omega_wrt` and `schlieren_wrt` require at least 2D (\f$n > 0\f$)
- 3D vorticity components (`omega_wrt(1)`, `omega_wrt(2)`) require \f$p > 0\f$
- Both require `fd_order` to be set

### 11.2 FFT Output

- Requires 3D (\f$n > 0\f$, \f$p > 0\f$)
- All boundaries must be periodic
- Global dimensions must be even
- Incompatible with cylindrical coordinates

### 11.3 Output Selection

At least one flow variable must be selected for post-processing output.

**Stage:** post_process | **Severity:** error

---

## References

The physics models and their constraints are described in detail in:

- \cite Wilfong26 — MFC 5.0: comprehensive description of all models
- \cite Bryngelson21 — MFC: An open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver
- \cite Allaire02 — Five-equation model for compressible two-phase flow
- \cite Saurel09 — Six-equation model with relaxation
- \cite Kapila01 — Two-phase modeling with interface mechanics

For the auto-generated constraint reference with compatibility tables and working examples, see @ref case_constraints "Case Creator Guide".
