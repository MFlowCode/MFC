@page equations Equations

# MFC: Comprehensive Equations Reference

This document catalogs every equation solved by MFC, organized by physical model.
Each section notes the input parameter(s) that activate the corresponding physics module and cross-references the relevant source files.

The models and algorithms described here are detailed in \cite Wilfong26 (MFC 5.0) and \cite Bryngelson21. Foundational references for each model are cited inline; see the \ref citelist "Bibliography" for full details.

For parameter details and allowed values, see @ref case "Case Files" and the @ref parameters "Case Parameters" reference.

---

## 1. Overview

MFC solves the compressible Navier-Stokes equations (or Euler equations when viscosity is off) in a finite volume framework. The general semi-discrete form is:

\f[\frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{q}) + \mathbf{h}(\mathbf{q})\,\nabla \cdot \mathbf{u} = \mathbf{s}(\mathbf{q})\f]

where:
- \f$\mathbf{q}\f$ is the conservative variable vector,
- \f$\mathbf{F}\f$ is the flux tensor,
- \f$\mathbf{h}(\mathbf{q})\,\nabla \cdot \mathbf{u}\f$ contains non-conservative terms (volume fraction advection),
- \f$\mathbf{s}(\mathbf{q})\f$ is the source vector (bubbles, body forces, chemistry, etc.).

The parameter `model_eqns` (1, 2, 3, or 4) selects the governing equation set.

**Key source files:** `src/simulation/m_rhs.fpp` (RHS evaluation), `src/common/m_variables_conversion.fpp` (EOS and variable conversion).

---

## 1b. Units, Dimensions, and Non-Dimensionalization {#sec-units-dimensions}

### General Users: Dimensional Handling {#sec-dimensional-handling}

#### Dimensions In = Dimensions Out {#sec-dimensions-in-out}

The main flow solver (Navier-Stokes equations, Riemann solvers, viscous stress, body forces, surface tension, etc.) is **unit-agnostic**: whatever units the user provides for the initial and boundary conditions, the solver preserves them throughout the computation. If the user inputs SI units, the outputs are in SI units. If the user inputs CGS, the outputs are in CGS. No internal non-dimensionalization is performed by the flow solver.

This means that for simulations **without** sub-grid bubble models, the user can work in any consistent unit system without additional effort.

#### Stored Parameter Conventions {#sec-stored-forms}

Several EOS and transport parameters use **transformed stored forms** that differ from the standard physical values. This is the most common source of input errors:

| Parameter | Physical quantity | What MFC expects (stored form) |
|---|---|---|
| `fluid_pp(i)%%gamma` | Heat capacity ratio \f$\gamma\f$ | \f$\Gamma = \frac{1}{\gamma - 1}\f$ |
| `fluid_pp(i)%%pi_inf` | Stiffness pressure \f$\pi_\infty\f$ [Pa] | \f$\Pi_\infty = \frac{\gamma\,\pi_\infty}{\gamma - 1}\f$ [Pa] |
| `fluid_pp(i)%%Re(1)` | Dynamic viscosity \f$\mu\f$ | \f$1/\mu\f$ (inverse viscosity) |
| `fluid_pp(i)%%Re(2)` | Bulk viscosity \f$\mu_b\f$ | \f$1/\mu_b\f$ (inverse bulk viscosity) |

These transformations arise because MFC internally solves the energy equation using the transformed variables \f$\Gamma\f$ and \f$\Pi_\infty\f$ (see Section 3.1), and the viscous stress is computed by dividing by `Re` rather than multiplying by \f$\mu\f$.

**Common mistake:** setting `fluid_pp(1)%%gamma = 1.4` for air. The correct value is `1.0 / (1.4 - 1.0) = 2.5`. Setting `gamma = 1.4` corresponds to a physical \f$\gamma \approx 1.71\f$, which is not a standard gas.

#### Common Material Values {#sec-material-values}

Pre-computed stored-form values for common fluids (SI units):

| Material | \f$\gamma\f$ | \f$\pi_\infty\f$ [Pa] | `gamma` (stored) | `pi_inf` (stored) [Pa] |
|---|---|---|---|---|
| Air | 1.4 | 0 | 2.5 | 0 |
| Helium | 5/3 | 0 | 1.5 | 0 |
| Water (Tait) | 4.4 | 6.0e8 | 0.2941 | 7.76e8 |
| Water (\cite LeMetayer04) | 6.12 | 3.43e8 | 0.1953 | 4.10e8 |

Example for an air-water simulation:

```python
# Air (fluid 1)
gam_a  = 1.4
"fluid_pp(1)%gamma":  1.0 / (gam_a - 1.0),   # = 2.5
"fluid_pp(1)%pi_inf": 0.0,

# Water (fluid 2)
gam_w  = 4.4
pi_w   = 6.0e8  # Pa
"fluid_pp(2)%gamma":  1.0 / (gam_w - 1.0),            # ≈ 0.294
"fluid_pp(2)%pi_inf": gam_w * pi_w / (gam_w - 1.0),   # ≈ 7.76e8
```

For viscous cases, provide the **reciprocal** of the dynamic viscosity:

```python
mu = 1.002e-3  # water viscosity [Pa·s]
"fluid_pp(1)%Re(1)": 1.0 / mu,  # ≈ 998
```

#### Unit Consistency {#sec-unit-consistency}

The solver does not check or convert units. All inputs must use the **same consistent unit system** (e.g., all SI or all CGS). Mixing units — for example, pressures in atmospheres with densities in kg/m³ — will produce silently incorrect results.

### Bubble Users: Non-Dimensional Framework {#sec-bubble-nondim}

#### Non-Dimensional Bubble Dynamics {#sec-nondim-bubble-dynamics}

The sub-grid bubble models (`bubbles_euler = .true.` or `bubbles_lagrange = .true.`) solve the bubble wall dynamics in **non-dimensional form**. The bubble wall pressure equation as implemented is:

\f[p_{bw} = \left(\text{Ca} + \frac{2}{\text{We}_b\,R_0}\right)\left(\frac{R_0}{R}\right)^{3\gamma} - \text{Ca} - \frac{4\,\text{Re}_{\text{inv}}\,\dot{R}}{R} - \frac{2}{R\,\text{We}_b}\f]

Here \f$R\f$ and \f$R_0\f$ are non-dimensional radii (scaled by \f$x_0\f$), and \f$\dot{R}\f$ is a non-dimensional wall speed (scaled by \f$u_0\f$); the entire bubble ODE is solved in non-dimensional variables.

The dimensionless groups are:

| Dimensionless group | Definition | Code variable | Computed from |
|---|---|---|---|
| \f$\text{Ca}\f$ (Cavitation number) | \f$p_{0,\text{ref}} - p_v\f$ | `Ca` | `bub_pp%%p0ref - bub_pp%%pv` |
| \f$\text{Eu}\f$ (Euler number) | \f$p_{0,\text{ref}}\f$ | `Eu` | `bub_pp%%p0ref` |
| \f$\text{We}_b\f$ (bubble Weber number) | \f$1/\sigma\f$ | `Web` | `1 / bub_pp%%ss` |
| \f$\text{Re}_{\text{inv}}\f$ (inverse bubble Reynolds number) | \f$\mu_l\f$ | `Re_inv` | `bub_pp%%mu_l` |

Because the bubble equations use these dimensionless numbers directly, all `bub_pp%%` inputs are interpreted by the code as **already non-dimensional**. The code does **not** non-dimensionalize bubble quantities internally. Therefore, when bubbles are enabled, the simulation must be run in a **fully non-dimensional** form: **all** inputs — flow ICs/BCs, EOS parameters, domain lengths, `dt`, and `bub_pp%%` values — must be scaled with the same \f$(x_0, p_0, \rho_0, u_0, t_0, T_0)\f$ reference quantities, or the coupled solution will be physically incorrect.

#### Reference Scales {#sec-reference-scales}

When using bubble models, the user must choose reference scales and non-dimensionalize **all** inputs (flow and bubble) consistently. The standard convention used in the MFC examples is:

| Reference quantity | Symbol | Typical choice |
|---|---|---|
| Length | \f$x_0\f$ | \f$R_{0,\text{ref}}\f$ (reference bubble radius) |
| Pressure | \f$p_0\f$ | \f$p_{0,\text{ref}}\f$ (reference bubble pressure) |
| Density | \f$\rho_0\f$ | \f$\rho_{0,\text{ref}}\f$ (reference liquid density) |
| Velocity | \f$u_0\f$ | \f$\sqrt{p_0 / \rho_0}\f$ (derived) |
| Time | \f$t_0\f$ | \f$x_0 / u_0\f$ (derived) |
| Temperature | \f$T_0\f$ | \f$T_{0,\text{ref}}\f$ (reference temperature) |

#### Non-Dimensionalization of Input Parameters {#sec-nondim-inputs}

The following table lists every `bub_pp%%` parameter and its required non-dimensionalization:

| Parameter | Physical meaning | Non-dimensional form |
|---|---|---|
| `bub_pp%%R0ref` | Reference bubble radius | \f$R_{0,\text{ref}} / x_0\f$ |
| `bub_pp%%p0ref` | Reference bubble pressure | \f$p_{0,\text{ref}} / p_0\f$ |
| `bub_pp%%rho0ref` | Reference liquid density | \f$\rho_{0,\text{ref}} / \rho_0\f$ |
| `bub_pp%%T0ref` | Reference temperature | \f$T_{0,\text{ref}} / T_0\f$ (typically 1) |
| `bub_pp%%ss` | Surface tension \f$\sigma\f$ | \f$\sigma / (\rho_0\,x_0\,u_0^2)\f$ |
| `bub_pp%%pv` | Vapor pressure | \f$p_v / p_0\f$ |
| `bub_pp%%mu_l` | Liquid dynamic viscosity | \f$\mu_l / (\rho_0\,x_0\,u_0)\f$ |
| `bub_pp%%mu_v` | Vapor dynamic viscosity | \f$\mu_v / (\rho_0\,x_0\,u_0)\f$ |
| `bub_pp%%mu_g` | Gas dynamic viscosity | \f$\mu_g / (\rho_0\,x_0\,u_0)\f$ |
| `bub_pp%%vd` | Vapor diffusivity | \f$D / (x_0\,u_0)\f$ |
| `bub_pp%%k_v` | Vapor thermal conductivity | \f$k_v\,T_0 / (x_0\,\rho_0\,u_0^3)\f$ |
| `bub_pp%%k_g` | Gas thermal conductivity | \f$k_g\,T_0 / (x_0\,\rho_0\,u_0^3)\f$ |
| `bub_pp%%cp_v` | Vapor specific heat | \f$c_{p,v}\,T_0 / u_0^2\f$ |
| `bub_pp%%cp_g` | Gas specific heat | \f$c_{p,g}\,T_0 / u_0^2\f$ |
| `bub_pp%%R_v` | Vapor gas constant | \f$R_v\,T_0 / u_0^2\f$ |
| `bub_pp%%R_g` | Gas gas constant | \f$R_g\,T_0 / u_0^2\f$ |
| `bub_pp%%gam_v` | Vapor heat capacity ratio | Already dimensionless (no scaling) |
| `bub_pp%%gam_g` | Gas heat capacity ratio | Already dimensionless (no scaling) |
| `bub_pp%%M_v` | Vapor molar mass | Consistent units; only ratios are used (no scaling needed) |
| `bub_pp%%M_g` | Gas molar mass | Consistent units; only ratios are used (no scaling needed) |

When the reference scales match the bubble reference values (e.g., \f$x_0 = R_{0,\text{ref}}\f$, \f$p_0 = p_{0,\text{ref}}\f$, \f$\rho_0 = \rho_{0,\text{ref}}\f$), the reference parameters simplify to unity: `bub_pp%%R0ref = 1`, `bub_pp%%p0ref = 1`, `bub_pp%%rho0ref = 1`.

#### Flow Parameters with Bubbles {#sec-flow-params-bubbles}

When bubbles are enabled, the flow-level parameters must also be non-dimensionalized with the same reference scales:

| Parameter | Non-dimensional form |
|---|---|
| `x_domain%%beg`, `x_domain%%end` | Domain bounds divided by \f$x_0\f$ |
| `patch_icpp(i)%%pres` | Pressure divided by \f$p_0\f$ |
| `patch_icpp(i)%%alpha_rho(j)` | Partial density divided by \f$\rho_0\f$ |
| `patch_icpp(i)%%vel(j)` | Velocity divided by \f$u_0\f$ |
| `fluid_pp(i)%%gamma` | \f$1/(\gamma_i - 1)\f$ (dimensionless, same as without bubbles) |
| `fluid_pp(i)%%pi_inf` | \f$\gamma_i\,\pi_{\infty,i} / [(\gamma_i - 1)\,p_0]\f$ (scaled by reference pressure) |
| `fluid_pp(i)%%Re(1)` | \f$\rho_0\,x_0\,u_0 / \mu_i\f$ (Reynolds number, inverse viscosity) |
| `dt` | Time step divided by \f$t_0\f$ |

#### Two Different Viscosity Parameters {#sec-two-viscosities}

MFC has two conceptually distinct viscosity-related parameters that serve different physical roles:

1. **`fluid_pp(i)%%Re(1)`** — Used for the **macroscopic flow viscous stress tensor** (Navier-Stokes equations). This is \f$1/\mu\f$ in dimensional simulations, or \f$\rho_0 x_0 u_0 / \mu\f$ (a Reynolds number) when non-dimensionalized. It appears as a **divisor** in the viscous stress computation:
   \f[\tau_{ij} \propto \frac{\nabla u}{\text{Re}}\f]
   Stored in the physical\_parameters derived type (`src/common/m_derived_types.fpp`).

2. **`bub_pp%%mu_l`** — Used for **microscale bubble wall viscous damping** (Rayleigh-Plesset / Keller-Miksis equations). This is the non-dimensional liquid viscosity \f$\mu_l / (\rho_0 x_0 u_0)\f$. It appears as a **multiplier** in the bubble wall pressure:
   \f[p_{bw} \ni -\frac{4\,\text{Re}_{\text{inv}}\,\dot{R}}{R}\f]
   Stored in the subgrid\_bubble\_physical\_parameters derived type (`src/common/m_derived_types.fpp`).

These two parameters represent viscous effects at fundamentally different scales — bulk flow dissipation vs. single-bubble-wall damping — and are stored in separate derived types with separate code paths. They are **not** interchangeable: `fluid_pp%%Re(1)` is an inverse viscosity while `bub_pp%%mu_l` is a viscosity (non-dimensionalized).

### Worked Examples {#sec-nondim-example}

#### Example: Non-Dimensionalizing a Bubble Case {#sec-bubble-example}

A typical bubble case setup in `case.py` follows this pattern:

```python
import math

# Physical properties (SI units)
rho_l  = 1.0e03     # liquid density [kg/m³]
mu_l   = 1.002e-03  # liquid viscosity [kg/(m·s)]
ss     = 0.07275    # surface tension [kg/s²]
pv     = 2.3388e03  # vapor pressure [Pa]
gam_l  = 7.15       # liquid stiffened gas gamma
pi_inf = 306.0e06   # liquid stiffened gas pi_inf [Pa]

# Bubble reference values (SI)
R0ref   = 10.0e-06  # reference bubble radius [m]
p0ref   = 112.9e03  # reference bubble pressure [Pa]
rho0ref = rho_l     # reference density [kg/m³]

# Derived reference scales
x0   = R0ref
p0   = p0ref
rho0 = rho0ref
u0   = math.sqrt(p0 / rho0)
t0   = x0 / u0

# Non-dimensional inputs
params = {
    "bub_pp%R0ref":   R0ref / x0,                     # = 1.0
    "bub_pp%p0ref":   p0ref / p0,                      # = 1.0
    "bub_pp%rho0ref": rho0ref / rho0,                  # = 1.0
    "bub_pp%ss":      ss / (rho0 * x0 * u0**2),       # surface tension
    "bub_pp%pv":      pv / p0,                         # vapor pressure
    "bub_pp%mu_l":    mu_l / (rho0 * x0 * u0),        # liquid viscosity

    "fluid_pp(1)%gamma":  1.0 / (gam_l - 1.0),
    "fluid_pp(1)%pi_inf": gam_l * (pi_inf / p0) / (gam_l - 1.0),
    "fluid_pp(1)%Re(1)":  rho0 * x0 * u0 / mu_l,     # flow Re (inverse!)
}
```

Note the inverse relationship: `fluid_pp%%Re(1) = 1 / bub_pp%%mu_l` when both use the same reference scales and the same physical viscosity. This is expected — they encode the same physical viscosity but in reciprocal forms for their respective equations.

---

## 2. Governing PDEs

### 2.1 Five-Equation Model (`model_eqns = 2`)

The primary workhorse model (\cite Allaire02; \cite Wilfong26 Sec. 2.1). The state vector is:

\f[\mathbf{q} = \bigl(\alpha_1 \rho_1,\;\alpha_2 \rho_2,\;\ldots,\;\rho u_1,\;\rho u_2,\;\rho u_3,\;\rho E,\;\alpha_1,\;\alpha_2,\;\ldots\bigr)^T\f]

**Continuity** (one per component):

\f[\frac{\partial (\alpha_i \rho_i)}{\partial t} + \nabla \cdot (\alpha_i \rho_i\,\mathbf{u}) = 0\f]

**Momentum:**

\f[\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot \bigl(\rho\,\mathbf{u} \otimes \mathbf{u} + p\,\mathbf{I} - \boldsymbol{\tau}^v\bigr) = 0\f]

**Energy:**

\f[\frac{\partial (\rho E)}{\partial t} + \nabla \cdot \bigl[(\rho E + p)\,\mathbf{u} - \boldsymbol{\tau}^v \cdot \mathbf{u}\bigr] = 0\f]

**Volume fraction advection:**

\f[\frac{\partial \alpha_i}{\partial t} + \mathbf{u} \cdot \nabla \alpha_i = K\,\nabla \cdot \mathbf{u}\f]

where the \f$K\f$ term enforces interface conditions via the Wood sound speed:

\f[K = \frac{\rho_2 c_2^2 - \rho_1 c_1^2}{\displaystyle\frac{\rho_1 c_1^2}{\alpha_1} + \displaystyle\frac{\rho_2 c_2^2}{\alpha_2}}\f]

Setting `alt_soundspeed = .true.` enables the \f$K\f$ correction (\cite Kapila01, with Wood sound speed). Setting `alt_soundspeed = .false.` uses the Allaire variant without the \f$K\f$ correction, which is conservative but does not strictly obey the second law of thermodynamics.

**Mixture rules:**

\f[1 = \sum_i \alpha_i, \qquad \rho = \sum_i \alpha_i \rho_i, \qquad \rho e = \sum_i \alpha_i \rho_i e_i\f]

### 2.2 Six-Equation Model (`model_eqns = 3`)

Allows pressure disequilibrium between phases (\cite Saurel09; \cite Wilfong26 Sec. 2.1).

**Continuity and momentum:** Same as the five-equation model.

**Separate phasic internal energy:**

\f[\frac{\partial (\alpha_i \rho_i e_i)}{\partial t} + \nabla \cdot (\alpha_i \rho_i e_i\,\mathbf{u}) + \alpha_i p_i\,\nabla \cdot \mathbf{u} = -\mu\,p_I\,(p_2 - p_1) - \alpha_i\,\boldsymbol{\tau}_i^v : \nabla \mathbf{u}\f]

**Volume fraction:**

\f[\frac{\partial \alpha_1}{\partial t} + \mathbf{u} \cdot \nabla \alpha_1 = \mu\,(p_1 - p_2)\f]

**Interfacial pressure:**

\f[p_I = \frac{z_2\,p_1 + z_1\,p_2}{z_1 + z_2}, \qquad z_i = \rho_i\,c_i\f]

Infinite pressure relaxation is applied at each Runge-Kutta stage to drive toward pressure equilibrium.

**Mixture speed of sound:**

\f[c^2 = \sum_k Y_k\,c_k^2\f]

With phase change (`relax = .true.`), additional source terms appear in the phasic energy and volume fraction equations:
- **Pressure relaxation:** \f$\mu\,\delta p\f$ where \f$\delta p = p_1 - p_2\f$
- **Thermal transfer:** \f$Q = \theta\,(T_2 - T_1)\f$
- **Mass transfer:** \f$\dot{m} = \nu\,(g_2 - g_1)\f$ (Gibbs free energy difference)

See Section 8 (Phase Change) below for details.

### 2.3 Other Model Variants

- `model_eqns = 1`: **Gamma/pi_inf model** — simplified single-fluid formulation using mixture \f$\gamma\f$ and \f$\pi_\infty\f$ directly without tracking individual volume fractions (\cite Johnsen08).
- `model_eqns = 4`: **Four-equation model** — reduced model from the six-equation system after full pressure-temperature equilibrium relaxation (Tait-like compressible liquid).

---

## 3. Equations of State

### 3.1 Stiffened Gas EOS (\cite Menikoff89; \cite LeMetayer04; \cite Wilfong26 Sec. 2.2)

The primary closure for each phase:

\f[p_k = (\gamma_k - 1)\,\rho_k\,e_k - \gamma_k\,\pi_{\infty,k}\f]

Equivalently:

\f[e_k = \frac{p_k + \gamma_k\,\pi_{\infty,k}}{(\gamma_k - 1)\,\rho_k}\f]

**Total energy relation:**

\f[\rho E = \Gamma\,p + \Pi_\infty + \frac{1}{2}\rho\,|\mathbf{u}|^2 + q_v\f]

where MFC internally tracks the transformed thermodynamic quantities:

\f[\Gamma_k = \frac{1}{\gamma_k - 1}, \qquad \Pi_{\infty,k} = \frac{\gamma_k\,\pi_{\infty,k}}{\gamma_k - 1}\f]

and the mixture rules are arithmetic averages of these transformed quantities:

\f[\Gamma = \sum_i \frac{\alpha_i}{\gamma_i - 1}, \qquad \Pi_\infty = \sum_i \frac{\alpha_i\,\gamma_i\,\pi_{\infty,i}}{\gamma_i - 1}, \qquad q_v = \sum_i \alpha_i\,\rho_i\,q_{v,i}\f]

The pressure is recovered from the total energy as:

\f[p = \frac{\rho E - \frac{1}{2}\rho\,|\mathbf{u}|^2 - \Pi_\infty - q_v}{\Gamma}\f]

**Phasic speed of sound:**

\f[c_k = \sqrt{\frac{\gamma_k\,(p + \pi_{\infty,k})}{\rho_k}}\f]

**Wood mixture sound speed:**

\f[\frac{1}{\rho\,c^2} = \sum_k \frac{\alpha_k}{\rho_k\,c_k^2}\f]

Input parameters per fluid: `gamma` (\f$\Gamma_k = 1/(\gamma_k - 1)\f$), `pi_inf` (\f$\Pi_{\infty,k} = \gamma_k\,\pi_{\infty,k}/(\gamma_k - 1)\f$), `cv` (\f$c_{v,k}\f$), `qv` (\f$q_{v,k}\f$), `qvp` (\f$q'_{v,k}\f$). Note that `gamma` and `pi_inf` are stored in transformed form, not as the raw physical values (see Section 1b).

### 3.2 Ideal Gas EOS (Chemistry, `chemistry = .true.`)

For reacting gas mixtures:

\f[p = \frac{\rho\,R_u\,T}{W}, \qquad W = \left(\sum_m \frac{Y_m}{W_m}\right)^{-1}\f]

Temperature is obtained from the internal energy by Newton iteration:

\f[e_g - \sum_m e_m(T)\,Y_m = 0\f]

**Species internal energy from enthalpy:**

\f[e_m(T) = \frac{\hat{h}_m(T) - R_u\,T}{W_m}\f]

**NASA polynomial enthalpies** (\cite McBride93):

\f[\frac{\hat{h}_m}{R_u\,T} = \frac{C_0}{T} + \sum_{r=1}^{5} \frac{C_r\,T^{r-1}}{r}\f]

---

## 4. Viscous Stress Tensor (`viscous = .true.`)

**Newtonian viscous stress (no bulk viscosity by default):**

\f[\boldsymbol{\tau}^v = 2\,\eta\left(\mathbf{D} - \frac{1}{3}\,\text{tr}(\mathbf{D})\,\mathbf{I}\right)\f]

where the strain rate tensor is:

\f[\mathbf{D} = \frac{1}{2}\bigl(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\bigr)\f]

**With bulk viscosity:**

\f[\tau_{ij} = \mu\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) + \left(\zeta - \frac{2\mu}{3}\right)\delta_{ij}\,\frac{\partial u_k}{\partial x_k}\f]

**Cartesian components:**

\f[\tau_{xx} = \mu\left(2\,\frac{\partial u}{\partial x} - \frac{2}{3}\nabla\cdot\mathbf{u}\right), \qquad \tau_{xy} = \mu\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)\f]

and similarly for all other components. Cylindrical coordinate formulations include additional \f$1/r\f$ terms.

**Viscosity averaging:**

\f[\frac{1}{\text{Re}_\text{mix}} = \sum_j \frac{\alpha_j}{\text{Re}_j}\f]

Input parameters: `Re_inv` (shear and volume Reynolds numbers per fluid).

---

## 5. Cylindrical Coordinates (`cyl_coord = .true.`) (\cite Wilfong26 Sec. 2.3)

Additional geometric source terms appear with \f$1/r\f$ factors in the continuity, momentum, and energy equations. Key modifications:

- **Radial momentum:** extra \f$p/r\f$ and \f$\tau_{\theta\theta}/r\f$ terms
- **Viscous stress:** \f$\tau_{yy}\f$ includes \f$v/r\f$ corrections:

\f[\tau_{yy} = \mu\left(\frac{4}{3}\frac{\partial v}{\partial r} - \frac{2}{3}\frac{\partial u}{\partial x} - \frac{2}{3}\frac{v}{r}\right)\f]

- **Axis singularity:** axis placed at cell boundary with spectral filtering in the azimuthal direction

---

## 6. Sub-Grid Bubble Dynamics (\cite Wilfong26 Sec. 4.1)

### 6.1 Euler-Euler Bubbles (`bubbles_euler = .true.`)

**Source:** `src/simulation/m_bubbles_EE.fpp`, `src/simulation/m_bubbles.fpp`

#### 6.1.1 Method of Classes (\cite Commander89; \cite Ando11)

**Modified mixture pressure:**

\f[p = (1 - \alpha)\,p_l + \alpha\left(\frac{R^3\,p_{bw}}{\bar{R}^3} - \frac{\rho\,R^3\,\dot{R}^2}{\bar{R}^3}\right)\f]

**Modified stiffened gas for the liquid phase:**

\f[\Gamma_l\,p_l + \Pi_{\infty,l} = \frac{1}{1 - \alpha}\left(E - \frac{1}{2}\rho\,|\mathbf{u}|^2\right)\f]

**Bubble wall pressure (polytropic):**

\f[p_{bw} = \left(p_0 + \frac{2\sigma}{R_0}\right)\left(\frac{R_0}{R}\right)^{3\gamma} - \frac{4\mu\,\dot{R}}{R} - \frac{2\sigma}{R}\f]

**Void fraction transport:**

\f[\frac{\partial \alpha}{\partial t} + \mathbf{u} \cdot \nabla \alpha = \frac{3\,\alpha\,\bar{R}^2\,\dot{R}}{\bar{R}^3}\f]

**Number density conservation:**

\f[\frac{\partial n_\text{bub}}{\partial t} + \nabla \cdot (n_\text{bub}\,\mathbf{u}) = 0\f]

where \f$n = \frac{3}{4\pi}\,\frac{\alpha}{\bar{R}^3}\f$.

**Polydispersity** (`polydisperse = .true.`): Log-normal PDF discretized into \f$N_\text{bin}\f$ equilibrium radii with standard deviation `poly_sigma`, integrated via Simpson's rule.

#### 6.1.2 Rayleigh-Plesset (`bubble_model = 3`) (\cite Rayleigh17; \cite Plesset49)

\f[R\,\ddot{R} + \frac{3}{2}\,\dot{R}^2 = \frac{p_{bw} - p_\infty}{\rho_l}\f]

#### 6.1.3 Keller-Miksis (`bubble_model = 2`) (\cite Keller80)

\f[R\,\ddot{R}\left(1 - \frac{\dot{R}}{c}\right) + \frac{3}{2}\,\dot{R}^2\left(1 - \frac{\dot{R}}{3c}\right) = \frac{p_{bw} - p_\infty}{\rho_l}\left(1 + \frac{\dot{R}}{c}\right) + \frac{R\,\dot{p}_{bw}}{\rho_l\,c}\f]

#### 6.1.4 Gilmore (`bubble_model = 1`) (\cite Gilmore52)

Enthalpy-based formulation with compressibility corrections via the Tait EOS:

\f[R\,\ddot{R}\left(1 - \frac{\dot{R}}{C}\right) + \frac{3}{2}\,\dot{R}^2\left(1 - \frac{\dot{R}}{3C}\right) = H\left(1 + \frac{\dot{R}}{C}\right) + \frac{R\,\dot{H}}{C}\left(1 - \frac{\dot{R}}{C}\right)\f]

where the enthalpy difference is:

\f[H = \frac{n_\text{tait}(1 + B)}{n_\text{tait} - 1}\left[\left(\frac{p_{bw}}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}} - \left(\frac{p_\infty}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}}\right]\f]

and the local liquid sound speed:

\f[C = \sqrt{n_\text{tait}(1+B)\left(\frac{p_\infty}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}} + (n_\text{tait} - 1)\,H}\f]

#### 6.1.5 Non-Polytropic Thermal Model (`polytropic = .false.`) (\cite Preston07)

**Internal bubble pressure ODE:**

\f[\dot{p}_b = \frac{3\gamma_b}{R}\left(-\dot{R}\,p_b + R_v\,T_{bw}\,\dot{m}_v + \frac{\gamma_b - 1}{\gamma_b}\,k_{bw}\left.\frac{\partial T}{\partial r}\right|_R\right)\f]

**Vapor mass flux:**

\f[\dot{m}_v = \frac{D\,\rho_{bw}}{1 - \chi_{vw}}\left.\frac{\partial \chi_v}{\partial r}\right|_R\f]

#### 6.1.6 QBMM Moment Transport (`qbmm = .true.`) (\cite Bryngelson20)

**Population balance equation:**

\f[\frac{\partial f}{\partial t} + \frac{\partial (f\,\dot{R})}{\partial R} + \frac{\partial (f\,\ddot{R})}{\partial \dot{R}} = 0\f]

**Moment transport:**

\f[\frac{\partial (n_\text{bub}\,\mu_i)}{\partial t} + \nabla \cdot (n_\text{bub}\,\mu_i\,\mathbf{u}) = n_\text{bub}\,\dot{\mu}_i\f]

where moments \f$\mu_{i_1,i_2} = \int R^{i_1}\,\dot{R}^{i_2}\,f\,dR\,d\dot{R}\f$.

**CHyQMOM inversion** recovers 4 quadrature nodes \f$(w_j, R_j, \dot{R}_j)\f$ from 6 moments via:

\f[\bar{u} = \frac{\mu_{10}}{\mu_{00}}, \quad \bar{v} = \frac{\mu_{01}}{\mu_{00}}, \quad c_{20} = \frac{\mu_{20}}{\mu_{00}} - \bar{u}^2, \quad c_{11} = \frac{\mu_{11}}{\mu_{00}} - \bar{u}\bar{v}, \quad c_{02} = \frac{\mu_{02}}{\mu_{00}} - \bar{v}^2\f]

### 6.2 Euler-Lagrange Bubbles (`bubbles_lagrange = .true.`) (\cite Maeda18)

**Source:** `src/simulation/m_bubbles_EL.fpp`

Volume-averaged carrier flow equations with bubble source terms:

**Continuity:**

\f[\frac{\partial \rho_l}{\partial t} + \nabla \cdot (\rho_l\,\mathbf{u}_l) = \frac{\rho_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right]\f]

**Momentum:**

\f[\frac{\partial (\rho_l\,\mathbf{u}_l)}{\partial t} + \nabla \cdot (\rho_l\,\mathbf{u}_l \otimes \mathbf{u}_l + p\,\mathbf{I} - \boldsymbol{\tau}_l) = \frac{\rho_l\,\mathbf{u}_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right] - \frac{\alpha}{1 - \alpha}\,\nabla \cdot (p\,\mathbf{I} - \boldsymbol{\tau}_l)\f]

**Energy:**

\f[\frac{\partial E_l}{\partial t} + \nabla \cdot \bigl[(E_l + p)\,\mathbf{u}_l - \boldsymbol{\tau}_l \cdot \mathbf{u}_l\bigr] = \frac{E_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right] - \frac{\alpha}{1 - \alpha}\,\nabla \cdot (p\,\mathbf{u}_l - \boldsymbol{\tau}_l \cdot \mathbf{u}_l)\f]

The left-hand side is the standard conservation law for the liquid phase; the right-hand side source terms capture the effect of the bubbles on the host liquid.

**Void fraction via regularization kernel:**

\f[\alpha(\mathbf{x}) = \sum_n V_n\,\delta_\sigma(\mathbf{x} - \mathbf{x}_n)\f]

where \f$\delta_\sigma\f$ is a Gaussian kernel:

\f[\delta_\sigma(\mathbf{r}) = \frac{1}{(2\pi\sigma^2)^{3/2}}\exp\!\left(-\frac{|\mathbf{r}|^2}{2\sigma^2}\right)\f]

with \f$\sigma = \varepsilon_b \max(\Delta x^{1/3}_\text{cell},\;R_\text{bubble})\f$.

Each bubble is tracked individually with Keller-Miksis dynamics and 4th-order adaptive Runge-Kutta time integration.

---

## 7. Fluid-Structure Interaction

### 7.1 Hypoelastic Model (`hypoelasticity = .true.`) (\cite Rodriguez19; \cite Wilfong26 Sec. 4.1.6)

**Source:** `src/simulation/m_hypoelastic.fpp`

**Cauchy stress decomposition:**

\f[\sigma_{ij} = -p\,\delta_{ij} + \tau_{ij}^{(v)} + \tau_{ij}^{(e)}\f]

**Elastic energy contribution to total energy:**

\f[E = e + \frac{|\mathbf{u}|^2}{2} + \frac{\boldsymbol{\tau}^e : \boldsymbol{\tau}^e}{4\,\rho\,G}\f]

**Elastic stress evolution:**

\f[\frac{\partial (\rho\,\boldsymbol{\tau}^e)}{\partial t} + \nabla \cdot (\rho\,\boldsymbol{\tau}^e \otimes \mathbf{u}) = \mathbf{S}^e\f]

**Source term:**

\f[\mathbf{S}^e = \rho\bigl(\mathbf{l} \cdot \boldsymbol{\tau}^e + \boldsymbol{\tau}^e \cdot \mathbf{l}^T - \boldsymbol{\tau}^e\,\text{tr}(\mathbf{D}) + 2G\,\mathbf{D}^d\bigr)\f]

where \f$\mathbf{l} = \nabla \mathbf{u}\f$ is the velocity gradient and \f$\mathbf{D}^d = \mathbf{D} - \frac{1}{3}\text{tr}(\mathbf{D})\,\mathbf{I}\f$ is the deviatoric strain rate.

**Lie objective temporal derivative (Kelvin-Voigt):**

\f[\hat{\boldsymbol{\tau}}^e = \frac{D\boldsymbol{\tau}^e}{Dt} - \mathbf{l} \cdot \boldsymbol{\tau}^e - \boldsymbol{\tau}^e \cdot \mathbf{l}^T + \boldsymbol{\tau}^e\,\text{tr}(\mathbf{D}) = 2G\,\mathbf{D}^d\f]

This adds 6 additional transport equations in 3D (symmetric stress tensor: \f$\tau_{xx}^e, \tau_{xy}^e, \tau_{yy}^e, \tau_{xz}^e, \tau_{yz}^e, \tau_{zz}^e\f$).

### 7.2 Hyperelastic Model (`hyperelasticity = .true.`) (\cite Kamrin12; \cite Wilfong26 Sec. 4.1.6)

**Source:** `src/simulation/m_hyperelastic.fpp`

**Reference map evolution:**

\f[\frac{\partial (\rho\,\boldsymbol{\xi})}{\partial t} + \nabla \cdot (\rho\,\boldsymbol{\xi} \otimes \mathbf{u}) = 0\f]

**Deformation gradient from reference map:**

\f[\mathbf{F} = (\nabla \boldsymbol{\xi})^{-1}\f]

**Left Cauchy-Green tensor:**

\f[\mathbf{b} = \mathbf{F}\,\mathbf{F}^T\f]

**Neo-Hookean Cauchy stress:**

\f[\boldsymbol{\tau}^e = \frac{G}{J}\left(\mathbf{b} - \frac{\text{tr}(\mathbf{b})}{3}\,\mathbf{I}\right)\f]

where \f$J = \det(\mathbf{F})\f$.

**Hyperelastic energy:**

\f[e^e = \frac{G}{2}\bigl(I_{\mathbf{b}} - 3\bigr), \qquad I_{\mathbf{b}} = \text{tr}(\mathbf{b})\f]

---

## 8. Phase Change (`relax = .true.`) (\cite Wilfong26 Sec. 4.1.3)

**Source:** `src/common/m_phase_change.fpp`

### 8.1 pT-Relaxation (`relax_model = 5`) (\cite Saurel08)

\f$N\f$-fluid pressure-temperature equilibrium. The equilibrium condition is:

\f[f(p) = \sum_i \alpha_i - 1 = 0\f]

**Temperature from energy conservation:**

\f[T = \frac{\rho e + p - \sum_i (\alpha_i \rho_i)\,q_{v,i}}{\sum_i (\alpha_i \rho_i)\,c_{v,i}\,\gamma_i}\f]

**Newton residual:**

\f[g(p) = \sum_i \frac{(\gamma_i - 1)\,(\alpha_i \rho_i)\,c_{v,i}}{(p + \pi_{\infty,i})} \cdot \frac{\rho e + p - \sum_j (\alpha_j \rho_j)\,q_{v,j}}{\sum_j (\alpha_j \rho_j)\,c_{v,j}\,\gamma_j}\f]

Solved via Newton's method for the equilibrium pressure.

### 8.2 pTg-Relaxation (`relax_model = 6`) (\cite Zein10)

Two coupled equations for \f$(\alpha_1 \rho_1,\;p)\f$:

**Gibbs free energy equilibrium (Clausius-Clapeyron):**

\f[F_1 = T\left[(c_{v,l}\gamma_l - c_{v,v}\gamma_v)(1 - \ln T) - (q'_l - q'_v) + c_{v,l}(\gamma_l - 1)\ln(p + \pi_{\infty,l}) - c_{v,v}(\gamma_v - 1)\ln(p + \pi_{\infty,v})\right] + q_{v,l} - q_{v,v} = 0\f]

**Energy conservation constraint:**

\f[F_2 = \rho e + p + m_l\,(q_{v,v} - q_{v,l}) - m_T\,q_{v,v} - m_{qD} + \frac{m_l\,(c_{v,v}\,\gamma_v - c_{v,l}\,\gamma_l) - m_T\,c_{v,v}\,\gamma_v - m_{cpD}}{m_l\left(\frac{c_{v,l}\,(\gamma_l - 1)}{p + \pi_{\infty,l}} - \frac{c_{v,v}\,(\gamma_v - 1)}{p + \pi_{\infty,v}}\right) + \frac{m_T\,c_{v,v}\,(\gamma_v - 1)}{p + \pi_{\infty,v}} + m_{cvgp}} = 0\f]

where \f$m_T\f$ is the total mass, \f$m_l = \alpha_l \rho_l\f$ is the liquid partial density, and \f$m_{qD}\f$, \f$m_{cpD}\f$, \f$m_{cvgp}\f$ are auxiliary thermodynamic sums over additional fluids (beyond the phase-changing pair).

Solved via 2D Newton-Raphson.

---

## 9. Chemistry and Combustion (`chemistry = .true.`) (\cite Wilfong26 Sec. 4.1.7)

**Source:** `src/common/m_chemistry.fpp`

**Species transport:**

\f[\frac{\partial (\rho_g\,Y_m)}{\partial t} + \frac{\partial (\rho_g\,u_i\,Y_m)}{\partial x_i} = W_m\,\dot{\omega}_m\f]

**Net production rate:**

\f[\dot{\omega}_m = \sum_n (\nu''_{mn} - \nu'_{mn})\,\mathcal{R}_n\f]

**Reaction rate (law of mass action):**

\f[\mathcal{R}_n = k_n(T)\left[\prod_j \left(\frac{\rho_g\,Y_j}{W_j}\right)^{\nu'_{jn}} - \frac{1}{K_n}\prod_k \left(\frac{\rho_g\,Y_k}{W_k}\right)^{\nu''_{kn}}\right]\f]

**Arrhenius rate:**

\f[k_n(T) = A_n\,T^{b_n}\exp\!\left(-\frac{T_{a,n}}{T}\right)\f]

**Molecular diffusion** (`transport_model`):
- **Mixture-average:** Species-specific diffusion coefficients \f$D_m^\text{mix}\f$, mass flux: \f$\dot{m}_k = \rho\,D_k^\text{mix}\,(W_k / W_\text{mix})\,\partial X_k / \partial x\f$
- **Unity Lewis number:** \f$D_m = \lambda / (\rho\,c_p)\f$

Enthalpy flux with diffusion:

\f[q_\text{diff} = \lambda\,\frac{\partial T}{\partial x} + \sum_k h_k\,\dot{m}_k\f]

Reaction mechanisms are code-generated via Pyrometheus (\cite Cisneros26), which provides symbolic abstractions for thermochemistry that enable portable GPU computation and automatic differentiation of chemical source terms.

---

## 10. Surface Tension (`surface_tension = .true.`) (\cite Schmidmayer17; \cite Wilfong26 Sec. 4.1.8)

**Source:** `src/simulation/m_surface_tension.fpp`, `src/simulation/include/inline_capillary.fpp`

**Color function advection:**

\f[\frac{\partial c}{\partial t} + \mathbf{u} \cdot \nabla c = 0\f]

**Capillary stress tensor (CSF model):**

\f[\boldsymbol{\Omega} = -\sigma\left(\|\nabla c\|\,\mathbf{I} - \frac{\nabla c \otimes \nabla c}{\|\nabla c\|}\right)\f]

In component form, with \f$\hat{w}_i = (\partial c / \partial x_i) / \|\nabla c\|\f$:

\f[\Omega_{xx} = -\sigma\,(\hat{w}_y^2 + \hat{w}_z^2)\,\|\nabla c\|, \qquad \Omega_{xy} = \sigma\,\hat{w}_x\,\hat{w}_y\,\|\nabla c\|\f]

The capillary stress divergence is added to the momentum and energy equations. The total energy equation becomes:

\f[\frac{\partial (\rho E + \varepsilon_0)}{\partial t} + \nabla \cdot \bigl[(\rho E + \varepsilon_0 + p)\,\mathbf{u} + (\boldsymbol{\Omega} - \boldsymbol{\tau}^v) \cdot \mathbf{u}\bigr] = 0\f]

**Capillary mixture energy:**

\f[\varepsilon_0 = \sigma\,\|\nabla c\|\f]

---

## 11. Magnetohydrodynamics

### 11.1 Ideal MHD (`mhd = .true.`) (\cite Wilfong26 Sec. 4.1.9)

**Continuity:**

\f[\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho\,\mathbf{u}) = 0\f]

**Momentum:**

\f[\frac{\partial (\rho\,\mathbf{u})}{\partial t} + \nabla \cdot \left[\rho\,\mathbf{u} \otimes \mathbf{u} + \left(p + \frac{|\mathbf{B}|^2}{2}\right)\mathbf{I} - \mathbf{B} \otimes \mathbf{B}\right] = 0\f]

**Energy:**

\f[\frac{\partial \mathcal{E}}{\partial t} + \nabla \cdot \left[\left(\mathcal{E} + p + \frac{|\mathbf{B}|^2}{2}\right)\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}\right] = 0\f]

**Induction:**

\f[\frac{\partial \mathbf{B}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{B} - \mathbf{B} \otimes \mathbf{u}) = 0\f]

**Total energy:**

\f[\mathcal{E} = \rho\,e + \frac{1}{2}\rho\,|\mathbf{u}|^2 + \frac{|\mathbf{B}|^2}{2}\f]

**Fast magnetosonic speed:**

\f[c_f = \sqrt{\frac{1}{2}\left(c_s^2 + v_A^2 + \sqrt{(c_s^2 + v_A^2)^2 - 4\,c_s^2\,v_A^2\cos^2\theta}\right)}\f]

**Alfven speed:**

\f[v_A = \sqrt{\frac{|\mathbf{B}|^2}{\rho}}\f]

Uses the HLLD Riemann solver (`riemann_solver = 4`). Hyperbolic divergence cleaning (`hyper_cleaning = .true.`) via the GLM method (\cite Dedner02).

### 11.2 Relativistic MHD (`relativity = .true.`) (\cite Wilfong26 Sec. 4.1.10)

**Conserved variables:**

\f[\mathbf{U} = (D,\;\mathbf{m},\;\tau,\;\mathbf{B})^T\f]

where:

\f[D = \Gamma\,\rho, \qquad \mathbf{m} = \Gamma^2\rho h\,\mathbf{u} + |\mathbf{B}|^2\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}\f]

\f[\tau = \Gamma^2\rho h - p + \frac{|\mathbf{B}|^2}{2} + \frac{|\mathbf{u}|^2|\mathbf{B}|^2 - (\mathbf{B} \cdot \mathbf{u})^2}{2} - \Gamma\,\rho\f]

Primitive recovery uses Newton-Raphson on the nonlinear conserved-to-primitive relation.

---

## 12. Information Geometric Regularization (`igr = .true.`) (\cite Wilfong25a)

**Source:** `src/simulation/m_igr.fpp`

**Modified momentum with entropic pressure** \f$\Sigma\f$**:**

\f[\frac{\partial (\rho\,\mathbf{u})}{\partial t} + \nabla \cdot \bigl[\rho\,\mathbf{u} \otimes \mathbf{u} + (p + \Sigma)\,\mathbf{I} - \boldsymbol{\tau}\bigr] = 0\f]

**Elliptic PDE for** \f$\Sigma\f$**:**

\f[\alpha\left[\text{tr}(\nabla \mathbf{u})^2 + \text{tr}^2(\nabla \mathbf{u})\right] = \frac{\Sigma}{\rho} - \alpha\,\nabla \cdot \left(\frac{\nabla \Sigma}{\rho}\right)\f]

where \f$\alpha \sim \Delta x^2\f$ (regularization strength proportional to mesh spacing squared):

\f[\alpha_\text{IGR} = \alpha_\text{factor} \cdot \max(\Delta x,\;\Delta y,\;\Delta z)^2\f]

**RHS strain-rate source (3D):**

\f[\text{RHS} = \alpha\left[2\left(\frac{\partial u}{\partial y}\frac{\partial v}{\partial x} + \frac{\partial u}{\partial z}\frac{\partial w}{\partial x} + \frac{\partial v}{\partial z}\frac{\partial w}{\partial y}\right) + \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial w}{\partial z}\right)^2 + (\nabla \cdot \mathbf{u})^2\right]\f]

**Iterative solver:** Jacobi (`igr_iter_solver = 1`) or Gauss-Seidel (`igr_iter_solver = 2`), up to `num_igr_iters` iterations (default 5).

Uses Lax-Friedrichs flux (replaces WENO + Riemann solver).

---

## 13. Body Forces (`bf_x`, `bf_y`, `bf_z`)

**Source:** `src/simulation/m_body_forces.fpp`

**Time-dependent acceleration:**

\f[a_i(t) = g_i + k_i\sin(\omega_i\,t - \phi_i)\f]

**Momentum source:**

\f[\frac{\partial (\rho\,u_i)}{\partial t}\bigg|_\text{bf} = \rho\,a_i(t)\f]

**Energy source:**

\f[\frac{\partial E}{\partial t}\bigg|_\text{bf} = \rho\,\mathbf{u} \cdot \mathbf{a}(t)\f]

---

## 14. Acoustic Sources (`acoustic_source = .true.`)

**Source:** `src/simulation/m_acoustic_src.fpp`

Source terms added to the RHS of the governing equations.

Formulation follows \cite Maeda17.

**Discrete delta function (spatial support):**

\f[\delta_h(r) = \frac{1}{(2\pi\sigma^2)^{d/2}}\exp\!\left(-\frac{r^2}{2\sigma^2}\right)\f]

**Forcing form (added to conservative variables):**

\f[\mathbf{s}_\text{ac} = \Omega_\Gamma\,f(t)\left[\frac{1}{c_0},\;\cos\theta,\;\sin\theta,\;\frac{c_0^2}{\gamma - 1}\right]^T\f]

**Temporal profiles:**
- **Sine** (`pulse = 1`): \f$S(t) = M\sin(\omega(t - t_\text{delay}))\f$
- **Gaussian** (`pulse = 2`): \f$S(t) = M\exp\!\bigl(-\frac{1}{2}((t - t_\text{delay})/\sigma_t)^2\bigr)\f$
- **Square** (`pulse = 3`): \f$S(t) = M\,\text{sign}(\sin(\omega(t - t_\text{delay})))\f$
- **Broadband** (`pulse = 4`): superposition of multiple frequencies across a bandwidth

**Spatial supports:** planar, spherical transducer, cylindrical transducer, transducer array (arcuate, annular, circular).

---

## 15. Numerical Methods

### 15.1 Spatial Reconstruction

#### WENO (`weno_order = 3, 5, 7`)

**Source:** `src/simulation/m_weno.fpp`

Weighted sum of candidate polynomials at cell interfaces:

\f[f_{i+1/2} = \sum_r \omega_r\,f_{i+1/2}^{(r)}\f]

**WENO-JS** (\cite Jiang96, default):

\f[\alpha_r = \frac{d_r}{(\beta_r + \varepsilon)^2}, \qquad \omega_r = \frac{\alpha_r}{\sum_s \alpha_s}\f]

where \f$d_r\f$ are ideal weights, \f$\beta_r\f$ are smoothness indicators, and \f$\varepsilon\f$ is a small regularization parameter (`weno_eps`).

**WENO-M** (`mapped_weno = .true.`): \cite Henrick05 mapped weights for improved accuracy at critical points:

\f[\omega_M^{(r)} = \frac{d_r\bigl(1 + d_r - 3\omega_0^{(r)} + (\omega_0^{(r)})^2\bigr)\,\omega_0^{(r)}}{d_r^2 + \omega_0^{(r)}(1 - 2d_r)}, \qquad \omega^{(r)} = \frac{\omega_M^{(r)}}{\sum_s \omega_M^{(s)}}\f]

**WENO-Z** (`wenoz = .true.`): \cite Borges08 improved weights with global smoothness measure:

\f[\alpha_r = d_r\left(1 + \left(\frac{\tau}{\beta_r + \varepsilon}\right)^q\right), \qquad \tau = |\beta_0 - \beta_{k-1}|\f]

The parameter \f$q\f$ controls the convergence rate at critical points (typically \f$q = 1\f$ for fifth-order reconstruction, as used in MFC).

**TENO** (`teno = .true.`): \cite Fu16 targeted ENO with smoothness threshold \f$C_T\f$ (`teno_CT`):

\f[\gamma_r = 1 + \frac{\tau}{\beta_r}, \qquad \xi_r = \frac{\gamma_r}{\sum_s \gamma_s}\f]

If \f$\xi_r < C_T\f$, set \f$\alpha_r = 0\f$ (stencil excluded).

Primitive variable reconstruction is used to avoid spurious oscillations at interfaces.

#### MUSCL (`muscl_order = 2`)

**Source:** `src/simulation/m_muscl.fpp`

\f[q_L(j) = q(j) - \frac{1}{2}\,\phi(r)\,\Delta q, \qquad q_R(j) = q(j) + \frac{1}{2}\,\phi(r)\,\Delta q\f]

**Five slope limiters** (`muscl_lim`):
1. **Minmod:** \f$\phi = \text{sign}(a)\,\min(|a|,\;|b|)\f$ if \f$ab > 0\f$, else \f$0\f$
2. **MC (Monotonized Central):** \f$\phi = \text{sign}(a)\,\min(2|a|,\;2|b|,\;\frac{1}{2}(|a|+|b|))\f$ if \f$ab > 0\f$
3. **OSPRE (Van Albada):** \f$\phi = (a^2 b + a b^2)/(a^2 + b^2)\f$
4. **Van Leer:** \f$\phi = 2ab/(a + b)\f$ if \f$ab > 0\f$
5. **Superbee:** \f$\phi = \max(\min(2|a|,|b|),\;\min(|a|,2|b|))\f$ if \f$ab > 0\f$

where \f$a\f$ and \f$b\f$ are the left and right slope differences.

**THINC interface compression** (`int_comp = .true.`): applies hyperbolic tangent compression near material interfaces:

\f[q_\text{THINC} = q_\min + \frac{q_\max}{2}\left(1 + \text{sign}(s)\,\frac{\tanh(\beta) + A}{1 + A\,\tanh(\beta)}\right)\f]

where \f$A = \frac{\exp(\text{sign}(s)\,\beta\,(2C - 1))\,/\,\cosh(\beta) - 1}{\tanh(\beta)}\f$ and \f$\beta\f$ controls compression steepness.

#### IGR Reconstruction

5th-order or 3rd-order polynomial interpolation without WENO nonlinear weights, using Lax-Friedrichs numerical flux. Stencil coefficients:

**5th-order:** \f$\text{coeff}_L = [-3/60,\;27/60,\;47/60,\;-13/60,\;2/60]\f$

**3rd-order:** \f$\text{coeff}_L = [2/6,\;5/6,\;-1/6]\f$

### 15.2 Riemann Solvers

**Source:** `src/simulation/m_riemann_solvers.fpp`

#### HLL (`riemann_solver = 1`) (\cite Harten83)

\f[\mathbf{F}^\text{HLL} = \frac{S_R\,\mathbf{F}_L - S_L\,\mathbf{F}_R + S_L\,S_R\,(\mathbf{q}_R - \mathbf{q}_L)}{S_R - S_L}\f]

with wave speed estimates \f$S_L = \min(u_L - c_L,\;u_R - c_R)\f$, \f$S_R = \max(u_R + c_R,\;u_L + c_L)\f$.

#### HLLC (`riemann_solver = 2`) (\cite Toro94)

Four-state solver resolving the contact discontinuity. Star-state satisfies:

\f[p_L^* = p_R^* = p^*, \qquad u_L^* = u_R^* = u^*\f]

#### Exact (`riemann_solver = 3`)

Iterative exact Riemann solver.

#### HLLD (`riemann_solver = 4`, MHD only)

Seven-state solver for ideal MHD resolving fast magnetosonic, Alfven, and contact waves (\cite Miyoshi05). The Riemann fan is divided by outer wave speeds \f$S_L\f$, \f$S_R\f$, Alfven speeds \f$S_L^*\f$, \f$S_R^*\f$, and a middle contact \f$S_M\f$:

\f[S_M = \frac{(S_R - u_R)\rho_R u_R - (S_L - u_L)\rho_L u_L - p_{T,R} + p_{T,L}}{(S_R - u_R)\rho_R - (S_L - u_L)\rho_L}\f]

\f[S_L^* = S_M - \frac{|B_x|}{\sqrt{\rho_L^*}}, \qquad S_R^* = S_M + \frac{|B_x|}{\sqrt{\rho_R^*}}\f]

where \f$p_T = p + |\mathbf{B}|^2/2\f$ is the total (thermal + magnetic) pressure. Continuity of normal velocity and total pressure is enforced across the Riemann fan.

### 15.3 Time Integration

**Source:** `src/simulation/m_time_steppers.fpp`

#### TVD Runge-Kutta (`time_stepper = 1, 2, 3`) (\cite Gottlieb98)

**RK1 (Forward Euler):**

\f[\mathbf{q}^{n+1} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)\f]

**RK2:**

\f[\mathbf{q}^{(1)} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)\f]

\f[\mathbf{q}^{n+1} = \frac{1}{2}\mathbf{q}^n + \frac{1}{2}\mathbf{q}^{(1)} + \frac{1}{2}\Delta t\,\mathbf{L}(\mathbf{q}^{(1)})\f]

**RK3 (SSP):**

\f[\mathbf{q}^{(1)} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)\f]

\f[\mathbf{q}^{(2)} = \frac{3}{4}\mathbf{q}^n + \frac{1}{4}\mathbf{q}^{(1)} + \frac{1}{4}\Delta t\,\mathbf{L}(\mathbf{q}^{(1)})\f]

\f[\mathbf{q}^{n+1} = \frac{1}{3}\mathbf{q}^n + \frac{2}{3}\mathbf{q}^{(2)} + \frac{2}{3}\Delta t\,\mathbf{L}(\mathbf{q}^{(2)})\f]

#### Adaptive Time Stepping (`adap_dt = .true.`)

Embedded RK pairs for error estimation with Hairer-Norsett-Wanner algorithm for initial step size.

#### Strang Splitting (\cite Strang68) (for stiff bubble dynamics)

For equations of the form \f$\partial \mathbf{q}/\partial t = -\nabla \cdot \mathbf{F}(\mathbf{q}) + \mathbf{s}(\mathbf{q})\f$, the Strang splitting scheme integrates three sub-equations per time step:

\f[\frac{\partial \mathbf{q}'}{\partial t} = \mathbf{s}(\mathbf{q}'), \quad t \in [0,\,\Delta t/2], \quad \mathbf{q}'(0) = \mathbf{q}^n\f]

\f[\frac{\partial \mathbf{q}''}{\partial t} = -\nabla \cdot \mathbf{F}(\mathbf{q}''), \quad t \in [0,\,\Delta t], \quad \mathbf{q}''(0) = \mathbf{q}'(\Delta t/2)\f]

\f[\frac{\partial \mathbf{q}'''}{\partial t} = \mathbf{s}(\mathbf{q}'''), \quad t \in [0,\,\Delta t/2], \quad \mathbf{q}^{n+1} = \mathbf{q}'''(\Delta t/2)\f]

The stiff bubble source terms are integrated using an adaptive embedded Runge-Kutta scheme with error control.

### 15.4 CFL Condition

**Inviscid:**

\f[\Delta t = \text{CFL} \cdot \min_{i,j,k}\left(\frac{\Delta x_\text{cell}}{|u| + |v| + |w| + c}\right)\f]

**Viscous:**

\f[\Delta t_v = \text{CFL} \cdot \min\left(\frac{\Delta x^2}{\nu}\right)\f]

where \f$c\f$ is the mixture sound speed and \f$\nu\f$ is the kinematic viscosity.

### 15.5 Finite Differences (`fd_order = 1, 2, 4`)

**Source:** `src/common/m_finite_differences.fpp`

Used for viscous fluxes and velocity gradients.

**1st-order:**

\f[f'_i = \frac{f_{i+1} - f_i}{x_{i+1} - x_i}\f]

**2nd-order centered:**

\f[f'_i = \frac{f_{i+1} - f_{i-1}}{x_{i+1} - x_{i-1}}\f]

**4th-order centered:**

\f[f'_i = \frac{-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}}{-x_{i+2} + 8x_{i+1} - 8x_{i-1} + x_{i-2}}\f]

**Boundary-adjusted one-sided stencils (3rd-order):**

\f[f'_i\big|_\text{left} = \frac{-3f_i + 4f_{i+1} - f_{i+2}}{x_{i+2} - x_i}\f]

\f[f'_i\big|_\text{right} = \frac{3f_i - 4f_{i-1} + f_{i-2}}{x_i - x_{i-2}}\f]

---

## 16. Boundary Conditions

**Source:** `src/simulation/m_cbc.fpp`, `src/simulation/m_compute_cbc.fpp`, `src/common/m_constants.fpp`

### 16.1 Simple BCs

| BC Code | Type |
|---------|------|
| `-1`  | Periodic |
| `-2`  | Reflective |
| `-3`  | Ghost cell extrapolation |
| `-4`  | Riemann extrapolation |
| `-14` | Axis symmetry |
| `-15` | Slip wall |
| `-16` | No-slip wall |

### 16.2 Characteristic BCs (\cite Thompson87, \cite Thompson90; `bc_x%%beg = -5` to `-12`)

**Characteristic decomposition:**

\f[\frac{\partial \mathbf{q}_c}{\partial t} + \mathbf{R}_x\,\boldsymbol{\Lambda}_x\,\mathbf{L}_x\,\frac{\partial \mathbf{q}_c}{\partial x} = 0\f]

**Characteristic amplitudes** \f$\mathcal{L}_1\f$ through \f$\mathcal{L}_5\f$ define wave interactions at boundaries:

\f[\mathcal{L}_1 = (u - c)\left(\frac{\partial p}{\partial x} - \rho\,c\,\frac{\partial u}{\partial x}\right), \qquad \mathcal{L}_2 = u\left(c^2\,\frac{\partial \rho}{\partial x} - \frac{\partial p}{\partial x}\right)\f]

\f[\mathcal{L}_3 = u\,\frac{\partial v}{\partial x}, \qquad \mathcal{L}_4 = u\,\frac{\partial w}{\partial x}, \qquad \mathcal{L}_5 = (u + c)\left(\frac{\partial p}{\partial x} + \rho\,c\,\frac{\partial u}{\partial x}\right)\f]

For non-reflecting boundaries, the incoming wave amplitudes are set to zero.

**GRCBC (incoming wave from ghost point):**

\f[\mathcal{L} = -\boldsymbol{\Lambda}^o\,\mathbf{L}_x\,\frac{\partial \mathbf{q}_c}{\partial x} + n_x\,\lambda^i\,\mathbf{L}_x\,\frac{\mathbf{P}\,(\mathbf{q}_p^{(g)} - \mathbf{q}_p)}{\Delta}\f]

**8 types:** slip wall (`-5`), non-reflecting buffer (`-6`), non-reflecting sub. inflow (`-7`), non-reflecting sub. outflow (`-8`), force-free sub. outflow (`-9`), constant-pressure sub. outflow (`-10`), supersonic inflow (`-11`), supersonic outflow (`-12`).

### 16.3 Immersed Boundary Method (`ib = .true.`) (\cite Tseng03; \cite Mittal05; \cite Wilfong26 Sec. 4.1.1)

**Source:** `src/simulation/m_ibm.fpp`

Ghost-cell IBM with level set field \f$\phi\f$.

**Image point:**

\f[\mathbf{x}_{ip} = \mathbf{x}_{gp} + 2\,\phi(\mathbf{x}_{gp})\,\hat{\mathbf{n}}\f]

**Velocity boundary conditions:**
- **No-slip:** \f$\mathbf{u}_{gp} = \mathbf{0}\f$
- **Slip:** \f$\mathbf{u}_{gp} = \mathbf{u}_{ip} - (\hat{\mathbf{n}} \cdot \mathbf{u}_{ip})\,\hat{\mathbf{n}}\f$

**Neumann BC** for pressure, density, and volume fractions.

Supports STL/OBJ geometry import with ray tracing for inside/outside determination.

---

## 17. Low Mach Number Corrections (\cite Wilfong26 Sec. 4.2.4)

**Chen correction** (`low_Mach = 1`, \cite Chen22): anti-dissipation pressure correction (APC) added to the HLLC flux:

\f[p_d = \frac{\rho_L\,\rho_R\,(S_R - u_R)(S_L - u_L)}{\rho_R(S_R - u_R) - \rho_L(S_L - u_L)}\f]

\f[\text{APC}_{p\mathbf{u}} = (z - 1)\,p_d\,\hat{\mathbf{n}}, \qquad \text{APC}_{pE} = (z - 1)\,p_d\,S_*\f]

The corrected HLLC flux in the star region becomes \f$\mathbf{F}^{*} + \text{APC}\f$.

**Thornber correction** (`low_Mach = 2`, \cite Thornber08): modifies the reconstructed velocities at cell interfaces:

\f[u'_L = \frac{u_L + u_R}{2} + z\,\frac{u_L - u_R}{2}, \qquad u'_R = \frac{u_L + u_R}{2} - z\,\frac{u_L - u_R}{2}\f]

\f[z = \min\!\left(\max\!\left(\frac{|u_L|}{c_L},\;\frac{|u_R|}{c_R}\right),\;1\right)\f]

This reduces numerical dissipation at low Mach numbers while recovering the standard scheme at high Mach numbers.

Additionally, the **artificial Mach number** technique (`pi_fac`) modifies \f$\pi_\infty\f$ to artificially increase the local Mach number.

---

## 18. Flux Limiting

**Volume fraction limiting:** enforces \f$0 \le \alpha_k \le 1\f$ and rescales as:

\f[\alpha_k \leftarrow \frac{\alpha_k}{\sum_j \alpha_j}\f]

**Advective flux limiting** based on local volume fraction gradient \f$\chi\f$ to prevent spurious oscillations near material interfaces.

