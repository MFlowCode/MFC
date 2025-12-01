# MFC Feature Compatibility Guide

> **Quick reference** for understanding which MFC features work together and configuration requirements.

> Auto-generated from validation rules in `case_validator.py`.

## 🚀 Common Configuration Patterns

Start with these proven combinations:


<details open>
<summary><b>💧 Multiphase Flow (Bubbles)</b></summary>

```python
'model_eqns': 2,              # 5-equation model
'num_fluids': 2,              # Two-phase
'bubbles_euler': 'T',         # Ensemble-averaged bubbles
'riemann_solver': 2,          # HLLC
'avg_state': 2,               # Arithmetic average
```
</details>

<details>
<summary><b>⚡ Magnetohydrodynamics (MHD)</b></summary>

```python
'model_eqns': 2,              # 5-equation model
'num_fluids': 1,              # Single component
'mhd': 'T',                   # Enable MHD
'riemann_solver': 1,          # HLL (or 4 for HLLD)
```
</details>

<details>
<summary><b>🌡️ Phase Change</b></summary>

```python
'model_eqns': 3,              # 6-equation model
'num_fluids': 2,              # Two-phase
'relax': 'T',                 # Phase change
'riemann_solver': 2,          # HLLC
```
</details>

<details>
<summary><b>💎 Elastic Materials</b></summary>

```python
'model_eqns': 2,              # 5-equation model
'hypoelasticity': 'T',        # Elastic solids
'riemann_solver': 1,          # HLL
```
</details>

## 📊 Feature Compatibility

What works together:


### Physics Models

| Feature | Requirements | Notes |
|---------|-------------|-------|
| Magnetohydrodynamics (MHD) | — | ✓ General use |
| Surface Tension Model | Specific model | ✓ General use |
| Hypoelasticity | — | ✓ General use |
| Hyperelasticity | — | ✓ General use |
| Phase Change (Relaxation) | — | ✓ General use |
| Viscosity | — | ✓ General use |
| Acoustic Sources | — | ✓ General use |


### Bubble Models

| Feature | Requirements | Notes |
|---------|-------------|-------|
| Euler–Euler Bubble Model | — | ✓ General use |
| Euler–Lagrange Bubble Model | — | ✓ General use |
| Quadrature-Based Moment Method (QBMM) | — | ✓ General use |
| Polydisperse Bubble Dynamics | — | ✓ General use |
| adv_n | — | ✓ General use |


### Numerics

| Feature | Requirements | Notes |
|---------|-------------|-------|
| Riemann Solver | Specific model | ✓ General use |
| WENO Order | — | ✓ General use |
| MUSCL Order | — | ✓ General use |


### Geometry

| Feature | Requirements | Notes |
|---------|-------------|-------|
| Immersed Boundaries | — | ✓ General use |
| Cylindrical Coordinates | — | ✓ General use |

## 🔢 Model Equations

Choose your governing equations:


<details>
<summary><b>Model 1: π-γ (Compressible Euler)</b></summary>

- **Use for:** Single-fluid compressible flow
- **Value:** `model_eqns = 1`
- **Note:** Cannot use `num_fluids`, bubbles, or certain WENO variants
</details>

<details>
<summary><b>Model 2: 5-Equation (Most versatile)</b></summary>

- **Use for:** Multiphase, bubbles, elastic materials, MHD
- **Value:** `model_eqns = 2`
- **Requirements:** Set `num_fluids`
- **Compatible with:** Most physics models
</details>

<details>
<summary><b>Model 3: 6-Equation (Phase change)</b></summary>

- **Use for:** Phase change, cavitation
- **Value:** `model_eqns = 3`
- **Requirements:** `riemann_solver = 2` (HLLC), `avg_state = 2`, `wave_speeds = 1`
- **Note:** Not compatible with bubbles or 3D cylindrical
</details>

<details>
<summary><b>Model 4: 4-Equation (Single component)</b></summary>

- **Use for:** Single-component flows with bubbles
- **Value:** `model_eqns = 4`
- **Requirements:** `num_fluids = 1`, set `rhoref` and `pref`
</details>

## ⚙️ Riemann Solvers

| Solver | `riemann_solver` | Best For | Requirements |
|--------|-----------------|----------|-------------|
| **HLL** | `1` | MHD, elastic materials | — |
| **HLLC** | `2` | Bubbles, phase change, multiphase | `avg_state=2` for bubbles |
| **Exact** | `3` | High accuracy (expensive) | — |
| **HLLD** | `4` | MHD (advanced) | MHD only, no relativity |
| **Lax-Friedrichs** | `5` | Robust fallback | Not with cylindrical+viscous |

## 💧 Bubble Models


<details>
<summary><b>Euler-Euler (`bubbles_euler`)</b></summary>

**Requirements:**
- `model_eqns = 2` or `4`
- `riemann_solver = 2` (HLLC)
- `avg_state = 2`
- Set `nb` (number of bins) ≥ 1

**Extensions:**
- `polydisperse = T`: Multiple bubble sizes (requires odd `nb > 1`)
- `qbmm = T`: Quadrature method (requires `nnode = 4`)
- `adv_n = T`: Number density advection (requires `num_fluids = 1`)
</details>

<details>
<summary><b>Euler-Lagrange (`bubbles_lagrange`)</b></summary>

**Requirements:**
- `n > 0` (2D or 3D only)
- `file_per_process = F`
- Not compatible with `model_eqns = 3`

**Note:** Tracks individual bubbles
</details>

## 📖 Quick Parameter Reference

Key parameters and their constraints:


<details>
<summary><b>MHD</b> (`mhd`)</summary>

**Requirements:**
- relativity requires mhd to be enabled
- Powell's method requires mhd to be enabled

**Incompatibilities:**
- Bx0 must not be set if MHD is not enabled
- Bx0 must not be set in 2D/3D MHD simulations

</details>


<details>
<summary><b>Surface Tension</b> (`surface_tension`)</summary>

**Requirements:**
- sigma must be set if surface_tension is enabled
- The surface tension model requires model_eqns = 2 or model_eqns = 3
- The surface tension model requires num_fluids = 2

**Valid values:**
- sigma must be greater than or equal to zero

</details>


<details>
<summary><b>Number of Fluids</b> (`num_fluids`)</summary>

**Requirements:**
- 5-equation model (model_eqns = 2) requires num_fluids to be set
- 6-equation model (model_eqns = 3) requires num_fluids to be set
- 4-equation model (model_eqns = 4) requires num_fluids to be set

**Incompatibilities:**
- num_fluids is not supported for model_eqns = 1
- num_fluids = 1 does not support mpp_lim

**Valid values:**
- num_fluids must be positive
- perturb_flow_fluid must be between 0 and num_fluids

</details>


<details>
<summary><b>Cylindrical Coordinates</b> (`cyl_coord`)</summary>

**Incompatibilities:**
- 6-equation model (model_eqns = 3) does not support cylindrical coordinates (cyl_coord = T and p != 0)
- Bubble models untested in cylindrical coordinates
- Acoustic source is not supported in 3D cylindrical simulations

**Valid values:**
- p must be odd for cylindrical coordinates

</details>


<details>
<summary><b>Immersed Boundaries</b> (`ib`)</summary>

**Requirements:**
- Immersed Boundaries do not work in 1D (requires n > 0)

**Incompatibilities:**
- output_partial_domain is incompatible with certain output flags

**Valid values:**
- num_ibs must be between 1 and num_patches_max (10)

</details>


---

💡 **Tip:** If you encounter a validation error, check the relevant section above or 
review [`case_validator.py`](../../toolchain/mfc/case_validator.py) for complete validation logic.

