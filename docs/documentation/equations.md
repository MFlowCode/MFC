@page equations Equations

# MFC: Comprehensive Equations Reference

This document catalogs every equation solved by MFC, organized by physical model.
Each section notes the input parameter(s) that activate the corresponding physics module, cross-references the relevant source files, and cites the originating papers.

**References:**
- **[bryngelson-CPC-20]** Bryngelson et al., "MFC: An open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver," *Computer Physics Communications* **266**, 107396 (2021).
- **[radhakrishnan-CPC-24]** Radhakrishnan et al., "Method for portable, scalable, and performant GPU-accelerated simulation of multiphase compressible flow," *Computer Physics Communications* **302**, 109238 (2024).
- **[wilfong-mfc-26]** Wilfong et al., "MFC 5.0: An exascale many-physics flow solver," *Computer Physics Communications* **322**, 110055 (2026).
- **[rodriguez-CAV-21]** Rodriguez et al., "Acoustically-induced bubble growth and phase change dynamics near compliant surfaces," *11th International Symposium on Cavitation* (2021).
- **[cisneros-cpc-25]** Cisneros-Garibay et al., "Pyrometheus: Symbolic abstractions for XPU and automatically differentiated computation of combustion kinetics and thermodynamics," *Computer Physics Communications* **320**, 109987 (2026).
- **[wilfong-igr-SC25]** Wilfong et al., "Information Geometric Regularization of the Barotropic Euler Equations," SC'25 (2025).
- **[wilfong-arxiv-25]** Wilfong et al., arxiv preprint (2025).

---

## 1. Overview

MFC solves the compressible Navier-Stokes equations (or Euler equations when viscosity is off) in a finite volume framework. The general semi-discrete form is:

$$\frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{q}) + \mathbf{h}(\mathbf{q})\,\nabla \cdot \mathbf{u} = \mathbf{s}(\mathbf{q})$$

where:
- $\mathbf{q}$ is the conservative variable vector,
- $\mathbf{F}$ is the flux tensor,
- $\mathbf{h}(\mathbf{q})\,\nabla \cdot \mathbf{u}$ contains non-conservative terms (volume fraction advection),
- $\mathbf{s}(\mathbf{q})$ is the source vector (bubbles, body forces, chemistry, etc.).

The parameter `model_eqns` (1, 2, 3, or 4) selects the governing equation set.

**Key source files:** `src/simulation/m_rhs.fpp` (RHS evaluation), `src/common/m_variables_conversion.fpp` (EOS and variable conversion).

---

## 2. Governing PDEs

### 2.1 Five-Equation Model (`model_eqns = 2`)

The primary workhorse model (Allaire et al., 2002). The state vector is:

$$\mathbf{q} = \bigl(\alpha_1 \rho_1,\;\alpha_2 \rho_2,\;\ldots,\;\rho u_1,\;\rho u_2,\;\rho u_3,\;\rho E,\;\alpha_1,\;\alpha_2,\;\ldots\bigr)^T$$

**Continuity** (one per component):

$$\frac{\partial (\alpha_i \rho_i)}{\partial t} + \nabla \cdot (\alpha_i \rho_i\,\mathbf{u}) = 0$$

**Momentum:**

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot \bigl(\rho\,\mathbf{u} \otimes \mathbf{u} + p\,\mathbf{I} - \boldsymbol{\tau}^v\bigr) = 0$$

**Energy:**

$$\frac{\partial (\rho E)}{\partial t} + \nabla \cdot \bigl[(\rho E + p)\,\mathbf{u} - \boldsymbol{\tau}^v \cdot \mathbf{u}\bigr] = 0$$

**Volume fraction advection:**

$$\frac{\partial \alpha_i}{\partial t} + \mathbf{u} \cdot \nabla \alpha_i = K\,\nabla \cdot \mathbf{u}$$

where the $K$ term enforces interface conditions via the Wood sound speed:

$$K = \frac{\rho_2 c_2^2 - \rho_1 c_1^2}{\dfrac{\rho_1 c_1^2}{\alpha_1} + \dfrac{\rho_2 c_2^2}{\alpha_2}}$$

Setting `alt_soundspeed = .true.` enables the $K$ correction (Kapila model with Wood sound speed). Setting `alt_soundspeed = .false.` uses the Allaire variant without the $K$ correction, which is conservative but does not strictly obey the second law of thermodynamics.

**Mixture rules:**

$$1 = \sum_i \alpha_i, \qquad \rho = \sum_i \alpha_i \rho_i, \qquad \rho e = \sum_i \alpha_i \rho_i e_i$$

**Ref:** [wilfong-mfc-26] Section 2.1.1, [bryngelson-CPC-20] Section 3.1, [radhakrishnan-CPC-24] Section 2.1 Eq. 2.

### 2.2 Six-Equation Model (`model_eqns = 3`)

Allows pressure disequilibrium between phases (Saurel et al., 2009).

**Continuity and momentum:** Same as the five-equation model.

**Separate phasic internal energy:**

$$\frac{\partial (\alpha_i \rho_i e_i)}{\partial t} + \nabla \cdot (\alpha_i \rho_i e_i\,\mathbf{u}) + \alpha_i p_i\,\nabla \cdot \mathbf{u} = -\mu\,p_I\,(p_2 - p_1) - \alpha_i\,\boldsymbol{\tau}_i^v : \nabla \mathbf{u}$$

**Volume fraction:**

$$\frac{\partial \alpha_1}{\partial t} + \mathbf{u} \cdot \nabla \alpha_1 = \mu\,(p_1 - p_2)$$

**Interfacial pressure:**

$$p_I = \frac{z_2\,p_1 + z_1\,p_2}{z_1 + z_2}, \qquad z_i = \rho_i\,c_i$$

Infinite pressure relaxation is applied at each Runge-Kutta stage to drive toward pressure equilibrium.

**Mixture speed of sound:**

$$c^2 = \sum_k Y_k\,c_k^2$$

With phase change (`relax = .true.`), additional source terms appear in the phasic energy and volume fraction equations:
- **Pressure relaxation:** $\mu\,\delta p$ where $\delta p = p_1 - p_2$
- **Thermal transfer:** $Q = \theta\,(T_2 - T_1)$
- **Mass transfer:** $\dot{m} = \nu\,(g_2 - g_1)$ (Gibbs free energy difference)

See [Section 8](#8-phase-change-relax--true) for details.

**Ref:** [wilfong-mfc-26] Section 2.1.2, [bryngelson-CPC-20] Section 3.2, [rodriguez-CAV-21] Section 2.1, `src/simulation/m_pressure_relaxation.fpp`.

### 2.3 Other Model Variants

- `model_eqns = 1`: **Gamma/pi_inf model** — simplified single-fluid formulation using mixture $\gamma$ and $\pi_\infty$ directly without tracking individual volume fractions (Johnsen, 2008).
- `model_eqns = 4`: **Four-equation model** — reduced model from the six-equation system after full pressure-temperature equilibrium relaxation (Tait-like compressible liquid).

---

## 3. Equations of State

### 3.1 Stiffened Gas EOS

The primary closure for each phase:

$$p_k = (\gamma_k - 1)\,\rho_k\,e_k - \gamma_k\,\pi_{\infty,k}$$

Equivalently:

$$e_k = \frac{p_k + \gamma_k\,\pi_{\infty,k}}{(\gamma_k - 1)\,\rho_k}$$

**Total energy relation:**

$$\rho E = \Gamma\,p + \Pi_\infty + \frac{1}{2}\rho\,|\mathbf{u}|^2 + q_v$$

where MFC internally tracks the transformed thermodynamic quantities:

$$\Gamma_k = \frac{1}{\gamma_k - 1}, \qquad \Pi_{\infty,k} = \frac{\gamma_k\,\pi_{\infty,k}}{\gamma_k - 1}$$

and the mixture rules are arithmetic averages of these transformed quantities:

$$\Gamma = \sum_i \frac{\alpha_i}{\gamma_i - 1}, \qquad \Pi_\infty = \sum_i \frac{\alpha_i\,\gamma_i\,\pi_{\infty,i}}{\gamma_i - 1}, \qquad q_v = \sum_i \alpha_i\,\rho_i\,q_{v,i}$$

The pressure is recovered from the total energy as:

$$p = \frac{\rho E - \frac{1}{2}\rho\,|\mathbf{u}|^2 - \Pi_\infty - q_v}{\Gamma}$$

**Phasic speed of sound:**

$$c_k = \sqrt{\frac{\gamma_k\,(p + \pi_{\infty,k})}{\rho_k}}$$

**Wood mixture sound speed:**

$$\frac{1}{\rho\,c^2} = \sum_k \frac{\alpha_k}{\rho_k\,c_k^2}$$

Input parameters per fluid: `gamma` ($\gamma_k$), `pi_inf` ($\pi_{\infty,k}$), `cv` ($c_{v,k}$), `qv` ($q_{v,k}$), `qvp` ($q'_{v,k}$).

**Ref:** [wilfong-mfc-26] Section 2.2, [radhakrishnan-CPC-24] Section 2.1 Eq. 1, `src/common/m_variables_conversion.fpp`.

### 3.2 Ideal Gas EOS (Chemistry, `chemistry = .true.`)

For reacting gas mixtures:

$$p = \frac{\rho\,R_u\,T}{W}, \qquad W = \left(\sum_m \frac{Y_m}{W_m}\right)^{-1}$$

Temperature is obtained from the internal energy by Newton iteration:

$$e_g - \sum_m e_m(T)\,Y_m = 0$$

**Species internal energy from enthalpy:**

$$e_m(T) = \frac{\hat{h}_m(T) - R_u\,T}{W_m}$$

**NASA polynomial enthalpies:**

$$\frac{\hat{h}_m}{R_u\,T} = \frac{C_0}{T} + \sum_{r=1}^{5} \frac{C_r\,T^{r-1}}{r}$$

**Ref:** [wilfong-mfc-26] Section 4.1.7 Eqs. 18-22.

---

## 4. Viscous Stress Tensor (`viscous = .true.`)

**Newtonian viscous stress (no bulk viscosity by default):**

$$\boldsymbol{\tau}^v = 2\,\eta\left(\mathbf{D} - \frac{1}{3}\,\text{tr}(\mathbf{D})\,\mathbf{I}\right)$$

where the strain rate tensor is:

$$\mathbf{D} = \frac{1}{2}\bigl(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\bigr)$$

**With bulk viscosity:**

$$\tau_{ij} = \mu\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) + \left(\zeta - \frac{2\mu}{3}\right)\delta_{ij}\,\frac{\partial u_k}{\partial x_k}$$

**Cartesian components:**

$$\tau_{xx} = \mu\left(2\,\frac{\partial u}{\partial x} - \frac{2}{3}\nabla\cdot\mathbf{u}\right), \qquad \tau_{xy} = \mu\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)$$

and similarly for all other components. Cylindrical coordinate formulations include additional $1/r$ terms.

**Viscosity averaging:**

$$\frac{1}{\text{Re}_\text{mix}} = \sum_j \frac{\alpha_j}{\text{Re}_j}$$

Input parameters: `Re_inv` (shear and volume Reynolds numbers per fluid).

**Ref:** [wilfong-mfc-26] Eqs. 1-2, `src/simulation/m_viscous.fpp`.

---

## 5. Cylindrical Coordinates (`cyl_coord = .true.`)

Additional geometric source terms appear with $1/r$ factors in the continuity, momentum, and energy equations. Key modifications:

- **Radial momentum:** extra $p/r$ and $\tau_{\theta\theta}/r$ terms
- **Viscous stress:** $\tau_{yy}$ includes $v/r$ corrections:

$$\tau_{yy} = \mu\left(\frac{4}{3}\frac{\partial v}{\partial r} - \frac{2}{3}\frac{\partial u}{\partial x} - \frac{2}{3}\frac{v}{r}\right)$$

- **Axis singularity:** axis placed at cell boundary with spectral filtering in the azimuthal direction

**Ref:** [wilfong-mfc-26] Section 2.3, [bryngelson-CPC-20] Section 4.3.

---

## 6. Sub-Grid Bubble Dynamics

### 6.1 Euler-Euler Bubbles (`bubbles_euler = .true.`)

**Source:** `src/simulation/m_bubbles_EE.fpp`, `src/simulation/m_bubbles.fpp`

#### 6.1.1 Method of Classes

**Modified mixture pressure:**

$$p = (1 - \alpha)\,p_l + \alpha\left(\frac{R^3\,p_{bw}}{\bar{R}^3} - \frac{\rho\,R^3\,\dot{R}^2}{\bar{R}^3}\right)$$

**Modified stiffened gas for the liquid phase:**

$$\Gamma_l\,p_l + \Pi_{\infty,l} = \frac{1}{1 - \alpha}\left(E - \frac{1}{2}\rho\,|\mathbf{u}|^2\right)$$

**Bubble wall pressure (polytropic):**

$$p_{bw} = \left(p_0 + \frac{2\sigma}{R_0}\right)\left(\frac{R_0}{R}\right)^{3\gamma} - \frac{4\mu\,\dot{R}}{R} - \frac{2\sigma}{R}$$

**Void fraction transport:**

$$\frac{\partial \alpha}{\partial t} + \mathbf{u} \cdot \nabla \alpha = \frac{3\,\alpha\,\bar{R}^2\,\dot{R}}{\bar{R}^3}$$

**Number density conservation:**

$$\frac{\partial n_\text{bub}}{\partial t} + \nabla \cdot (n_\text{bub}\,\mathbf{u}) = 0$$

where $n = \frac{3}{4\pi}\,\frac{\alpha}{\bar{R}^3}$.

**Polydispersity** (`polydisperse = .true.`): Log-normal PDF discretized into $N_\text{bin}$ equilibrium radii with standard deviation `poly_sigma`, integrated via Simpson's rule.

#### 6.1.2 Rayleigh-Plesset (`bubble_model = 3`)

$$R\,\ddot{R} + \frac{3}{2}\,\dot{R}^2 = \frac{p_{bw} - p_\infty}{\rho_l}$$

#### 6.1.3 Keller-Miksis (`bubble_model = 2`)

$$R\,\ddot{R}\left(1 - \frac{\dot{R}}{c}\right) + \frac{3}{2}\,\dot{R}^2\left(1 - \frac{\dot{R}}{3c}\right) = \frac{p_{bw} - p_\infty}{\rho_l}\left(1 + \frac{\dot{R}}{c}\right) + \frac{R\,\dot{p}_{bw}}{\rho_l\,c}$$

#### 6.1.4 Gilmore (`bubble_model = 1`)

Enthalpy-based formulation with compressibility corrections via the Tait EOS:

$$R\,\ddot{R}\left(1 - \frac{\dot{R}}{C}\right) + \frac{3}{2}\,\dot{R}^2\left(1 - \frac{\dot{R}}{3C}\right) = H\left(1 + \frac{\dot{R}}{C}\right) + \frac{R\,\dot{H}}{C}\left(1 - \frac{\dot{R}}{C}\right)$$

where the enthalpy difference is:

$$H = \frac{n_\text{tait}(1 + B)}{n_\text{tait} - 1}\left[\left(\frac{p_{bw}}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}} - \left(\frac{p_\infty}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}}\right]$$

and the local liquid sound speed:

$$C = \sqrt{n_\text{tait}(1+B)\left(\frac{p_\infty}{1+B} + 1\right)^{(n_\text{tait}-1)/n_\text{tait}} + (n_\text{tait} - 1)\,H}$$

#### 6.1.5 Non-Polytropic Thermal Model (`polytropic = .false.`)

**Internal bubble pressure ODE:**

$$\dot{p}_b = \frac{3\gamma_b}{R}\left(-\dot{R}\,p_b + R_v\,T_{bw}\,\dot{m}_v + \frac{\gamma_b - 1}{\gamma_b}\,k_{bw}\left.\frac{\partial T}{\partial r}\right|_R\right)$$

**Vapor mass flux:**

$$\dot{m}_v = \frac{D\,\rho_{bw}}{1 - \chi_{vw}}\left.\frac{\partial \chi_v}{\partial r}\right|_R$$

#### 6.1.6 QBMM Moment Transport (`qbmm = .true.`)

**Population balance equation:**

$$\frac{\partial f}{\partial t} + \frac{\partial (f\,\dot{R})}{\partial R} + \frac{\partial (f\,\ddot{R})}{\partial \dot{R}} = 0$$

**Moment transport:**

$$\frac{\partial (n_\text{bub}\,\mu_i)}{\partial t} + \nabla \cdot (n_\text{bub}\,\mu_i\,\mathbf{u}) = n_\text{bub}\,\dot{\mu}_i$$

where moments $\mu_{i_1,i_2} = \int R^{i_1}\,\dot{R}^{i_2}\,f\,dR\,d\dot{R}$.

**CHyQMOM inversion** recovers 4 quadrature nodes $(w_j, R_j, \dot{R}_j)$ from 6 moments via:

$$\bar{u} = \frac{\mu_{10}}{\mu_{00}}, \quad \bar{v} = \frac{\mu_{01}}{\mu_{00}}, \quad c_{20} = \frac{\mu_{20}}{\mu_{00}} - \bar{u}^2, \quad c_{11} = \frac{\mu_{11}}{\mu_{00}} - \bar{u}\bar{v}, \quad c_{02} = \frac{\mu_{02}}{\mu_{00}} - \bar{v}^2$$

**Ref:** [wilfong-mfc-26] Section 4.1.4, `src/simulation/m_qbmm.fpp`.

### 6.2 Euler-Lagrange Bubbles (`bubbles_lagrange = .true.`)

**Source:** `src/simulation/m_bubbles_EL.fpp`

Volume-averaged carrier flow equations with bubble source terms (Eq. 12 of [wilfong-mfc-26]):

**Continuity:**

$$\frac{\partial \rho_l}{\partial t} + \nabla \cdot (\rho_l\,\mathbf{u}_l) = \frac{\rho_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right]$$

**Momentum:**

$$\frac{\partial (\rho_l\,\mathbf{u}_l)}{\partial t} + \nabla \cdot (\rho_l\,\mathbf{u}_l \otimes \mathbf{u}_l + p\,\mathbf{I} - \boldsymbol{\tau}_l) = \frac{\rho_l\,\mathbf{u}_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right] - \frac{\alpha}{1 - \alpha}\,\nabla \cdot (p\,\mathbf{I} - \boldsymbol{\tau}_l)$$

**Energy:**

$$\frac{\partial E_l}{\partial t} + \nabla \cdot \bigl[(E_l + p)\,\mathbf{u}_l - \boldsymbol{\tau}_l \cdot \mathbf{u}_l\bigr] = \frac{E_l}{1 - \alpha}\left[\frac{\partial \alpha}{\partial t} + \mathbf{u}_l \cdot \nabla \alpha\right] - \frac{\alpha}{1 - \alpha}\,\nabla \cdot (p\,\mathbf{u}_l - \boldsymbol{\tau}_l \cdot \mathbf{u}_l)$$

The left-hand side is the standard conservation law for the liquid phase; the right-hand side source terms capture the effect of the bubbles on the host liquid.

**Void fraction via regularization kernel:**

$$\alpha(\mathbf{x}) = \sum_n V_n\,\delta_\sigma(\mathbf{x} - \mathbf{x}_n)$$

where $\delta_\sigma$ is a Gaussian kernel:

$$\delta_\sigma(\mathbf{r}) = \frac{1}{(2\pi\sigma^2)^{3/2}}\exp\!\left(-\frac{|\mathbf{r}|^2}{2\sigma^2}\right)$$

with $\sigma = \varepsilon_b \max(\Delta x^{1/3}_\text{cell},\;R_\text{bubble})$.

Each bubble is tracked individually with Keller-Miksis dynamics and 4th-order adaptive Runge-Kutta time integration.

**Ref:** [wilfong-mfc-26] Section 4.1.5 Eqs. 13-14.

---

## 7. Fluid-Structure Interaction

### 7.1 Hypoelastic Model (`hypoelasticity = .true.`)

**Source:** `src/simulation/m_hypoelastic.fpp`

**Cauchy stress decomposition:**

$$\sigma_{ij} = -p\,\delta_{ij} + \tau_{ij}^{(v)} + \tau_{ij}^{(e)}$$

**Elastic energy contribution to total energy:**

$$E = e + \frac{|\mathbf{u}|^2}{2} + \frac{\boldsymbol{\tau}^e : \boldsymbol{\tau}^e}{4\,\rho\,G}$$

**Elastic stress evolution:**

$$\frac{\partial (\rho\,\boldsymbol{\tau}^e)}{\partial t} + \nabla \cdot (\rho\,\boldsymbol{\tau}^e \otimes \mathbf{u}) = \mathbf{S}^e$$

**Source term:**

$$\mathbf{S}^e = \rho\bigl(\mathbf{l} \cdot \boldsymbol{\tau}^e + \boldsymbol{\tau}^e \cdot \mathbf{l}^T - \boldsymbol{\tau}^e\,\text{tr}(\mathbf{D}) + 2G\,\mathbf{D}^d\bigr)$$

where $\mathbf{l} = \nabla \mathbf{u}$ is the velocity gradient and $\mathbf{D}^d = \mathbf{D} - \frac{1}{3}\text{tr}(\mathbf{D})\,\mathbf{I}$ is the deviatoric strain rate.

**Lie objective temporal derivative (Kelvin-Voigt):**

$$\hat{\boldsymbol{\tau}}^e = \frac{D\boldsymbol{\tau}^e}{Dt} - \mathbf{l} \cdot \boldsymbol{\tau}^e - \boldsymbol{\tau}^e \cdot \mathbf{l}^T + \boldsymbol{\tau}^e\,\text{tr}(\mathbf{D}) = 2G\,\mathbf{D}^d$$

This adds 6 additional transport equations in 3D (symmetric stress tensor: $\tau_{xx}^e, \tau_{xy}^e, \tau_{yy}^e, \tau_{xz}^e, \tau_{yz}^e, \tau_{zz}^e$).

**Ref:** [wilfong-mfc-26] Section 4.1.6, [radhakrishnan-CPC-24] Section 2.2 Eqs. 3-5.

### 7.2 Hyperelastic Model (`hyperelasticity = .true.`)

**Source:** `src/simulation/m_hyperelastic.fpp`

**Reference map evolution:**

$$\frac{\partial (\rho\,\boldsymbol{\xi})}{\partial t} + \nabla \cdot (\rho\,\boldsymbol{\xi} \otimes \mathbf{u}) = 0$$

**Deformation gradient from reference map:**

$$\mathbf{F} = (\nabla \boldsymbol{\xi})^{-1}$$

**Left Cauchy-Green tensor:**

$$\mathbf{b} = \mathbf{F}\,\mathbf{F}^T$$

**Neo-Hookean Cauchy stress:**

$$\boldsymbol{\tau}^e = \frac{G}{J}\left(\mathbf{b} - \frac{\text{tr}(\mathbf{b})}{3}\,\mathbf{I}\right)$$

where $J = \det(\mathbf{F})$.

**Hyperelastic energy:**

$$e^e = \frac{G}{2}\bigl(I_{\mathbf{b}} - 3\bigr), \qquad I_{\mathbf{b}} = \text{tr}(\mathbf{b})$$

**Ref:** [wilfong-mfc-26] Section 4.1.6.

---

## 8. Phase Change (`relax = .true.`)

**Source:** `src/common/m_phase_change.fpp`

### 8.1 pT-Relaxation (`relax_model = 5`)

$N$-fluid pressure-temperature equilibrium. The equilibrium condition is:

$$f(p) = \sum_i \alpha_i - 1 = 0$$

**Temperature from energy conservation:**

$$T = \frac{\rho e + p - \sum_i (\alpha_i \rho_i)\,q_{v,i}}{\sum_i (\alpha_i \rho_i)\,c_{v,i}\,\gamma_i}$$

**Newton residual:**

$$g(p) = \sum_i \frac{(\gamma_i - 1)\,(\alpha_i \rho_i)\,c_{v,i}}{(p + \pi_{\infty,i})} \cdot \frac{\rho e + p - \sum_j (\alpha_j \rho_j)\,q_{v,j}}{\sum_j (\alpha_j \rho_j)\,c_{v,j}\,\gamma_j}$$

Solved via Newton's method for the equilibrium pressure.

### 8.2 pTg-Relaxation (`relax_model = 6`)

Two coupled equations for $(\alpha_1 \rho_1,\;p)$:

**Gibbs free energy equilibrium (Clausius-Clapeyron):**

$$F_1 = T\left[(c_{v,l}\gamma_l - c_{v,v}\gamma_v)(1 - \ln T) - (q'_l - q'_v) + c_{v,l}(\gamma_l - 1)\ln(p + \pi_{\infty,l}) - c_{v,v}(\gamma_v - 1)\ln(p + \pi_{\infty,v})\right] + q_{v,l} - q_{v,v} = 0$$

**Energy conservation constraint:**

$$F_2 = \rho e + p + m_l\,(q_{v,v} - q_{v,l}) - m_T\,q_{v,v} - m_{qD} + \frac{m_l\,(c_{v,v}\,\gamma_v - c_{v,l}\,\gamma_l) - m_T\,c_{v,v}\,\gamma_v - m_{cpD}}{m_l\left(\frac{c_{v,l}\,(\gamma_l - 1)}{p + \pi_{\infty,l}} - \frac{c_{v,v}\,(\gamma_v - 1)}{p + \pi_{\infty,v}}\right) + \frac{m_T\,c_{v,v}\,(\gamma_v - 1)}{p + \pi_{\infty,v}} + m_{cvgp}} = 0$$

where $m_T$ is the total mass, $m_l = \alpha_l \rho_l$ is the liquid partial density, and $m_{qD}$, $m_{cpD}$, $m_{cvgp}$ are auxiliary thermodynamic sums over additional fluids (beyond the phase-changing pair).

Solved via 2D Newton-Raphson.

**Ref:** [wilfong-mfc-26] Section 4.1.3, [rodriguez-CAV-21] Section 2.

---

## 9. Chemistry and Combustion (`chemistry = .true.`)

**Source:** `src/common/m_chemistry.fpp`

**Species transport:**

$$\frac{\partial (\rho_g\,Y_m)}{\partial t} + \frac{\partial (\rho_g\,u_i\,Y_m)}{\partial x_i} = W_m\,\dot{\omega}_m$$

**Net production rate:**

$$\dot{\omega}_m = \sum_n (\nu''_{mn} - \nu'_{mn})\,\mathcal{R}_n$$

**Reaction rate (law of mass action):**

$$\mathcal{R}_n = k_n(T)\left[\prod_j \left(\frac{\rho_g\,Y_j}{W_j}\right)^{\nu'_{jn}} - \frac{1}{K_n}\prod_k \left(\frac{\rho_g\,Y_k}{W_k}\right)^{\nu''_{kn}}\right]$$

**Arrhenius rate:**

$$k_n(T) = A_n\,T^{b_n}\exp\!\left(-\frac{T_{a,n}}{T}\right)$$

**Molecular diffusion** (`transport_model`):
- **Mixture-average:** Species-specific diffusion coefficients $D_m^\text{mix}$, mass flux: $\dot{m}_k = \rho\,D_k^\text{mix}\,(W_k / W_\text{mix})\,\partial X_k / \partial x$
- **Unity Lewis number:** $D_m = \lambda / (\rho\,c_p)$

Enthalpy flux with diffusion:

$$q_\text{diff} = \lambda\,\frac{\partial T}{\partial x} + \sum_k h_k\,\dot{m}_k$$

Reaction mechanisms are code-generated via Pyrometheus, which provides symbolic abstractions for thermochemistry that enable portable GPU computation and automatic differentiation of chemical source terms.

**Ref:** [wilfong-mfc-26] Section 4.1.7 Eqs. 15-22, [cisneros-cpc-25] Section 5.

---

## 10. Surface Tension (`surface_tension = .true.`)

**Source:** `src/simulation/m_surface_tension.fpp`, `src/simulation/include/inline_capillary.fpp`

**Color function advection:**

$$\frac{\partial c}{\partial t} + \mathbf{u} \cdot \nabla c = 0$$

**Capillary stress tensor (CSF model):**

$$\boldsymbol{\Omega} = -\sigma\left(\|\nabla c\|\,\mathbf{I} - \frac{\nabla c \otimes \nabla c}{\|\nabla c\|}\right)$$

In component form, with $\hat{w}_i = (\partial c / \partial x_i) / \|\nabla c\|$:

$$\Omega_{xx} = -\sigma\,(\hat{w}_y^2 + \hat{w}_z^2)\,\|\nabla c\|, \qquad \Omega_{xy} = \sigma\,\hat{w}_x\,\hat{w}_y\,\|\nabla c\|$$

The capillary stress divergence is added to the momentum and energy equations. The total energy equation becomes:

$$\frac{\partial (\rho E + \varepsilon_0)}{\partial t} + \nabla \cdot \bigl[(\rho E + \varepsilon_0 + p)\,\mathbf{u} + (\boldsymbol{\Omega} - \boldsymbol{\tau}^v) \cdot \mathbf{u}\bigr] = 0$$

**Capillary mixture energy:**

$$\varepsilon_0 = \sigma\,\|\nabla c\|$$

**Ref:** [wilfong-mfc-26] Section 4.1.8.

---

## 11. Magnetohydrodynamics

### 11.1 Ideal MHD (`mhd = .true.`)

**Continuity:**

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho\,\mathbf{u}) = 0$$

**Momentum:**

$$\frac{\partial (\rho\,\mathbf{u})}{\partial t} + \nabla \cdot \left[\rho\,\mathbf{u} \otimes \mathbf{u} + \left(p + \frac{|\mathbf{B}|^2}{2}\right)\mathbf{I} - \mathbf{B} \otimes \mathbf{B}\right] = 0$$

**Energy:**

$$\frac{\partial \mathcal{E}}{\partial t} + \nabla \cdot \left[\left(\mathcal{E} + p + \frac{|\mathbf{B}|^2}{2}\right)\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}\right] = 0$$

**Induction:**

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{B} - \mathbf{B} \otimes \mathbf{u}) = 0$$

**Total energy:**

$$\mathcal{E} = \rho\,e + \frac{1}{2}\rho\,|\mathbf{u}|^2 + \frac{|\mathbf{B}|^2}{2}$$

**Fast magnetosonic speed:**

$$c_f = \sqrt{\frac{1}{2}\left(c_s^2 + v_A^2 + \sqrt{(c_s^2 + v_A^2)^2 - 4\,c_s^2\,v_A^2\cos^2\theta}\right)}$$

**Alfven speed:**

$$v_A = \sqrt{\frac{|\mathbf{B}|^2}{\rho}}$$

Uses the HLLD Riemann solver (`riemann_solver = 3`). Hyperbolic divergence cleaning (`hyper_cleaning = .true.`) via the GLM method (Dedner et al., 2002).

**Ref:** [wilfong-mfc-26] Section 4.1.9.

### 11.2 Relativistic MHD (`relativity = .true.`)

**Conserved variables:**

$$\mathbf{U} = (D,\;\mathbf{m},\;\tau,\;\mathbf{B})^T$$

where:

$$D = \Gamma\,\rho, \qquad \mathbf{m} = \Gamma^2\rho h\,\mathbf{u} + |\mathbf{B}|^2\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}$$

$$\tau = \Gamma^2\rho h - p + \frac{|\mathbf{B}|^2}{2} + \frac{|\mathbf{u}|^2|\mathbf{B}|^2 - (\mathbf{B} \cdot \mathbf{u})^2}{2} - \Gamma\,\rho$$

Primitive recovery uses Newton-Raphson on the nonlinear conserved-to-primitive relation.

**Ref:** [wilfong-mfc-26] Section 4.1.10.

---

## 12. Information Geometric Regularization (`igr = .true.`)

**Source:** `src/simulation/m_igr.fpp`

**Modified momentum with entropic pressure $\Sigma$:**

$$\frac{\partial (\rho\,\mathbf{u})}{\partial t} + \nabla \cdot \bigl[\rho\,\mathbf{u} \otimes \mathbf{u} + (p + \Sigma)\,\mathbf{I} - \boldsymbol{\tau}\bigr] = 0$$

**Elliptic PDE for $\Sigma$:**

$$\alpha\left[\text{tr}(\nabla \mathbf{u})^2 + \text{tr}^2(\nabla \mathbf{u})\right] = \frac{\Sigma}{\rho} - \alpha\,\nabla \cdot \left(\frac{\nabla \Sigma}{\rho}\right)$$

where $\alpha \sim \Delta x^2$ (regularization strength proportional to mesh spacing squared):

$$\alpha_\text{IGR} = \alpha_\text{factor} \cdot \max(\Delta x,\;\Delta y,\;\Delta z)^2$$

**RHS strain-rate source (3D):**

$$\text{RHS} = \alpha\left[2\left(\frac{\partial u}{\partial y}\frac{\partial v}{\partial x} + \frac{\partial u}{\partial z}\frac{\partial w}{\partial x} + \frac{\partial v}{\partial z}\frac{\partial w}{\partial y}\right) + \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + \left(\frac{\partial w}{\partial z}\right)^2 + (\nabla \cdot \mathbf{u})^2\right]$$

**Iterative solver:** Jacobi (`igr_iter_solver = 1`) or Gauss-Seidel (`igr_iter_solver = 2`), up to `num_igr_iters` iterations (default 5).

Uses Lax-Friedrichs flux (replaces WENO + Riemann solver).

**Ref:** [wilfong-arxiv-25] Eqs. 6-9.

---

## 13. Body Forces (`bodyforces = .true.`)

**Source:** `src/simulation/m_body_forces.fpp`

**Time-dependent acceleration:**

$$a_i(t) = g_i + k_i\sin(\omega_i\,t - \phi_i)$$

**Momentum source:**

$$\frac{\partial (\rho\,u_i)}{\partial t}\bigg|_\text{bf} = \rho\,a_i(t)$$

**Energy source:**

$$\frac{\partial E}{\partial t}\bigg|_\text{bf} = \rho\,\mathbf{u} \cdot \mathbf{a}(t)$$

---

## 14. Acoustic Sources (`acoustic_source = .true.`)

**Source:** `src/simulation/m_acoustic_src.fpp`

Source terms added to the RHS of the governing equations.

**Discrete delta function (spatial support):**

$$\delta_h(r) = \frac{1}{(2\pi\sigma^2)^{d/2}}\exp\!\left(-\frac{r^2}{2\sigma^2}\right)$$

**Forcing form (added to conservative variables):**

$$\mathbf{s}_\text{ac} = \Omega_\Gamma\,f(t)\left[\frac{1}{c_0},\;\cos\theta,\;\sin\theta,\;\frac{c_0^2}{\gamma - 1}\right]^T$$

**Temporal profiles:**
- **Sine** (`pulse = 1`): $S(t) = M\sin(\omega(t - t_\text{delay}))$
- **Gaussian** (`pulse = 2`): $S(t) = M\exp\!\bigl(-\tfrac{1}{2}((t - t_\text{delay})/\sigma_t)^2\bigr)$
- **Square** (`pulse = 3`): $S(t) = M\,\text{sign}(\sin(\omega(t - t_\text{delay})))$
- **Broadband** (`pulse = 4`): superposition of multiple frequencies across a bandwidth

**Spatial supports:** planar, spherical transducer, cylindrical transducer, transducer array (arcuate, annular, circular).

**Ref:** [bryngelson-CPC-20] Section 4.4 Eqs. 35-38.

---

## 15. Numerical Methods

### 15.1 Spatial Reconstruction

#### WENO (`weno_order = 3, 5, 7`)

**Source:** `src/simulation/m_weno.fpp`

Weighted sum of candidate polynomials at cell interfaces:

$$f_{i+1/2} = \sum_r \omega_r\,f_{i+1/2}^{(r)}$$

**WENO-JS** (default):

$$\alpha_r = \frac{d_r}{(\beta_r + \varepsilon)^2}, \qquad \omega_r = \frac{\alpha_r}{\sum_s \alpha_s}$$

where $d_r$ are ideal weights, $\beta_r$ are smoothness indicators, and $\varepsilon$ is a small regularization parameter (`weno_eps`).

**WENO-M** (`mapped_weno = .true.`): Henrick et al. (2005) mapped weights for improved accuracy at critical points:

$$\omega_M^{(r)} = \frac{d_r\bigl(1 + d_r - 3\omega_0^{(r)} + (\omega_0^{(r)})^2\bigr)\,\omega_0^{(r)}}{d_r^2 + \omega_0^{(r)}(1 - 2d_r)}, \qquad \omega^{(r)} = \frac{\omega_M^{(r)}}{\sum_s \omega_M^{(s)}}$$

**WENO-Z** (`wenoz = .true.`): Borges et al. (2008) improved weights with global smoothness measure:

$$\alpha_r = d_r\left(1 + \left(\frac{\tau}{\beta_r + \varepsilon}\right)^q\right), \qquad \tau = |\beta_0 - \beta_{k-1}|$$

The parameter $q$ controls the convergence rate at critical points (typically $q = 1$ for fifth-order reconstruction, as used in MFC).

**TENO** (`teno = .true.`): Fu et al. (2016) targeted ENO with smoothness threshold $C_T$ (`teno_CT`):

$$\gamma_r = 1 + \frac{\tau}{\beta_r}, \qquad \xi_r = \frac{\gamma_r}{\sum_s \gamma_s}$$

If $\xi_r < C_T$, set $\alpha_r = 0$ (stencil excluded).

Primitive variable reconstruction is used to avoid spurious oscillations at interfaces.

**Ref:** [wilfong-mfc-26] Sections 3.2.1, 4.2.2.

#### MUSCL (`muscl_order = 2`)

**Source:** `src/simulation/m_muscl.fpp`

$$q_L(j) = q(j) - \frac{1}{2}\,\phi(r)\,\Delta q, \qquad q_R(j) = q(j) + \frac{1}{2}\,\phi(r)\,\Delta q$$

**Five slope limiters** (`muscl_lim`):
1. **Minmod:** $\phi = \text{sign}(a)\,\min(|a|,\;|b|)$ if $ab > 0$, else $0$
2. **MC (Monotonized Central):** $\phi = \text{sign}(a)\,\min(2|a|,\;2|b|,\;\tfrac{1}{2}(|a|+|b|))$ if $ab > 0$
3. **OSPRE (Van Albada):** $\phi = (a^2 b + a b^2)/(a^2 + b^2)$
4. **Van Leer:** $\phi = 2ab/(a + b)$ if $ab > 0$
5. **Superbee:** $\phi = \max(\min(2|a|,|b|),\;\min(|a|,2|b|))$ if $ab > 0$

where $a$ and $b$ are the left and right slope differences.

**THINC interface compression** (`interface_compression > 0`): applies hyperbolic tangent compression near material interfaces:

$$q_\text{THINC} = q_\min + \frac{q_\max}{2}\left(1 + \text{sign}(s)\,\frac{\tanh(\beta) + A}{1 + A\,\tanh(\beta)}\right)$$

where $A = \exp(\text{sign}(s)\,\beta\,(2C - 1)) / \cosh(\beta) - 1) / \tanh(\beta)$ and $\beta$ controls compression steepness.

#### IGR Reconstruction

5th-order or 3rd-order polynomial interpolation without WENO nonlinear weights, using Lax-Friedrichs numerical flux. Stencil coefficients:

**5th-order:** $\text{coeff}_L = [-3/60,\;27/60,\;47/60,\;-13/60,\;2/60]$

**3rd-order:** $\text{coeff}_L = [2/6,\;5/6,\;-1/6]$

**Ref:** [wilfong-arxiv-25] Algorithm 1.

### 15.2 Riemann Solvers

**Source:** `src/simulation/m_riemann_solvers.fpp`

#### HLL (`riemann_solver = 1`)

$$\mathbf{F}^\text{HLL} = \frac{S_R\,\mathbf{F}_L - S_L\,\mathbf{F}_R + S_L\,S_R\,(\mathbf{q}_R - \mathbf{q}_L)}{S_R - S_L}$$

with wave speed estimates $S_L = \min(u_L - c_L,\;u_R - c_R)$, $S_R = \max(u_R + c_R,\;u_L + c_L)$.

#### HLLC (`riemann_solver = 2`)

Four-state solver resolving the contact discontinuity. Star-state satisfies:

$$p_L^* = p_R^* = p^*, \qquad u_L^* = u_R^* = u^*$$

#### HLLD (`riemann_solver = 3`, MHD only)

Seven-state solver for ideal MHD resolving fast magnetosonic, Alfven, and contact waves (Miyoshi and Kusano, 2005). The Riemann fan is divided by outer wave speeds $S_L$, $S_R$, Alfven speeds $S_L^*$, $S_R^*$, and a middle contact $S_M$:

$$S_M = \frac{(S_R - u_R)\rho_R u_R - (S_L - u_L)\rho_L u_L - p_{T,R} + p_{T,L}}{(S_R - u_R)\rho_R - (S_L - u_L)\rho_L}$$

$$S_L^* = S_M - \frac{|B_x|}{\sqrt{\rho_L^*}}, \qquad S_R^* = S_M + \frac{|B_x|}{\sqrt{\rho_R^*}}$$

where $p_T = p + |\mathbf{B}|^2/2$ is the total (thermal + magnetic) pressure. Continuity of normal velocity and total pressure is enforced across the Riemann fan.

#### Exact (`riemann_solver = 4`)

Iterative exact Riemann solver.

### 15.3 Time Integration

**Source:** `src/simulation/m_time_steppers.fpp`

#### TVD Runge-Kutta (`time_stepper = 1, 2, 3`)

**RK1 (Forward Euler):**

$$\mathbf{q}^{n+1} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)$$

**RK2:**

$$\mathbf{q}^{(1)} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)$$

$$\mathbf{q}^{n+1} = \frac{1}{2}\mathbf{q}^n + \frac{1}{2}\mathbf{q}^{(1)} + \frac{1}{2}\Delta t\,\mathbf{L}(\mathbf{q}^{(1)})$$

**RK3 (SSP):**

$$\mathbf{q}^{(1)} = \mathbf{q}^n + \Delta t\,\mathbf{L}(\mathbf{q}^n)$$

$$\mathbf{q}^{(2)} = \frac{3}{4}\mathbf{q}^n + \frac{1}{4}\mathbf{q}^{(1)} + \frac{1}{4}\Delta t\,\mathbf{L}(\mathbf{q}^{(1)})$$

$$\mathbf{q}^{n+1} = \frac{1}{3}\mathbf{q}^n + \frac{2}{3}\mathbf{q}^{(2)} + \frac{2}{3}\Delta t\,\mathbf{L}(\mathbf{q}^{(2)})$$

#### Adaptive Time Stepping (`adapt_dt = .true.`)

Embedded RK pairs for error estimation with Hairer-Norsett-Wanner algorithm for initial step size.

#### Strang Splitting (for stiff bubble dynamics)

For equations of the form $\partial \mathbf{q}/\partial t = -\nabla \cdot \mathbf{F}(\mathbf{q}) + \mathbf{s}(\mathbf{q})$, the Strang splitting scheme integrates three sub-equations per time step:

$$\frac{\partial \mathbf{q}'}{\partial t} = \mathbf{s}(\mathbf{q}'), \quad t \in [0,\,\Delta t/2], \quad \mathbf{q}'(0) = \mathbf{q}^n$$

$$\frac{\partial \mathbf{q}''}{\partial t} = -\nabla \cdot \mathbf{F}(\mathbf{q}''), \quad t \in [0,\,\Delta t], \quad \mathbf{q}''(0) = \mathbf{q}'(\Delta t/2)$$

$$\frac{\partial \mathbf{q}'''}{\partial t} = \mathbf{s}(\mathbf{q}'''), \quad t \in [0,\,\Delta t/2], \quad \mathbf{q}^{n+1} = \mathbf{q}'''(\Delta t/2)$$

The stiff bubble source terms are integrated using an adaptive embedded Runge-Kutta scheme with error control.

### 15.4 CFL Condition

**Inviscid:**

$$\Delta t = \text{CFL} \cdot \min_{i,j,k}\left(\frac{\Delta x_\text{cell}}{|u| + |v| + |w| + c}\right)$$

**Viscous:**

$$\Delta t_v = \text{CFL} \cdot \min\left(\frac{\Delta x^2}{\nu}\right)$$

where $c$ is the mixture sound speed and $\nu$ is the kinematic viscosity.

### 15.5 Finite Differences (`fd_order = 1, 2, 4`)

**Source:** `src/common/m_finite_differences.fpp`

Used for viscous fluxes and velocity gradients.

**1st-order:**

$$f'_i = \frac{f_{i+1} - f_i}{x_{i+1} - x_i}$$

**2nd-order centered:**

$$f'_i = \frac{f_{i+1} - f_{i-1}}{x_{i+1} - x_{i-1}}$$

**4th-order centered:**

$$f'_i = \frac{-f_{i+2} + 8f_{i+1} - 8f_{i-1} + f_{i-2}}{-x_{i+2} + 8x_{i+1} - 8x_{i-1} + x_{i-2}}$$

**Boundary-adjusted one-sided stencils (3rd-order):**

$$f'_i\big|_\text{left} = \frac{-3f_i + 4f_{i+1} - f_{i+2}}{x_{i+2} - x_i}$$

$$f'_i\big|_\text{right} = \frac{3f_i - 4f_{i-1} + f_{i-2}}{x_i - x_{i-2}}$$

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

### 16.2 Characteristic BCs (Thompson/LODI, `bc_x%beg = -5` to `-12`)

**Characteristic decomposition:**

$$\frac{\partial \mathbf{q}_c}{\partial t} + \mathbf{R}_x\,\boldsymbol{\Lambda}_x\,\mathbf{L}_x\,\frac{\partial \mathbf{q}_c}{\partial x} = 0$$

**Characteristic amplitudes** $\mathcal{L}_1$ through $\mathcal{L}_5$ define wave interactions at boundaries:

$$\mathcal{L}_1 = (u - c)\left(\frac{\partial p}{\partial x} - \rho\,c\,\frac{\partial u}{\partial x}\right), \qquad \mathcal{L}_2 = a\left(\frac{\partial \rho}{\partial x} + \frac{2}{\partial x}\frac{\partial p}{\partial x}\right)$$

$$\mathcal{L}_3 = u\,\frac{\partial v}{\partial x}, \qquad \mathcal{L}_4 = u\,\frac{\partial w}{\partial x}, \qquad \mathcal{L}_5 = (u + c)\left(\frac{\partial p}{\partial x} + \rho\,c\,\frac{\partial u}{\partial x}\right)$$

For non-reflecting boundaries, the incoming wave amplitudes are set to zero.

**GRCBC (incoming wave from ghost point):**

$$\mathcal{L} = -\boldsymbol{\Lambda}^o\,\mathbf{L}_x\,\frac{\partial \mathbf{q}_c}{\partial x} + n_x\,\lambda^i\,\mathbf{L}_x\,\frac{\mathbf{P}\,(\mathbf{q}_p^{(g)} - \mathbf{q}_p)}{\Delta}$$

**8 types:** slip wall (`-5`), non-reflecting buffer (`-6`), non-reflecting sub. inflow (`-7`), non-reflecting sub. outflow (`-8`), force-free sub. outflow (`-9`), constant-pressure sub. outflow (`-10`), supersonic inflow (`-11`), supersonic outflow (`-12`).

**Ref:** [wilfong-mfc-26] Section 4.2.1.

### 16.3 Immersed Boundary Method (`ib = .true.`)

**Source:** `src/simulation/m_ibm.fpp`

Ghost-cell IBM with level set field $\phi$.

**Image point:**

$$\mathbf{x}_{ip} = \mathbf{x}_{gp} + 2\,\phi(\mathbf{x}_{gp})\,\hat{\mathbf{n}}$$

**Velocity boundary conditions:**
- **No-slip:** $\mathbf{u}_{gp} = \mathbf{0}$
- **Slip:** $\mathbf{u}_{gp} = \mathbf{u}_{ip} - (\hat{\mathbf{n}} \cdot \mathbf{u}_{ip})\,\hat{\mathbf{n}}$

**Neumann BC** for pressure, density, and volume fractions.

Supports STL/OBJ geometry import with ray tracing for inside/outside determination.

**Ref:** [wilfong-mfc-26] Section 4.1.1.

---

## 17. Low Mach Number Corrections

**Thornber correction** (`low_Mach = 1`): modifies the reconstructed velocities at cell interfaces:

$$u'_L = \frac{u_L + u_R}{2} + z\,\frac{u_L - u_R}{2}, \qquad u'_R = \frac{u_L + u_R}{2} - z\,\frac{u_L - u_R}{2}$$

$$z = \min\!\left(\max\!\left(\frac{|u_L|}{c_L},\;\frac{|u_R|}{c_R}\right),\;1\right)$$

This reduces numerical dissipation at low Mach numbers while recovering the standard scheme at high Mach numbers.

**Chen correction** (`low_Mach = 2`): anti-dissipation pressure correction (APC) added to the HLLC flux:

$$p_d = \frac{\rho_L\,\rho_R\,(S_R - u_R)(S_L - u_L)}{\rho_R(S_R - u_R) - \rho_L(S_L - u_L)}$$

$$\text{APC}_{p\mathbf{u}} = (z - 1)\,p_d\,\hat{\mathbf{n}}, \qquad \text{APC}_{pE} = (z - 1)\,p_d\,S_*$$

The corrected HLLC flux in the star region becomes $\mathbf{F}^{*} + \text{APC}$.

Additionally, the **artificial Mach number** technique (`pi_fac`) modifies $\pi_\infty$ to artificially increase the local Mach number.

**Ref:** [wilfong-mfc-26] Section 4.2.4.

---

## 18. Flux Limiting

**Volume fraction limiting:** enforces $0 \le \alpha_k \le 1$ and rescales as:

$$\alpha_k \leftarrow \frac{\alpha_k}{\sum_j \alpha_j}$$

**Advective flux limiting** based on local volume fraction gradient $\chi$ to prevent spurious oscillations near material interfaces.

**Ref:** [bryngelson-CPC-20] Section 4.2.
