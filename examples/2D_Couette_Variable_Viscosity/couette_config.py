#!/usr/bin/env python3
"""Single source of truth for the variable-viscosity Couette validation.

case.py (the MFC input), reference.py (the exact solution), and validate.py
(the comparison) all import these constants so they can never drift apart.

Physics: planar Couette flow of a single fluid between two no-slip walls a
distance H apart. The bottom wall (y = 0) is held fixed at temperature T0; the
top wall (y = H) slides at speed U and is held at temperature T1. The fluid has
a temperature-dependent (Arrhenius) dynamic viscosity

    mu(T) = exp(C + D / T)            with D > 0  ->  mu falls as T rises

so the hot side near the top wall is thinner than the cold side near the
bottom. In steady state the shear stress mu(T) du/dy is uniform across the gap,
which forces du/dy to vary inversely with mu -- the velocity profile is *curved*
rather than the straight line a constant-viscosity code would give. That
curvature is the signature of mu(T), and it is what this case validates against
an exact solution.

MFC is a compressible solver, so we deliberately soften the sound speed (small
cv) to keep the flow at low Mach number; that makes the incompressible Couette
solution in reference.py a faithful reference and keeps the acoustic time step
affordable. Temperature is not a stored field: it follows from the stiffened-gas
EOS T = (p + p_inf)/((gam - 1) rho cv), so a target temperature profile is
imposed through a density profile rho(T) at uniform pressure.
"""

import numpy as np

# Gas / EOS
gam = 1.4  # ratio of specific heats
p_inf = 0.0  # ideal gas (no stiffening needed: the soft cv already softens c)

# Geometry and drive (the streamwise x-extent is set per-grid in case.py to keep
# cells square; x is periodic with no gradient, so its size is physically free).
H = 1.0  # wall-normal gap (y extent)
U = 1.0  # top-wall speed (x-component)

# Wall temperatures (bottom fixed at T0, top moving at T1)
T0 = 300.0
T1 = 400.0
T_ref = 0.5 * (T0 + T1)  # mid-gap reference temperature

# Soft sound speed: c = sqrt(gam*(gam-1)*cv*T). Pick c_ref at T_ref, back out cv.
c_ref = 10.0
cv = c_ref**2 / (gam * (gam - 1.0) * T_ref)
cp = gam * cv

# Uniform pressure giving rho = rho_ref at T_ref (mechanical equilibrium: no body
# force, so a y-varying density at uniform pressure is a valid steady state).
rho_ref = 1.0
p0 = (gam - 1.0) * rho_ref * cv * T_ref

# Arrhenius viscosity mu(T) = exp(C + D/T).
mu_contrast = 3.5  # target mu(T0)/mu(T1) across the gap
D = np.log(mu_contrast) / (1.0 / T0 - 1.0 / T1)
mu_ref = 0.05  # dynamic viscosity at T_ref  ->  Re = rho_ref*U*H/mu_ref
C = np.log(mu_ref) - D / T_ref

# Thermal conductivity from a chosen Prandtl number Pr = mu*cp/k.
Pr = 1.0
k_therm = mu_ref * cp / Pr


def mu_of_T(T):
    """Arrhenius dynamic viscosity, same formula MFC evaluates (visc_model = 1)."""
    return np.exp(C + D / T)


def T_linear(y):
    """Pure-conduction (Brinkman -> 0) temperature profile across the gap."""
    return T0 + (T1 - T0) * (y / H)


def rho_of_T(T):
    """Density that encodes temperature T at the uniform pressure p0 via the EOS."""
    return (p0 + p_inf) / ((gam - 1.0) * cv * T)


# Derived dimensionless diagnostics (reported, not tuned).
Re = rho_ref * U * H / mu_ref  # Reynolds number at T_ref
Ma = U / c_ref  # Mach number at T_ref
Br = mu_ref * U**2 / (k_therm * (T1 - T0))  # Brinkman: viscous-heating strength

# Slowest viscous relaxation mode sets how long to run to reach steady state.
nu_min = mu_of_T(T1) / rho_of_T(T1)  # smallest kinematic viscosity (hot wall)
nu_max = mu_of_T(T0) / rho_of_T(T0)  # largest kinematic viscosity (cold wall)
tau_relax = H**2 / (np.pi**2 * nu_min)  # first-mode decay time
t_end = 10.0 * tau_relax  # run ~10 decay times: steady transient well below grid error


if __name__ == "__main__":
    print(f"cv        = {cv:.6g}   cp = {cp:.6g}")
    print(f"p0        = {p0:.6g}")
    print(f"C         = {C:.6g}   D = {D:.6g}")
    print(f"k_therm   = {k_therm:.6g}")
    print(f"mu(T0)    = {mu_of_T(T0):.6g}   mu(T1) = {mu_of_T(T1):.6g}   contrast = {mu_of_T(T0) / mu_of_T(T1):.3f}")
    print(f"Re        = {Re:.4g}   Ma = {Ma:.4g}   Br = {Br:.4g}   Pr = {Pr:.4g}")
    print(f"nu_min    = {nu_min:.4g}   nu_max = {nu_max:.4g}")
    print(f"tau_relax = {tau_relax:.4g}   t_end = {t_end:.4g}")
