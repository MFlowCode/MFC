#!/usr/bin/env python3
# Thermocapillary migration of a 2D drop in the near-zero-Marangoni limit -- Samareh, Mostaghimi &
# Moreau, Int. J. Heat Mass Transfer 73 (2014) 616-626, Sec. 4.1.1 / Fig. 5. The temperature field is
# imposed through density (rho = rho_coeff/T at uniform p, recovered from the EOS); a large thermal
# conductivity (small Marangoni number Ma) + isothermal gradient walls then hold T near the imposed
# linear profile against the drop's own advection. As Ma -> 0 this approaches Samareh's invariant-T
# limit (v_t / v_YGB ~ 0.80). Ma is very small here (deep limit), so dt is conduction-limited
# (dt ~ Ma) and the run is long -- coarse grids only are practical.
#
# MFC is compressible, so a temperature gradient IS a density gradient: at uniform pressure the EOS
# gives T = (p+p_inf)/((gam-1)*rho*cv), so rho(y) = rho_coeff/T(y) encodes the imposed field. Both
# fluids are identical, so one analytic patch (shared smooth circle eta(x,y) driving color, the
# Laplace pressure jump, and density together) makes the recovered T exactly linear everywhere --
# drop included. T is shifted up by T0 = 10 (Samareh's T0 = 0) to keep rho positive; only gradT and
# sigma_T (which set v_YGB) are physical.

import json

# Marangoni number: the small-Ma realization of Samareh's Ma = 0 limit (smaller -> closer to invariant T)
Ma = 0.001

# Geometry: D = 1 drop, 5D wide x 7.5D tall box; drop 1.5D above the cold floor, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_drop = -Ly / 2.0 + 1.5 * D

Nx = 64  # cells across the box width (Samareh used 64, 128, 256)
dx = W / Nx
Ny = round(Ly / dx)

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1  # dynamic viscosity of both phases (MFC takes Re = 1/mu)

# Imposed linear field T(y) = T0 + gradT*y, encoded as rho(y) = rho_coeff/T(y) at uniform p0
T0 = 10.0
gradT = 2.0 / 15.0  # |dT/dy| = 0.1333 (Samareh)
sigma0 = 0.1
sigma_T = -0.1  # dsigma/dT
G = abs(sigma_T * gradT)  # Marangoni stress scale = 0.013333
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2 (density at the drop's reference temperature T0)
rho_coeff = rho_drop * T0
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(T0) = rho_drop
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's initial position

# Conduction properties from Ma: smaller Ma = larger diffusivity = T held closer to the gradient.
# k = alpha*rho*cp at the reference state; the local diffusivity k/(rho*cp) varies mildly with rho.
alpha_b = G * r**2 / (mu * Ma)  # bulk thermal diffusivity
cp = gam * cv
k_b = alpha_b * rho_drop * cp  # bulk conductivity

eps = 1.0e-9  # trace volume fraction of the (identical) second phase

# One smooth circle profile eta(x,y) -- ~1 in the drop, ~0 outside -- drives the color field, the
# Laplace pressure jump, AND the density together, so the EOS-recovered T is the imposed linear field
# everywhere (drop included): rho tracks p so the (p+p_inf) factors cancel. A single patch => no
# analytic-density leak across patches. Hardcode the drop center/radius as decimals (`r` is a reserved
# IC-expression name).
xc_d, yc_d, r_d = 0.0, y_drop, r
w_if = 0.75 * dx  # interface half-width (~3-cell transition)
laplace = sigma0 / r  # Laplace pressure jump sigma/r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
cf_expr = f"({eta})"
pres_expr = f"{p0:.9f} + {laplace:.9f}*({eta})"
rho_num = (1.0 - eps) * rho_coeff
rho_expr = f"{rho_num:.9f}*({p0 + p_inf:.9f} + {laplace:.9f}*({eta}))/({p0 + p_inf:.9f}*({T0} + {gradT:.9f}*y))"


def T_of_y(y):
    return T0 + gradT * y


# Time stepping: min(acoustic CFL, explicit-conduction limit). Small Ma => conduction-limited dt.
rho_min = rho_coeff / (T0 + gradT * Ly / 2.0)  # hot wall: lowest density, max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
t_r = mu / G  # capillary-thermal time = 7.5
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_b))  # 2D explicit diffusion number
t_step_stop = round(2.0 * t_r / mydt)  # 2 capillary-thermal times (the conduction-sweep window)
t_step_save = max(1, t_step_stop // 80)

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation Algorithm
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Free-slip box (Samareh): reflective side walls, slip walls on y (isothermal walls require -15/-16)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -15,
    "bc_y%end": -15,
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension sigma(T); k_therm > 0 enables conduction
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (continuous phase)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (identical properties)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_b,
    # Single patch over the whole box. The drop lives entirely in the analytic profiles: cf, pres,
    # and alpha_rho all share the one eta(x,y) above, so rho tracks the (jump-carrying) pressure and
    # the EOS-recovered T is exactly the imposed linear field everywhere -- drop included, no bump.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Ly / 2),
    "bc_y%Twall_out": T_of_y(Ly / 2),
}

print(json.dumps(data))
