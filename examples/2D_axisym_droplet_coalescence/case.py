#!/usr/bin/env python3
# Phase 1 drainage validation case — acoustic_demulsification project
#   Geometry:  2D axisymmetric, head-on collision of two equal droplets
#   Fluids  :  tetradecane drops in nitrogen ambient
#   Target  :  resolve gas-film drainage; compare h_min(t) to Klaseboer 2000
#   Stop    :  when h_min reaches ~100–200 nm (before numerical rupture)
# Roadmap reference: +Projects/acoustic_demulsification/simulation_roadmap.md
# (Phase 1 — Drainage validation)
# Run:
#   ./mfc.sh run case.py -n <N>          # single node
#   ./mfc.sh run case.py --gpu -n <N>    # GPU build
# Knobs at top of file:
#   SMOKE_TEST     — True: fast end-to-end smoke test (no stretch, fixed dt, 100 steps)
#                    False: production (stretched grid, CFL-adaptive dt, t_stop = 200 us)
#   SCHEME         — 'wenoz' or 'muscl_thinc'
#   DX_FINE        — cell size in drainage film (m)
#   DX_BULK        — cell size in bulk (m)
#   FILM_HALF      — axial half-width of the uniform fine zone (m)
#   WE             — target Weber number (default 30)
# Known issues:
#   - With SMOKE_TEST = False and stretching enabled, simulation hits NaN at
#     corner cell (m, n) after ~10 steps. Cause under investigation. The
#     smoke-test path (uniform grid) runs end-to-end cleanly.

import json
import math

# User-facing knobs

SMOKE_TEST = True  # True: coarse + uniform + 100 fixed-dt steps
SCHEME = "wenoz"  # 'wenoz' | 'muscl_thinc'
FILM_HALF = 2e-6  # ± half-width of the uniform fine zone on the axis (m)
WE = 30.0  # Weber number (coalescence regime)

# Droplet diameter defined here so DX_FINE can reference it; also used below.
D = 300e-6  # droplet diameter, 300 um

# DX_FINE = D/1000 = 300 nm for both modes. Sets the UNIFORM PRE-STRETCH cell
# size domain/N (see Nx/Ny derivation below). MFC's tanh stretch keeps cells
# inside [x_a, x_b] at ~this size and grows cells outward; effective bulk
# cell size is emergent, not a user knob.
DX_FINE = D / 1000.0

# Optional explicit overrides for (Nx, Ny). If set, DX_FINE becomes diagnostic.
NX_OVERRIDE = None
NY_OVERRIDE = None

# Physical constants (tetradecane / nitrogen at ~25 C)

# Tetradecane (liquid)
rho_l = 762.0  # kg/m^3
mu_l = 2.1e-3  # Pa*s
sigma = 0.027  # N/m (≈ 27 mN/m)

# Stiffened-gas fit for tetradecane: gamma = 4.4, p_inf = 3.0e8 Pa
# => c_l = sqrt(gamma*(p+p_inf)/rho) ≈ 1316 m/s at 1 atm — matches bulk data.
gamma_l = 4.4
pi_inf_l = 3.0e8

# Nitrogen (ideal gas)
rho_g = 1.165  # kg/m^3 at 1 atm, 20 C
gamma_g = 1.4
pi_inf_g = 0.0

p_amb = 1.01325e5  # 1 atm

# Droplet geometry and kinematics

R = D / 2.0  # radius, 150 um (D = 300 um defined at top)

# Weber number: We = rho_l * U_rel^2 * D / sigma,  U_rel = 2U per Qian & Law.
# U_rel   = 1.8823 m/s  ->  U = 0.9412 m/s per droplet
U_rel = math.sqrt(WE * sigma / (rho_l * D))
U = 0.5 * U_rel

# Initial axial gap between droplet surfaces (pre-drainage), and centers.
gap0 = 20e-6  # 20 um initial gap
x_center = R + 0.5 * gap0  # droplet centers at x = +/- x_center

# Domain & grid

Lx_half = 5.0 * R  # = 750 um axial half-domain
Ly_max = 3.0 * R  # = 450 um radial extent


def _post_stretch_extent(lo, hi, a, xa, xb, loops):
    """Mimic MFC's hyperbolic-tangent stretching in pre_process/m_grid.f90 to
    compute the post-stretch domain extent. The formula uses normalized coords:
      x <- x / length; xa <- xa / length; xb <- xb / length;
      x <- x/a * (a + log(cosh(a*(x-xa))) + log(cosh(a*(x-xb))) - 2*log(cosh(a*(xb-xa)/2)))
      x <- x * length
    This expands the domain beyond [lo, hi]; patches must cover the result or
    the outer cells are left uninitialized and the sim NaNs at step 1.
    """
    length = abs(hi - lo)
    x0, x1 = lo / length, hi / length
    xa_n, xb_n = xa / length, xb / length
    for _ in range(loops):
        const = 2.0 * math.log(math.cosh(a * (xb_n - xa_n) / 2.0))

        def f(x):
            return x / a * (a + math.log(math.cosh(a * (x - xa_n))) + math.log(math.cosh(a * (x - xb_n))) - const)

        x0, x1 = f(x0), f(x1)
    return x0 * length, x1 * length


# Cell count to hit DX_FINE in the fine zone.
# MFC's tanh stretch (m_grid.f90 s_generate_serial_grid) starts with N uniform
# cells of size dx = domain/N, then remaps their boundaries. The remap
# derivative is ~1 inside [x_a, x_b] (cells stay at ~dx) and grows outside
# (cells expand). So to achieve DX_FINE in the fine zone:
#     N = ceil(domain / DX_FINE)
# Bulk cell sizes after stretch are set by (a_*, loops_*) and are reported in
# the stderr diagnostic at the bottom of this file.
Nx = NX_OVERRIDE if NX_OVERRIDE is not None else int(math.ceil(2.0 * Lx_half / DX_FINE))
Ny = NY_OVERRIDE if NY_OVERRIDE is not None else int(math.ceil(Ly_max / DX_FINE))

# Numerical scheme switches

if SCHEME == "wenoz":
    recon = {
        "recon_type": 1,
        "weno_order": 5,
        "weno_eps": 1e-6,
        "mapped_weno": "F",
        "wenoz": "T",
        "mp_weno": "T",
        "teno": "F",
    }
elif SCHEME == "muscl_thinc":
    recon = {
        "recon_type": 2,
        "muscl_order": 2,
        "muscl_lim": 4,  # Van Leer
        "int_comp": "T",  # THINC interface compression
        "ic_eps": 1e-4,
        "ic_beta": 1.6,
    }
else:
    raise ValueError(f"Unknown SCHEME: {SCHEME!r}")

# Time stepping & grid stretching (depends on SMOKE_TEST)

if SMOKE_TEST:
    # Inching roadmap (edit dt / t_step_stop to advance stages):
    #   stage 0 (old): dt=1e-12, n=100,    save=10    -> 100 ps  total, 13 s/run
    #   stage 1 (now): dt=1e-11, n=10000,  save=1000  -> 100 ns  total, ~21 min @-n1
    #   stage 2:       dt=1e-11, n=100000, save=10000 -> 1 µs    total, ~3.5 hr @-n1
    #   stage 3:       switch to cfl_adap_dt + t_stop=5e-6 (see production block)
    #   stage 4:       refine DX_FINE to 20 nm -> 10 nm
    #   stage 5:       enable stretch_y (currently NaNs at step ~10 — debug first)
    # CFL on 40-nm grid with c_liquid=1316 m/s: dt_max(CFL=0.3)=9.1e-12 s.
    # dt=1e-11 sits at CFL~0.33 — safely inside stability.
    time_block = {
        "dt": 1.0e-11,
        "t_step_start": 0,
        "t_step_stop": 10000,
        "t_step_save": 1000,
    }
    stretch_block = {
        "stretch_x": "T",
        "a_x": 2.0,
        "x_a": -FILM_HALF,
        "x_b": +FILM_HALF,
        "loops_x": 1,
        "stretch_y": "F",
    }
else:
    # Production: CFL-adaptive dt, stretched grid.
    # t_stop = 2e-4 s (~1.25 * tau_inertia = 1.59e-4 s), 40 frames.
    time_block = {
        "cfl_adap_dt": "T",
        "cfl_target": 0.3,
        "t_stop": 2.0e-4,
        "t_save": 5.0e-6,
        "n_start": 0,
    }
    # Production: uniform 300 nm cells everywhere (Nx/Ny already large enough
    # to hit DX_FINE without tanh concentration). Disabling stretch keeps the
    # domain at the intended [-Lx_half, +Lx_half] x [0, Ly_max]; with a=4
    # loops=2 the tanh remap was expanding to ±2892 um x 7 mm and wasting
    # cells on empty bulk. Re-enable if you need non-uniform spacing later.
    stretch_block = {
        "stretch_x": "F",
        "stretch_y": "F",
    }

# Compute post-stretch domain extent so patches cover the actual cell range.
# Add 10% margin on top, to account for cell centers vs boundaries.
if stretch_block["stretch_x"] == "T":
    _xlo, _xhi = _post_stretch_extent(-Lx_half, +Lx_half, stretch_block["a_x"], stretch_block["x_a"], stretch_block["x_b"], stretch_block["loops_x"])
    patch1_length_x = 1.10 * (_xhi - _xlo)
else:
    patch1_length_x = 2.0 * Lx_half

if stretch_block["stretch_y"] == "T":
    _ylo, _yhi = _post_stretch_extent(0.0, Ly_max, stretch_block["a_y"], stretch_block["y_a"], stretch_block["y_b"], stretch_block["loops_y"])
    patch1_length_y = 1.10 * (_yhi - _ylo)
else:
    patch1_length_y = Ly_max

# MFC input deck

eps = 1e-9  # avoid exactly 0/1 volume fractions

case = {
    # grid
    "x_domain%beg": -Lx_half,
    "x_domain%end": +Lx_half,
    "y_domain%beg": 0.0,
    "y_domain%end": +Ly_max,
    "m": Nx - 1,  # MFC uses zero-based cell-count (m = Nx-1)
    "n": Ny - 1,
    "p": 0,  # 2D
    "cyl_coord": "T",  # 2D axisymmetric: x is axial, y is radial
    # boundary conditions
    # Axial (x): outflow/extrapolation on both ends.
    "bc_x%beg": -6,  # TRIAGE: char non-reflecting (was -3)
    "bc_x%end": -6,
    # Radial (y): axisymmetric at y=0, outflow at y_max.
    "bc_y%beg": -2,  # axis
    "bc_y%end": -6,  # TRIAGE: char non-reflecting (was -3)
    # time integration
    "time_stepper": 3,  # TVD-RK3
    # physical model
    "model_eqns": 2,  # 5-equation (Allaire) diffuse interface
    "num_fluids": 2,
    "num_patches": 3,
    "mixture_err": "T",
    "mpp_lim": "T",
    "avg_state": 2,
    "riemann_solver": 2,  # HLLC
    "wave_speeds": 1,
    # surface tension
    "surface_tension": "T",
    "sigma": sigma,
    # fluids (index 1 = tetradecane, 2 = nitrogen)
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),  # 1/(4.4-1) = 0.2941
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),  # 4.4*3e8/3.4 ≈ 3.882e8
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),  # 1/(1.4-1) = 2.5
    "fluid_pp(2)%pi_inf": 0.0,
    # patch 1: nitrogen background filling the whole domain
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": patch1_length_y / 2.0,
    "patch_icpp(1)%length_x": patch1_length_x,
    "patch_icpp(1)%length_y": patch1_length_y,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_amb,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0.0,
    # patch 2: left droplet, moving in +x
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.95,
    "patch_icpp(2)%geometry": 2,  # circle (axisymmetric disk becomes a sphere)
    "patch_icpp(2)%x_centroid": -x_center,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%vel(1)": +U,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_amb,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
    # patch 3: right droplet, moving in -x
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.95,
    "patch_icpp(3)%geometry": 2,
    "patch_icpp(3)%x_centroid": +x_center,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%vel(1)": -U,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": p_amb,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val": 1.0,
    # post processing
    "alpha_wrt(1)": "T",  # liquid volume fraction (the field h_min lives in)
    "alpha_wrt(2)": "T",  # gas volume fraction
    "rho_wrt": "T",
    "pres_wrt": "T",
    "mom_wrt(1)": "T",
    "mom_wrt(2)": "T",
    "cf_wrt": "T",  # color function (requires surface_tension)
    "parallel_io": "T",
    "format": 1,  # 1 = silo-HDF5; 2 = binary
}

case.update(recon)
case.update(time_block)
case.update(stretch_block)

# Diagnostic: report ACHIEVED cell sizes after MFC's tanh stretch (stderr so
# it doesn't pollute the JSON on stdout).
import sys as _sys


def _stretched_cb(m, beg, end, stretch, a, xa, xb, loops):
    cb = [beg + (end - beg) * i / (m + 1) for i in range(-1, m + 1)]  # m+2 bnds
    if stretch:
        length = abs(cb[-1] - cb[0])
        v = [x / length for x in cb]
        xa_n, xb_n = xa / length, xb / length
        for _ in range(loops):
            c = 2.0 * math.log(math.cosh(a * (xb_n - xa_n) / 2.0))
            v = [x / a * (a + math.log(math.cosh(a * (x - xa_n))) + math.log(math.cosh(a * (x - xb_n))) - c) for x in v]
        cb = [x * length for x in v]
    return cb


def _report(axis, m, beg, end):
    sk_stretch = case.get(f"stretch_{axis}", "F") == "T"
    cb = _stretched_cb(m, beg, end, sk_stretch, case.get(f"a_{axis}", 1.0), case.get(f"{axis}_a", 0.0), case.get(f"{axis}_b", 0.0), case.get(f"loops_{axis}", 0))
    cells = [cb[i + 1] - cb[i] for i in range(m + 1)]
    _sys.stderr.write(f"[case.py] {axis}: N={m + 1:6d}  min={min(cells) * 1e9:>9.1f} nm  max={max(cells) * 1e6:>8.2f} um  span=[{cb[0] * 1e6:8.1f},{cb[-1] * 1e6:8.1f}] um\n")


_sys.stderr.write(f"[case.py] SMOKE_TEST={SMOKE_TEST}  DX_FINE target={DX_FINE * 1e9:.0f} nm\n")
_report("x", Nx - 1, -Lx_half, +Lx_half)
_report("y", Ny - 1, 0.0, +Ly_max)
_sys.stderr.write(f"[case.py] total cells = {Nx * Ny:,}\n")
_sys.stderr.flush()

print(json.dumps(case))
