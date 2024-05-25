from enum import Enum

from ..state import ARG

class ParamType(Enum):
    INT = 0
    REAL = 1
    LOG = 2
    STR = 3

# Currently doing this manually -- the first automated solution I came up with involved grep and awk, but it seemed quite brittle.
INTEGERS = {
    "m", "n", "p", 'thermal', 't_step_old', 't_step_start', 't_step_stop',
    't_step_save', 'num_patches', 'model_eqns', 'num_fluids', "weno_order",
    "precision", "relax_model", 'perturb_flow_fluid', 'perturb_sph_fluid',
    'dist_type', 'R0_type', 'num_ibs', 't_step_print', 'time_stepper',
    'riemann_solver', 'wave_speeds', 'avg_state', 'fd_order', 'num_probes',
    'bubble_model', 'num_mono', 'R0_type', 'num_integrals', 'format', 
    'flux_lim',
}
REALS = {
    "pref", "Web", "poly_sigma", 'x_domain%beg', 'x_domain%end', 'y_domain%beg',
    'y_domain%end', 'dt', "nb", "rhoref", "R0ref", "Re_inv", "Ca",
    'perturb_flow_mag', 'sigR', 'sigV', 'rhoRV', "palpha_eps", "ptgalpha_eps",
    'pi_fac', 'weno_eps', 'schlieren_alpha', 'fluid_rho',
}
LOGICALS = {
    "hypoelasticity", "cyl_coord", "parallel_io", "run_time_info", "bubbles",
    "polytropic", "polydisperse", "file_per_process", "relax", "mpp_lim", 
    "adv_alphan",  "adv_n", 'old_grid', 'old_ic', 'instability_wave',
    'perturb_flow', 'perturb_sph', 'qbmm', 'ib', 'mapped_weno', 'mp_weno',
    'weno_avg', 'weno_Re_flux', 'prim_vars_wrt', 'alt_soundspeed', 
    'null_weights', 'probe_wrt', 'Monopole', 'integral_wrt', 'cu_mpi', 
    'adap_dt', 'schlieren_wrt', 'alpha_rho_wrt', 'rho_wrt', 'mom_wrt', 
    'vel_wrt', 'flux_wrt', 'E_wrt', 'pres_wrt', 'alpha_wrt', 'kappa_wrt',
    'gamma_wrt', 'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt', 
    'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt', 'qm_wrt',
    'vel_profile', 'mixture_err',
}
STRINGS = {"case_dir"}


def set_type(param_name: str, ty: ParamType):
    match ty:
        case ParamType.INT: 
            INTEGERS.add(param_name)
        case ParamType.REAL:
            REALS.add(param_name)
        case ParamType.LOG:
            LOGICALS.add(param_name)
        case ParamType.STR:
            STRINGS.add(param_name)


COMMON = [
    "hypoelasticity", "cyl_coord", "pref", "p", "parallel_io",
    "Web", "poly_sigma", "case_dir", "thermal", "polytropic",
    "m", "mpp_lim", "R0ref", "adv_alphan", "num_fluids", "model_eqns",
    "nb", "weno_order", "rhoref", "bubbles", "Re_inv", "n", "precision",
    "Ca", "polydisperse", "file_per_process", "relax", "relax_model",
    "adv_n"
]


PRE_PROCESS = COMMON + [
    'old_grid', 'old_ic', 't_step_old', 't_step_start', 'vel_profile',
    'instability_wave', 'perturb_flow', 'perturb_flow_fluid', 'perturb_flow_mag',
    'perturb_sph', 'perturb_sph_fluid', 'fluid_rho', 'num_patches', 'qbmm',
    'dist_type', 'R0_type', 'sigR', 'sigV', 'rhoRV', "palpha_eps", "ptgalpha_eps",
    'pi_fac', 'ib', 'num_ibs',
]

for ib_id in range(1, 10+1):
    for real_attr, ty in [("geometry", ParamType.INT), ("radius", ParamType.REAL),
                          ("theta", ParamType.REAL), ("slip", ParamType.LOG),
                          ("c", ParamType.REAL), ("p", ParamType.REAL), 
                          ("t", ParamType.REAL), ("m", ParamType.REAL)]:
        PRE_PROCESS.append(f"patch_ib({ib_id})%{real_attr}")
        set_type(f"patch_ib({ib_id})%{real_attr}", ty)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_ib({ib_id})%{cmp}_centroid')
        set_type(f'patch_ib({ib_id})%{cmp}_centroid', ParamType.REAL)
        PRE_PROCESS.append(f'patch_ib({ib_id})%length_{cmp}')
        set_type(f'patch_ib({ib_id})%length_{cmp}', ParamType.REAL)

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS.append(f"{cmp}_{prepend}")
        set_type(f"{cmp}_{prepend}", ParamType.REAL)

    for append, ty in [("stretch", ParamType.LOG), ("a", ParamType.REAL), ("loops", ParamType.INT)]:
        PRE_PROCESS.append(f"{append}_{cmp}")
        set_type(f"{append}_{cmp}", ty)

    PRE_PROCESS.append(f"bc_{cmp}%beg")
    PRE_PROCESS.append(f"bc_{cmp}%end")
    set_type(f"bc_{cmp}%beg", ParamType.INT)
    set_type(f"bc_{cmp}%end", ParamType.INT)

for f_id in range(1, 10+1):
    PRE_PROCESS.append(f'fluid_rho({f_id})')
    set_type(f'fluid_rho({f_id})', ParamType.REAL)

    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{real_attr}")
        set_type(f"fluid_pp({f_id})%{real_attr}", ParamType.REAL)

for p_id in range(1, 10+1):
    for real_attr, ty in [("geometry", ParamType.INT), ("smoothen", ParamType.LOG),
                      ("smooth_patch_id", ParamType.INT), ("hcid", ParamType.INT)]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{real_attr}")
        set_type(f"patch_icpp({p_id})%{real_attr}", ty)

    for real_attr in ["radius",  "radii", "epsilon", "beta", "normal", "alpha_rho", 
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0", "cv", "qv", "qvp"]: 
        PRE_PROCESS.append(f"patch_icpp({p_id})%{real_attr}")
        set_type(f"patch_icpp({p_id})%{real_attr}", ParamType.REAL)

    # (cameron): This parameter has since been removed.
    # for i in range(100):
    #     PRE_PROCESS.append(f"patch_icpp({p_id})%Y({i})")

    PRE_PROCESS.append(f"patch_icpp({p_id})%model%filepath")
    set_type(f"patch_icpp({p_id})%model%filepath", ParamType.STR)

    for real_attr in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS.append(f"patch_icpp({p_id})%model%{real_attr}({j})")
            set_type(f"patch_icpp({p_id})%model%{real_attr}({j})", ParamType.REAL)

    PRE_PROCESS.append(f"patch_icpp({p_id})%model%spc")
    set_type(f"patch_icpp({p_id})%model%spc", ParamType.INT)
    PRE_PROCESS.append(f"patch_icpp({p_id})%model%threshold")
    set_type(f"patch_icpp({p_id})%model%threshold", ParamType.REAL)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_icpp({p_id})%{cmp}_centroid')
        set_type(f'patch_icpp({p_id})%{cmp}_centroid', ParamType.REAL)
        PRE_PROCESS.append(f'patch_icpp({p_id})%length_{cmp}')
        set_type(f'patch_icpp({p_id})%length_{cmp}', ParamType.REAL)

        for append in ["radii", "normal", "vel"]:
            PRE_PROCESS.append(f'patch_icpp({p_id})%{append}({cmp_id})')
            set_type(f'patch_icpp({p_id})%{append}({cmp_id})', ParamType.REAL)

    for arho_id in range(1, 10+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha({arho_id})')
        set_type(f'patch_icpp({p_id})%alpha({arho_id})', ParamType.REAL)
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha_rho({arho_id})')
        set_type(f'patch_icpp({p_id})%alpha_rho({arho_id})', ParamType.REAL)

    for taue_id in range(1, 6+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%tau_e({taue_id})')
        set_type(f'patch_icpp({p_id})%tau_e({taue_id})', ParamType.REAL)

    if p_id >= 2:
        PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch')
        set_type(f'patch_icpp({p_id})%alter_patch', ParamType.LOG)

        for alter_id in range(1, p_id):
            PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch({alter_id})')
            set_type(f'patch_icpp({p_id})%alter_patch({alter_id})', ParamType.LOG)

# NOTE: Currently unused.
# for f_id in range(1, 10+1):
#     PRE_PROCESS.append(f"spec_pp({f_id})")


# Removed: 't_tol', 'alt_crv', 'regularization', 'lsq_deriv',
# Feel free to put them back if they are needed once more.
# Be sure to add them to the correct type set at the top of the file too!
SIMULATION = COMMON + [
    'run_time_info', 't_step_old', 'dt', 't_step_start',
    't_step_stop', 't_step_save', 't_step_print', 'time_stepper', 'weno_eps',
    'mapped_weno', 'mp_weno', 'weno_avg', 'weno_Re_flux',
    'riemann_solver', 'wave_speeds', 'avg_state', 'prim_vars_wrt',
    'alt_soundspeed', 'null_weights',
    'mixture_err', 'fd_order', 'num_probes', 'probe_wrt', 
    'bubble_model', 'Monopole', 'num_mono', 'qbmm', 'R0_type', 'integral_wrt', 
    'num_integrals', 'cu_mpi', 'palpha_eps', 'ptgalpha_eps', 
    'pi_fac', 'adap_dt', 'ib', 'num_ibs'
]

# NOTE: Not currently present
# for var in [ 'advection', 'diffusion', 'reactions' ]:
#     SIMULATION.append(f'chem_params%{var}')

for ib_id in range(1, 10+1):
    for real_attr, ty in [("geometry", ParamType.INT), ("radius", ParamType.REAL),
                          ("theta", ParamType.REAL), ("slip", ParamType.LOG),
                          ("c", ParamType.REAL), ("p", ParamType.REAL), 
                          ("t", ParamType.REAL), ("m", ParamType.REAL)]:
        SIMULATION.append(f"patch_ib({ib_id})%{real_attr}")
        set_type(f"patch_ib({ib_id})%{real_attr}", ty)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        SIMULATION.append(f'patch_ib({ib_id})%{cmp}_centroid')
        set_type(f'patch_ib({ib_id})%{cmp}_centroid', ParamType.REAL)
        SIMULATION.append(f'patch_ib({ib_id})%length_{cmp}')
        set_type(f'patch_ib({ib_id})%length_{cmp}', ParamType.REAL)

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg') # Type already set to int
    SIMULATION.append(f'bc_{cmp}%end') # Type already set to int
    SIMULATION.append(f'bc_{cmp}%vb1')
    SIMULATION.append(f'bc_{cmp}%vb2')
    SIMULATION.append(f'bc_{cmp}%vb3')
    SIMULATION.append(f'bc_{cmp}%ve1')
    SIMULATION.append(f'bc_{cmp}%ve2')
    SIMULATION.append(f'bc_{cmp}%ve3')
    set_type(f'bc_{cmp}%vb1', ParamType.REAL)
    set_type(f'bc_{cmp}%vb2', ParamType.REAL)
    set_type(f'bc_{cmp}%vb3', ParamType.REAL)
    set_type(f'bc_{cmp}%ve1', ParamType.REAL)
    set_type(f'bc_{cmp}%ve2', ParamType.REAL)
    set_type(f'bc_{cmp}%ve3', ParamType.REAL)

    for prepend in ["domain%beg", "domain%end"]:
        SIMULATION.append(f"{cmp}_{prepend}")
        set_type(f"{cmp}_{prepend}", ParamType.REAL)

# NOTE: This is now just "probe_wrt"
# for wrt_id in range(1,10+1):
#    for cmp in ["x", "y", "z"]:
#        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')
#        set_type(f'probe_wrt({wrt_id})%{cmp}', ParamType.LOG)

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')
        set_type(f'probe({probe_id})%{cmp}', ParamType.REAL)

for f_id in range(1,10+1):
    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        SIMULATION.append(f"fluid_pp({f_id})%{real_attr}")
        set_type(f"fluid_pp({f_id})%{real_attr}", ParamType.REAL)

    for re_id in [1, 2]:
        SIMULATION.append(f"fluid_pp({f_id})%Re({re_id})")
        set_type(f"fluid_pp({f_id})%Re({re_id})", ParamType.REAL)

    for mono_id in range(1,4+1):
        for int_attr in ["pulse", "support"]:
            SIMULATION.append(f"Mono({mono_id})%{int_attr}")
            set_type(f"Mono({mono_id})%{int_attr}", ParamType.INT)

        for real_attr in ["mag", "length", "dir", "npulse", "delay", 
                          "foc_length", "aperture", "support_width"]:
            SIMULATION.append(f"Mono({mono_id})%{real_attr}")
            set_type(f"Mono({mono_id})%{real_attr}", ParamType.REAL)

        for cmp_id in range(1,3+1):
            SIMULATION.append(f"Mono({mono_id})%loc({cmp_id})")
            set_type(f"Mono({mono_id})%loc({cmp_id})", ParamType.REAL)

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION.append(f"integral({int_id})%{cmp}min")
            SIMULATION.append(f"integral({int_id})%{cmp}max")
            set_type(f"integral({int_id})%{cmp}min", ParamType.REAL)
            set_type(f"integral({int_id})%{cmp}max", ParamType.REAL)


# Removed: 'fourier_modes%beg', 'fourier_modes%end', 'chem_wrt'
# Feel free to return them if they are needed once more.
POST_PROCESS = COMMON + [
    't_step_start', 't_step_stop', 't_step_save', 'alt_soundspeed',
    'mixture_err', 'format', 'schlieren_wrt', 'schlieren_alpha', 'fd_order',
    'alpha_rho_wrt', 'rho_wrt',
    'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt', 'pres_wrt',
    'alpha_wrt', 'kappa_wrt', 'gamma_wrt', 'heat_ratio_wrt', 'pi_inf_wrt',
    'pres_inf_wrt', 'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt','qbmm',
    'qm_wrt'
]

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg') # Type already set to int
    POST_PROCESS.append(f'bc_{cmp}%end') # Type already set to int

    for real_attr in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{real_attr}({cmp_id})')
        set_type(f'{real_attr}({cmp_id})', ParamType.LOG)

# NOTE: `chem_wrt` is missing
# for cmp_id in range(100):
#     POST_PROCESS.append(f'chem_wrt({cmp_id})')

for fl_id in range(1,10+1):
    for append, ty in [("schlieren_alpha", ParamType.REAL), ("alpha_rho_wrt", ParamType.LOG),
                       ("alpha_wrt", ParamType.LOG), ("kappa_wrt", ParamType.LOG)]:
        POST_PROCESS.append(f'{append}({fl_id})')
        set_type(f'{append}({fl_id})', ty)

    for real_attr in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0",
                      "cv", "qv", "qvp" ]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{real_attr}")
        set_type(f"fluid_pp({fl_id})%{real_attr}", ParamType.REAL)

ALL = list(set(PRE_PROCESS + SIMULATION + POST_PROCESS))

CASE_OPTIMIZATION = [ "nb", "weno_order", "num_fluids" ]

TYPED_SET = INTEGERS | REALS | LOGICALS | STRINGS
assert TYPED_SET - set(ALL) == set(), "Found parameter without stage (COMMON, PREPROCESS, SIMULATION, POSTPROCESS)"
assert set(ALL) - TYPED_SET  == set(), "Found untyped parameter!"

_properties = {}
for param in ALL:
    if param in INTEGERS:
        _properties[param] = {"type": "integer"}
    elif param in REALS:
        _properties[param] = {"type": "number"}
    elif param in LOGICALS:
        _properties[param] = {"enum": ["T", "F"]}
    elif param in STRINGS:
        _properties[param] = {"type": "string"}
    else:
        print(f'WARNING: Found parameter {param} without type!')

SCHEMA = {
    "type": "object",
    "properties": _properties
    }


def get_input_dict_keys(target_name: str) -> list:
    result = {
        "pre_process"  : PRE_PROCESS,
        "simulation"   : SIMULATION,
        "post_process" : POST_PROCESS
    }.get(target_name, {}).copy()

    if not ARG("case_optimization") or target_name != "simulation":
        return result

    return [ x for x in result if x not in CASE_OPTIMIZATION ]

