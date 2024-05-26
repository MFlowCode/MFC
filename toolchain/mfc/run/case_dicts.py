from ..state import ARG

# Currently doing this manually -- the first automated solution I came up with involved grep and awk, but it seemed quite brittle.
INTEGERS = {
    "m", "n", "p", 'thermal', 't_step_old', 't_step_start', 't_step_stop',
    't_step_save', 'num_patches', 'model_eqns', 'num_fluids', "weno_order",
    "precision", "relax_model", 'perturb_flow_fluid', 'perturb_sph_fluid',
    'dist_type', 'R0_type', 'num_ibs', 't_step_print', 'time_stepper',
    'riemann_solver', 'wave_speeds', 'avg_state', 'fd_order', 'num_probes',
    'bubble_model', 'num_mono', 'num_integrals', 'format', 'flux_lim',
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
    'cons_vars_wrt', 'c_wrt', 'omega_wrt', 'qm_wrt', 'vel_profile',
    'mixture_err',
}
STRINGS = {"case_dir"}


def add_param(param_name: str, stage_list: list[str], type_set: set[str]):
    '''Arguments:
    stage_list -- List containing parameters used for a particular stage. 
                  One of COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS.
    type_set   -- Set containing parameters of a certain type. Currently
                  one of INTEGERS, REALS, LOGICALS, STRINGS.'''
    stage_list.append(param_name)
    type_set.add(param_name)


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
    for real_attr, ty in [("geometry", INTEGERS), ("radius", REALS),
                          ("theta", REALS), ("slip", LOGICALS),
                          ("c", REALS), ("p", REALS),
                          ("t", REALS), ("m", REALS)]:
        add_param(f"patch_ib({ib_id})%{real_attr}", PRE_PROCESS, ty)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        add_param(f'patch_ib({ib_id})%{cmp}_centroid', PRE_PROCESS, REALS)
        add_param(f'patch_ib({ib_id})%length_{cmp}', PRE_PROCESS, REALS)

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        add_param(f"{cmp}_{prepend}", PRE_PROCESS, REALS)

    for append, ty in [("stretch", LOGICALS), ("a", REALS), ("loops", INTEGERS)]:
        add_param(f"{append}_{cmp}", PRE_PROCESS, ty)

    add_param(f"bc_{cmp}%beg", PRE_PROCESS, INTEGERS)
    add_param(f"bc_{cmp}%end", PRE_PROCESS, INTEGERS)

for f_id in range(1, 10+1):
    add_param(f'fluid_rho({f_id})', PRE_PROCESS, REALS)

    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        add_param(f"fluid_pp({f_id})%{real_attr}", PRE_PROCESS, REALS)

for p_id in range(1, 10+1):
    for attribute, ty in [("geometry", INTEGERS), ("smoothen", LOGICALS),
                      ("smooth_patch_id", INTEGERS), ("hcid", INTEGERS)]:
        add_param(f"patch_icpp({p_id})%{attribute}", PRE_PROCESS, ty)

    for real_attr in ["radius",  "radii", "epsilon", "beta", "normal", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0", "cv", "qv", "qvp"]: 
        add_param(f"patch_icpp({p_id})%{real_attr}", PRE_PROCESS, REALS)

    # (cameron): This parameter has since been removed.
    # for i in range(100):
    #     PRE_PROCESS.append(f"patch_icpp({p_id})%Y({i})")

    add_param(f"patch_icpp({p_id})%model%filepath", PRE_PROCESS, STRINGS)

    for real_attr in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            add_param(f"patch_icpp({p_id})%model%{real_attr}({j})", PRE_PROCESS, REALS)

    add_param(f"patch_icpp({p_id})%model%spc", PRE_PROCESS, INTEGERS)
    add_param(f"patch_icpp({p_id})%model%threshold", PRE_PROCESS, REALS)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        add_param(f'patch_icpp({p_id})%{cmp}_centroid', PRE_PROCESS, REALS)
        add_param(f'patch_icpp({p_id})%length_{cmp}', PRE_PROCESS, REALS)

        for append in ["radii", "normal", "vel"]:
            add_param(f'patch_icpp({p_id})%{append}({cmp_id})', PRE_PROCESS, REALS)

    for arho_id in range(1, 10+1):
        add_param(f'patch_icpp({p_id})%alpha({arho_id})', PRE_PROCESS, REALS)
        add_param(f'patch_icpp({p_id})%alpha_rho({arho_id})', PRE_PROCESS, REALS)

    for taue_id in range(1, 6+1):
        add_param(f'patch_icpp({p_id})%tau_e({taue_id})', PRE_PROCESS, REALS)

    if p_id >= 2:
        add_param(f'patch_icpp({p_id})%alter_patch', PRE_PROCESS, LOGICALS)

        for alter_id in range(1, p_id):
            add_param(f'patch_icpp({p_id})%alter_patch({alter_id})', PRE_PROCESS, LOGICALS)

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
    for real_attr, ty in [("geometry", INTEGERS), ("radius", REALS),
                          ("theta", REALS), ("slip", LOGICALS),
                          ("c", REALS), ("p", REALS),
                          ("t", REALS), ("m", REALS)]:
        add_param(f"patch_ib({ib_id})%{real_attr}", SIMULATION, ty)

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        add_param(f'patch_ib({ib_id})%{cmp}_centroid', SIMULATION, REALS)
        add_param(f'patch_ib({ib_id})%length_{cmp}', SIMULATION, REALS)

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg') # Type already set to int
    SIMULATION.append(f'bc_{cmp}%end') # Type already set to int
    add_param(f'bc_{cmp}%vb1', SIMULATION, REALS)
    add_param(f'bc_{cmp}%vb2', SIMULATION, REALS)
    add_param(f'bc_{cmp}%vb3', SIMULATION, REALS)
    add_param(f'bc_{cmp}%ve1', SIMULATION, REALS)
    add_param(f'bc_{cmp}%ve2', SIMULATION, REALS)
    add_param(f'bc_{cmp}%ve3', SIMULATION, REALS)

    for prepend in ["domain%beg", "domain%end"]:
        add_param(f"{cmp}_{prepend}", SIMULATION, REALS)

# NOTE: This is now just "probe_wrt"
# for wrt_id in range(1,10+1):
#    for cmp in ["x", "y", "z"]:
#        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')
#        set_type(f'probe_wrt({wrt_id})%{cmp}', LOGICALS)

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        add_param(f'probe({probe_id})%{cmp}', SIMULATION, REALS)

for f_id in range(1,10+1):
    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        add_param(f"fluid_pp({f_id})%{real_attr}", SIMULATION, REALS)

    for re_id in [1, 2]:
        add_param(f"fluid_pp({f_id})%Re({re_id})", SIMULATION, REALS)

    for mono_id in range(1,4+1):
        for int_attr in ["pulse", "support"]:
            add_param(f"Mono({mono_id})%{int_attr}", SIMULATION, INTEGERS)

        for real_attr in ["mag", "length", "dir", "npulse", "delay",
                          "foc_length", "aperture", "support_width"]:
            add_param(f"Mono({mono_id})%{real_attr}", SIMULATION, REALS)

        for cmp_id in range(1,3+1):
            add_param(f"Mono({mono_id})%loc({cmp_id})", SIMULATION, REALS)

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            add_param(f"integral({int_id})%{cmp}min", SIMULATION, REALS)
            add_param(f"integral({int_id})%{cmp}max", SIMULATION, REALS)


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
        add_param(f'{real_attr}({cmp_id})', POST_PROCESS, LOGICALS)

# NOTE: `chem_wrt` is missing
# for cmp_id in range(100):
#     POST_PROCESS.append(f'chem_wrt({cmp_id})')

for fl_id in range(1,10+1):
    for append, ty in [("schlieren_alpha", REALS), ("alpha_rho_wrt", LOGICALS),
                       ("alpha_wrt", LOGICALS), ("kappa_wrt", LOGICALS)]:
        add_param(f'{append}({fl_id})', POST_PROCESS, ty)

    for real_attr in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0",
                      "cv", "qv", "qvp" ]:
        add_param(f"fluid_pp({fl_id})%{real_attr}", POST_PROCESS, REALS)

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
