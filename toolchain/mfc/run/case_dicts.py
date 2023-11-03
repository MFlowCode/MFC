from ..      import common
from ..state import ARG


COMMON = [
    "hypoelasticity", "cyl_coord", "pref", "p", "parallel_io",
    "Web", "poly_sigma", "case_dir", "thermal", "polytropic",
    "m", "mpp_lim", "R0ref", "adv_alphan", "num_fluids", "model_eqns",
    "nb", "weno_order", "rhoref", "bubbles", "Re_inv", "n", "precision",
    "Ca", "polydisperse"
]


PRE_PROCESS = COMMON + [
    'old_grid', 'old_ic', 't_step_old', 't_step_start', 'vel_profile',
    'instability_wave', 'perturb_flow', 'perturb_flow_fluid',
    'perturb_sph', 'perturb_sph_fluid', 'fluid_rho', 'num_patches', 'qbmm',
    'dist_type', 'R0_type', 'sigR', 'sigV', 'rhoRV'
]

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS.append(f"{cmp}_{prepend}")

    for append in ["stretch", "a", "loops"]:
        PRE_PROCESS.append(f"{append}_{cmp}")

    PRE_PROCESS.append(f"bc_{cmp}%beg")
    PRE_PROCESS.append(f"bc_{cmp}%end")

for f_id in range(1, 10+1):
    PRE_PROCESS.append(f'fluid_rho({f_id})')

    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G"]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{attribute}")

for p_id in range(1, 10+1):
    for attribute in ["geometry", "radius", "radii", "epsilon", "beta",
                      "normal", "smoothen", "smooth_patch_id", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0"]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{attribute}")

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_icpp({p_id})%{cmp}_centroid')
        PRE_PROCESS.append(f'patch_icpp({p_id})%length_{cmp}')

        for append in ["radii", "normal", "vel"]:
            PRE_PROCESS.append(f'patch_icpp({p_id})%{append}({cmp_id})')

    for arho_id in range(1, 10+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha({arho_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha_rho({arho_id})')

    for taue_id in range(1, 6+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%tau_e({taue_id})')

    if p_id >= 2:
        PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch')

        for alter_id in range(1, p_id):
            PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch({alter_id})')


SIMULATION = COMMON + [
    'run_time_info', 't_step_old', 't_tol', 'dt', 't_step_start',
    't_step_stop', 't_step_save', 'time_stepper', 'weno_eps',
    'mapped_weno', 'mp_weno', 'weno_avg', 'weno_Re_flux',
    'riemann_solver', 'wave_speeds', 'avg_state', 'prim_vars_wrt',
    'alt_crv', 'alt_soundspeed', 'regularization', 'null_weights',
    'mixture_err', 'lsq_deriv', 'fd_order', 'num_probes', 'probe_wrt', 
    'bubble_model', 'Monopole', 'num_mono', 'qbmm', 'R0_type', 'integral_wrt', 
    'num_integrals', 'cu_mpi'
]

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')

for wrt_id in range(1,10+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for f_id in range(1,10+1):
    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G"]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for re_id in [1, 2]:
        SIMULATION.append(f"fluid_pp({f_id})%Re({re_id})")

    for mono_id in range(1,4+1):
        for attribute in ["mag", "length", "dir", "npulse", "pulse", "support",
                          "delay", "foc_length", "aperture"]:
            SIMULATION.append(f"Mono({mono_id})%{attribute}")

        for cmp_id in range(1,3+1):
            SIMULATION.append(f"Mono({mono_id})%loc({cmp_id})")

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION.append(f"integral({int_id})%{cmp}min")
            SIMULATION.append(f"integral({int_id})%{cmp}max")


POST_PROCESS = COMMON + [
    't_step_start', 't_step_stop', 't_step_save', 'alt_soundspeed',
    'mixture_err', 'format', 'schlieren_wrt', 'schlieren_alpha', 'fd_order',
    'fourier_modes%beg', 'fourier_modes%end', 'alpha_rho_wrt', 'rho_wrt',
    'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt', 'pres_wrt',
    'alpha_wrt', 'kappa_wrt', 'gamma_wrt', 'heat_ratio_wrt', 'pi_inf_wrt',
    'pres_inf_wrt', 'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt',
    'qm_wrt'
]

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg')
    POST_PROCESS.append(f'bc_{cmp}%end')

    for attribute in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{attribute}({cmp_id})')

for fl_id in range(1,10+1):
    for append in ["schlieren_alpha", "alpha_rho_wrt", "alpha_wrt", "kappa_wrt"]:
        POST_PROCESS.append(f'{append}({fl_id})')

    for attribute in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0"]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{attribute}")


CASE_OPTIMIZATION = [ "nb", "weno_order", "num_fluids" ]


def get_input_dict_keys(target_name: str) -> list:
    result = None
    if target_name == "pre_process":  result = PRE_PROCESS.copy()
    if target_name == "simulation":   result = SIMULATION.copy()
    if target_name == "post_process": result = POST_PROCESS.copy()

    if result == None:
        raise common.MFCException(f"[INPUT DICTS] Target {target_name} doesn't have an input dict.")

    if not ARG("case_optimization") or target_name != "simulation":
        return result
    
    return [ x for x in result if x not in CASE_OPTIMIZATION ]
