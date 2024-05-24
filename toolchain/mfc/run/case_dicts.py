from ..state import ARG


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
    for attribute in ["geometry", "radius", "theta", "slip", "c", "p", "t", "m"]:
        PRE_PROCESS.append(f"patch_ib({ib_id})%{attribute}")

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_ib({ib_id})%{cmp}_centroid')
        PRE_PROCESS.append(f'patch_ib({ib_id})%length_{cmp}')

    for attribute in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS.append(f"patch_ib({ib_id})%model%{attribute}({j})")
    
    PRE_PROCESS.append(f"patch_ib({ib_id})%model%filepath")
    PRE_PROCESS.append(f"patch_ib({ib_id})%model%spc")
    PRE_PROCESS.append(f"patch_ib({ib_id})%model%threshold")


    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_scale({1})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_scale({2})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_scale({3})")

    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_translate({1})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_translate({2})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_translate({3})")

    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_rotate({1})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_rotate({2})")
    # PRE_PROCESS.append(f"patch_ib({ib_id})%model_rotate({3})")

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
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{attribute}")

for p_id in range(1, 10+1):
    for attribute in ["geometry", "radius", "radii", "epsilon", "beta",
                      "normal", "smoothen", "smooth_patch_id", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0", "hcid", "cv", "qv", "qvp" ]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{attribute}")

    for i in range(100):
        PRE_PROCESS.append(f"patch_icpp({p_id})%Y({i})")

    PRE_PROCESS.append(f"patch_icpp({p_id})%model%filepath")

    for attribute in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS.append(f"patch_icpp({p_id})%model%{attribute}({j})")

    PRE_PROCESS.append(f"patch_icpp({p_id})%model%spc")
    PRE_PROCESS.append(f"patch_icpp({p_id})%model%threshold")

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


for f_id in range(1, 10+1):
    PRE_PROCESS.append(f"spec_pp({f_id})")


SIMULATION = COMMON + [
    'run_time_info', 't_step_old', 't_tol', 'dt', 't_step_start',
    't_step_stop', 't_step_save', 't_step_print', 'time_stepper', 'weno_eps',
    'mapped_weno', 'mp_weno', 'weno_avg', 'weno_Re_flux',
    'riemann_solver', 'wave_speeds', 'avg_state', 'prim_vars_wrt',
    'alt_crv', 'alt_soundspeed', 'regularization', 'null_weights',
    'mixture_err', 'lsq_deriv', 'fd_order', 'num_probes', 'probe_wrt', 
    'bubble_model', 'Monopole', 'num_mono', 'qbmm', 'R0_type', 'integral_wrt', 
    'num_integrals', 'cu_mpi', 'palpha_eps', 'ptgalpha_eps', 
    'pi_fac', 'adap_dt', 'ib', 'num_ibs'
]

for var in [ 'advection', 'diffusion', 'reactions' ]:
    SIMULATION.append(f'chem_params%{var}')

for ib_id in range(1, 10+1):
    for attribute in ["geometry", "radius", "theta","slip", "c", "m", "t", "p"]:
        SIMULATION.append(f"patch_ib({ib_id})%{attribute}")

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        SIMULATION.append(f'patch_ib({ib_id})%{cmp}_centroid')
        SIMULATION.append(f'patch_ib({ib_id})%length_{cmp}')

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')
    SIMULATION.append(f'bc_{cmp}%vb1')
    SIMULATION.append(f'bc_{cmp}%vb2')
    SIMULATION.append(f'bc_{cmp}%vb3')
    SIMULATION.append(f'bc_{cmp}%ve1')
    SIMULATION.append(f'bc_{cmp}%ve2')
    SIMULATION.append(f'bc_{cmp}%ve3')
    for prepend in ["domain%beg", "domain%end"]:
        SIMULATION.append(f"{cmp}_{prepend}")

for wrt_id in range(1,10+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for f_id in range(1,10+1):
    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for re_id in [1, 2]:
        SIMULATION.append(f"fluid_pp({f_id})%Re({re_id})")

    for mono_id in range(1,4+1):
        for attribute in ["mag", "length", "dir", "npulse", "pulse", "support",
                          "delay", "foc_length", "aperture", "support_width"]:
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
    'pres_inf_wrt', 'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt','qbmm',
    'qm_wrt', 'chem_wrt'
]

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg')
    POST_PROCESS.append(f'bc_{cmp}%end')

    for attribute in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{attribute}({cmp_id})')

for cmp_id in range(100):
    POST_PROCESS.append(f'chem_wrt({cmp_id})')

for fl_id in range(1,10+1):
    for append in ["schlieren_alpha", "alpha_rho_wrt", "alpha_wrt", "kappa_wrt"]:
        POST_PROCESS.append(f'{append}({fl_id})')

    for attribute in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0",
                      "cv", "qv", "qvp" ]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{attribute}")

ALL = list(set(PRE_PROCESS + SIMULATION + POST_PROCESS))

CASE_OPTIMIZATION = [ "nb", "weno_order", "num_fluids"]


def get_input_dict_keys(target_name: str) -> list:
    result = {
        "pre_process"  : PRE_PROCESS,
        "simulation"   : SIMULATION,
        "post_process" : POST_PROCESS
    }.get(target_name, {}).copy()

    if not ARG("case_optimization") or target_name != "simulation":
        return result

    return [ x for x in result if x not in CASE_OPTIMIZATION ]
