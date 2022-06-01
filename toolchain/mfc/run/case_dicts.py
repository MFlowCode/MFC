import common

PRE_PROCESS = ['case_dir', 'old_grid', 'old_ic', 't_step_old', 'm', 'n', 'p',
               'cyl_coord', 'model_eqns', 'num_fluids', 'adv_alphan', 'mpp_lim',
               'weno_order', 'precision', 'parallel_io', 'perturb_flow',
               'perturb_flow_fluid', 'perturb_sph', 'perturb_sph_fluid',
               'fluid_rho', 'hypoelasticity', 'num_patches', 'Ca', 'Web',
               'Re_inv', 'pref', 'rhoref', 'bubbles' , 'polytropic',
               'polydisperse', 'poly_sigma', 'thermal', 'nb', 'R0ref', 'qbmm',
               'dist_type', 'R0_type', 'nnode', 'sigR', 'sigV', 'rhoRV']

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


SIMULATION = ['case_dir', 'run_time_info', 't_step_old', 't_tol', 'debug', 'm',
              'n', 'p', 'cyl_coord', 'dt', 't_step_start', 't_step_stop',
              't_step_save', 'model_eqns', 'num_fluids', 'adv_alphan',
              'mpp_lim', 'time_stepper', 'weno_vars', 'weno_order', 'weno_eps',
              'char_decomp', 'mapped_weno', 'mp_weno', 'weno_avg',
              'weno_Re_flux', 'riemann_solver', 'wave_speeds', 'avg_state',
              'commute_err', 'split_err', 'alt_crv', 'alt_soundspeed',
              'regularization', 'reg_eps', 'null_weights', 'mixture_err',
              'tvd_riemann_flux', 'tvd_rhs_flux', 'tvd_wave_speeds', 'flux_lim',
              'We_riemann_flux', 'We_rhs_flux', 'We_src', 'We_wave_speeds',
              'lsq_deriv', 'parallel_io', 'precision', 'hypoelasticity',
              'fd_order' , 'com_wrt', 'num_probes', 'probe_wrt', 'cb_wrt',
              'threshold_mf', 'moment_order', 'pref', 'rhoref', 'polydisperse',
              'poly_sigma', 'bubbles', 'bubble_model', 'polytropic', 'thermal',
              'R0ref', 'Ca', 'Web', 'Re_inv', 'nb', 'Monopole', 'num_mono',
              'qbmm', 'R0_type', 'nnode', 'integral_wrt', 'num_integrals',
              "cu_mpi"]

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')

for wrt_id in range(1,10+1):
    SIMULATION.append(f'com_wrt({wrt_id})')
    SIMULATION.append(f'cb_wrt({wrt_id})')

    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for mf_id in range(1,5+1):
    SIMULATION.append(f'threshold_mf({mf_id})')

for order_id in range(1,5+1):
    SIMULATION.append(f'moment_order({order_id})')

for f_id in range(1,10+1):
    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G"]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for re_id in [1, 2]:
        SIMULATION.append(f"fluid_pp({f_id})%Re({re_id})")

    for mono_id in range(1,4+1):
        for attribute in ["mag", "length", "dir", "npulse", "pulse", "support",
                          "delay"]:
            SIMULATION.append(f"Mono({mono_id})%{attribute}")

        for cmp_id in range(1,3+1):
            SIMULATION.append(f"Mono({mono_id})%loc({cmp_id})")

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION.append(f"integral({int_id})%{cmp}min")
            SIMULATION.append(f"integral({int_id})%{cmp}max")


POST_PROCESS = ['case_dir', 'cyl_coord', 'm', 'n', 'p', 't_step_start',
                't_step_stop', 't_step_save', 'model_eqns', 'num_fluids',
                'adv_alphan', 'mpp_lim', 'weno_order', 'alt_soundspeed',
                'mixture_err', 'parallel_io', 'hypoelasticity',
                'polydisperse', 'poly_sigma', 'polytropic', 'thermal',
                'pref', 'Ca', 'Web', 'Re_inv', 'rhoref', 'bubbles',
                'R0ref', 'nb', 'format', 'precision', 'coarsen_silo',
                'fourier_decomp', 'fourier_modes%beg',
                'fourier_modes%end', 'alpha_rho_wrt', 'rho_wrt',
                'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt',
                'pres_wrt', 'alpha_wrt', 'kappa_wrt', 'gamma_wrt',
                'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt',
                'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt',
                'schlieren_wrt', 'schlieren_alpha', 'fd_order']

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


def get_input_dict_keys(target_name: str) -> list:
    if target_name == "pre_process":  return PRE_PROCESS.copy()
    if target_name == "simulation":   return SIMULATION.copy()
    if target_name == "post_process": return POST_PROCESS.copy()

    raise common.MFCException(f"[INPUT DICTS] Target {target_name} doesn't have an input dict.")
