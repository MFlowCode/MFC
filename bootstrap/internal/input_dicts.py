from internal.common import MFCException


COMPONENTS =["x", "y", "z"]
FLUID_COUNT=10
TAU_E_COUNT=6
PROBE_COUNT=3
MF_COUNT=5
MOMENT_ORDER_COUNT=5
MONO_COUNT=4
INTEGRAL_COUNT=5

PRE_PROCESS: list = ['case_dir', 'old_grid', 'old_ic', 't_step_old', 'm', 'n', 'p', 'cyl_coord', 'model_eqns', 'num_fluids', 'adv_alphan', 'mpp_lim' , 'weno_order', 'precision' , 'parallel_io' , 'perturb_flow', 'perturb_flow_fluid', 'perturb_sph', 'perturb_sph_fluid' , 'fluid_rho', 'hypoelasticity', 'num_patches', 'Ca', 'Web' , 'Re_inv', 'pref', 'rhoref', 'bubbles' , 'polytropic', 'polydisperse', 'poly_sigma', 'thermal', 'nb', 'R0ref' , 'qbmm', 'dist_type' , 'R0_type' , 'nnode' , 'sigR', 'sigV', 'rhoRV']

for cmp in COMPONENTS:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS.append(f"{cmp}_{prepend}")
    
    for append in ["stretch", "a", "loops"]:
        PRE_PROCESS.append(f"{append}_{cmp}")
    
    PRE_PROCESS.append(f"bc_{cmp}%beg")
    PRE_PROCESS.append(f"bc_{cmp}%end")

for f_id in range(1, FLUID_COUNT+1):
    PRE_PROCESS.append(f'fluid_rho({f_id})')

    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G"]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{attribute}")

for p_id in range(1, FLUID_COUNT+1):
    for attribute in ["geometry", "radius", "radii", "epsilon", "beta", "normal", "smoothen", "smooth_patch_id", "alpha_rho", "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma", "pi_inf", "r0", "v0", "p0", "m0"]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{attribute}")

    for cmp_id, cmp in enumerate(COMPONENTS):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_icpp({p_id})%{cmp}_centroid')
        PRE_PROCESS.append(f'patch_icpp({p_id})%length_{cmp}')
        PRE_PROCESS.append(f'patch_icpp({p_id})%radii({cmp_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%normal({cmp_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%vel({cmp_id})')

    for arho_id in range(1, FLUID_COUNT+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha({arho_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha_rho({arho_id})')
    
    for taue_id in range(1, TAU_E_COUNT+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%tau_e({arho_id})')

SIMULATION: list = ['case_dir', 'run_time_info', 't_step_old', 't_tol', 'debug', 'm', 'n', 'p', 'cyl_coord', 'dt', 't_step_start', 't_step_stop', 't_step_save', 'model_eqns', 'num_fluids', 'adv_alphan', 'mpp_lim', 'time_stepper', 'weno_vars', 'weno_order', 'weno_eps', 'char_decomp', 'mapped_weno', 'mp_weno', 'weno_avg', 'weno_Re_flux', 'riemann_solver', 'wave_speeds', 'avg_state', 'commute_err', 'split_err', 'alt_crv', 'alt_soundspeed', 'regularization', 'reg_eps', 'null_weights', 'mixture_err', 'tvd_riemann_flux', 'tvd_rhs_flux', 'tvd_wave_speeds', 'flux_lim', 'We_riemann_flux', 'We_rhs_flux', 'We_src', 'We_wave_speeds', 'lsq_deriv', 'parallel_io', 'precision', 'hypoelasticity', 'fd_order' , 'com_wrt', 'num_probes', 'probe_wrt', 'cb_wrt', 'threshold_mf', 'moment_order', 'pref', 'rhoref', 'polydisperse', 'poly_sigma', 'bubbles', 'bubble_model', 'polytropic', 'thermal', 'R0ref', 'Ca', 'Web', 'Re_inv', 'nb', 'Monopole', 'num_mono', 'qbmm', 'R0_type', 'nnode', 'integral_wrt', 'num_integrals']

for cmp in COMPONENTS:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')

for wrt_id in range(1,FLUID_COUNT+1):
    SIMULATION.append(f'com_wrt({wrt_id})')
    SIMULATION.append(f'cb_wrt({wrt_id})')

    for cmp in COMPONENTS:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,PROBE_COUNT+1):
    for cmp in COMPONENTS:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for mf_id in range(1,MF_COUNT+1):
    SIMULATION.append(f'threshold_mf({mf_id})')

for order_id in range(1,MOMENT_ORDER_COUNT+1):
    SIMULATION.append(f'moment_order({order_id})')

for f_id in range(1,FLUID_COUNT+1):
    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G"]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for mono_id in range(1,MONO_COUNT+1):
        for attribute in ["mag", "length", "dir", "npulse", "pulse", "support", "delay"]:
            SIMULATION.append(f"Mono({mono_id})%{attribute}")
        
        for cmp_id in range(1,len(COMPONENTS)+1):
            SIMULATION.append(f"Mono({mono_id})%loc({cmp_id})")

    for int_id in range(1,INTEGRAL_COUNT+1):
        for cmp in COMPONENTS:
            SIMULATION.append(f"integral({int_id})%{cmp}min")
            SIMULATION.append(f"integral({int_id})%{cmp}max")

    for r_id in range(1,FLUID_COUNT+1):
        # FIXME: Add all fluid_pp({f_id})%...(...)
        # I don't get how some indices are selected

        print("FIX THIS")

POST_PROCESS: list = ['case_dir', 'cyl_coord', 'm', 'n', 'p', 't_step_start', 't_step_stop', 't_step_save', 'model_eqns', 'num_fluids', 'adv_alphan', 'mpp_lim', 'weno_order', 'alt_soundspeed', 'mixture_err', 'parallel_io', 'hypoelasticity', 'polydisperse', 'poly_sigma', 'polytropic', 'thermal', 'pref', 'Ca', 'Web', 'Re_inv', 'rhoref', 'bubbles', 'R0ref', 'nb', 'format', 'precision', 'coarsen_silo', 'fourier_decomp', 'fourier_modes%beg', 'fourier_modes%end', 'alpha_rho_wrt', 'rho_wrt', 'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt', 'pres_wrt', 'alpha_wrt', 'kappa_wrt', 'gamma_wrt', 'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt', 'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt', 'schlieren_wrt', 'schlieren_alpha', 'fd_order']

for cmp_id in range(1,len(COMPONENTS)+1):
    cmp = COMPONENTS[cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg')
    POST_PROCESS.append(f'bc_{cmp}%end')

    for attribute in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{attribute}({cmp_id})')

for fl_id in range(1,FLUID_COUNT+1):
    POST_PROCESS.append(f'schlieren_alpha({fl_id})')
    POST_PROCESS.append(f'alpha_rho_wrt({fl_id})')
    POST_PROCESS.append(f'alpha_wrt({fl_id})')
    POST_PROCESS.append(f'kappa_wrt({fl_id})')

    for attribute in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0"]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{attribute}")

def get_keys(target_name: str) -> list:
    if target_name == "pre_process":
        return [x for x in PRE_PROCESS]
    if target_name == "simulation":
        return [x for x in SIMULATION]
    if target_name == "post_process":
        return [x for x in POST_PROCESS]
    
    raise MFCException(f"[INPUT DICTS] Target {target_name} doesn't have an input dict.")
