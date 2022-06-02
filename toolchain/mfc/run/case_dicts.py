import typing
import dataclasses

DFLT_REAL = -1e6
DFLT_INT  = -100

class CaseAttribute:
    pass

@dataclasses.dataclass
class CaseAttribute:
    name:    str
    default: typing.Any # If None, it must be set by the user
    inner:   typing.List[CaseAttribute] = None

# Shorthand for CaseAttribute
MFCCA = CaseAttribute

PARAMS: typing.List[CaseAttribute] = [
    MFCCA('case_dir',           '.'),
    MFCCA('old_grid',           None),
    MFCCA('old_ic',             None),
    MFCCA('t_step_old',         None),
    MFCCA('m',                  None),
    MFCCA('n',                  0),
    MFCCA('p',                  0),
    MFCCA('cyl_coord',          None),
    MFCCA('model_eqns',         None),
    MFCCA('adv_alphan',         None),
    MFCCA('mpp_lim',            None),
    MFCCA('weno_order',         None),
    MFCCA('precision',          None),
    MFCCA('parallel_io',        None),
    MFCCA('perturb_flow',       None),
    MFCCA('perturb_flow_fluid', None),
    MFCCA('perturb_sph',        None),
    MFCCA('perturb_sph_fluid',  None),
    MFCCA('fluid_rho',          None),
    MFCCA('hypoelasticity',     None),
    MFCCA('num_patches',        None),
    MFCCA('Ca',                 None),
    MFCCA('Web',                None),
    MFCCA('Re_inv',             None),
    MFCCA('pref',         None),
    MFCCA('rhoref',       None),
    MFCCA('bubbles',      None),
    MFCCA('polytropic',   None),
    MFCCA('polydisperse', None),
    MFCCA('poly_sigma', None),
    MFCCA('thermal', None),
    MFCCA('nb', None),
    MFCCA('R0ref', None),
    MFCCA('qbmm', None),
    MFCCA('dist_type', None),
    MFCCA('R0_type', None),
    MFCCA('nnode', None),
    MFCCA('sigR', None),
    MFCCA('sigV', None),
    MFCCA('rhoRV', None),
    MFCCA('run_time_info', None),
    MFCCA('t_tol', None),
    MFCCA('debug', None),
    MFCCA('dt', None),
    MFCCA('t_step_start', None),
    MFCCA('t_step_stop', None),
    MFCCA('t_step_save', None),
    MFCCA('num_fluids', None),
    MFCCA('time_stepper', None),
    MFCCA('weno_vars', None),
    MFCCA('weno_eps', None),
    MFCCA('char_decomp', None),
    MFCCA('mapped_weno', None),
    MFCCA('mp_weno', None),
    MFCCA('weno_avg', None),
    MFCCA('weno_Re_flux', None),
    MFCCA('riemann_solver', None),
    MFCCA('wave_speeds', None),
    MFCCA('avg_state', None),
    MFCCA('commute_err', None),
    MFCCA('split_err', None),
    MFCCA('alt_crv', None),
    MFCCA('alt_soundspeed', None),
    MFCCA('regularization', None),
    MFCCA('reg_eps', None),
    MFCCA('null_weights', None),
    MFCCA('mixture_err', None),
    MFCCA('tvd_riemann_flux', None),
    MFCCA('tvd_rhs_flux', None),
    MFCCA('tvd_wave_speeds', None),
    MFCCA('flux_lim', None),
    MFCCA('We_riemann_flux', None),
    MFCCA('We_rhs_flux', None),
    MFCCA('We_src', None),
    MFCCA('We_wave_speeds', None),
    MFCCA('lsq_deriv', None),
    MFCCA('fd_order', None),
    MFCCA('com_wrt', None),
    MFCCA('num_probes', None),
    MFCCA('probe_wrt', None),
    MFCCA('cb_wrt', None),
    MFCCA('threshold_mf', None),
    MFCCA('moment_order', None),
    MFCCA('bubble_model', None),
    MFCCA('Monopole', None),
    MFCCA('num_mono', None),
    MFCCA('integral_wrt', None),
    MFCCA('num_integrals', None),
    MFCCA('cu_mpi', None),
    MFCCA('format', None),
    MFCCA('coarsen_silo', None),
    MFCCA('fourier_decomp', None),
    MFCCA('fourier_modes%beg', None),
    MFCCA('fourier_modes%end', None),
    MFCCA('alpha_rho_wrt', None),
    MFCCA('rho_wrt', None),
    MFCCA('mom_wrt', None),
    MFCCA('vel_wrt', None),
    MFCCA('flux_wrt', None),
    MFCCA('E_wrt', None),
    MFCCA('pres_wrt', None),
    MFCCA('alpha_wrt', None),
    MFCCA('kappa_wrt', None),
    MFCCA('gamma_wrt', None),
    MFCCA('heat_ratio_wrt', None),
    MFCCA('pi_inf_wrt', None),
    MFCCA('pres_inf_wrt', None),
    MFCCA('cons_vars_wrt', None),
    MFCCA('prim_vars_wrt', None),
    MFCCA('c_wrt', None),
    MFCCA('omega_wrt', None),
    MFCCA('schlieren_wrt', None),
    MFCCA('schlieren_alpha', None),
    MFCCA('fluids',  [], [
        MFCCA("gamma", None),
        MFCCA("pi_inf", None),
        MFCCA("mul0", None),
        MFCCA("ss", None),
        MFCCA("pv", None),
        MFCCA("gamma_v", None),
        MFCCA("M_v", None),
        MFCCA("mu_v", None),
        MFCCA("k_v", None),
        MFCCA("G", None),
        MFCCA("Re(1)", None),
        MFCCA("Re(2)", None),
    ]),
    MFCCA('patches', [], [
        MFCCA("geometry", None),
        MFCCA("radius", None),
        MFCCA("radii", None),
        MFCCA("epsilon", None),
        MFCCA("beta", None),
        MFCCA("normal", None),
        MFCCA("smoothen", None),
        MFCCA("smooth_patch_id", None),
        MFCCA("alpha_rho", None),
        MFCCA("smooth_coeff", None),
        MFCCA("rho", None),
        MFCCA("vel", None),
        MFCCA("pres", None),
        MFCCA("alpha", None),
        MFCCA("gamma", None),
        MFCCA("pi_inf", None),
        MFCCA("r0", None),
        MFCCA("v0", None),
        MFCCA("p0", None),
        MFCCA("m0", None),
    ]),
]

for wrt_id in range(1,10+1):
    PARAMS.append(MFCCA(f'com_wrt({wrt_id})', None))
    PARAMS.append(MFCCA(f'cb_wrt({wrt_id})', None))

    for cmp in ["x", "y", "z"]:
        PARAMS.append(MFCCA(f'probe_wrt({wrt_id})%{cmp}', None))

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        PARAMS.append(MFCCA(f'probe({probe_id})%{cmp}', None))

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PARAMS.append(MFCCA(f"{cmp}_{prepend}", DFLT_REAL))

    for append in [("stretch", False), ("a", DFLT_REAL), ("loops", 1)]:
        PARAMS.append(MFCCA(f"{append[0]}_{cmp}", append[1]))

    PARAMS.append(MFCCA(f"bc_{cmp}%beg", DFLT_REAL))
    PARAMS.append(MFCCA(f"bc_{cmp}%end", DFLT_REAL))

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PARAMS["patches"].inner.append(MFCCA(f'{cmp}_centroid'))
        PARAMS["patches"].inner.append(MFCCA(f'length_{cmp}'))

        for append in ["radii", "normal", "vel"]:
            PARAMS["patches"].inner.append(MFCCA(f'{append}({cmp_id})'))

def generate(user: dict, master: dict = None):
    if master is None:
        master = PARAMS.copy()

    result: dict = {}

    if not set(user.keys()).issubset(set([ e.name for e in master ])):
        print(set(user.keys()) - set([ e.name for e in master ]))
        print(set(user.keys()))
        print(set([ e.name for e in master ]))
        raise "case_dicts::generate invalid."

    for param in PARAMS:
        entry: dict = {
            param.name: user[param.name] if param.name in user else param.default
        }

        if param.inner is not None:
            for entry in entry[param.name]:
                entry = generate(entry, param.inner)

        result.update(entry.copy())

    return result
