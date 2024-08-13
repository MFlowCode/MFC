from enum import Enum
from ..state import ARG

class ParamType(Enum):
    INT = {"type": "integer"}
    REAL = {"type": "number"}
    LOG = {"enum": ["T", "F"]}
    STR = {"type": "string"}

    _ANALYTIC_INT = {"type": ["integer", "string"]}
    _ANALYTIC_REAL = {"type": ["number", "string"]}

    def analytic(self):
        if self == self.INT:
            return self._ANALYTIC_INT
        if self == self.REAL:
            return self._ANALYTIC_REAL
        return self.STR

COMMON = {
    'hypoelasticity': ParamType.LOG,
    'cyl_coord': ParamType.LOG,
    'pref': ParamType.REAL,
    'p': ParamType.INT,
    'parallel_io': ParamType.LOG,
    'Web': ParamType.REAL,
    'poly_sigma': ParamType.REAL,
    'case_dir': ParamType.STR,
    'thermal': ParamType.INT,
    'polytropic': ParamType.LOG,
    'm': ParamType.INT,
    'mpp_lim': ParamType.LOG,
    'R0ref': ParamType.REAL,
    'num_fluids': ParamType.INT,
    'model_eqns': ParamType.INT,
    'nb': ParamType.REAL,
    'weno_order': ParamType.INT,
    'rhoref': ParamType.REAL,
    'bubbles': ParamType.LOG,
    'Re_inv': ParamType.REAL,
    'n': ParamType.INT,
    'precision': ParamType.INT,
    'Ca': ParamType.REAL,
    'polydisperse': ParamType.LOG,
    'file_per_process': ParamType.LOG,
    'relax': ParamType.LOG,
    'relax_model': ParamType.INT,
    'sigma': ParamType.REAL,
    'adv_n': ParamType.LOG,
}

PRE_PROCESS = COMMON.copy()
PRE_PROCESS.update({
    'old_grid': ParamType.LOG,
    'old_ic': ParamType.LOG,
    't_step_old': ParamType.INT,
    't_step_start': ParamType.INT,
    'mixlayer_vel_profile': ParamType.LOG,
    'mixlayer_vel_coef': ParamType.REAL,
    'mixlayer_domain': ParamType.REAL,
    'mixlayer_perturb': ParamType.LOG,
    'perturb_flow': ParamType.LOG,
    'perturb_flow_fluid': ParamType.INT,
    'perturb_flow_mag': ParamType.REAL,
    'perturb_sph': ParamType.LOG,
    'perturb_sph_fluid': ParamType.INT,
    'fluid_rho': ParamType.REAL,
    'num_patches': ParamType.INT,
    'qbmm': ParamType.LOG,
    'dist_type': ParamType.INT,
    'R0_type': ParamType.INT,
    'sigR': ParamType.REAL,
    'sigV': ParamType.REAL,
    'rhoRV': ParamType.REAL,
    'palpha_eps': ParamType.REAL,
    'ptgalpha_eps': ParamType.REAL,
    'pi_fac': ParamType.REAL,
    'ib': ParamType.LOG,
    'num_ibs': ParamType.INT,
})

for ib_id in range(1, 10+1):
    for real_attr, ty in [("geometry", ParamType.INT), ("radius", ParamType.REAL),
                          ("theta", ParamType.REAL), ("slip", ParamType.LOG),
                          ("c", ParamType.REAL), ("p", ParamType.REAL),
                          ("t", ParamType.REAL), ("m", ParamType.REAL)]:
        PRE_PROCESS[f"patch_ib({ib_id})%{real_attr}"] = ty

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS[f'patch_ib({ib_id})%{cmp}_centroid'] = ParamType.REAL
        PRE_PROCESS[f'patch_ib({ib_id})%length_{cmp}'] = ParamType.REAL

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS[f"{cmp}_{prepend}"] = ParamType.REAL

    for append, ty in [("stretch", ParamType.LOG), ("a", ParamType.REAL),
                       ("loops", ParamType.INT)]:
        PRE_PROCESS[f"{append}_{cmp}"] = ty

    PRE_PROCESS[f"bc_{cmp}%beg"] = ParamType.INT
    PRE_PROCESS[f"bc_{cmp}%end"] = ParamType.INT

for f_id in range(1, 10+1):
    PRE_PROCESS[f'fluid_rho({f_id})'] = ParamType.REAL

    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        PRE_PROCESS[f"fluid_pp({f_id})%{real_attr}"] = ParamType.REAL

for p_id in range(1, 10+1):
    for attribute, ty in [("geometry", ParamType.INT), ("smoothen", ParamType.LOG),
                      ("smooth_patch_id", ParamType.INT), ("hcid", ParamType.INT)]:
        PRE_PROCESS[f"patch_icpp({p_id})%{attribute}"] = ty

    for real_attr in ["radius",  "radii", "epsilon", "beta", "normal", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0", "cv", "qv", "qvp", "cf_val"]: 
        PRE_PROCESS[f"patch_icpp({p_id})%{real_attr}"] = ParamType.REAL
    PRE_PROCESS[f"patch_icpp({p_id})%pres"] = ParamType.REAL.analytic()

    # (cameron): This parameter has since been removed.
    # for i in range(100):
    #     PRE_PROCESS.append(f"patch_icpp({p_id})%Y({i})")

    PRE_PROCESS[f"patch_icpp({p_id})%model%filepath"] = ParamType.STR

    for real_attr in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS[f"patch_icpp({p_id})%model%{real_attr}({j})"] = ParamType.REAL

    PRE_PROCESS[f"patch_icpp({p_id})%model%spc"] = ParamType.INT
    PRE_PROCESS[f"patch_icpp({p_id})%model%threshold"] = ParamType.REAL

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS[f'patch_icpp({p_id})%{cmp}_centroid'] = ParamType.REAL
        PRE_PROCESS[f'patch_icpp({p_id})%length_{cmp}'] = ParamType.REAL

        for append in ["radii", "normal"]:
            PRE_PROCESS[f'patch_icpp({p_id})%{append}({cmp_id})'] = ParamType.REAL
        PRE_PROCESS[f'patch_icpp({p_id})%vel({cmp_id})'] = ParamType.REAL.analytic()

    for arho_id in range(1, 10+1):
        PRE_PROCESS[f'patch_icpp({p_id})%alpha({arho_id})'] = ParamType.REAL.analytic()
        PRE_PROCESS[f'patch_icpp({p_id})%alpha_rho({arho_id})'] = ParamType.REAL.analytic()

    for taue_id in range(1, 6+1):
        PRE_PROCESS[f'patch_icpp({p_id})%tau_e({taue_id})'] = ParamType.REAL.analytic()

    if p_id >= 2:
        PRE_PROCESS[f'patch_icpp({p_id})%alter_patch'] = ParamType.LOG

        for alter_id in range(1, p_id):
            PRE_PROCESS[f'patch_icpp({p_id})%alter_patch({alter_id})'] = ParamType.LOG

# NOTE: Currently unused.
# for f_id in range(1, 10+1):
#     PRE_PROCESS.append(f"spec_pp({f_id})")


# Removed: 't_tol', 'alt_crv', 'regularization', 'lsq_deriv',
# Feel free to put them back if they are needed once more.
# Be sure to add them to the correct type set at the top of the file too!
SIMULATION = COMMON.copy()
SIMULATION.update({
    'run_time_info': ParamType.LOG,
    't_step_old': ParamType.INT,
    'dt': ParamType.REAL,
    't_step_start': ParamType.INT,
    't_step_stop': ParamType.INT,
    't_step_save': ParamType.INT,
    't_step_print': ParamType.INT,
    'time_stepper': ParamType.INT,
    'weno_eps': ParamType.REAL,
    'teno_CT': ParamType.REAL,
    'mapped_weno': ParamType.LOG,
    'wenoz': ParamType.LOG,
    'teno': ParamType.LOG,
    'mp_weno': ParamType.LOG,
    'weno_avg': ParamType.LOG,
    'weno_Re_flux': ParamType.LOG,
    'riemann_solver': ParamType.INT,
    'wave_speeds': ParamType.INT,
    'avg_state': ParamType.INT,
    'prim_vars_wrt': ParamType.LOG,
    'alt_soundspeed': ParamType.LOG,
    'null_weights': ParamType.LOG,
    'mixture_err': ParamType.LOG,
    'fd_order': ParamType.INT,
    'num_probes': ParamType.INT,
    'probe_wrt': ParamType.LOG,
    'bubble_model': ParamType.INT,
    'acoustic_source': ParamType.LOG,
    'num_source': ParamType.INT,
    'qbmm': ParamType.LOG,
    'R0_type': ParamType.INT,
    'integral_wrt': ParamType.LOG,
    'num_integrals': ParamType.INT,
    'rdma_mpi': ParamType.LOG,
    'palpha_eps': ParamType.REAL,
    'ptgalpha_eps': ParamType.REAL,
    'pi_fac': ParamType.REAL,
    'adap_dt': ParamType.LOG,
    'ib': ParamType.LOG,
    'num_ibs': ParamType.INT,
    'low_Mach': ParamType.INT,
})

# NOTE: Not currently present
# for var in [ 'advection', 'diffusion', 'reactions' ]:
#     SIMULATION.append(f'chem_params%{var}')

for ib_id in range(1, 10+1):
    for real_attr, ty in [("geometry", ParamType.INT), ("radius", ParamType.REAL),
                          ("theta", ParamType.REAL), ("slip", ParamType.LOG),
                          ("c", ParamType.REAL), ("p", ParamType.REAL),
                          ("t", ParamType.REAL), ("m", ParamType.REAL)]:
        SIMULATION[f"patch_ib({ib_id})%{real_attr}"] = ty

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        SIMULATION[f'patch_ib({ib_id})%{cmp}_centroid'] = ParamType.REAL
        SIMULATION[f'patch_ib({ib_id})%length_{cmp}'] = ParamType.REAL

for cmp in ["x", "y", "z"]:
    SIMULATION[f'bc_{cmp}%beg'] = ParamType.INT
    SIMULATION[f'bc_{cmp}%end'] = ParamType.INT
    SIMULATION[f'bc_{cmp}%vb1'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%vb2'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%vb3'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%ve1'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%ve2'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%ve3'] = ParamType.REAL

    for var in ["k", "w", "p", "g"]:
        SIMULATION[f'{var}_{cmp}'] = ParamType.REAL
    SIMULATION[f'bf_{cmp}'] = ParamType.LOG


    for prepend in ["domain%beg", "domain%end"]:
        SIMULATION[f"{cmp}_{prepend}"] = ParamType.REAL

# NOTE: This is now just "probe_wrt"
# for wrt_id in range(1,10+1):
#    for cmp in ["x", "y", "z"]:
#        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')
#        set_type(f'probe_wrt({wrt_id})%{cmp}', ParamType.LOG)

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION[f'probe({probe_id})%{cmp}'] = ParamType.REAL

for f_id in range(1,10+1):
    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G", "cv", "qv", "qvp" ]:
        SIMULATION[f"fluid_pp({f_id})%{real_attr}"] = ParamType.REAL

    for re_id in [1, 2]:
        SIMULATION[f"fluid_pp({f_id})%Re({re_id})"] = ParamType.REAL

    for mono_id in range(1,4+1):
        for int_attr in ["pulse", "support", "num_elements", "element_on"]:
            SIMULATION[f"acoustic({mono_id})%{int_attr}"] = ParamType.INT

        SIMULATION[f"acoustic({mono_id})%dipole"] = ParamType.LOG

        for real_attr in ["mag", "length", "height", "wavelength", "frequency",
                          "gauss_sigma_dist", "gauss_sigma_time", "npulse",
                          "dir", "delay", "foc_length", "aperture",
                          "element_spacing_angle", "element_polygon_ratio",
                          "rotate_angle"]:
            SIMULATION[f"acoustic({mono_id})%{real_attr}"] = ParamType.REAL

        for cmp_id in range(1,3+1):
            SIMULATION[f"acoustic({mono_id})%loc({cmp_id})"] = ParamType.REAL

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION[f"integral({int_id})%{cmp}min"] = ParamType.REAL
            SIMULATION[f"integral({int_id})%{cmp}max"] = ParamType.REAL


# Removed: 'fourier_modes%beg', 'fourier_modes%end', 'chem_wrt'
# Feel free to return them if they are needed once more.
POST_PROCESS = COMMON.copy()
POST_PROCESS.update({
    't_step_start': ParamType.INT,
    't_step_stop': ParamType.INT,
    't_step_save': ParamType.INT,
    'alt_soundspeed': ParamType.LOG,
    'mixture_err': ParamType.LOG,
    'format': ParamType.INT,
    'schlieren_wrt': ParamType.LOG,
    'schlieren_alpha': ParamType.REAL,
    'fd_order': ParamType.INT,
    'alpha_rho_wrt': ParamType.LOG,
    'rho_wrt': ParamType.LOG,
    'mom_wrt': ParamType.LOG,
    'vel_wrt': ParamType.LOG,
    'flux_lim': ParamType.INT,
    'flux_wrt': ParamType.LOG,
    'E_wrt': ParamType.LOG,
    'pres_wrt': ParamType.LOG,
    'alpha_wrt': ParamType.LOG,
    'kappa_wrt': ParamType.LOG,
    'gamma_wrt': ParamType.LOG,
    'heat_ratio_wrt': ParamType.LOG,
    'pi_inf_wrt': ParamType.LOG,
    'pres_inf_wrt': ParamType.LOG,
    'cons_vars_wrt': ParamType.LOG,
    'prim_vars_wrt': ParamType.LOG,
    'c_wrt': ParamType.LOG,
    'omega_wrt': ParamType.LOG,
    'qbmm': ParamType.LOG,
    'qm_wrt': ParamType.LOG,
    'cf_wrt': ParamType.LOG,
    'ib': ParamType.LOG
})

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS[f'bc_{cmp}%beg'] = ParamType.INT
    POST_PROCESS[f'bc_{cmp}%end'] = ParamType.INT

    for real_attr in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS[f'{real_attr}({cmp_id})'] = ParamType.LOG

# NOTE: `chem_wrt` is missing
# for cmp_id in range(100):
#     POST_PROCESS.append(f'chem_wrt({cmp_id})')

for fl_id in range(1,10+1):
    for append, ty in [("schlieren_alpha", ParamType.REAL),
                       ("alpha_rho_wrt", ParamType.LOG),
                       ("alpha_wrt", ParamType.LOG), ("kappa_wrt", ParamType.LOG)]:
        POST_PROCESS[f'{append}({fl_id})'] = ty

    for real_attr in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0",
                      "cv", "qv", "qvp" ]:
        POST_PROCESS[f"fluid_pp({fl_id})%{real_attr}"] = ParamType.REAL


ALL = COMMON.copy()
ALL.update(PRE_PROCESS)
ALL.update(SIMULATION)
ALL.update(POST_PROCESS)

CASE_OPTIMIZATION = [ "mapped_weno", "wenoz", "teno", "nb", "weno_order", "num_fluids" ]

_properties = { k: v.value for k, v in ALL.items() }

SCHEMA = {
    "type": "object",
    "properties": _properties
}


def get_input_dict_keys(target_name: str) -> list:
    result = {
        "pre_process"  : PRE_PROCESS,
        "simulation"   : SIMULATION,
        "post_process" : POST_PROCESS
    }.get(target_name, {}).keys()

    if not ARG("case_optimization") or target_name != "simulation":
        return result

    return [ x for x in result if x not in CASE_OPTIMIZATION ]
