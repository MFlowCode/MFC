import fastjsonschema

from enum import Enum
from ..state import ARG
from functools import cache


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
    'mhd': ParamType.LOG,
    'hypoelasticity': ParamType.LOG,
    'hyperelasticity': ParamType.LOG,
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
    'bubbles_euler': ParamType.LOG,
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
    'cfl_adap_dt': ParamType.LOG,
    'cfl_const_dt': ParamType.LOG,
    'chemistry': ParamType.LOG,
    'cantera_file': ParamType.STR,
    'Bx0': ParamType.REAL,
    'relativity': ParamType.LOG,
    'cont_damage': ParamType.LOG,
    'num_bc_patches': ParamType.INT,
    'igr': ParamType.LOG,
    'igr_order': ParamType.INT,
    'recon_type': ParamType.INT,
    'muscl_order': ParamType.INT,
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
    'mixlayer_perturb_nk': ParamType.INT,
    'mixlayer_perturb_k0': ParamType.REAL,
    'perturb_flow': ParamType.LOG,
    'perturb_flow_fluid': ParamType.INT,
    'perturb_flow_mag': ParamType.REAL,
    'perturb_sph': ParamType.LOG,
    'perturb_sph_fluid': ParamType.INT,
    'fluid_rho': ParamType.REAL,
    'num_patches': ParamType.INT,
    'qbmm': ParamType.LOG,
    'dist_type': ParamType.INT,
    'sigR': ParamType.REAL,
    'sigV': ParamType.REAL,
    'rhoRV': ParamType.REAL,
    'palpha_eps': ParamType.REAL,
    'ptgalpha_eps': ParamType.REAL,
    'pi_fac': ParamType.REAL,
    'ib': ParamType.LOG,
    'num_ibs': ParamType.INT,
    'pre_stress': ParamType.LOG,
    'cfl_dt': ParamType.LOG,
    'n_start': ParamType.INT,
    'n_start_old': ParamType.INT,
    'surface_tension': ParamType.LOG,
    'elliptic_smoothing': ParamType.LOG,
    'elliptic_smoothing_iters': ParamType.INT,
    'viscous': ParamType.LOG,
    'bubbles_lagrange': ParamType.LOG,
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

    for real_attr_stl, ty_stl in [("filepath", ParamType.STR), ("spc", ParamType.INT),
                          ("threshold", ParamType.REAL)]:
        PRE_PROCESS[f"patch_ib({ib_id})%model_{real_attr_stl}"] = ty_stl

    for real_attr_stl2 in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS[f"patch_ib({ib_id})%model_{real_attr_stl2}({j})"] = ParamType.REAL

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
                      "mu_v", "k_v", "cp_v", "G", "cv", "qv", "qvp" ]:
        PRE_PROCESS[f"fluid_pp({f_id})%{real_attr}"] = ParamType.REAL

for bc_p_id in range(1, 10+1):
    for attribute in ["geometry","type","dir","loc"]:
        PRE_PROCESS[f"patch_bc({bc_p_id})%{attribute}"] = ParamType.INT

    for attribute in ["centroid","length"]:
        for d_id in range(1, 3+1):
            PRE_PROCESS[f"patch_bc({bc_p_id})%{attribute}({d_id})"] = ParamType.REAL

    PRE_PROCESS[f"patch_bc({bc_p_id})%radius"] = ParamType.REAL

for p_id in range(1, 10+1):
    for attribute, ty in [("geometry", ParamType.INT), ("smoothen", ParamType.LOG),
                      ("smooth_patch_id", ParamType.INT), ("hcid", ParamType.INT)]:
        PRE_PROCESS[f"patch_icpp({p_id})%{attribute}"] = ty

    for real_attr in ["radius",  "radii", "epsilon", "beta", "normal", "alpha_rho",
                      'non_axis_sym', "normal", "smooth_coeff", "rho", "vel",
                      "alpha", "gamma", "pi_inf", "r0", "v0", "p0", "m0", "cv",
                      "qv", "qvp"]:
        PRE_PROCESS[f"patch_icpp({p_id})%{real_attr}"] = ParamType.REAL

    for real_attr in range(2, 9+1):
        PRE_PROCESS[f"patch_icpp({p_id})%a({real_attr})"] = ParamType.REAL

    PRE_PROCESS[f"patch_icpp({p_id})%pres"] = ParamType.REAL.analytic()

    PRE_PROCESS[f"patch_icpp({p_id})%Bx"] = ParamType.REAL.analytic()
    PRE_PROCESS[f"patch_icpp({p_id})%By"] = ParamType.REAL.analytic()
    PRE_PROCESS[f"patch_icpp({p_id})%Bz"] = ParamType.REAL.analytic()

    for i in range(100):
        PRE_PROCESS[f"patch_icpp({p_id})%Y({i})"] = ParamType.REAL.analytic()

    PRE_PROCESS[f"patch_icpp({p_id})%model_filepath"] = ParamType.STR

    for real_attr in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            PRE_PROCESS[f"patch_icpp({p_id})%model_{real_attr}({j})"] = ParamType.REAL

    PRE_PROCESS[f"patch_icpp({p_id})%model_spc"] = ParamType.INT
    PRE_PROCESS[f"patch_icpp({p_id})%model_threshold"] = ParamType.REAL

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

    PRE_PROCESS[f'patch_icpp({p_id})%cf_val'] = ParamType.REAL.analytic()

    if p_id >= 2:
        PRE_PROCESS[f'patch_icpp({p_id})%alter_patch'] = ParamType.LOG

        for alter_id in range(1, p_id):
            PRE_PROCESS[f'patch_icpp({p_id})%alter_patch({alter_id})'] = ParamType.LOG

    PRE_PROCESS[f'patch_icpp({p_id})%cf_val'] = ParamType.REAL.analytic()

    for cmp in ["x", "y", "z"]:
        PRE_PROCESS[f'bc_{cmp}%beg'] = ParamType.INT
        PRE_PROCESS[f'bc_{cmp}%end'] = ParamType.INT
        PRE_PROCESS[f'bc_{cmp}%vb1'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%vb2'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%vb3'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%ve1'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%ve2'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%ve3'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%pres_in'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%pres_out'] = ParamType.REAL
        PRE_PROCESS[f'bc_{cmp}%grcbc_in'] = ParamType.LOG
        PRE_PROCESS[f'bc_{cmp}%grcbc_out'] = ParamType.LOG
        PRE_PROCESS[f'bc_{cmp}%grcbc_vel_out'] = ParamType.LOG

        for int_id in range(1, 10+1):
            PRE_PROCESS[f"bc_{cmp}%alpha_rho_in({int_id})"] = ParamType.REAL
            PRE_PROCESS[f"bc_{cmp}%alpha_in({int_id})"] = ParamType.REAL

        for int_id in range(1, 3+1):
            PRE_PROCESS[f"bc_{cmp}%vel_in({int_id})"] = ParamType.REAL
            PRE_PROCESS[f"bc_{cmp}%vel_out({int_id})"] = ParamType.REAL

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
    'wenoz_q': ParamType.REAL,
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
    'integral_wrt': ParamType.LOG,
    'num_integrals': ParamType.INT,
    'rdma_mpi': ParamType.LOG,
    'palpha_eps': ParamType.REAL,
    'ptgalpha_eps': ParamType.REAL,
    'pi_fac': ParamType.REAL,
    'adap_dt': ParamType.LOG,
    'adap_dt_tol': ParamType.REAL,
    'adap_dt_max_iters': ParamType.INT,
    'ib': ParamType.LOG,
    'num_ibs': ParamType.INT,
    'n_start': ParamType.INT,
    't_stop': ParamType.REAL,
    't_save': ParamType.REAL,
    'cfl_target': ParamType.REAL,
    'low_Mach': ParamType.INT,
    'surface_tension': ParamType.LOG,
    'viscous': ParamType.LOG,
    'bubbles_lagrange': ParamType.LOG,
    'num_bc_patches': ParamType.INT,
    'powell': ParamType.LOG,
    'tau_star': ParamType.REAL,
    'cont_damage_s': ParamType.REAL,
    'alpha_bar': ParamType.REAL,
    'num_igr_iters': ParamType.INT,
    'num_igr_warm_start_iters': ParamType.INT,
    'alf_factor': ParamType.REAL,
    'igr_iter_solver': ParamType.INT,
    'igr_pres_lim': ParamType.LOG,
    'recon_type': ParamType.INT,
    'muscl_order': ParamType.INT,
    'muscl_lim': ParamType.INT,
    'int_comp': ParamType.LOG,
    'ic_eps': ParamType.REAL,
    'ic_beta': ParamType.REAL,
})

for var in [ 'heatTransfer_model', 'massTransfer_model', 'pressure_corrector',
             'write_bubbles', 'write_bubbles_stats' ]:
    SIMULATION[f'lag_params%{var}'] = ParamType.LOG

for var in [ 'solver_approach', 'cluster_type', 'smooth_type', 'nBubs_glb']:
    SIMULATION[f'lag_params%{var}'] = ParamType.INT

for var in [ 'epsilonb', 'valmaxvoid', 'charwidth', 'diffcoefvap',
            'c0', 'rho0', 'T0', 'x0', 'Thost' ]:
    SIMULATION[f'lag_params%{var}'] = ParamType.REAL

for var in [ 'diffusion', 'reactions' ]:
    SIMULATION[f'chem_params%{var}'] = ParamType.LOG

for var in [ 'gamma_method' ]:
    SIMULATION[f'chem_params%{var}'] = ParamType.INT

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
    SIMULATION[f'bc_{cmp}%pres_in'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%pres_out'] = ParamType.REAL
    SIMULATION[f'bc_{cmp}%grcbc_in'] = ParamType.LOG
    SIMULATION[f'bc_{cmp}%grcbc_out'] = ParamType.LOG
    SIMULATION[f'bc_{cmp}%grcbc_vel_out'] = ParamType.LOG

    for int_id in range(1, 10+1):
        SIMULATION[f"bc_{cmp}%alpha_rho_in({int_id})"] = ParamType.REAL
        SIMULATION[f"bc_{cmp}%alpha_in({int_id})"] = ParamType.REAL

    for int_id in range(1, 3+1):
        SIMULATION[f"bc_{cmp}%vel_in({int_id})"] = ParamType.REAL
        SIMULATION[f"bc_{cmp}%vel_out({int_id})"] = ParamType.REAL

    for var in ["k", "w", "p", "g"]:
        SIMULATION[f'{var}_{cmp}'] = ParamType.REAL
    SIMULATION[f'bf_{cmp}'] = ParamType.LOG


    for prepend in ["domain%beg", "domain%end"]:
        SIMULATION[f"{cmp}_{prepend}"] = ParamType.REAL

for probe_id in range(1,10+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION[f'probe({probe_id})%{cmp}'] = ParamType.REAL

for f_id in range(1,10+1):
    for real_attr in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "cp_v", "G", "cv", "qv", "qvp" ]:
        SIMULATION[f"fluid_pp({f_id})%{real_attr}"] = ParamType.REAL

    for re_id in [1, 2]:
        SIMULATION[f"fluid_pp({f_id})%Re({re_id})"] = ParamType.REAL

    for mono_id in range(1,4+1):
        for int_attr in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            SIMULATION[f"acoustic({mono_id})%{int_attr}"] = ParamType.INT

        SIMULATION[f"acoustic({mono_id})%dipole"] = ParamType.LOG

        for real_attr in ["mag", "length", "height", "wavelength", "frequency",
                          "gauss_sigma_dist", "gauss_sigma_time", "npulse",
                          "dir", "delay", "foc_length", "aperture",
                          "element_spacing_angle", "element_polygon_ratio",
                          "rotate_angle", "bb_bandwidth", "bb_lowest_freq"]:
            SIMULATION[f"acoustic({mono_id})%{real_attr}"] = ParamType.REAL

        for cmp_id in range(1,3+1):
            SIMULATION[f"acoustic({mono_id})%loc({cmp_id})"] = ParamType.REAL

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION[f"integral({int_id})%{cmp}min"] = ParamType.REAL
            SIMULATION[f"integral({int_id})%{cmp}max"] = ParamType.REAL

# Removed: 'fourier_modes%beg', 'fourier_modes%end'.
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
    'sim_data': ParamType.LOG,
    'ib': ParamType.LOG,
    'num_ibs': ParamType.INT,
    'cfl_target': ParamType.REAL,
    't_save': ParamType.REAL,
    't_stop': ParamType.REAL,
    'n_start': ParamType.INT,
    'surface_tension': ParamType.LOG,
    'output_partial_domain': ParamType.LOG,
    'bubbles_lagrange': ParamType.LOG,
})

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS[f'bc_{cmp}%beg'] = ParamType.INT
    POST_PROCESS[f'bc_{cmp}%end'] = ParamType.INT

    POST_PROCESS[f'{cmp}_output%beg'] = ParamType.REAL
    POST_PROCESS[f'{cmp}_output%end'] = ParamType.REAL

    for real_attr in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS[f'{real_attr}({cmp_id})'] = ParamType.LOG

for cmp_id in range(100):
    POST_PROCESS[f'chem_wrt_Y({cmp_id})'] = ParamType.LOG
POST_PROCESS['chem_wrt_T'] = ParamType.LOG

for fl_id in range(1,10+1):
    for append, ty in [("schlieren_alpha", ParamType.REAL),
                       ("alpha_rho_wrt", ParamType.LOG),
                       ("alpha_wrt", ParamType.LOG), ("kappa_wrt", ParamType.LOG)]:
        POST_PROCESS[f'{append}({fl_id})'] = ty

    for real_attr in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "cp_v",
                      "G", "mul0", "cv", "qv", "qvp" ]:
        POST_PROCESS[f"fluid_pp({fl_id})%{real_attr}"] = ParamType.REAL

IGNORE = ["cantera_file", "chemistry"]

ALL = COMMON.copy()
ALL.update(PRE_PROCESS)
ALL.update(SIMULATION)
ALL.update(POST_PROCESS)

CASE_OPTIMIZATION = [ "mapped_weno", "wenoz", "teno", "wenoz_q", "nb", "weno_order",
                     "num_fluids", "mhd", "relativity", "igr_order", "viscous",
                     "igr_iter_solver", "igr", "igr_pres_lim", "recon_type", "muscl_order", "muscl_lim" ]

_properties = { k: v.value for k, v in ALL.items() }

SCHEMA = {
    "type": "object",
    "properties": _properties,
    "additionalProperties": False
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


@cache
def get_validator():
    return fastjsonschema.compile(SCHEMA)
