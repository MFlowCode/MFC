"""
Post-process and remaining parameter definitions.

Defines parameters for post-processing output, chemistry species output,
and remaining miscellaneous parameters.
"""

from ..schema import ParamDef, ParamType, Stage
from ..registry import REGISTRY


NUM_FLUIDS = 10


def _reg(name: str, ptype: ParamType, stages: set, category: str = None):
    """Helper to register a parameter."""
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        stages=stages,
        category=category,
    ))


def register_post_process_wrt_params():
    """Register post-process output control parameters."""

    # Per-fluid output controls
    for fl_id in range(1, NUM_FLUIDS + 1):
        for attr, ptype in [("schlieren_alpha", ParamType.REAL),
                            ("alpha_rho_wrt", ParamType.LOG),
                            ("alpha_wrt", ParamType.LOG),
                            ("kappa_wrt", ParamType.LOG),
                            ("alpha_rho_e_wrt", ParamType.LOG)]:
            _reg(f"{attr}({fl_id})", ptype, {Stage.POST_PROCESS}, "post_process_wrt")

    # Component-wise output controls (1-3)
    for cmp_id in range(1, 4):
        for attr in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
            _reg(f"{attr}({cmp_id})", ParamType.LOG, {Stage.POST_PROCESS}, "post_process_wrt")

    # Partial domain output
    for cmp in ["x", "y", "z"]:
        _reg(f"{cmp}_output%beg", ParamType.REAL, {Stage.POST_PROCESS}, "post_process_output")
        _reg(f"{cmp}_output%end", ParamType.REAL, {Stage.POST_PROCESS}, "post_process_output")


def register_chem_wrt_params():
    """Register chemistry species output parameters."""

    # Y species output (0-99)
    for cmp_id in range(100):
        _reg(f"chem_wrt_Y({cmp_id})", ParamType.LOG, {Stage.POST_PROCESS}, "chem_wrt")

    # Temperature output
    _reg("chem_wrt_T", ParamType.LOG, {Stage.POST_PROCESS}, "chem_wrt")


def register_fluid_rho_params():
    """Register fluid_rho parameters for pre-process."""

    for f_id in range(1, NUM_FLUIDS + 1):
        _reg(f"fluid_rho({f_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "fluid_rho")


def register_remaining_common_params():
    """Register remaining COMMON parameters not yet defined."""

    # MHD and physics models
    _reg("mhd", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("hypoelasticity", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("hyperelasticity", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("pref", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("parallel_io", ParamType.LOG, {Stage.COMMON}, "io")
    _reg("poly_sigma", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("case_dir", ParamType.STR, {Stage.COMMON}, "io")
    _reg("thermal", ParamType.INT, {Stage.COMMON}, "physics")
    _reg("polytropic", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("mpp_lim", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("R0ref", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("nb", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("rhoref", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("bubbles_euler", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("precision", ParamType.INT, {Stage.COMMON}, "numerics")
    _reg("polydisperse", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("file_per_process", ParamType.LOG, {Stage.COMMON}, "io")
    _reg("relax", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("relax_model", ParamType.INT, {Stage.COMMON}, "physics")
    _reg("sigma", ParamType.REAL, {Stage.COMMON}, "physics")
    _reg("adv_n", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("cfl_adap_dt", ParamType.LOG, {Stage.COMMON}, "time_stepping")
    _reg("cfl_const_dt", ParamType.LOG, {Stage.COMMON}, "time_stepping")
    _reg("chemistry", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("cantera_file", ParamType.STR, {Stage.COMMON}, "physics")
    _reg("Bx0", ParamType.REAL, {Stage.COMMON}, "mhd")
    _reg("relativity", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("cont_damage", ParamType.LOG, {Stage.COMMON}, "physics")
    _reg("num_bc_patches", ParamType.INT, {Stage.COMMON}, "boundary_conditions")
    _reg("igr", ParamType.LOG, {Stage.COMMON}, "numerics")
    _reg("igr_order", ParamType.INT, {Stage.COMMON}, "numerics")
    _reg("down_sample", ParamType.LOG, {Stage.COMMON}, "io")
    _reg("recon_type", ParamType.INT, {Stage.COMMON}, "numerics")
    _reg("muscl_order", ParamType.INT, {Stage.COMMON}, "numerics")


def register_remaining_pre_process_params():
    """Register remaining PRE_PROCESS-only parameters."""

    _reg("old_grid", ParamType.LOG, {Stage.PRE_PROCESS}, "restart")
    _reg("old_ic", ParamType.LOG, {Stage.PRE_PROCESS}, "restart")
    _reg("t_step_old", ParamType.INT, {Stage.PRE_PROCESS}, "restart")
    _reg("mixlayer_vel_profile", ParamType.LOG, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("mixlayer_vel_coef", ParamType.REAL, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("mixlayer_domain", ParamType.REAL, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("mixlayer_perturb", ParamType.LOG, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("mixlayer_perturb_nk", ParamType.INT, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("mixlayer_perturb_k0", ParamType.REAL, {Stage.PRE_PROCESS}, "mixlayer")
    _reg("perturb_flow", ParamType.LOG, {Stage.PRE_PROCESS}, "perturb")
    _reg("perturb_flow_fluid", ParamType.INT, {Stage.PRE_PROCESS}, "perturb")
    _reg("perturb_flow_mag", ParamType.REAL, {Stage.PRE_PROCESS}, "perturb")
    _reg("perturb_sph", ParamType.LOG, {Stage.PRE_PROCESS}, "perturb")
    _reg("perturb_sph_fluid", ParamType.INT, {Stage.PRE_PROCESS}, "perturb")
    _reg("fluid_rho", ParamType.REAL, {Stage.PRE_PROCESS}, "physics")
    _reg("qbmm", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "physics")
    _reg("dist_type", ParamType.INT, {Stage.PRE_PROCESS}, "physics")
    _reg("sigR", ParamType.REAL, {Stage.PRE_PROCESS}, "physics")
    _reg("sigV", ParamType.REAL, {Stage.PRE_PROCESS}, "physics")
    _reg("rhoRV", ParamType.REAL, {Stage.PRE_PROCESS}, "physics")
    _reg("palpha_eps", ParamType.REAL, {Stage.PRE_PROCESS, Stage.SIMULATION}, "physics")
    _reg("ptgalpha_eps", ParamType.REAL, {Stage.PRE_PROCESS, Stage.SIMULATION}, "physics")
    _reg("pi_fac", ParamType.REAL, {Stage.PRE_PROCESS, Stage.SIMULATION}, "physics")
    _reg("ib", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "ibm")
    _reg("num_ibs", ParamType.INT, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "ibm")
    _reg("pre_stress", ParamType.LOG, {Stage.PRE_PROCESS}, "physics")
    _reg("cfl_dt", ParamType.LOG, {Stage.PRE_PROCESS}, "time_stepping")
    _reg("n_start", ParamType.INT, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "time_stepping")
    _reg("n_start_old", ParamType.INT, {Stage.PRE_PROCESS}, "restart")
    _reg("surface_tension", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "physics")
    _reg("elliptic_smoothing", ParamType.LOG, {Stage.PRE_PROCESS}, "mesh")
    _reg("elliptic_smoothing_iters", ParamType.INT, {Stage.PRE_PROCESS}, "mesh")
    _reg("viscous", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION}, "physics")
    _reg("bubbles_lagrange", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "physics")
    _reg("simplex_perturb", ParamType.LOG, {Stage.PRE_PROCESS}, "perturb")
    _reg("fft_wrt", ParamType.LOG, {Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS}, "io")

    # Grid stretching parameters (x_a, x_b, y_a, etc.)
    for cmp in ["x", "y", "z"]:
        for prepend in ["a", "b"]:
            _reg(f"{cmp}_{prepend}", ParamType.REAL, {Stage.PRE_PROCESS}, "grid_stretch")
        _reg(f"a_{cmp}", ParamType.REAL, {Stage.PRE_PROCESS}, "grid_stretch")
        _reg(f"loops_{cmp}", ParamType.INT, {Stage.PRE_PROCESS}, "grid_stretch")


def register_remaining_simulation_params():
    """Register remaining SIMULATION-only parameters."""

    _reg("run_time_info", ParamType.LOG, {Stage.SIMULATION}, "io")
    _reg("t_step_old", ParamType.INT, {Stage.SIMULATION}, "restart")
    _reg("t_step_print", ParamType.INT, {Stage.SIMULATION}, "io")
    _reg("teno_CT", ParamType.REAL, {Stage.SIMULATION}, "weno")
    _reg("wenoz_q", ParamType.REAL, {Stage.SIMULATION}, "weno")
    _reg("mp_weno", ParamType.LOG, {Stage.SIMULATION}, "weno")
    _reg("riemann_solver", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("wave_speeds", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("avg_state", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("prim_vars_wrt", ParamType.LOG, {Stage.SIMULATION, Stage.POST_PROCESS}, "io")
    _reg("alt_soundspeed", ParamType.LOG, {Stage.SIMULATION, Stage.POST_PROCESS}, "physics")
    _reg("null_weights", ParamType.LOG, {Stage.SIMULATION}, "numerics")
    _reg("mixture_err", ParamType.LOG, {Stage.SIMULATION, Stage.POST_PROCESS}, "physics")
    _reg("fd_order", ParamType.INT, {Stage.SIMULATION, Stage.POST_PROCESS}, "numerics")
    _reg("num_probes", ParamType.INT, {Stage.SIMULATION}, "io")
    _reg("probe_wrt", ParamType.LOG, {Stage.SIMULATION}, "io")
    _reg("bubble_model", ParamType.INT, {Stage.SIMULATION}, "physics")
    _reg("acoustic_source", ParamType.LOG, {Stage.SIMULATION}, "physics")
    _reg("num_source", ParamType.INT, {Stage.SIMULATION}, "physics")
    _reg("integral_wrt", ParamType.LOG, {Stage.SIMULATION}, "io")
    _reg("num_integrals", ParamType.INT, {Stage.SIMULATION}, "io")
    _reg("rdma_mpi", ParamType.LOG, {Stage.SIMULATION}, "parallel")
    _reg("adap_dt", ParamType.LOG, {Stage.SIMULATION}, "time_stepping")
    _reg("adap_dt_tol", ParamType.REAL, {Stage.SIMULATION}, "time_stepping")
    _reg("adap_dt_max_iters", ParamType.INT, {Stage.SIMULATION}, "time_stepping")
    _reg("t_stop", ParamType.REAL, {Stage.SIMULATION, Stage.POST_PROCESS}, "time_stepping")
    _reg("t_save", ParamType.REAL, {Stage.SIMULATION, Stage.POST_PROCESS}, "time_stepping")
    _reg("low_Mach", ParamType.INT, {Stage.SIMULATION}, "physics")
    _reg("powell", ParamType.LOG, {Stage.SIMULATION}, "mhd")
    _reg("tau_star", ParamType.REAL, {Stage.SIMULATION}, "physics")
    _reg("cont_damage_s", ParamType.REAL, {Stage.SIMULATION}, "physics")
    _reg("alpha_bar", ParamType.REAL, {Stage.SIMULATION}, "physics")
    _reg("num_igr_iters", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("num_igr_warm_start_iters", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("alf_factor", ParamType.REAL, {Stage.SIMULATION}, "numerics")
    _reg("igr_iter_solver", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("igr_pres_lim", ParamType.LOG, {Stage.SIMULATION}, "numerics")
    _reg("muscl_lim", ParamType.INT, {Stage.SIMULATION}, "numerics")
    _reg("int_comp", ParamType.LOG, {Stage.SIMULATION}, "physics")
    _reg("ic_eps", ParamType.REAL, {Stage.SIMULATION}, "physics")
    _reg("ic_beta", ParamType.REAL, {Stage.SIMULATION}, "physics")
    _reg("nv_uvm_out_of_core", ParamType.LOG, {Stage.SIMULATION}, "gpu")
    _reg("nv_uvm_igr_temps_on_gpu", ParamType.INT, {Stage.SIMULATION}, "gpu")
    _reg("nv_uvm_pref_gpu", ParamType.LOG, {Stage.SIMULATION}, "gpu")


def register_remaining_post_process_params():
    """Register remaining POST_PROCESS-only parameters."""

    _reg("format", ParamType.INT, {Stage.POST_PROCESS}, "io")
    _reg("schlieren_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("schlieren_alpha", ParamType.REAL, {Stage.POST_PROCESS}, "io")
    _reg("alpha_rho_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("rho_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("mom_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("vel_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("flux_lim", ParamType.INT, {Stage.POST_PROCESS}, "io")
    _reg("flux_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("E_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("pres_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("alpha_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("kappa_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("gamma_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("heat_ratio_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("pi_inf_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("pres_inf_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("cons_vars_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("c_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("omega_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("qm_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("liutex_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("cf_wrt", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("sim_data", ParamType.LOG, {Stage.POST_PROCESS}, "io")
    _reg("output_partial_domain", ParamType.LOG, {Stage.POST_PROCESS}, "io")

    # Lagrangian bubble post-process output
    for var in ["lag_header", "lag_txt_wrt", "lag_db_wrt", "lag_id_wrt", "lag_pos_wrt",
                "lag_pos_prev_wrt", "lag_vel_wrt", "lag_rad_wrt", "lag_rvel_wrt",
                "lag_r0_wrt", "lag_rmax_wrt", "lag_rmin_wrt", "lag_dphidt_wrt",
                "lag_pres_wrt", "lag_mv_wrt", "lag_mg_wrt", "lag_betaT_wrt", "lag_betaC_wrt"]:
        _reg(var, ParamType.LOG, {Stage.POST_PROCESS}, "lag_post")


# Auto-register when module is imported
register_post_process_wrt_params()
register_chem_wrt_params()
register_fluid_rho_params()
register_remaining_common_params()
register_remaining_pre_process_params()
register_remaining_simulation_params()
register_remaining_post_process_params()
