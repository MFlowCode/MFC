"""
Parse Fortran namelist definitions to extract valid parameters for each target.

This module reads the Fortran source files and extracts the parameter names
from each target's namelist definition. This ensures the Python toolchain
stays in sync with what the Fortran code actually accepts.

When Fortran sources are unavailable (e.g. Homebrew installs), a built-in
fallback parameter set is used instead.
"""

import re
from pathlib import Path
from typing import Dict, Set


# Fallback parameters for when Fortran source files are not available.
# Generated from the namelist definitions in src/*/m_start_up.fpp.
# To regenerate: python3 toolchain/mfc/params/namelist_parser.py
_FALLBACK_PARAMS = {
    'pre_process': {
        'Bx0', 'Ca', 'R0ref', 'Re_inv', 'Web',
        'a_x', 'a_y', 'a_z', 'adv_n', 'bc_x',
        'bc_y', 'bc_z', 'bub_pp', 'bubbles_euler', 'bubbles_lagrange',
        'case_dir', 'cfl_adap_dt', 'cfl_const_dt', 'cont_damage', 'cyl_coord',
        'dist_type', 'down_sample', 'elliptic_smoothing', 'elliptic_smoothing_iters', 'fft_wrt',
        'file_per_process', 'fluid_pp', 'fluid_rho', 'hyper_cleaning', 'hyperelasticity',
        'hypoelasticity', 'ib', 'igr', 'igr_order', 'loops_x',
        'loops_y', 'loops_z', 'm', 'mhd', 'mixlayer_perturb',
        'mixlayer_perturb_k0', 'mixlayer_perturb_nk', 'mixlayer_vel_coef', 'mixlayer_vel_profile', 'model_eqns',
        'mpp_lim', 'muscl_order', 'n', 'n_start', 'n_start_old',
        'nb', 'num_bc_patches', 'num_fluids', 'num_ibs', 'num_patches',
        'old_grid', 'old_ic', 'p', 'palpha_eps', 'parallel_io',
        'patch_bc', 'patch_ib', 'patch_icpp', 'perturb_flow', 'perturb_flow_fluid',
        'perturb_flow_mag', 'perturb_sph', 'perturb_sph_fluid', 'pi_fac', 'poly_sigma',
        'polydisperse', 'polytropic', 'pre_stress', 'precision', 'pref',
        'ptgalpha_eps', 'qbmm', 'recon_type', 'relativity', 'relax',
        'relax_model', 'rhoRV', 'rhoref', 'sigR', 'sigV',
        'sigma', 'simplex_params', 'simplex_perturb', 'stretch_x', 'stretch_y',
        'stretch_z', 'surface_tension', 't_step_old', 't_step_start', 'thermal',
        'viscous', 'weno_order', 'x_a', 'x_b', 'x_domain',
        'y_a', 'y_b', 'y_domain', 'z_a', 'z_b',
        'z_domain',
    },
    'simulation': {
        'Bx0', 'Ca', 'R0ref', 'Re_inv', 'Web',
        'acoustic', 'acoustic_source', 'adap_dt', 'adap_dt_max_iters', 'adap_dt_tol',
        'adv_n', 'alf_factor', 'alpha_bar', 'alt_soundspeed', 'avg_state',
        'bc_x', 'bc_y', 'bc_z', 'bf_x', 'bf_y',
        'bf_z', 'bub_pp', 'bubble_model', 'bubbles_euler', 'bubbles_lagrange',
        'case_dir', 'cfl_adap_dt', 'cfl_const_dt', 'cfl_target', 'chem_params',
        'cont_damage', 'cont_damage_s', 'cyl_coord', 'down_sample', 'dt',
        'fd_order', 'fft_wrt', 'file_per_process', 'fluid_pp', 'g_x',
        'g_y', 'g_z', 'hyper_cleaning', 'hyper_cleaning_speed', 'hyper_cleaning_tau',
        'hyperelasticity', 'hypoelasticity', 'ib', 'ic_beta', 'ic_eps',
        'igr', 'igr_iter_solver', 'igr_order', 'igr_pres_lim', 'int_comp',
        'integral', 'integral_wrt', 'k_x', 'k_y', 'k_z',
        'lag_params', 'low_Mach', 'm', 'mapped_weno', 'mhd',
        'mixture_err', 'model_eqns', 'mp_weno', 'mpp_lim', 'muscl_lim',
        'muscl_order', 'n', 'n_start', 'nb', 'null_weights',
        'num_bc_patches', 'num_fluids', 'num_ibs', 'num_igr_iters', 'num_igr_warm_start_iters',
        'num_integrals', 'num_probes', 'num_source', 'nv_uvm_igr_temps_on_gpu', 'nv_uvm_out_of_core',
        'nv_uvm_pref_gpu', 'p', 'p_x', 'p_y', 'p_z',
        'palpha_eps', 'parallel_io', 'patch_ib', 'pi_fac', 'poly_sigma',
        'polydisperse', 'polytropic', 'precision', 'pref', 'prim_vars_wrt',
        'probe', 'probe_wrt', 'ptgalpha_eps', 'qbmm', 'rdma_mpi',
        'recon_type', 'relativity', 'relax', 'relax_model', 'rhoref',
        'riemann_solver', 'run_time_info', 'sigma', 'surface_tension', 't_save',
        't_step_old', 't_step_print', 't_step_save', 't_step_start', 't_step_stop',
        't_stop', 'tau_star', 'teno', 'teno_CT', 'thermal',
        'time_stepper', 'viscous', 'w_x', 'w_y', 'w_z',
        'wave_speeds', 'weno_Re_flux', 'weno_avg', 'weno_eps', 'weno_order',
        'wenoz', 'wenoz_q', 'x_a', 'x_b', 'x_domain',
        'y_a', 'y_b', 'y_domain', 'z_a', 'z_b',
        'z_domain',
    },
    'post_process': {
        'Bx0', 'Ca', 'E_wrt', 'G', 'R0ref',
        'Re_inv', 'Web', 'adv_n', 'alpha_rho_e_wrt', 'alpha_rho_wrt',
        'alpha_wrt', 'alt_soundspeed', 'avg_state', 'bc_x', 'bc_y',
        'bc_z', 'bub_pp', 'bubbles_euler', 'bubbles_lagrange', 'c_wrt',
        'case_dir', 'cf_wrt', 'cfl_adap_dt', 'cfl_const_dt', 'cfl_target',
        'chem_wrt_T', 'chem_wrt_Y', 'cons_vars_wrt', 'cont_damage', 'cyl_coord',
        'down_sample', 'fd_order', 'fft_wrt', 'file_per_process', 'fluid_pp',
        'flux_lim', 'flux_wrt', 'format', 'gamma_wrt', 'heat_ratio_wrt',
        'hyper_cleaning', 'hyperelasticity', 'hypoelasticity', 'ib', 'igr',
        'igr_order', 'lag_betaC_wrt', 'lag_betaT_wrt', 'lag_db_wrt', 'lag_dphidt_wrt',
        'lag_header', 'lag_id_wrt', 'lag_mg_wrt', 'lag_mv_wrt', 'lag_pos_prev_wrt',
        'lag_pos_wrt', 'lag_pres_wrt', 'lag_r0_wrt', 'lag_rad_wrt', 'lag_rmax_wrt',
        'lag_rmin_wrt', 'lag_rvel_wrt', 'lag_txt_wrt', 'lag_vel_wrt', 'liutex_wrt',
        'm', 'mhd', 'mixture_err', 'model_eqns', 'mom_wrt',
        'mpp_lim', 'muscl_order', 'n', 'n_start', 'nb',
        'num_bc_patches', 'num_fluids', 'num_ibs', 'omega_wrt', 'output_partial_domain',
        'p', 'parallel_io', 'pi_inf_wrt', 'poly_sigma', 'polydisperse',
        'polytropic', 'precision', 'pref', 'pres_inf_wrt', 'pres_wrt',
        'prim_vars_wrt', 'qbmm', 'qm_wrt', 'recon_type', 'relativity',
        'relax', 'relax_model', 'rho_wrt', 'rhoref', 'schlieren_alpha',
        'schlieren_wrt', 'sigR', 'sigma', 'sim_data', 'surface_tension',
        't_save', 't_step_save', 't_step_start', 't_step_stop', 't_stop',
        'thermal', 'vel_wrt', 'weno_order', 'x_output', 'y_output',
        'z_output',
    },
}


def parse_namelist_from_file(filepath: Path) -> Set[str]:
    """
    Parse a Fortran file and extract parameter names from the namelist definition.

    Args:
        filepath: Path to the Fortran source file (m_start_up.fpp)

    Returns:
        Set of parameter names found in the namelist
    """
    content = filepath.read_text()

    # Find the namelist block - starts with "namelist /user_inputs/"
    # and continues until a line without continuation (&) or a blank line
    namelist_match = re.search(
        r'namelist\s+/user_inputs/\s*(.+?)(?=\n\s*\n|\n\s*!(?!\s*&)|\n\s*[a-zA-Z_]+\s*=)',
        content,
        re.DOTALL | re.IGNORECASE
    )

    if not namelist_match:
        raise ValueError(f"Could not find namelist /user_inputs/ in {filepath}")

    namelist_text = namelist_match.group(1)

    # Remove Fortran line continuations (&) and join lines
    namelist_text = re.sub(r'&\s*\n\s*', ' ', namelist_text)

    # Remove preprocessor directives (#:if, #:endif, etc.)
    namelist_text = re.sub(r'#:.*', '', namelist_text)

    # Remove comments (! to end of line, but not inside strings)
    namelist_text = re.sub(r'!.*', '', namelist_text)

    # Extract parameter names - they're comma-separated identifiers
    # Parameter names are alphanumeric with underscores
    found_params = set()
    for match in re.finditer(r'\b([a-zA-Z_][a-zA-Z0-9_]*)\b', namelist_text):
        name = match.group(1)
        # Skip Fortran keywords that might appear
        if name.lower() not in {'namelist', 'user_inputs', 'if', 'endif', 'not'}:
            found_params.add(name)

    return found_params


def parse_all_namelists(mfc_root: Path) -> Dict[str, Set[str]]:
    """
    Parse namelist definitions from all MFC targets.

    Args:
        mfc_root: Path to MFC root directory

    Returns:
        Dict mapping target name to set of valid parameter names.
        Falls back to built-in parameter sets when sources are unavailable.
    """
    targets = {
        'pre_process': mfc_root / 'src' / 'pre_process' / 'm_start_up.fpp',
        'simulation': mfc_root / 'src' / 'simulation' / 'm_start_up.fpp',
        'post_process': mfc_root / 'src' / 'post_process' / 'm_start_up.fpp',
    }

    result = {}
    for target_name, filepath in targets.items():
        if not filepath.exists():
            # Source files not available (e.g. Homebrew install).
            # Use built-in fallback parameters.
            return dict(_FALLBACK_PARAMS)
        result[target_name] = parse_namelist_from_file(filepath)

    return result


def get_mfc_root() -> Path:
    """Get the MFC root directory from this file's location."""
    # This file is at toolchain/mfc/params/namelist_parser.py
    # MFC root is 4 levels up
    return Path(__file__).resolve().parent.parent.parent.parent


# Module-level cache for parsed target params
_TARGET_PARAMS_CACHE: Dict[str, Set[str]] = {}


def get_target_params() -> Dict[str, Set[str]]:
    """
    Get the valid parameters for each target, parsing Fortran if needed.

    Returns:
        Dict mapping target name to set of valid parameter names
    """
    if not _TARGET_PARAMS_CACHE:
        _TARGET_PARAMS_CACHE.update(parse_all_namelists(get_mfc_root()))
    return _TARGET_PARAMS_CACHE


def is_param_valid_for_target(param_name: str, target_name: str) -> bool:
    """
    Check if a parameter is valid for a given target.

    This handles both scalar params (like "m") and indexed params
    (like "patch_icpp(1)%geometry") by checking the base name.

    Args:
        param_name: The parameter name (may include indices like "(1)%attr")
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        True if the parameter is valid for the target
    """
    valid_params = get_target_params().get(target_name, set())

    # Extract base parameter name (before any index or attribute)
    # e.g., "patch_icpp(1)%geometry" -> "patch_icpp"
    # e.g., "fluid_pp(2)%gamma" -> "fluid_pp"
    base_match = re.match(r'^([a-zA-Z_][a-zA-Z0-9_]*)', param_name)
    if base_match:
        return base_match.group(1) in valid_params

    return param_name in valid_params


if __name__ == '__main__':
    # Test the parser
    import sys

    try:
        parsed_targets = parse_all_namelists(get_mfc_root())

        print("Parsed namelist parameters:\n")
        for tgt, tgt_params in sorted(parsed_targets.items()):
            print(f"{tgt}: {len(tgt_params)} parameters")
            # Print first 10 as sample
            sorted_list = sorted(tgt_params)
            for param in sorted_list[:10]:
                print(f"  - {param}")
            if len(tgt_params) > 10:
                print(f"  ... and {len(tgt_params) - 10} more")
            print()

        # Show params unique to each target
        print("Parameters unique to each target:\n")
        all_param_names = set.union(*parsed_targets.values())
        for tgt, tgt_params in sorted(parsed_targets.items()):
            other = set.union(*[p for t, p in parsed_targets.items() if t != tgt])
            unique = tgt_params - other
            print(f"{tgt} only ({len(unique)}): {sorted(unique)[:15]}...")
            print()

        # Show params in all targets
        common = set.intersection(*parsed_targets.values())
        print(f"Parameters in ALL targets ({len(common)}): {sorted(common)[:20]}...")

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
