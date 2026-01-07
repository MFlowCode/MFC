"""
MFC Case Parameter Constraint Validator

Validates inter-parameter dependencies and constraints before running
the Fortran codes. This catches configuration errors early and provides
user-friendly error messages.

Based on the constraints enforced in:
- src/common/m_checker_common.fpp
- src/pre_process/m_checker.fpp
- src/simulation/m_checker.fpp
- src/post_process/m_checker.fpp
"""
# pylint: disable=too-many-lines
# Justification: Comprehensive validator covering all MFC parameter constraints

from typing import Dict, Any, List
from .common import MFCException


class CaseConstraintError(MFCException):
    """Exception raised when case parameters violate constraints"""


class CaseValidator:  # pylint: disable=too-many-public-methods
    """Validates MFC case parameter constraints"""

    def __init__(self, params: Dict[str, Any]):
        self.params = params
        self.errors: List[str] = []

    def get(self, key: str, default=None):
        """Get parameter value with default"""
        return self.params.get(key, default)

    def is_set(self, key: str) -> bool:
        """Check if parameter is set (not None and present)"""
        return key in self.params and self.params[key] is not None

    def prohibit(self, condition: bool, message: str):
        """Assert that condition is False, otherwise add error"""
        if condition:
            self.errors.append(message)

    # ===================================================================
    # Common Checks (All Stages)
    # ===================================================================

    def check_simulation_domain(self):
        """Checks constraints on dimensionality and number of cells"""
        m = self.get('m')
        n = self.get('n', 0)
        p = self.get('p', 0)
        cyl_coord = self.get('cyl_coord', 'F') == 'T'

        self.prohibit(m is None, "m must be set")
        self.prohibit(m is not None and m <= 0, "m must be positive")
        self.prohibit(n is not None and n < 0, "n must be non-negative")
        self.prohibit(p is not None and p < 0, "p must be non-negative")
        self.prohibit(cyl_coord and p is not None and p > 0 and p % 2 == 0,
                     "p must be odd for cylindrical coordinates")
        self.prohibit(n is not None and p is not None and n == 0 and p > 0,
                     "p must be 0 if n = 0")

    def check_model_eqns_and_num_fluids(self):
        """Checks constraints on model equations and number of fluids"""
        model_eqns = self.get('model_eqns')
        num_fluids = self.get('num_fluids')
        mpp_lim = self.get('mpp_lim', 'F') == 'T'
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        p = self.get('p', 0)

        self.prohibit(model_eqns is not None and model_eqns not in [1, 2, 3, 4],
                     "model_eqns must be 1, 2, 3, or 4")
        self.prohibit(num_fluids is not None and num_fluids < 1,
                     "num_fluids must be positive")
        self.prohibit(model_eqns == 1 and num_fluids is not None,
                     "num_fluids is not supported for model_eqns = 1")
        self.prohibit(model_eqns == 2 and num_fluids is None,
                     "5-equation model (model_eqns = 2) requires num_fluids to be set")
        self.prohibit(model_eqns == 3 and num_fluids is None,
                     "6-equation model (model_eqns = 3) requires num_fluids to be set")
        self.prohibit(model_eqns == 4 and num_fluids is None,
                     "4-equation model (model_eqns = 4) requires num_fluids to be set")
        self.prohibit(model_eqns == 1 and mpp_lim,
                     "model_eqns = 1 does not support mpp_lim")
        self.prohibit(num_fluids == 1 and mpp_lim,
                     "num_fluids = 1 does not support mpp_lim")
        self.prohibit(model_eqns == 3 and cyl_coord and p != 0,
                     "6-equation model (model_eqns = 3) does not support cylindrical coordinates (cyl_coord = T and p != 0)")

    def check_igr(self):
        """Checks constraints regarding IGR order"""
        igr = self.get('igr', 'F') == 'T'

        if not igr:
            return

        igr_order = self.get('igr_order')
        m = self.get('m', 0)
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(igr_order not in [None, 3, 5],
                     "igr_order must be 3 or 5")
        if igr_order:
            self.prohibit(m + 1 < igr_order,
                         f"m must be at least igr_order - 1 (= {igr_order - 1})")
            self.prohibit(n is not None and n > 0 and n + 1 < igr_order,
                         f"n must be at least igr_order - 1 (= {igr_order - 1})")
            self.prohibit(p is not None and p > 0 and p + 1 < igr_order,
                         f"p must be at least igr_order - 1 (= {igr_order - 1})")

    def check_weno(self):
        """Checks constraints regarding WENO order"""
        recon_type = self.get('recon_type', 1)

        # WENO_TYPE = 1
        if recon_type != 1:
            return

        weno_order = self.get('weno_order')
        m = self.get('m', 0)
        n = self.get('n', 0)
        p = self.get('p', 0)

        if weno_order is None:
            return

        self.prohibit(weno_order not in [1, 3, 5, 7],
                     "weno_order must be 1, 3, 5, or 7")
        self.prohibit(m + 1 < weno_order,
                     f"m must be at least weno_order - 1 (= {weno_order - 1})")
        self.prohibit(n is not None and n > 0 and n + 1 < weno_order,
                     f"For 2D simulation, n must be at least weno_order - 1 (= {weno_order - 1})")
        self.prohibit(p is not None and p > 0 and p + 1 < weno_order,
                     f"For 3D simulation, p must be at least weno_order - 1 (= {weno_order - 1})")

    def check_muscl(self):
        """Check constraints regarding MUSCL order"""
        recon_type = self.get('recon_type', 1)

        # MUSCL_TYPE = 2
        if recon_type != 2:
            return

        muscl_order = self.get('muscl_order')
        m = self.get('m', 0)
        n = self.get('n', 0)
        p = self.get('p', 0)

        if muscl_order is None:
            return

        self.prohibit(muscl_order not in [1, 2],
                     "muscl_order must be 1 or 2")
        self.prohibit(m + 1 < muscl_order,
                     f"m must be at least muscl_order - 1 (= {muscl_order - 1})")
        self.prohibit(n is not None and n > 0 and n + 1 < muscl_order,
                     f"For 2D simulation, n must be at least muscl_order - 1 (= {muscl_order - 1})")
        self.prohibit(p is not None and p > 0 and p + 1 < muscl_order,
                     f"For 3D simulation, p must be at least muscl_order - 1 (= {muscl_order - 1})")

    def check_boundary_conditions(self):  # pylint: disable=too-many-locals
        """Checks constraints on boundary conditions"""
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        m = self.get('m', 0)
        n = self.get('n', 0)
        p = self.get('p', 0)

        for dir, var in [('x', 'm'), ('y', 'n'), ('z', 'p')]:
            var_val = {'m': m, 'n': n, 'p': p}[var]

            for bound in ['beg', 'end']:
                bc_key = f'bc_{dir}%{bound}'
                bc_val = self.get(bc_key)

                self.prohibit(var_val is not None and var_val == 0 and bc_val is not None,
                             f"{bc_key} is not supported for {var} = 0")
                self.prohibit(var_val is not None and var_val > 0 and bc_val is None,
                             f"{var} != 0 but {bc_key} is not set")

            # Check periodicity matches
            beg_bc = self.get(f'bc_{dir}%beg')
            end_bc = self.get(f'bc_{dir}%end')
            if beg_bc is not None and end_bc is not None:
                self.prohibit((beg_bc == -1 and end_bc != -1) or (end_bc == -1 and beg_bc != -1),
                             f"bc_{dir}%beg and bc_{dir}%end must be both periodic (= -1) or both non-periodic")

            # Range check (skip for cylindrical y/z)
            skip_check = cyl_coord and dir in ['y', 'z']
            for bound in ['beg', 'end']:
                bc_key = f'bc_{dir}%{bound}'
                bc_val = self.get(bc_key)

                if not skip_check and bc_val is not None:
                    self.prohibit(bc_val > -1 or bc_val < -17,
                                 f"{bc_key} must be between -1 and -17")
                    self.prohibit(bc_val == -14 and not cyl_coord,
                                 f"{bc_key} must not be -14 (BC_AXIS) for non-cylindrical coordinates")

        # Check BC_NULL is not used
        for dir in ['x', 'y', 'z']:
            for bound in ['beg', 'end']:
                bc_val = self.get(f'bc_{dir}%{bound}')
                self.prohibit(bc_val == -13,
                             "Boundary condition -13 (BC_NULL) is not supported")

        # Cylindrical specific checks
        if cyl_coord:
            self.prohibit(n is not None and n == 0, "n must be positive (2D or 3D) for cylindrical coordinates")
            bc_y_beg = self.get('bc_y%beg')
            bc_y_end = self.get('bc_y%end')
            bc_z_beg = self.get('bc_z%beg')
            bc_z_end = self.get('bc_z%end')

            self.prohibit(p is not None and p == 0 and bc_y_beg != -2,
                         "bc_y%beg must be -2 (BC_REFLECTIVE) for 2D cylindrical coordinates (p = 0)")
            self.prohibit(p is not None and p > 0 and bc_y_beg != -14,
                         "bc_y%beg must be -14 (BC_AXIS) for 3D cylindrical coordinates (p > 0)")

            if bc_y_end is not None:
                self.prohibit(bc_y_end > -1 or bc_y_end < -17,
                             "bc_y%end must be between -1 and -17")
                self.prohibit(bc_y_end == -14,
                             "bc_y%end must not be -14 (BC_AXIS)")

            # 3D cylindrical
            if p is not None and p > 0:
                self.prohibit(bc_z_beg is not None and bc_z_beg not in [-1, -2],
                             "bc_z%beg must be -1 (periodic) or -2 (reflective) for 3D cylindrical coordinates")
                self.prohibit(bc_z_end is not None and bc_z_end not in [-1, -2],
                             "bc_z%end must be -1 (periodic) or -2 (reflective) for 3D cylindrical coordinates")

    def check_bubbles_euler(self):
        """Checks constraints on bubble parameters"""
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'

        if not bubbles_euler:
            return

        nb = self.get('nb')
        polydisperse = self.get('polydisperse', 'F') == 'T'
        thermal = self.get('thermal')
        model_eqns = self.get('model_eqns')
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        rhoref = self.get('rhoref')
        pref = self.get('pref')
        num_fluids = self.get('num_fluids')

        self.prohibit(nb is None or nb < 1,
                     "The Ensemble-Averaged Bubble Model requires nb >= 1")
        self.prohibit(polydisperse and nb == 1,
                     "Polydisperse bubble dynamics requires nb > 1")
        self.prohibit(polydisperse and nb is not None and nb % 2 == 0,
                     "nb must be odd for polydisperse bubbles")
        self.prohibit(thermal is not None and thermal > 3,
                     "thermal must be <= 3")
        self.prohibit(model_eqns == 3,
                     "Bubble models untested with 6-equation model (model_eqns = 3)")
        self.prohibit(model_eqns == 1,
                     "Bubble models untested with pi-gamma model (model_eqns = 1)")
        self.prohibit(model_eqns == 4 and rhoref is None,
                     "rhoref must be set if using bubbles_euler with model_eqns = 4")
        self.prohibit(model_eqns == 4 and pref is None,
                     "pref must be set if using bubbles_euler with model_eqns = 4")
        self.prohibit(model_eqns == 4 and num_fluids != 1,
                     "4-equation model (model_eqns = 4) is single-component and requires num_fluids = 1")
        self.prohibit(cyl_coord,
                     "Bubble models untested in cylindrical coordinates")

    def check_qbmm_and_polydisperse(self):
        """Checks constraints on QBMM and polydisperse bubble parameters"""
        polydisperse = self.get('polydisperse', 'F') == 'T'
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'
        poly_sigma = self.get('poly_sigma')
        qbmm = self.get('qbmm', 'F') == 'T'
        nnode = self.get('nnode')

        self.prohibit(polydisperse and not bubbles_euler,
                     "Polydisperse bubble modeling requires the bubbles_euler flag to be set")
        self.prohibit(polydisperse and poly_sigma is None,
                     "Polydisperse bubble modeling requires poly_sigma to be set")
        self.prohibit(polydisperse and poly_sigma is not None and poly_sigma <= 0,
                     "poly_sigma must be positive")
        self.prohibit(qbmm and not bubbles_euler,
                     "QBMM requires the bubbles_euler flag to be set")
        self.prohibit(qbmm and nnode is not None and nnode != 4,
                     "QBMM requires nnode = 4")

    def check_adv_n(self):
        """Checks constraints on adv_n flag"""
        adv_n = self.get('adv_n', 'F') == 'T'
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'
        num_fluids = self.get('num_fluids')
        qbmm = self.get('qbmm', 'F') == 'T'

        if not adv_n:
            return

        self.prohibit(not bubbles_euler,
                     "adv_n requires bubbles_euler to be enabled")
        self.prohibit(num_fluids != 1,
                     "adv_n requires num_fluids = 1")
        self.prohibit(qbmm,
                     "adv_n is not compatible with qbmm")

    def check_hypoelasticity(self):
        """Checks constraints on hypoelasticity parameters"""
        hypoelasticity = self.get('hypoelasticity', 'F') == 'T'
        model_eqns = self.get('model_eqns')
        riemann_solver = self.get('riemann_solver')

        if not hypoelasticity:
            return

        self.prohibit(model_eqns is not None and model_eqns != 2,
                     "hypoelasticity requires model_eqns = 2")
        self.prohibit(riemann_solver is not None and riemann_solver != 1,
                     "hypoelasticity requires HLL Riemann solver (riemann_solver = 1)")

    def check_phase_change(self):
        """Checks constraints on phase change parameters"""
        relax = self.get('relax', 'F') == 'T'
        relax_model = self.get('relax_model')
        model_eqns = self.get('model_eqns')
        palpha_eps = self.get('palpha_eps')
        ptgalpha_eps = self.get('ptgalpha_eps')

        if not relax:
            return

        self.prohibit(model_eqns is not None and model_eqns != 3,
                     "phase change (relax) requires model_eqns = 3")
        self.prohibit(relax_model is not None and (relax_model < 0 or relax_model > 6),
                     "relax_model must be between 0 and 6")
        self.prohibit(palpha_eps is not None and palpha_eps <= 0,
                     "palpha_eps must be positive")
        self.prohibit(palpha_eps is not None and palpha_eps >= 1,
                     "palpha_eps must be less than 1")
        self.prohibit(ptgalpha_eps is not None and ptgalpha_eps <= 0,
                     "ptgalpha_eps must be positive")
        self.prohibit(ptgalpha_eps is not None and ptgalpha_eps >= 1,
                     "ptgalpha_eps must be less than 1")

    def check_ibm(self):
        """Checks constraints on Immersed Boundaries parameters"""
        ib = self.get('ib', 'F') == 'T'
        n = self.get('n', 0)
        num_ibs = self.get('num_ibs', 0)

        self.prohibit(ib and n <= 0,
                     "Immersed Boundaries do not work in 1D (requires n > 0)")
        self.prohibit(ib and (num_ibs <= 0 or num_ibs > 10),
                     "num_ibs must be between 1 and num_patches_max (10)")
        self.prohibit(not ib and num_ibs > 0,
                     "num_ibs is set, but ib is not enabled")

    def check_stiffened_eos(self):
        """Checks constraints on stiffened equation of state fluids parameters"""
        num_fluids = self.get('num_fluids')
        model_eqns = self.get('model_eqns')
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'

        if num_fluids is None:
            return

        # Allow one extra fluid property slot when using bubbles_euler
        bub_fac = 1 if (bubbles_euler) else 0

        for i in range(1, num_fluids + 1 + bub_fac):
            gamma = self.get(f'fluid_pp({i})%gamma')
            pi_inf = self.get(f'fluid_pp({i})%pi_inf')
            cv = self.get(f'fluid_pp({i})%cv')

            # Positivity checks
            if gamma is not None:
                self.prohibit(gamma <= 0,
                             f"fluid_pp({i})%gamma must be positive")
            if pi_inf is not None:
                self.prohibit(pi_inf < 0,
                             f"fluid_pp({i})%pi_inf must be non-negative")
            if cv is not None:
                self.prohibit(cv < 0,
                             f"fluid_pp({i})%cv must be positive")

            # Model-specific support
            if model_eqns == 1:
                self.prohibit(gamma is not None,
                             f"model_eqns = 1 does not support fluid_pp({i})%gamma")
                self.prohibit(pi_inf is not None,
                             f"model_eqns = 1 does not support fluid_pp({i})%pi_inf")

    def check_surface_tension(self):
        """Checks constraints on surface tension"""
        surface_tension = self.get('surface_tension', 'F') == 'T'
        sigma = self.get('sigma')
        model_eqns = self.get('model_eqns')
        num_fluids = self.get('num_fluids')

        if not surface_tension and sigma is None:
            return

        self.prohibit(surface_tension and sigma is None,
                     "sigma must be set if surface_tension is enabled")
        self.prohibit(surface_tension and sigma is not None and sigma < 0,
                     "sigma must be greater than or equal to zero")
        self.prohibit(sigma is not None and not surface_tension,
                     "sigma is set but surface_tension is not enabled")
        self.prohibit(surface_tension and model_eqns not in [2, 3],
                     "The surface tension model requires model_eqns = 2 or model_eqns = 3")
        self.prohibit(surface_tension and num_fluids != 2,
                     "The surface tension model requires num_fluids = 2")

    def check_mhd(self):
        """Checks constraints on MHD parameters"""
        mhd = self.get('mhd', 'F') == 'T'
        num_fluids = self.get('num_fluids')
        model_eqns = self.get('model_eqns')
        relativity = self.get('relativity', 'F') == 'T'
        Bx0 = self.get('Bx0')
        n = self.get('n', 0)

        self.prohibit(mhd and num_fluids != 1,
                     "MHD is only available for single-component flows (num_fluids = 1)")
        self.prohibit(mhd and model_eqns != 2,
                     "MHD is only available for the 5-equation model (model_eqns = 2)")
        self.prohibit(relativity and not mhd,
                     "relativity requires mhd to be enabled")
        self.prohibit(Bx0 is not None and not mhd,
                     "Bx0 must not be set if MHD is not enabled")
        self.prohibit(mhd and n is not None and n == 0 and Bx0 is None,
                     "Bx0 must be set in 1D MHD simulations")
        self.prohibit(mhd and n is not None and n > 0 and Bx0 is not None,
                     "Bx0 must not be set in 2D/3D MHD simulations")

    # ===================================================================
    # Simulation-Specific Checks
    # ===================================================================

    def check_riemann_solver(self):
        """Checks constraints on Riemann solver (simulation only)"""
        riemann_solver = self.get('riemann_solver')
        model_eqns = self.get('model_eqns')
        wave_speeds = self.get('wave_speeds')
        avg_state = self.get('avg_state')
        low_Mach = self.get('low_Mach', 0)
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        viscous = self.get('viscous', 'F') == 'T'

        if riemann_solver is None:
            return

        self.prohibit(riemann_solver < 1 or riemann_solver > 5,
                     "riemann_solver must be 1, 2, 3, 4 or 5")
        self.prohibit(riemann_solver != 2 and model_eqns == 3,
                     "6-equation model (model_eqns = 3) requires riemann_solver = 2 (HLLC)")
        self.prohibit(wave_speeds is not None and wave_speeds not in [1, 2],
                     "wave_speeds must be 1 or 2")
        self.prohibit(riemann_solver == 3 and wave_speeds is not None,
                     "Exact Riemann (riemann_solver = 3) does not support wave_speeds")
        self.prohibit(avg_state is not None and avg_state not in [1, 2],
                     "avg_state must be 1 or 2")
        self.prohibit(riemann_solver not in [3, 5] and wave_speeds is None,
                     "wave_speeds must be set if riemann_solver != 3,5")
        self.prohibit(riemann_solver not in [3, 5] and avg_state is None,
                     "avg_state must be set if riemann_solver != 3,5")
        self.prohibit(low_Mach not in [0, 1, 2],
                     "low_Mach must be 0, 1, or 2")
        self.prohibit(riemann_solver != 2 and low_Mach == 2,
                     "low_Mach = 2 requires riemann_solver = 2")
        self.prohibit(low_Mach != 0 and model_eqns not in [2, 3],
                     "low_Mach = 1 or 2 requires model_eqns = 2 or 3")
        self.prohibit(riemann_solver == 5 and cyl_coord and viscous,
                     "Lax Friedrichs with cylindrical viscous flux not supported")

    def check_time_stepping(self):
        """Checks time stepping parameters (simulation/post-process)"""
        cfl_dt = self.get('cfl_dt', 'F') == 'T'
        time_stepper = self.get('time_stepper')

        # Check time_stepper bounds
        self.prohibit(time_stepper is not None and (time_stepper < 1 or time_stepper > 3),
                     "time_stepper must be 1, 2, or 3")

        if cfl_dt:
            cfl_target = self.get('cfl_target')
            t_stop = self.get('t_stop')
            t_save = self.get('t_save')
            n_start = self.get('n_start')

            self.prohibit(cfl_target is not None and (cfl_target < 0 or cfl_target > 1),
                         "cfl_target must be between 0 and 1")
            self.prohibit(t_stop is not None and t_stop <= 0,
                         "t_stop must be positive")
            self.prohibit(t_save is not None and t_save <= 0,
                         "t_save must be positive")
            self.prohibit(t_save is not None and t_stop is not None and t_save > t_stop,
                         "t_save must be <= t_stop")
            self.prohibit(n_start is not None and n_start < 0,
                         "n_start must be non-negative")
        else:
            t_step_start = self.get('t_step_start')
            t_step_stop = self.get('t_step_stop')
            t_step_save = self.get('t_step_save')
            dt = self.get('dt')

            self.prohibit(t_step_start is not None and t_step_start < 0,
                         "t_step_start must be non-negative")
            self.prohibit(t_step_stop is not None and t_step_start is not None and t_step_stop <= t_step_start,
                         "t_step_stop must be > t_step_start")
            self.prohibit(t_step_save is not None and t_step_stop is not None and t_step_start is not None and
                         t_step_save > t_step_stop - t_step_start,
                         "t_step_save must be <= (t_step_stop - t_step_start)")
            self.prohibit(dt is not None and dt <= 0,
                         "dt must be positive")

    def check_finite_difference(self):
        """Checks constraints on finite difference parameters"""
        fd_order = self.get('fd_order')

        if fd_order is None:
            return

        self.prohibit(fd_order not in [1, 2, 4],
                     "fd_order must be 1, 2, or 4")

    def check_weno_simulation(self):
        """Checks WENO-specific constraints for simulation"""
        weno_order = self.get('weno_order')
        weno_eps = self.get('weno_eps')
        wenoz = self.get('wenoz', 'F') == 'T'
        wenoz_q = self.get('wenoz_q')
        teno = self.get('teno', 'F') == 'T'
        teno_CT = self.get('teno_CT')
        mapped_weno = self.get('mapped_weno', 'F') == 'T'
        mp_weno = self.get('mp_weno', 'F') == 'T'
        weno_avg = self.get('weno_avg', 'F') == 'T'
        model_eqns = self.get('model_eqns')

        # Check for multiple WENO schemes (regardless of weno_order being set)
        num_schemes = sum([mapped_weno, wenoz, teno])
        self.prohibit(num_schemes >= 2,
                     "Only one of mapped_weno, wenoz, or teno can be set to true")

        # Early return if weno_order not set (other checks need it)
        if weno_order is None:
            return

        self.prohibit(weno_order != 1 and weno_eps is None,
                     "weno_order != 1 requires weno_eps to be set. A typical value is 1e-6")
        self.prohibit(weno_eps is not None and weno_eps <= 0,
                     "weno_eps must be positive. A typical value is 1e-6")
        self.prohibit(wenoz and weno_order == 7 and wenoz_q is None,
                     "wenoz at 7th order requires wenoz_q to be set (should be 2, 3, or 4)")
        self.prohibit(wenoz and weno_order == 7 and wenoz_q is not None and wenoz_q not in [2, 3, 4],
                     "wenoz_q must be either 2, 3, or 4)")
        self.prohibit(teno and teno_CT is None,
                     "teno requires teno_CT to be set. A typical value is 1e-6")
        self.prohibit(teno and teno_CT is not None and teno_CT <= 0,
                     "teno_CT must be positive. A typical value is 1e-6")

        self.prohibit(weno_order == 1 and mapped_weno,
                     "mapped_weno is not compatible with weno_order = 1")
        self.prohibit(weno_order == 1 and wenoz,
                     "wenoz is not compatible with weno_order = 1")
        self.prohibit(weno_order in [1, 3] and teno,
                     "teno requires weno_order = 5 or 7")
        self.prohibit(weno_order != 5 and mp_weno,
                     "mp_weno requires weno_order = 5")
        self.prohibit(model_eqns == 1 and weno_avg,
                     "weno_avg is not compatible with model_eqns = 1")

    def check_muscl_simulation(self):
        """Checks MUSCL-specific constraints for simulation"""
        muscl_order = self.get('muscl_order')
        muscl_lim = self.get('muscl_lim')

        if muscl_order is None:
            return

        self.prohibit(muscl_order == 2 and muscl_lim is None,
                     "muscl_lim must be defined if using muscl_order = 2")
        self.prohibit(muscl_lim is not None and (muscl_lim < 1 or muscl_lim > 5),
                     "muscl_lim must be 1, 2, 3, 4, or 5")

    def check_model_eqns_simulation(self):
        """Checks model equation constraints specific to simulation"""
        model_eqns = self.get('model_eqns')
        avg_state = self.get('avg_state')
        wave_speeds = self.get('wave_speeds')

        if model_eqns != 3:
            return

        self.prohibit(avg_state is not None and avg_state != 2,
                     "6-equation model (model_eqns = 3) requires avg_state = 2")
        self.prohibit(wave_speeds is not None and wave_speeds != 1,
                     "6-equation model (model_eqns = 3) requires wave_speeds = 1")

    def check_bubbles_euler_simulation(self):
        """Checks bubble constraints specific to simulation"""
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'
        bubbles_lagrange = self.get('bubbles_lagrange', 'F') == 'T'
        riemann_solver = self.get('riemann_solver')
        avg_state = self.get('avg_state')
        model_eqns = self.get('model_eqns')
        bubble_model = self.get('bubble_model')

        self.prohibit(bubbles_euler and bubbles_lagrange,
                     "Activate only one of the bubble subgrid models (bubbles_euler or bubbles_lagrange)")

        if not bubbles_euler:
            return

        self.prohibit(riemann_solver is not None and riemann_solver != 2,
                     "Bubble modeling requires HLLC Riemann solver (riemann_solver = 2)")
        self.prohibit(avg_state is not None and avg_state != 2,
                     "Bubble modeling requires arithmetic average (avg_state = 2)")
        self.prohibit(model_eqns == 2 and bubble_model == 1,
                     "The 5-equation bubbly flow model does not support bubble_model = 1 (Gilmore)")

    def check_body_forces(self):
        """Checks constraints on body forces parameters"""
        for dir in ['x', 'y', 'z']:
            bf = self.get(f'bf_{dir}', 'F') == 'T'

            if not bf:
                continue

            self.prohibit(self.get(f'k_{dir}') is None,
                         f"k_{dir} must be specified if bf_{dir} is true")
            self.prohibit(self.get(f'w_{dir}') is None,
                         f"w_{dir} must be specified if bf_{dir} is true")
            self.prohibit(self.get(f'p_{dir}') is None,
                         f"p_{dir} must be specified if bf_{dir} is true")
            self.prohibit(self.get(f'g_{dir}') is None,
                         f"g_{dir} must be specified if bf_{dir} is true")

    def check_viscosity(self):
        """Checks constraints on viscosity parameters"""
        viscous = self.get('viscous', 'F') == 'T'
        num_fluids = self.get('num_fluids')
        model_eqns = self.get('model_eqns')
        weno_order = self.get('weno_order')
        weno_avg = self.get('weno_avg', 'F') == 'T'
        igr = self.get('igr', 'F') == 'T'

        # If num_fluids is not set, check at least fluid 1 (for model_eqns=1)
        if num_fluids is None:
            num_fluids = 1

        for i in range(1, num_fluids + 1):
            Re1 = self.get(f'fluid_pp({i})%Re(1)')
            Re2 = self.get(f'fluid_pp({i})%Re(2)')

            for j, Re_val in [(1, Re1), (2, Re2)]:
                if Re_val is not None:
                    self.prohibit(Re_val <= 0,
                                 f"fluid_pp({i})%Re({j}) must be positive")
                    self.prohibit(model_eqns == 1,
                                 f"model_eqns = 1 does not support fluid_pp({i})%Re({j})")
                    if not igr:
                        self.prohibit(weno_order == 1 and not weno_avg,
                                     f"weno_order = 1 without weno_avg does not support fluid_pp({i})%Re({j})")
                    self.prohibit(not viscous,
                                 f"Re({j}) is specified, but viscous is not set to true")

            # Check Re(1) requirement
            self.prohibit(Re1 is None and viscous,
                         f"viscous is set to true, but fluid_pp({i})%Re(1) is not specified")

    def check_mhd_simulation(self):
        """Checks MHD constraints specific to simulation"""
        mhd = self.get('mhd', 'F') == 'T'
        riemann_solver = self.get('riemann_solver')
        relativity = self.get('relativity', 'F') == 'T'
        powell = self.get('powell', 'F') == 'T'
        n = self.get('n', 0)
        fd_order = self.get('fd_order')

        self.prohibit(mhd and riemann_solver is not None and riemann_solver not in [1, 4],
                     "MHD simulations require riemann_solver = 1 (HLL) or riemann_solver = 4 (HLLD)")
        self.prohibit(riemann_solver == 4 and not mhd,
                     "HLLD (riemann_solver = 4) is only available for MHD simulations")
        self.prohibit(riemann_solver == 4 and relativity,
                     "HLLD is not available for RMHD (relativity)")
        self.prohibit(powell and not mhd,
                     "Powell's method requires mhd to be enabled")
        self.prohibit(powell and n is not None and n == 0,
                     "Powell's method is not supported for 1D simulations")
        self.prohibit(powell and fd_order is None,
                     "fd_order must be set if Powell's method is enabled")


    def check_igr_simulation(self):  # pylint: disable=too-many-locals
        """Checks IGR constraints specific to simulation"""
        igr = self.get('igr', 'F') == 'T'

        if not igr:
            return

        num_igr_iters = self.get('num_igr_iters')
        num_igr_warm_start_iters = self.get('num_igr_warm_start_iters')
        igr_iter_solver = self.get('igr_iter_solver')
        alf_factor = self.get('alf_factor')
        model_eqns = self.get('model_eqns')
        ib = self.get('ib', 'F') == 'T'
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'
        bubbles_lagrange = self.get('bubbles_lagrange', 'F') == 'T'
        alt_soundspeed = self.get('alt_soundspeed', 'F') == 'T'
        surface_tension = self.get('surface_tension', 'F') == 'T'
        hypoelasticity = self.get('hypoelasticity', 'F') == 'T'
        acoustic_source = self.get('acoustic_source', 'F') == 'T'
        relax = self.get('relax', 'F') == 'T'
        mhd = self.get('mhd', 'F') == 'T'
        hyperelasticity = self.get('hyperelasticity', 'F') == 'T'
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        probe_wrt = self.get('probe_wrt', 'F') == 'T'

        self.prohibit(num_igr_iters is not None and num_igr_iters < 0,
                     "num_igr_iters must be greater than or equal to 0")
        self.prohibit(num_igr_warm_start_iters is not None and num_igr_warm_start_iters < 0,
                     "num_igr_warm_start_iters must be greater than or equal to 0")
        self.prohibit(igr_iter_solver is not None and igr_iter_solver not in [1, 2],
                     "igr_iter_solver must be 1 or 2")
        self.prohibit(alf_factor is not None and alf_factor < 0,
                     "alf_factor must be non-negative")
        self.prohibit(model_eqns is not None and model_eqns != 2,
                     "IGR only supports model_eqns = 2")
        self.prohibit(ib,
                     "IGR does not support the immersed boundary method")
        self.prohibit(bubbles_euler,
                     "IGR does not support Euler-Euler bubble models")
        self.prohibit(bubbles_lagrange,
                     "IGR does not support Euler-Lagrange bubble models")
        self.prohibit(alt_soundspeed,
                     "IGR does not support alt_soundspeed = T")
        self.prohibit(surface_tension,
                     "IGR does not support surface tension")
        self.prohibit(hypoelasticity,
                     "IGR does not support hypoelasticity")
        self.prohibit(acoustic_source,
                     "IGR does not support acoustic sources")
        self.prohibit(relax,
                     "IGR does not support phase change")
        self.prohibit(mhd,
                     "IGR does not support magnetohydrodynamics")
        self.prohibit(hyperelasticity,
                     "IGR does not support hyperelasticity")
        self.prohibit(cyl_coord,
                     "IGR does not support cylindrical or axisymmetric coordinates")
        self.prohibit(probe_wrt,
                     "IGR does not support probe writes")

        # Check BCs - IGR does not support characteristic BCs
        # Characteristic BCs are BC_CHAR_SLIP_WALL (-5) through BC_CHAR_SUP_OUTFLOW (-12)
        for dir in ['x', 'y', 'z']:
            for bound in ['beg', 'end']:
                bc = self.get(f'bc_{dir}%{bound}')
                if bc is not None:
                    self.prohibit(-12 <= bc <= -5,
                                 f"Characteristic boundary condition bc_{dir}%{bound} is not compatible with IGR")

    def check_acoustic_source(self):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
        """Checks acoustic source parameters (simulation)"""
        acoustic_source = self.get('acoustic_source', 'F') == 'T'

        if not acoustic_source:
            return

        num_source = self.get('num_source')
        n = self.get('n', 0)
        p = self.get('p', 0)
        cyl_coord = self.get('cyl_coord', 'F') == 'T'

        # Determine dimensionality
        if n is not None and n == 0:
            dim = 1
        elif p is not None and p == 0:
            dim = 2
        else:
            dim = 3

        self.prohibit(num_source is None,
                     "num_source must be specified for acoustic_source")
        self.prohibit(num_source is not None and num_source < 0,
                     "num_source must be non-negative")

        if num_source is None or num_source <= 0:
            return

        # Check each acoustic source
        for j in range(1, num_source + 1):
            jstr = str(j)

            support = self.get(f'acoustic({j})%support')
            loc = [self.get(f'acoustic({j})%loc({i})') for i in range(1, 4)]
            mag = self.get(f'acoustic({j})%mag')
            pulse = self.get(f'acoustic({j})%pulse')
            frequency = self.get(f'acoustic({j})%frequency')
            wavelength = self.get(f'acoustic({j})%wavelength')
            gauss_sigma_time = self.get(f'acoustic({j})%gauss_sigma_time')
            gauss_sigma_dist = self.get(f'acoustic({j})%gauss_sigma_dist')
            bb_num_freq = self.get(f'acoustic({j})%bb_num_freq')
            bb_bandwidth = self.get(f'acoustic({j})%bb_bandwidth')
            bb_lowest_freq = self.get(f'acoustic({j})%bb_lowest_freq')
            npulse = self.get(f'acoustic({j})%npulse')
            dipole = self.get(f'acoustic({j})%dipole', 'F') == 'T'
            dir_val = self.get(f'acoustic({j})%dir')
            delay = self.get(f'acoustic({j})%delay')
            length = self.get(f'acoustic({j})%length')
            height = self.get(f'acoustic({j})%height')
            foc_length = self.get(f'acoustic({j})%foc_length')
            aperture = self.get(f'acoustic({j})%aperture')
            num_elements = self.get(f'acoustic({j})%num_elements')
            element_on = self.get(f'acoustic({j})%element_on')
            element_spacing_angle = self.get(f'acoustic({j})%element_spacing_angle')
            element_polygon_ratio = self.get(f'acoustic({j})%element_polygon_ratio')

            self.prohibit(support is None,
                         f"acoustic({jstr})%support must be specified for acoustic_source")

            # Dimension-specific support checks (only if support was specified)
            if support is not None:
                if dim == 1:
                    self.prohibit(support != 1,
                                 f"Only acoustic({jstr})%support = 1 is allowed for 1D simulations")
                    self.prohibit(support == 1 and loc[0] is None,
                                 f"acoustic({jstr})%loc(1) must be specified for support = 1")

                elif dim == 2:
                    if cyl_coord:
                        self.prohibit(support not in [2, 6, 10],
                                     f"Only acoustic({jstr})%support = 2, 6, or 10 is allowed for 2D axisymmetric")
                    else:
                        self.prohibit(support not in [2, 5, 6, 9, 10],
                                     f"Only acoustic({jstr})%support = 2, 5, 6, 9, or 10 is allowed for 2D")

                    if support in [2, 5, 6, 9, 10]:
                        self.prohibit(loc[0] is None or loc[1] is None,
                                     f"acoustic({jstr})%loc(1:2) must be specified for support = {support}")

                elif dim == 3:
                    self.prohibit(support not in [3, 7, 11],
                                 f"Only acoustic({jstr})%support = 3, 7, or 11 is allowed for 3D")
                    self.prohibit(cyl_coord,
                                 "Acoustic source is not supported in 3D cylindrical simulations")

                    if support == 3:
                        self.prohibit(loc[0] is None or loc[1] is None,
                                     f"acoustic({jstr})%loc(1:2) must be specified for support = 3")
                    elif support in [7, 11]:
                        self.prohibit(loc[0] is None or loc[1] is None or loc[2] is None,
                                     f"acoustic({jstr})%loc(1:3) must be specified for support = {support}")

            # Pulse parameters
            self.prohibit(mag is None,
                         f"acoustic({jstr})%mag must be specified")
            self.prohibit(pulse is None,
                         f"acoustic({jstr})%pulse must be specified")
            self.prohibit(pulse is not None and pulse not in [1, 2, 3, 4],
                         f"Only acoustic({jstr})%pulse = 1, 2, 3, or 4 is allowed")

            # Pulse-specific requirements
            if pulse in [1, 3]:
                freq_set = frequency is not None
                wave_set = wavelength is not None
                self.prohibit(freq_set == wave_set,
                             f"One and only one of acoustic({jstr})%frequency or wavelength must be specified for pulse = {pulse}")

            if pulse == 2:
                time_set = gauss_sigma_time is not None
                dist_set = gauss_sigma_dist is not None
                self.prohibit(time_set == dist_set,
                             f"One and only one of acoustic({jstr})%gauss_sigma_time or gauss_sigma_dist must be specified for pulse = 2")
                self.prohibit(delay is None,
                             f"acoustic({jstr})%delay must be specified for pulse = 2 (Gaussian)")

            if pulse == 4:
                self.prohibit(bb_num_freq is None,
                             f"acoustic({jstr})%bb_num_freq must be specified for pulse = 4")
                self.prohibit(bb_bandwidth is None,
                             f"acoustic({jstr})%bb_bandwidth must be specified for pulse = 4")
                self.prohibit(bb_lowest_freq is None,
                             f"acoustic({jstr})%bb_lowest_freq must be specified for pulse = 4")

            # npulse checks
            self.prohibit(npulse is None,
                         f"acoustic({jstr})%npulse must be specified")
            self.prohibit(support is not None and support >= 5 and npulse is not None and not isinstance(npulse, int),
                         f"acoustic({jstr})%npulse must be an integer for support >= 5 (non-planar)")
            self.prohibit(npulse is not None and npulse >= 5 and dipole,
                         f"acoustic({jstr})%dipole is not supported for npulse >= 5")
            self.prohibit(support is not None and support < 5 and dir_val is None,
                         f"acoustic({jstr})%dir must be specified for support < 5 (planar)")
            self.prohibit(support == 1 and dir_val is not None and dir_val == 0,
                         f"acoustic({jstr})%dir must be non-zero for support = 1")
            self.prohibit(pulse == 3 and support is not None and support >= 5,
                         f"acoustic({jstr})%support >= 5 is not allowed for pulse = 3 (square wave)")

            # Geometry checks
            if support in [2, 3]:
                self.prohibit(length is None,
                             f"acoustic({jstr})%length must be specified for support = {support}")
                self.prohibit(length is not None and length <= 0,
                             f"acoustic({jstr})%length must be positive for support = {support}")

            if support == 3:
                self.prohibit(height is None,
                             f"acoustic({jstr})%height must be specified for support = 3")
                self.prohibit(height is not None and height <= 0,
                             f"acoustic({jstr})%height must be positive for support = 3")

            if support is not None and support >= 5:
                self.prohibit(foc_length is None,
                             f"acoustic({jstr})%foc_length must be specified for support >= 5 (non-planar)")
                self.prohibit(foc_length is not None and foc_length <= 0,
                             f"acoustic({jstr})%foc_length must be positive for support >= 5")
                self.prohibit(aperture is None,
                             f"acoustic({jstr})%aperture must be specified for support >= 5 (non-planar)")
                self.prohibit(aperture is not None and aperture <= 0,
                             f"acoustic({jstr})%aperture must be positive for support >= 5")

            # Transducer array checks
            if support in [9, 10, 11]:
                self.prohibit(num_elements is None,
                             f"acoustic({jstr})%num_elements must be specified for support = {support} (transducer array)")
                self.prohibit(num_elements is not None and num_elements <= 0,
                             f"acoustic({jstr})%num_elements must be positive for support = {support}")
                self.prohibit(element_on is not None and element_on < 0,
                             f"acoustic({jstr})%element_on must be non-negative for support = {support}")
                self.prohibit(element_on is not None and num_elements is not None and element_on > num_elements,
                             f"acoustic({jstr})%element_on must be <= num_elements for support = {support}")

            if support in [9, 10]:
                self.prohibit(element_spacing_angle is None,
                             f"acoustic({jstr})%element_spacing_angle must be specified for support = {support} (2D transducer)")
                self.prohibit(element_spacing_angle is not None and element_spacing_angle < 0,
                             f"acoustic({jstr})%element_spacing_angle must be non-negative for support = {support}")

            if support == 11:
                self.prohibit(element_polygon_ratio is None,
                             f"acoustic({jstr})%element_polygon_ratio must be specified for support = 11 (3D transducer)")
                self.prohibit(element_polygon_ratio is not None and element_polygon_ratio <= 0,
                             f"acoustic({jstr})%element_polygon_ratio must be positive for support = 11")

    def check_adaptive_time_stepping(self):
        """Checks adaptive time stepping parameters (simulation)"""
        adap_dt = self.get('adap_dt', 'F') == 'T'

        if not adap_dt:
            return

        time_stepper = self.get('time_stepper')
        model_eqns = self.get('model_eqns')
        polytropic = self.get('polytropic', 'F') == 'T'
        bubbles_lagrange = self.get('bubbles_lagrange', 'F') == 'T'
        qbmm = self.get('qbmm', 'F') == 'T'
        adv_n = self.get('adv_n', 'F') == 'T'

        self.prohibit(time_stepper is not None and time_stepper != 3,
                     "adap_dt requires Runge-Kutta 3 (time_stepper = 3)")
        self.prohibit(model_eqns == 1,
                     "adap_dt is not supported for model_eqns = 1")
        self.prohibit(qbmm,
                     "adap_dt is not compatible with qbmm")
        self.prohibit(not polytropic and not bubbles_lagrange,
                     "adap_dt requires polytropic = T or bubbles_lagrange = T")
        self.prohibit(not adv_n and not bubbles_lagrange,
                     "adap_dt requires adv_n = T or bubbles_lagrange = T")

    def check_alt_soundspeed(self):
        """Checks alternative sound speed parameters (simulation)"""
        alt_soundspeed = self.get('alt_soundspeed', 'F') == 'T'

        if not alt_soundspeed:
            return

        model_eqns = self.get('model_eqns')
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'
        avg_state = self.get('avg_state')
        riemann_solver = self.get('riemann_solver')
        num_fluids = self.get('num_fluids')

        self.prohibit(model_eqns is not None and model_eqns != 2,
                     "5-equation model (model_eqns = 2) is required for alt_soundspeed")
        self.prohibit(bubbles_euler,
                     "alt_soundspeed is not compatible with bubbles_euler")
        self.prohibit(avg_state is not None and avg_state != 2,
                     "alt_soundspeed requires avg_state = 2")
        self.prohibit(riemann_solver is not None and riemann_solver != 2,
                     "alt_soundspeed requires HLLC Riemann solver (riemann_solver = 2)")
        self.prohibit(num_fluids is not None and num_fluids not in [2, 3],
                     "alt_soundspeed requires num_fluids = 2 or 3")

    def check_bubbles_lagrange(self):
        """Checks Lagrangian bubble parameters (simulation)"""
        bubbles_lagrange = self.get('bubbles_lagrange', 'F') == 'T'

        if not bubbles_lagrange:
            return

        n = self.get('n', 0)
        file_per_process = self.get('file_per_process', 'F') == 'T'
        model_eqns = self.get('model_eqns')
        cluster_type = self.get('lag_params%cluster_type')
        smooth_type = self.get('lag_params%smooth_type')
        polytropic = self.get('polytropic', 'F') == 'T'
        thermal = self.get('thermal')

        self.prohibit(n is not None and n == 0,
                     "bubbles_lagrange accepts 2D and 3D simulations only")
        self.prohibit(file_per_process,
                     "file_per_process must be false for bubbles_lagrange")
        self.prohibit(model_eqns == 3,
                     "The 6-equation flow model does not support bubbles_lagrange")
        self.prohibit(polytropic,
                     "bubbles_lagrange requires polytropic = F")
        self.prohibit(thermal is not None and thermal != 3,
                     "bubbles_lagrange requires thermal = 3")
        self.prohibit(cluster_type is not None and cluster_type >= 2 and smooth_type != 1,
                     "cluster_type >= 2 requires smooth_type = 1")

    def check_continuum_damage(self):
        """Checks continuum damage model parameters (simulation)"""
        cont_damage = self.get('cont_damage', 'F') == 'T'

        if not cont_damage:
            return

        tau_star = self.get('tau_star')
        cont_damage_s = self.get('cont_damage_s')
        alpha_bar = self.get('alpha_bar')
        model_eqns = self.get('model_eqns')

        self.prohibit(tau_star is None,
                     "tau_star must be specified for cont_damage")
        self.prohibit(cont_damage_s is None,
                     "cont_damage_s must be specified for cont_damage")
        self.prohibit(alpha_bar is None,
                     "alpha_bar must be specified for cont_damage")
        self.prohibit(model_eqns is not None and model_eqns != 2,
                     "cont_damage requires model_eqns = 2")

    def check_grcbc(self):
        """Checks Generalized Relaxation Characteristics BC (simulation)"""
        for dir in ['x', 'y', 'z']:
            grcbc_in = self.get(f'bc_{dir}%grcbc_in', 'F') == 'T'
            grcbc_out = self.get(f'bc_{dir}%grcbc_out', 'F') == 'T'
            grcbc_vel_out = self.get(f'bc_{dir}%grcbc_vel_out', 'F') == 'T'
            bc_beg = self.get(f'bc_{dir}%beg')
            bc_end = self.get(f'bc_{dir}%end')

            if grcbc_in:
                # Check if EITHER beg OR end is set to -7
                self.prohibit(bc_beg != -7 and bc_end != -7,
                             f"Subsonic Inflow (grcbc_in) requires bc_{dir}%beg = -7 or bc_{dir}%end = -7")
            if grcbc_out:
                # Check if EITHER beg OR end is set to -8
                self.prohibit(bc_beg != -8 and bc_end != -8,
                             f"Subsonic Outflow (grcbc_out) requires bc_{dir}%beg = -8 or bc_{dir}%end = -8")
            if grcbc_vel_out:
                self.prohibit(bc_beg != -8 and bc_end != -8,
                             f"Subsonic Outflow Velocity (grcbc_vel_out) requires bc_{dir}%beg = -8 or bc_{dir}%end = -8")

    def check_probe_integral_output(self):
        """Checks probe and integral output requirements (simulation)"""
        probe_wrt = self.get('probe_wrt', 'F') == 'T'
        integral_wrt = self.get('integral_wrt', 'F') == 'T'
        fd_order = self.get('fd_order')
        bubbles_euler = self.get('bubbles_euler', 'F') == 'T'

        self.prohibit(probe_wrt and fd_order is None,
                     "fd_order must be specified for probe_wrt")
        self.prohibit(integral_wrt and fd_order is None,
                     "fd_order must be specified for integral_wrt")
        self.prohibit(integral_wrt and not bubbles_euler,
                     "integral_wrt requires bubbles_euler to be enabled")

    def check_hyperelasticity(self):
        """Checks hyperelasticity constraints"""
        hyperelasticity = self.get('hyperelasticity', 'F') == 'T'

        if not hyperelasticity:
            return

        model_eqns = self.get('model_eqns')

        self.prohibit(model_eqns == 1,
                     "hyperelasticity is not supported for model_eqns = 1")
        self.prohibit(model_eqns is not None and model_eqns > 3,
                     "hyperelasticity is not supported for model_eqns > 3")

    # ===================================================================
    # Pre-Process Specific Checks
    # ===================================================================

    def check_restart(self):
        """Checks constraints on restart parameters (pre-process)"""
        old_grid = self.get('old_grid', 'F') == 'T'
        old_ic = self.get('old_ic', 'F') == 'T'
        t_step_old = self.get('t_step_old')
        num_patches = self.get('num_patches', 0)

        self.prohibit(not old_grid and old_ic,
                     "old_ic can only be enabled with old_grid enabled")
        self.prohibit(old_grid and t_step_old is None,
                     "old_grid requires t_step_old to be set")
        self.prohibit(num_patches < 0,
                     "num_patches must be non-negative")
        self.prohibit(num_patches == 0 and t_step_old is None,
                     "num_patches must be positive for the non-restart case")

    def check_qbmm_pre_process(self):
        """Checks QBMM constraints for pre-process"""
        qbmm = self.get('qbmm', 'F') == 'T'
        dist_type = self.get('dist_type')
        rhoRV = self.get('rhoRV')

        if not qbmm:
            return

        self.prohibit(dist_type is None,
                     "dist_type must be set if using QBMM")
        self.prohibit(dist_type is not None and dist_type != 1 and rhoRV is not None and rhoRV > 0,
                     "rhoRV cannot be used with dist_type != 1")

    def check_parallel_io_pre_process(self):
        """Checks parallel I/O constraints (pre-process)"""
        parallel_io = self.get('parallel_io', 'F') == 'T'
        down_sample = self.get('down_sample', 'F') == 'T'
        igr = self.get('igr', 'F') == 'T'
        p = self.get('p', 0)
        file_per_process = self.get('file_per_process', 'F') == 'T'
        m = self.get('m', 0)
        n = self.get('n', 0)

        if down_sample:
            self.prohibit(not parallel_io,
                         "down sample requires parallel_io = T")
            self.prohibit(not igr,
                         "down sample requires igr = T")
            self.prohibit(p == 0,
                         "down sample requires 3D (p > 0)")
            self.prohibit(not file_per_process,
                         "down sample requires file_per_process = T")
            if m is not None and m >= 0:
                self.prohibit((m + 1) % 3 != 0,
                             "down sample requires m divisible by 3")
            if n is not None and n >= 0:
                self.prohibit((n + 1) % 3 != 0,
                             "down sample requires n divisible by 3")
            if p is not None and p >= 0:
                self.prohibit((p + 1) % 3 != 0,
                             "down sample requires p divisible by 3")

    def check_grid_stretching(self):  # pylint: disable=too-many-branches
        """Checks grid stretching constraints (pre-process)"""
        loops_x = self.get('loops_x', 1)
        loops_y = self.get('loops_y', 1)
        stretch_y = self.get('stretch_y', 'F') == 'T'
        stretch_z = self.get('stretch_z', 'F') == 'T'
        old_grid = self.get('old_grid', 'F') == 'T'
        n = self.get('n', 0)
        p = self.get('p', 0)
        cyl_coord = self.get('cyl_coord', 'F') == 'T'

        self.prohibit(loops_x < 1,
                     "loops_x must be at least 1")
        self.prohibit(loops_y < 1,
                     "loops_y must be at least 1")
        self.prohibit(stretch_y and n == 0,
                     "stretch_y requires n > 0")
        self.prohibit(stretch_z and p == 0,
                     "stretch_z requires p > 0")
        self.prohibit(stretch_z and cyl_coord,
                     "stretch_z is not compatible with cylindrical coordinates")

        for direction in ['x', 'y', 'z']:
            stretch = self.get(f'stretch_{direction}', 'F') == 'T'
            if not stretch:
                continue

            a = self.get(f'a_{direction}')
            coord_a = self.get(f'{direction}_a')
            coord_b = self.get(f'{direction}_b')

            self.prohibit(old_grid,
                         f"old_grid and stretch_{direction} are incompatible")
            self.prohibit(a is None,
                         f"a_{direction} must be set with stretch_{direction} enabled")
            self.prohibit(coord_a is None,
                         f"{direction}_a must be set with stretch_{direction} enabled")
            self.prohibit(coord_b is None,
                         f"{direction}_b must be set with stretch_{direction} enabled")
            if coord_a is not None and coord_b is not None:
                self.prohibit(coord_a >= coord_b,
                             f"{direction}_a must be less than {direction}_b with stretch_{direction} enabled")

    def check_perturb_density(self):
        """Checks initial partial density perturbation constraints (pre-process)"""
        perturb_flow = self.get('perturb_flow', 'F') == 'T'
        perturb_flow_fluid = self.get('perturb_flow_fluid')
        perturb_flow_mag = self.get('perturb_flow_mag')
        perturb_sph = self.get('perturb_sph', 'F') == 'T'
        perturb_sph_fluid = self.get('perturb_sph_fluid')
        num_fluids = self.get('num_fluids')

        if perturb_flow:
            self.prohibit(perturb_flow_fluid is None or perturb_flow_mag is None,
                         "perturb_flow_fluid and perturb_flow_mag must be set with perturb_flow = T")
        else:
            self.prohibit(perturb_flow_fluid is not None or perturb_flow_mag is not None,
                         "perturb_flow_fluid and perturb_flow_mag must not be set with perturb_flow = F")

        if num_fluids is not None and perturb_flow_fluid is not None:
            self.prohibit(perturb_flow_fluid > num_fluids or perturb_flow_fluid < 0,
                         "perturb_flow_fluid must be between 0 and num_fluids")

        if perturb_sph:
            self.prohibit(perturb_sph_fluid is None,
                         "perturb_sph_fluid must be set with perturb_sph = T")
        else:
            self.prohibit(perturb_sph_fluid is not None,
                         "perturb_sph_fluid must not be set with perturb_sph = F")

        if num_fluids is not None and perturb_sph_fluid is not None:
            self.prohibit(perturb_sph_fluid > num_fluids or perturb_sph_fluid < 0,
                         "perturb_sph_fluid must be between 0 and num_fluids")

    def check_chemistry(self):
        """Checks chemistry constraints (pre-process)
        
        Note: num_species is set automatically by Cantera at runtime when cantera_file
        is provided. No static validation is performed here - chemistry will fail at
        runtime if misconfigured.
        """

    def check_misc_pre_process(self):
        """Checks miscellaneous pre-process constraints"""
        mixlayer_vel_profile = self.get('mixlayer_vel_profile', 'F') == 'T'
        mixlayer_perturb = self.get('mixlayer_perturb', 'F') == 'T'
        elliptic_smoothing = self.get('elliptic_smoothing', 'F') == 'T'
        elliptic_smoothing_iters = self.get('elliptic_smoothing_iters')
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(mixlayer_vel_profile and n == 0,
                     "mixlayer_vel_profile requires n > 0")
        self.prohibit(mixlayer_perturb and p == 0,
                     "mixlayer_perturb requires p > 0")
        if elliptic_smoothing and elliptic_smoothing_iters is not None:
            self.prohibit(elliptic_smoothing_iters < 1,
                         "elliptic_smoothing_iters must be positive")

    def check_bc_patches(self):  # pylint: disable=too-many-branches,too-many-statements
        """Checks boundary condition patch geometry (pre-process)"""
        num_bc_patches = self.get('num_bc_patches', 0)
        num_bc_patches_max = self.get('num_bc_patches_max', 10)

        if num_bc_patches <= 0:
            return

        self.prohibit(num_bc_patches > num_bc_patches_max,
                     f"num_bc_patches must be <= {num_bc_patches_max}")

        for i in range(1, num_bc_patches + 1):
            geometry = self.get(f'patch_bc({i})%geometry')
            bc_type = self.get(f'patch_bc({i})%type')
            direction = self.get(f'patch_bc({i})%dir')
            radius = self.get(f'patch_bc({i})%radius')
            centroid = [self.get(f'patch_bc({i})%centroid({j})') for j in range(1, 4)]
            length = [self.get(f'patch_bc({i})%length({j})') for j in range(1, 4)]

            if geometry is None:
                continue

            # Line Segment BC (geometry = 1)
            if geometry == 1:
                self.prohibit(radius is not None,
                             f"Line Segment Patch {i} can't have radius defined")
                if direction in [1, 2]:
                    self.prohibit(centroid[direction - 1] is not None or centroid[2] is not None,
                                 f"Line Segment Patch {i} of Dir {direction} can't have centroid in Dir {direction} or 3")
                    self.prohibit(length[direction - 1] is not None or length[2] is not None,
                                 f"Line Segment Patch {i} of Dir {direction} can't have length in Dir {direction} or 3")

            # Circle BC (geometry = 2)
            elif geometry == 2:
                self.prohibit(radius is None,
                             f"Circle Patch {i} must have radius defined")
                self.prohibit(any(length_val is not None for length_val in length),
                             f"Circle Patch {i} can't have lengths defined")
                if direction in [1, 2, 3]:
                    self.prohibit(centroid[direction - 1] is not None,
                                 f"Circle Patch {i} of Dir {direction} can't have centroid in Dir {direction}")

            # Rectangle BC (geometry = 3)
            elif geometry == 3:
                self.prohibit(radius is not None,
                             f"Rectangle Patch {i} can't have radius defined")
                if direction in [1, 2, 3]:
                    self.prohibit(centroid[direction - 1] is not None,
                                 f"Rectangle Patch {i} of Dir {direction} can't have centroid in Dir {direction}")
                    self.prohibit(length[direction - 1] is not None,
                                 f"Rectangle Patch {i} of Dir {direction} can't have length in Dir {direction}")

            # Check for incompatible BC types
            if bc_type is not None:
                # BC types -14 to -4, -1 (periodic), or < -17 (dirichlet) are incompatible with patches
                self.prohibit((-14 <= bc_type <= -4) or bc_type == -1 or bc_type < -17,
                             f"Incompatible BC type for boundary condition patch {i}")

    # ===================================================================
    # Post-Process Specific Checks
    # ===================================================================

    def check_output_format(self):
        """Checks output format parameters (post-process)"""
        format = self.get('format')
        precision = self.get('precision')

        if format is not None:
            self.prohibit(format not in [1, 2],
                         "format must be 1 or 2")

        if precision is not None:
            self.prohibit(precision not in [1, 2],
                         "precision must be 1 or 2")

    def check_vorticity(self):
        """Checks vorticity parameters (post-process)"""
        omega_wrt = [self.get(f'omega_wrt({i})', 'F') == 'T' for i in range(1, 4)]
        n = self.get('n', 0)
        p = self.get('p', 0)
        fd_order = self.get('fd_order')

        self.prohibit(n is not None and n == 0 and any(omega_wrt),
                     "omega_wrt requires n > 0 (at least 2D)")
        self.prohibit(p is not None and p == 0 and (omega_wrt[0] or omega_wrt[1]),
                     "omega_wrt(1) and omega_wrt(2) require p > 0 (3D)")
        self.prohibit(any(omega_wrt) and fd_order is None,
                     "fd_order must be set for omega_wrt")

    def check_schlieren(self):
        """Checks schlieren parameters (post-process)"""
        schlieren_wrt = self.get('schlieren_wrt', 'F') == 'T'
        n = self.get('n', 0)
        fd_order = self.get('fd_order')
        num_fluids = self.get('num_fluids')

        self.prohibit(n is not None and n == 0 and schlieren_wrt,
                     "schlieren_wrt requires n > 0 (at least 2D)")
        self.prohibit(schlieren_wrt and fd_order is None,
                     "fd_order must be set for schlieren_wrt")

        if num_fluids is not None:
            for i in range(1, num_fluids + 1):
                schlieren_alpha = self.get(f'schlieren_alpha({i})')
                if schlieren_alpha is not None:
                    self.prohibit(schlieren_alpha <= 0,
                                 f"schlieren_alpha({i}) must be greater than zero")
                    self.prohibit(not schlieren_wrt,
                                 f"schlieren_alpha({i}) should be set only with schlieren_wrt enabled")

    def check_partial_domain(self):  # pylint: disable=too-many-locals
        """Checks partial domain output constraints (post-process)"""
        output_partial_domain = self.get('output_partial_domain', 'F') == 'T'

        if not output_partial_domain:
            return

        format_val = self.get('format')
        precision = self.get('precision')
        flux_wrt = self.get('flux_wrt', 'F') == 'T'
        heat_ratio_wrt = self.get('heat_ratio_wrt', 'F') == 'T'
        pres_inf_wrt = self.get('pres_inf_wrt', 'F') == 'T'
        c_wrt = self.get('c_wrt', 'F') == 'T'
        schlieren_wrt = self.get('schlieren_wrt', 'F') == 'T'
        qm_wrt = self.get('qm_wrt', 'F') == 'T'
        liutex_wrt = self.get('liutex_wrt', 'F') == 'T'
        ib = self.get('ib', 'F') == 'T'
        omega_wrt = [self.get(f'omega_wrt({i})', 'F') == 'T' for i in range(1, 4)]
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(format_val == 1,
                     "output_partial_domain requires format = 2")
        self.prohibit(precision == 1,
                     "output_partial_domain requires precision = 2")
        self.prohibit(flux_wrt or heat_ratio_wrt or pres_inf_wrt or c_wrt or
                     schlieren_wrt or qm_wrt or liutex_wrt or ib or any(omega_wrt),
                     "output_partial_domain is incompatible with certain output flags")

        x_output_beg = self.get('x_output%beg')
        x_output_end = self.get('x_output%end')
        self.prohibit(x_output_beg is None or x_output_end is None,
                     "x_output%beg and x_output%end must be set for output_partial_domain")

        if n is not None and n != 0:
            y_output_beg = self.get('y_output%beg')
            y_output_end = self.get('y_output%end')
            self.prohibit(y_output_beg is None or y_output_end is None,
                         "y_output%beg and y_output%end must be set for output_partial_domain with n > 0")

        if p is not None and p != 0:
            z_output_beg = self.get('z_output%beg')
            z_output_end = self.get('z_output%end')
            self.prohibit(z_output_beg is None or z_output_end is None,
                         "z_output%beg and z_output%end must be set for output_partial_domain with p > 0")

        for direction in ['x', 'y', 'z']:
            beg = self.get(f'{direction}_output%beg')
            end = self.get(f'{direction}_output%end')
            if beg is not None and end is not None:
                self.prohibit(beg > end,
                             f"{direction}_output%beg must be <= {direction}_output%end")

    def check_partial_density(self):
        """Checks partial density output constraints (post-process)"""
        num_fluids = self.get('num_fluids')
        model_eqns = self.get('model_eqns')

        if num_fluids is None:
            return

        for i in range(1, num_fluids + 1):
            alpha_rho_wrt = self.get(f'alpha_rho_wrt({i})', 'F') == 'T'
            if alpha_rho_wrt:
                self.prohibit(model_eqns == 1,
                             f"alpha_rho_wrt({i}) is not supported for model_eqns = 1")

    def check_momentum_post(self):
        """Checks momentum output constraints (post-process)"""
        mom_wrt = [self.get(f'mom_wrt({i})', 'F') == 'T' for i in range(1, 4)]
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(n == 0 and mom_wrt[1],
                     "mom_wrt(2) requires n > 0")
        self.prohibit(p == 0 and mom_wrt[2],
                     "mom_wrt(3) requires p > 0")

    def check_velocity_post(self):
        """Checks velocity output constraints (post-process)"""
        vel_wrt = [self.get(f'vel_wrt({i})', 'F') == 'T' for i in range(1, 4)]
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(n == 0 and vel_wrt[1],
                     "vel_wrt(2) requires n > 0")
        self.prohibit(p == 0 and vel_wrt[2],
                     "vel_wrt(3) requires p > 0")

    def check_flux_limiter(self):
        """Checks flux limiter constraints (post-process)"""
        flux_wrt = [self.get(f'flux_wrt({i})', 'F') == 'T' for i in range(1, 4)]
        flux_lim = self.get('flux_lim')
        n = self.get('n', 0)
        p = self.get('p', 0)

        self.prohibit(n == 0 and flux_wrt[1],
                     "flux_wrt(2) requires n > 0")
        self.prohibit(p == 0 and flux_wrt[2],
                     "flux_wrt(3) requires p > 0")

        if flux_lim is not None:
            self.prohibit(flux_lim not in [1, 2, 3, 4, 5, 6, 7],
                         "flux_lim must be between 1 and 7")

    def check_volume_fraction(self):
        """Checks volume fraction output constraints (post-process)"""
        num_fluids = self.get('num_fluids')
        model_eqns = self.get('model_eqns')

        if num_fluids is None:
            return

        for i in range(1, num_fluids + 1):
            alpha_wrt = self.get(f'alpha_wrt({i})', 'F') == 'T'
            if alpha_wrt:
                self.prohibit(model_eqns == 1,
                             f"alpha_wrt({i}) is not supported for model_eqns = 1")

    def check_fft(self):
        """Checks FFT output constraints (post-process)"""
        fft_wrt = self.get('fft_wrt', 'F') == 'T'

        if not fft_wrt:
            return

        n = self.get('n', 0)
        p = self.get('p', 0)
        cyl_coord = self.get('cyl_coord', 'F') == 'T'
        m_glb = self.get('m_glb')
        n_glb = self.get('n_glb')
        p_glb = self.get('p_glb')

        self.prohibit(n == 0 or p == 0,
                     "FFT WRT only supported in 3D")
        self.prohibit(cyl_coord,
                     "FFT WRT incompatible with cylindrical coordinates")

        if m_glb is not None and n_glb is not None and p_glb is not None:
            self.prohibit((m_glb + 1) % 2 != 0 or (n_glb + 1) % 2 != 0 or (p_glb + 1) % 2 != 0,
                         "FFT WRT requires global dimensions divisible by 2")

        # BC checks: all boundaries must be periodic (-1)
        for direction in ['x', 'y', 'z']:
            for end in ['beg', 'end']:
                bc_val = self.get(f'bc_{direction}%{end}')
                if bc_val is not None:
                    self.prohibit(bc_val != -1,
                                 "FFT WRT requires periodic BCs (all BCs should be -1)")

    def check_qm(self):
        """Checks Q-criterion output constraints (post-process)"""
        qm_wrt = self.get('qm_wrt', 'F') == 'T'
        n = self.get('n', 0)

        self.prohibit(n == 0 and qm_wrt,
                     "qm_wrt requires n > 0 (at least 2D)")

    def check_liutex_post(self):
        """Checks liutex output constraints (post-process)"""
        liutex_wrt = self.get('liutex_wrt', 'F') == 'T'
        n = self.get('n', 0)

        self.prohibit(n == 0 and liutex_wrt,
                     "liutex_wrt requires n > 0 (at least 2D)")

    def check_surface_tension_post(self):
        """Checks surface tension output constraints (post-process)"""
        cf_wrt = self.get('cf_wrt', 'F') == 'T'
        surface_tension = self.get('surface_tension', 'F') == 'T'

        self.prohibit(cf_wrt and not surface_tension,
                     "cf_wrt can only be enabled if surface_tension is enabled")

    def check_no_flow_variables(self):  # pylint: disable=too-many-locals
        """Checks that at least one flow variable is selected (post-process)"""
        rho_wrt = self.get('rho_wrt', 'F') == 'T'
        E_wrt = self.get('E_wrt', 'F') == 'T'
        pres_wrt = self.get('pres_wrt', 'F') == 'T'
        gamma_wrt = self.get('gamma_wrt', 'F') == 'T'
        heat_ratio_wrt = self.get('heat_ratio_wrt', 'F') == 'T'
        pi_inf_wrt = self.get('pi_inf_wrt', 'F') == 'T'
        pres_inf_wrt = self.get('pres_inf_wrt', 'F') == 'T'
        cons_vars_wrt = self.get('cons_vars_wrt', 'F') == 'T'
        prim_vars_wrt = self.get('prim_vars_wrt', 'F') == 'T'
        c_wrt = self.get('c_wrt', 'F') == 'T'
        schlieren_wrt = self.get('schlieren_wrt', 'F') == 'T'

        # Check array variables
        num_fluids = self.get('num_fluids')
        if num_fluids is None:
            num_fluids = 1
        alpha_rho_wrt_any = any(self.get(f'alpha_rho_wrt({i})', 'F') == 'T' for i in range(1, num_fluids + 1))
        mom_wrt_any = any(self.get(f'mom_wrt({i})', 'F') == 'T' for i in range(1, 4))
        vel_wrt_any = any(self.get(f'vel_wrt({i})', 'F') == 'T' for i in range(1, 4))
        flux_wrt_any = any(self.get(f'flux_wrt({i})', 'F') == 'T' for i in range(1, 4))
        alpha_wrt_any = any(self.get(f'alpha_wrt({i})', 'F') == 'T' for i in range(1, num_fluids + 1))
        omega_wrt_any = any(self.get(f'omega_wrt({i})', 'F') == 'T' for i in range(1, 4))

        has_output = (rho_wrt or E_wrt or pres_wrt or gamma_wrt or heat_ratio_wrt or
                     pi_inf_wrt or pres_inf_wrt or cons_vars_wrt or prim_vars_wrt or
                     c_wrt or schlieren_wrt or alpha_rho_wrt_any or mom_wrt_any or
                     vel_wrt_any or flux_wrt_any or alpha_wrt_any or omega_wrt_any)

        self.prohibit(not has_output,
                     "None of the flow variables have been selected for post-process")

    # ===================================================================
    # Main Validation Entry Points
    # ===================================================================

    def validate_common(self):
        """Validate parameters common to all stages"""
        self.check_simulation_domain()
        self.check_model_eqns_and_num_fluids()
        self.check_igr()
        self.check_weno()
        self.check_muscl()
        self.check_boundary_conditions()
        self.check_bubbles_euler()
        self.check_qbmm_and_polydisperse()
        self.check_adv_n()
        self.check_hypoelasticity()
        self.check_hyperelasticity()
        self.check_phase_change()
        self.check_ibm()
        self.check_stiffened_eos()
        self.check_surface_tension()
        self.check_mhd()

    def validate_simulation(self):
        """Validate simulation-specific parameters"""
        self.validate_common()
        self.check_finite_difference()
        self.check_time_stepping()
        self.check_riemann_solver()
        self.check_weno_simulation()
        self.check_muscl_simulation()
        self.check_model_eqns_simulation()
        self.check_bubbles_euler_simulation()
        self.check_body_forces()
        self.check_viscosity()
        self.check_mhd_simulation()
        self.check_igr_simulation()
        self.check_acoustic_source()
        self.check_adaptive_time_stepping()
        self.check_alt_soundspeed()
        self.check_bubbles_lagrange()
        self.check_continuum_damage()
        self.check_grcbc()
        self.check_probe_integral_output()

    def validate_pre_process(self):
        """Validate pre-process-specific parameters"""
        self.validate_common()
        self.check_restart()
        self.check_qbmm_pre_process()
        self.check_parallel_io_pre_process()
        self.check_grid_stretching()
        self.check_perturb_density()
        self.check_chemistry()
        self.check_misc_pre_process()
        self.check_bc_patches()

    def validate_post_process(self):
        """Validate post-process-specific parameters"""
        self.validate_common()
        self.check_finite_difference()
        self.check_time_stepping()
        self.check_output_format()
        self.check_partial_domain()
        self.check_partial_density()
        self.check_momentum_post()
        self.check_velocity_post()
        self.check_flux_limiter()
        self.check_volume_fraction()
        self.check_vorticity()
        self.check_fft()
        self.check_qm()
        self.check_liutex_post()
        self.check_schlieren()
        self.check_surface_tension_post()
        self.check_no_flow_variables()

    def validate(self, stage: str = 'simulation'):
        """Main validation method

        Args:
            stage: One of 'simulation', 'pre_process', or 'post_process'
                   Other stages (like 'syscheck') have no case constraints and are skipped

        Raises:
            CaseConstraintError: If any constraint violations are found
        """
        self.errors = []

        if stage == 'simulation':
            self.validate_simulation()
        elif stage == 'pre_process':
            self.validate_pre_process()
        elif stage == 'post_process':
            self.validate_post_process()
        else:
            # No stage-specific constraints for auxiliary targets like 'syscheck'.
            # Silently skip validation rather than treating this as an error.
            return

        if self.errors:
            error_msg = "Case parameter constraint violations:\n" + "\n".join(f"  • {err}" for err in self.errors)
            raise CaseConstraintError(error_msg)


def validate_case_constraints(params: Dict[str, Any], stage: str = 'simulation'):
    """Convenience function to validate case parameters

    Args:
        params: Dictionary of case parameters
        stage: One of 'simulation', 'pre_process', or 'post_process'

    Raises:
        CaseConstraintError: If any constraint violations are found
    """
    validator = CaseValidator(params)
    validator.validate(stage)
