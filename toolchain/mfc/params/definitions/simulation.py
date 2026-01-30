"""
Simulation-specific Parameter Definitions.

Defines parameters for patch_ib (immersed boundaries), acoustic sources,
probes, integrals, and other simulation-specific indexed parameters.
"""

from ..schema import ParamDef, ParamType, Stage
from ..registry import REGISTRY


NUM_IBS = 10
NUM_PROBES = 10
NUM_ACOUSTIC = 4
NUM_INTEGRALS = 5
NUM_FLUIDS = 10


def _reg(name: str, ptype: ParamType, stages: set, category: str = None):
    """Helper to register a parameter."""
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        stages=stages,
        category=category,
    ))


def register_patch_ib_params():
    """Register patch_ib (immersed boundary) parameters."""

    for ib_id in range(1, NUM_IBS + 1):
        prefix = f"patch_ib({ib_id})%"

        # Basic attributes - PRE_PROCESS
        for attr, ptype in [("geometry", ParamType.INT),
                            ("radius", ParamType.REAL),
                            ("theta", ParamType.REAL),
                            ("slip", ParamType.LOG),
                            ("c", ParamType.REAL),
                            ("p", ParamType.REAL),
                            ("t", ParamType.REAL),
                            ("m", ParamType.REAL),
                            ("moving_ibm", ParamType.INT),
                            ("mass", ParamType.REAL)]:
            _reg(f"{prefix}{attr}", ptype, {Stage.PRE_PROCESS, Stage.SIMULATION}, "patch_ib")

        # Angles - PRE_PROCESS and SIMULATION
        for dir_id in range(1, 4):
            _reg(f"{prefix}angles({dir_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "patch_ib")

        # Centroid and length - PRE_PROCESS and SIMULATION
        for cmp in ["x", "y", "z"]:
            _reg(f"{prefix}{cmp}_centroid", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "patch_ib")
            _reg(f"{prefix}length_{cmp}", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "patch_ib")

        # Model STL - PRE_PROCESS only
        for attr, ptype in [("model_filepath", ParamType.STR),
                            ("model_spc", ParamType.INT),
                            ("model_threshold", ParamType.REAL)]:
            _reg(f"{prefix}{attr}", ptype, {Stage.PRE_PROCESS}, "patch_ib")

        for transform in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _reg(f"{prefix}model_{transform}({j})", ParamType.REAL,
                     {Stage.PRE_PROCESS}, "patch_ib")

        # SIMULATION-only: vel and angular_vel (analytic)
        for dir_id in range(1, 4):
            _reg(f"{prefix}vel({dir_id})", ParamType.ANALYTIC_REAL,
                 {Stage.SIMULATION}, "patch_ib")
            _reg(f"{prefix}angular_vel({dir_id})", ParamType.ANALYTIC_REAL,
                 {Stage.SIMULATION}, "patch_ib")


def register_acoustic_params():
    """Register acoustic source parameters."""

    for mono_id in range(1, NUM_ACOUSTIC + 1):
        prefix = f"acoustic({mono_id})%"

        # Integer attributes
        for attr in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            _reg(f"{prefix}{attr}", ParamType.INT, {Stage.SIMULATION}, "acoustic")

        # Logical attribute
        _reg(f"{prefix}dipole", ParamType.LOG, {Stage.SIMULATION}, "acoustic")

        # Real attributes
        for attr in ["mag", "length", "height", "wavelength", "frequency",
                     "gauss_sigma_dist", "gauss_sigma_time", "npulse",
                     "dir", "delay", "foc_length", "aperture",
                     "element_spacing_angle", "element_polygon_ratio",
                     "rotate_angle", "bb_bandwidth", "bb_lowest_freq"]:
            _reg(f"{prefix}{attr}", ParamType.REAL, {Stage.SIMULATION}, "acoustic")

        # Location (3D)
        for cmp_id in range(1, 4):
            _reg(f"{prefix}loc({cmp_id})", ParamType.REAL, {Stage.SIMULATION}, "acoustic")


def register_probe_params():
    """Register probe parameters."""

    for probe_id in range(1, NUM_PROBES + 1):
        for cmp in ["x", "y", "z"]:
            _reg(f"probe({probe_id})%{cmp}", ParamType.REAL, {Stage.SIMULATION}, "probe")


def register_integral_params():
    """Register integral region parameters."""

    for int_id in range(1, NUM_INTEGRALS + 1):
        for cmp in ["x", "y", "z"]:
            _reg(f"integral({int_id})%{cmp}min", ParamType.REAL, {Stage.SIMULATION}, "integral")
            _reg(f"integral({int_id})%{cmp}max", ParamType.REAL, {Stage.SIMULATION}, "integral")


def register_bc_extended_params():
    """Register extended boundary condition parameters."""

    for cmp in ["x", "y", "z"]:
        prefix = f"bc_{cmp}%"

        # Velocity BCs - PRE_PROCESS and SIMULATION
        for attr in ["vb1", "vb2", "vb3", "ve1", "ve2", "ve3"]:
            _reg(f"{prefix}{attr}", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")

        # Pressure BCs
        for attr in ["pres_in", "pres_out"]:
            _reg(f"{prefix}{attr}", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")

        # GRCBC (Ghost-cell Riemann Characteristic BC)
        for attr in ["grcbc_in", "grcbc_out", "grcbc_vel_out"]:
            _reg(f"{prefix}{attr}", ParamType.LOG,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")

        # Indexed: alpha_rho_in, alpha_in (1-10)
        for int_id in range(1, NUM_FLUIDS + 1):
            _reg(f"{prefix}alpha_rho_in({int_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")
            _reg(f"{prefix}alpha_in({int_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")

        # vel_in, vel_out (1-3)
        for int_id in range(1, 4):
            _reg(f"{prefix}vel_in({int_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")
            _reg(f"{prefix}vel_out({int_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS, Stage.SIMULATION}, "boundary_conditions")

        # Body force parameters - SIMULATION only
        for var in ["k", "w", "p", "g"]:
            _reg(f"{var}_{cmp}", ParamType.REAL, {Stage.SIMULATION}, "body_force")
        _reg(f"bf_{cmp}", ParamType.LOG, {Stage.SIMULATION}, "body_force")


def register_patch_bc_params():
    """Register patch_bc (boundary condition patch) parameters."""

    for bc_p_id in range(1, NUM_IBS + 1):
        prefix = f"patch_bc({bc_p_id})%"

        # Integer attributes
        for attr in ["geometry", "type", "dir", "loc"]:
            _reg(f"{prefix}{attr}", ParamType.INT, {Stage.PRE_PROCESS}, "patch_bc")

        # Centroid and length (indexed 1-3)
        for d_id in range(1, 4):
            _reg(f"{prefix}centroid({d_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "patch_bc")
            _reg(f"{prefix}length({d_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "patch_bc")

        _reg(f"{prefix}radius", ParamType.REAL, {Stage.PRE_PROCESS}, "patch_bc")


def register_simplex_params():
    """Register simplex perturbation parameters."""

    # perturb_dens for each fluid
    for f_id in range(1, NUM_FLUIDS + 1):
        _reg(f"simplex_params%perturb_dens({f_id})", ParamType.LOG, {Stage.PRE_PROCESS}, "simplex")
        _reg(f"simplex_params%perturb_dens_freq({f_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "simplex")
        _reg(f"simplex_params%perturb_dens_scale({f_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "simplex")

        # perturb_dens_offset(f_id, dir)
        for dir_id in range(1, 4):
            _reg(f"simplex_params%perturb_dens_offset({f_id}, {dir_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS}, "simplex")

    # perturb_vel for each direction
    for d_id in range(1, 4):
        _reg(f"simplex_params%perturb_vel({d_id})", ParamType.LOG, {Stage.PRE_PROCESS}, "simplex")
        _reg(f"simplex_params%perturb_vel_freq({d_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "simplex")
        _reg(f"simplex_params%perturb_vel_scale({d_id})", ParamType.REAL, {Stage.PRE_PROCESS}, "simplex")

        # perturb_vel_offset(d_id, dir)
        for dir_id in range(1, 4):
            _reg(f"simplex_params%perturb_vel_offset({d_id},{dir_id})", ParamType.REAL,
                 {Stage.PRE_PROCESS}, "simplex")


def register_lag_params():
    """Register Lagrangian bubble parameters."""

    # Logical params
    for var in ['heatTransfer_model', 'massTransfer_model', 'pressure_corrector',
                'write_bubbles', 'write_bubbles_stats']:
        _reg(f"lag_params%{var}", ParamType.LOG, {Stage.SIMULATION}, "lag_params")

    # Integer params
    for var in ['solver_approach', 'cluster_type', 'smooth_type', 'nBubs_glb']:
        _reg(f"lag_params%{var}", ParamType.INT, {Stage.SIMULATION}, "lag_params")

    # Real params
    for var in ['epsilonb', 'valmaxvoid', 'charwidth']:
        _reg(f"lag_params%{var}", ParamType.REAL, {Stage.SIMULATION}, "lag_params")


def register_chem_params():
    """Register chemistry parameters."""

    # Logical params
    for var in ['diffusion', 'reactions']:
        _reg(f"chem_params%{var}", ParamType.LOG, {Stage.SIMULATION}, "chem_params")

    # Integer params
    for var in ['gamma_method', 'transport_model']:
        _reg(f"chem_params%{var}", ParamType.INT, {Stage.SIMULATION}, "chem_params")


# Auto-register when module is imported
register_patch_ib_params()
register_acoustic_params()
register_probe_params()
register_integral_params()
register_bc_extended_params()
register_patch_bc_params()
register_simplex_params()
register_lag_params()
register_chem_params()
