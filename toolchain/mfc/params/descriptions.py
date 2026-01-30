"""
Parameter Descriptions.

Manual descriptions for simple parameters + pattern-based auto-generation for indexed params.
Descriptions extracted from docs/documentation/case.md where available.
"""

import re

# Manual descriptions for simple parameters (from docs/documentation/case.md)
DESCRIPTIONS = {
    # Computational domain
    "m": "Number of grid cells in the x-direction",
    "n": "Number of grid cells in the y-direction",
    "p": "Number of grid cells in the z-direction",
    "x_domain%beg": "Beginning of the x-direction domain",
    "x_domain%end": "End of the x-direction domain",
    "y_domain%beg": "Beginning of the y-direction domain",
    "y_domain%end": "End of the y-direction domain",
    "z_domain%beg": "Beginning of the z-direction domain",
    "z_domain%end": "End of the z-direction domain",
    "cyl_coord": "Enable cylindrical coordinates (2D: axisymmetric, 3D: cylindrical)",
    "stretch_x": "Enable grid stretching in the x-direction",
    "stretch_y": "Enable grid stretching in the y-direction",
    "stretch_z": "Enable grid stretching in the z-direction",
    "a_x": "Rate of grid stretching in the x-direction",
    "a_y": "Rate of grid stretching in the y-direction",
    "a_z": "Rate of grid stretching in the z-direction",
    "x_a": "Start of stretching in negative x-direction",
    "x_b": "Start of stretching in positive x-direction",
    "y_a": "Start of stretching in negative y-direction",
    "y_b": "Start of stretching in positive y-direction",
    "z_a": "Start of stretching in negative z-direction",
    "z_b": "Start of stretching in positive z-direction",
    "loops_x": "Number of times to apply grid stretching in x",
    "loops_y": "Number of times to apply grid stretching in y",
    "loops_z": "Number of times to apply grid stretching in z",

    # Time stepping
    "dt": "Time step size",
    "t_step_start": "Starting time step index",
    "t_step_stop": "Ending time step index",
    "t_step_save": "Time step interval for saving data",
    "t_step_print": "Time step interval for printing info",
    "t_stop": "Simulation stop time",
    "t_save": "Time interval for saving data",
    "time_stepper": "Time integration scheme (1=Euler, 2=TVD-RK2, 3=TVD-RK3, 4=RK4, 5=RK5)",
    "cfl_adap_dt": "Enable adaptive time stepping based on CFL",
    "cfl_const_dt": "Use constant CFL for time stepping",
    "cfl_target": "Target CFL number for adaptive time stepping",
    "cfl_max": "Maximum allowed CFL number",
    "adap_dt": "Enable adaptive time stepping",
    "adap_dt_tol": "Tolerance for adaptive time stepping",
    "adap_dt_max_iters": "Maximum iterations for adaptive time stepping",

    # Model equations
    "model_eqns": "Model equations (1=gamma-law, 2=5-eq, 3=6-eq, 4=4-eq)",
    "num_fluids": "Number of fluid components",
    "num_patches": "Number of initial condition patches",
    "mpp_lim": "Enable mixture pressure positivity limiter",
    "mixture_err": "Enable mixture error checking",
    "alt_soundspeed": "Use alternative sound speed formulation",

    # WENO reconstruction
    "weno_order": "Order of WENO reconstruction (1, 3, 5, or 7)",
    "weno_eps": "WENO epsilon parameter for smoothness",
    "mapped_weno": "Enable mapped WENO scheme",
    "wenoz": "Enable WENO-Z scheme",
    "wenoz_q": "WENO-Z power parameter",
    "teno": "Enable TENO scheme",
    "teno_CT": "TENO cutoff parameter",
    "mp_weno": "Enable monotonicity-preserving WENO",
    "weno_Re_flux": "Enable WENO for viscous fluxes",
    "weno_avg": "Enable WENO averaging",
    "null_weights": "Allow null WENO weights",

    # MUSCL reconstruction
    "recon_type": "Reconstruction type (1=WENO, 2=MUSCL)",
    "muscl_order": "Order of MUSCL reconstruction",
    "muscl_lim": "MUSCL limiter type",

    # Riemann solver
    "riemann_solver": "Riemann solver (1=HLL, 2=HLLC, 3=exact)",
    "wave_speeds": "Wave speed estimates (1=direct, 2=pressure)",
    "avg_state": "Average state for Riemann solver (1=Roe, 2=arithmetic)",
    "low_Mach": "Low Mach number correction",

    # Boundary conditions
    "bc_x%beg": "Boundary condition at x-begin (-1=periodic, -2=reflective, -3=symmetric, etc.)",
    "bc_x%end": "Boundary condition at x-end",
    "bc_y%beg": "Boundary condition at y-begin",
    "bc_y%end": "Boundary condition at y-end",
    "bc_z%beg": "Boundary condition at z-begin",
    "bc_z%end": "Boundary condition at z-end",
    "num_bc_patches": "Number of boundary condition patches",

    # Physics models
    "bubbles_euler": "Enable Euler-Euler bubble model",
    "bubbles_lagrange": "Enable Lagrangian bubble tracking",
    "bubble_model": "Bubble dynamics model (1=Gilmore, 2=Keller-Miksis, 3=Rayleigh-Plesset)",
    "polytropic": "Enable polytropic gas behavior for bubbles",
    "polydisperse": "Enable polydisperse bubble distribution",
    "nb": "Number of bubble bins for polydisperse model",
    "qbmm": "Enable quadrature-based moment method",
    "R0ref": "Reference bubble radius",
    "viscous": "Enable viscous effects",
    "hypoelasticity": "Enable hypoelastic model",
    "hyperelasticity": "Enable hyperelastic model",
    "surface_tension": "Enable surface tension effects",
    "chemistry": "Enable chemical reactions",
    "mhd": "Enable magnetohydrodynamics",
    "relativity": "Enable special relativity",

    # Output
    "run_time_info": "Output run-time information",
    "prim_vars_wrt": "Write primitive variables",
    "cons_vars_wrt": "Write conservative variables",
    "probe_wrt": "Write probe data",
    "integral_wrt": "Write integral data",
    "parallel_io": "Enable parallel I/O",
    "file_per_process": "Write separate file per MPI process",
    "format": "Output format (1=Silo, 2=binary)",
    "precision": "Output precision (1=single, 2=double)",
    "schlieren_wrt": "Write schlieren images",
    "rho_wrt": "Write density field",
    "pres_wrt": "Write pressure field",
    "vel_wrt": "Write velocity field",
    "E_wrt": "Write energy field",
    "gamma_wrt": "Write gamma field",
    "alpha_wrt": "Write volume fraction field",
    "alpha_rho_wrt": "Write partial density field",
    "c_wrt": "Write sound speed field",
    "omega_wrt": "Write vorticity field",
    "cf_wrt": "Write color function field",

    # Immersed boundaries
    "ib": "Enable immersed boundary method",
    "num_ibs": "Number of immersed boundary patches",

    # Acoustic sources
    "acoustic_source": "Enable acoustic source terms",
    "num_source": "Number of acoustic sources",

    # Probes and integrals
    "num_probes": "Number of probe points",
    "num_integrals": "Number of integral regions",

    # MPI/GPU
    "rdma_mpi": "Enable RDMA for MPI communication (GPUs)",

    # Misc
    "case_dir": "Case directory path",
    "cantera_file": "Cantera mechanism file for chemistry",
    "old_grid": "Use grid from previous simulation",
    "old_ic": "Use initial conditions from previous simulation",
    "t_step_old": "Time step to restart from",
    "fd_order": "Finite difference order for gradients",
}

# Patterns for auto-generating descriptions of indexed parameters
PATTERNS = [
    # patch_icpp patterns
    (r"patch_icpp\((\d+)\)%geometry", "Geometry type for initial condition patch {0}"),
    (r"patch_icpp\((\d+)\)%x_centroid", "X-coordinate of centroid for patch {0}"),
    (r"patch_icpp\((\d+)\)%y_centroid", "Y-coordinate of centroid for patch {0}"),
    (r"patch_icpp\((\d+)\)%z_centroid", "Z-coordinate of centroid for patch {0}"),
    (r"patch_icpp\((\d+)\)%length_x", "X-dimension length for patch {0}"),
    (r"patch_icpp\((\d+)\)%length_y", "Y-dimension length for patch {0}"),
    (r"patch_icpp\((\d+)\)%length_z", "Z-dimension length for patch {0}"),
    (r"patch_icpp\((\d+)\)%radius", "Radius for patch {0}"),
    (r"patch_icpp\((\d+)\)%radii\((\d+)\)", "Radius component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%normal\((\d+)\)", "Normal component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%vel\((\d+)\)", "Velocity component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%alpha\((\d+)\)", "Volume fraction of fluid {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%alpha_rho\((\d+)\)", "Partial density of fluid {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%pres", "Pressure for patch {0}"),
    (r"patch_icpp\((\d+)\)%rho", "Density for patch {0}"),
    (r"patch_icpp\((\d+)\)%gamma", "Specific heat ratio for patch {0}"),
    (r"patch_icpp\((\d+)\)%pi_inf", "Stiffness pressure for patch {0}"),
    (r"patch_icpp\((\d+)\)%smoothen", "Enable smoothing for patch {0}"),
    (r"patch_icpp\((\d+)\)%smooth_patch_id", "Patch ID to smooth against for patch {0}"),
    (r"patch_icpp\((\d+)\)%smooth_coeff", "Smoothing coefficient for patch {0}"),
    (r"patch_icpp\((\d+)\)%alter_patch\((\d+)\)", "Alter patch {1} with patch {0}"),
    (r"patch_icpp\((\d+)\)%alter_patch", "Enable patch alteration for patch {0}"),
    (r"patch_icpp\((\d+)\)%Y\((\d+)\)", "Mass fraction of species {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%tau_e\((\d+)\)", "Elastic stress component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%Bx", "X-component of magnetic field for patch {0}"),
    (r"patch_icpp\((\d+)\)%By", "Y-component of magnetic field for patch {0}"),
    (r"patch_icpp\((\d+)\)%Bz", "Z-component of magnetic field for patch {0}"),
    (r"patch_icpp\((\d+)\)%model_filepath", "STL model file path for patch {0}"),
    (r"patch_icpp\((\d+)\)%model_translate\((\d+)\)", "Model translation component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%model_scale\((\d+)\)", "Model scale component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%model_rotate\((\d+)\)", "Model rotation component {1} for patch {0}"),
    (r"patch_icpp\((\d+)\)%epsilon", "Interface thickness for patch {0}"),
    (r"patch_icpp\((\d+)\)%beta", "Shape parameter beta for patch {0}"),

    # fluid_pp patterns
    (r"fluid_pp\((\d+)\)%gamma", "Specific heat ratio for fluid {0}"),
    (r"fluid_pp\((\d+)\)%pi_inf", "Stiffness pressure for fluid {0}"),
    (r"fluid_pp\((\d+)\)%G", "Shear modulus for fluid {0}"),
    (r"fluid_pp\((\d+)\)%cv", "Specific heat at constant volume for fluid {0}"),
    (r"fluid_pp\((\d+)\)%qv", "Heat of formation for fluid {0}"),
    (r"fluid_pp\((\d+)\)%qvp", "Heat of formation prime for fluid {0}"),
    (r"fluid_pp\((\d+)\)%Re\((\d+)\)", "Reynolds number component {1} for fluid {0}"),

    # patch_ib patterns
    (r"patch_ib\((\d+)\)%geometry", "Geometry type for immersed boundary {0}"),
    (r"patch_ib\((\d+)\)%x_centroid", "X-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%y_centroid", "Y-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%z_centroid", "Z-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_x", "X-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_y", "Y-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_z", "Z-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%radius", "Radius for IB patch {0}"),
    (r"patch_ib\((\d+)\)%vel\((\d+)\)", "Velocity component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%angular_vel\((\d+)\)", "Angular velocity component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%angles\((\d+)\)", "Orientation angle {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%slip", "Enable slip condition for IB patch {0}"),
    (r"patch_ib\((\d+)\)%moving_ibm", "Enable moving boundary for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_filepath", "STL model file path for IB patch {0}"),

    # bc patterns
    (r"bc_([xyz])%vel_in\((\d+)\)", "Inlet velocity component {1} at {0}-boundary"),
    (r"bc_([xyz])%vel_out\((\d+)\)", "Outlet velocity component {1} at {0}-boundary"),
    (r"bc_([xyz])%alpha_rho_in\((\d+)\)", "Inlet partial density of fluid {1} at {0}-boundary"),
    (r"bc_([xyz])%alpha_in\((\d+)\)", "Inlet volume fraction of fluid {1} at {0}-boundary"),
    (r"bc_([xyz])%pres_in", "Inlet pressure at {0}-boundary"),
    (r"bc_([xyz])%pres_out", "Outlet pressure at {0}-boundary"),
    (r"bc_([xyz])%vb(\d+)", "Boundary velocity component {1} at {0}-begin"),
    (r"bc_([xyz])%ve(\d+)", "Boundary velocity component {1} at {0}-end"),
    (r"bc_([xyz])%grcbc_in", "Enable GRCBC at {0}-inlet"),
    (r"bc_([xyz])%grcbc_out", "Enable GRCBC at {0}-outlet"),
    (r"bc_([xyz])%grcbc_vel_out", "Enable GRCBC velocity at {0}-outlet"),

    # patch_bc patterns
    (r"patch_bc\((\d+)\)%geometry", "Geometry type for BC patch {0}"),
    (r"patch_bc\((\d+)\)%type", "BC type for patch {0}"),
    (r"patch_bc\((\d+)\)%dir", "Direction for BC patch {0}"),
    (r"patch_bc\((\d+)\)%loc", "Location for BC patch {0}"),
    (r"patch_bc\((\d+)\)%centroid\((\d+)\)", "Centroid component {1} for BC patch {0}"),
    (r"patch_bc\((\d+)\)%length\((\d+)\)", "Length component {1} for BC patch {0}"),
    (r"patch_bc\((\d+)\)%radius", "Radius for BC patch {0}"),

    # acoustic patterns
    (r"acoustic\((\d+)\)%loc\((\d+)\)", "Location component {1} for acoustic source {0}"),
    (r"acoustic\((\d+)\)%mag", "Magnitude for acoustic source {0}"),
    (r"acoustic\((\d+)\)%pulse", "Pulse type for acoustic source {0}"),
    (r"acoustic\((\d+)\)%support", "Support type for acoustic source {0}"),
    (r"acoustic\((\d+)\)%frequency", "Frequency for acoustic source {0}"),
    (r"acoustic\((\d+)\)%wavelength", "Wavelength for acoustic source {0}"),
    (r"acoustic\((\d+)\)%length", "Length for acoustic source {0}"),
    (r"acoustic\((\d+)\)%height", "Height for acoustic source {0}"),
    (r"acoustic\((\d+)\)%delay", "Delay for acoustic source {0}"),
    (r"acoustic\((\d+)\)%dipole", "Enable dipole for acoustic source {0}"),

    # probe patterns
    (r"probe\((\d+)\)%x", "X-coordinate of probe {0}"),
    (r"probe\((\d+)\)%y", "Y-coordinate of probe {0}"),
    (r"probe\((\d+)\)%z", "Z-coordinate of probe {0}"),

    # integral patterns
    (r"integral\((\d+)\)%xmin", "X-min of integral region {0}"),
    (r"integral\((\d+)\)%xmax", "X-max of integral region {0}"),
    (r"integral\((\d+)\)%ymin", "Y-min of integral region {0}"),
    (r"integral\((\d+)\)%ymax", "Y-max of integral region {0}"),
    (r"integral\((\d+)\)%zmin", "Z-min of integral region {0}"),
    (r"integral\((\d+)\)%zmax", "Z-max of integral region {0}"),

    # bub_pp patterns
    (r"bub_pp%(\w+)", "Bubble parameter: {0}"),

    # Output array patterns
    (r"schlieren_alpha\((\d+)\)", "Schlieren coefficient for fluid {0}"),
    (r"alpha_rho_wrt\((\d+)\)", "Write partial density for fluid {0}"),
    (r"alpha_wrt\((\d+)\)", "Write volume fraction for fluid {0}"),
    (r"mom_wrt\((\d+)\)", "Write momentum component {0}"),
    (r"vel_wrt\((\d+)\)", "Write velocity component {0}"),
    (r"flux_wrt\((\d+)\)", "Write flux component {0}"),
    (r"omega_wrt\((\d+)\)", "Write vorticity component {0}"),
    (r"chem_wrt_Y\((\d+)\)", "Write mass fraction of species {0}"),

    # Body force patterns
    (r"([kgwp])_([xyz])", "Body force parameter {0} in {1}-direction"),
    (r"bf_([xyz])", "Enable body force in {0}-direction"),

    # simplex patterns
    (r"simplex_params%perturb_dens\((\d+)\)", "Enable density perturbation for fluid {0}"),
    (r"simplex_params%perturb_dens_freq\((\d+)\)", "Density perturbation frequency for fluid {0}"),
    (r"simplex_params%perturb_dens_scale\((\d+)\)", "Density perturbation scale for fluid {0}"),
    (r"simplex_params%perturb_vel\((\d+)\)", "Enable velocity perturbation for direction {0}"),
    (r"simplex_params%perturb_vel_freq\((\d+)\)", "Velocity perturbation frequency for direction {0}"),
    (r"simplex_params%perturb_vel_scale\((\d+)\)", "Velocity perturbation scale for direction {0}"),

    # lag_params patterns
    (r"lag_params%(\w+)", "Lagrangian tracking parameter: {0}"),

    # chem_params patterns
    (r"chem_params%(\w+)", "Chemistry parameter: {0}"),

    # fluid_rho patterns
    (r"fluid_rho\((\d+)\)", "Reference density for fluid {0}"),
]


def get_description(param_name: str) -> str:
    """Get description for a parameter, using manual or auto-generated."""
    # Check manual descriptions first
    if param_name in DESCRIPTIONS:
        return DESCRIPTIONS[param_name]

    # Try pattern matching for indexed params
    for pattern, template in PATTERNS:
        match = re.match(pattern, param_name)
        if match:
            return template.format(*match.groups())

    return ""


def get_pattern_description(pattern_name: str) -> str:
    """Get description for a collapsed pattern like patch_icpp(N)%geometry."""
    # Convert pattern back to example: patch_icpp(N)%geometry -> patch_icpp(1)%geometry
    # Use different placeholder values so we can distinguish them later
    example = pattern_name.replace("(N)", "(1)").replace("(M)", "(2)").replace("(K)", "(3)")
    desc = get_description(example)

    if desc:
        # Replace specific index values with generic labels
        # First, handle the secondary index (2 -> M)
        desc = re.sub(r"(species|fluid|component|direction) 2", r"\1 M", desc)
        desc = re.sub(r"component 2", "component M", desc)
        # Then handle primary index (1 -> N)
        desc = re.sub(r"(patch|fluid|IB patch|source|probe|region|species|direction|component) 1", r"\1 N", desc)
        # Generic fallback
        desc = re.sub(r" 1([,\s]|$)", r" N\1", desc)
        desc = re.sub(r" 2([,\s]|$)", r" M\1", desc)
        desc = re.sub(r" 3([,\s]|$)", r" K\1", desc)

    return desc
