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
    "time_stepper": "Time integration scheme",
    "cfl_adap_dt": "Enable adaptive time stepping based on CFL",
    "cfl_const_dt": "Use constant CFL for time stepping",
    "cfl_target": "Target CFL number for adaptive time stepping",
    "cfl_max": "Maximum allowed CFL number",
    "adap_dt": "Enable adaptive time stepping",
    "adap_dt_tol": "Tolerance for adaptive time stepping",
    "adap_dt_max_iters": "Maximum iterations for adaptive time stepping",

    # Model equations
    "model_eqns": "Model equations",
    "num_fluids": "Number of fluid components",
    "num_patches": "Number of initial condition patches",
    "mpp_lim": "Enable mixture pressure positivity limiter",
    "mixture_err": "Enable mixture error checking",
    "alt_soundspeed": "Use alternative sound speed formulation",

    # WENO reconstruction
    "weno_order": "Order of WENO reconstruction",
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
    "recon_type": "Reconstruction type",
    "muscl_order": "Order of MUSCL reconstruction",
    "muscl_lim": "MUSCL limiter type",

    # Riemann solver
    "riemann_solver": "Riemann solver",
    "wave_speeds": "Wave speed estimates",
    "avg_state": "Average state for Riemann solver",
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
    "bubble_model": "Bubble dynamics model",
    "polytropic": "Enable polytropic gas behavior for bubbles",
    "polydisperse": "Enable polydisperse bubble distribution",
    "nb": "Number of bubble bins for polydisperse model",
    "qbmm": "Enable quadrature-based moment method",
    "R0ref": "Reference bubble radius",
    "Ca": "Cavitation number",
    "Web": "Weber number",
    "Re_inv": "Inverse Reynolds number",
    "viscous": "Enable viscous effects",
    "hypoelasticity": "Enable hypoelastic model",
    "hyperelasticity": "Enable hyperelastic model",
    "surface_tension": "Enable surface tension effects",
    "chemistry": "Enable chemical reactions",
    "mhd": "Enable magnetohydrodynamics",
    "hyper_cleaning": "Enable hyperbolic divergence cleaning for MHD",
    "hyper_cleaning_speed": "Wave speed for hyperbolic divergence cleaning",
    "hyper_cleaning_tau": "Damping time constant for hyperbolic divergence cleaning",
    "relativity": "Enable special relativity",

    # Output
    "run_time_info": "Output run-time information",
    "prim_vars_wrt": "Write primitive variables",
    "cons_vars_wrt": "Write conservative variables",
    "probe_wrt": "Write probe data",
    "integral_wrt": "Write integral data",
    "parallel_io": "Enable parallel I/O",
    "file_per_process": "Write separate file per MPI process",
    "format": "Output format",
    "precision": "Output precision",
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

    # Additional simple params
    "thermal": "Thermal model selection",
    "relax_model": "Relaxation model type",
    "igr_order": "Implicit gradient reconstruction order",
    "pref": "Reference pressure",
    "poly_sigma": "Polydisperse distribution standard deviation",
    "rhoref": "Reference density",
    "sigma": "Surface tension coefficient",
    "Bx0": "Background magnetic field in x-direction",
    "relax": "Enable relaxation terms",
    "adv_n": "Enable advection of number density",
    "cont_damage": "Enable continuum damage model",
    "igr": "Enable implicit gradient reconstruction",
    "down_sample": "Enable output downsampling",
    "perturb_flow_fluid": "Fluid index for flow perturbation",
    "perturb_sph_fluid": "Fluid index for spherical perturbation",
    "dist_type": "Distribution type for polydisperse bubbles",
    "mixlayer_perturb_nk": "Number of perturbation modes for mixing layer",
    "elliptic_smoothing_iters": "Number of elliptic smoothing iterations",
    "mixlayer_vel_coef": "Velocity coefficient for mixing layer",
    "mixlayer_domain": "Mixing layer domain size",
    "mixlayer_perturb_k0": "Base wavenumber for mixing layer perturbation",
    "perturb_flow_mag": "Magnitude of flow perturbation",
    "fluid_rho": "Reference fluid density",
    "sigR": "Bubble radius standard deviation",
    "sigV": "Bubble velocity standard deviation",
    "rhoRV": "Bubble radius-velocity correlation",
    "mixlayer_vel_profile": "Enable mixing layer velocity profile",
    "mixlayer_perturb": "Enable mixing layer perturbation",
    "perturb_flow": "Enable flow perturbation",
    "perturb_sph": "Enable spherical perturbation",
    "cfl_dt": "Enable CFL-based time stepping",
    "pre_stress": "Enable pre-stress initialization",
    "elliptic_smoothing": "Enable elliptic smoothing",
    "simplex_perturb": "Enable simplex noise perturbation",
    "n_start_old": "Starting index from previous simulation",
    "palpha_eps": "Volume fraction epsilon for pressure relaxation",
    "ptgalpha_eps": "Volume fraction epsilon for PTG relaxation",
    "pi_fac": "Pi infinity factor",
    "n_start": "Starting time step index",
    "tau_star": "Non-dimensional relaxation time",
    "cont_damage_s": "Continuum damage shape parameter",
    "alpha_bar": "Average volume fraction",
    "alf_factor": "Artificial viscosity factor",
    "ic_eps": "Interface compression epsilon",
    "ic_beta": "Interface compression beta",
    "powell": "Enable Powell source terms for MHD",
    "igr_pres_lim": "Enable IGR pressure limiting",
    "int_comp": "Enable interface compression",
    "nv_uvm_out_of_core": "Enable NVIDIA UVM out-of-core",
    "nv_uvm_pref_gpu": "Prefer GPU for NVIDIA UVM",
    "nv_uvm_igr_temps_on_gpu": "Store IGR temporaries on GPU",
    "num_igr_iters": "Number of IGR iterations",
    "num_igr_warm_start_iters": "Number of IGR warm-start iterations",
    "igr_iter_solver": "IGR iterative solver type",
    "schlieren_alpha": "Schlieren alpha coefficient",
    "t_tol": "Time tolerance",
    "flux_lim": "Flux limiter type",
    "heat_ratio_wrt": "Write heat capacity ratio field",
    "pi_inf_wrt": "Write pi_inf field",
    "pres_inf_wrt": "Write reference pressure field",
    "qm_wrt": "Write Q-criterion field",
    "liutex_wrt": "Write Liutex vortex field",
    "sim_data": "Enable simulation data output",
    "output_partial_domain": "Enable partial domain output",
    "fft_wrt": "Enable FFT output",
    "kappa_wrt": "Write curvature field",
    "lag_header": "Enable Lagrangian output header",
    "chem_wrt_T": "Write temperature field for chemistry",
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
    (r"patch_icpp\((\d+)\)%model_threshold", "Model threshold for patch {0}"),
    (r"patch_icpp\((\d+)\)%epsilon", "Interface thickness for patch {0}"),
    (r"patch_icpp\((\d+)\)%beta", "Shape parameter beta for patch {0}"),
    (r"patch_icpp\((\d+)\)%a\((\d+)\)", "Shape coefficient a({1}) for patch {0}"),
    (r"patch_icpp\((\d+)\)%cf_val", "Color function value for patch {0}"),
    (r"patch_icpp\((\d+)\)%cv", "Specific heat at constant volume for patch {0}"),
    (r"patch_icpp\((\d+)\)%qv", "Heat of formation for patch {0}"),
    (r"patch_icpp\((\d+)\)%qvp", "Heat of formation prime for patch {0}"),
    (r"patch_icpp\((\d+)\)%hcid", "Hard-coded patch ID for patch {0}"),
    (r"patch_icpp\((\d+)\)%model_spc", "Model spacing for patch {0}"),
    (r"patch_icpp\((\d+)\)%non_axis_sym", "Non-axisymmetric parameter for patch {0}"),
    (r"patch_icpp\((\d+)\)%r0", "Initial bubble radius for patch {0}"),
    (r"patch_icpp\((\d+)\)%v0", "Initial bubble velocity for patch {0}"),
    (r"patch_icpp\((\d+)\)%p0", "Initial bubble pressure for patch {0}"),
    (r"patch_icpp\((\d+)\)%m0", "Initial bubble mass for patch {0}"),
    (r"patch_icpp\((\d+)\)%vel", "Velocity magnitude for patch {0}"),
    (r"patch_icpp\((\d+)\)%alpha", "Volume fraction for patch {0}"),
    (r"patch_icpp\((\d+)\)%alpha_rho", "Partial density for patch {0}"),
    (r"patch_icpp\((\d+)\)%radii", "Radii for patch {0}"),
    (r"patch_icpp\((\d+)\)%normal", "Normal direction for patch {0}"),

    # fluid_pp patterns
    (r"fluid_pp\((\d+)\)%gamma", "Specific heat ratio for fluid {0}"),
    (r"fluid_pp\((\d+)\)%pi_inf", "Stiffness pressure for fluid {0}"),
    (r"fluid_pp\((\d+)\)%G", "Shear modulus for fluid {0}"),
    (r"fluid_pp\((\d+)\)%cv", "Specific heat at constant volume for fluid {0}"),
    (r"fluid_pp\((\d+)\)%qv", "Heat of formation for fluid {0}"),
    (r"fluid_pp\((\d+)\)%qvp", "Heat of formation prime for fluid {0}"),
    (r"fluid_pp\((\d+)\)%Re\((\d+)\)", "Reynolds number component {1} for fluid {0}"),
    (r"fluid_pp\((\d+)\)%mul0", "Reference liquid viscosity for fluid {0}"),
    (r"fluid_pp\((\d+)\)%ss", "Surface tension for fluid {0}"),
    (r"fluid_pp\((\d+)\)%pv", "Vapor pressure for fluid {0}"),
    (r"fluid_pp\((\d+)\)%gamma_v", "Specific heat ratio of vapor phase for fluid {0}"),
    (r"fluid_pp\((\d+)\)%M_v", "Molecular weight of vapor phase for fluid {0}"),
    (r"fluid_pp\((\d+)\)%mu_v", "Viscosity of vapor phase for fluid {0}"),
    (r"fluid_pp\((\d+)\)%k_v", "Thermal conductivity of vapor phase for fluid {0}"),
    (r"fluid_pp\((\d+)\)%cp_v", "Specific heat capacity (const. pressure) of vapor for fluid {0}"),
    (r"fluid_pp\((\d+)\)%D_v", "Vapor mass diffusivity for fluid {0}"),

    # patch_ib patterns
    (r"patch_ib\((\d+)\)%geometry", "Geometry type for immersed boundary {0}"),
    (r"patch_ib\((\d+)\)%x_centroid", "X-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%y_centroid", "Y-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%z_centroid", "Z-coordinate of centroid for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_x", "X-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_y", "Y-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%length_z", "Z-dimension length for IB patch {0}"),
    (r"patch_ib\((\d+)\)%radius", "Radius for IB patch {0}"),
    (r"patch_ib\((\d+)\)%theta", "Theta angle for IB patch {0}"),
    (r"patch_ib\((\d+)\)%c", "Shape parameter c for IB patch {0}"),
    (r"patch_ib\((\d+)\)%p", "Shape parameter p for IB patch {0}"),
    (r"patch_ib\((\d+)\)%t", "Shape parameter t for IB patch {0}"),
    (r"patch_ib\((\d+)\)%m", "Shape parameter m for IB patch {0}"),
    (r"patch_ib\((\d+)\)%mass", "Mass for IB patch {0}"),
    (r"patch_ib\((\d+)\)%vel\((\d+)\)", "Velocity component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%angular_vel\((\d+)\)", "Angular velocity component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%angles\((\d+)\)", "Orientation angle {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%slip", "Enable slip condition for IB patch {0}"),
    (r"patch_ib\((\d+)\)%moving_ibm", "Enable moving boundary for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_filepath", "STL model file path for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_spc", "Model spacing for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_threshold", "Model threshold for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_translate\((\d+)\)", "Model translation component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_scale\((\d+)\)", "Model scale component {1} for IB patch {0}"),
    (r"patch_ib\((\d+)\)%model_rotate\((\d+)\)", "Model rotation component {1} for IB patch {0}"),

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
    (r"acoustic\((\d+)\)%dir", "Direction for acoustic source {0}"),
    (r"acoustic\((\d+)\)%npulse", "Number of pulses for acoustic source {0}"),
    (r"acoustic\((\d+)\)%gauss_sigma_dist", "Gaussian spatial width for acoustic source {0}"),
    (r"acoustic\((\d+)\)%gauss_sigma_time", "Gaussian temporal width for acoustic source {0}"),
    (r"acoustic\((\d+)\)%num_elements", "Number of array elements for acoustic source {0}"),
    (r"acoustic\((\d+)\)%element_on", "Active element index for acoustic source {0}"),
    (r"acoustic\((\d+)\)%element_spacing_angle", "Element spacing angle for acoustic source {0}"),
    (r"acoustic\((\d+)\)%element_polygon_ratio", "Element polygon ratio for acoustic source {0}"),
    (r"acoustic\((\d+)\)%foc_length", "Focal length for acoustic source {0}"),
    (r"acoustic\((\d+)\)%aperture", "Aperture for acoustic source {0}"),
    (r"acoustic\((\d+)\)%rotate_angle", "Rotation angle for acoustic source {0}"),
    (r"acoustic\((\d+)\)%bb_num_freq", "Number of broadband frequencies for source {0}"),
    (r"acoustic\((\d+)\)%bb_bandwidth", "Broadband bandwidth for acoustic source {0}"),
    (r"acoustic\((\d+)\)%bb_lowest_freq", "Lowest broadband frequency for source {0}"),

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
    (r"bub_pp%R0ref", "Reference bubble radius"),
    (r"bub_pp%p0ref", "Reference pressure for bubbles"),
    (r"bub_pp%rho0ref", "Reference density for bubbles"),
    (r"bub_pp%T0ref", "Reference temperature for bubbles"),
    (r"bub_pp%ss", "Surface tension between host and gas (bubble)"),
    (r"bub_pp%pv", "Vapor pressure of host fluid"),
    (r"bub_pp%vd", "Vapor diffusion coefficient"),
    (r"bub_pp%mu_l", "Viscosity of host in liquid state"),
    (r"bub_pp%mu_v", "Viscosity of host in vapor state"),
    (r"bub_pp%mu_g", "Viscosity of gas (bubble)"),
    (r"bub_pp%gam_v", "Specific heat ratio of host in vapor state"),
    (r"bub_pp%gam_g", "Specific heat ratio of gas (bubble)"),
    (r"bub_pp%M_v", "Molecular weight of host vapor"),
    (r"bub_pp%M_g", "Molecular weight of gas (bubble)"),
    (r"bub_pp%k_v", "Thermal conductivity of host in vapor state"),
    (r"bub_pp%k_g", "Thermal conductivity of gas (bubble)"),
    (r"bub_pp%cp_v", "Specific heat (const. pressure) of host vapor"),
    (r"bub_pp%cp_g", "Specific heat (const. pressure) of gas (bubble)"),
    (r"bub_pp%R_v", "Gas constant of host in vapor state"),
    (r"bub_pp%R_g", "Gas constant of gas (bubble)"),
    (r"bub_pp%(\w+)", "Bubble parameter: {0}"),

    # Output array patterns
    (r"schlieren_alpha\((\d+)\)", "Schlieren coefficient for fluid {0}"),
    (r"alpha_rho_wrt\((\d+)\)", "Write partial density for fluid {0}"),
    (r"alpha_wrt\((\d+)\)", "Write volume fraction for fluid {0}"),
    (r"alpha_rho_e_wrt\((\d+)\)", "Write partial energy for fluid {0}"),
    (r"kappa_wrt\((\d+)\)", "Write curvature for fluid {0}"),
    (r"mom_wrt\((\d+)\)", "Write momentum component {0}"),
    (r"vel_wrt\((\d+)\)", "Write velocity component {0}"),
    (r"flux_wrt\((\d+)\)", "Write flux component {0}"),
    (r"omega_wrt\((\d+)\)", "Write vorticity component {0}"),
    (r"chem_wrt_Y\((\d+)\)", "Write mass fraction of species {0}"),

    # Lagrangian output patterns - specific fields first
    (r"lag_pos_wrt", "Write Lagrangian bubble position"),
    (r"lag_pos_prev_wrt", "Write Lagrangian bubble previous position"),
    (r"lag_vel_wrt", "Write Lagrangian bubble velocity"),
    (r"lag_rvel_wrt", "Write Lagrangian bubble radial velocity"),
    (r"lag_rad_wrt", "Write Lagrangian bubble radius"),
    (r"lag_r0_wrt", "Write Lagrangian initial bubble radius"),
    (r"lag_rmax_wrt", "Write Lagrangian max bubble radius"),
    (r"lag_rmin_wrt", "Write Lagrangian min bubble radius"),
    (r"lag_pres_wrt", "Write Lagrangian bubble pressure"),
    (r"lag_mv_wrt", "Write Lagrangian vapor mass"),
    (r"lag_mg_wrt", "Write Lagrangian gas mass"),
    (r"lag_db_wrt", "Write Lagrangian bubble diameter"),
    (r"lag_dphidt_wrt", "Write Lagrangian void fraction time derivative"),
    (r"lag_betaT_wrt", "Write Lagrangian thermal beta coefficient"),
    (r"lag_betaC_wrt", "Write Lagrangian concentration beta coefficient"),
    (r"lag_id_wrt", "Write Lagrangian bubble ID"),
    (r"lag_txt_wrt", "Write Lagrangian data to text files"),
    (r"lag_(\w+)_wrt", "Write Lagrangian {0} field"),

    # Body force patterns
    (r"([kgwp])_([xyz])", "Body force parameter {0} in {1}-direction"),
    (r"bf_([xyz])", "Enable body force in {0}-direction"),

    # simplex patterns
    (r"simplex_params%perturb_dens\((\d+)\)", "Enable density perturbation for fluid {0}"),
    (r"simplex_params%perturb_dens_freq\((\d+)\)", "Density perturbation frequency for fluid {0}"),
    (r"simplex_params%perturb_dens_scale\((\d+)\)", "Density perturbation scale for fluid {0}"),
    (r"simplex_params%perturb_dens_offset\((\d+),\s*(\d+)\)", "Density perturbation offset ({1}) for fluid {0}"),
    (r"simplex_params%perturb_vel\((\d+)\)", "Enable velocity perturbation for direction {0}"),
    (r"simplex_params%perturb_vel_freq\((\d+)\)", "Velocity perturbation frequency for direction {0}"),
    (r"simplex_params%perturb_vel_scale\((\d+)\)", "Velocity perturbation scale for direction {0}"),
    (r"simplex_params%perturb_vel_offset\((\d+),(\d+)\)", "Velocity perturbation offset ({1}) for direction {0}"),

    # lag_params patterns - specific fields first
    (r"lag_params%solver_approach", "Lagrangian solver approach (1=one-way, 2=two-way coupling)"),
    (r"lag_params%cluster_type", "Cluster model for pressure at infinity"),
    (r"lag_params%pressure_corrector", "Enable cell pressure correction for Lagrangian bubbles"),
    (r"lag_params%smooth_type", "Smoothing function type (1=Gaussian, 2=Delta 3x3)"),
    (r"lag_params%heatTransfer_model", "Enable heat transfer at bubble-liquid interface"),
    (r"lag_params%massTransfer_model", "Enable mass transfer at bubble-liquid interface"),
    (r"lag_params%write_bubbles", "Write bubble evolution data each time step"),
    (r"lag_params%write_bubbles_stats", "Write max/min radius statistics for bubbles"),
    (r"lag_params%nBubs_glb", "Global number of Lagrangian bubbles"),
    (r"lag_params%epsilonb", "Standard deviation scaling for Gaussian smoothing"),
    (r"lag_params%charwidth", "Domain virtual depth for 2D simulations"),
    (r"lag_params%valmaxvoid", "Maximum permitted void fraction"),
    (r"lag_params%T0", "Initial bubble temperature"),
    (r"lag_params%Thost", "Host fluid temperature"),
    (r"lag_params%c0", "Initial sound speed"),
    (r"lag_params%rho0", "Initial density"),
    (r"lag_params%x0", "Initial bubble position"),
    (r"lag_params%(\w+)", "Lagrangian tracking parameter: {0}"),

    # chem_params patterns - specific fields first
    (r"chem_params%diffusion", "Enable species diffusion for chemistry"),
    (r"chem_params%reactions", "Enable chemical reactions"),
    (r"chem_params%gamma_method", "Gamma calculation method (1=formulation, 2=cp/cv ratio)"),
    (r"chem_params%transport_model", "Transport model selection for chemistry"),
    (r"chem_params%(\w+)", "Chemistry parameter: {0}"),

    # fluid_rho patterns
    (r"fluid_rho\((\d+)\)", "Reference density for fluid {0}"),
]


def get_description(param_name: str) -> str:
    """Get description for a parameter from hand-curated or auto-generated sources.

    Priority: hand-curated DESCRIPTIONS > PATTERNS > auto-generated param.description.
    """
    # 1. Hand-curated descriptions (highest quality)
    if param_name in DESCRIPTIONS:
        return DESCRIPTIONS[param_name]

    # 2. Pattern matching for indexed params (hand-curated templates)
    for pattern, template in PATTERNS:
        match = re.fullmatch(pattern, param_name)
        if match:
            return template.format(*match.groups())

    # 3. Auto-generated description from registry (set by _auto_describe at registration)
    from . import REGISTRY  # pylint: disable=import-outside-toplevel
    param = REGISTRY.all_params.get(param_name)
    if param and param.description:
        return param.description

    # 4. Last resort: naming convention inference
    return _infer_from_naming(param_name)


def _infer_from_naming(param_name: str) -> str:  # pylint: disable=too-many-return-statements,too-many-branches
    """Infer description from naming conventions."""
    name = param_name

    # Handle nested params (e.g., simplex_params%perturb_dens_offset)
    if "%" in name:
        parts = name.split("%")
        prefix = parts[0]
        suffix = parts[1]

        # Extract prefix context
        prefix_map = {
            "simplex_params": "Simplex noise",
            "lag_params": "Lagrangian tracking",
            "chem_params": "Chemistry",
            "bub_pp": "Bubble",
            "x_output": "X-direction output",
            "y_output": "Y-direction output",
            "z_output": "Z-direction output",
        }

        # Remove index from prefix if present
        base_prefix = re.sub(r"\(\d+\)", "", prefix)
        context = prefix_map.get(base_prefix, "")

        # Handle common suffix patterns
        if suffix.endswith("_wrt"):
            field = suffix[:-4].replace("_", " ")
            return f"Write {field} output" + (f" ({context})" if context else "")

        if "offset" in suffix:
            return "Offset parameter" + (f" for {context}" if context else "")

        if context:
            # Clean up suffix for display
            clean_suffix = re.sub(r"\(\d+\)", "", suffix).replace("_", " ")
            return f"{context} {clean_suffix} parameter"

    # Handle *_wrt patterns (write flags)
    if name.endswith("_wrt"):
        field = name[:-4].replace("_", " ")
        return f"Write {field} to output"

    # Handle num_* patterns
    if name.startswith("num_"):
        thing = name[4:].replace("_", " ")
        return f"Number of {thing}"

    # Handle *_order patterns
    if name.endswith("_order"):
        thing = name[:-6].replace("_", " ")
        return f"Order of {thing}"

    # Handle *_model patterns
    if name.endswith("_model"):
        thing = name[:-6].replace("_", " ")
        return f"{thing.title()} model selection"

    # Handle *_tol patterns
    if name.endswith("_tol"):
        thing = name[:-4].replace("_", " ")
        return f"Tolerance for {thing}"

    # Handle *_eps patterns
    if name.endswith("_eps"):
        thing = name[:-4].replace("_", " ")
        return f"Epsilon parameter for {thing}"

    # Handle *_coef or *_coeff patterns
    if name.endswith("_coef") or name.endswith("_coeff"):
        thing = name.rsplit("_", 1)[0].replace("_", " ")
        return f"Coefficient for {thing}"

    # Handle *_max / *_min patterns
    if name.endswith("_max"):
        thing = name[:-4].replace("_", " ")
        return f"Maximum {thing}"
    if name.endswith("_min"):
        thing = name[:-4].replace("_", " ")
        return f"Minimum {thing}"

    # Handle *%beg / *%end patterns
    if name.endswith("%beg"):
        thing = name[:-4].replace("_", " ").replace("%", " ")
        return f"Beginning value for {thing}"
    if name.endswith("%end"):
        thing = name[:-4].replace("_", " ").replace("%", " ")
        return f"End value for {thing}"

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
        # Generic fallback for space-separated indices
        desc = re.sub(r" 1([,\s]|$)", r" N\1", desc)
        desc = re.sub(r" 2([,\s]|$)", r" M\1", desc)
        desc = re.sub(r" 3([,\s]|$)", r" K\1", desc)
        # Handle parenthesized indices (e.g., "a(2)" -> "a(M)")
        desc = re.sub(r"\(1\)", "(N)", desc)
        desc = re.sub(r"\(2\)", "(M)", desc)
        desc = re.sub(r"\(3\)", "(K)", desc)

    return desc


def get_math_symbol(param_name: str) -> str:
    """Get the LaTeX math symbol for a parameter, if one is defined.

    Looks up the math_symbol field from the parameter registry (single source of truth).
    Symbols are defined via math= in the _r() calls in definitions.py.
    """
    from . import REGISTRY  # pylint: disable=import-outside-toplevel
    param = REGISTRY.all_params.get(param_name)
    return param.math_symbol if param else ""


# Feature group descriptions (for display purposes)
# The actual parameter-to-tag mapping is in definitions.py (single source of truth)
FEATURE_DESCRIPTIONS = {
    "mhd": "Magnetohydrodynamics parameters",
    "bubbles": "Bubble dynamics and cavitation",
    "viscosity": "Viscous flow parameters",
    "weno": "WENO reconstruction scheme",
    "time": "Time stepping and integration",
    "output": "Output and visualization",
    "chemistry": "Chemical reactions and species transport",
    "elasticity": "Elastic and hyperelastic materials",
    "acoustic": "Acoustic sources and wave generation",
    "ib": "Immersed boundary method",
    "grid": "Computational grid and domain",
    "bc": "Boundary conditions",
    "riemann": "Riemann solver settings",
    "probes": "Probe points and integral regions",
    "surface_tension": "Surface tension and interface",
    "relativity": "Special relativity",
}
