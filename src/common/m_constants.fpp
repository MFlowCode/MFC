!>
!! @file m_constants.f90
!! @brief Contains constant values used throughout the code(s).

module m_constants

    use m_precision_select

    character, parameter :: dflt_char = ' ' !< Default string value

    real(wp), parameter :: dflt_real = -1e6_wp                !< Default real value
    real(wp), parameter :: sgm_eps = 1e-16_wp               !< Segmentation tolerance
    real(wp), parameter :: small_alf = 1e-11_wp                !< Small alf tolerance
    real(wp), parameter :: pi = 3.141592653589793_wp !< Pi
    real(wp), parameter :: verysmall = 1.e-12_wp              !< Very small number

    integer, parameter :: num_stcls_min = 5                       !< Minimum # of stencils
    integer, parameter :: path_len = 400                          !< Maximum path length
    integer, parameter :: name_len = 50                           !< Maximum name length
    integer, parameter :: dflt_int = -100                         !< Default integer value
    integer, parameter :: fourier_rings = 5                       !< Fourier filter ring limit
    integer, parameter :: num_fluids_max = 10                     !< Maximum number of fluids in the simulation
    integer, parameter :: num_probes_max = 10                     !< Maximum number of flow probes in the simulation
    integer, parameter :: num_patches_max = 10
    integer, parameter :: pathlen_max = 400
    integer, parameter :: nnode = 4    !< Number of QBMM nodes
    integer, parameter :: gp_layers = 3 !< Number of ghost point layers for IBM
    real(wp), parameter :: capillary_cutoff = 1e-6 !< color function gradient magnitude at which to apply the surface tension fluxes
    real(wp), parameter :: acoustic_spatial_support_width = 2.5_wp !< Spatial support width of acoustic source, used in s_source_spatial
    real(wp), parameter :: dflt_vcfl_dt = 100._wp !< value of vcfl_dt when viscosity is off for computing adaptive timestep size
    real(wp), parameter :: broadband_spectral_level_constant = 20._wp !< The constant to scale the spectral level at the lower frequency bound
    real(wp), parameter :: broadband_spectral_level_growth_rate = 10._wp !< The spectral level constant to correct the magnitude at each frqeuency to ensure the source is overall broadband

    ! Chemistry
    real(wp), parameter :: dflt_T_guess = 1200._wp ! Default guess for temperature (when a previous value is not available)

    ! IBM+STL interpolation constants
    integer, parameter :: Ifactor_2D = 50 !< Multiple factor of the ratio (edge to cell width) for interpolation along edges for 2D models
    integer, parameter :: Ifactor_3D = 5 !< Multiple factor of the ratio (edge to cell width) for interpolation along edges for 3D models
    integer, parameter :: Ifactor_bary_3D = 20 !< Multiple factor of the ratio (triangle area to cell face area) for interpolation on triangle facets for 3D models
    integer, parameter :: num_ray = 20 !< Default number of rays traced per cell
    real(wp), parameter :: ray_tracing_threshold = 0.9_wp !< Threshold above which the cell is marked as the model patch
    real(wp), parameter :: threshold_vector_zero = 1e-10 !< Threshold to treat the component of a vector to be zero
    real(wp), parameter :: threshold_edge_zero = 1e-10 !< Threshold to treat two edges to be overlapped
    real(wp), parameter :: threshold_bary = 1e-1 !< Threshold to interpolate a barycentric facet
    real(wp), parameter :: initial_distance_buffer = 1e12_wp !< Initialized levelset distance for the shortest path pair algorithm

    ! Lagrange bubbles constants
    integer, parameter :: mapCells = 3 !< Number of cells around the bubble where the smoothening function will have effect
    real(wp), parameter :: R_uni = 8314._wp ! Universal gas constant - J/kmol/K

    ! Strang Splitting constants
    real(wp), parameter :: dflt_adap_dt_tol = 1e-4_wp !< Default tolerance for adaptive step size
    integer, parameter :: adap_dt_max_iters = 100 !< Maximum number of iterations

    ! Relativity
    integer, parameter :: relativity_cons_to_prim_max_iter = 100

    ! System constants
    integer, parameter :: CASE_FILE_ERROR_CODE = 22

end module m_constants
