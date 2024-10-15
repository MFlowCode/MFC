!>
!! @file m_constants.f90
!! @brief Contains constant values used throughout the code(s).

module m_constants

    character, parameter :: dflt_char = ' ' !< Default string value

    real(kind(0d0)), parameter :: dflt_real = -1d6                !< Default real value
    real(kind(0d0)), parameter :: sgm_eps = 1d-16               !< Segmentation tolerance
    real(kind(0d0)), parameter :: small_alf = 1d-11                !< Small alf tolerance
    real(kind(0d0)), parameter :: pi = 3.141592653589793d0 !< Pi
    real(kind(0d0)), parameter :: verysmall = 1.d-12              !< Very small number

    integer, parameter :: num_stcls_min = 5    !< Minimum # of stencils
    integer, parameter :: path_len = 400  !< Maximum path length
    integer, parameter :: name_len = 50   !< Maximum name length
    integer, parameter :: dflt_int = -100 !< Default integer value
    integer, parameter :: fourier_rings = 5    !< Fourier filter ring limit
    integer, parameter :: num_fluids_max = 10   !< Maximum number of fluids in the simulation
    integer, parameter :: num_probes_max = 10   !< Maximum number of flow probes in the simulation
    integer, parameter :: num_patches_max = 10
    integer, parameter :: num_bc_patches_max = 20
    integer, parameter :: pathlen_max = 400
    integer, parameter :: nnode = 4    !< Number of QBMM nodes
    integer, parameter :: gp_layers = 3 !< Number of ghost point layers
    integer, parameter :: num_bcs_max = 100 !< Maximum number of boundary condition types.

    real(kind(0d0)), parameter :: capillary_cutoff = 1e-6 !< color function gradient magnitude at which to apply the surface tension fluxes
    real(kind(0d0)), parameter :: acoustic_spatial_support_width = 2.5d0 !< Spatial support width of acoustic source, used in s_source_spatial
    real(kind(0d0)), parameter :: dflt_vcfl_dt = 100d0 !< value of vcfl_dt when viscosity is off for computing adaptive timestep size
    real(kind(0d0)), parameter :: broadband_spectral_level_constant = 20d0 !< The constant to scale the spectral level at the lower frequency bound
    real(kind(0d0)), parameter :: broadband_spectral_level_growth_rate = 10d0 !< The spectral level constant to correct the magnitude at each frqeuency to ensure the source is overall broadband

    ! IBM+STL interpolation constants
    integer, parameter :: Ifactor_2D = 50 !< Multiple factor of the ratio (edge to cell width) for interpolation along edges for 2D models
    integer, parameter :: Ifactor_3D = 5 !< Multiple factor of the ratio (edge to cell width) for interpolation along edges for 3D models
    integer, parameter :: Ifactor_bary_3D = 20 !< Multiple factor of the ratio (triangle area to cell face area) for interpolation on triangle facets for 3D models
    integer, parameter :: num_ray = 20 !< Default number of rays traced per cell
    real(kind(0d0)), parameter :: ray_tracing_threshold = 0.9d0 !< Threshold above which the cell is marked as the model patch
    real(kind(0d0)), parameter :: threshold_vector_zero = 1d-10 !< Threshold to treat the component of a vector to be zero
    real(kind(0d0)), parameter :: threshold_edge_zero = 1d-10 !< Threshold to treat two edges to be overlapped
    real(kind(0d0)), parameter :: threshold_bary = 1d-1 !< Threshold to interpolate a barycentric facet
    real(kind(0d0)), parameter :: initial_distance_buffer = 1d12 !< Initialized levelset distance for the shortest path pair algorithm

    ! Characteristic Boundary Conditions
    integer, parameter :: cbc_range_beg = -13
    integer, parameter :: cbc_range_end =  -5

end module m_constants
