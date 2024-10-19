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

    integer, parameter :: num_stcls_min = 5    !< Minimum # of stencils
    integer, parameter :: path_len = 400  !< Maximum path length
    integer, parameter :: name_len = 50   !< Maximum name length
    integer, parameter :: dflt_int = -100 !< Default integer value
    integer, parameter :: fourier_rings = 5    !< Fourier filter ring limit
    integer, parameter :: num_fluids_max = 10   !< Maximum number of fluids in the simulation
    integer, parameter :: num_probes_max = 10   !< Maximum number of flow probes in the simulation
    integer, parameter :: num_patches_max = 10
    integer, parameter :: pathlen_max = 400
    integer, parameter :: nnode = 4    !< Number of QBMM nodes
    real(wp), parameter :: capillary_cutoff = 1e-6_wp !< color function gradient magnitude at which to apply the surface tension fluxes
    real(wp), parameter :: acoustic_spatial_support_width = 2.5_wp !< Spatial support width of acoustic source, used in s_source_spatial
    real(wp), parameter :: dflt_vcfl_dt = 100._wp !< value of vcfl_dt when viscosity is off for computing adaptive timestep size

end module m_constants
