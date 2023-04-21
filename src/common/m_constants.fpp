!>
!! @file m_constants.f90
!! @brief Contains constant values used throughout the code(s).

module m_constants

    use m_precision_select

    character, parameter :: dflt_char = ' ' !< Default string value
    
    real(wp), parameter :: dflt_real = (-1._wp * (10._wp ** 6))                !< Default real value
    real(wp), parameter :: sgm_eps   = (1._wp * (10._wp ** -(16)))               !< Segmentation tolerance
    real(wp), parameter :: small_alf = (1._wp * (10._wp ** -(7)))                !< Small alf tolerance
    real(wp), parameter :: pi        = 3.141592653589793_wp !< Pi
    real(wp), parameter :: verysmall = (1._wp * (10._wp ** -(12)))              !< Very small number
    
    integer, parameter :: num_stcls_min   = 5    !< Mininum # of stencils
    integer, parameter :: path_len        = 400  !< Maximum path length
    integer, parameter :: name_len        = 50   !< Maximum name length
    integer, parameter :: dflt_int        = -100 !< Default integer value
    integer, parameter :: fourier_rings   = 5    !< Fourier filter ring limit
    integer, parameter :: num_fluids_max  = 10   !< Maximum number of fluids in the simulation
    integer, parameter :: num_probes_max  = 10   !< Maximum number of flow probes in the simulation
    integer, parameter :: num_patches_max = 10
    integer, parameter :: pathlen_max     = 400
    integer, parameter :: nnode           = 4    !< Number of QBMM nodes

end module m_constants
