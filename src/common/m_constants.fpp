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
    integer, parameter :: pathlen_max = 400
    integer, parameter :: nnode = 4    !< Number of QBMM nodes

end module m_constants
