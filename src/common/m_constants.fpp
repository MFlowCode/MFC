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
    real(kind(0d0)), parameter :: capillary_cutoff = 1e-6 !< color function gradient magnitude at which to apply the surface tension fluxes
    real(kind(0d0)), parameter :: acoustic_spatial_support_width = 2.5d0 !< Spatial support width of acoustic source, used in s_source_spatial

    integer, parameter :: kLineSegmentPatch = 1
    integer, parameter :: circularPatch = 2
    integer, parameter :: rectangularPatch = 3
    integer, parameter :: sweptLinePatch = 4
    integer, parameter :: ellipticalPatch = 5
    integer, parameter :: unimplementedPatch = 6
    integer, parameter :: funcPatch = 7
    integer, parameter :: spherePatch = 8
    integer, parameter :: cuboidPatch = 9
    integer, parameter :: cylinderPatch = 10
    integer, parameter :: sweptPlanePatch = 11
    integer, parameter :: ellipsoidPatch = 12
    integer, parameter :: analyticalFunctionPatch = 13
    integer, parameter :: sphericalHarmonicPatch = 14
    integer, parameter :: analytical1DPatch = 15
    integer, parameter :: bubble_pulse1D = 16
    integer, parameter :: spiralPatch = 17
    integer, parameter :: modifiedCircularPatch = 18
    integer, parameter :: modifiedCircular3DPatch = 19
    integer, parameter :: taylorGreenVortexPatch = 20
    integer, parameter :: stlPatch = 21


    integer, parameter :: gamma/pi_inf_model = 1
    integer, parameter :: five_eqn_model = 2
    integer, parameter :: six_eqn_model = 3
    integer, parameter :: four_eqn_model = 4

    integer, parameter :: grid_3Dcylindrical = 3
end module m_constants
