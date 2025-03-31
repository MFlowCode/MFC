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

    ! RKCK constants
    integer, parameter :: num_ts_rkck = 6 !< Number of time-stages in the RKCK stepper
    ! RKCK 4th/5th time stepper coefficients based on Cash J. and Karp A. (1990)
    real(wp), parameter :: rkck_c1 = 0._wp, rkck_c2 = 0.2_wp, rkck_c3 = 0.3_wp, rkck_c4 = 0.6_wp, &     ! c1 c2 c3 c4 c5 c6
                           rkck_c5 = 1._wp, rkck_c6 = 0.875_wp
    real(wp), dimension(6), parameter :: rkck_coef1 = (/0.2_wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp/)     ! a21
    real(wp), dimension(6), parameter :: rkck_coef2 = (/3._wp/40._wp, 9._wp/40._wp, 0._wp, 0._wp, &     ! a31 a32
                                                        0._wp, 0._wp/)
    real(wp), dimension(6), parameter :: rkck_coef3 = (/0.3_wp, -0.9_wp, 1.2_wp, 0._wp, 0._wp, 0._wp/)  ! a41 a42 a43
    real(wp), dimension(6), parameter :: rkck_coef4 = (/-11._wp/54._wp, 2.5_wp, -70._wp/27._wp, &       ! a51 a52 a53 a54
                                                        35._wp/27._wp, 0._wp, 0._wp/)
    real(wp), dimension(6), parameter :: rkck_coef5 = (/1631._wp/55296._wp, 175._wp/512._wp, &          ! a61 a62 a63 a64 a65
                                                        575._wp/13824._wp, 44275._wp/110592._wp, &
                                                        253._wp/4096._wp, 0._wp/)
    real(wp), dimension(6), parameter :: rkck_coef6 = (/37._wp/378._wp, 0._wp, 250._wp/621._wp, &       ! b1 b2 b3 b4 b5 b6
                                                        125._wp/594._wp, 0._wp, 512._wp/1771._wp/)
    real(wp), dimension(6), parameter :: rkck_coefE = (/37._wp/378._wp - 2825._wp/27648._wp, 0._wp, &   ! er1 er2 er3 er4 er5 er6 (4th/5th error)
                                                        250._wp/621._wp - 18575._wp/48384._wp, &
                                                        125._wp/594._wp - 13525._wp/55296._wp, &
                                                        -277._wp/14336._wp, 512._wp/1771._wp - 0.25_wp/)
    ! Adaptive rkck constants
    real(wp), parameter :: verysmall_dt = 1e-14_wp !< Very small dt, stop stepper
    real(wp), parameter :: SAFETY = 0.9_wp !< Next dt will be maximum 0.9*dt if truncation error is above tolerance.
    real(wp), parameter :: RNDDEC = 1e8_wp !< Round calculated dt (16th digit) to avoid the inclusion of random decimals
    real(wp), parameter :: PSHRNK = -0.25_wp !< Factor to reduce dt when truncation error above tolerance
    real(wp), parameter :: SHRNKDT = 0.5_wp !< Factor to reduce dt due to negative bubble radius
    real(wp), parameter :: ERRCON = 1.89e-4_wp !< Limit to slightly increase dt when truncation error is between ERRCON and 1
    real(wp), parameter :: PGROW = -0.2_wp !< Factor to increase dt when truncation error is between ERRCON and 1

    ! System constants
    integer, parameter :: CASE_FILE_ERROR_CODE = 22

end module m_constants
