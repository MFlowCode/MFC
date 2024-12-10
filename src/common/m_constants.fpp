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

    ! Lagrange bubbles constants
    integer, parameter :: mapCells = 3 !< Number of cells around the bubble where the smoothening function will have effect
    real(kind(0d0)), parameter :: R_uni = 8314 ! Universal gas constant - J/kmol/K

    ! RKCK constants
    integer, parameter :: num_ts_rkck = 6 !< Number of time-stages in the RKCK stepper
    ! RKCK 4th/5th time stepper coefficients based on Cash J. and Karp A. (1990)
    real(kind(0d0)), parameter :: rkck_c1 = 0.0d0, rkck_c2 = 0.2d0, rkck_c3 = 0.3d0, rkck_c4 = 0.6d0, &  ! c1 c2 c3 c4 c5 c6
                                  rkck_c5 = 1.0d0, rkck_c6 = 0.875d0
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef1 = (/0.2d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)   ! a21
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef2 = (/3.0d0/40.0d0, 9.0d0/40.0d0, 0.0d0, 0.0d0, &  ! a31 a32
                                                                0.0d0, 0.0d0/)
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef3 = (/0.3d0, -0.9d0, 1.2d0, 0.0d0, 0.0d0, 0.0d0/)  ! a41 a42 a43
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef4 = (/-11.0d0/54.0d0, 2.5d0, -70.0d0/27.0d0, &  ! a51 a52 a53 a54
                                                                35.d0/27.d0, 0.0d0, 0.0d0/)
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef5 = (/1631.0d0/55296.0d0, 175.0d0/512.0d0, &  ! a61 a62 a63 a64 a65
                                                                575.d0/13824.d0, 44275.d0/110592.d0, 253.d0/4096.d0, 0.0d0/)
    real(kind(0.d0)), dimension(6), parameter :: rkck_coef6 = (/37.d0/378.d0, 0.0d0, 250.d0/621.d0, &  ! b1 b2 b3 b4 b5 b6
                                                                125.0d0/594.0d0, 0.0d0, 512.0d0/1771.0d0/)
    real(kind(0.d0)), dimension(6), parameter :: rkck_coefE = (/37.d0/378.d0 - 2825.0d0/27648.0d0, 0.0d0, &  ! er1 er2 er3 er4 er5 er6 (4th/5th error)
                                                                250.d0/621.d0 - 18575.0d0/48384.0d0, 125.0d0/594.0d0 - 13525.0d0/55296.0d0, &
                                                                -277.0d0/14336.0d0, 512.0d0/1771.0d0 - 0.25d0/)
    ! Adaptive rkck constants
    real(kind(0d0)), parameter :: verysmall_dt = 1.d-14 !< Very small dt, stop stepper
    real(kind(0d0)), parameter :: SAFETY = 0.9d0 !< Next dt will be maximum 0.9*dt if truncation error is above tolerance.
    real(kind(0d0)), parameter :: RNDDEC = 1.0d05 !< Need to round the relative truncation error (5th decimal) to avoid the inclusion of random decimals when dividing by a very small number (rkck_tolerance)
    real(kind(0d0)), parameter :: PSHRNK = -0.25d0 !< Factor to reduce dt when truncation error above tolerance
    real(kind(0d0)), parameter :: SHRNKDT = 0.5d0 !< Factor to reduce dt due to negative bubble radius
    real(kind(0d0)), parameter :: ERRCON = 1.89d-4 !< Limit to slightly increase dt when truncation error is between ERRCON and 1
    real(kind(0d0)), parameter :: PGROW = -0.2d0 !< Factor to increase dt when truncation error is between ERRCON and 1

end module m_constants
