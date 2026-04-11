!>
!! @file
!! @brief Contains module m_derived_types

#:include "macros.fpp"

!> @brief Shared derived types for field data, patch geometry, bubble dynamics, and MPI I/O structures
module m_derived_types

    use m_constants
    use m_precision_select
    use m_thermochem, only: num_species

    implicit none

    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(stp), allocatable, dimension(:,:,:) :: fp  !< Field position
    end type field_position

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(stp), pointer, dimension(:,:,:) :: sf => null()
    end type scalar_field

    !> Derived type for bubble variables pb and mv at quadrature nodes (qbmm)
    type pres_field
        real(stp), pointer, dimension(:,:,:,:,:) :: sf => null()
    end type pres_field

    !> Derived type annexing an integer scalar field (SF)
    type integer_field
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), pointer, dimension(:,:,:) :: sf => null()
#else
        integer, pointer, dimension(:,:,:) :: sf => null()
#endif
    end type integer_field

    !> Derived type for levelset
    type levelset_field
        real(stp), pointer, dimension(:,:,:,:) :: sf => null()
    end type levelset_field

    !> Derived type for levelset norm
    type levelset_norm_field
        real(stp), pointer, dimension(:,:,:,:,:) :: sf => null()
    end type levelset_norm_field

    type mpi_io_var
        integer, allocatable, dimension(:)            :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    type mpi_io_ib_var
        integer             :: view
        type(integer_field) :: var
    end type mpi_io_ib_var

    type mpi_io_levelset_var
        integer              :: view
        type(levelset_field) :: var
    end type mpi_io_levelset_var

    type mpi_io_levelset_norm_var
        integer                   :: view
        type(levelset_norm_field) :: var
    end type mpi_io_levelset_norm_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf  !< Vector field
    end type vector_field

    !> Generic 3-component vector (e.g., spatial coordinates or field components) Named _dt (derived types: x,y,z) to differentiate
    !! from t_vec3 (3-component vector)
    type vec3_dt  ! dt for derived types
        real(wp) :: x
        real(wp) :: y
        real(wp) :: z
    end type vec3_dt

    !> Left and right Riemann states
    type riemann_states
        real(wp) :: L
        real(wp) :: R
    end type riemann_states

    !> Left and right Riemann states for 3-component vectors
    type riemann_states_vec3
        real(wp) :: L(3)
        real(wp) :: R(3)
    end type riemann_states_vec3

    !> Lightweight beg/end pair for equation index ranges (no BC payload).
    type idx_bounds_info
        integer :: beg
        integer :: end
    end type idx_bounds_info

    !> Integer bounds for variables
    type int_bounds_info
        integer                             :: beg
        integer                             :: end
        real(wp)                            :: vb1
        real(wp)                            :: vb2
        real(wp)                            :: vb3
        real(wp)                            :: ve1
        real(wp)                            :: ve2
        real(wp)                            :: ve3
        real(wp)                            :: pres_in, pres_out
        real(wp), dimension(3)              :: vel_in, vel_out
        real(wp), dimension(num_fluids_max) :: alpha_rho_in, alpha_in
        logical                             :: grcbc_in, grcbc_out, grcbc_vel_out
        logical                             :: isothermal_in, isothermal_out
        real(wp)                            :: Twall_in, Twall_out
    end type int_bounds_info

    !> Groups the x, y, z boundary condition begin/end codes for passing as a single argument.
    type bc_xyz_info
        type(int_bounds_info) :: x, y, z
    end type bc_xyz_info

    !> QBMM moment index mappings - separate from bub beg/end so eqn_idx contains no allocatables.
    type qbmm_idx_info
        integer, dimension(:), allocatable     :: rs       !< R moment indices per bubble bin
        integer, dimension(:), allocatable     :: vs       !< V moment indices per bubble bin
        integer, dimension(:), allocatable     :: ps       !< Pressure moment indices per bubble bin
        integer, dimension(:), allocatable     :: ms       !< Mass moment indices per bubble bin
        integer, dimension(:,:), allocatable   :: moms     !< Moment indices for qbmm
        integer, dimension(:,:,:), allocatable :: fullmom  !< Full moment indices for qbmm
    end type qbmm_idx_info

    !> All conserved-variable equation indices, computed at startup from model_eqns and enabled features.
    !> Range indices (beg/end) use int_bounds_info; scalar indices are plain integers (0 = inactive).
    !> Contains no allocatable members - safe for GPU_DECLARE as a single struct.
    type eqn_idx_info
        type(idx_bounds_info) :: cont     !< Partial densities (continuity equations)
        type(idx_bounds_info) :: mom      !< Momentum components
        type(idx_bounds_info) :: adv      !< Volume fractions (advection equations)
        type(idx_bounds_info) :: bub      !< Bubble equation range (beg/end only)
        type(idx_bounds_info) :: stress   !< Stress tensor components
        type(idx_bounds_info) :: xi       !< Reference map equations
        type(idx_bounds_info) :: B        !< Magnetic field components
        type(idx_bounds_info) :: int_en   !< Internal energy equations
        type(idx_bounds_info) :: species  !< Chemistry species equations
        integer               :: E        !< Energy/pressure equation
        integer               :: n        !< Number density equation
        integer               :: alf      !< Void fraction (scalar, model_eqns=4)
        integer               :: gamma    !< Specific heat ratio function (model_eqns=1)
        integer               :: pi_inf   !< Liquid stiffness function (model_eqns=1)
        integer               :: c        !< Color function equation
        integer               :: damage   !< Damage variable equation
        integer               :: psi      !< Psi variable equation
    end type eqn_idx_info

    type bc_patch_parameters
        integer                :: geometry
        integer                :: type
        integer                :: dir
        integer                :: loc
        real(wp), dimension(3) :: centroid
        real(wp), dimension(3) :: length
        real(wp)               :: radius
    end type bc_patch_parameters

    !> Derived type adding beginning (beg) and end bounds info as attributes
    type bounds_info
        real(wp) :: beg
        real(wp) :: end
    end type bounds_info

    !> Defines parameters for a Model Patch
    type ic_model_parameters
        character(LEN=pathlen_max) :: filepath   !< Path the STL file relative to case_dir.
        real(wp), dimension(1:3)   :: translate  !< Translation of the STL object.
        real(wp), dimension(1:3)   :: scale      !< Scale factor for the STL object.
        real(wp), dimension(1:3)   :: rotate     !< Angle to rotate the STL object along each cartesian coordinate axis, in radians.
        integer                    :: spc        !< Number of samples per cell to use when discretizing the STL object.
        real(wp)                   :: threshold  !< Threshold to turn on smoothen STL patch.
    end type ic_model_parameters

    type :: t_triangle
        real(wp), dimension(1:3,1:3) :: v  !< Vertices of the triangle
        real(wp), dimension(1:3)     :: n  !< Normal vector
    end type t_triangle

    type :: t_ray
        real(wp), dimension(1:3) :: o  !< Origin
        real(wp), dimension(1:3) :: d  !< Direction
    end type t_ray

    type :: t_bbox
        real(wp), dimension(1:3) :: min  !< Minimum coordinates
        real(wp), dimension(1:3) :: max  !< Maximum coordinates
    end type t_bbox

    type :: t_model
        integer                       :: ntrs    !< Number of triangles
        type(t_triangle), allocatable :: trs(:)  !< Triangles
    end type t_model

    type :: t_model_array
        ! Original CPU-side fields (unchanged)
        type(t_model), allocatable              :: model                    !< STL/OBJ geometry model
        real(wp), allocatable, dimension(:,:,:) :: boundary_v               !< Boundary vertices
        real(wp), allocatable, dimension(:,:)   :: interpolated_boundary_v  !< Interpolated boundary vertices
        integer                                 :: boundary_edge_count      !< Number of boundary edges
        integer                                 :: total_vertices           !< Total vertex count
        integer                                 :: interpolate              !< Interpolation flag

        ! GPU-friendly flattened arrays
        integer                                 :: ntrs   !< Copy of model%ntrs
        real(wp), allocatable, dimension(:,:,:) :: trs_v  !< Triangle vertices (3, 3, ntrs)
        real(wp), allocatable, dimension(:,:)   :: trs_n  !< Triangle normals (3, ntrs)
    end type t_model_array

    !> Derived type adding initial condition (ic) patch parameters as attributes NOTE: The requirements for the specification of the
    !! above parameters are strongly dependent on both the choice of the multicomponent flow model as well as the choice of the
    !! patch geometry.
    type ic_patch_parameters

        integer :: geometry  !< Type of geometry for the patch
        real(wp) :: x_centroid, y_centroid, z_centroid  !< Geometric center coordinates of the patch
        real(wp) :: length_x, length_y, length_z  !< Dimensions of the patch. x,y,z Lengths.
        real(wp) :: radius  !< Dimensions of the patch. radius.
        real(wp), dimension(3) :: radii  !< Elliptical/ellipsoidal patch radii in x, y, z
        real(wp) :: epsilon, beta  !< The isentropic vortex parameters for the amplitude of the disturbance and domain of influence.
        real(wp), dimension(2:9) :: a  !< Used by hardcoded IC and as temporary variables.
        logical :: non_axis_sym

        ! Geometry 13 (2D modal Fourier): fourier_cos(n), fourier_sin(n) for mode n
        real(wp), dimension(1:max_2d_fourier_modes) :: fourier_cos, fourier_sin
        logical :: modal_clip_r_to_min  !< When true, clip boundary radius: R(theta) = max(R(theta), modal_r_min) (Non-exp form only)
        real(wp) :: modal_r_min  !< Minimum boundary radius when modal_clip_r_to_min is true (Non-exp form only)
        logical :: modal_use_exp_form  !< When true, boundary = radius*exp(Fourier series)
        ! Geometry 14 (3D spherical harmonic): sph_har_coeff(l,m) for real Y_lm
        real(wp), dimension(0:max_sph_harm_degree,-max_sph_harm_degree:max_sph_harm_degree) :: sph_har_coeff
        real(wp), dimension(3) :: normal  !< Patch orientation normal vector (x, y, z)
        logical, dimension(0:num_patches_max - 1) :: alter_patch  !< Overwrite permissions for preceding patches
        logical :: smoothen  !< Whether patch boundaries are smoothed across cells
        integer :: smooth_patch_id  !< Identity (id) of the patch with which current patch is to get smoothed
        real(wp) :: smooth_coeff  !< Smoothing stencil size coefficient
        real(wp), dimension(num_fluids_max) :: alpha_rho
        real(wp) :: rho
        real(wp), dimension(3) :: vel
        real(wp) :: pres
        real(wp), dimension(num_fluids_max) :: alpha
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: cv
        real(wp) :: qv
        real(wp) :: qvp  !< Reference entropy per unit mass (SGEOS)
        real(wp) :: Bx, By, Bz  !< Magnetic field components; B%x is not used for 1D
        real(wp), dimension(6) :: tau_e  !< Elastic stresses added to primitive variables if hypoelasticity = True
        real(wp) :: R0  !< Bubble size
        real(wp) :: V0  !< Bubble velocity
        real(wp) :: p0  !< Bubble size
        real(wp) :: m0  !< Bubble velocity
        integer :: hcid  !< Hardcoded initial condition ID
        real(wp) :: cf_val  !< Color function value
        real(wp) :: Y(1:num_species)  !< Species mass fractions

        ! STL or OBJ model input parameter
        character(LEN=pathlen_max) :: model_filepath   !< Path the STL file relative to case_dir.
        real(wp), dimension(1:3)   :: model_translate  !< Translation of the STL object.
        real(wp), dimension(1:3)   :: model_scale      !< Scale factor for the STL object.
        !> Angle to rotate the STL object along each cartesian coordinate axis, in radians.
        real(wp), dimension(1:3) :: model_rotate
        integer                  :: model_spc        !< Number of samples per cell to use when discretizing the STL object.
        real(wp)                 :: model_threshold  !< Threshold to turn on smoothen STL patch.
    end type ic_patch_parameters

    type ib_patch_parameters

        integer  :: geometry                            !< Type of geometry for the patch
        real(wp) :: x_centroid, y_centroid, z_centroid  !< Geometric center coordinates of the patch
        !> Centroid locations of intermediate steps in the time_stepper module
        real(wp)                 :: step_x_centroid, step_y_centroid, step_z_centroid
        real(wp), dimension(1:3) :: centroid_offset  !< offset of center of mass from computed cell center for odd-shaped IBs
        real(wp), dimension(1:3) :: angles
        real(wp), dimension(1:3) :: step_angles
        !> matrix that converts from IB reference frame to fluid reference frame
        real(wp), dimension(1:3,1:3) :: rotation_matrix
        !> matrix that converts from fluid reference frame to IB reference frame
        real(wp), dimension(1:3,1:3) :: rotation_matrix_inverse
        real(wp)                     :: c, p, t, m
        real(wp)                     :: length_x, length_y, length_z  !< Dimensions of the patch. x,y,z Lengths.
        real(wp)                     :: radius                        !< Dimensions of the patch. radius.
        real(wp)                     :: theta
        logical                      :: slip

        ! STL or OBJ model input parameter
        character(LEN=pathlen_max) :: model_filepath   !< Path the STL file relative to case_dir.
        real(wp), dimension(1:3)   :: model_translate  !< Translation of the STL object.
        real(wp), dimension(1:3)   :: model_scale      !< Scale factor for the STL object.
        !> Angle to rotate the STL object along each cartesian coordinate axis, in radians.
        real(wp), dimension(1:3) :: model_rotate
        integer :: model_spc  !< Number of samples per cell to use when discretizing the STL object.
        real(wp) :: model_threshold  !< Threshold to turn on smoothen STL patch. Patch conditions for moving imersed boundaries
        integer :: moving_ibm  !< 0 for no moving, 1 for moving, 2 for moving on forced path
        real(wp) :: mass, moment  !< mass and moment of inertia of object used to compute forces in 2-way coupling
        real(wp), dimension(1:3) :: force, torque  !< vectors for the computed force and torque values applied to an IB
        real(wp), dimension(1:3) :: vel
        real(wp), dimension(1:3) :: step_vel  !< velocity array used to store intermediate steps in the time_stepper module
        real(wp), dimension(1:3) :: angular_vel
        real(wp), dimension(1:3) :: step_angular_vel  !< velocity array used to store intermediate steps in the time_stepper module
    end type ib_patch_parameters

    !> Derived type annexing the physical parameters (PP) of the fluids. These include the specific heat ratio function and liquid
    !! stiffness function.
    type physical_parameters
        real(wp)               :: gamma   !< Sp. heat ratio
        real(wp)               :: pi_inf  !< Liquid stiffness
        real(wp), dimension(2) :: Re      !< Reynolds number
        real(wp)               :: cv      !< heat capacity
        real(wp)               :: qv      !< reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
        real(wp)               :: qvp     !< reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
        real(wp)               :: G
    end type physical_parameters

    !> Derived type annexing the physical parameters required for sub-grid bubble models
    type subgrid_bubble_physical_parameters
        real(wp) :: R0ref    !< reference bubble radius
        real(wp) :: p0ref    !< reference pressure
        real(wp) :: rho0ref  !< reference density
        real(wp) :: T0ref    !< reference temperature
        real(wp) :: ss       !< surface tension between host and gas (bubble)
        real(wp) :: pv       !< vapor pressure of host
        real(wp) :: vd       !< vapor diffusivity in gas (bubble)
        real(wp) :: mu_l     !< viscosity of host in liquid state
        real(wp) :: mu_v     !< viscosity of host in vapor state
        real(wp) :: mu_g     !< viscosity of gas (bubble)
        real(wp) :: gam_v    !< specific heat ratio of host in vapor state
        real(wp) :: gam_g    !< specific heat ratio of gas (bubble)
        real(wp) :: M_v      !< Molecular weight of host
        real(wp) :: M_g      !< Molecular weight of gas (bubble)
        real(wp) :: k_v      !< thermal conductivity of host in vapor state
        real(wp) :: k_g      !< thermal conductivity of gas (bubble)
        real(wp) :: cp_v     !< specific heat capacity in constant pressure of host in vapor state
        real(wp) :: cp_g     !< specific heat capacity in constant pressure of gas (bubble)
        real(wp) :: R_v      !< gas constant of host in vapor state
        real(wp) :: R_g      !< gas constant of gas (bubble)
    end type subgrid_bubble_physical_parameters

    type mpi_io_airfoil_ib_var
        integer, dimension(2)                    :: view
        type(vec3_dt), allocatable, dimension(:) :: var
    end type mpi_io_airfoil_ib_var

    !> Derived type annexing integral regions
    type integral_parameters
        real(wp) :: xmin  !< Min. boundary first coordinate direction
        real(wp) :: xmax  !< Max. boundary first coordinate direction
        real(wp) :: ymin  !< Min. boundary second coordinate direction
        real(wp) :: ymax  !< Max. boundary second coordinate direction
        real(wp) :: zmin  !< Min. boundary third coordinate direction
        real(wp) :: zmax  !< Max. boundary third coordinate direction
    end type integral_parameters

    !> Acoustic source parameters
    type acoustic_parameters
        integer                :: pulse                  !< Type of pulse
        integer                :: support                !< Type of support
        logical                :: dipole                 !< Whether the source is a dipole or monopole
        real(wp), dimension(3) :: loc                    !< Physical location of acoustic source
        real(wp)               :: mag                    !< Acoustic pulse magnitude
        real(wp)               :: length                 !< Length of planar source (2D/3D)
        real(wp)               :: height                 !< Height of planar source (3D)
        real(wp)               :: wavelength             !< Wave length of pulse
        real(wp)               :: frequency              !< Frequency of pulse
        real(wp)               :: gauss_sigma_dist       !< sigma of Gaussian pulse multiplied by speed of sound
        real(wp)               :: gauss_sigma_time       !< sigma of Gaussian pulse
        real(wp)               :: npulse                 !< Number of cycles of pulse
        real(wp)               :: dir                    !< Direction of pulse
        real(wp)               :: delay                  !< Time-delay of pulse start
        real(wp)               :: foc_length             !< Focal length of transducer
        real(wp)               :: aperture               !< Aperture diameter of transducer
        real(wp)               :: element_spacing_angle  !< Spacing between aperture elements in 2D acoustic array
        !> Ratio of aperture element diameter to side length of polygon connecting their centers, in 3D acoustic array
        real(wp) :: element_polygon_ratio
        real(wp) :: rotate_angle    !< Angle of rotation of the entire circular 3D acoustic array
        real(wp) :: bb_bandwidth    !< Bandwidth of each frequency in broadband wave
        real(wp) :: bb_lowest_freq  !< The lower frequency bound of broadband wave
        integer  :: num_elements    !< Number of elements in the acoustic array
        integer  :: element_on      !< Element in the acoustic array to turn on
        integer  :: bb_num_freq     !< Number of frequencies in the broadband wave
    end type acoustic_parameters

    !> Acoustic source source_spatial pre-calculated values
    type source_spatial_type
        integer, pointer, dimension(:,:)  :: coord => null()  !< List of grid points indices with non-zero source_spatial values
        real(wp), pointer, dimension(:)   :: val => null()  !< List of non-zero source_spatial values
        real(wp), pointer, dimension(:)   :: angle => null()  !< List of angles with x-axis for mom source term vector
        real(wp), pointer, dimension(:,:) :: xyz_to_r_ratios => null()  !< List of [xyz]/r for mom source term vector
    end type source_spatial_type

    !> Ghost Point for Immersed Boundaries
    type ghost_point
        integer, dimension(3)        :: loc            !< Physical location of the ghost point
        real(wp), dimension(3)       :: ip_loc         !< Physical location of the image point
        integer, dimension(3)        :: ip_grid        !< Top left grid point of IP
        real(wp), dimension(2, 2, 2) :: interp_coeffs  !< Interpolation Coefficients of image point
        integer                      :: ib_patch_id    !< ID of the IB Patch the ghost point is part of
        real(wp)                     :: levelset
        real(wp), dimension(1:3)     :: levelset_norm
        logical                      :: slip
        integer, dimension(3)        :: DB
        integer                      :: x_periodicity, y_periodicity, z_periodicity
    end type ghost_point

    !> Species parameters
    type species_parameters
        character(LEN=name_len) :: name  !< Name of species
    end type species_parameters

    !> Chemistry parameters
    type chemistry_parameters
        character(LEN=name_len) :: cantera_file  !< Path to Cantera file
        logical                 :: diffusion
        logical                 :: reactions

        !> Method of determining gamma.
        !> gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
        !> gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
        integer :: gamma_method
        integer :: transport_model
    end type chemistry_parameters

    !> Lagrangian bubble parameters
    type bubbles_lagrange_parameters

        integer  :: solver_approach      !< 1: One-way coupling, 2: two-way coupling
        integer  :: cluster_type         !< Cluster model to find p_inf
        logical  :: pressure_corrector   !< Cell pressure correction term
        integer  :: smooth_type          !< Smoothing function. 1: Gaussian, 2:Delta 3x3
        logical  :: heatTransfer_model   !< Activate HEAT transfer model at the bubble-liquid interface
        logical  :: massTransfer_model   !< Activate MASS transfer model at the bubble-liquid interface
        logical  :: write_bubbles        !< Write files to track the bubble evolution each time step
        logical  :: write_bubbles_stats  !< Write the maximum and minimum radius of each bubble
        integer  :: nBubs_glb            !< Global number of bubbles
        real(wp) :: epsilonb             !< Standard deviation scaling for the gaussian function
        real(wp) :: charwidth            !< Domain virtual depth (z direction, for 2D simulations)
        real(wp) :: valmaxvoid           !< Maximum void fraction permitted
    end type bubbles_lagrange_parameters

    !> Max and min number of cells in a direction of each combination of x-,y-, and z-
    type cell_num_bounds
        integer :: mn_max, np_max, mp_max, mnp_max
        integer :: mn_min, np_min, mp_min, mnp_min
    end type cell_num_bounds

    type simplex_noise_params
        logical, dimension(3)                   :: perturb_vel
        real(wp), dimension(3)                  :: perturb_vel_freq
        real(wp), dimension(3)                  :: perturb_vel_scale
        real(wp), dimension(3, 3)               :: perturb_vel_offset
        logical, dimension(1:num_fluids_max)    :: perturb_dens
        real(wp), dimension(1:num_fluids_max)   :: perturb_dens_freq
        real(wp), dimension(1:num_fluids_max)   :: perturb_dens_scale
        real(wp), dimension(1:num_fluids_max,3) :: perturb_dens_offset
    end type simplex_noise_params
end module m_derived_types
