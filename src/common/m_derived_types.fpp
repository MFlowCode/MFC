!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
! New line at end of file is required for FYPP# 2 "/storage/home/hcoda1/5/tprathi3/OSPO/MFC-prathi.git/openmp/src/common/include/parallel_macros.fpp" 2
! New line at end of file is required for FYPP# 2 "/storage/home/hcoda1/5/tprathi3/OSPO/MFC-prathi.git/openmp/src/common/include/omp_macros.fpp" 2
! New line at end of file is required for FYPP# 3 "/storage/home/hcoda1/5/tprathi3/OSPO/MFC-prathi.git/openmp/src/common/include/parallel_macros.fpp" 2
! New line at end of file is required for FYPP# 2 "/storage/home/hcoda1/5/tprathi3/OSPO/MFC-prathi.git/openmp/src/common/include/acc_macros.fpp" 2
! New line at end of file is required for FYPP# 4 "/storage/home/hcoda1/5/tprathi3/OSPO/MFC-prathi.git/openmp/src/common/include/parallel_macros.fpp" 2
! New line at end of file is required for FYPP
! New line at end of file is required for FYPP
!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
module m_derived_types
    use m_constants  !< Constants
    use m_precision_select
    use m_thermochem, only: num_species
    implicit none
    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(wp), allocatable, dimension(:, :, :) :: fp !< Field position
    end type field_position
!$omp declare mapper (field_position::x) map ( &
!$omp  x%fp &
!$omp )
    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(wp), pointer, dimension(:, :, :) :: sf => null()
    end type scalar_field
!$omp declare mapper (scalar_field::x) map ( &
!$omp  x%sf &
!$omp )
    !> Derived type for bubble variables pb and mv at quadrature nodes (qbmm)
    type pres_field
        real(wp), pointer, dimension(:, :, :, :, :) :: sf => null()
    end type pres_field
!$omp declare mapper (pres_field::x) map ( &
!$omp  x%sf &
!$omp )
    !> Derived type annexing an integer scalar field (SF)
    type integer_field
        integer, pointer, dimension(:, :, :) :: sf => null()
    end type integer_field
!$omp declare mapper (integer_field::x) map ( &
!$omp  x%sf &
!$omp )
    !> Derived type for levelset
    type levelset_field
        real(wp), pointer, dimension(:, :, :, :) :: sf => null()
    end type levelset_field
!$omp declare mapper (levelset_field::x) map ( &
!$omp  x%sf &
!$omp )
    !> Derived type for levelset norm
    type levelset_norm_field
        real(wp), pointer, dimension(:, :, :, :, :) :: sf => null()
    end type levelset_norm_field
!$omp declare mapper (levelset_norm_field::x) map ( &
!$omp  x%sf &
!$omp )
    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var
!$omp declare mapper (mpi_io_var::x) map ( &
!$omp  x%view &
!$omp , x%var &
!$omp )
    type mpi_io_ib_var
        integer :: view
        type(integer_field) :: var
    end type mpi_io_ib_var
!$omp declare mapper (mpi_io_ib_var::x) map ( &
!$omp  x%view &
!$omp , x%var &
!$omp )
    type mpi_io_levelset_var
        integer :: view
        type(levelset_field) :: var
    end type mpi_io_levelset_var
!$omp declare mapper (mpi_io_levelset_var::x) map ( &
!$omp  x%view &
!$omp , x%var &
!$omp )
    type mpi_io_levelset_norm_var
        integer :: view
        type(levelset_norm_field) :: var
    end type mpi_io_levelset_norm_var
!$omp declare mapper (mpi_io_levelset_norm_var::x) map ( &
!$omp  x%view &
!$omp , x%var &
!$omp )
    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf !< Vector field
    end type vector_field
!$omp declare mapper (vector_field::x) map ( &
!$omp  x%vf &
!$omp )
    !> Generic 3-component vector (e.g., spatial coordinates or field components)
    !! Named _dt (derived types: x,y,z) to differentiate from t_vec3 (3-component vector)
    type vec3_dt ! dt for derived types
        real(wp) :: x
        real(wp) :: y
        real(wp) :: z
    end type vec3_dt
!$omp declare mapper (vec3_dt::x) map ( &
!$omp  x%x &
!$omp , x%y &
!$omp , x%z &
!$omp )
    !> Left and right Riemann states
    type riemann_states
        real(wp) :: L
        real(wp) :: R
    end type riemann_states
!$omp declare mapper (riemann_states::x) map ( &
!$omp  x%L &
!$omp , x%R &
!$omp )
    !> Left and right Riemann states for 3-component vectors
    type riemann_states_vec3
        real(wp) :: L(3)
        real(wp) :: R(3)
    end type riemann_states_vec3
!$omp declare mapper (riemann_states_vec3::x) map ( &
!$omp  x%L(3) &
!$omp , x%R(3) &
!$omp )
    !> Integer bounds for variables
    type int_bounds_info
        integer :: beg
        integer :: end
        real(wp) :: vb1
        real(wp) :: vb2
        real(wp) :: vb3
        real(wp) :: ve1
        real(wp) :: ve2
        real(wp) :: ve3
        real(wp) :: pres_in, pres_out
        real(wp), dimension(3) :: vel_in, vel_out
        real(wp), dimension(num_fluids_max) :: alpha_rho_in, alpha_in
        logical :: grcbc_in, grcbc_out, grcbc_vel_out
    end type int_bounds_info
!$omp declare mapper (int_bounds_info::x) map ( &
!$omp  x%beg &
!$omp , x%end &
!$omp , x%vb1 &
!$omp , x%vb2 &
!$omp , x%vb3 &
!$omp , x%ve1 &
!$omp , x%ve2 &
!$omp , x%ve3 &
!$omp , x%pres_in &
!$omp , x%pres_out &
!$omp , x%vel_in &
!$omp , x%vel_out &
!$omp , x%alpha_rho_in &
!$omp , x%alpha_in &
!$omp , x%grcbc_in &
!$omp , x%grcbc_out &
!$omp , x%grcbc_vel_out &
!$omp )
    type bc_patch_parameters
        integer :: geometry
        integer :: type
        integer :: dir
        integer :: loc
        real(wp), dimension(3) :: centroid
        real(wp), dimension(3) :: length
        real(wp) :: radius
    end type bc_patch_parameters
!$omp declare mapper (bc_patch_parameters::x) map ( &
!$omp  x%geometry &
!$omp , x%type &
!$omp , x%dir &
!$omp , x%loc &
!$omp , x%centroid &
!$omp , x%length &
!$omp , x%radius &
!$omp )
    !> Derived type adding beginning (beg) and end bounds info as attributes
    type bounds_info
        real(wp) :: beg
        real(wp) :: end
    end type bounds_info
!$omp declare mapper (bounds_info::x) map ( &
!$omp  x%beg &
!$omp , x%end &
!$omp )
    !> bounds for the bubble dynamic variables
    type bub_bounds_info
        integer :: beg
        integer :: end
        integer, dimension(:), allocatable :: rs
        integer, dimension(:), allocatable :: vs
        integer, dimension(:), allocatable :: ps
        integer, dimension(:), allocatable :: ms
        integer, dimension(:, :), allocatable :: moms !< Moment indices for qbmm
        integer, dimension(:, :, :), allocatable :: fullmom !< Moment indices for qbmm
    end type bub_bounds_info
!$omp declare mapper (bub_bounds_info::x) map ( &
!$omp  x%beg &
!$omp , x%end &
!$omp , x%rs &
!$omp , x%vs &
!$omp , x%ps &
!$omp , x%ms &
!$omp , x%moms &
!$omp , x%fullmom &
!$omp )
    !> Defines parameters for a Model Patch
    type ic_model_parameters
        character(LEN=pathlen_max) :: filepath !<
        !! Path the STL file relative to case_dir.
        real(wp), dimension(1:3) :: translate !<
        !! Translation of the STL object.
        real(wp), dimension(1:3) :: scale !<
        !! Scale factor for the STL object.
        real(wp), dimension(1:3) :: rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.
        integer :: spc !<
        !! Number of samples per cell to use when discretizing the STL object.
        real(wp) :: threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ic_model_parameters
!$omp declare mapper (ic_model_parameters::x) map ( &
!$omp  x%filepath &
!$omp , x%translate &
!$omp , x%scale &
!$omp , x%rotate &
!$omp , x%spc &
!$omp , x%threshold &
!$omp )
    type :: t_triangle
        real(wp), dimension(1:3, 1:3) :: v ! Vertices of the triangle
        real(wp), dimension(1:3) :: n ! Normal vector
    end type t_triangle
!$omp declare mapper (t_triangle::x) map ( &
!$omp  x%v &
!$omp , x%n &
!$omp )
    type :: t_ray
        real(wp), dimension(1:3) :: o ! Origin
        real(wp), dimension(1:3) :: d ! Direction
    end type t_ray
!$omp declare mapper (t_ray::x) map ( &
!$omp  x%o &
!$omp , x%d &
!$omp )
    type :: t_bbox
        real(wp), dimension(1:3) :: min ! Minimum coordinates
        real(wp), dimension(1:3) :: max ! Maximum coordinates
    end type t_bbox
!$omp declare mapper (t_bbox::x) map ( &
!$omp  x%min &
!$omp , x%max &
!$omp )
    type :: t_model
        integer :: ntrs   ! Number of triangles
        type(t_triangle), allocatable :: trs(:) ! Triangles
    end type t_model
!$omp declare mapper (t_model::x) map ( &
!$omp  x%ntrs &
!$omp , x%trs(:) &
!$omp )
    !> Derived type adding initial condition (ic) patch parameters as attributes
    !! NOTE: The requirements for the specification of the above parameters
    !! are strongly dependent on both the choice of the multicomponent flow
    !! model as well as the choice of the patch geometry.
    type ic_patch_parameters
        integer :: geometry !< Type of geometry for the patch
        real(wp) :: x_centroid, y_centroid, z_centroid !<
        !! Location of the geometric center, i.e. the centroid, of the patch. It
        !! is specified through its x-, y- and z-coordinates, respectively.
        real(wp) :: length_x, length_y, length_z !< Dimensions of the patch. x,y,z Lengths.
        real(wp) :: radius !< Dimensions of the patch. radius.
        real(wp), dimension(3) :: radii !<
        !! Vector indicating the various radii for the elliptical and ellipsoidal
        !! patch geometries. It is specified through its x-, y-, and z-components
        !! respectively.
        real(wp) :: epsilon, beta !<
        !! The isentropic vortex parameters for the amplitude of the disturbance and
        !! domain of influence.
        real(wp), dimension(2:9) :: a !<
        !! The parameters needed for the spherical harmonic patch
        logical :: non_axis_sym
        real(wp), dimension(3) :: normal !<
        !! Normal vector indicating the orientation of the patch. It is specified
        !! through its x-, y- and z-components, respectively.
        logical, dimension(0:num_patches_max - 1) :: alter_patch !<
        !! List of permissions that indicate to the current patch which preceding
        !! patches it is allowed to overwrite when it is in process of being laid
        !! out in the domain
        logical :: smoothen !<
        !! Permission indicating to the current patch whether its boundaries will
        !! be smoothed out across a few cells or whether they are to remain sharp
        integer :: smooth_patch_id !<
        !! Identity (id) of the patch with which current patch is to get smoothed
        real(wp) :: smooth_coeff !<
        !! Smoothing coefficient (coeff) for the size of the stencil of
        !! cells across which boundaries of the current patch will be smeared out
        real(wp), dimension(num_fluids_max) :: alpha_rho
        real(wp) :: rho
        real(wp), dimension(3) :: vel
        real(wp) :: pres
        real(wp), dimension(num_fluids_max) :: alpha
        real(wp) :: gamma
        real(wp) :: pi_inf !<
        real(wp) :: cv !<
        real(wp) :: qv !<
        real(wp) :: qvp !<
        !! Primitive variables associated with the patch. In order, these include
        !! the partial densities, density, velocity, pressure, volume fractions,
        !! specific heat ratio function and the liquid stiffness function.
        real(wp) :: Bx, By, Bz !<
        !! Magnetic field components; B%x is not used for 1D
        real(wp), dimension(6) :: tau_e !<
        !! Elastic stresses added to primitive variables if hypoelasticity = True
        real(wp) :: R0 !< Bubble size
        real(wp) :: V0 !< Bubble velocity
        real(wp) :: p0 !< Bubble size
        real(wp) :: m0 !< Bubble velocity
        integer :: hcid
        !! id for hard coded initial condition
        real(wp) :: cf_val !! color function value
        real(wp) :: Y(1:num_species)
        !! STL or OBJ model input parameter
        character(LEN=pathlen_max) :: model_filepath !<
        !! Path the STL file relative to case_dir.
        real(wp), dimension(1:3) :: model_translate !<
        !! Translation of the STL object.
        real(wp), dimension(1:3) :: model_scale !<
        !! Scale factor for the STL object.
        real(wp), dimension(1:3) :: model_rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.
        integer :: model_spc !<
        !! Number of samples per cell to use when discretizing the STL object.
        real(wp) :: model_threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ic_patch_parameters
!$omp declare mapper (ic_patch_parameters::x) map ( &
!$omp  x%geometry &
!$omp , x%x_centroid &
!$omp , x%y_centroid &
!$omp , x%z_centroid &
!$omp , x%length_x &
!$omp , x%length_y &
!$omp , x%length_z &
!$omp , x%radius &
!$omp , x%radii &
!$omp , x%epsilon &
!$omp , x%beta &
!$omp , x%a &
!$omp , x%non_axis_sym &
!$omp , x%normal &
!$omp , x%alter_patch &
!$omp , x%smoothen &
!$omp , x%smooth_patch_id &
!$omp , x%smooth_coeff &
!$omp , x%alpha_rho &
!$omp , x%rho &
!$omp , x%vel &
!$omp , x%pres &
!$omp , x%alpha &
!$omp , x%gamma &
!$omp , x%pi_inf &
!$omp , x%cv &
!$omp , x%qv &
!$omp , x%qvp &
!$omp , x%Bx &
!$omp , x%By &
!$omp , x%Bz &
!$omp , x%tau_e &
!$omp , x%R0 &
!$omp , x%V0 &
!$omp , x%p0 &
!$omp , x%m0 &
!$omp , x%hcid &
!$omp , x%cf_val &
!$omp , x%Y(1:num_species) &
!$omp , x%model_filepath &
!$omp , x%model_translate &
!$omp , x%model_scale &
!$omp , x%model_rotate &
!$omp , x%model_spc &
!$omp , x%model_threshold &
!$omp )
    type ib_patch_parameters
        integer :: geometry !< Type of geometry for the patch
        real(wp) :: x_centroid, y_centroid, z_centroid !<
        !! Location of the geometric center, i.e. the centroid, of the patch. It
        !! is specified through its x-, y- and z-coordinates, respectively.
        real(wp) :: c, p, t, m
        real(wp) :: length_x, length_y, length_z !< Dimensions of the patch. x,y,z Lengths.
        real(wp) :: radius !< Dimensions of the patch. radius.
        real(wp) :: theta
        logical :: slip
        !! STL or OBJ model input parameter
        character(LEN=pathlen_max) :: model_filepath !<
        !! Path the STL file relative to case_dir.
        real(wp), dimension(1:3) :: model_translate !<
        !! Translation of the STL object.
        real(wp), dimension(1:3) :: model_scale !<
        !! Scale factor for the STL object.
        real(wp), dimension(1:3) :: model_rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.
        integer :: model_spc !<
        !! Number of samples per cell to use when discretizing the STL object.
        real(wp) :: model_threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ib_patch_parameters
!$omp declare mapper (ib_patch_parameters::x) map ( &
!$omp  x%geometry &
!$omp , x%x_centroid &
!$omp , x%y_centroid &
!$omp , x%z_centroid &
!$omp , x%c &
!$omp , x%p &
!$omp , x%t &
!$omp , x%m &
!$omp , x%length_x &
!$omp , x%length_y &
!$omp , x%length_z &
!$omp , x%radius &
!$omp , x%theta &
!$omp , x%slip &
!$omp , x%model_filepath &
!$omp , x%model_translate &
!$omp , x%model_scale &
!$omp , x%model_rotate &
!$omp , x%model_spc &
!$omp , x%model_threshold &
!$omp )
    !> Derived type annexing the physical parameters (PP) of the fluids. These
    !! include the specific heat ratio function and liquid stiffness function.
    type physical_parameters
        real(wp) :: gamma   !< Sp. heat ratio
        real(wp) :: pi_inf  !< Liquid stiffness
        real(wp), dimension(2) :: Re      !< Reynolds number
        real(wp) :: cv      !< heat capacity
        real(wp) :: qv      !< reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
        real(wp) :: qvp     !< reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
        real(wp) :: mul0    !< Bubble viscosity
        real(wp) :: ss      !< Bubble surface tension
        real(wp) :: pv      !< Bubble vapour pressure
        real(wp) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: M_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: mu_v    !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: k_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: cp_v
        real(wp) :: G
    end type physical_parameters
!$omp declare mapper (physical_parameters::x) map ( &
!$omp  x%gamma &
!$omp , x%pi_inf &
!$omp , x%Re &
!$omp , x%cv &
!$omp , x%qv &
!$omp , x%qvp &
!$omp , x%mul0 &
!$omp , x%ss &
!$omp , x%pv &
!$omp , x%gamma_v &
!$omp , x%M_v &
!$omp , x%mu_v &
!$omp , x%k_v &
!$omp , x%cp_v &
!$omp , x%G &
!$omp )
    type mpi_io_airfoil_ib_var
        integer, dimension(2) :: view
        type(vec3_dt), allocatable, dimension(:) :: var
    end type mpi_io_airfoil_ib_var
!$omp declare mapper (mpi_io_airfoil_ib_var::x) map ( &
!$omp  x%view &
!$omp , x%var &
!$omp )
    !> Derived type annexing integral regions
    type integral_parameters
        real(wp) :: xmin !< Min. boundary first coordinate direction
        real(wp) :: xmax !< Max. boundary first coordinate direction
        real(wp) :: ymin !< Min. boundary second coordinate direction
        real(wp) :: ymax !< Max. boundary second coordinate direction
        real(wp) :: zmin !< Min. boundary third coordinate direction
        real(wp) :: zmax !< Max. boundary third coordinate direction
    end type integral_parameters
!$omp declare mapper (integral_parameters::x) map ( &
!$omp  x%xmin &
!$omp , x%xmax &
!$omp , x%ymin &
!$omp , x%ymax &
!$omp , x%zmin &
!$omp , x%zmax &
!$omp )
    !> Acoustic source parameters
    type acoustic_parameters
        integer :: pulse !< Type of pulse
        integer :: support !< Type of support
        logical :: dipole !< Whether the source is a dipole or monopole
        real(wp), dimension(3) :: loc !< Physical location of acoustic source
        real(wp) :: mag !< Acoustic pulse magnitude
        real(wp) :: length !< Length of planar source (2D/3D)
        real(wp) :: height !< Height of planar source (3D)
        real(wp) :: wavelength !< Wave length of pulse
        real(wp) :: frequency !< Frequency of pulse
        real(wp) :: gauss_sigma_dist !< sigma of Gaussian pulse multiplied by speed of sound
        real(wp) :: gauss_sigma_time !< sigma of Gaussian pulse
        real(wp) :: npulse !< Number of cycles of pulse
        real(wp) :: dir !< Direction of pulse
        real(wp) :: delay !< Time-delay of pulse start
        real(wp) :: foc_length ! < Focal length of transducer
        real(wp) :: aperture ! < Aperture diameter of transducer
        real(wp) :: element_spacing_angle !< Spacing between aperture elements in 2D acoustic array
        real(wp) :: element_polygon_ratio !< Ratio of aperture element diameter to side length of polygon connecting their centers, in 3D acoustic array
        real(wp) :: rotate_angle !< Angle of rotation of the entire circular 3D acoustic array
        real(wp) :: bb_bandwidth !< Bandwidth of each frequency in broadband wave
        real(wp) :: bb_lowest_freq !< The lower frequency bound of broadband wave
        integer :: num_elements !< Number of elements in the acoustic array
        integer :: element_on !< Element in the acoustic array to turn on
        integer :: bb_num_freq !< Number of frequencies in the broadband wave
    end type acoustic_parameters
!$omp declare mapper (acoustic_parameters::x) map ( &
!$omp  x%pulse &
!$omp , x%support &
!$omp , x%dipole &
!$omp , x%loc &
!$omp , x%mag &
!$omp , x%length &
!$omp , x%height &
!$omp , x%wavelength &
!$omp , x%frequency &
!$omp , x%gauss_sigma_dist &
!$omp , x%gauss_sigma_time &
!$omp , x%npulse &
!$omp , x%dir &
!$omp , x%delay &
!$omp , x%foc_length &
!$omp , x%aperture &
!$omp , x%element_spacing_angle &
!$omp , x%element_polygon_ratio &
!$omp , x%rotate_angle &
!$omp , x%bb_bandwidth &
!$omp , x%bb_lowest_freq &
!$omp , x%num_elements &
!$omp , x%element_on &
!$omp , x%bb_num_freq &
!$omp )
    !> Acoustic source source_spatial pre-calculated values
    type source_spatial_type
        integer, dimension(:, :), allocatable :: coord !< List of grid points indices with non-zero source_spatial values
        real(wp), dimension(:), allocatable :: val !< List of non-zero source_spatial values
        real(wp), dimension(:), allocatable :: angle !< List of angles with x-axis for mom source term vector
        real(wp), dimension(:, :), allocatable :: xyz_to_r_ratios !< List of [xyz]/r for mom source term vector
    end type source_spatial_type
!$omp declare mapper (source_spatial_type::x) map ( &
!$omp  x%coord &
!$omp , x%val &
!$omp , x%angle &
!$omp , x%xyz_to_r_ratios &
!$omp )
    !> Ghost Point for Immersed Boundaries
    type ghost_point
        integer, dimension(3) :: loc !< Physical location of the ghost point
        real(wp), dimension(3) :: ip_loc !< Physical location of the image point
        integer, dimension(3) :: ip_grid !< Top left grid point of IP
        real(wp), dimension(2, 2, 2) :: interp_coeffs !< Interpolation Coefficients of image point
        integer :: ib_patch_id !< ID of the IB Patch the ghost point is part of
        logical :: slip
        integer, dimension(3) :: DB
    end type ghost_point
!$omp declare mapper (ghost_point::x) map ( &
!$omp  x%loc &
!$omp , x%ip_loc &
!$omp , x%ip_grid &
!$omp , x%interp_coeffs &
!$omp , x%ib_patch_id &
!$omp , x%slip &
!$omp , x%DB &
!$omp )
    !> Species parameters
    type species_parameters
        character(LEN=name_len) :: name !< Name of species
    end type species_parameters
!$omp declare mapper (species_parameters::x) map ( &
!$omp  x%name &
!$omp )
    !> Chemistry parameters
    type chemistry_parameters
        character(LEN=name_len) :: cantera_file !< Path to Cantera file
        logical :: diffusion
        logical :: reactions
        !> Method of determining gamma.
        !> gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
        !> gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
        integer :: gamma_method
    end type chemistry_parameters
!$omp declare mapper (chemistry_parameters::x) map ( &
!$omp  x%cantera_file &
!$omp , x%diffusion &
!$omp , x%reactions &
!$omp , x%gamma_method &
!$omp )
    !> Lagrangian bubble parameters
    type bubbles_lagrange_parameters
        integer :: solver_approach          !< 1: One-way coupling, 2: two-way coupling
        integer :: cluster_type             !< Cluster model to find p_inf
        logical :: pressure_corrector       !< Cell pressure correction term
        integer :: smooth_type              !< Smoothing function. 1: Gaussian, 2:Delta 3x3
        logical :: heatTransfer_model       !< Activate HEAT transfer model at the bubble-liquid interface
        logical :: massTransfer_model       !< Activate MASS transfer model at the bubble-liquid interface
        logical :: write_bubbles            !< Write files to track the bubble evolution each time step
        logical :: write_bubbles_stats      !< Write the maximum and minimum radius of each bubble
        integer :: nBubs_glb                !< Global number of bubbles
        real(wp) :: epsilonb         !< Standard deviation scaling for the gaussian function
        real(wp) :: charwidth        !< Domain virtual depth (z direction, for 2D simulations)
        real(wp) :: valmaxvoid       !< Maximum void fraction permitted
        real(wp) :: c0               !< Reference speed
        real(wp) :: rho0             !< Reference density
        real(wp) :: T0, Thost        !< Reference temperature and host temperature
        real(wp) :: x0               !< Reference length
        real(wp) :: diffcoefvap      !< Vapor diffusivity in the gas
    end type bubbles_lagrange_parameters
!$omp declare mapper (bubbles_lagrange_parameters::x) map ( &
!$omp  x%solver_approach &
!$omp , x%cluster_type &
!$omp , x%pressure_corrector &
!$omp , x%smooth_type &
!$omp , x%heatTransfer_model &
!$omp , x%massTransfer_model &
!$omp , x%write_bubbles &
!$omp , x%write_bubbles_stats &
!$omp , x%nBubs_glb &
!$omp , x%epsilonb &
!$omp , x%charwidth &
!$omp , x%valmaxvoid &
!$omp , x%c0 &
!$omp , x%rho0 &
!$omp , x%T0 &
!$omp , x%Thost &
!$omp , x%x0 &
!$omp , x%diffcoefvap &
!$omp )
    !> Max and min number of cells in a direction of each combination of x-,y-, and z-
    type cell_num_bounds
        integer :: mn_max, np_max, mp_max, mnp_max
        integer :: mn_min, np_min, mp_min, mnp_min
    end type cell_num_bounds
!$omp declare mapper (cell_num_bounds::x) map ( &
!$omp  x%mn_max &
!$omp , x%np_max &
!$omp , x%mp_max &
!$omp , x%mnp_max &
!$omp , x%mn_min &
!$omp , x%np_min &
!$omp , x%mp_min &
!$omp , x%mnp_min &
!$omp )
end module m_derived_types

! Code was translated using: /media/shared/Documents/GitHub/OSPO/intel-application-migration-tool-for-openacc-to-openmp/simulation/src/intel-application-migration-tool-for-openacc-to-openmp -keep-binding-clauses=all simulation/p_main.fpp.f90
