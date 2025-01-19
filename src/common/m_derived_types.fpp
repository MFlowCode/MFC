!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types

#:include "macros.fpp"

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

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(wp), pointer, dimension(:, :, :) :: sf => null()
    end type scalar_field

    !> Derived type for bubble variables pb and mv at quadrature nodes (qbmm)
    type pres_field
        real(wp), pointer, dimension(:, :, :, :, :) :: sf => null()
    end type pres_field

    !> Derived type annexing an integer scalar field (SF)
    type integer_field
        integer, pointer, dimension(:, :, :) :: sf => null()
    end type integer_field

    !> Derived type for levelset
    type levelset_field
        real(wp), pointer, dimension(:, :, :, :) :: sf => null()
    end type levelset_field

    !> Derived type for levelset norm
    type levelset_norm_field
        real(wp), pointer, dimension(:, :, :, :, :) :: sf => null()
    end type levelset_norm_field

    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    type mpi_io_ib_var
        integer :: view
        type(integer_field) :: var
    end type mpi_io_ib_var

    type mpi_io_levelset_var
        integer :: view
        type(levelset_field) :: var
    end type mpi_io_levelset_var

    type mpi_io_levelset_norm_var
        integer :: view
        type(levelset_norm_field) :: var
    end type mpi_io_levelset_norm_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf !< Vector field
    end type vector_field

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

    !> Derived type adding beginning (beg) and end bounds info as attributes
    type bounds_info
        real(wp) :: beg
        real(wp) :: end
    end type bounds_info

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

    !> Defines parameters for a Model Patch
    type ic_model_parameters
        character(LEN=pathlen_max) :: filepath !<
        !! Path the STL file relative to case_dir.

        t_vec3 :: translate !<
        !! Translation of the STL object.

        t_vec3 :: scale !<
        !! Scale factor for the STL object.

        t_vec3 :: rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.

        integer :: spc !<
        !! Number of samples per cell to use when discretizing the STL object.

        real(wp) :: threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ic_model_parameters

    type :: t_triangle
        real(wp), dimension(1:3, 1:3) :: v ! Vertices of the triangle
        t_vec3 :: n ! Normal vector
    end type t_triangle

    type :: t_ray
        t_vec3 :: o ! Origin
        t_vec3 :: d ! Direction
    end type t_ray

    type :: t_bbox
        t_vec3 :: min ! Minimum coordinates
        t_vec3 :: max ! Maximum coordinates
    end type t_bbox

    type :: t_model
        integer :: ntrs   ! Number of triangles
        type(t_triangle), allocatable :: trs(:) ! Triangles
    end type t_model

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

        t_vec3 :: model_translate !<
        !! Translation of the STL object.

        t_vec3 :: model_scale !<
        !! Scale factor for the STL object.

        t_vec3 :: model_rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.

        integer :: model_spc !<
        !! Number of samples per cell to use when discretizing the STL object.

        real(wp) :: model_threshold !<
        !! Threshold to turn on smoothen STL patch.

    end type ic_patch_parameters

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

        t_vec3 :: model_translate !<
        !! Translation of the STL object.

        t_vec3 :: model_scale !<
        !! Scale factor for the STL object.

        t_vec3 :: model_rotate !<
        !! Angle to rotate the STL object along each cartesian coordinate axis,
        !! in radians.

        integer :: model_spc !<
        !! Number of samples per cell to use when discretizing the STL object.

        real(wp) :: model_threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ib_patch_parameters

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

    !> Derived type annexing the flow probe location
    type probe_parameters
        real(wp) :: x !< First coordinate location
        real(wp) :: y !< Second coordinate location
        real(wp) :: z !< Third coordinate location
    end type probe_parameters

    type mpi_io_airfoil_ib_var
        integer, dimension(2) :: view
        type(probe_parameters), allocatable, dimension(:) :: var
    end type mpi_io_airfoil_ib_var

    !> Derived type annexing integral regions
    type integral_parameters
        real(wp) :: xmin !< Min. boundary first coordinate direction
        real(wp) :: xmax !< Max. boundary first coordinate direction
        real(wp) :: ymin !< Min. boundary second coordinate direction
        real(wp) :: ymax !< Max. boundary second coordinate direction
        real(wp) :: zmin !< Min. boundary third coordinate direction
        real(wp) :: zmax !< Max. boundary third coordinate direction
    end type integral_parameters

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

    !> Acoustic source source_spatial pre-calculated values
    type source_spatial_type
        integer, dimension(:, :), allocatable :: coord !< List of grid points indices with non-zero source_spatial values
        real(wp), dimension(:), allocatable :: val !< List of non-zero source_spatial values
        real(wp), dimension(:), allocatable :: angle !< List of angles with x-axis for mom source term vector
        real(wp), dimension(:, :), allocatable :: xyz_to_r_ratios !< List of [xyz]/r for mom source term vector
    end type source_spatial_type

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

    !> Species parameters
    type species_parameters
        character(LEN=name_len) :: name !< Name of species
    end type species_parameters

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

end module m_derived_types
