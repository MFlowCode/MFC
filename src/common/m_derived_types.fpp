!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types

#:include "macros.fpp"

!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
module m_derived_types

    use m_constants !< Constants

    implicit none

    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(kind(0d0)), allocatable, dimension(:, :, :) :: fp !< Field position
    end type field_position

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(kind(0d0)), pointer, dimension(:, :, :) :: sf => null()
    end type scalar_field

    !> Derived type for bubble variables pb and mv at quadrature nodes (qbmm)
    type pres_field
        real(kind(0d0)), pointer, dimension(:, :, :, :, :) :: sf => null()
    end type pres_field

    !> Derived type annexing an integer scalar field (SF)
    type integer_field
        integer, pointer, dimension(:, :, :) :: sf => null()
    end type integer_field

    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    type mpi_io_ib_var
        integer :: view
        type(integer_field) :: var
    end type mpi_io_ib_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf !< Vector field
    end type vector_field

    !> Integer bounds for variables
    type int_bounds_info
        integer :: beg
        integer :: end
        real(kind(0d0)) :: vb1
        real(kind(0d0)) :: vb2
        real(kind(0d0)) :: vb3
        real(kind(0d0)) :: ve1
        real(kind(0d0)) :: ve2
        real(kind(0d0)) :: ve3
    end type int_bounds_info

    !> Derived type adding beginning (beg) and end bounds info as attributes
    type bounds_info
        real(kind(0d0)) :: beg
        real(kind(0d0)) :: end
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
    type :: ic_model_parameters
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

        real(kind(0d0)) :: threshold !<
        !! Threshold to turn on smoothen STL patch.
    end type ic_model_parameters

    type :: t_triangle
        real(kind(0d0)), dimension(1:3, 1:3) :: v ! Vertices of the triangle
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

        real(kind(0d0)) :: x_centroid, y_centroid, z_centroid !<
        !! Location of the geometric center, i.e. the centroid, of the patch. It
        !! is specified through its x-, y- and z-coordinates, respectively.

        real(kind(0d0)) :: length_x, length_y, length_z !< Dimensions of the patch. x,y,z Lengths.
        real(kind(0d0)) :: radius !< Dimensions of the patch. radius.

        real(kind(0d0)), dimension(3) :: radii !<
        !! Vector indicating the various radii for the elliptical and ellipsoidal
        !! patch geometries. It is specified through its x-, y-, and z-components
        !! respectively.

        type(ic_model_parameters) :: model !< Model parameters

        real(kind(0d0)) :: epsilon, beta !<
        !! The spherical harmonics eccentricity parameters.

        real(kind(0d0)), dimension(3) :: normal !<
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

        real(kind(0d0)) :: smooth_coeff !<
        !! Smoothing coefficient (coeff) adminstrating the size of the stencil of
        !! cells across which boundaries of the current patch will be smeared out

        real(kind(0d0)), dimension(num_fluids_max) :: alpha_rho
        real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(3) :: vel
        real(kind(0d0)) :: pres
        real(kind(0d0)), dimension(num_fluids_max) :: alpha
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf !<
        real(kind(0d0)) :: cv !<
        real(kind(0d0)) :: qv !<
        real(kind(0d0)) :: qvp !<

        !! Primitive variables associated with the patch. In order, these include
        !! the partial densities, density, velocity, pressure, volume fractions,
        !! specific heat ratio function and the liquid stiffness function.

        real(kind(0d0)), dimension(6) :: tau_e
        !! Elastic stresses added to primitive variables if hypoelasticity = True

        real(kind(0d0)) :: R0 !< Bubble size
        real(kind(0d0)) :: V0 !< Bubble velocity

        real(kind(0d0)) :: p0 !< Bubble size
        real(kind(0d0)) :: m0 !< Bubble velocity

        integer :: hcid
        !! id for hard coded initial condition

        real(kind(0d0)) :: cf_val !! color function value

    end type ic_patch_parameters

    type ib_patch_parameters

        integer :: geometry !< Type of geometry for the patch

        real(kind(0d0)) :: x_centroid, y_centroid, z_centroid !<
        !! Location of the geometric center, i.e. the centroid, of the patch. It
        !! is specified through its x-, y- and z-coordinates, respectively.

        real(kind(0d0)) :: c, p, t, m

        real(kind(0d0)) :: length_x, length_y, length_z !< Dimensions of the patch. x,y,z Lengths.
        real(kind(0d0)) :: radius !< Dimensions of the patch. radius.
        real(kind(0d0)) :: theta

        logical :: slip

    end type ib_patch_parameters

    !> Derived type annexing the physical parameters (PP) of the fluids. These
    !! include the specific heat ratio function and liquid stiffness function.
    type physical_parameters
        real(kind(0d0)) :: gamma   !< Sp. heat ratio
        real(kind(0d0)) :: pi_inf  !< Liquid stiffness
        real(kind(0d0)), dimension(2) :: Re      !< Reynolds number
        real(kind(0d0)) :: cv      !< heat capacity
        real(kind(0d0)) :: qv      !< reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
        real(kind(0d0)) :: qvp     !< reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
        real(kind(0d0)) :: mul0    !< Bubble viscosity
        real(kind(0d0)) :: ss      !< Bubble surface tension
        real(kind(0d0)) :: pv      !< Bubble vapour pressure
        real(kind(0d0)) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: M_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: mu_v    !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: k_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: G
    end type physical_parameters

    !> Derived type annexing the flow probe location
    type probe_parameters
        real(kind(0d0)) :: x !< First coordinate location
        real(kind(0d0)) :: y !< Second coordinate location
        real(kind(0d0)) :: z !< Third coordinate location
    end type probe_parameters

    type mpi_io_airfoil_ib_var
        integer, dimension(2) :: view
        type(probe_parameters), allocatable, dimension(:) :: var
    end type mpi_io_airfoil_ib_var

    !> Derived type annexing integral regions
    type integral_parameters
        real(kind(0d0)) :: xmin !< Min. boundary first coordinate direction
        real(kind(0d0)) :: xmax !< Max. boundary first coordinate direction
        real(kind(0d0)) :: ymin !< Min. boundary second coordinate direction
        real(kind(0d0)) :: ymax !< Max. boundary second coordinate direction
        real(kind(0d0)) :: zmin !< Min. boundary third coordinate direction
        real(kind(0d0)) :: zmax !< Max. boundary third coordinate direction
    end type integral_parameters

    !> Acoustic source parameters
    type acoustic_parameters
        integer :: pulse !< Type of pulse
        integer :: support !< Type of support
        logical :: dipole !< Whether the source is a dipole or monopole
        real(kind(0d0)), dimension(3) :: loc !< Physical location of acoustic source
        real(kind(0d0)) :: mag !< Acoustic pulse magnitude
        real(kind(0d0)) :: length !< Length of planar source (2D/3D)
        real(kind(0d0)) :: height !< Height of planar source (3D)
        real(kind(0d0)) :: wavelength !< Wave length of pulse
        real(kind(0d0)) :: frequency !< Frequency of pulse
        real(kind(0d0)) :: gauss_sigma_dist !< sigma of Gaussian pulse multiplied by speed of sound
        real(kind(0d0)) :: gauss_sigma_time !< sigma of Gaussian pulse
        real(kind(0d0)) :: npulse !< Number of cycles of pulse
        real(kind(0d0)) :: dir !< Direction of pulse
        real(kind(0d0)) :: delay !< Time-delay of pulse start
        real(kind(0d0)) :: foc_length ! < Focal length of transducer
        real(kind(0d0)) :: aperture ! < Aperture diameter of transducer
        real(kind(0d0)) :: element_spacing_angle !< Spacing between aperture elements in 2D acoustic array
        real(kind(0d0)) :: element_polygon_ratio !< Ratio of aperture element diameter to side length of polygon connecting their centers, in 3D acoustic array
        real(kind(0d0)) :: rotate_angle !< Angle of rotation of the entire circular 3D acoustic array
        integer :: num_elements !< Number of elements in the acoustic array
        integer :: element_on !< Element in the acoustic array to turn on
    end type acoustic_parameters

    !> Acoustic source source_spatial pre-calculated values
    type source_spatial_type
        integer, dimension(:, :), allocatable :: coord !< List of grid points indices with non-zero source_spatial values
        real(kind(0d0)), dimension(:), allocatable :: val !< List of non-zero source_spatial values
        real(kind(0d0)), dimension(:), allocatable :: angle !< List of angles with x-axis for mom source term vector
        real(kind(0d0)), dimension(:, :), allocatable :: xyz_to_r_ratios !< List of [xyz]/r for mom source term vector
    end type source_spatial_type

    !> Ghost Point for Immersed Boundaries
    type ghost_point

        real(kind(0d0)), dimension(3) :: loc !< Physical location of the ghost point
        real(kind(0d0)), dimension(3) :: ip_loc !< Physical location of the image point
        integer, dimension(3) :: ip_grid !< Top left grid point of IP
        real(kind(0d0)), dimension(2, 2, 2) :: interp_coeffs !< Interpolation Coefficients of image point
        integer :: ib_patch_id !< ID of the IB Patch the ghost point is part of
        logical :: slip
        integer, dimension(3) :: DB

    end type ghost_point

end module m_derived_types
