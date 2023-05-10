!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types

!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
module m_derived_types

    use m_constants !< Constants
    use m_precision_select
    
    implicit none

    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(wp), allocatable, dimension(:, :, :) :: fp !< Field position
    end type field_position

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(wp), pointer, dimension(:, :, :) :: sf => null()
    end type scalar_field

    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf !< Vector field
    end type vector_field

    !> Integer boounds for variables
    type int_bounds_info
        integer :: beg
        integer :: end
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
        !! The isentropic vortex parameters administrating, respectively, both
        !! the amplitude of the disturbance as well as its domain of influence.

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
        !! Smoothing coefficient (coeff) adminstrating the size of the stencil of
        !! cells across which boundaries of the current patch will be smeared out

        real(wp), dimension(num_fluids_max) :: alpha_rho
        real(wp) :: rho
        real(wp), dimension(3) :: vel
        real(wp) :: pres
        real(wp), dimension(num_fluids_max) :: alpha
        real(wp) :: gamma
        real(wp) :: pi_inf !<


        !! Primitive variables associated with the patch. In order, these include
        !! the partial densities, density, velocity, pressure, volume fractions,
        !! specific heat ratio function and the liquid stiffness function.

        real(wp), dimension(6) :: tau_e
        !! Elastic stresses added to primitive variables if hypoelasticity = True

        real(wp) :: R0 !< Bubble size
        real(wp) :: V0 !< Bubble velocity

        real(wp) :: p0 !< Bubble size
        real(wp) :: m0 !< Bubble velocity

    end type ic_patch_parameters

    !> Derived type annexing the physical parameters (PP) of the fluids. These
    !! include the specific heat ratio function and liquid stiffness function.
    type physical_parameters
        real(wp) :: gamma   !< Sp. heat ratio
        real(wp) :: pi_inf  !< Liquid stiffness
        real(wp), dimension(2) :: Re      !< Reynolds number
        real(wp) :: mul0    !< Bubble viscosity
        real(wp) :: ss      !< Bubble surface tension
        real(wp) :: pv      !< Bubble vapour pressure
        real(wp) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: M_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: mu_v    !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: k_v     !< Bubble constants (see Preston (2007), Ando (2010))
        real(wp) :: G
    end type physical_parameters

    !> Derived type annexing the flow probe location
    type probe_parameters
        real(wp) :: x !< First coordinate location
        real(wp) :: y !< Second coordinate location
        real(wp) :: z !< Third coordinate location
    end type probe_parameters

    !> Derived type annexing integral regions
    type integral_parameters
        real(wp) :: xmin !< Min. boundary first coordinate direction
        real(wp) :: xmax !< Max. boundary first coordinate direction
        real(wp) :: ymin !< Min. boundary second coordinate direction
        real(wp) :: ymax !< Max. boundary second coordinate direction
        real(wp) :: zmin !< Min. boundary third coordinate direction
        real(wp) :: zmax !< Max. boundary third coordinate direction
    end type integral_parameters

    !> Monopole acoustic source parameters
    type mono_parameters
        real(wp), dimension(3) :: loc !< Physical location of acoustic source
        real(wp) :: mag !< Magnitude
        real(wp) :: length !< Length of line source
        real(wp) :: npulse !< Number of cycles of pulse
        real(wp) :: dir !< Direction of pulse
        real(wp) :: delay !< Time-delay of pulse start
        integer :: pulse
        integer :: support
        real(wp) :: aperture
        real(wp) :: foc_length
    end type mono_parameters

end module m_derived_types
