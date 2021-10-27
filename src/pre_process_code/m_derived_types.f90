!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
module m_derived_types

    implicit none

    integer, parameter :: num_patches_max = 10 !<
    !! Maximum number of patches allowed

    integer, parameter :: num_fluids_max = 10 !<
    !! Maximum number of fluids allowed

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(kind(0d0)), pointer, dimension(:, :, :) :: sf => null()
    end type scalar_field

    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    !> Integer boounds for variables
    type int_bounds_info
        integer :: beg
        integer :: end
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

        real(kind(0d0)) :: epsilon, beta !<
        !! The isentropic vortex parameters administrating, respectively, both
        !! the amplitude of the disturbance as well as its domain of influence.

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
        real(kind(0d0))                            :: rho
        real(kind(0d0)), dimension(3)              :: vel
        real(kind(0d0))                            :: pres
        real(kind(0d0)), dimension(num_fluids_max) :: alpha
        real(kind(0d0))                            :: gamma
        real(kind(0d0))                            :: pi_inf !<
        !! Primitive variables associated with the patch. In order, these include
        !! the partial densities, density, velocity, pressure, volume fractions,
        !! specific heat ratio function and the liquid stiffness function.

        real(kind(0d0))    :: R0 !< Bubble size
        real(kind(0d0))    :: V0 !< Bubble velocity

        real(kind(0d0))    :: p0 !< Bubble size
        real(kind(0d0))    :: m0 !< Bubble velocity

    end type ic_patch_parameters

    !> Derived type annexing the physical parameters (PP) of the fluids. These
    !! include the specific heat ratio function and liquid stiffness function.
    type physical_parameters
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: mul0
        real(kind(0d0)) :: ss
        real(kind(0d0)) :: pv
        real(kind(0d0)) :: gamma_v
        real(kind(0d0)) :: M_v
        real(kind(0d0)) :: mu_v
        real(kind(0d0)) :: k_v
    end type physical_parameters

end module m_derived_types
