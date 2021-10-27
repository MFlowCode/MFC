!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features the definitions of all the custom-defined
!!        derived types that are utilized throughout the post_process code.
module m_derived_types

    implicit none

    integer, parameter :: num_fluids_max = 10  !<
    !! Maximum number of fluids allowed

    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(kind(0d0)), allocatable, dimension(:, :, :) :: fp !< Field position
    end type field_position

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(kind(0d0)), pointer, dimension(:, :, :) :: sf !< Scalar field
    end type scalar_field

    !> MPI pass-throughs
    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    !> Derived type attaching information about the beginning and the end bounds
    !! as attributes (beg - beginning)
    type bounds_info
        integer :: beg  !< Beginning
        integer :: end  !< End
    end type bounds_info

    !> Variables for bubble dynamic variables
    type bub_bounds_info
        integer :: beg !< Beginning of bub. variables
        integer :: end !< End of bub. variables
        integer, dimension(:), allocatable :: rs !< Bubble radii
        integer, dimension(:), allocatable :: vs !< Bubble radial velocities
        integer, dimension(:), allocatable :: ps !< Bubble pressures
        integer, dimension(:), allocatable :: ms !< Bubble mass fluxes
    end type bub_bounds_info

    !> Derived type annexing the physical parameters (PP) of fluids.
    type physical_parameters
        real(kind(0d0)) :: gamma !< Sp. heat ratio
        real(kind(0d0)) :: pi_inf !< Liquid stiffness
        real(kind(0d0)) :: mul0 !< Bubble viscosity
        real(kind(0d0)) :: ss   !< Bubble surface tension
        real(kind(0d0)) :: pv   !< Bubble vapour pressure
        real(kind(0d0)) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: M_v  !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: mu_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: k_v  !< Bubble constants (see Preston (2007), Ando (2010))
    end type physical_parameters

end module m_derived_types
