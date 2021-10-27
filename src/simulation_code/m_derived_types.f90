!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features the definitions of all the custom-defined
!!        derived types that are utilized throughout the simulation code.
module m_derived_types

    implicit none

    integer, parameter :: num_fluids_max = 10 !<
    !! Maximum number of fluids in the simulation

    integer, parameter :: num_probes_max = 10 !<
    !! Maximum number of flow probes in the simulation

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(kind(0d0)), pointer, dimension(:, :, :) :: sf => null() !< Scalar field
    end type scalar_field

    type mpi_io_var
        integer, allocatable, dimension(:) :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf !< Vector field
    end type vector_field

    !> Derived type annexing beginning and ending bounds information
    type bounds_info
        integer :: beg !< Begin boundary
        integer :: end !< End boundary
    end type bounds_info

    !> Variables for bubble dynamic variables
    type bub_bounds_info
        integer :: beg !< Beginning of bub. variables
        integer :: end !< End of bub. variables
        integer, dimension(:), allocatable :: rs !< Bubble radii
        integer, dimension(:), allocatable :: vs !< Bubble radial velocities
        integer, dimension(:), allocatable :: ps !< Bubble pressures
        integer, dimension(:), allocatable :: ms !< Bubble mass fluxes

        integer, dimension(:, :), allocatable :: moms !< Moment indices for qbmm
    end type bub_bounds_info

    !> Derived type annexing the physical parameters (PP) of fluids.
    type physical_parameters
        real(kind(0d0))                            :: gamma  !< Sp. heat ratio
        real(kind(0d0))                            :: pi_inf !< Liquid stiffness
        real(kind(0d0)), dimension(2)              :: Re    !< Reynolds number
        real(kind(0d0)) :: mul0 !< Bubble viscosity
        real(kind(0d0)) :: ss   !< Bubble surface tension
        real(kind(0d0)) :: pv   !< Bubble vapour pressure
        real(kind(0d0)) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: M_v  !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: mu_v !< Bubble constants (see Preston (2007), Ando (2010))
        real(kind(0d0)) :: k_v  !< Bubble constants (see Preston (2007), Ando (2010))
    end type physical_parameters

    !> Derived type annexing the flow probe location
    type probe_parameters
        real(kind(0d0)) :: x !< First coordinate location
        real(kind(0d0)) :: y !< Second coordinate location
        real(kind(0d0)) :: z !< Third coordinate location
    end type probe_parameters

    !> Derived type annexing integral regions
    type integral_parameters
        real(kind(0d0)) :: xmin !< Min. boundary first coordinate direction
        real(kind(0d0)) :: xmax !< Max. boundary first coordinate direction
        real(kind(0d0)) :: ymin !< Min. boundary second coordinate direction
        real(kind(0d0)) :: ymax !< Max. boundary second coordinate direction
        real(kind(0d0)) :: zmin !< Min. boundary third coordinate direction
        real(kind(0d0)) :: zmax !< Max. boundary third coordinate direction
    end type integral_parameters

    !> Monopole acoustic source parameters
    type mono_parameters
        real(kind(0d0)), dimension(3) :: loc !< Physical location of acoustic source
        real(kind(0d0)) :: mag !< Magnitude
        real(kind(0d0)) :: length !< Length of line source
        real(kind(0d0)) :: npulse !< Number of cycles of pulse
        real(kind(0d0)) :: dir !< Direction of pulse
        real(kind(0d0)) :: delay !< Time-delay of pulse start
        integer :: pulse
        integer :: support
    end type mono_parameters

end module m_derived_types
