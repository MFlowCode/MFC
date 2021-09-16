!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features the definitions of all the custom-defined
!!        derived types that are utilized throughout the post_process code.
MODULE m_derived_types
    
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: num_fluids_max = 10  !<
    !! Maximum number of fluids allowed
    
    
    !> Derived type adding the field position (fp) as an attribute
    TYPE field_position
        REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: fp !< Field position
    END TYPE field_position

    !> Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf !< Scalar field
    END TYPE scalar_field

    !> MPI pass-throughs
    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var
    
    !> Derived type attaching information about the beginning and the end bounds
    !! as attributes (beg - beginning)
    TYPE bounds_info
        INTEGER :: beg  !< Beginning
        INTEGER :: end  !< End
    END TYPE bounds_info
    
    !> Variables for bubble dynamic variables
    TYPE bub_bounds_info
        integer :: beg !< Beginning of bub. variables
        integer :: end !< End of bub. variables
        integer, dimension(:), allocatable :: rs !< Bubble radii
        integer, dimension(:), allocatable :: vs !< Bubble radial velocities
        integer, dimension(:), allocatable :: ps !< Bubble pressures
        integer, dimension(:), allocatable :: ms !< Bubble mass fluxes
    END TYPE bub_bounds_info   

    !> Derived type annexing the physical parameters (PP) of fluids. 
    TYPE physical_parameters
        REAL(KIND(0d0)) :: gamma !< Sp. heat ratio
        REAL(KIND(0d0)) :: pi_inf !< Liquid stiffness
        REAL(KIND(0d0)) :: mul0 !< Bubble viscosity
        REAL(KIND(0d0)) :: ss   !< Bubble surface tension 
        REAL(KIND(0d0)) :: pv   !< Bubble vapour pressure 
        REAL(KIND(0d0)) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: M_v  !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: mu_v !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: k_v  !< Bubble constants (see Preston (2007), Ando (2010))
    END TYPE physical_parameters
    
    
END MODULE m_derived_types
