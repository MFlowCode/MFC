! MFC v3.0 - Post-process Code: m_derived_types.f90
! Description: This file contains the definitions of all of the custom-defined
!              types used in the post-process code.
! Author: Vedran Coralic
! Date: 06/08/12


MODULE m_derived_types
    
    
    IMPLICIT NONE
    
    
    ! Maximum number of fluids allowed
    INTEGER, PARAMETER :: num_fluids_max = 10
    
    
    ! Derived type adding the field position (fp) as an attribute
    TYPE field_position
        REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: fp
    END TYPE field_position

    ! Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf
    END TYPE scalar_field

    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var
    
    
    ! Derived type attaching information about the beginning and the end bounds
    ! as attributes (beg - beginning)
    TYPE bounds_info
        INTEGER :: beg
        INTEGER :: end
    END TYPE bounds_info
    
    TYPE bub_bounds_info
        integer :: beg
        integer :: end
        integer, dimension(:), allocatable :: rs
        integer, dimension(:), allocatable :: vs
        integer, dimension(:), allocatable :: ps
        integer, dimension(:), allocatable :: ms
    END TYPE bub_bounds_info
    
    ! Derived type annexing the physical parameters (PP) of the fluids. These
    ! include the specific heat ratio function and liquid stiffness function.
    TYPE physical_parameters
        REAL(KIND(0d0)) :: gamma
        REAL(KIND(0d0)) :: pi_inf
        REAL(KIND(0d0)) :: mul0
        REAL(KIND(0d0)) :: ss
        REAL(KIND(0d0)) :: pv
        REAL(KIND(0d0)) :: gamma_v
        REAL(KIND(0d0)) :: M_v
        REAL(KIND(0d0)) :: mu_v
        REAL(KIND(0d0)) :: k_v
    END TYPE physical_parameters
    
    
END MODULE m_derived_types
