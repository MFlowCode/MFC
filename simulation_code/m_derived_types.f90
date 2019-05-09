! MFC v3.0 - Simulation Code: m_derived_types.f90
! Description: This module features the definitions of all the custom-defined
!              derived types that are utilized throughout the simulation code.
! Author: Vedran Coralic
! Date: 06/08/12


MODULE m_derived_types
    
    
    IMPLICIT NONE
    
    
    ! Maximum number of fluids in the simulation
    INTEGER, PARAMETER :: num_fluids_max = 10
    
    ! Maximum number of flow probes in the simulation
    INTEGER, PARAMETER :: num_probes_max = 10

    ! Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf => NULL()
    END TYPE scalar_field
    
    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var
    
    ! Derived type annexing a vector field (VF)
    TYPE vector_field
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: vf
    END TYPE vector_field
    
    
    ! Derived type annexing beginning and ending bounds information
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
    
    ! Derived type annexing the physical parameters (PP) of fluids. They include
    ! stiffened equation of state parameters, the Reynolds (Re) numbers, and the
    ! Weber (We) numbers, respectively.
    TYPE physical_parameters
        REAL(KIND(0d0))                            :: gamma
        REAL(KIND(0d0))                            :: pi_inf
        REAL(KIND(0d0)), DIMENSION(2)              :: Re
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: We
        REAL(KIND(0d0)) :: mul0
        REAL(KIND(0d0)) :: ss
        REAL(KIND(0d0)) :: pv
        REAL(KIND(0d0)) :: gamma_v
        REAL(KIND(0d0)) :: M_v
        REAL(KIND(0d0)) :: mu_v
        REAL(KIND(0d0)) :: k_v
    END TYPE physical_parameters
    
    
    ! Derived type annexing the flow probe location
    TYPE probe_parameters
        REAL(KIND(0d0)) :: x
        REAL(KIND(0d0)) :: y
        REAL(KIND(0d0)) :: z
    END TYPE probe_parameters

    TYPE integral_parameters
        REAL(KIND(0d0)) :: xmin
        REAL(KIND(0d0)) :: xmax
        REAL(KIND(0d0)) :: ymin
        REAL(KIND(0d0)) :: ymax
        REAL(KIND(0d0)) :: zmin
        REAL(KIND(0d0)) :: zmax
    END TYPE integral_parameters




    TYPE mono_parameters
        REAL(KIND(0d0)),dimension(3) :: loc
        REAL(KIND(0d0)) :: mag
        REAL(KIND(0d0)) :: length
        REAL(KIND(0d0)) :: npulse
        REAL(KIND(0d0)) :: dir
        integer :: pulse
        integer :: support
    END TYPE mono_parameters

END MODULE m_derived_types
