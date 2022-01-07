!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types

!> @brief This module features the definitions of all the custom-defined
!!        derived types that are utilized throughout the simulation code.
MODULE m_derived_types
    
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: num_fluids_max = 10 !<
    !! Maximum number of fluids in the simulation
    
    INTEGER, PARAMETER :: num_probes_max = 10 !<
    !! Maximum number of flow probes in the simulation

    !> Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf => NULL() !< Scalar field
    END TYPE scalar_field
    
    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var
    
    !> Derived type annexing a vector field (VF)
    TYPE vector_field
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: vf !< Vector field
    END TYPE vector_field
    
    
    !> Derived type annexing beginning and ending bounds information
    TYPE bounds_info
        INTEGER :: beg !< Begin boundary
        INTEGER :: end !< End boundary
    END TYPE bounds_info
    
    !> Variables for bubble dynamic variables
    TYPE bub_bounds_info
        INTEGER :: beg !< Beginning of bub. variables
        INTEGER :: end !< End of bub. variables
        INTEGER, DIMENSION(:), ALLOCATABLE :: rs !< Bubble radii
        INTEGER, DIMENSION(:), ALLOCATABLE :: vs !< Bubble radial velocities
        INTEGER, DIMENSION(:), ALLOCATABLE :: ps !< Bubble pressures
        INTEGER, DIMENSION(:), ALLOCATABLE :: ms !< Bubble mass fluxes

        INTEGER, DIMENSION(:,:), ALLOCATABLE :: moms !< Moment indices for qbmm
    END TYPE bub_bounds_info    
    

    !> Derived type annexing the physical parameters (PP) of fluids. 
    TYPE physical_parameters
        REAL(KIND(0d0))                            :: gamma  !< Sp. heat ratio
        REAL(KIND(0d0))                            :: pi_inf !< Liquid stiffness
        REAL(KIND(0d0)), DIMENSION(2)              :: Re    !< Reynolds number
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: We    !< Weber number
        REAL(KIND(0d0)) :: mul0 !< Bubble viscosity
        REAL(KIND(0d0)) :: ss   !< Bubble surface tension 
        REAL(KIND(0d0)) :: pv   !< Bubble vapour pressure 
        REAL(KIND(0d0)) :: gamma_v !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: M_v  !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: mu_v !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: k_v  !< Bubble constants (see Preston (2007), Ando (2010))
        REAL(KIND(0d0)) :: G    !< Shear Modulus
    END TYPE physical_parameters
    
    
    !> Derived type annexing the flow probe location
    TYPE probe_parameters
        REAL(KIND(0d0)) :: x !< First coordinate location
        REAL(KIND(0d0)) :: y !< Second coordinate location
        REAL(KIND(0d0)) :: z !< Third coordinate location
    END TYPE probe_parameters

    !> Derived type annexing integral regions
    TYPE integral_parameters
        REAL(KIND(0d0)) :: xmin !< Min. boundary first coordinate direction
        REAL(KIND(0d0)) :: xmax !< Max. boundary first coordinate direction
        REAL(KIND(0d0)) :: ymin !< Min. boundary second coordinate direction
        REAL(KIND(0d0)) :: ymax !< Max. boundary second coordinate direction
        REAL(KIND(0d0)) :: zmin !< Min. boundary third coordinate direction
        REAL(KIND(0d0)) :: zmax !< Max. boundary third coordinate direction
    END TYPE integral_parameters


    !> Monopole acoustic source parameters
    TYPE mono_parameters
        REAL(KIND(0d0)), DIMENSION(3) :: loc !< Physical location of acoustic source
        REAL(KIND(0d0)) :: mag !< Magnitude
        REAL(KIND(0d0)) :: length !< Length of line source
        REAL(KIND(0d0)) :: npulse !< Number of cycles of pulse
        REAL(KIND(0d0)) :: dir !< Direction of pulse
        REAL(KIND(0d0)) :: delay !< Time-delay of pulse start
        INTEGER :: pulse 
        INTEGER :: support
    END TYPE mono_parameters

END MODULE m_derived_types
