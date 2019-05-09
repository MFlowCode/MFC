! MFC v3.0 - Pre-process Code: m_derived_types.f90
! Description: This file contains the definitions of all of the custom-defined
!              types used in the pre-process code.
! Author: Vedran Coralic
! Date: 06/09/12


MODULE m_derived_types
    
    
    IMPLICIT NONE
    
    
    ! Maximum number of patches allowed
    INTEGER, PARAMETER :: num_patches_max   = 10
    
    ! maximum number of fluids allowed
    integer, parameter :: num_fluids_max    = 10
     
    ! Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf => NULL()
    END TYPE scalar_field

    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var

    TYPE int_bounds_info
        INTEGER :: beg
        INTEGER :: end
    END TYPE int_bounds_info
    
    
    ! Derived type adding beginning (beg) and end bounds info as attributes
    TYPE bounds_info
        REAL(KIND(0d0)) :: beg
        REAL(KIND(0d0)) :: end
    END TYPE bounds_info
   
    TYPE bub_bounds_info
        integer :: beg
        integer :: end
        integer, dimension(:), allocatable :: rs
        integer, dimension(:), allocatable :: vs
        integer, dimension(:), allocatable :: ps
        integer, dimension(:), allocatable :: ms
    END TYPE bub_bounds_info
    
    ! Derived type adding initial condition (ic) patch parameters as attributes
    TYPE ic_patch_parameters
        
        ! Type of geometry for the patch
        INTEGER :: geometry
        
        ! Location of the geometric center, i.e. the centroid, of the patch. It
        ! is specified through its x-, y- and z-coordinates, respectively.
        REAL(KIND(0d0)) :: x_centroid, y_centroid, z_centroid
        
        ! Dimensions of the patch. These include, respectively, the lengths in
        ! the x-, y- and z-coordinate directions, as well as the radius.
        REAL(KIND(0d0)) :: length_x, length_y, length_z
        REAL(KIND(0d0)) :: radius

        ! Vector indicating the various radii for the elliptical and ellipsoidal
        ! patch geometries. It is specified through its x-, y-, and z-components
        ! respectively.
        REAL(KIND(0d0)), DIMENSION(3) :: radii
        
        ! The isentropic vortex parameters administrating, respectively, both
        ! the amplitude of the disturbance as well as its domain of influence.
        REAL(KIND(0d0)) :: epsilon, beta
        
        ! Normal vector indicating the orientation of the patch. It is specified
        ! through its x-, y- and z-components, respectively.
        REAL(KIND(0d0)), DIMENSION(3) :: normal
        
        ! List of permissions that indicate to the current patch which preceding
        ! patches it is allowed to overwrite when it is in process of being laid
        ! out in the domain
        LOGICAL, DIMENSION(0:num_patches_max-1) :: alter_patch
        
        ! Permission indicating to the current patch whether its boundaries will
        ! be smoothed out across a few cells or whether they are to remain sharp
        LOGICAL :: smoothen
        
        ! Identity (id) of the patch with which current patch is to get smoothed
        INTEGER :: smooth_patch_id
        
        ! Smoothing coefficient (coeff) adminstrating the size of the stencil of
        ! cells across which boundaries of the current patch will be smeared out
        REAL(KIND(0d0)) :: smooth_coeff
        
        ! Primitive variables associated with the patch. In order, these include
        ! the partial densities, density, velocity, pressure, volume fractions,
        ! specific heat ratio function and the liquid stiffness function.
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: alpha_rho
        REAL(KIND(0d0))                            :: rho
        REAL(KIND(0d0)), DIMENSION(3)              :: vel
        REAL(KIND(0d0))                            :: pres
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: alpha
        REAL(KIND(0d0))                            :: gamma
        REAL(KIND(0d0))                            :: pi_inf
        
        ! SHB: For bubbles
        REAL(KIND(0d0))    :: R0
        REAL(KIND(0d0))    :: V0
        
        ! NOTE: The requirements for the specification of the above parameters
        ! are strongly dependent on both the choice of the multicomponent flow
        ! model as well as the choice of the patch geometry.
        
    END TYPE ic_patch_parameters
    
    
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
