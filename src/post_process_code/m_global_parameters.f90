!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains all of the parameters characterizing the
!!      computational domain, simulation algorithm, stiffened equation of
!!      state and finally, the formatted database file(s) structure.
MODULE m_global_parameters
    
    
    ! Dependencies =============================================================
    USE mpi                     !< Message passing interface (MPI) module
    
    USE m_derived_types         !< Definitions of the derived types
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    !> @name Logistics
    !> @{ 
    INTEGER                    :: num_procs            !< Number of processors
    INTEGER        , PARAMETER :: num_stcls_min = 5    !< Mininum # of stencils
    INTEGER        , PARAMETER :: path_len      = 400  !< Maximum path length
    INTEGER        , PARAMETER :: name_len      = 50   !< Maximum name length
    REAL(KIND(0d0)), PARAMETER :: dflt_real     = -1d6 !< Default real value
    INTEGER        , PARAMETER :: dflt_int      = -100 !< Default integer value
    REAL(KIND(0d0)), PARAMETER :: sgm_eps       = 1d-16 !< Segmentation tolerance
    CHARACTER(LEN = path_len)  :: case_dir             !< Case folder location
    !> @} 
    
    ! Computational Domain Parameters ==========================================
    
    INTEGER :: proc_rank !< Rank of the local processor
    
    !> @name Number of cells in the x-, y- and z-coordinate directions
    !> @{
    INTEGER :: m, m_root
    INTEGER :: n
    INTEGER :: p
    !> @}

    !> @name Cylindrical coordinates (either axisymmetric or full 3D)
    !> @{
    LOGICAL :: cyl_coord
    INTEGER :: grid_geometry
    !> @}

    !> @name Global number of cells in each direction
    !> @{
    INTEGER :: m_glb, n_glb, p_glb
    !> @}

    INTEGER :: num_dims !< Number of spatial dimensions

    !> @name Cell-boundary locations in the x-, y- and z-coordinate directions
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cb, x_root_cb, y_cb, z_cb
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: coarse_x_cb, coarse_y_cb, coarse_z_cb
    !> @}

    !> @name Cell-center locations in the x-, y- and z-coordinate directions
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cc, x_root_cc, y_cc, z_cc
    !> @}

    !> Cell-width distributions in the x-, y- and z-coordinate directions
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: dx, dy, dz
    !> @}

    INTEGER :: buff_size !<
    !! Number of cells in buffer region. For the variables which feature a buffer
    !! region, this region is used to store information outside the computational
    !! domain based on the boundary conditions.
    
    INTEGER :: t_step_start  !< First time-step directory
    INTEGER :: t_step_stop   !< Last time-step directory
    INTEGER :: t_step_save   !< Interval between consecutive time-step directory
    
    ! NOTE: The variables m_root, x_root_cb and x_root_cc contain the grid data
    ! of the defragmented computational domain. They are only used in 1D. For
    ! serial simulations, they are equal to m, x_cb and x_cc, respectively.
    
    ! ==========================================================================
    
    
    !> @name Simulation Algorithm Parameters
    !> @{
    INTEGER :: model_eqns      !< Multicomponent flow model
    INTEGER :: num_fluids      !< Number of different fluids present in the flow
    LOGICAL :: adv_alphan      !< Advection of the last volume fraction
    LOGICAL :: mpp_lim         !< Maximum volume fraction limiter
    INTEGER :: sys_size        !< Number of unknowns in the system of equations
    INTEGER :: weno_order      !< Order of accuracy for the WENO reconstruction
    LOGICAL :: mixture_err     !< Mixture error limiter
    LOGICAL :: alt_soundspeed  !< Alternate sound speed
    !> @}

    !> @name Annotations of the structure, i.e. the organization, of the state vectors
    !> @{
    TYPE(bounds_info) :: cont_idx                  !< Indexes of first & last continuity eqns.
    TYPE(bounds_info) :: mom_idx                   !< Indexes of first & last momentum eqns.
    INTEGER           :: E_idx                     !< Index of energy equation
    TYPE(bounds_info) :: adv_idx                   !< Indexes of first & last advection eqns.
    TYPE(bounds_info) :: internalEnergies_idx      !< Indexes of first & last internal energy eqns.
    TYPE(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    INTEGER           :: gamma_idx                 !< Index of specific heat ratio func. eqn.
    INTEGER           :: alf_idx                   !< Index of specific heat ratio func. eqn.
    INTEGER           :: pi_inf_idx                !< Index of liquid stiffness func. eqn.
    !> @}

    !> @name Boundary conditions in the x-, y- and z-coordinate directions
    !> @{
    TYPE(bounds_info) :: bc_x, bc_y, bc_z
    !> @}

    LOGICAL :: parallel_io    !< Format of the data files


    INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_coords !< 
    !! Processor coordinates in MPI_CART_COMM

    INTEGER, ALLOCATABLE, DIMENSION(:) :: start_idx !< 
    !! Starting cell-center index of local processor in global grid

    TYPE(mpi_io_var), PUBLIC :: MPI_IO_DATA

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    CHARACTER(LEN = name_len) :: mpiiofs
    INTEGER :: mpi_info_int
    !> @}

    INTEGER, PRIVATE :: ierr
    ! ==========================================================================
    
    
    TYPE(physical_parameters), DIMENSION(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.
    
    ! ==========================================================================
    
    
    ! Formatted Database File(s) Structure Parameters ==========================
    

    INTEGER :: format !< Format of the database file(s)

    LOGICAL :: coarsen_silo
    

    INTEGER :: precision !< Floating point precision of the database file(s)
    
    !> @name Size of the ghost zone layer in the x-, y- and z-coordinate directions.
    !! The definition of the ghost zone layers is only necessary when using the
    !! Silo database file format in multidimensions. These zones provide VisIt
    !! with the subdomain connectivity information that it requires in order to
    !! produce smooth plots.
    !> @{
    TYPE(bounds_info) :: offset_x, offset_y, offset_z
    !> @}

    !> @name The list of all possible flow variables that may be written to a database
    !! file. It includes partial densities, density, momentum, velocity, energy,
    !! pressure, volume fraction(s), specific heat ratio function, specific heat
    !! ratio, liquid stiffness function, liquid stiffness, primitive variables,
    !! conservative variables, speed of sound, the vorticity, 
    !! and the numerical Schlieren function.
    !> @{
    LOGICAL, DIMENSION(num_fluids_max) :: alpha_rho_wrt
    LOGICAL                            :: rho_wrt
    LOGICAL, DIMENSION(3)              :: mom_wrt
    LOGICAL, DIMENSION(3)              :: vel_wrt
    INTEGER                            :: flux_lim
    LOGICAL, DIMENSION(3)              :: flux_wrt
    LOGICAL, DIMENSION(num_fluids_max) :: kappa_wrt
    LOGICAL                            :: E_wrt
    LOGICAL                            :: pres_wrt
    LOGICAL, DIMENSION(num_fluids_max) :: alpha_wrt
    LOGICAL                            :: gamma_wrt
    LOGICAL                            :: heat_ratio_wrt
    LOGICAL                            :: pi_inf_wrt
    LOGICAL                            :: pres_inf_wrt
    LOGICAL                            :: prim_vars_wrt
    LOGICAL                            :: cons_vars_wrt
    LOGICAL                            :: c_wrt
    LOGICAL, DIMENSION(3)              :: omega_wrt
    LOGICAL                            :: schlieren_wrt
    !> @}

    !> @name Options for Fourier decomposition in the azimuthal direction if 3D
    !! cylindrical coordinates are used
    !> @{
    LOGICAL                            :: fourier_decomp
    TYPE(bounds_info)                  :: fourier_modes
    !> @}

    REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: schlieren_alpha    !<
    !! Amplitude coefficients of the numerical Schlieren function that are used
    !! to adjust the intensity of numerical Schlieren renderings for individual
    !! fluids. This enables waves and interfaces of varying strenghts and in all
    !! of the fluids to be made simulatenously visible on a single plot.
    
    INTEGER :: fd_order !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when vorticity and/or the numerical
    !! Schlieren function are to be outputted to the formatted database file(s).
    
    INTEGER :: fd_number !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.
    
    ! ==========================================================================
    
    !> @name Reference parameters for Tait EOS
    !> @{
    REAL(KIND(0d0)) :: rhoref, pref
    !> @}

    !> @name Bubble modeling variables and parameters
    !> @{
    INTEGER         :: nb
    REAL(KIND(0d0)) :: R0ref
    REAL(KIND(0d0)) :: Ca, Web, Re_inv
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: weight, R0, V0
    LOGICAL         :: bubbles
    LOGICAL         :: polytropic
    LOGICAL         :: polydisperse
    INTEGER         :: thermal  !< 1 = adiabatic, 2 = isotherm, 3 = transfer
    REAL(KIND(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T 
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN 
    REAL(KIND(0d0)) :: poly_sigma
    !> @}

    ! Mathematical and Physical Constants ======================================
    REAL(KIND(0d0)), PARAMETER :: pi = 3.141592653589793d0
    ! ==========================================================================


    CONTAINS
        
        
        
        
        !> Assigns default values to user inputs prior to reading
        !!      them in. This allows for an easier consistency check of
        !!      these parameters once they are read from the input file.        
        SUBROUTINE s_assign_default_values_to_user_inputs() ! ------------------

            INTEGER :: i !< Generic loop iterator
            
            
            ! Logistics
            case_dir = ' '
            
            
            ! Computational domain parameters
            m      = dflt_int
            m_root = dflt_int
            n      = dflt_int
            p      = dflt_int
            
            cyl_coord = .FALSE.
            
            t_step_start = dflt_int
            t_step_stop  = dflt_int
            t_step_save  = dflt_int
            
            
            ! Simulation algorithm parameters
            model_eqns = dflt_int
            num_fluids = dflt_int
            adv_alphan = .FALSE.
            weno_order = dflt_int
            mixture_err = .FALSE.
            alt_soundspeed = .FALSE.
            
            bc_x%beg = dflt_int
            bc_x%end = dflt_int
            bc_y%beg = dflt_int
            bc_y%end = dflt_int
            bc_z%beg = dflt_int
            bc_z%end = dflt_int
            
            
            ! Fluids physical parameters
            DO i = 1, num_fluids_max
                fluid_pp(i)%gamma  = dflt_real
                fluid_pp(i)%pi_inf = dflt_real
            END DO
            
            
            ! Formatted database file(s) structure parameters
            format = dflt_int
            
            precision = dflt_int
            
            coarsen_silo  = .FALSE.

            alpha_rho_wrt  = .FALSE.
            rho_wrt        = .FALSE.
            mom_wrt        = .FALSE.
            vel_wrt        = .FALSE.
            flux_lim  = dflt_int
            flux_wrt = .FALSE.
            parallel_io = .FALSE.
            kappa_wrt      = .FALSE.
            E_wrt          = .FALSE.
            pres_wrt       = .FALSE.
            alpha_wrt      = .FALSE.
            gamma_wrt      = .FALSE.
            heat_ratio_wrt = .FALSE.
            pi_inf_wrt     = .FALSE.
            pres_inf_wrt   = .FALSE.
            prim_vars_wrt  = .FALSE.
            cons_vars_wrt  = .FALSE.
            c_wrt          = .FALSE.
            omega_wrt      = .FALSE.
            schlieren_wrt  = .FALSE.
            
            schlieren_alpha = dflt_real

            fourier_decomp = .FALSE.
            fourier_modes%beg = dflt_int
            fourier_modes%end = dflt_int
            
            fd_order = dflt_int
 
            ! Tait EOS
            rhoref  = dflt_real
            pref    = dflt_real

            ! Bubble modeling
            bubbles = .FALSE.
            R0ref   = dflt_real 
            nb      = dflt_int
            polydisperse= .FALSE.
            poly_sigma = dflt_real
             
        END SUBROUTINE s_assign_default_values_to_user_inputs ! ----------------
        
        
        
        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module        
        SUBROUTINE s_initialize_global_parameters_module() ! ----------------------

            INTEGER :: i, fac
            
            ! Setting m_root equal to m in the case of a 1D serial simulation
            IF(num_procs == 1 .AND. n == 0) m_root = m
            
            
            ! Gamma/Pi_inf Model ===============================================
            IF(model_eqns == 1) THEN
                
                ! Setting number of fluids
                num_fluids = 1
                
                ! Annotating structure of the state and flux vectors belonging
                ! to the system of equations defined by the selected number of
                ! spatial dimensions and the gamma/pi_inf model
                cont_idx%beg = 1
                cont_idx%end = cont_idx%beg
                mom_idx%beg  = cont_idx%end + 1
                mom_idx%end  = cont_idx%end + num_dims
                E_idx        = mom_idx%end + 1
                adv_idx%beg  = E_idx + 1
                adv_idx%end  = adv_idx%beg + 1
                gamma_idx    = adv_idx%beg
                pi_inf_idx   = adv_idx%end
                sys_size     = adv_idx%end
                
            ! ==================================================================
                
                
            ! Volume Fraction Model (5-equation model) =========================
            ELSE IF(model_eqns == 2) THEN
                
                ! Annotating structure of the state and flux vectors belonging
                ! to the system of equations defined by the selected number of
                ! spatial dimensions and the volume fraction model
                cont_idx%beg = 1
                cont_idx%end = num_fluids
                mom_idx%beg  = cont_idx%end + 1
                mom_idx%end  = cont_idx%end + num_dims
                E_idx        = mom_idx%end + 1
                adv_idx%beg  = E_idx + 1
                adv_idx%end  = E_idx + num_fluids
                
                IF( (adv_alphan .NEQV. .TRUE.) .AND. &
                    (num_fluids > 1)) adv_idx%end = adv_idx%end - 1
                
                sys_size     = adv_idx%end
                
                IF (bubbles) THEN
                    alf_idx  = adv_idx%end
                ELSE
                    alf_idx = 0
                END IF
                
                IF (bubbles) THEN
                    bub_idx%beg = sys_size+1
                    bub_idx%end = sys_size+2*nb
                    IF (polytropic .NEQV. .TRUE.) THEN
                        bub_idx%end = sys_size+4*nb
                    END IF
                    sys_size = bub_idx%end

                    ALLOCATE( bub_idx%rs(nb), bub_idx%vs(nb) )
                    ALLOCATE( bub_idx%ps(nb), bub_idx%ms(nb) )
                    ALLOCATE( weight(nb),R0(nb),V0(nb) )

                    DO i = 1, nb
                        IF (polytropic .NEQV. .TRUE.) THEN
                            fac = 4
                        ELSE
                            fac = 2
                        END IF
                        
                        bub_idx%rs(i) = bub_idx%beg+(i-1)*fac
                        bub_idx%vs(i) = bub_idx%rs(i)+1

                        IF (polytropic .NEQV. .TRUE.) THEN
                            bub_idx%ps(i) = bub_idx%vs(i)+1
                            bub_idx%ms(i) = bub_idx%ps(i)+1
                        END IF
                    END DO

                    IF (nb == 1) THEN
                        weight(:)   = 1d0
                        R0(:)       = 1d0
                        V0(:)       = 0d0
                    ELSE IF (nb > 1) THEN
                        CALL s_simpson(nb)
                        V0(:)       = 0d0
                    ELSE
                        STOP 'Invalid value of nb'
                    END IF

                    IF (polytropic .NEQV. .TRUE.) THEN
                        CALL s_initialize_nonpoly
                    ELSE
                        rhoref  = 1.d0
                        pref    = 1.d0
                    END IF
                END IF          
                
            ! ==================================================================

            ! Volume Fraction Model (6-equation model) =========================
            ELSE IF(model_eqns == 3) THEN

                ! Annotating structure of the state and flux vectors belonging
                ! to the system of equations defined by the selected number of
                ! spatial dimensions and the volume fraction model
                cont_idx%beg = 1
                cont_idx%end = num_fluids
                mom_idx%beg  = cont_idx%end + 1
                mom_idx%end  = cont_idx%end + num_dims
                E_idx        = mom_idx%end + 1
                adv_idx%beg  = E_idx + 1
                adv_idx%end  = E_idx + num_fluids
                IF(adv_alphan .NEQV. .TRUE.) adv_idx%end = adv_idx%end - 1
                internalEnergies_idx%beg  = adv_idx%end + 1
                internalEnergies_idx%end  = adv_idx%end + num_fluids
                sys_size     = internalEnergies_idx%end

            ELSE IF(model_eqns == 4) THEN
                cont_idx%beg = 1 ! one continuity equation
                cont_idx%end = 1 !num_fluids
                mom_idx%beg  = cont_idx%end + 1 ! one momentum equation in each
                mom_idx%end  = cont_idx%end + num_dims
                E_idx        = mom_idx%end + 1 ! one energy equation
                adv_idx%beg  = E_idx + 1
                adv_idx%end  = adv_idx%beg !one volume advection equation
                alf_idx      = adv_idx%end
                sys_size     = alf_idx !adv_idx%end
                
                IF (bubbles) THEN
                    bub_idx%beg = sys_size+1
                    bub_idx%end = sys_size+2*nb
                    IF (polytropic .NEQV. .TRUE.) THEN
                        bub_idx%end = sys_size+4*nb
                    END IF
                    sys_size = bub_idx%end

                    ALLOCATE( bub_idx%rs(nb), bub_idx%vs(nb) )
                    ALLOCATE( bub_idx%ps(nb), bub_idx%ms(nb) )
                    ALLOCATE( weight(nb),R0(nb),V0(nb) )

                    DO i = 1, nb
                        IF (polytropic .NEQV. .TRUE.) THEN
                            fac = 4
                        ELSE
                            fac = 2
                        END IF
                        
                        bub_idx%rs(i) = bub_idx%beg+(i-1)*fac
                        bub_idx%vs(i) = bub_idx%rs(i)+1

                        IF (polytropic .NEQV. .TRUE.) THEN
                            bub_idx%ps(i) = bub_idx%vs(i)+1
                            bub_idx%ms(i) = bub_idx%ps(i)+1
                        END IF
                    END DO

                    IF (nb == 1) THEN
                        weight(:)   = 1d0
                        R0(:)       = 1d0
                        V0(:)       = 0d0
                    ELSE IF (nb > 1) THEN
                        CALL s_simpson(nb)
                        V0(:)       = 0d0
                    ELSE
                        STOP 'Invalid value of nb'
                    END IF

                    IF (polytropic .NEQV. .TRUE.) THEN
                        CALL s_initialize_nonpoly
                    ELSE
                        rhoref  = 1.d0
                        pref    = 1.d0
                    END IF
                END IF
             END IF
            ! ==================================================================
            
            ALLOCATE(MPI_IO_DATA%view(1:sys_size))
            ALLOCATE(MPI_IO_DATA%var(1:sys_size))

            DO i = 1, sys_size
                ALLOCATE(MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => NULL()
            END DO

            
            ! Size of the ghost zone layer is non-zero only when post-processing
            ! the raw simulation data of a parallel multidimensional computation
            ! in the Silo-HDF5 format. If this is the case, one must also verify
            ! whether the raw simulation data is 2D or 3D. In the 2D case, size
            ! of the z-coordinate direction ghost zone layer must be zeroed out.
            IF(num_procs == 1 .OR. format /= 1 .OR. n == 0) THEN
                
                offset_x%beg = 0
                offset_x%end = 0
                offset_y%beg = 0
                offset_y%end = 0
                offset_z%beg = 0
                offset_z%end = 0
                
            ELSEIF(p == 0) THEN
                
                offset_z%beg = 0
                offset_z%end = 0
                
            END IF
            
            
            ! Determining the finite-difference number and the buffer size. Note
            ! that the size of the buffer is unrelated to the order of the WENO
            ! scheme. Rather, it is directly dependent on maximum size of ghost
            ! zone layers and possibly the order of the finite difference scheme
            ! used for the computation of vorticity and/or numerical Schlieren
            ! function.
            buff_size = MAX( offset_x%beg, offset_x%end, offset_y%beg, &
                               offset_y%end, offset_z%beg, offset_z%end  )
            
            IF(ANY(omega_wrt) .OR. schlieren_wrt) THEN
                fd_number = MAX(1, fd_order/2)
                buff_size = buff_size + fd_number
            END IF
            
            
            ! Allocating the grid variables in the x-coordinate direction
            ALLOCATE(x_cb(-1-offset_x%beg : m + offset_x%end))
            ALLOCATE(x_cc(    -buff_size  : m + buff_size   ))
            ALLOCATE(  dx(    -buff_size  : m + buff_size   ))
            
            
            ! Allocating grid variables in the y- and z-coordinate directions
            IF(n > 0) THEN
                
                ALLOCATE(y_cb(-1-offset_y%beg : n + offset_y%end))
                ALLOCATE(y_cc(     -buff_size : n + buff_size   ))
                ALLOCATE(  dy(     -buff_size : n + buff_size   ))
                
                IF(p > 0) THEN
                    ALLOCATE(z_cb(-1-offset_z%beg : p + offset_z%end))
                    ALLOCATE(z_cc(    -buff_size  : p + buff_size   ))
                    ALLOCATE(  dz(    -buff_size  : p + buff_size   ))
                END IF
                
            ! Allocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
            ELSE
                
                ALLOCATE(x_root_cb(-1:m_root))
                ALLOCATE(x_root_cc( 0:m_root))
                
            END IF

            IF (coarsen_silo) THEN
                ALLOCATE(coarse_x_cb(-1-offset_x%beg : (m/2) + offset_x%end))
                IF (n > 0) THEN
                    ALLOCATE(coarse_y_cb(-1-offset_y%beg : (n/2) + offset_y%end))
                    IF (p > 0) ALLOCATE(coarse_z_cb(-1-offset_z%beg : (p/2) + offset_z%end))
                END IF
            END IF            
            
            IF (cyl_coord .NEQV. .TRUE.) THEN ! Cartesian grid
                grid_geometry = 1
            ELSEIF (cyl_coord .AND. p == 0) THEN ! Axisymmetric cylindrical grid
                grid_geometry = 2
            ELSE ! Fully 3D cylindrical grid
                grid_geometry = 3
            END IF

        END SUBROUTINE s_initialize_global_parameters_module ! --------------------

        !> Subroutine to initialize variable for non-polytropic gas modeling processes
        SUBROUTINE s_initialize_nonpoly

            INTEGER :: ir
            REAL(KIND(0.D0)) :: rhol0
            REAL(KIND(0.D0)) :: pl0
            REAL(KIND(0.D0)) :: uu
            REAL(KIND(0.D0)) :: D_m
            REAL(KIND(0.D0)) :: temp
            REAL(KIND(0.D0)) :: omega_ref
            REAL(KIND(0.D0)), DIMENSION(Nb) :: chi_vw0
            REAL(KIND(0.D0)), DIMENSION(Nb) :: cp_m0
            REAL(KIND(0.D0)), DIMENSION(Nb) :: k_m0
            REAL(KIND(0.D0)), DIMENSION(Nb) :: rho_m0
            REAL(KIND(0.D0)), DIMENSION(Nb) :: x_vw
            ! polytropic index used to compute isothermal natural frequency
            REAL(KIND(0.D0)), PARAMETER :: k_poly = 1.D0
            ! universal gas constant
            REAL(KIND(0.D0)), PARAMETER :: Ru = 8314.D0


            ! liquid physical properties
            REAL(KIND(0.D0)) :: mul0, ss, pv, gamma_v, M_v, mu_v

            ! gas physical properties
            REAL(KIND(0.D0)) :: gamma_m, gamma_n, M_n, mu_n

            rhol0 = rhoref
            pl0   = pref
            
            ALLOCATE( pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb) )
            ALLOCATE( k_n(nb), k_v(nb), omegaN(nb) )
            ALLOCATE( Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb) ) 
            
            pb0(:)      = dflt_real
            mass_n0(:)  = dflt_real
            mass_v0(:)  = dflt_real
            Pe_T(:)     = dflt_real
            omegaN(:)   = dflt_real

            mul0    = fluid_pp(1)%mul0
            ss      = fluid_pp(1)%ss
            pv      = fluid_pp(1)%pv
            gamma_v = fluid_pp(1)%gamma_v
            M_v     = fluid_pp(1)%M_v
            mu_v    = fluid_pp(1)%mu_v
            k_v(:)  = fluid_pp(1)%k_v

            gamma_n = fluid_pp(2)%gamma_v
            M_n     = fluid_pp(2)%M_v
            mu_n    = fluid_pp(2)%mu_v
            k_n(:)  = fluid_pp(2)%k_v

            gamma_m = gamma_n
            IF (thermal==2 ) gamma_m = 1.d0 !isothermal

            temp = 293.15D0
            D_m  = 0.242D-4
            uu = DSQRT( pl0/rhol0 )

            omega_ref = 3.D0*k_poly*Ca + 2.D0*( 3.D0*k_poly-1.D0 )/Web
           
            !!! thermal properties !!!
            ! gas constants
            R_n = Ru/M_n
            R_v = Ru/M_v
            ! phi_vn & phi_nv (phi_nn = phi_vv = 1)
            phi_vn = ( 1.D0+DSQRT(mu_v/mu_n)*(M_n/M_v)**(0.25D0) )**2 &
                / ( DSQRT(8.D0)*DSQRT(1.D0+M_v/M_n) )
            phi_nv = ( 1.D0+DSQRT(mu_n/mu_v)*(M_v/M_n)**(0.25D0) )**2 &
                / ( DSQRT(8.D0)*DSQRT(1.D0+M_n/M_v) )
            ! internal bubble pressure
            pb0 = pl0 + 2.D0*ss/( R0ref*R0 )            
            
            ! mass fraction of vapor
            chi_vw0 = 1.D0/( 1.D0+R_v/R_n*(pb0/pv-1.D0) )
            ! specific heat for gas/vapor mixture
            cp_m0 = chi_vw0*R_v*gamma_v/( gamma_v-1.D0 ) &
                  + ( 1.D0-chi_vw0 )*R_n*gamma_n/( gamma_n-1.D0 )
            ! mole fraction of vapor
            x_vw = M_n*chi_vw0/( M_v+(M_n-M_v)*chi_vw0 )
            ! thermal conductivity for gas/vapor mixture
            k_m0 = x_vw*k_v/( x_vw+(1.D0-x_vw)*phi_vn ) &
                 + ( 1.D0-x_vw )*k_n/( x_vw*phi_nv+1.D0-x_vw )
            ! mixture density
            rho_m0 = pv/( chi_vw0*R_v*temp )

            ! mass of gas/vapor computed using dimensional quantities
            mass_n0 = 4.D0*( pb0-pv )*pi/( 3.D0*R_n*temp*rhol0 )*R0**3
            mass_v0 = 4.D0*pv*pi/( 3.D0*R_v*temp*rhol0 )*R0**3
            ! Peclet numbers
            Pe_T = rho_m0*cp_m0*uu*R0ref/k_m0
            Pe_c = uu*R0ref/D_m
            ! nondimensional properties
            R_n = rhol0*R_n*temp/pl0
            R_v = rhol0*R_v*temp/pl0
            k_n = k_n/k_m0
            k_v = k_v/k_m0
            pb0 = pb0/pl0
            pv = pv/pl0

            ! bubble wall temperature, normalized by T0, in the liquid
            ! keeps a constant (cold liquid assumption)
            Tw = 1.D0
            ! natural frequencies
            omegaN = DSQRT( 3.D0*k_poly*Ca+2.D0*(3.D0*k_poly-1.D0)/(Web*R0) )/R0

            pl0 = 1.d0
            DO ir = 1,Nb
                CALL s_transcoeff( omegaN(ir)*R0(ir),Pe_T(ir)*R0(ir), &
                          Re_trans_T(ir),Im_trans_T(ir) )
                CALL s_transcoeff( omegaN(ir)*R0(ir),Pe_c*R0(ir), &
                          Re_trans_c(ir),Im_trans_c(ir) )
            END DO     
            Im_trans_T = 0d0
            Im_trans_c = 0d0

            rhoref = 1.d0
            pref = 1.d0
        END SUBROUTINE s_initialize_nonpoly
       
        !> Subroutine to compute the transfer coefficient for non-polytropic gas modeling
        SUBROUTINE s_transcoeff( omega,peclet,Re_trans,Im_trans )

            REAL(KIND(0.D0)), INTENT(IN) :: omega
            REAL(KIND(0.D0)), INTENT(IN) :: peclet
            REAL(KIND(0.D0)), INTENT(OUT) :: Re_trans
            REAL(KIND(0.D0)), INTENT(OUT) :: Im_trans
            COMPLEX :: trans, c1, c2, c3
            COMPLEX :: imag = ( 0.,1. )
            REAL(KIND(0.D0)) :: f_transcoeff

            c1 = imag*omega*peclet
            c2 = CSQRT( c1 )
            c3 = ( CEXP(c2)-CEXP(-c2) )/( CEXP(c2)+CEXP(-c2) ) ! TANH(c2)
            trans = ( (c2/c3-1.D0)**(-1)-3.D0/c1 )**( -1 ) ! transfer function

            Re_trans = DBLE( trans )
            Im_trans = AIMAG( trans )

        END SUBROUTINE s_transcoeff

        
        !> Subroutine to initialize parallel infrastructure
        SUBROUTINE s_initialize_parallel_io() ! --------------------------------

            num_dims = 1 + MIN(1,n) + MIN(1,p)

            ALLOCATE(proc_coords(1:num_dims))

            IF (parallel_io .NEQV. .TRUE.) RETURN

            ! Option for Lustre file system (Darter/Comet/Stampede)
            WRITE(mpiiofs, '(A)') '/lustre_'
            mpiiofs = TRIM(mpiiofs)
            CALL MPI_INFO_CREATE(mpi_info_int, ierr)
            CALL MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

            ! Option for UNIX file system (Hooke/Thomson)
            ! WRITE(mpiiofs, '(A)') '/ufs_'
            ! mpiiofs = TRIM(mpiiofs)
            ! mpi_info_int = MPI_INFO_NULL

            ALLOCATE(start_idx(1:num_dims))

        END SUBROUTINE s_initialize_parallel_io ! ------------------------------
        
        
        
        !> Deallocation procedures for the module
        SUBROUTINE s_finalize_global_parameters_module() ! -------------------
            
            INTEGER :: i
            
            ! Deallocating the grid variables for the x-coordinate direction
            DEALLOCATE(x_cb, x_cc, dx)
            
            
            ! Deallocating grid variables for the y- and z-coordinate directions
            IF(n > 0) THEN
                
                DEALLOCATE(y_cb, y_cc, dy)
                
                IF(p > 0) DEALLOCATE(z_cb, z_cc, dz)
                
            ! Deallocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
            ELSE
                
                DEALLOCATE(x_root_cb, x_root_cc)
                
            END IF

            IF (coarsen_silo) THEN
                DEALLOCATE(coarse_x_cb)
                IF (n > 0) THEN
                    DEALLOCATE(coarse_y_cb)
                    IF (p > 0) DEALLOCATE(coarse_z_cb)
                END IF
            END IF            
            
            DEALLOCATE(proc_coords)
            IF (parallel_io) THEN
                DEALLOCATE(start_idx)
                DO i = 1, sys_size
                    MPI_IO_DATA%var(i)%sf => NULL()
                END DO

                DEALLOCATE(MPI_IO_DATA%var)
                DEALLOCATE(MPI_IO_DATA%view)
            END IF
            
        END SUBROUTINE s_finalize_global_parameters_module ! -----------------
        
        
        !> Computes the bubble number density n from the conservative variables
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
        SUBROUTINE s_comp_n_from_cons( vftmp,nRtmp,ntmp )

            REAL(KIND(0.D0)), INTENT(IN) :: vftmp
            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: nRtmp
            REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
            REAL(KIND(0.D0)) :: nR3

            CALL s_quad( nRtmp**3,nR3 )  !returns itself if NR0 = 1
            ntmp = DSQRT( (4.d0*pi/3.d0)*nR3/vftmp )

        END SUBROUTINE s_comp_n_from_cons


        !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
        SUBROUTINE s_comp_n_from_prim( vftmp,Rtmp,ntmp )
            
            REAL(KIND(0.D0)), INTENT(IN) :: vftmp
            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: Rtmp
            REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
            REAL(KIND(0.D0)) :: R3

            CALL s_quad( Rtmp**3,R3 )  !returns itself if NR0 = 1
            ntmp = (3.d0/(4.d0*pi)) * vftmp/R3

        END SUBROUTINE s_comp_n_from_prim


        !> Computes the quadrature for polydisperse bubble populations
        !! @param func is the bubble dynamic variables for each bin
        !! @param mom is the computed moment
        SUBROUTINE s_quad( func,mom )

            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: func
            REAL(KIND(0.D0)), INTENT(OUT) :: mom

            mom = DOT_PRODUCT( weight,func )

        END SUBROUTINE s_quad


        !> Computes the Simpson weights for quadrature
        !! @param Npt is the number of bins that represent the polydisperse bubble population
        SUBROUTINE s_simpson( Npt )

            INTEGER, INTENT(IN) :: Npt
            INTEGER :: ir
            REAL(KIND(0.D0)) :: R0mn
            REAL(KIND(0.D0)) :: R0mx
            REAL(KIND(0.D0)) :: dphi
            REAL(KIND(0.D0)) :: tmp
            REAL(KIND(0.D0)) :: sd
            REAL(KIND(0.D0)), DIMENSION(Npt) :: phi

            ! nondiml. min. & max. initial radii for numerical quadrature
            !sd   = 0.05D0
            !R0mn = 0.75D0
            !R0mx = 1.3D0

            sd = poly_sigma
            R0mn = 0.8D0*DEXP(-2.8D0 * sd)
            R0mx = 0.2D0*DEXP( 9.5D0 * sd) + 1.D0

            ! phi = ln( R0 ) & return R0
            DO ir = 1,Npt
                phi(ir) = DLOG( R0mn ) &
                    + DBLE( ir-1 )*DLOG( R0mx/R0mn )/DBLE( Npt-1 )
                R0(ir) = DEXP( phi(ir) )
            END DO
            dphi = phi(2) - phi(1)

            ! weights for quadrature using Simpson's rule
            DO ir = 2,Npt-1
                ! Gaussian
                tmp = DEXP( -0.5D0*(phi(ir)/sd)**2 )/DSQRT( 2.D0*pi )/sd
                IF ( MOD(ir,2)==0 ) THEN
                    weight(ir) = tmp*4.D0*dphi/3.D0
                ELSE
                    weight(ir) = tmp*2.D0*dphi/3.D0
                END IF
            END DO
            
            tmp = DEXP( -0.5D0*(phi(1)/sd)**2 )/DSQRT( 2.D0*pi )/sd
            weight(1) = tmp*dphi/3.D0
            tmp = DEXP( -0.5D0*(phi(Npt)/sd)**2 )/DSQRT( 2.D0*pi )/sd
            weight(Npt) = tmp*dphi/3.D0

        END SUBROUTINE s_simpson
        
END MODULE m_global_parameters
