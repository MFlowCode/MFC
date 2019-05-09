! MFC v3.0 - Simulation Code: m_global_parameters.f90
! Description: The module contains all of the parameters describing the program
!              logistics, the computational domain and the simulation algorithm.
!              Additionally, for the volume fraction model, physical parameters
!              of each of the fluids present in the flow are located here. They
!              include stiffened gas equation of state parameters, the Reynolds
!              numbers and the Weber numbers.
! Author: Vedran Coralic
! Date: 06/27/12


MODULE m_global_parameters
    
    
    ! Dependencies =============================================================
    USE mpi                    ! Message passing interface (MPI) module
    
    USE m_derived_types        ! Definitions of the derived types
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    ! Logistics ================================================================
    INTEGER                    :: num_procs             ! Number of processors
    INTEGER        , PARAMETER :: num_stcls_min = 5     ! Mininum # of stencils
    INTEGER        , PARAMETER :: path_len      = 400   ! Maximum path length
    INTEGER        , PARAMETER :: name_len      = 50    ! Maximum name length
    CHARACTER      , PARAMETER :: dflt_char     = ' '   ! Default string value
    REAL(KIND(0d0)), PARAMETER :: dflt_real     = -1d6  ! Default real value
    INTEGER        , PARAMETER :: dflt_int      = -100  ! Default integer value
    REAL(KIND(0d0)), PARAMETER :: sgm_eps       = 1d-16 ! Segmentation tolerance
    INTEGER        , PARAMETER :: fourier_rings = 5     ! Fourier filter ring limit
    CHARACTER(LEN = path_len)  :: case_dir              ! Case folder location
    LOGICAL                    :: run_time_info         ! Run-time output flag
    INTEGER                    :: t_step_old            ! Existing IC/grid folder
    ! ==========================================================================
    
    
    ! Computational Domain Parameters ==========================================
    
    ! Rank of the local processor
    INTEGER :: proc_rank
    
    ! Number of cells in the x-, y- and z-directions, respectively
    INTEGER :: m, n, p
    
    ! Global number of cells in each direction
    INTEGER :: m_glb, n_glb, p_glb

    ! Cylindrical coordinates (either axisymmetric or full 3D)
    LOGICAL :: cyl_coord
    INTEGER :: grid_geometry

    ! Cell-boundary (CB) locations in the x-, y- and z-directions, respectively
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:) :: x_cb, y_cb, z_cb
    
    ! Cell-center (CC) locations in the x-, y- and z-directions, respectively
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:) :: x_cc, y_cc, z_cc
    
    ! Cell-width distributions in the x-, y- and z-directions, respectively
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:) :: dx, dy, dz
    
    ! Size of the time-step
    REAL(KIND(0d0)) :: dt
    
    ! Starting time-step iteration, stopping time-step iteration and the number
    ! of time-step iterations between successive solution backups, respectively
    INTEGER :: t_step_start, t_step_stop, t_step_save
    
    ! ==========================================================================
    
    
    ! Simulation Algorithm Parameters ==========================================
    INTEGER         :: model_eqns     ! Multicomponent flow model
    INTEGER         :: num_dims       ! Number of spatial dimensions
    INTEGER         :: num_fluids     ! Number of fluids in the flow
    LOGICAL         :: adv_alphan     ! Advection of the last volume fraction
    LOGICAL         :: mpp_lim        ! Mixture physical parameters (MPP) limits
    INTEGER         :: time_stepper   ! Time-stepper algorithm
    INTEGER         :: weno_vars      ! WENO-reconstructed state variables type
    INTEGER         :: weno_order     ! Order of the WENO reconstruction
    INTEGER         :: weno_polyn     ! Degree of the WENO polynomials (polyn)
    REAL(KIND(0d0)) :: weno_eps       ! Binding for the WENO nonlinear weights
    LOGICAL         :: char_decomp    ! Characteristic decomposition
    LOGICAL         :: mapped_weno    ! WENO with mapping of nonlinear weights
    LOGICAL         :: mp_weno        ! Monotonicity preserving (MP) WENO
    LOGICAL         :: weno_avg       ! Average left/right cell-boundary states
    LOGICAL         :: weno_Re_flux   ! WENO reconstruct velocity gradients for viscous stress tensor
    INTEGER         :: riemann_solver ! Riemann solver algorithm
    INTEGER         :: wave_speeds    ! Wave speeds estimation method
    INTEGER         :: avg_state      ! Average state evaluation method
    LOGICAL         :: commute_err    ! Commutative error correction
    LOGICAL         :: split_err      ! Dimensional splitting error correction
    LOGICAL         :: alt_crv        ! Alternate curvature definition
    LOGICAL         :: alt_soundspeed ! Alternate mixture sound speed
    LOGICAL         :: regularization ! Regularization terms of Tiwari (2013)
    REAL(KIND(0d0)) :: reg_eps        ! User-defined interface thickness parameter for regularization terms
    LOGICAL         :: null_weights   ! Null undesired WENO weights
    LOGICAL         :: mixture_err    ! Mixture properties correction
    LOGICAL         :: tvd_riemann_flux ! Apply TVD flux limiter to left and right states inside Riemann solver
    LOGICAL         :: tvd_rhs_flux   ! Apply TVD flux limiter to to intercell fluxes outside Riemann solver
    LOGICAL         :: tvd_wave_speeds! Use TVD wavespeeds when computing fluxes inside Riemann solver
    INTEGER         :: flux_lim       ! Choice of flux limiter
    LOGICAL         :: We_riemann_flux ! Account for capillary effects in the Riemann solver
    LOGICAL         :: We_rhs_flux    ! Account for capillary effects using the conservative formulation in RHS
    LOGICAL         :: We_src         ! Account for capillary effects in non-conservative formulation in RHS
    LOGICAL         :: We_wave_speeds ! Account for capillary effects when computing the contact wave speed
    LOGICAL         :: lsq_deriv      ! Use linear least squares to calculate normals and curvatures

    LOGICAL         :: anti_diffusion ! Use anti-diffusion step of Shyue
    
    INTEGER         :: cpu_start, cpu_end, cpu_rate
    
    ! Boundary conditions (BC) in the x-, y- and z-directions, respectively
    TYPE(bounds_info) :: bc_x, bc_y, bc_z

    TYPE(bounds_info) :: bc_x_glb, bc_y_glb, bc_z_glb

    ! Format of the data files
    LOGICAL :: parallel_io

    ! Processor coordinates in MPI_CART_COMM
    INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_coords
    ! Starting cell-center index of local processor in global grid
    INTEGER, ALLOCATABLE, DIMENSION(:) :: start_idx

    TYPE(mpi_io_var), PUBLIC :: MPI_IO_DATA

    ! MPI info for parallel IO with Lustre file systems
    CHARACTER(LEN = name_len) :: mpiiofs
    INTEGER :: mpi_info_int

    INTEGER, PRIVATE :: ierr
    
    ! Annotations of the structure of the state and flux vectors in terms of the
    ! size and the configuration of the system of equations to which they belong
    INTEGER           :: sys_size                  ! Number of unknowns in system of eqns.
    TYPE(bounds_info) :: cont_idx                  ! Indexes of first & last continuity eqns.
    TYPE(bounds_info) :: mom_idx                   ! Indexes of first & last momentum eqns.
    INTEGER           :: E_idx                     ! Index of energy equation
    TYPE(bounds_info) :: adv_idx                   ! Indexes of first & last advection eqns.
    TYPE(bounds_info) :: internalEnergies_idx      ! Indexes of first & last internal energy eqns.
    TYPE(bub_bounds_info) :: bub_idx                    ! Indexes of first & last bubble variable eqns.
    INTEGER               :: alf_idx                    ! Index of void fraction
    INTEGER           :: gamma_idx                 ! Index of specific heat ratio func. eqn.
    INTEGER           :: pi_inf_idx                ! Index of liquid stiffness func. eqn.
    
    ! The number of fluids, along with their identifying indexes, respectively,
    ! for which viscous effects, e.g. the shear and/or the volume Reynolds (Re)
    ! numbers, will be non-negligible.
    INTEGER,              DIMENSION(2)   :: Re_size
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Re_idx
    
    ! The number of idiosyncratic material interfaces, along with the indexes of
    ! the fluids comprising them, respectively, for which the effects of surface
    ! tension, e.g. the Weber (We) numbers, will be non-negligible.
    INTEGER                              :: We_size
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: We_idx
    
    ! The number of fluids, along with their identifying indexes, respectively,
    ! for which the physical and the geometric curvatures (crv) of the material
    ! interfaces will be calculated.
    INTEGER                            :: crv_size
    INTEGER, ALLOCATABLE, DIMENSION(:) :: crv_idx
    
    ! The coordinate direction indexes and flags (flg), respectively, for which
    ! the configurations will be determined with respect to a working direction
    ! and that will be used to isolate the contributions, in that direction, in
    ! the dimensionally split system of equations.
    INTEGER        , DIMENSION(3) :: dir_idx
    REAL(KIND(0d0)), DIMENSION(3) :: dir_flg
    
    ! The WENO average (WA) flag regulates whether the calculation of any cell-
    ! average spatial derivatives is carried out in each cell by utilizing the
    ! arithmetic mean of the left and right, WENO-reconstructed, cell-boundary
    ! values or simply, the unaltered left and right, WENO-reconstructed, cell-
    ! boundary values.
    REAL(KIND(0d0)) :: wa_flg
    
    ! The number of cells that are necessary to be able to store enough boundary
    ! conditions data to march the solution in the physical computational domain
    ! to the next time-step.
    INTEGER :: buff_size
    ! END: Simulation Algorithm Parameters =====================================
    
    
    ! Fluids Physical Parameters ===============================================
    
    ! Database of the physical parameters of each of the fluids that is present
    ! in the flow. These include the stiffened gas equation of state parameters,
    ! the Reynolds numbers and the Weber numbers.
    TYPE(physical_parameters), DIMENSION(num_fluids_max) :: fluid_pp
    
    ! ==========================================================================
    
    
    ! The order of the finite-difference (fd) approximations of the first-order
    ! derivatives that need to be evaluated when the CoM or flow probe data
    ! files are to be written at each time step
    INTEGER :: fd_order

    ! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    ! it is a measure of the half-size of the finite-difference stencil for the
    ! selected order of accuracy.
    INTEGER :: fd_number

    LOGICAL, DIMENSION(num_fluids_max) :: com_wrt, cb_wrt
    LOGICAL :: probe_wrt
    LOGICAL :: integral_wrt
    INTEGER :: num_probes
    INTEGER :: num_integrals
    TYPE(probe_parameters), DIMENSION(num_probes_max) :: probe
    TYPE(integral_parameters), DIMENSION(num_probes_max) :: integral
    REAL(KIND(0d0)), DIMENSION(5) :: threshold_mf
    INTEGER, DIMENSION(5) :: moment_order

    ! SHB: Reference parameters for Tait EOS
    REAL(KIND(0d0)) :: rhoref, pref
    ! SHB: For bubble modeling
    INTEGER         :: nb
    REAL(KIND(0d0)) :: R0ref
    REAL(KIND(0d0)) :: Ca, Web, Re_inv
    REAL(KIND(0d0)), dimension(:), allocatable :: weight, R0, V0
    LOGICAL         :: bubbles
    LOGICAL         :: bubble_switch
    INTEGER         :: bubble_model 


    ! SHB: For non-polytropic modeling
    LOGICAL         :: polytropic
    INTEGER         :: thermal  !1 = adiabatic, 2 = isotherm, 3 = transfer
    REAL(kind(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v
    REAL(kind(0d0)), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T 
    REAL(kind(0d0)), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN 
    ! liquid physical properties
    REAL(KIND(0.D0)) :: mul0, ss, gamma_v, mu_v
    ! gas physical properties
    REAL(KIND(0.D0)) :: gamma_m, gamma_n, mu_n


    ! SHB: For acoustic monopole
    LOGICAL         :: monopole
    TYPE(mono_parameters), dimension(num_probes_max) :: mono
    INTEGER         :: num_mono

    !weights for switches
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: a_coeff, b_coeff
    
    !Source for bubble transfer
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: s_term

    ! SHB: hold onto p_tilde
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: ptil 
    ! ======================================================================

    ! Mathematical and Physical Constants ======================================
    REAL(KIND(0d0)), PARAMETER :: pi = 3.141592653589793d0
    ! ==========================================================================
  
    
    CONTAINS
        
        SUBROUTINE s_assign_default_values_to_user_inputs() ! ------------------
        ! Description: Assigns default values to the user inputs before reading
        !              them in. This enables for an easier consistency check of
        !              these parameters once they are read from the input file.
            
            
            ! Generic loop iterator
            INTEGER :: i,j
            
            
            ! Logistics
            case_dir      = dflt_char
            run_time_info = .FALSE.
            t_step_old = dflt_int
            
            
            ! Computational domain parameters
            m = dflt_int; n = dflt_int; p = dflt_int
            
            cyl_coord = .FALSE.

            dt = dflt_real
            
            t_step_start = dflt_int
            t_step_stop  = dflt_int
            t_step_save  = dflt_int
            
            
            ! Simulation algorithm parameters
            model_eqns       = dflt_int
            num_fluids       = dflt_int
            adv_alphan       = .FALSE.
            mpp_lim          = .FALSE.
            time_stepper     = dflt_int
            weno_vars        = dflt_int
            weno_order       = dflt_int
            weno_eps         = dflt_real
            char_decomp      = .FALSE.
            mapped_weno      = .FALSE.
            mp_weno          = .FALSE.
            weno_avg         = .FALSE.
            weno_Re_flux     = .FALSE.
            riemann_solver   = dflt_int
            wave_speeds      = dflt_int
            avg_state        = dflt_int
            commute_err      = .FALSE.
            split_err        = .FALSE.
            alt_crv          = .FALSE.
            alt_soundspeed   = .FALSE.
            regularization   = .FALSE.
            reg_eps          = dflt_real
            null_weights     = .FALSE.
            mixture_err      = .FALSE.
            tvd_riemann_flux = .FALSE.
            tvd_rhs_flux     = .FALSE.
            tvd_wave_speeds  = .FALSE.
            flux_lim         = dflt_int
            We_riemann_flux  = .FALSE.
            We_rhs_flux      = .FALSE.
            We_src           = .FALSE.
            We_wave_speeds   = .FALSE.
            lsq_deriv        = .FALSE.
            parallel_io      = .FALSE.
            anti_diffusion   = .FALSE.
            
            bc_x%beg = dflt_int; bc_x%end = dflt_int
            bc_y%beg = dflt_int; bc_y%end = dflt_int
            bc_z%beg = dflt_int; bc_z%end = dflt_int
            
            
            
            ! Fluids physical parameters
            DO i = 1, num_fluids_max
                fluid_pp(i)%gamma   = dflt_real
                fluid_pp(i)%pi_inf  = dflt_real
                fluid_pp(i)%Re(:)   = dflt_real
                fluid_pp(i)%We(:)   = dflt_real
                fluid_pp(i)%mul0    = dflt_real
                fluid_pp(i)%ss      = dflt_real
                fluid_pp(i)%pv      = dflt_real
                fluid_pp(i)%gamma_v = dflt_real
                fluid_pp(i)%M_v     = dflt_real
                fluid_pp(i)%mu_v    = dflt_real
                fluid_pp(i)%k_v     = dflt_real
            END DO

            !SHB Tait EOS
            rhoref  = dflt_real
            pref    = dflt_real

            !SHB Bubble modeling
            bubbles     = .FALSE.
            bubble_switch = .FALSE.
            bubble_model  = 1
            polytropic  = .TRUE.
            thermal     = dflt_int
            R0ref       = dflt_real
            nb          = dflt_int
            
            Ca      = dflt_real
            Re_inv  = dflt_real
            Web     = dflt_real
            
            !SHB monopole source
            monopole = .FALSE.
            num_mono = 1

            do j = 1,num_probes_max
                do i = 1,3
                    mono(j)%loc(i) = dflt_real
                end do
                mono(j)%mag    = dflt_real
                mono(j)%length = dflt_real
                mono(j)%dir    = dflt_real
                mono(j)%npulse = 1.d0
                mono(j)%pulse = 1
                mono(j)%support = 1
            end do

            fd_order = dflt_int
            com_wrt = .FALSE.
            cb_wrt = .FALSE.
            probe_wrt = .FALSE.
            integral_wrt = .FALSE.
            num_probes = dflt_int
            num_integrals = dflt_int

            DO i = 1, num_probes_max
                probe(i)%x = dflt_real
                probe(i)%y = dflt_real
                probe(i)%z = dflt_real
            END DO

            DO i = 1, num_probes_max
                integral(i)%xmin = dflt_real
                integral(i)%xmax = dflt_real
                integral(i)%ymin = dflt_real
                integral(i)%ymax = dflt_real
                integral(i)%ymin = dflt_real
                integral(i)%ymax = dflt_real
            END DO

            DO i = 1, 5
                threshold_mf(i) = dflt_real
                moment_order(i) = dflt_int
            END DO
            
        END SUBROUTINE s_assign_default_values_to_user_inputs ! ----------------
        
        
        
        SUBROUTINE s_initialize_global_parameters_module() ! -------------------
        ! Description: The computation of parameters, the allocation of memory,
        !              the association of pointers and/or the execution of any
        !              other procedures that are necessary to setup the module.
            
            
            ! Temporary indexes storage
            INTEGER :: tmp_idx
            
            ! Generic loop iterators
            INTEGER :: i,j
            
            ! Generic counter
            INTEGER :: k

            INTEGER :: fac
            
            TYPE(bounds_info) :: ix,iy,iz
            
            
            ! Determining the degree of the WENO polynomials
            weno_polyn = (weno_order - 1) / 2
            
            
            ! Initializing the number of fluids for which viscous effects will
            ! be non-negligible, the number of distinctive material interfaces
            ! for which surface tension will be important and also, the number
            ! of fluids for which the physical and geometric curvatures of the
            ! interfaces will be computed
            Re_size = 0; We_size = 0; crv_size = 0
            
            
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
                
                
            ! Volume Fraction Model ============================================
            ELSE
                
                ! Annotating structure of the state and flux vectors belonging
                ! to the system of equations defined by the selected number of
                ! spatial dimensions and the volume fraction model
                IF(model_eqns == 2) THEN
                    cont_idx%beg = 1
                    cont_idx%end = num_fluids
                    mom_idx%beg  = cont_idx%end + 1
                    mom_idx%end  = cont_idx%end + num_dims
                    E_idx        = mom_idx%end + 1
                    adv_idx%beg  = E_idx + 1
                    adv_idx%end  = E_idx + num_fluids                    
                    
                    IF( (adv_alphan .NEQV. .TRUE.) .and. &
                        (num_fluids > 1)) adv_idx%end = adv_idx%end - 1
 
                    sys_size = adv_idx%end

                    IF (bubbles) THEN
                        alf_idx  = adv_idx%end
                    ELSE
                        alf_idx = 0
                    END IF
                    
                    if (bubbles) then
                        bub_idx%beg = sys_size+1
                        bub_idx%end = sys_size+2*nb
                        if (polytropic .neqv. .TRUE.) then
                            bub_idx%end = sys_size+4*nb
                        end if
                        sys_size = bub_idx%end

                        allocate( bub_idx%rs(nb), bub_idx%vs(nb) )
                        allocate( bub_idx%ps(nb), bub_idx%ms(nb) )
                        allocate( weight(nb),R0(nb),V0(nb) )

                        do i = 1, nb
                            if (polytropic .neqv. .TRUE.) then
                                fac = 4
                            else
                                fac = 2
                            end if

                            bub_idx%rs(i) = bub_idx%beg+(i-1)*fac
                            bub_idx%vs(i) = bub_idx%rs(i)+1
                            
                            if (polytropic .neqv. .TRUE.) then
                                bub_idx%ps(i) = bub_idx%vs(i)+1
                                bub_idx%ms(i) = bub_idx%ps(i)+1
                            end if
                        end do

                        if (nb == 1) then
                            weight(:)   = 1d0
                            R0(:)       = 1d0
                            V0(:)       = 0d0
                        else if (nb > 1) then
                            call s_simpson(nb)
                            V0(:)       = 0d0
                        else
                            stop 'Invalid value of nb'
                        end if

                        if (polytropic .neqv. .TRUE.) then
                            call s_initialize_nonpoly
                        else
                            rhoref  = 1.d0
                            pref    = 1.d0
                        end if
                    end if
                ELSE IF(model_eqns == 3) THEN
                    cont_idx%beg = 1
                    cont_idx%end = num_fluids
                    mom_idx%beg  = cont_idx%end + 1
                    mom_idx%end  = cont_idx%end + num_dims
                    E_idx        = mom_idx%end + 1
                    adv_idx%beg  = E_idx + 1
                    adv_idx%end  = E_idx + num_fluids
                    internalEnergies_idx%beg  = adv_idx%end + 1
                    internalEnergies_idx%end  = adv_idx%end + num_fluids
                    sys_size     = internalEnergies_idx%end
                ELSE IF(model_eqns == 4) THEN
                    !SHB: 4 equation model with subgrid bubbles
                    !only works with one fluid
                    cont_idx%beg = 1 ! one continuity equation
                    cont_idx%end = 1 !num_fluids
                    mom_idx%beg  = cont_idx%end + 1 ! one momentum equation in each direction
                    mom_idx%end  = cont_idx%end + num_dims
                    E_idx        = mom_idx%end + 1 ! one energy equation
                    adv_idx%beg  = E_idx + 1
                    adv_idx%end  = adv_idx%beg !one volume advection equation
                    alf_idx      = adv_idx%end
                    sys_size     = alf_idx !adv_idx%end
                    
                    if (bubbles) then
                        bub_idx%beg = sys_size+1
                        bub_idx%end = sys_size+2*nb
                        if (polytropic .neqv. .TRUE.) then
                            bub_idx%end = sys_size+4*nb
                        end if
                        sys_size = bub_idx%end


                        allocate( bub_idx%rs(nb), bub_idx%vs(nb) )
                        allocate( bub_idx%ps(nb), bub_idx%ms(nb) )
                        allocate( weight(nb),R0(nb),V0(nb) )

                        do i = 1, nb
                            if (polytropic .neqv. .TRUE.) then
                                fac = 4
                            else
                                fac = 2
                            end if

                            bub_idx%rs(i) = bub_idx%beg+(i-1)*fac
                            bub_idx%vs(i) = bub_idx%rs(i)+1

                            if (polytropic .neqv. .TRUE.) then
                                bub_idx%ps(i) = bub_idx%vs(i)+1
                                bub_idx%ms(i) = bub_idx%ps(i)+1
                            end if
                        end do

                        if (nb == 1) then
                            weight(:)   = 1d0
                            R0(:)       = 1d0
                            V0(:)       = 0d0
                        else if (nb > 1) then
                            call s_simpson(nb)
                            V0(:)       = 0d0
                        else
                            stop 'Invalid value of nb'
                        end if

                        if (polytropic .neqv. .TRUE.) then
                            call s_initialize_nonpoly
                        else
                            rhoref  = 1.d0
                            pref    = 1.d0
                        end if
                    end if 
                END IF
            
                ! Determining the number of fluids for which the shear and the
                ! volume Reynolds numbers, e.g. viscous effects, are important
                DO i = 1, num_fluids
                    IF(fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
                    IF(fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
                END DO
                
                ! Determining the number of distinct material interfaces for
                ! which the Weber number, e.g. the capillarity, is important
                DO i = 1, num_fluids-1
                    DO j = i+1, num_fluids
                        IF(fluid_pp(i)%We(j) > 0) We_size = We_size + 1
                    END DO
                END DO
                
                ! Constraining the utilization of WENO averaging so that it may
                ! only be applied when the presence of viscous and/or capillary
                ! effects requires that the cell-average spatial derivatives be
                ! determined
                IF(   proc_rank == 0  .AND.   weno_avg   &
                                      .AND.              &
                    ALL(Re_size == 0) .AND. We_size == 0 ) THEN
                    PRINT '(A)', 'Unsupported combination of values '  // &
                                 'of weno_avg, fluid_pp(:)%Re(:) and ' // &
                                 'fluid_pp(:)%We(:). Exiting ...'
                    STOP
                END IF
                
                ! Bookkeeping the indexes of any viscous fluids and any pairs of
                ! fluids whose interface will support effects of surface tension
                IF(ANY(Re_size > 0)) THEN
                    
                    ALLOCATE(Re_idx(1:2,1:MAXVAL(Re_size)))
                    
                    k = 0
                    DO i = 1, num_fluids
                        IF(fluid_pp(i)%Re(1) > 0) THEN
                            k = k + 1; Re_idx(1,k) = i
                        END IF
                    END DO
                    
                    k = 0
                    DO i = 1, num_fluids
                        IF(fluid_pp(i)%Re(2) > 0) THEN
                            k = k + 1; Re_idx(2,k) = i
                        END IF
                    END DO
                    
                END IF
                
                IF(We_size > 0) THEN
                    
                    ALLOCATE(We_idx(1:We_size,1:2))
                    
                    k = 0
                    DO i = 1, num_fluids-1
                        DO j = i+1, num_fluids
                            IF(fluid_pp(i)%We(j) > 0) THEN
                                k = k + 1; We_idx(k,:) = (/i,j/)
                            END IF
                        END DO
                    END DO
                    
                    crv_size = 2
                    DO i = 2, We_size
                        IF(We_idx(i,1) > We_idx(i-1,1)) THEN
                            crv_size = crv_size + 1
                        END IF
                        IF(ALL(We_idx(i,2) /= We_idx(1:i-1,2))) THEN
                            crv_size = crv_size + 1
                        END IF
                    END DO
                    
                    ALLOCATE(crv_idx(1:crv_size))
                    
                    k = 2; crv_idx(1:2) = We_idx(1,:)
                    DO i = 2, We_size
                        IF(We_idx(i,1) > We_idx(i-1,1)) THEN
                            k = k + 1; crv_idx(k) = We_idx(i,1)
                        END IF
                        IF(ALL(We_idx(i,2) /= We_idx(1:i-1,2))) THEN
                            k = k + 1; crv_idx(k) = We_idx(i,2)
                        END IF
                    END DO
                    
                    DO j = 1, crv_size
                        DO i = 1, crv_size-1
                            IF(crv_idx(i) > crv_idx(i+1)) THEN
                                tmp_idx      = crv_idx(i+1)
                                crv_idx(i+1) = crv_idx( i )
                                crv_idx( i ) = tmp_idx
                            END IF
                        END DO
                    END DO
                    
                END IF
                
            END IF
            ! END: Volume Fraction Model =======================================
            
            ALLOCATE(MPI_IO_DATA%view(1:sys_size))
            ALLOCATE(MPI_IO_DATA%var(1:sys_size))

            DO i = 1, sys_size
                ALLOCATE(MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => NULL()
            END DO
            
            ! Configuring the WENO average flag that will be used to regulate
            ! whether any spatial derivatives are to computed in each cell by
            ! using the arithmetic mean of left and right, WENO-reconstructed,
            ! cell-boundary values or otherwise, the unaltered left and right,
            ! WENO-reconstructed, cell-boundary values
            wa_flg = 0d0; IF(weno_avg) wa_flg = 1d0
            
            
            ! Determining the number of cells that are needed in order to store
            ! sufficient boundary conditions data as to iterate the solution in
            ! the physical computational domain from one time-step iteration to
            ! the next one
            IF(commute_err .AND. ANY(Re_size > 0) .AND. We_size > 0) THEN
                buff_size = 4*weno_polyn + 3
            ELSEIF(We_size > 0) THEN
                buff_size = 3*weno_polyn + 3
            ELSEIF(commute_err .AND. ANY(Re_size > 0)) THEN
                buff_size = 3*weno_polyn + 2
            ELSEIF(ANY(Re_size > 0)) THEN
                buff_size = 2*weno_polyn + 2
            ELSEIF(commute_err) THEN
                buff_size = 2*weno_polyn + 1
            ELSE
                buff_size =   weno_polyn + 2
            END IF
            
             ! Configuring Coordinate Direction Indexes =========================
            if (bubbles) then
                ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
                IF(n > 0) iy%beg = -buff_size; IF(p > 0) iz%beg = -buff_size
                ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
                allocate(   a_coeff(    ix%beg:ix%end,   &
                                        iy%beg:iy%end,   &
                                        iz%beg:iz%end ), &
                            b_coeff(    ix%beg:ix%end,   &
                                        iy%beg:iy%end,   &
                                        iz%beg:iz%end ) )
                allocate(   s_term(     ix%beg:ix%end,   &
                                        iy%beg:iy%end,   &
                                        iz%beg:iz%end ) )
                allocate(   ptil   (    ix%beg:ix%end, &
                                        iy%beg:iy%end, &
                                        iz%beg:iz%end ) )
            end if           
 
            IF (probe_wrt) THEN
                fd_number = MAX(1, fd_order/2)
                buff_size = buff_size + fd_number
            END IF

            IF (cyl_coord .NEQV. .TRUE.) THEN ! Cartesian grid
                grid_geometry = 1
            ELSEIF (cyl_coord .AND. p == 0) THEN ! Axisymmetric cylindrical grid
                grid_geometry = 2
            ELSE ! Fully 3D cylindrical grid
                grid_geometry = 3
            END IF

            ! Allocating grid variables for the x-, y- and z-directions
            ALLOCATE(x_cb(-1-buff_size:m+buff_size))
            ALLOCATE(x_cc(  -buff_size:m+buff_size))
            ALLOCATE(  dx(  -buff_size:m+buff_size))
            
            IF(n == 0) RETURN; ALLOCATE(y_cb(-1-buff_size:n+buff_size))
                               ALLOCATE(y_cc(  -buff_size:n+buff_size))
                               ALLOCATE(  dy(  -buff_size:n+buff_size))
            
            IF(p == 0) RETURN; ALLOCATE(z_cb(-1-buff_size:p+buff_size))
                               ALLOCATE(z_cc(  -buff_size:p+buff_size))
                               ALLOCATE(  dz(  -buff_size:p+buff_size))


            
        END SUBROUTINE s_initialize_global_parameters_module ! -----------------


        subroutine s_initialize_nonpoly
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


            rhol0 = rhoref
            pl0   = pref
            
            allocate( pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb) )
            allocate( k_n(nb), k_v(nb), omegaN(nb) )
            allocate( Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb) ) 
            
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
            if (thermal==2) gamma_m = 1.d0
            
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
        end subroutine s_initialize_nonpoly

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
        
        
        
        
        
        SUBROUTINE s_finalize_global_parameters_module() ! ---------------------
        ! Description: Module deallocation and/or disassociation procedures
            
            INTEGER :: i
            
            ! Deallocating the variables bookkeeping the indexes of any viscous
            ! fluids and any pairs of fluids whose interfaces supported effects
            ! of surface tension
            IF(ANY(Re_size > 0)) DEALLOCATE(Re_idx)
            IF(We_size > 0) DEALLOCATE(We_idx, crv_idx)
            
            
            ! Deallocating grid variables for the x-, y- and z-directions
            DEALLOCATE(x_cb, x_cc, dx)
            
            IF(n == 0) RETURN; DEALLOCATE(y_cb, y_cc, dy)
            
            IF(p == 0) RETURN; DEALLOCATE(z_cb, z_cc, dz)
            
            DEALLOCATE(proc_coords)
            IF (parallel_io) THEN
                DEALLOCATE(start_idx)
                DO i = 1, sys_size
                    MPI_IO_DATA%var(i)%sf => NULL()
                END DO

                DEALLOCATE(MPI_IO_DATA%var)
                DEALLOCATE(MPI_IO_DATA%view)
            END IF
            
        END SUBROUTINE s_finalize_global_parameters_module ! -------------------
        
        
        !bubble routines
        SUBROUTINE s_comp_n_from_cons( vftmp,nRtmp,ntmp )
            ! vftemp is \alpha, nRtemp is nR=n*R, ntmp is n
            ! compute n from \alpha and nR(:)
            REAL(KIND(0.D0)), INTENT(IN) :: vftmp
            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: nRtmp
            REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
            REAL(KIND(0.D0)) :: nR3

            !if (nRtmp(1) < 0d0) nRtmp = 1.d-12
            CALL s_quad( nRtmp**3,nR3 )  !returns itself if NR0 = 1
            if ( nR3 < 0.d0 ) stop 'nR3 is negative'

            if (vftmp < 0.d0) stop 'vf negative'
            ntmp = DSQRT( (4.d0*pi/3.d0)*nR3/vftmp )

        END SUBROUTINE s_comp_n_from_cons

        SUBROUTINE s_comp_n_from_prim( vftmp,Rtmp,ntmp )
            ! vftemp is \alpha, Rtemp is R, ntmp is n
            ! compute n from \alpha and R(:)
            REAL(KIND(0.D0)), INTENT(IN) :: vftmp
            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: Rtmp
            REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
            REAL(KIND(0.D0)) :: R3
            integer :: i

            CALL s_quad( Rtmp**3,R3 )  !returns itself if NR0 = 1
            do i = 1,nb
                if (Rtmp(i) < 0d0) then
                    print*, Rtmp(i)
                    stop 'R neg'
                end if
            end do
            ntmp = (3.d0/(4.d0*pi)) * vftmp/R3

        END SUBROUTINE s_comp_n_from_prim

        SUBROUTINE s_quad( func,mom )

            REAL(KIND(0.D0)), DIMENSION(nb), INTENT(IN) :: func
            REAL(KIND(0.D0)), INTENT(OUT) :: mom

            
            mom = DOT_PRODUCT( weight,func )

        END SUBROUTINE s_quad

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

            !sd   = 0.3D0
            !R0mn = 0.3D0
            !R0mx = 6.D0
            
            sd   = 0.7D0
            R0mn = 0.12D0
            R0mx = 150.D0

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
            !NR0beg = 1

            print*, 'weights = ', weight(:)
            print*, 'R0s = ', R0(:)

        END SUBROUTINE s_simpson

        SUBROUTINE s_get_coeffs(a,b,alf)
            !eventually this should weight mixture vs bubble terms according to \alf
            real(kind(0d0)), intent(out)    :: a,b
            real(kind(0d0)), intent(in)     :: alf

            real(kind(0d0)) :: alf_p, alf_m

            alf_m = 0.025; alf_p = 0.045;
            !a -> mixture part; b -> bubble part
            if (bubble_switch) then
                if (alf > alf_p) then
                    a = 1d0; b = 0d0 !only mixture terms
                else if (alf < alf_m) then 
                    a = 0d0; b = 1d0 !only bubble terms
                else
                    !linear relationship between a, b and \alpha
                    a = (alf - alf_m)/(alf_p - alf_m)
                    b = (alf_p - alf)/(alf_p - alf_m) 
                end if
            else
                !a = 0d0; b = 1d0 !only bubble terms
                a = 1d0; b = 1d0 !all terms
                !a = 1d0; b = 0d0 !no bubble terms
            end if

        END SUBROUTINE s_get_coeffs

        SUBROUTINE s_get_scoeffs(a,b,s,alfgr,alfgs)
            !eventually this should weight mixture vs bubble terms according to \alf
            real(kind(0d0)), intent(out)    :: a,b,s
            real(kind(0d0)), intent(in)     :: alfgr, alfgs

            real(kind(0d0)) :: alf_p, alf_m, alf_g, amp

            !a -> mixture part; b -> bubble part
            a = 1d0; b = 1d0 !no weighting on source terms
            
            alf_g = alfgr + alfgs
            alf_m = 1.d-5; alf_p = 1.d-3;
            
            !amp = 0.1d0
            amp = 0.005d0/dt

            if (bubble_switch) then
                !s = c( H[alf_g - alf_-] alf_gr - 
                !       H[alf_+ - alf_g] alf_gs )
                if (alf_g < alf_m) then
                    s = amp*alfgr
                else if (alf_g > alf_p) then
                    s = -1.d0*amp*alfgs
                end if
            else
                s = 0d0
            end if

        END SUBROUTINE s_get_scoeffs

END MODULE m_global_parameters
