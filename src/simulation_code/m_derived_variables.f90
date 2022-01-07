!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!              numerous flow variables from the conservative and primitive ones.
!!              Currently, the available derived variables include the unadvected
!!              volume fraction, specific heat ratio, liquid stiffness, speed of
!!              sound, vorticity and the numerical Schlieren function.
MODULE m_derived_variables
    
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy

    USE m_data_output           !< Data output module

    USE m_time_steppers         !< Time-stepping algorithms
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_derived_variables_module   , &
                       s_initialize_derived_variables          , &
                       s_compute_derived_variables             , &
                       s_finalize_derived_variables_module

    !> @name Finite-difference coefficients
    !! Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: fd_coeff_x
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: fd_coeff_y
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: fd_coeff_z
    !> @}
    
    CONTAINS
        
        
        
        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module        
        SUBROUTINE s_initialize_derived_variables_module() ! ----------------------

            
            ! Allocating the variables which will store the coefficients of the
            ! centered family of finite-difference schemes. Note that sufficient
            ! space is allocated so that the coefficients up to any chosen order
            ! of accuracy may be bookkept. However, if higher than fourth-order
            ! accuracy coefficients are wanted, the formulae required to compute
            ! these coefficients will have to be implemented in the subroutine
            ! s_compute_finite_difference_coefficients.
            
            ! Allocating centered finite-difference coefficients
            IF (probe_wrt) THEN
                ALLOCATE( fd_coeff_x(-fd_number : fd_number, 0:m))
                IF (n > 0) THEN
                    ALLOCATE(fd_coeff_y(-fd_number : fd_number, 0:n))
                    IF (p > 0) THEN
                        ALLOCATE(fd_coeff_z(-fd_number : fd_number, 0:p))
                    END IF
                END IF
            END IF
            
        END SUBROUTINE s_initialize_derived_variables_module ! --------------------
        
        


        !> Allocate and open derived variables. Computing FD coefficients.
        SUBROUTINE s_initialize_derived_variables() ! -----------------------------

            ! Opening and writing header of CoM and flow probe files
            IF (proc_rank == 0) THEN
                IF (ANY(com_wrt)) THEN
                    CALL s_open_com_files()
                END IF
                IF (ANY(cb_wrt)) THEN
                    CALL s_open_cb_files()
                END IF
                IF (probe_wrt) THEN
                    CALL s_open_probe_files()
                END IF
            END IF
        
        
            ! Computing centered finite difference coefficients
            IF (probe_wrt) THEN
                CALL s_compute_finite_difference_coefficients(m,x_cc,fd_coeff_x)
                IF (n > 0) THEN
                    CALL s_compute_finite_difference_coefficients(n,y_cc,fd_coeff_y)
                    IF (p > 0) THEN
                        CALL s_compute_finite_difference_coefficients(p,z_cc,fd_coeff_z)
                    END IF
                END IF
            END IF

        END SUBROUTINE s_initialize_derived_variables ! -----------------------------




        !> Writes coherent body information, communication files, and probes.
        !!  @param t_step Current time-step
        SUBROUTINE s_compute_derived_variables(t_step) ! -----------------------

            INTEGER, INTENT(IN) :: t_step


            INTEGER :: i,j,k !< Generic loop iterators

            ! IF ((ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) .AND. (t_step > t_step_start + 2)) THEN
            IF ((ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) ) THEN
                IF (ANY(com_wrt)) THEN
                    CALL s_derive_center_of_mass(q_prim_ts(0)%vf, &
                                                 q_prim_ts(1)%vf, &
                                                 q_prim_ts(2)%vf, &
                                                 q_prim_ts(3)%vf, &
                                                 q_com)
                    CALL s_derive_higher_moments(q_prim_ts(0)%vf, moments)
                    CALL s_write_com_files(t_step,q_com,moments)
                END IF
        
                IF (ANY(cb_wrt)) THEN
                    CALL s_derive_fluid_bounds(q_prim_ts(0)%vf, bounds)
                    CALL s_derive_coherent_body(q_prim_ts(0)%vf, cb_mass)
                    CALL s_derive_centerline(q_prim_ts(0)%vf, cntrline)
                    CALL s_write_cb_files(t_step,cb_mass,bounds,cntrline)
                END IF
        
                IF (probe_wrt) THEN
                    CALL s_derive_acceleration_component(1, q_prim_ts(0)%vf, &
                                                            q_prim_ts(1)%vf, &
                                                            q_prim_ts(2)%vf, &
                                                            q_prim_ts(3)%vf, &
                                                            x_accel)
                    IF (n > 0) THEN
                        CALL s_derive_acceleration_component(2, q_prim_ts(0)%vf, &
                                                                q_prim_ts(1)%vf, &
                                                                q_prim_ts(2)%vf, &
                                                                q_prim_ts(3)%vf, &
                                                                y_accel)
                        IF (p > 0) THEN
                            CALL s_derive_acceleration_component(3, q_prim_ts(0)%vf, &
                                                                    q_prim_ts(1)%vf, &
                                                                    q_prim_ts(2)%vf, &
                                                                    q_prim_ts(3)%vf, &
                                                                    z_accel)
                        END IF
                    END IF
        
                    DO k = 0, p
                        DO j = 0, n
                            DO i = 0, m
                                IF (p > 0) THEN
                                    accel_mag(i,j,k) = SQRT(x_accel(i,j,k)**2d0 + &
                                                            y_accel(i,j,k)**2d0 + &
                                                            z_accel(i,j,k)**2d0)
                                ELSEIF (n > 0) THEN
                                    accel_mag(i,j,k) = SQRT(x_accel(i,j,k)**2d0 + &
                                                            y_accel(i,j,k)**2d0)
                                ELSE
                                    accel_mag(i,j,k) = x_accel(i,j,k)
                                END IF
                            END DO
                        END DO
                    END DO
        
                    CALL s_write_probe_files(t_step,q_cons_ts(1)%vf,accel_mag)
                END IF
            END IF

        END SUBROUTINE s_compute_derived_variables ! ---------------------------




        !>  The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
        SUBROUTINE s_compute_finite_difference_coefficients(    q,s_cc, fd_coeff_s  )
           
            INTEGER, INTENT(IN) :: q
            
            REAL(KIND(0d0)), &
            DIMENSION(-buff_size:q+buff_size), &
            INTENT(IN) :: s_cc
            
            REAL(KIND(0d0)), &
            DIMENSION(-fd_number:fd_number,0:q), &
            INTENT(INOUT) :: fd_coeff_s
            

            INTEGER :: i !< Generic loop iterator
            
            
            ! Computing the 1st order finite-difference coefficients
            IF(fd_order == 1) THEN
                DO i = 0, q
                    fd_coeff_s(-1,i) =  0d0
                    fd_coeff_s( 0,i) = -1d0 / (s_cc(i+1) - s_cc(i))
                    fd_coeff_s( 1,i) = -fd_coeff_s(0,i)
                END DO
                
            ! Computing the 2nd order finite-difference coefficients
            ELSEIF(fd_order == 2) THEN
                DO i = 0, q
                    fd_coeff_s(-1,i) = -1d0 / (s_cc(i+1) - s_cc(i-1))
                    fd_coeff_s( 0,i) =  0d0
                    fd_coeff_s( 1,i) = -fd_coeff_s(-1,i)
                END DO
                
            ! Computing the 4th order finite-difference coefficients
            ELSE
                DO i = 0, q
                    fd_coeff_s(-2,i) =  1d0 / (s_cc(i-2) - 8d0*s_cc(i-1) - s_cc(i+2) + 8d0*s_cc(i+1))
                    fd_coeff_s(-1,i) = -8d0 * fd_coeff_s(-2,i)
                    fd_coeff_s( 0,i) =  0d0
                    fd_coeff_s( 1,i) = -fd_coeff_s(-1,i)
                    fd_coeff_s( 2,i) = -fd_coeff_s(-2,i)
                END DO
                
            END IF
            
            
        END SUBROUTINE s_compute_finite_difference_coefficients ! --------------
        
        
        
        !> This subroutine receives as inputs the indicator of the
        !!      component of the acceleration that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired acceleration component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Acceleration component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_prim_vf1 Primitive variables
        !!  @param q_prim_vf2 Primitive variables
        !!  @param q_prim_vf3 Primitive variables
        !!  @param q_sf Acceleration component
        SUBROUTINE s_derive_acceleration_component(i, q_prim_vf, q_prim_vf1, &
                            q_prim_vf2, q_prim_vf3, q_sf) ! ----------
           
            
            INTEGER, INTENT(IN) :: i

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf1
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf2
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf3

            REAL(KIND(0d0)), DIMENSION(0:m,0:n,0:p), INTENT(OUT) :: q_sf
            

            INTEGER :: j,k,l,r !< Generic loop iterators
            
            ! Computing the acceleration component in the x-coordinate direction
            IF(i == 1) THEN
                DO l = 0, p
                    DO k = 0, n
                        DO j = 0, m

                            q_sf(j,k,l) = (11d0*q_prim_vf(mom_idx%beg)%sf(j,k,l) & 
                                        - 18d0*q_prim_vf1(mom_idx%beg)%sf(j,k,l) &
                                        +  9d0*q_prim_vf2(mom_idx%beg)%sf(j,k,l) &
                                        -  2d0*q_prim_vf3(mom_idx%beg)%sf(j,k,l) ) / (6d0*dt)
                            
                            DO r = -fd_number, fd_number
                                IF ( n == 0) THEN ! 1D simulation
                                    q_sf(j,k,l) = q_sf(j,k,l) &
                                        + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                        q_prim_vf(mom_idx%beg)%sf(r+j, k , l )
                                ELSEIF ( p == 0) THEN ! 2D simulation
                                    q_sf(j,k,l) = q_sf(j,k,l) &
                                        + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                        q_prim_vf(mom_idx%beg)%sf(r+j, k , l ) &
                                        + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                        q_prim_vf(mom_idx%beg)%sf( j ,r+k, l )
                                ELSE ! 3D simulation
                                    IF (grid_geometry == 3) THEN
                                        q_sf(j,k,l) = q_sf(j,k,l) &
                                            + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                            q_prim_vf(mom_idx%beg)%sf(r+j, k , l ) &
                                            + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                            q_prim_vf(mom_idx%beg)%sf( j ,r+k, l ) &
                                            + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * &
                                            q_prim_vf(mom_idx%beg)%sf( j , k ,r+l)/y_cc(k)
                                    ELSE
                                        q_sf(j,k,l) = q_sf(j,k,l) &
                                            + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                            q_prim_vf(mom_idx%beg)%sf(r+j, k , l ) &
                                            + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                            q_prim_vf(mom_idx%beg)%sf( j ,r+k, l ) &
                                            + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * &
                                            q_prim_vf(mom_idx%beg)%sf( j , k ,r+l)
                                    END IF
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
                
                
            ! Computing the acceleration component in the y-coordinate direction
            ELSEIF(i == 2) THEN
                DO l = 0, p
                    DO k = 0, n
                        DO j = 0, m
                            
                            q_sf(j,k,l) = (11d0*q_prim_vf(mom_idx%beg+1)%sf(j,k,l) & 
                                        - 18d0*q_prim_vf1(mom_idx%beg+1)%sf(j,k,l) &
                                        +  9d0*q_prim_vf2(mom_idx%beg+1)%sf(j,k,l) &
                                        -  2d0*q_prim_vf3(mom_idx%beg+1)%sf(j,k,l) ) / (6d0*dt)
                            
                            DO r = -fd_number, fd_number
                                IF ( p == 0) THEN ! 2D simulation
                                    q_sf(j,k,l) = q_sf(j,k,l) &
                                        + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                        q_prim_vf(mom_idx%beg+1)%sf(r+j, k , l ) &
                                        + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                        q_prim_vf(mom_idx%beg+1)%sf( j ,r+k, l )
                                ELSE ! 3D simulation
                                    IF (grid_geometry == 3) THEN
                                        q_sf(j,k,l) = q_sf(j,k,l) &
                                            + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                            q_prim_vf(mom_idx%beg+1)%sf(r+j, k , l ) &
                                            + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * & 
                                            q_prim_vf(mom_idx%beg+1)%sf( j ,r+k, l ) &
                                            + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * & 
                                            q_prim_vf(mom_idx%beg+1)%sf( j , k ,r+l)/y_cc(k) &
                                            -(q_prim_vf( mom_idx%end )%sf(j,k,l)**2d0)/y_cc(k)
                                    ELSE
                                        q_sf(j,k,l) = q_sf(j,k,l) &
                                            + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                            q_prim_vf(mom_idx%beg+1)%sf(r+j, k , l ) &
                                            + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                            q_prim_vf(mom_idx%beg+1)%sf( j ,r+k, l ) &
                                            + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * &
                                            q_prim_vf(mom_idx%beg+1)%sf( j , k ,r+l)
                                    END IF
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
                
                
            ! Computing the acceleration component in the z-coordinate direction
            ELSE
                DO l = 0, p
                    DO k = 0, n
                        DO j = 0, m
                            q_sf(j,k,l) = (11d0*q_prim_vf(mom_idx%end)%sf(j,k,l) & 
                                        - 18d0*q_prim_vf1(mom_idx%end)%sf(j,k,l) &
                                        +  9d0*q_prim_vf2(mom_idx%end)%sf(j,k,l) &
                                        -  2d0*q_prim_vf3(mom_idx%end)%sf(j,k,l) ) / (6d0*dt)
                            
                            DO r = -fd_number, fd_number
                                IF (grid_geometry == 3) THEN
                                    q_sf(j,k,l) = q_sf(j,k,l) &
                                        + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                        q_prim_vf(mom_idx%end)%sf(r+j, k , l ) &
                                        + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                        q_prim_vf(mom_idx%end)%sf( j ,r+k, l ) &
                                        + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * &
                                        q_prim_vf(mom_idx%end)%sf( j , k ,r+l)/y_cc(k) &
                                        +(q_prim_vf( mom_idx%end )%sf(j,k,l) * &
                                        q_prim_vf(mom_idx%beg+1)%sf(j,k,l))/y_cc(k)
                                ELSE
                                    q_sf(j,k,l) = q_sf(j,k,l) &
                                        + q_prim_vf( mom_idx%beg )%sf(j,k,l)*fd_coeff_x(r,j) * &
                                        q_prim_vf(mom_idx%end)%sf(r+j, k , l ) &
                                        + q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*fd_coeff_y(r,k) * &
                                        q_prim_vf(mom_idx%end)%sf( j ,r+k, l ) &
                                        + q_prim_vf( mom_idx%end )%sf(j,k,l)*fd_coeff_z(r,l) * &
                                        q_prim_vf(mom_idx%end)%sf( j , k ,r+l)
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END IF
            
            
        END SUBROUTINE s_derive_acceleration_component ! --------------------------




        !> This subroutine is used together with the volume fraction
        !!      model and when called upon, it computes the location of
        !!      of the center of mass for each fluid from the inputted 
        !!      primitive variables, q_prim_vf. The computed location
        !!      is then written to a formatted data file by the root
        !!      process.
        !!  @param q_prim_vf Primitive variables
        !!  @param q_prim_vf1 Primitive variables
        !!  @param q_prim_vf2 Primitive variables
        !!  @param q_prim_vf3 Primitive variables
        !!  @param q_com Mass,x-location,y-location,z-location,x-velocity,y-velocity,z-velocity,
        !!  x-acceleration, y-acceleration, z-acceleration, weighted 
        SUBROUTINE s_derive_center_of_mass(q_prim_vf,q_prim_vf1,q_prim_vf2,q_prim_vf3,q_com)
           
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf1
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf2
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf3
            REAL(KIND(0d0)), DIMENSION(num_fluids,10), INTENT(INOUT) :: q_com

            REAL(KIND(0d0)) :: xbeg,xend,ybeg,yend,zbeg,zend !<
            !! Maximum and minimum values of cell boundaries in each direction used in check for
            !! reflective BC in computation of center of mass
 

            INTEGER :: i,j,k,l !< Generic loop iterators


            REAL(KIND(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce


            REAL(KIND(0d0)) :: dV !< Discrete cell volume


            REAL(KIND(0d0)) :: cart_u_x,  cart_u_y, &
                               cart_u_x1, cart_u_y1,&
                               cart_u_x2, cart_u_y2,&
                               cart_u_x3, cart_u_y3 !<
            !! Cartesian velocities


            IF (n == 0)  THEN !1D simulation

                DO i = 1,num_fluids !Loop over individual fluids
                    IF (com_wrt(i)) THEN
                        q_com(i,:) = 0d0
                        DO l = 0, p !Loop over grid
                            DO k = 0, n
                                DO j = 0, m
                                    
                                    dV = dx(j)

                                    ! Mass
                                    q_com(i,1) = q_com(i,1) + q_prim_vf(i)%sf(j,k,l)*dV
                                    ! x-location weighted
                                    q_com(i,2) = q_com(i,2) + q_prim_vf(i)%sf(j,k,l)*dV*x_cc(j)
                                    ! x-velocity weighted
                                    q_com(i,5) = q_com(i,5) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg)%sf(j,k,l)
                                    ! x-acceleration weighted
                                    q_com(i,8) = q_com(i,8) + dV*(  11d0*( q_prim_vf(i)%sf(j,k,l) &
                                                    * q_prim_vf(mom_idx%beg)%sf(j,k,l)) &
                                                    - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg)%sf(j,k,l)) &
                                                    +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg)%sf(j,k,l)) &
                                                    -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg)%sf(j,k,l)))/(6d0*dt)
                                END DO
                            END DO
                        END DO
                        ! Sum all components across all processors using MPI_ALLREDUCE
                        IF (num_procs > 1) THEN
                            tmp = q_com(i,1)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,1))
                            tmp = q_com(i,2)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,2))
                            tmp = q_com(i,5)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,5))
                            tmp = q_com(i,8)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,8))
                        END IF

                        ! Compute quotients
                        q_com(i,2) = q_com(i,2)/q_com(i,1)
                        q_com(i,5) = q_com(i,5)/q_com(i,1)
                        q_com(i,8) = q_com(i,8)/q_com(i,1)
                    END IF
                END DO

            ELSEIF (p == 0) THEN !2D simulation

                DO i = 1,num_fluids !Loop over individual fluids
                    IF (com_wrt(i)) THEN
                        q_com(i,:) = 0d0
                        DO l = 0, p !Loop over grid
                            DO k = 0, n
                                DO j = 0, m

                                    dV = dx(j)*dy(k)

                                    ! Mass
                                    q_com(i,1) = q_com(i,1) + q_prim_vf(i)%sf(j,k,l)*dV
                                    ! x-location weighted
                                    q_com(i,2) = q_com(i,2) + q_prim_vf(i)%sf(j,k,l)*dV*x_cc(j)
                                    ! y-location weighted
                                    q_com(i,3) = q_com(i,3) + q_prim_vf(i)%sf(j,k,l)*dV*y_cc(k)
                                    ! x-velocity weighted
                                    q_com(i,5) = q_com(i,5) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg)%sf(j,k,l)
                                    ! y-velocity weighted
                                    q_com(i,6) = q_com(i,6) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg+1)%sf(j,k,l)
                                    ! x-acceleration weighted
                                    q_com(i,8) = q_com(i,8) + dV* &
                                        (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%beg)%sf(j,k,l)) &
                                        - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg)%sf(j,k,l)) &
                                        +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg)%sf(j,k,l)) &
                                        -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg)%sf(j,k,l)))/(6d0*dt)
                                    ! y-acceleration weighted
                                    q_com(i,9) = q_com(i,9) + dV * &
                                        (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%beg+1)%sf(j,k,l)) &
                                        - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg+1)%sf(j,k,l)) &
                                        +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg+1)%sf(j,k,l)) &
                                        -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg+1)%sf(j,k,l)))/(6d0*dt)
                                END DO
                            END DO
                        END DO
                        ! Sum all components across all processors using MPI_ALLREDUCE
                        IF (num_procs > 1) THEN
                            tmp = q_com(i,1)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,1))
                            tmp = q_com(i,2)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,2))
                            tmp = q_com(i,3)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,3))
                            tmp = q_com(i,5)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,5))
                            tmp = q_com(i,6)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,6))
                            tmp = q_com(i,8)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,8))
                            tmp = q_com(i,9)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,9))
                        END IF

                        ! Compute quotients
                        q_com(i,2) = q_com(i,2)/q_com(i,1)
                        q_com(i,3) = q_com(i,3)/q_com(i,1)
                        q_com(i,5) = q_com(i,5)/q_com(i,1)
                        q_com(i,6) = q_com(i,6)/q_com(i,1)
                        q_com(i,8) = q_com(i,8)/q_com(i,1)
                        q_com(i,9) = q_com(i,9)/q_com(i,1)
                    END IF
                END DO

            ELSE !3D simulation

                DO i = 1,num_fluids !Loop over individual fluids
                    IF (com_wrt(i)) THEN
                        q_com(i,:) = 0d0
                        DO l = 0, p !Loop over grid
                            DO k = 0, n
                                DO j = 0, m
                                    IF (grid_geometry == 3) THEN

                                        dV = (2d0*y_cb(k-1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                        cart_u_x  =  q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*COS(z_cc(l)) - & 
                                            q_prim_vf(mom_idx%end)%sf(j,k,l)*SIN(z_cc(l))
                                        cart_u_y  =  q_prim_vf(mom_idx%beg+1)%sf(j,k,l)*SIN(z_cc(l)) + &  
                                            q_prim_vf(mom_idx%end)%sf(j,k,l)*COS(z_cc(l))
                                        cart_u_x1 = q_prim_vf1(mom_idx%beg+1)%sf(j,k,l)*COS(z_cc(l)) - & 
                                            q_prim_vf1(mom_idx%end)%sf(j,k,l)*SIN(z_cc(l))
                                        cart_u_y1 = q_prim_vf1(mom_idx%beg+1)%sf(j,k,l)*SIN(z_cc(l)) + & 
                                            q_prim_vf1(mom_idx%end)%sf(j,k,l)*COS(z_cc(l))
                                        cart_u_x2 = q_prim_vf2(mom_idx%beg+1)%sf(j,k,l)*COS(z_cc(l)) - & 
                                            q_prim_vf2(mom_idx%end)%sf(j,k,l)*SIN(z_cc(l))
                                        cart_u_y2 = q_prim_vf2(mom_idx%beg+1)%sf(j,k,l)*SIN(z_cc(l)) + & 
                                            q_prim_vf2(mom_idx%end)%sf(j,k,l)*COS(z_cc(l))
                                        cart_u_x3 = q_prim_vf3(mom_idx%beg+1)%sf(j,k,l)*COS(z_cc(l)) - & 
                                            q_prim_vf3(mom_idx%end)%sf(j,k,l)*SIN(z_cc(l))
                                        cart_u_y3 = q_prim_vf3(mom_idx%beg+1)%sf(j,k,l)*SIN(z_cc(l)) + & 
                                            q_prim_vf3(mom_idx%end)%sf(j,k,l)*COS(z_cc(l))

                                        ! Mass
                                        q_com(i,1) = q_com(i,1) + q_prim_vf(i)%sf(j,k,l)*dV
                                        ! x-location weighted
                                        q_com(i,2) = q_com(i,2) + q_prim_vf(i)%sf(j,k,l)*dV*y_cc(k)*COS(z_cc(l))
                                        ! y-location weighted
                                        q_com(i,3) = q_com(i,3) + q_prim_vf(i)%sf(j,k,l)*dV*y_cc(k)*SIN(z_cc(l))
                                        ! z-location weighted
                                        q_com(i,4) = q_com(i,4) + q_prim_vf(i)%sf(j,k,l)*dV*x_cc(j)
                                        ! x-velocity weighted
                                        q_com(i,5) = q_com(i,5) + q_prim_vf(i)%sf(j,k,l)*dV*cart_u_x
                                        ! y-velocity weighted
                                        q_com(i,6) = q_com(i,6) + q_prim_vf(i)%sf(j,k,l)*dV*cart_u_y
                                        ! z-velocity weighted
                                        q_com(i,7) = q_com(i,7) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg)%sf(j,k,l)
                                        ! x-acceleration weighted
                                        q_com(i,8) = q_com(i,8) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)*cart_u_x ) &
                                            - 18d0*(q_prim_vf1(i)%sf(j,k,l)*cart_u_x1) &
                                            +  9d0*(q_prim_vf2(i)%sf(j,k,l)*cart_u_x2) &
                                            -  2d0*(q_prim_vf3(i)%sf(j,k,l)*cart_u_x3))/(6d0*dt)
                                        ! y-acceleration weighted
                                        q_com(i,9) = q_com(i,9) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)*cart_u_y ) &
                                                - 18d0*(q_prim_vf1(i)%sf(j,k,l)*cart_u_y1) &
                                                +  9d0*(q_prim_vf2(i)%sf(j,k,l)*cart_u_y2) &
                                                -  2d0*(q_prim_vf3(i)%sf(j,k,l)*cart_u_y3))/(6d0*dt)
                                        ! z-acceleration weighted
                                        q_com(i,10) = q_com(i,10) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%beg)%sf(j,k,l)) &
                                            - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg)%sf(j,k,l)) &
                                            +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg)%sf(j,k,l)) &
                                            -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg)%sf(j,k,l)))/(6d0*dt)
                                    ELSE

                                        dV = dx(j)*dy(k)*dz(l)

                                        ! Mass
                                        q_com(i,1) = q_com(i,1) + q_prim_vf(i)%sf(j,k,l)*dV
                                        ! x-location weighted
                                        q_com(i,2) = q_com(i,2) + q_prim_vf(i)%sf(j,k,l)*dV*x_cc(j)
                                        ! y-location weighted
                                        q_com(i,3) = q_com(i,3) + q_prim_vf(i)%sf(j,k,l)*dV*y_cc(k)
                                        ! z-location weighted
                                        q_com(i,4) = q_com(i,4) + q_prim_vf(i)%sf(j,k,l)*dV*z_cc(l)
                                        ! x-velocity weighted
                                        q_com(i,5) = q_com(i,5) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg)%sf(j,k,l)
                                        ! y-velocity weighted
                                        q_com(i,6) = q_com(i,6) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%beg+1)%sf(j,k,l)
                                        ! z-velocity weighted
                                        q_com(i,7) = q_com(i,7) + q_prim_vf(i)%sf(j,k,l)*dV*q_prim_vf(mom_idx%end)%sf(j,k,l)
                                        ! x-acceleration weighted
                                        q_com(i,8) = q_com(i,8) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%beg)%sf(j,k,l)) &
                                            - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg)%sf(j,k,l)) &
                                            +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg)%sf(j,k,l)) &
                                            -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg)%sf(j,k,l)))/(6d0*dt)
                                        ! y-acceleration weighted
                                        q_com(i,9) = q_com(i,9) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%beg+1)%sf(j,k,l)) &
                                            - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%beg+1)%sf(j,k,l)) &
                                            +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%beg+1)%sf(j,k,l)) &
                                            -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%beg+1)%sf(j,k,l)))/(6d0*dt)
                                        ! z-acceleration weighted
                                        q_com(i,10) = q_com(i,10) + dV * &
                                            (  11d0*( q_prim_vf(i)%sf(j,k,l)* q_prim_vf(mom_idx%end)%sf(j,k,l)) &
                                            - 18d0*(q_prim_vf1(i)%sf(j,k,l)*q_prim_vf1(mom_idx%end)%sf(j,k,l)) &
                                            +  9d0*(q_prim_vf2(i)%sf(j,k,l)*q_prim_vf2(mom_idx%end)%sf(j,k,l)) &
                                            -  2d0*(q_prim_vf3(i)%sf(j,k,l)*q_prim_vf3(mom_idx%end)%sf(j,k,l)))/(6d0*dt)
                                    END IF
                                END DO
                            END DO
                        END DO
                        ! Sum all components across all processors using MPI_ALLREDUCE
                        IF (num_procs > 1) THEN
                            tmp = q_com(i,1)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,1))
                            tmp = q_com(i,2)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,2))
                            tmp = q_com(i,3)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,3))
                            tmp = q_com(i,4)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,4))
                            tmp = q_com(i,5)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,5))
                            tmp = q_com(i,6)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,6))
                            tmp = q_com(i,7)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,7))
                            tmp = q_com(i,8)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,8))
                            tmp = q_com(i,9)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,9))
                            tmp = q_com(i,10)
                            CALL s_mpi_allreduce_sum(tmp,q_com(i,10))
                        END IF

                        ! Compute quotients
                        q_com(i,2) = q_com(i,2)/q_com(i,1)
                        q_com(i,3) = q_com(i,3)/q_com(i,1)
                        q_com(i,4) = q_com(i,4)/q_com(i,1)
                        q_com(i,5) = q_com(i,5)/q_com(i,1)
                        q_com(i,6) = q_com(i,6)/q_com(i,1)
                        q_com(i,7) = q_com(i,7)/q_com(i,1)
                        q_com(i,8) = q_com(i,8)/q_com(i,1)
                        q_com(i,9) = q_com(i,9)/q_com(i,1)
                        q_com(i,10) = q_com(i,10)/q_com(i,1)
                    END IF
                END DO

            END IF

            ! Find computational domain boundaries
            IF (num_procs > 1) THEN
                CALL s_mpi_allreduce_min(MINVAL(x_cb(-1:m)),xbeg)
                CALL s_mpi_allreduce_max(MAXVAL(x_cb(-1:m)),xend)
                IF (n > 0) THEN
                    CALL s_mpi_allreduce_min(MINVAL(y_cb(-1:n)),ybeg)
                    CALL s_mpi_allreduce_max(MAXVAL(y_cb(-1:n)),yend)
                    IF (p > 0) THEN
                        CALL s_mpi_allreduce_min(MINVAL(z_cb(-1:p)),zbeg)
                        CALL s_mpi_allreduce_max(MAXVAL(z_cb(-1:p)),zend)
                    END IF
                END IF
            ELSE
                xbeg = MINVAL(x_cb(-1:m))
                xend = MAXVAL(x_cb(-1:m))
                IF (n > 0) THEN
                    ybeg = MINVAL(y_cb(-1:n))
                    yend = MAXVAL(y_cb(-1:n))
                    IF (p > 0) THEN
                        zbeg = MINVAL(z_cb(-1:p))
                        zend = MAXVAL(z_cb(-1:p))
                    END IF
                END IF
            END IF

            DO i = 1, num_fluids
                IF (com_wrt(i)) THEN
                    ! Check for reflective BC in x-direction
                    IF (bc_x_glb%beg == -2) THEN
                        q_com(i,1) = q_com(i,1)*2d0
                        q_com(i,2) = xbeg
                        q_com(i,5) = 0d0
                        q_com(i,8) = 0d0
                    ELSEIF (bc_x_glb%end == -2) THEN
                        q_com(i,1) = q_com(i,1)*2d0
                        q_com(i,2) = xend
                        q_com(i,5) = 0d0
                        q_com(i,8) = 0d0
                    END IF
                    IF ( n > 0 ) THEN
                        ! Check for reflective BC in y-direction
                        IF (bc_y_glb%beg == -2) THEN
                            q_com(i,1) = q_com(i,1)*2d0
                            q_com(i,3) = ybeg
                            q_com(i,6) = 0d0
                            q_com(i,9) = 0d0
                        ELSEIF (bc_y_glb%end == -2) THEN
                            q_com(i,1) = q_com(i,1)*2d0
                            q_com(i,3) = yend
                            q_com(i,6) = 0d0
                            q_com(i,9) = 0d0
                        END IF
                        IF ( p > 0 ) THEN
                            ! Check for reflective BC in z-direction
                            IF (bc_z_glb%beg == -2) THEN
                                q_com(i,1) = q_com(i,1)*2d0
                                q_com(i,4) = zbeg
                                q_com(i,7) = 0d0
                                q_com(i,10) = 0d0
                            ELSEIF (bc_z_glb%end == -2) THEN
                                q_com(i,1) = q_com(i,1)*2d0
                                q_com(i,4) = zend
                                q_com(i,7) = 0d0
                                q_com(i,10) = 0d0
                            END IF
                            
                        END IF
                    END IF
                END IF
            END DO
                        
        END SUBROUTINE s_derive_center_of_mass ! ----------------------------------




        !>  Subroutine to compute the higher moments in an attempt to find
        !!      the maximal size of the droplet
        !!  @param q_prim_vf Primitive variables
        !!  @param moments Higher moments (2 lateral directions, 5 moment orders)
        SUBROUTINE s_derive_higher_moments(q_prim_vf, moments)

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,2,5), INTENT(INOUT) :: moments


            INTEGER :: i,r !< Generic loop iterators

            ! Using the global boundary conditions, determine method of computing
            ! higher moments for y-direction
            IF (n > 0) THEN 
                IF ((bc_y_glb%beg /= -2) .AND. (bc_y_glb%end /= -2)) THEN
                    ! Non-symmetric moments
                    CALL s_non_symmetric_moments(q_prim_vf, moments, 1)
                ELSEIF (((bc_y_glb%beg == -2) .AND. (bc_y_glb%end == -2)) & 
                            .OR. &
                    ((bc_y_glb%beg == -1) .AND. (bc_y_glb%end == -1))) THEN
                    PRINT '(A)', 'Periodic boundary conditions in y-direction. ' // &
                        'Cannot compute higher moments. Exiting...'
                    CALL s_mpi_abort()
                ELSE
                    CALL s_symmetric_moments(q_prim_vf, moments, 1)
                END IF

                IF (p > 0) THEN
                    IF ((bc_z_glb%beg /= -2) .AND. (bc_z_glb%end /= -2)) THEN
                        ! Non-symmetric moments
                        CALL s_non_symmetric_moments(q_prim_vf, moments, 2)
                    ELSEIF (((bc_z_glb%beg == -2) .AND. (bc_z_glb%end == -2)) & 
                                .OR. & 
                        ((bc_z_glb%beg == -1) .AND. (bc_z_glb%end == -1))) THEN
                        PRINT '(A)', 'Periodic boundary conditions in z-direction. ' // &
                            'Cannot compute higher moments. Exiting...'
                        CALL s_mpi_abort()
                    ELSE
                        CALL s_symmetric_moments(q_prim_vf, moments, 2)
                    END IF
                END IF
            ELSE !1D simulation
                DO i = 1,num_fluids !Loop over individual fluids
                    IF (com_wrt(i)) THEN
                        DO r = 1, 5
                            IF (moment_order(r) /= dflt_int) THEN
                                moments(i,:,r) = 0d0
                            ELSE
                                moments(i,:,r) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            END IF

        END SUBROUTINE s_derive_higher_moments ! -----------------------------------------




        !> Compute non-symmetric moments
        !! @param q_prim_vf Primitive variables
        !! @param moments Higher moments(2 lateral directions, 5 moment orders)
        !! @param dir Current lateral direction
        SUBROUTINE s_non_symmetric_moments(q_prim_vf, moments, dir) ! ---------------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,2,5), INTENT(INOUT) :: moments
            INTEGER, INTENT(IN) :: dir

            REAL(KIND(0d0)), DIMENSION(num_fluids,5) :: pos_numer, neg_numer, pos_denom, neg_denom !<
            !! Numerator and denominator place holders for computation

            REAL(KIND(0d0)) :: numer_weight     !< Numerator weight
            REAL(KIND(0d0)) :: main_term        !< Constant term in both numerator and denominator
            REAL(KIND(0d0)) :: dV               !< Discrete cell volume
            REAL(KIND(0d0)) :: cart_x, cart_y   !< Cartesian x- and y-locations
            INTEGER :: i,j,k,l,r    !< Generic loop iterators
            REAL(KIND(0d0)) :: tmp  !< Temporary variable to store quantity for mpi_allreduce

            DO i = 1, num_fluids
                IF (com_wrt(i)) THEN
                pos_numer(i,:) = 0d0
                neg_numer(i,:) = 0d0
                pos_denom(i,:) = 0d0
                neg_denom(i,:) = 0d0
                    DO r = 1, 5
                        IF (moment_order(r) /= dflt_int) THEN
                            DO l = 0, p
                                DO k = 0, n
                                    DO j = 0, m
                                        IF (q_prim_vf(i)%sf(j,k,l) > 5d-1) THEN
                                            IF (grid_geometry == 3) THEN
                                                dV = (2d0*y_cb(k-1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                                cart_x = y_cc(k)*COS(z_cc(l))
                                                cart_y = y_cc(k)*SIN(z_cc(l))
                                                main_term =       q_prim_vf(i+E_idx)%sf(j,k,l)  * &
                                                           (1d0 - q_prim_vf(i+E_idx)%sf(j,k,l)) * dV
                                                IF ((dir == 1) .AND. (cart_x >= 0d0)) THEN
                                                    numer_weight = cart_x**moment_order(r)
    
                                                    pos_numer(i,r) = pos_numer(i,r) + numer_weight * main_term
                                                    pos_denom(i,r) = pos_denom(i,r) + main_term
                                                ELSEIF ((dir == 1) .AND. (cart_x < 0d0)) THEN
                                                    numer_weight = cart_x**moment_order(r)
    
                                                    neg_numer(i,r) = neg_numer(i,r) + numer_weight * main_term
                                                    neg_denom(i,r) = neg_denom(i,r) + main_term
                                                ELSEIF ((dir == 2) .AND. (cart_y >= 0d0)) THEN
                                                    numer_weight = cart_y**moment_order(r)
    
                                                    pos_numer(i,r) = pos_numer(i,r) + numer_weight * main_term
                                                    pos_denom(i,r) = pos_denom(i,r) + main_term
                                                ELSEIF ((dir == 2) .AND. (cart_y < 0d0)) THEN
                                                    numer_weight = cart_y**moment_order(r)
    
                                                    neg_numer(i,r) = neg_numer(i,r) + numer_weight * main_term
                                                    neg_denom(i,r) = neg_denom(i,r) + main_term
                                                END IF
                                            ELSE
                                                IF (n > 0) THEN
                                                    main_term =       q_prim_vf(i+E_idx)%sf(j,k,l)  * &
                                                               (1d0 - q_prim_vf(i+E_idx)%sf(j,k,l)) * &
                                                               dx(j)*dy(k)
                                                    IF (p > 0) THEN
                                                        main_term = main_term * dz(l)
                                                    END IF
                                                END IF
                                                IF ((dir == 1) .AND. (y_cc(k) >= 0d0)) THEN
                                                    numer_weight = y_cc(k)**moment_order(r)
    
                                                    pos_numer(i,r) = pos_numer(i,r) + numer_weight * main_term
                                                    pos_denom(i,r) = pos_denom(i,r) + main_term
                                                ELSEIF ((dir == 1) .AND. (y_cc(k) < 0d0)) THEN
                                                    numer_weight = y_cc(k)**moment_order(r)
    
                                                    neg_numer(i,r) = neg_numer(i,r) + numer_weight * main_term
                                                    neg_denom(i,r) = neg_denom(i,r) + main_term
                                                ELSEIF ((dir == 2) .AND. (z_cc(l) >= 0d0)) THEN
                                                    numer_weight = z_cc(l)**moment_order(r)
    
                                                    pos_numer(i,r) = pos_numer(i,r) + numer_weight * main_term
                                                    pos_denom(i,r) = pos_denom(i,r) + main_term
                                                ELSEIF ((dir == 2) .AND. (z_cc(l) < 0d0)) THEN
                                                    numer_weight = z_cc(l)**moment_order(r)
    
                                                    neg_numer(i,r) = neg_numer(i,r) + numer_weight * main_term
                                                    neg_denom(i,r) = neg_denom(i,r) + main_term
                                                END IF
                                            END IF
                                        END IF
                                    END DO
                                END DO
                            END DO
                            ! Sum all components across all procs using MPI_ALLREDUCE
                            IF (num_procs > 1) THEN
                                tmp = pos_numer(i,r)
                                CALL s_mpi_allreduce_sum(tmp,pos_numer(i,r))
                                tmp = neg_numer(i,r)
                                CALL s_mpi_allreduce_sum(tmp,neg_numer(i,r))
                                tmp = pos_denom(i,r)
                                CALL s_mpi_allreduce_sum(tmp,pos_denom(i,r))
                                tmp = neg_denom(i,r)
                                CALL s_mpi_allreduce_sum(tmp,neg_denom(i,r))
                            END IF
                            ! Compute quotients and sum to get total moment
                            moments(i,dir,r) = (pos_numer(i,r)/pos_denom(i,r))**(1d0/moment_order(r)) + &
                                       (neg_numer(i,r)/neg_denom(i,r))**(1d0/moment_order(r))
                        ELSE
                            moments(i,dir,r) = dflt_real
                        END IF
                    END DO
                END IF
            END DO

        END SUBROUTINE s_non_symmetric_moments ! ------------------------------------------ 




        !> Compute symmetric moments
        !! @param q_prim_vf Primitive variables
        !! @param moments Higher moments(2 lateral directions, 5 moment orders)
        !! @param dir Current lateral direction
        SUBROUTINE s_symmetric_moments(q_prim_vf, moments, dir) ! ---------------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,2,5), INTENT(INOUT) :: moments
            INTEGER, INTENT(IN) :: dir


            REAL(KIND(0d0)), DIMENSION(num_fluids,5) :: numer, denom !<
            !! Numerator and denominator place holders for computation

            REAL(KIND(0d0)) :: numer_weight     !< Numerator weight
            REAL(KIND(0d0)) :: main_term        !< Constant term in both numerator and denominator
            REAL(KIND(0d0)) :: dV               !< Discrete cell volume
            REAL(KIND(0d0)) :: cart_x, cart_y   !< Cartesian x- and y-locations
            REAL(KIND(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

            INTEGER :: i,j,k,l,r !< Generic loop iterators


            DO i = 1, num_fluids
                IF (com_wrt(i)) THEN
                numer(i,:) = 0d0
                denom(i,:) = 0d0
                    DO r = 1, 5
                        IF (moment_order(r) /= dflt_int) THEN
                            DO l = 0, p
                                DO k = 0, n
                                    DO j = 0, m
                                        IF (q_prim_vf(i)%sf(j,k,l) > 5d-1) THEN
                                            IF (grid_geometry == 3) THEN
                                                dV = (2d0*y_cb(k-1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                                cart_x = y_cc(k)*COS(z_cc(l))
                                                cart_y = y_cc(k)*SIN(z_cc(l))
                                                main_term =       q_prim_vf(i+E_idx)%sf(j,k,l)  * &
                                                           (1d0 - q_prim_vf(i+E_idx)%sf(j,k,l)) * dV
                                                IF (dir == 1) THEN
                                                    numer_weight = cart_x**moment_order(r)
        
                                                    numer(i,r) = numer(i,r) + numer_weight * main_term
                                                    denom(i,r) = denom(i,r) + main_term
                                                ELSEIF (dir == 2) THEN
                                                    numer_weight = cart_y**moment_order(r)
        
                                                    numer(i,r) = numer(i,r) + numer_weight * main_term
                                                    denom(i,r) = denom(i,r) + main_term
                                                END IF
                                            ELSE
                                                IF (n > 0) THEN
                                                    main_term =       q_prim_vf(i+E_idx)%sf(j,k,l)  * &
                                                               (1d0 - q_prim_vf(i+E_idx)%sf(j,k,l)) * &
                                                               dx(j)*dy(k)
                                                    IF (p > 0) THEN
                                                        main_term = main_term * dz(l)
                                                    END IF
                                                END IF
                                                IF (dir == 1) THEN
                                                    numer_weight = y_cc(k)**moment_order(r)
        
                                                    numer(i,r) = numer(i,r) + numer_weight * main_term
                                                    denom(i,r) = denom(i,r) + main_term
                                                ELSEIF (dir == 2) THEN
                                                    numer_weight = z_cc(l)**moment_order(r)
        
                                                    numer(i,r) = numer(i,r) + numer_weight * main_term
                                                    denom(i,r) = denom(i,r) + main_term
                                                END IF
                                            END IF
                                        END IF
                                    END DO
                                END DO
                            END DO
                            ! Sum all components across all procs using MPI_ALLREDUCE
                            IF (num_procs > 1) THEN
                                tmp = numer(i,r)
                                CALL s_mpi_allreduce_sum(tmp,numer(i,r))
                                tmp = denom(i,r)
                                CALL s_mpi_allreduce_sum(tmp,denom(i,r))
                            END IF
                            ! Compute quotients and sum to get total moment
                            moments(i,dir,r) = (numer(i,r)/denom(i,r))**(1d0/moment_order(r))
                        ELSE
                            moments(i,dir,r) = dflt_real
                        END IF
                    END DO
                END IF
            END DO

        END SUBROUTINE s_symmetric_moments ! ------------------------------------------ 




        !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the min and max bounds  of the 
        !!      fluid in each direction in the domain. 
        !!  @param q_prim_vf Primitive variables
        !!  @param bounds Variables storing the min and max bounds  of the fluids
        SUBROUTINE s_derive_fluid_bounds(q_prim_vf, bounds) ! -----------------------


            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,5,6), INTENT(INOUT) :: bounds


            REAL(KIND(0d0)) :: cart_x, cart_y, cart_z !< Cartesian x,y,z-locations
            REAL(KIND(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

            INTEGER :: i,j,k,l,r !< Generic loop iterators


            IF (n == 0) THEN ! 1D simulation
                DO i = 1,num_fluids !Loop over individual fluids
                    IF (cb_wrt(i)) THEN
                        bounds(i,:,1) = -1d0*dflt_real ! 1d6
                        bounds(i,:,2) = dflt_real ! -1d6
                        DO r = 1, 5
                            IF (threshold_mf(r) /= dflt_real) THEN
                                DO l = 0, p !Loop over grid
                                    DO k = 0, n
                                        DO j = 0, m
                                            IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                .AND. (x_cb(j-1) <= bounds(i,r,1))) THEN
                                                bounds(i,r,1) = x_cb(j-1)
                                            ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                    .AND. (x_cb(j) >= bounds(i,r,2))) THEN
                                                bounds(i,r,2) = x_cb(j)
                                            END IF
                                        END DO
                                    END DO
                                END DO

                                IF (num_procs > 1) THEN
                                    tmp = bounds(i,r,1)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,1))
                                    tmp = bounds(i,r,2)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,2))
                                END IF
                            ELSE
                                bounds(i,r,1) = dflt_real
                                bounds(i,r,2) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            ELSEIF (p == 0) THEN ! 2D simulation
                DO i = 1,num_fluids !Loop over individual fluids
                    IF (cb_wrt(i)) THEN
                        bounds(i,:,1) = -1d0*dflt_real ! 1d6
                        bounds(i,:,2) = dflt_real ! -1d6
                        bounds(i,:,3) = -1d0*dflt_real ! 1d6
                        bounds(i,:,4) = dflt_real ! -1d6
                        DO r = 1, 5
                            IF (threshold_mf(r) /= dflt_real) THEN 
                                DO l = 0, p ! Loop over grid
                                    DO k = 0, n
                                        DO j = 0, m 
                                            IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                .AND. (x_cb(j-1) <= bounds(i,r,1))) THEN
                                                bounds(i,r,1) = x_cb(j-1)
                                            ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                    .AND. (x_cb(j) >= bounds(i,r,2))) THEN
                                                bounds(i,r,2) = x_cb(j)
                                            END IF
                                            IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                .AND. (y_cb(k-1) <= bounds(i,r,3))) THEN
                                                bounds(i,r,3) = y_cb(k-1)
                                            ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                    .AND. (y_cb(k) >= bounds(i,r,4))) THEN
                                                bounds(i,r,4) = y_cb(k)
                                            END IF
                                        END DO
                                    END DO
                                END DO

                                IF (num_procs > 1) THEN
                                    tmp = bounds(i,r,1)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,1))
                                    tmp = bounds(i,r,2)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,2))
                                    tmp = bounds(i,r,3)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,3))
                                    tmp = bounds(i,r,4)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,4))
                                END IF
                            ELSE
                                bounds(i,r,1) = dflt_real
                                bounds(i,r,2) = dflt_real
                                bounds(i,r,3) = dflt_real
                                bounds(i,r,4) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            ELSE ! 3D simulation
                DO i = 1,num_fluids !Loop over individual fluids
                    IF (cb_wrt(i)) THEN
                        bounds(i,:,1) = -1d0*dflt_real ! 1d6
                        bounds(i,:,2) = dflt_real ! -1d6
                        bounds(i,:,3) = -1d0*dflt_real ! 1d6
                        bounds(i,:,4) = dflt_real ! -1d6
                        bounds(i,:,5) = -1d0*dflt_real ! 1d6
                        bounds(i,:,6) = dflt_real ! -1d6
                        DO r = 1, 5 
                            IF (threshold_mf(r) /= dflt_real) THEN
                                DO l = 0, p ! Loop over grid
                                    DO k = 0, n
                                        DO j = 0, m 
                                            IF (grid_geometry == 3) THEN
                                                cart_x = y_cc(k)*COS(z_cc(l))
                                                cart_y = y_cc(k)*SIN(z_cc(l))
                                                cart_z = x_cc(j)
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &
                                                    .AND. (cart_x <= bounds(i,r,1))) THEN
                                                    bounds(i,r,1) = cart_x
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (cart_x >= bounds(i,r,2))) THEN
                                                    bounds(i,r,2) = cart_x
                                                END IF
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (cart_y <= bounds(i,r,3))) THEN
                                                    bounds(i,r,3) = cart_y
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &  
                                                    .AND. (cart_y >= bounds(i,r,4))) THEN
                                                    bounds(i,r,4) = cart_y
                                                END IF
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (cart_z <= bounds(i,r,5))) THEN
                                                    bounds(i,r,5) = cart_z
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &  
                                                    .AND. (cart_z >= bounds(i,r,6))) THEN
                                                    bounds(i,r,6) = cart_z
                                                END IF
                                            ELSE
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &  
                                                    .AND. (x_cb(j-1) <= bounds(i,r,1))) THEN
                                                    bounds(i,r,1) = x_cb(j-1)
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (x_cb(j) >= bounds(i,r,2))) THEN
                                                    bounds(i,r,2) = x_cb(j)
                                                END IF
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) &  
                                                    .AND. (y_cb(k-1) <= bounds(i,r,3))) THEN
                                                    bounds(i,r,3) = y_cb(k-1)
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (y_cb(k) >= bounds(i,r,4))) THEN
                                                    bounds(i,r,4) = y_cb(k)
                                                END IF
                                                IF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (z_cb(l-1) <= bounds(i,r,5))) THEN
                                                    bounds(i,r,5) = z_cb(l-1)
                                                ELSEIF ((q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) & 
                                                    .AND. (z_cb(l) >= bounds(i,r,6))) THEN
                                                    bounds(i,r,6) = z_cb(l)
                                                END IF
                                            END IF
                                        END DO
                                    END DO
                                END DO

                                IF (num_procs > 1) THEN
                                    tmp = bounds(i,r,1)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,1))
                                    tmp = bounds(i,r,2)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,2))
                                    tmp = bounds(i,r,3)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,3))
                                    tmp = bounds(i,r,4)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,4))
                                    tmp = bounds(i,r,5)
                                    CALL s_mpi_allreduce_min(tmp,bounds(i,r,5))
                                    tmp = bounds(i,r,6)
                                    CALL s_mpi_allreduce_max(tmp,bounds(i,r,6))
                                END IF
                            ELSE
                                bounds(i,r,1) = dflt_real
                                bounds(i,r,2) = dflt_real
                                bounds(i,r,3) = dflt_real
                                bounds(i,r,4) = dflt_real
                                bounds(i,r,5) = dflt_real
                                bounds(i,r,6) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            END IF

        END SUBROUTINE s_derive_fluid_bounds ! -------------------------------------------



        !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the total mass of a fluid in the 
        !!      entire domain for which the volume fraction is greater than a
        !!      threshold value. This gives the mass of the coherent body.
        !!  @param q_prim_vf Primitive variables
        !!  @param cb_mass Coherent body mass
        SUBROUTINE s_derive_coherent_body(q_prim_vf, cb_mass) ! --------------------------


            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,10), INTENT(INOUT) :: cb_mass


            REAL(KIND(0d0)) :: dV  !< Discrete cell volume
            INTEGER :: i,j,k,l,r   !< Generic loop iterators
            REAL(KIND(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

            DO i = 1,num_fluids !Loop over individual fluids
                IF (cb_wrt(i)) THEN
                    cb_mass(i,:) = 0d0
                    DO r = 1, 5 ! Volume fraction threshold values
                        IF (threshold_mf(r) /= dflt_real) THEN
                            DO l = 0, p !Loop over grid
                                DO k = 0, n
                                    DO j = 0, m
                                        IF (q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) THEN
                                            IF (n == 0) THEN
                                                dV = dx(j)
                                            ELSEIF (p == 0) THEN
                                                dV = dx(j)*dy(k)
                                            ELSE
                                                IF (grid_geometry == 3) THEN
                                                    dV = (2d0*y_cb(k-1)*dy(k) + dy(k)**2d0)/2d0*dx(j)*dz(l)
                                                ELSE
                                                    dV = dx(j)*dy(k)*dz(l)
                                                END IF
                                            END IF
                                            cb_mass(i,r) = cb_mass(i,r) + q_prim_vf(i)%sf(j,k,l)*dV ! Mass
                                            cb_mass(i,r+5) = cb_mass(i,r+5) + dV ! Volume
                                        END IF
                                    END DO
                                END DO
                            END DO

                            IF (num_procs > 1) THEN
                                tmp = cb_mass(i,r)
                                CALL s_mpi_allreduce_sum(tmp,cb_mass(i,r))
                                tmp = cb_mass(i,r+5)
                                CALL s_mpi_allreduce_sum(tmp,cb_mass(i,r+5))
                            END IF
                        ELSE
                            cb_mass(i,r) = dflt_real
                            cb_mass(i,r+5) = dflt_real
                        END IF
                    END DO
                END IF
            END DO

            DO i = 1, num_fluids
                IF (cb_wrt(i)) THEN
                    DO r = 1, 5
                        IF (threshold_mf(r) /= dflt_real) THEN
                            ! Check for reflective BC in x-direction
                            IF (bc_x_glb%beg == -2) THEN
                                cb_mass(i,r) = cb_mass(i,r)*2d0
                                cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                            ELSEIF (bc_x_glb%end == -2) THEN
                                cb_mass(i,r) = cb_mass(i,r)*2d0
                                cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                            END IF
                            IF ( n > 0 ) THEN
                                ! Check for reflective BC in y-direction
                                IF (bc_y_glb%beg == -2) THEN
                                    cb_mass(i,r) = cb_mass(i,r)*2d0
                                    cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                                ELSEIF (bc_y_glb%end == -2) THEN
                                    cb_mass(i,r) = cb_mass(i,r)*2d0
                                    cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                                END IF
                                IF ( p > 0 ) THEN
                                    ! Check for reflective BC in z-direction
                                    IF (bc_z_glb%beg == -2) THEN
                                        cb_mass(i,r) = cb_mass(i,r)*2d0
                                        cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                                    ELSEIF (bc_z_glb%end == -2) THEN
                                        cb_mass(i,r) = cb_mass(i,r)*2d0
                                        cb_mass(i,r+5) = cb_mass(i,r+5)*2d0
                                    END IF
                                    
                                END IF
                            END IF
                        END IF
                    END DO
                END IF
            END DO


        END SUBROUTINE s_derive_coherent_body ! ------------------------------



        !>  This subroutine is used together with the volume fraction model
        !!      and when called upon, it computes the centerline length of the
        !!      fluid.
        !!  @param q_prim_vf Primitive variables
        !!  @param cntrline Variables storing the centerline length of the fluids
        SUBROUTINE s_derive_centerline(q_prim_vf, cntrline) ! --------------------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            REAL(KIND(0d0)), DIMENSION(num_fluids,5), INTENT(INOUT) :: cntrline


            REAL(KIND(0d0)), DIMENSION(5) :: cntrmin, cntrmax !< Placeholders
            REAL(KIND(0d0)) :: tmp !< Temporary variable to store quantity for mpi_allreduce

            INTEGER :: i,j,k,l,r !< Generic loop iterators

            IF (n == 0) THEN ! 1D simulation
                DO i = 1, num_fluids
                    IF (cb_wrt(i)) THEN
                        DO r = 1, 5
                            IF (threshold_mf(r) /= dflt_real) THEN
                                cntrline(i,r) = 0d0
                            ELSE
                                cntrline(i,r) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            ELSEIF ((p == 0) .OR. (grid_geometry == 3)) THEN ! 2D simulation
                DO i = 1, num_fluids
                    IF (cb_wrt(i)) THEN
                        cntrmin(:) = -1d0*dflt_real
                        cntrmax(:) = dflt_real
                        DO r = 1, 5
                            IF (threshold_mf(r) /= dflt_real) THEN
                                DO l = 0, p
                                    DO k = 0, n
                                        DO j = 0, m
                                            IF (    (y_cb(k-1) == 0d0) .AND. &
                                                (q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) .AND. &
                                                (x_cb(j-1) <= cntrmin(r)) ) THEN
                                                    cntrmin(r) = x_cb(j-1)
                                            ELSEIF ((y_cb(k-1) == 0d0) .AND. &
                                                (q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) .AND. &
                                                (x_cb(j) >= cntrmax(r)) ) THEN
                                                    cntrmax(r) = x_cb(j)
                                            END IF
                                        END DO
                                    END DO
                                END DO

                                IF (num_procs > 1) THEN
                                    tmp = cntrmin(r)
                                    CALL s_mpi_allreduce_min(tmp,cntrmin(r))
                                    tmp = cntrmax(r)
                                    CALL s_mpi_allreduce_max(tmp,cntrmax(r))
                                END IF

                                cntrline(i,r) = cntrmax(r) - cntrmin(r)
                            ELSE
                                cntrline(i,r) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            ELSE ! 3D simulation
                DO i = 1, num_fluids
                    IF (cb_wrt(i)) THEN
                        cntrmin(:) = -1d0*dflt_real
                        cntrmax(:) = dflt_real
                        DO r = 1, 5
                            IF (threshold_mf(r) /= dflt_real) THEN
                                DO l = 0, p
                                    DO k = 0, n
                                        DO j = 0, m
                                            IF (    (y_cb(k-1) == 0d0) .AND. &
                                                    (z_cb(l-1) == 0d0) .AND. &
                                                (q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) .AND. &
                                                (x_cb(j-1) <= cntrmin(r)) ) THEN
                                                    cntrmin(r) = x_cb(j-1)
                                            ELSEIF ((y_cb(k-1) == 0d0) .AND. &
                                                    (z_cb(l-1) == 0d0) .AND. &
                                                (q_prim_vf(i+E_idx)%sf(j,k,l) >= threshold_mf(r)) .AND. &
                                                (x_cb(j) >= cntrmax(r)) ) THEN
                                                    cntrmax(r) = x_cb(j)
                                            END IF
                                        END DO
                                    END DO
                                END DO

                                IF (num_procs > 1) THEN
                                    tmp = cntrmin(r)
                                    CALL s_mpi_allreduce_min(tmp,cntrmin(r))
                                    tmp = cntrmax(r)
                                    CALL s_mpi_allreduce_max(tmp,cntrmax(r))
                                END IF

                                cntrline(i,r) = cntrmax(r) - cntrmin(r)
                            ELSE
                                cntrline(i,r) = dflt_real
                            END IF
                        END DO
                    END IF
                END DO
            END IF

        END SUBROUTINE s_derive_centerline ! ----------------------------------





        !> Deallocation procedures for the module
        SUBROUTINE s_finalize_derived_variables_module() ! -------------------

            ! Closing CoM and flow probe files
            IF (proc_rank == 0) THEN
                IF (ANY(com_wrt)) THEN
                    CALL s_close_com_files()
                END IF
                IF (ANY(cb_wrt)) THEN
                    CALL s_close_cb_files()
                END IF
                IF (probe_wrt) THEN
                    CALL s_close_probe_files()
                END IF
            END IF
        
            ! Deallocating the variables that might have been used to bookkeep
            ! the finite-difference coefficients in the x-, y- and z-directions
            IF(ALLOCATED(fd_coeff_x)) DEALLOCATE(fd_coeff_x)
            IF(ALLOCATED(fd_coeff_y)) DEALLOCATE(fd_coeff_y)
            IF(ALLOCATED(fd_coeff_z)) DEALLOCATE(fd_coeff_z)
            
            
        END SUBROUTINE s_finalize_derived_variables_module ! -----------------
        
        
        
        
        
END MODULE m_derived_variables
