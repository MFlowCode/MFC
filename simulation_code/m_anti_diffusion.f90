!>
!! @file m_anti_diffusion.f90
!! @brief This module is used to compute the anti-diffusion equation
!! @author spencer
!! @version 1.1
!! @date 1/1/1
MODULE m_anti_diffusion


    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types

    USE m_global_parameters    !< Definitions of the global parameters

    USE m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    USE m_data_output          !< Run-time info & solution data output procedures

    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy

    USE m_time_steppers        !< Time-stepping algorithms

    USE m_variables_conversion !< Variable conversion
    ! ==========================================================================

    IMPLICIT NONE


    TYPE(scalar_field), PRIVATE, ALLOCATABLE, DIMENSION(:) :: ad_rhs_vf !<
    !! Cell-average RHS variables at the current time-stage


    TYPE(scalar_field), PRIVATE, ALLOCATABLE, DIMENSION(:) :: ad_prim_vf !<
    !! Cell-average RHS variables at the current time-stage


    !> @name Variables to compute grad alpha
    !> @{
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:) :: grad_alpha
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:) :: grad_alpharho1
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:) :: grad_alpharho2
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:) :: grad_rho
    !> @}


    !> @name Temporary array variables
    !> @{
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:) :: ad_a, ad_b, ad_c
    !> @}


    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:) :: face_flux !< Face flux

    !> @name Diffusion coefficients
    !> @{
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:) :: diff_coeff_x,diff_coeff_y,diff_coeff_z
    !> @}


    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp !<
    !!Temp variables to store primitives


    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: gm_alpha !<
    !!gamma_alpha_qp -- Dummy variable in convert from conservative to primitive

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    TYPE(bounds_info) :: ad_ix, ad_iy, ad_iz
    !> @}

    CONTAINS



        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
        SUBROUTINE s_initialize_anti_diffusion_module()
 
            INTEGER :: i,j

            ! Allocating the cell-average RHS variables
            ALLOCATE(ad_rhs_vf(1:sys_size))

            DO i = 1, sys_size
                ALLOCATE(ad_rhs_vf(i)%sf(0:m,0:n,0:p))
            END DO

            ! Setting the indical bounds in the x-, y- and z-directions
            ad_ix%beg = -buff_size; ad_ix%end = m + buff_size

            IF(n > 0) THEN

                ad_iy%beg = -buff_size; ad_iy%end = n + buff_size

                IF(p > 0) THEN
                    ad_iz%beg = -buff_size; ad_iz%end = p + buff_size
                ELSE
                    ad_iz%beg = 0; ad_iz%end = 0
                END IF

            ELSE

                ad_iy%beg = 0; ad_iy%end = 0
                ad_iz%beg = 0; ad_iz%end = 0

            END IF

            ! Allocating the cell-average RHS variables
            ALLOCATE(ad_prim_vf(1:sys_size))

            DO i = mom_idx%beg, E_idx
                ALLOCATE(ad_prim_vf(i)%sf(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
            END DO

            !Allocataing gm_alpha_qp
            ALLOCATE(gm_alpha(1:num_fluids))
            DO i = 1,num_fluids
                ALLOCATE(gm_alpha(i)%sf(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
            END DO

            !Allocating grad alpha
            ALLOCATE(grad_alpha(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
            ALLOCATE(grad_alpharho1(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
            ALLOCATE(grad_alpharho2(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
            ALLOCATE(grad_rho(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))

            !Allocating face_flux
            ALLOCATE(face_flux(1:4,-1:m,-1:n,-1:p))

            !Allocating diffusion coefficent
            ALLOCATE(diff_coeff_x(0:m,0:n,0:p))
            ALLOCATE(diff_coeff_y(0:m,0:n,0:p))
            ALLOCATE(diff_coeff_z(0:m,0:n,0:p))

            !Allocating temp variable
            ALLOCATE(temp(1:7,0:m,0:n,0:p))

        END SUBROUTINE s_initialize_anti_diffusion_module



        !> Subroutine to compute anti-diffusion equation
        SUBROUTINE s_compute_anti_diffusion()

            INTEGER :: i, j, k, r, no_of_itr
            REAL(KIND(0d0)) :: dtau, maxD_loc, maxD_glb
            REAL(KIND(0d0)), DIMENSION(num_dims) :: vel
            REAL(KIND(0d0)), DIMENSION(0:m,0:n,0:p) :: Lheaviside, U0, velmag
            REAL(KIND(0d0)) :: U0_loc, U0_glb

            no_of_itr = 1

            IF(num_fluids /= 2) THEN 
                PRINT '(A)', 'Error in number of species for anti-diffusion...exiting'
                CALL s_mpi_abort()
            END IF

            !convert conservative to primitive
            DO i = 1, cont_idx%end
                ad_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            DO i = adv_idx%beg, adv_idx%end
                ad_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            CALL s_convert_conservative_to_primitive_variables(q_cons_ts(1)%vf, ad_prim_vf, gm_alpha,ad_ix,ad_iy,ad_iz)

            !Setting diffusion coefficient to velocity magnitude
            DO i = 0, m
                DO j = 0, n
                    DO k = 0, p
                        diff_coeff_x(i,j,k) = ABS(ad_prim_vf(mom_idx%beg)%sf(i,j,k))
                        diff_coeff_y(i,j,k) = ABS(ad_prim_vf(mom_idx%beg+1)%sf(i,j,k))
                        diff_coeff_z(i,j,k) = ABS(ad_prim_vf(mom_idx%end)%sf(i,j,k))
                    END DO
                END DO
            END DO
            maxD_loc = MAX(MAXVAL(diff_coeff_x),MAXVAL(diff_coeff_y),MAXVAL(diff_coeff_z))
            IF (num_procs > 1) THEN
                CALL s_mpi_allreduce_max(maxD_loc,maxD_glb)
            ELSE
                maxD_glb = maxD_loc
            END IF

            IF(n.eq.0) THEN
                dtau = MIN(MINVAL(dx)**2d0/(2d0*maxD_glb),dt)
            ELSE IF(p.eq.0) THEN
                IF(cyl_coord.NEQV..TRUE.) THEN
                    dtau = MIN((MIN((minval(dx))**2d0,(minval(dy))**2d0))/(4d0*maxd_glb),dt)
                ELSE
                    dtau = 1.0d1*MIN((MIN((minval(dx))**2d0,(minval(dy))**2d0))/(6d0*maxd_glb),dt)
                END IF
            ELSE
                dtau = MIN((minval(dx))**2d0,(minval(dy))**2d0)
                dtau = MIN((MIN( dtau ,(minval(dz))**2d0))/(6d0*maxd_glb),dt)
            END IF

            DO i = 0, m
                DO j = 0, n
                    DO k = 0, p
                        IF ((ad_prim_vf(adv_idx%beg)%sf(i,j,k) > 1d-6) &
                                         .AND.                         &
                            (ad_prim_vf(adv_idx%beg)%sf(i,j,k) < (1d0-1d-6))) THEN
                            Lheaviside(i,j,k) = 1d0
                        ELSE
                            Lheaviside(i,j,k) = 0d0
                        END IF

                        DO r = 1, num_dims
                            vel(r) = ad_prim_vf(cont_idx%end+r)%sf(i,j,k)
                        END DO

                        velmag(i,j,k) = SQRT(DOT_PRODUCT(vel,vel))

                        U0(i,j,k) = 4d0*ad_prim_vf(adv_idx%beg)%sf(i,j,k)   * &
                                    (1d0-ad_prim_vf(adv_idx%beg)%sf(i,j,k)) * &
                                    velmag(i,j,k)
                    END DO
                END DO
            END DO

            U0_loc = MAXVAL(U0)
            IF (num_procs > 1) THEN
                CALL s_mpi_allreduce_max(U0_loc,U0_glb)
            ELSE
                U0_glb = U0_loc
            END IF

            DO j = 1, no_of_itr
                !convert conservative to primitive
                CALL s_convert_conservative_to_primitive_variables(q_cons_ts(1)%vf, ad_prim_vf, gm_alpha,ad_ix,ad_iy,ad_iz)

                !Compute anti-diffusive rhs
                CALL s_compute_anti_diffusive_rhs(q_cons_ts(1)%vf, ad_prim_vf, ad_rhs_vf)

                DO i = 1, adv_idx%end
                    q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                                    q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                                    + dtau*Lheaviside(0:m,0:n,0:p)*U0_glb*ad_rhs_vf(i)%sf(0:m,0:n,0:p)
                END DO
            END DO

            DO i = 1, cont_idx%end
                ad_prim_vf(i)%sf => NULL()
            END DO

            DO i = adv_idx%beg, adv_idx%end
                ad_prim_vf(i)%sf => NULL()
            END DO

        END SUBROUTINE s_compute_anti_diffusion




        !> Compute anti-diffusive rhs terms of Shyue
        !!      Details follow the imcompressible anit-diffusion method of So and Adams
        SUBROUTINE s_compute_anti_diffusive_rhs(ad_cons_vf, ad_prim_vf, ad_rhs_vf)
            

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: ad_cons_vf !<
            !! Cell-average conservative variables


            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: ad_prim_vf !<
            !! Cell-average conservative variables


            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: ad_rhs_vf !<
            !! Cell-average RHS variables


            INTEGER :: h,i,j,k,l !< Generic loop variables

            INTEGER :: ALLOC_ERR !< Error variable

            !Call buffer subroutine to update bcs and interprocess bcs
            CALL s_populate_variables_buffers(ad_cons_vf)

            temp(3,0:m,0:n,0:p) = ad_prim_vf(mom_idx%beg)%sf(0:m,0:n,0:p)
            IF(n > 0) THEN
                temp(4,0:m,0:n,0:p) = ad_prim_vf(mom_idx%beg+1)%sf(0:m,0:n,0:p)
                IF(p > 0) THEN
                    temp(5,0:m,0:n,0:p) = ad_prim_vf(mom_idx%end)%sf(0:m,0:n,0:p)
                ELSE
                    temp(5,0:m,0:n,0:p) = 0d0
                END IF
            ELSE
                temp(4,0:m,0:n,0:p) = 0d0
                temp(5,0:m,0:n,0:p) = 0d0
            END IF
            temp(6,0:m,0:n,0:p) = ad_prim_vf(E_idx)%sf(0:m,0:n,0:p)*fluid_pp(1)%gamma + fluid_pp(1)%pi_inf !rho1 e1
            temp(7,0:m,0:n,0:p) = ad_prim_vf(E_idx)%sf(0:m,0:n,0:p)*fluid_pp(2)%gamma + fluid_pp(2)%pi_inf !rho2 e2

            DO i = 1, num_dims
                IF (i == 1) THEN

                    ALLOCATE(ad_a(ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
                    ALLOCATE(ad_b(ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))
                    ALLOCATE(ad_c(ad_iy%beg:ad_iy%end,ad_iz%beg:ad_iz%end))

                    DO j = 0, m
                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(j+1,:,:) - &
                            ad_cons_vf(cont_idx%beg)%sf( j ,:,:))/(x_cc(j+1) - x_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf( j ,:,:) - &
                            ad_cons_vf(cont_idx%beg)%sf(j-1,:,:))/(x_cc( j ) - x_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho1(j,:,:) = ad_c(:,:) 

                        ad_a(:,:) = (ad_cons_vf(cont_idx%end)%sf(j+1,:,:) - &
                            ad_cons_vf(cont_idx%end)%sf( j ,:,:))/(x_cc(j+1) - x_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%end)%sf( j ,:,:) - &
                            ad_cons_vf(cont_idx%end)%sf(j-1,:,:))/(x_cc( j ) - x_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho2(j,:,:) = ad_c(:,:)

                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(j+1,:,:) + &
                                    ad_cons_vf(cont_idx%end)%sf(j+1,:,:) &
                                    -ad_cons_vf(cont_idx%beg)%sf( j ,:,:) - &
                                    ad_cons_vf(cont_idx%end)%sf( j ,:,:))/(x_cc(j+1) - x_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf( j ,:,:) + &
                                    ad_cons_vf(cont_idx%end)%sf( j ,:,:) &
                                    -ad_cons_vf(cont_idx%beg)%sf(j-1,:,:) - &
                                    ad_cons_vf(cont_idx%end)%sf(j-1,:,:))/(x_cc( j ) - x_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_rho(j,:,:) = ad_c(:,:)

                        ad_a(:,:) = (ad_cons_vf(adv_idx%beg)%sf(j+1,:,:) - &
                                ad_cons_vf(adv_idx%beg)%sf( j ,:,:))/(x_cc(j+1) - x_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(adv_idx%beg)%sf( j ,:,:) - &
                                ad_cons_vf(adv_idx%beg)%sf(j-1,:,:))/(x_cc( j ) - x_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpha(j,:,:) = ad_c(:,:) 
                    END DO

                    !Setting buffer values to zero ( since alpha at buffers is same as interior--need to check if this is right)
                    grad_alpharho1(-buff_size:-1,:,:)  = 0d0
                    grad_alpharho1(m+1:m+buff_size,:,:) = 0d0
                    grad_alpharho2(-buff_size:-1,:,:)  = 0d0
                    grad_alpharho2(m+1:m+buff_size,:,:) = 0d0
                    grad_rho(-buff_size:-1,:,:)  = 0d0
                    grad_rho(m+1:m+buff_size,:,:) = 0d0
                    grad_alpha(-buff_size:-1,:,:)  = 0d0
                    grad_alpha(m+1:m+buff_size,:,:) = 0d0

                    !Compute face flux
                    DO l = 0, p
                        DO k = 0, n
                            DO j = -1, m
                                IF (ABS(grad_alpharho1(j+1,k,l)) > ABS(grad_alpharho1(j,k,l))) THEN
                                    face_flux(1,j,k,l) = grad_alpharho1(j,k,l) 
                                ELSE
                                    face_flux(1,j,k,l) = grad_alpharho1(j+1,k,l) 
                                END IF

                                IF (ABS(grad_alpharho2(j+1,k,l)) > ABS(grad_alpharho2(j,k,l))) THEN
                                    face_flux(2,j,k,l) = grad_alpharho2(j,k,l) 
                                ELSE
                                    face_flux(2,j,k,l) = grad_alpharho2(j+1,k,l) 
                                END IF

                                IF (ABS(grad_rho(j+1,k,l)) > ABS(grad_rho(j,k,l))) THEN
                                    face_flux(3,j,k,l) = grad_rho(j,k,l) 
                                ELSE
                                    face_flux(3,j,k,l) = grad_rho(j+1,k,l) 
                                END IF

                                IF (ABS(grad_alpha(j+1,k,l)) > ABS(grad_alpha(j,k,l))) THEN
                                    face_flux(4,j,k,l) = grad_alpha(j,k,l) 
                                ELSE
                                    face_flux(4,j,k,l) = grad_alpha(j+1,k,l) 
                                END IF
                            END DO
                        END DO
                    END DO

                    DO l = 0, p
                        DO k = 0, n
                            DO j = 0, m
                                ad_rhs_vf(cont_idx%beg)%sf(j,k,l) = -1d0/dx(j)*diff_coeff_x(j,k,l)*(face_flux(1,j,k,l) &
                                    - face_flux(1,j-1,k,l))
                                ad_rhs_vf(cont_idx%end)%sf(j,k,l) = -1d0/dx(j)*diff_coeff_x(j,k,l)*(face_flux(2,j,k,l) &
                                    - face_flux(2,j-1,k,l))
                                DO h = 1, num_dims
                                    ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) = &
                                        -1d0/dx(j)*diff_coeff_x(j,k,l)*ad_prim_vf(cont_idx%end+h)%sf(j,k,l) * &
                                        (face_flux(3,j,k,l) - face_flux(3,j-1,k,l))
                                END DO
                                ad_rhs_vf(E_idx)%sf(j,k,l) = -1d0/dx(j)*diff_coeff_x(j,k,l)*( &
                                                            (5d-1*(temp(3,j,k,l)**2d0 + temp(4,j,k,l)**2d0 + temp(5,j,k,l)**2d0)) * &
                                                            (face_flux(3,j,k,l) - face_flux(3,j-1,k,l)) + &
                                                            (temp(6,j,k,l) - temp(7,j,k,l))*(face_flux(4,j,k,l)-face_flux(4,j-1,k,l)))
                                ad_rhs_vf(adv_idx%beg)%sf(j,k,l) = -1d0/dx(j)*diff_coeff_x(j,k,l)*(face_flux(4,j,k,l)-face_flux(4,j-1,k,l))
                            END DO
                        END DO
                    END DO

                    DEALLOCATE(ad_a,ad_b,ad_c, STAT = ALLOC_ERR)

                ELSEIF (i == 2) THEN

                    ALLOCATE(ad_a(ad_ix%beg:ad_ix%end,ad_iz%beg:ad_iz%end))
                    ALLOCATE(ad_b(ad_ix%beg:ad_ix%end,ad_iz%beg:ad_iz%end))
                    ALLOCATE(ad_c(ad_ix%beg:ad_ix%end,ad_iz%beg:ad_iz%end))

                    DO j = 0, n
                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,j+1,:) - ad_cons_vf(cont_idx%beg)%sf(:, j ,:))/(y_cc(j+1) - y_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:, j ,:) - ad_cons_vf(cont_idx%beg)%sf(:,j-1,:))/(y_cc( j ) - y_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho1(:,j,:) = ad_c(:,:) 

                        ad_a(:,:) = (ad_cons_vf(cont_idx%end)%sf(:,j+1,:) - ad_cons_vf(cont_idx%end)%sf(:, j ,:))/(y_cc(j+1) - y_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%end)%sf(:, j ,:) - ad_cons_vf(cont_idx%end)%sf(:,j-1,:))/(y_cc( j ) - y_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho2(:,j,:) = ad_c(:,:)
          
                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,j+1,:) + ad_cons_vf(cont_idx%end)%sf(:,j+1,:) &
                                    -ad_cons_vf(cont_idx%beg)%sf(:, j ,:) - ad_cons_vf(cont_idx%end)%sf(:, j ,:))/(y_cc(j+1) - y_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:, j ,:) + ad_cons_vf(cont_idx%end)%sf(:, j ,:) &
                                    -ad_cons_vf(cont_idx%beg)%sf(:,j-1,:) - ad_cons_vf(cont_idx%end)%sf(:,j-1,:))/(y_cc( j ) - y_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_rho(:,j,:) = ad_c(:,:)
          
                        ad_a(:,:) = (ad_cons_vf(adv_idx%beg)%sf(:,j+1,:) - ad_cons_vf(adv_idx%beg)%sf(:, j ,:))/(y_cc(j+1) - y_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(adv_idx%beg)%sf(:, j ,:) - ad_cons_vf(adv_idx%beg)%sf(:,j-1,:))/(y_cc( j ) - y_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpha(:,j,:) = ad_c(:,:) 
                    END DO
          
                    !Setting buffer values to zero ( since alpha at buffers is same as interior--need to check if this is right)
                    grad_alpharho1(:,-buff_size:-1,:)  = 0d0
                    grad_alpharho1(:,n+1:n+buff_size,:) = 0d0
                    grad_alpharho2(:,-buff_size:-1,:)  = 0d0
                    grad_alpharho2(:,n+1:n+buff_size,:) = 0d0
                    grad_rho(:,-buff_size:-1,:)  = 0d0
                    grad_rho(:,n+1:n+buff_size,:) = 0d0
                    grad_alpha(:,-buff_size:-1,:)  = 0d0
                    grad_alpha(:,n+1:n+buff_size,:) = 0d0

                    !Compute face flux
                    DO l = 0, p
                        DO k = -1, n
                            DO j = 0, m
                                IF (ABS(grad_alpharho1(j,k+1,l)) > ABS(grad_alpharho1(j,k,l))) THEN
                                    face_flux(1,j,k,l) = grad_alpharho1(j,k,l) 
                                ELSE
                                    face_flux(1,j,k,l) = grad_alpharho1(j,k+1,l) 
                                END IF
             
                                IF (ABS(grad_alpharho2(j,k+1,l)) > ABS(grad_alpharho2(j,k,l))) THEN
                                    face_flux(2,j,k,l) = grad_alpharho2(j,k,l) 
                                ELSE
                                    face_flux(2,j,k,l) = grad_alpharho2(j,k+1,l) 
                                END IF
             
                                IF (ABS(grad_rho(j,k+1,l)) > ABS(grad_rho(j,k,l))) THEN
                                    face_flux(3,j,k,l) = grad_rho(j,k,l) 
                                ELSE
                                    face_flux(3,j,k,l) = grad_rho(j,k+1,l) 
                                END IF
             
                                IF (ABS(grad_alpha(j,k+1,l)) > ABS(grad_alpha(j,k,l))) THEN
                                    face_flux(4,j,k,l) = grad_alpha(j,k,l) 
                                ELSE
                                    face_flux(4,j,k,l) = grad_alpha(j,k+1,l) 
                                END IF
                            END DO
                        END DO
                    END DO
                    
                    DO l = 0, p
                        DO k = 0, n
                            DO j = 0, m
                                ad_rhs_vf(cont_idx%beg)%sf(j,k,l) = &
                                    ad_rhs_vf(cont_idx%beg)%sf(j,k,l)-1d0/dy(k)*diff_coeff_y(j,k,l) * &
                                    (face_flux(1,j,k,l) - face_flux(1,j,k-1,l))
                                ad_rhs_vf(cont_idx%end)%sf(j,k,l) = &
                                    ad_rhs_vf(cont_idx%end)%sf(j,k,l)-1d0/dy(k)*diff_coeff_y(j,k,l) * &
                                    (face_flux(2,j,k,l) - face_flux(2,j,k-1,l))
                                DO h = 1, num_dims
                                    ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) = ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) - &
                                        1d0/dy(k)*diff_coeff_y(j,k,l)*ad_prim_vf(cont_idx%end+h)%sf(j,k,l) * &
                                        (face_flux(3,j,k,l) - face_flux(3,j,k-1,l))
                                END DO
                                ad_rhs_vf(E_idx)%sf(j,k,l) = &
                                    ad_rhs_vf(E_idx)%sf(j,k,l)-1d0/dy(k)*diff_coeff_y(j,k,l)*( &
                                    (5d-1*(temp(3,j,k,l)**2d0 + temp(4,j,k,l)**2d0 + temp(5,j,k,l)**2d0)) * &
                                    (face_flux(3,j,k,l) - face_flux(3,j,k-1,l)) + &
                                    (temp(6,j,k,l) - temp(7,j,k,l))*(face_flux(4,j,k,l)-face_flux(4,j,k-1,l)))
                                ad_rhs_vf(adv_idx%beg)%sf(j,k,l) = &
                                    ad_rhs_vf(adv_idx%beg)%sf(j,k,l) - &
                                    1d0/dy(k)*diff_coeff_y(j,k,l)*(face_flux(4,j,k,l)-face_flux(4,j,k-1,l))
  
                                IF(cyl_coord.AND.p.eq.0) THEN
                                    ad_rhs_vf(cont_idx%beg)%sf(j,k,l) = &
                                        ad_rhs_vf(cont_idx%beg)%sf(j,k,l)-5d-1/y_cc(k)*diff_coeff_y(j,k,l) * &
                                        (face_flux(1,j,k,l) + face_flux(1,j,k-1,l))
                                    ad_rhs_vf(cont_idx%end)%sf(j,k,l) = &
                                        ad_rhs_vf(cont_idx%end)%sf(j,k,l)-5d-1/y_cc(k)*diff_coeff_y(j,k,l) * &
                                        (face_flux(2,j,k,l) + face_flux(2,j,k-1,l))
                                    DO h = 1, num_dims
                                        ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) = &
                                            ad_rhs_vf(cont_idx%end+h)%sf(j,k,l)-5d-1/y_cc(k)*diff_coeff_y(j,k,l) * &
                                            ad_prim_vf(cont_idx%end+h)%sf(j,k,l)*(face_flux(3,j,k,l) + face_flux(3,j,k-1,l))
                                    END DO
                                    ad_rhs_vf(E_idx)%sf(j,k,l) = &
                                        ad_rhs_vf(E_idx)%sf(j,k,l)-5d-1/y_cc(k)*diff_coeff_y(j,k,l)*( &
                                        (5d-1*(temp(3,j,k,l)**2d0 + temp(4,j,k,l)**2d0 + temp(5,j,k,l)**2d0)) * &
                                        (face_flux(3,j,k,l) + face_flux(3,j,k-1,l)) + &
                                        (temp(6,j,k,l) - temp(7,j,k,l))*(face_flux(4,j,k,l)+face_flux(4,j,k-1,l)))
                                    ad_rhs_vf(adv_idx%beg)%sf(j,k,l) = ad_rhs_vf(adv_idx%beg)%sf(j,k,l) - &
                                        5d-1/y_cc(k)*diff_coeff_y(j,k,l)*(face_flux(4,j,k,l)+face_flux(4,j,k-1,l))
                                END IF
                            END DO
                        END DO
                    END DO
          
                    DEALLOCATE(ad_a,ad_b,ad_c, STAT = ALLOC_ERR)

                ELSE

                    ALLOCATE(ad_a(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end))
                    ALLOCATE(ad_b(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end))
                    ALLOCATE(ad_c(ad_ix%beg:ad_ix%end,ad_iy%beg:ad_iy%end))

                    DO j = 0, p
                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,:,j+1) - ad_cons_vf(cont_idx%beg)%sf(:,:, j ))/(z_cc(j+1) - z_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,:, j ) - ad_cons_vf(cont_idx%beg)%sf(:,:,j-1))/(z_cc( j ) - z_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho1(:,:,j) = ad_c(:,:) 

                        ad_a(:,:) = (ad_cons_vf(cont_idx%end)%sf(:,:,j+1) - ad_cons_vf(cont_idx%end)%sf(:,:, j ))/(z_cc(j+1) - z_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%end)%sf(:,:, j ) - ad_cons_vf(cont_idx%end)%sf(:,:,j-1))/(z_cc( j ) - z_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpharho2(:,:,j) = ad_c(:,:)

                        ad_a(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,:,j+1) + ad_cons_vf(cont_idx%end)%sf(:,:,j+1) &
                                    -ad_cons_vf(cont_idx%beg)%sf(:,:, j ) - ad_cons_vf(cont_idx%end)%sf(:,:, j ))/(z_cc(j+1) - z_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(cont_idx%beg)%sf(:,:, j ) + ad_cons_vf(cont_idx%end)%sf(:,:, j ) &
                                    -ad_cons_vf(cont_idx%beg)%sf(:,:,j-1) - ad_cons_vf(cont_idx%end)%sf(:,:,j-1))/(z_cc( j ) - z_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_rho(:,:,j) = ad_c(:,:)

                        ad_a(:,:) = (ad_cons_vf(adv_idx%beg)%sf(:,:,j+1) - ad_cons_vf(adv_idx%beg)%sf(:,:, j ))/(z_cc(j+1) - z_cc( j ))
                        ad_b(:,:) = (ad_cons_vf(adv_idx%beg)%sf(:,:, j ) - ad_cons_vf(adv_idx%beg)%sf(:,:,j-1))/(z_cc( j ) - z_cc(j-1))
                        CALL s_minmod(ad_a,ad_b,ad_c)
                        grad_alpha(:,:,j) = ad_c(:,:) 
                    END DO

                    !Setting buffer values to zero ( since alpha at buffers is same as interior--need to check if this is right)
                    grad_alpharho1(:,:,-buff_size:-1)  = 0d0
                    grad_alpharho1(:,:,p+1:p+buff_size) = 0d0
                    grad_alpharho2(:,:,-buff_size:-1)  = 0d0
                    grad_alpharho2(:,:,p+1:p+buff_size) = 0d0
                    grad_rho(:,:,-buff_size:-1)  = 0d0
                    grad_rho(:,:,p+1:p+buff_size) = 0d0
                    grad_alpha(:,:,-buff_size:-1)  = 0d0
                    grad_alpha(:,:,p+1:p+buff_size) = 0d0

                    !Compute face flux
                    DO l = -1, p
                        DO k = 0, n
                            DO j = 0, m
                                IF (ABS(grad_alpharho1(j,k,l+1)) > ABS(grad_alpharho1(j,k,l))) THEN
                                    face_flux(1,j,k,l) = grad_alpharho1(j,k,l) 
                                ELSE
                                    face_flux(1,j,k,l) = grad_alpharho1(j,k,l+1) 
                                END IF

                                IF (ABS(grad_alpharho2(j,k,l+1)) > ABS(grad_alpharho2(j,k,l))) THEN
                                    face_flux(2,j,k,l) = grad_alpharho2(j,k,l) 
                                ELSE
                                    face_flux(2,j,k,l) = grad_alpharho2(j,k,l+1) 
                                END IF

                                IF (ABS(grad_rho(j,k,l+1)) > ABS(grad_rho(j,k,l))) THEN
                                    face_flux(3,j,k,l) = grad_rho(j,k,l) 
                                ELSE
                                    face_flux(3,j,k,l) = grad_rho(j,k,l+1) 
                                END IF

                                IF (ABS(grad_alpha(j,k,l+1)) > ABS(grad_alpha(j,k,l))) THEN
                                    face_flux(4,j,k,l) = grad_alpha(j,k,l) 
                                ELSE
                                    face_flux(4,j,k,l) = grad_alpha(j,k,l+1) 
                                END IF
                            END DO
                        END DO
                    END DO

                    DO l = 0, p
                        DO k = 0, n
                            DO j = 0, m
                                ad_rhs_vf(cont_idx%beg)%sf(j,k,l) = ad_rhs_vf(cont_idx%beg)%sf(j,k,l) - &
                                    1d0/dz(l)*diff_coeff_z(j,k,l)*(face_flux(1,j,k,l) - face_flux(1,j,k,l-1))
                                ad_rhs_vf(cont_idx%end)%sf(j,k,l) = ad_rhs_vf(cont_idx%end)%sf(j,k,l) - &
                                    1d0/dz(l)*diff_coeff_z(j,k,l)*(face_flux(2,j,k,l) - face_flux(2,j,k,l-1))
                                DO h = 1, num_dims
                                    ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) = ad_rhs_vf(cont_idx%end+h)%sf(j,k,l) - &
                                        1d0/dz(l)*diff_coeff_z(j,k,l)*ad_prim_vf(cont_idx%end+h)%sf(j,k,l) * &
                                        (face_flux(3,j,k,l) - face_flux(3,j,k,l-1))
                                END DO
                                ad_rhs_vf(E_idx)%sf(j,k,l) = ad_rhs_vf(E_idx)%sf(j,k,l) - &
                                    1d0/dz(l)*diff_coeff_z(j,k,l)*( &
                                    (5d-1*(temp(3,j,k,l)**2d0 + temp(4,j,k,l)**2d0 + temp(5,j,k,l)**2d0)) * &
                                    (face_flux(3,j,k,l) - face_flux(3,j,k,l-1)) + &
                                    (temp(6,j,k,l) - temp(7,j,k,l))*(face_flux(4,j,k,l)-face_flux(4,j,k,l-1)))
                                ad_rhs_vf(adv_idx%beg)%sf(j,k,l) = ad_rhs_vf(adv_idx%beg)%sf(j,k,l) - &
                                    1d0/dz(l)*diff_coeff_z(j,k,l)*(face_flux(4,j,k,l)-face_flux(4,j,k,l-1))
                            END DO
                        END DO
                    END DO

                    DEALLOCATE(ad_a,ad_b,ad_c, STAT = ALLOC_ERR)

                END IF
            END DO

        END SUBROUTINE s_compute_anti_diffusive_rhs




        !> Computes minmod of two variables
        !> @param ad_a Variable 1
        !> @param ad_b Variable 2
        !> @param ad_c Minmod result
        SUBROUTINE s_minmod(ad_a,ad_b,ad_c)


            REAL(KIND(0d0)), INTENT(IN) :: ad_a(:,:), ad_b(:,:)
            REAL(KIND(0d0)), INTENT(INOUT) :: ad_c(:,:)

            INTEGER :: array_size(2) ! To find size of 2d arrays a and b
            INTEGER :: i,j

            array_size = shape(ad_a) ! a and b should have the same size

            DO j = 1, array_size(2)
                DO i = 1, array_size(1)
                    IF((ad_a(i,j)*ad_b(i,j) > 0d0) .AND. (ABS(ad_a(i,j)) <= ABS(ad_b(i,j)))) THEN
                        ad_c(i,j) = ad_a(i,j)
                    ELSE IF((ad_a(i,j)*ad_b(i,j) > 0d0) .AND. (ABS(ad_b(i,j)) <= ABS(ad_a(i,j)))) THEN
                        ad_c(i,j) = ad_b(i,j)
                    ELSE
                        ad_c(i,j) = 0d0
                    END IF
                END DO
            END DO

        END SUBROUTINE s_minmod




        !> Finalizes the anti diffusion module via deallocation
        SUBROUTINE s_finalize_anti_diffusion_module()

            INTEGER :: i

            DO i = 1, sys_size
                DEALLOCATE(ad_rhs_vf(i)%sf)
            END DO
            DO i = mom_idx%beg, E_idx
                DEALLOCATE(ad_prim_vf(i)%sf)
            END DO

            DEALLOCATE(ad_rhs_vf, ad_prim_vf)

            DO i = 1, num_fluids
                DEALLOCATE(gm_alpha(i)%sf)
            END DO

            DEALLOCATE(grad_alpha, grad_alpharho1, grad_alpharho2, grad_rho)
            DEALLOCATE(gm_alpha, face_flux, diff_coeff_x, diff_coeff_y, diff_coeff_z, temp)

        END SUBROUTINE s_finalize_anti_diffusion_module





END MODULE m_anti_diffusion
