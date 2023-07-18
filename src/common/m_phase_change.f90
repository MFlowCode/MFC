!>
!! @file m_phase_change.f90
!! @brief Contains module m_phasechange

!> @brief This module is used to relax pressure, temperature, and/or gibbs free energies of the species towards 
!> equilibrium 
module m_phase_change

#ifndef MFC_POST_PROCESS

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    ! ==========================================================================

    IMPLICIT NONE

    PRIVATE; PUBLIC :: s_initialize_phasechange_module, &
                        s_relaxation_solver,             & 
                        s_relaxation_finite_solver,      & 
                        s_finite_ptg_relaxation,         &
                        s_infinite_p_relaxation,         &
                        s_infinite_p_relaxation_k,       &
                        s_infinite_pt_relaxation,        &
                        s_infinite_relaxation_k,         &
                        s_infinite_ptg_relaxation,       &
                        s_finalize_relaxation_solver_module
                        
    !> @name Abstract interface for creating function pointers
    !> @{
    ABSTRACT INTERFACE

        !> @name Abstract subroutine for the finite relaxation solver 
        !> @{
        SUBROUTINE s_abstract_relaxation_solver(q_cons_vf) ! -------
            IMPORT :: scalar_field, sys_size          
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
        END SUBROUTINE

        !> @name Abstract subroutine for the finite relaxation solver 
        !> @{
        SUBROUTINE s_abstract_relaxation_finite_solver(q_cons_vf, rhs_vf) ! -------
            IMPORT :: scalar_field, sys_size
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: rhs_vf
        END SUBROUTINE

    END INTERFACE

    !> @name Parameters for the phase change part of the code
    !> @{
    INTEGER,         PARAMETER :: newton_iter       = 50        !< p_relaxk \alpha iter,                set to 25
    REAL(KIND(0d0)), PARAMETER :: pknewton_eps      = 1.d-15    !< p_relaxk \alpha threshold,           set to 1E-15
    REAL(KIND(0d0)), PARAMETER :: pTsatnewton_eps   = 1.d-10    !< Saturation temperature tol,          set to 1E-10
    REAL(KIND(0d0)), PARAMETER :: ptgnewton_eps     = 1.d-8     !< Saturation p-T-mu tolerance,         set to 1.d-10
    REAL(KIND(0d0)), PARAMETER :: pres_crit         = 22.09D6   !< Critical water pressure              set to 22.06d6
    REAL(KIND(0d0)), PARAMETER :: T_crit            = 647.29D0  !< Critical water temperature           set to 648
    REAL(KIND(0d0)), PARAMETER :: TsatHv            = 1000.d0   !< Saturation temperature threshold,    set to 900
    REAL(KIND(0d0)), PARAMETER :: TsatLv            = 250.d0    !< Saturation temperature threshold,    set to 250
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    REAL(KIND(0d0)) :: n1, n2, pinf1, pinf2
    REAL(KIND(0d0)) :: gibbsA, gibbsB, gibbsC, gibbsD
    !> @}

    PROCEDURE(s_abstract_relaxation_solver), & 
    POINTER :: s_relaxation_solver => NULL()

    PROCEDURE(s_abstract_relaxation_finite_solver), &
    POINTER :: s_relaxation_finite_solver => NULL()


    CONTAINS

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_initialize_phasechange_module()
            n1    = 1.d0/fluid_pp(1)%gamma + 1.d0
            pinf1 = fluid_pp(1)%pi_inf/(1.d0 + fluid_pp(1)%gamma)
            n2    = 1.d0/fluid_pp(2)%gamma + 1.d0
            pinf2 = fluid_pp(2)%pi_inf/(1.d0 + fluid_pp(2)%gamma)
            gibbsA = (n1*fluid_pp(1)%cv - n2*fluid_pp(2)%cv +  & 
                        fluid_pp(2)%qvp - fluid_pp(1)%qvp) / &
                        (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsB = (fluid_pp(1)%qv - fluid_pp(2)%qv) / &
                        (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsC = (n2*fluid_pp(2)%cv - n1*fluid_pp(1)%cv) / &
                        (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsD = (n1*fluid_pp(1)%cv - fluid_pp(1)%cv) / & 
                        (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)

            ! Associating procedural pointer to the subroutine that will be
            ! utilized to calculate the solution of a given Riemann problem
        
            ! Opening and writing the header of the run-time information file
            IF(relax_model == 0 .OR. relax_model == 1) THEN
                s_relaxation_solver => s_infinite_p_relaxation
            ELSEIF (relax_model == 2) THEN
                s_relaxation_solver => s_infinite_pt_relaxation
            ELSEIF (relax_model == 3) THEN
                s_relaxation_solver => s_infinite_ptg_relaxation
            ELSEIF (relax_model == 4) THEN
                s_relaxation_solver => s_infinite_p_relaxation_k
            ELSEIF ( ( relax_model == 5 ) .OR. ( relax_model == 6 ) ) THEN
                s_relaxation_solver => s_infinite_relaxation_k 
            ELSE
                PRINT '(A)', 'relaxation solver was not set!'
                CALL s_mpi_abort()
            END IF

            !IF (relax_model == 1) THEN
            !    s_relaxation_finite_solver => s_finite_ptg_relaxation
            !END IF      

        END SUBROUTINE s_initialize_phasechange_module !-------------------------------

        !> The purpose of this procedure is to employ the inputted
        !!      cell-average conservative variables in order to compute
        !!      the cell-average RHS variables of the semidiscrete form
        !!      of the governing equations by utilizing the appropriate
        !!      Riemann solver.        
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param rhs_vf Cell-average RHS variables
        SUBROUTINE s_finite_ptg_relaxation(q_cons_vf, rhs_vf) ! -------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: rhs_vf

            REAL(KIND(0d0))                                   ::  sum_alpha, Tsat
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::  p_k, T_k, g_k, Z_k
            REAL(KIND(0d0))                                   ::  n_k, pinf_k
            REAL(KIND(0d0))                                   ::  rho_k, rhoeq_k, rhoe
            REAL(KIND(0d0))                                   ::  e_k, phi, psi
            REAL(KIND(0d0))                                   ::  f1, f2, f3, f4
            REAL(KIND(0d0))                                   ::  A1, B1, A2, B2
            REAL(KIND(0d0))                                   ::  rho_I, Q, kappa
            REAL(KIND(0d0))                                   ::  deltap, p_I, mu 
            REAL(KIND(0d0))                                   ::  e_I, mdot, nu, theta
            REAL(KIND(0d0))                                   ::  mdotalpha, mdotrhoe

            INTEGER :: i,j,k,l,r !< Generic loop iterators

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                        END IF

                        IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps .OR. &
                            q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps) THEN
        
                        ! computing n_k, pinf_k, p_k, T_k, and g_k for finite relaxation
                        phi = 0.d0; psi = 0.d0; f1 = 0.d0; f2 = 0.d0; f3 = 0.d0; f4 = 0.d0;
                        !rhoe = 0.d0;

                        DO i = 1, num_fluids
                            n_k    = 1.d0/fluid_pp(i)%gamma + 1.d0
                            pinf_k = fluid_pp(i)%pi_inf/(1.d0 + fluid_pp(i)%gamma)
                            rho_k  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                            rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                            !rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            e_k = q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) &
                                /q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)

                            p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            T_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf/(1.d0+fluid_pp(i)%gamma)) &
                                /(rho_k*fluid_pp(i)%cv)
                            g_k(i) = (n_k*fluid_pp(i)%cv-fluid_pp(i)%qvp)*T_k(i) &
                                -fluid_pp(i)%cv*T_k(i)*log(T_k(i)**(n_k) &
                                /((p_k(i)+pinf_k)**(n_k-1.d0)))+fluid_pp(i)%qv
                            Z_k(i) = n_k*(p_k(i)+pinf_k);

                            phi = phi + 1.d0/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                            psi = psi + 1.d0/(fluid_pp(i)%gamma*q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))

                            f1 = f1 + (p_k(i)+n_k*pinf_k)/q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                            f2 = f2 + pinf_k/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                            f3 = f3 + fluid_pp(i)%qv/(fluid_pp(i)%gamma & 
                                    *q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))
                            f4 = f4 + (e_k-pinf_k/rho_k) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                        END DO
                        !TODO IMPROVE THIS APPROACH
                        Tsat = f_Tsat(p_k(1))
                        !PRINT *,'prelax =',p_k(1),'Trelax = ',T_k(1),', Tsat = ',Tsat

                        IF ( ieee_is_nan(p_k(1)) .OR. p_k(1) < 0.d0 ) THEN 
                            PRINT *, 'code crashed' 
                            CALL s_mpi_abort()
                        END IF

                        kappa = f1/psi
                        rho_I = (phi*f1-psi*f2)/(psi*f4-phi*f3)
                        p_I = (Z_k(2)*p_k(1)+Z_k(1)*p_k(2))/(Z_k(1)+Z_k(2));
                        e_I = f4/phi + f2/(rho_I*phi)

                        !mu = 1.d8
                        mu = 0.d0;
                        theta = 1.d8 
                        nu = 0.d0
                        IF (T_k(1) .GT. Tsat) nu = 1d-3
    
                        deltap = mu*(p_k(1)-p_k(2))
                        Q = theta*(T_k(2)-T_k(1))
                        mdot = nu*(g_k(2)-g_k(1))
                        mdotalpha = mdot/rho_I
                        mdotrhoe = mdot*e_I

                        rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) = & 
                                rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) + deltap + Q/kappa + mdotalpha
                        rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) = & 
                                rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) - deltap - Q/kappa - mdotalpha
                        rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) = & 
                                rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) + mdot
                        rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) = & 
                                rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) - mdot
                        rhs_vf(1+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                rhs_vf(1+internalEnergies_idx%beg-1)%sf(j,k,l) - p_I*deltap + Q + mdotrhoe
                        rhs_vf(2+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                rhs_vf(2+internalEnergies_idx%beg-1)%sf(j,k,l) + p_I*deltap - Q - mdotrhoe
                        END IF
                    END DO
                END DO
            END DO
        END SUBROUTINE s_finite_ptg_relaxation ! --------------------------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_p_relaxation(q_cons_vf) ! ----------------
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            REAL(KIND(0d0))                                   ::      rhoeq_k
            REAL(KIND(0d0))                                   ::           a1
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: p_k, alpha_k
            !> @}
            INTEGER :: i,j,k,l           !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ==================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                    alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                            -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                            /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                            CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        END IF
                    END DO
                END DO
            END DO

        END SUBROUTINE s_infinite_p_relaxation ! ----------------

        !> Description: The purpose of this procedure is to infinitely relax
        !!              the pressures from the internal-energy equations to a
        !!              unique pressure, from which the corresponding volume
        !!              fraction of each phase are recomputed. For conservation
        !!              purpose, this pressure is finally corrected using the
        !!              mixture-total-energy equation.
        SUBROUTINE s_infinite_p_relaxation_k(q_cons_vf) ! ----------------        
            ! Relaxed pressure, initial partial pressures, function f(p) and its partial
            ! derivative df(p), isentropic partial density, sum of volume fractions,
            ! mixture density, dynamic pressure, surface energy, specific heat ratio
            ! function, liquid stiffness function (two variations of the last two
            ! ones), shear and volume Reynolds numbers and the Weber numbers
            ! Cell-average conservative variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf         
            REAL(KIND(0d0)), DIMENSION(num_fluids) ::     rho_K_s, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids) ::   gamma_min, pres_inf
            ! Generic loop iterators
            INTEGER :: i,j,k,l
            ! Relaxation procedure determination variable
            LOGICAL :: relax
            DO i = 1, num_fluids
                gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
                pres_inf(i)  = fluid_pp(i)%pi_inf / (1d0+fluid_pp(i)%gamma)
            END DO
            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        ! Pressures relaxation procedure ===================================
                        relax = .TRUE.
                        DO i = 1, num_fluids
                            IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. (1.d0-palpha_eps)) relax = .FALSE.
                        END DO
                        IF (relax) THEN
                            ! Calculating the initial pressure
                            DO i = 1, num_fluids
                                pres_K_init(i) = ((q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                    -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                    -fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                                IF (pres_K_init(i) .LE. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) & 
                                    pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-0
                            END DO
                            CALL s_compute_p_relax_k(rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                            ! Cell update of the volume fraction
                            DO i = 1, num_fluids
                                IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps) & 
                                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                            END DO
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_p_relaxation_k ! -----------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_pt_relaxation(q_cons_vf) ! ----------------
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: p_k, alpha_k 
            REAL(KIND(0d0))                                   :: rhoeq_k, rhoe, a1
            REAL(KIND(0d0))                                   :: rhoalpha1, rhoalpha2
            !> @}
            INTEGER :: i, j, k, l        !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            !< Computing the constant saturation properties 

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ==================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                    alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                            -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                            /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )

                        ! PT RELAXATION ==================================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            END DO
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_pt_relaxation ! -----------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_ptg_relaxation(q_cons_vf) ! ----------------
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0))                                   :: pres_relax, Trelax
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: p_k, alpha_k, Tk
            REAL(KIND(0d0))                                   :: rhoalpha1, rhoalpha2
            REAL(KIND(0d0))                                   :: rho, rhoe, rhoeq_k
            REAL(KIND(0d0))                                   :: rho1, rho2
            REAL(KIND(0d0))                                   :: a1, a2
            REAL(KIND(0d0))                                   :: gamma, pi_inf, p_infk
            REAL(KIND(0d0))                                   :: qv
            REAL(KIND(0d0))                                   :: Tsat

            !> @}
            INTEGER :: i, j, k, l        !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            !< Computing the constant saturation properties 

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ========================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, j, k, l, rho, &
                                                            gamma, pi_inf, qv )

                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                    alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                            -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                            /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                    p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! PT RELAXATION ==========================================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, j, k, l, rho, &
                                                            gamma, pi_inf, qv )
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            END DO
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! CHECKING IF PTG RELAXATION IS NEEDED  =====================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, j, k, l, rho, &
                                                                gamma, pi_inf, qv )
                        IF ((q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. ptgalpha_eps ) .AND. &
                                q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-ptgalpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) 
                            END DO                   
                            pres_relax = (rhoe - pi_inf)/gamma
                            Tsat = f_Tsat(pres_relax)
                            DO i = 1, num_fluids
                                Tk(i) = ((q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                    -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                    -fluid_pp(i)%pi_inf & 
                                    /(1.d0+fluid_pp(i)%gamma)) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)) 
                            END DO
                            IF (Tk(1) .LT. Tsat) relax = .FALSE.
                        END IF
                        ! PTG RELAXATION PROCEDURE ===========================
                        IF (relax) THEN
                            CALL s_compute_ptg_pTrelax(pres_relax,Trelax,rho,rhoe)
                            p_infk = fluid_pp(1)%pi_inf/(1.d0+fluid_pp(1)%gamma)
                            rho1 = (pres_relax + p_infk)*fluid_pp(1)%gamma /& 
                                    (fluid_pp(1)%cv*Trelax)
                            p_infk = fluid_pp(2)%pi_inf/(1.d0+fluid_pp(2)%gamma)
                            rho2 = (pres_relax + p_infk)*fluid_pp(2)%gamma /& 
                                    (fluid_pp(2)%cv*Trelax)
                            ! Calculate vapor and liquid volume fractions
                            a1 = (rho-rho2)/(rho1-rho2)
                            a2 = 1.d0 - a1
                            ! Cell update of the volume fraction
                            q_cons_vf(cont_idx%beg)%sf(j,k,l)   = rho1*a1
                            q_cons_vf(1+cont_idx%beg)%sf(j,k,l) = rho2*a2
                            q_cons_vf(adv_idx%beg)%sf(j,k,l)    = a1
                            q_cons_vf(1+adv_idx%beg)%sf(j,k,l)  = a2
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_ptg_relaxation ! -----------------------

        ! This subroutine is created to activate either the pT- or the pTg-equilibrium
        ! model, with changes such that mass depletion is taken into consideration
        SUBROUTINE s_infinite_relaxation_k(q_cons_vf) ! ----------------
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)  :: q_cons_vf 
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: g_min, p_inf, p_infA, pDumb, cv, qv, qvp
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: pk, sk, hk, gk
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: ek, rhok
            REAL(KIND(0.0D0))                                       :: A, B, C, D
            REAL(KIND(0.0D0))                                       :: pS, pSOV, pSSL
            REAL(KIND(0.0D0))                                       :: TS, TSOV, TSatOV, TSatSL, TSSL
            REAL(KIND(0.0D0))                                       :: rhoe, Ewe, mCP, mQ, TvF
            REAL(KIND(0.0D0))                                       :: rho, rM, m1, m2, MCT, mixM, rhoT, rMT 
            REAL(KIND(0.0D0))                                       :: TvFT, rhos
            REAL(KIND(0.0D0))                                       :: dynE, VFT
            
            !< Generic loop iterators
            INTEGER :: i, j, k, l, nf
            INTEGER :: lp, vp, ns

            ! indices for the reacting liquid and gas
            lp = 1 
            vp = 2

            ! renaming fluid properties JUST FOR EASINESS
            cv( 1 : num_fluids ) = fluid_pp( 1 : num_fluids )%cv
            qv( 1 : num_fluids ) = fluid_pp( 1 : num_fluids )%qv
            qvp( 1 : num_fluids ) = fluid_pp( 1 : num_fluids )%qvp
            
            ! gamma_min
            g_min( 1 : num_fluids ) = 1.0D0 / fluid_pp( 1 : num_fluids )%gamma + 1.0D0
            
            ! p_inf
            p_inf( 1 : num_fluids )  = fluid_pp( 1 : num_fluids )%pi_inf &
            / ( 1.0D0 + fluid_pp( 1 : num_fluids )%gamma )
            
            ! variables used in the calculation of the saturation curves for fluids 1 and 2
            A = ( g_min( lp ) * cv( lp )  - g_min( vp ) * cv( vp ) &
            + qv( vp ) - qv( lp ) ) / ( ( g_min( vp ) - 1.0D0 ) * cv( vp ) )
            
            B = ( qv( lp ) - qv( vp ) ) / ( ( g_min( vp ) - 1.0D0 ) * cv( vp ) )
            
            C = ( g_min( vp ) * cv( vp ) - g_min( lp ) * cv( lp ) ) & 
            / ( ( g_min( vp ) - 1.0D0 ) * cv( vp ) )
            
            D = ( ( g_min( lp ) - 1.0D0 ) * cv( lp ) ) &
            / ( ( g_min( vp ) - 1.0D0 ) * cv( vp ) )

            ! threshold for 'mixture cell'. If Y < mixM, the cell is not considered a mixture for pTg purposes
            mixM = 1.0D-08

            ! starting equilibrium solver
            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p

                        rho = 0.0D0
                        rhos = 0.0D0
                        TvF = 0.0D0
                        DO i = 1, num_fluids
                        
                            ! pressure that comes from the homogeneous solver at that cell
                            pk( i ) =   ( g_min( i ) - 1 ) &
                            * ( q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) & 
                            -   q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * qv( i ) ) &
                            /   q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )              &
                            - g_min( i ) * p_inf ( i )
                            
                            ! Mixture density
                            rho = rho + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )
                            
                            ! total entropy
                            rhos = rhos + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * sk( i )
                            
                            ! Total Volume Fraction
                            TvF = TvF + q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                        END DO
                        
                        ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                        ! throughout the phase-change process.
                        m1 = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

                        m2 = q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

                        rM  = m1 + m2

                        ! correcting negative volume and mass fraction values in case they happen
                        ! IF ( mpp_lim ) THEN

                        CALL s_correct_partial_densities(g_min, MCT, mixM, pK, p_inf, q_cons_vf, 0.0D0, rM, j, k, l, lp, vp )

                        ! END IF

                        ! kinetic energy so as to calculate the total internal energy as the total energy minus the 
                        ! kinetic energy.rMT
                        dynE = 0.0D0
                        DO i = mom_idx%beg, mom_idx%end

                            dynE = dynE + 5.0D-1 * q_cons_vf( i )%sf( j, k, l ) * q_cons_vf( i )%sf( j, k, l ) / rho

                        END DO

                        ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                        ! at each of the cells. Note I calculate UE as TE - KE due to numerical reasons (stability at the dicontinuities)
                        rhoe = q_cons_vf( E_idx )%sf( j, k, l ) - dynE
                        
                        ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                        CALL s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'pT', pk, pS, p_inf &
                        , p_infA, qv, q_cons_vf, rho, rhoe, TS )
                        
                        ! check if pTg-equilibrium is required
                        ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                        ! given the outputs from the pT- and pTg-equilibrium solvers are just p, T (byproduct)
                        ! , and the partial masses (pTg- case)
                        IF ( ( relax_model .EQ. 6 ) .AND. ( ( q_cons_vf( lp + cont_idx%beg - 1  )%sf( j, k, l ) &
                        .GT. MixM * rM ) .OR. ( q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l )        &
                        .GT. MixM * rM ) ) .AND. ( pS .LT. pres_crit ) .AND. ( TS .LT. T_crit ) ) THEN

                            ! Checking if phase change can happen.
                            ! overheated vapor
                            ! correcting the liquid partial density
                            q_cons_vf( lp + cont_idx%beg - 1  )%sf( j, k, l ) = mixM * rM

                            ! correcting the vapor partial density
                            q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = ( 1.0D0 - mixM ) * rM

                            ! calling pT-equilibrium for overheated vapor
                            CALL s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'OV', pk, pSOV, p_inf &
                            , pDumb, qv, q_cons_vf, rho, rhoe, TSOV )

                            ! calculating Saturation temperature
                            CALL s_TSat(A, B, C, cv, D, g_min, lp, pSOV, p_inf, qv, qvp, TSatOV, TSOV, vp)
                                 
                            ! subcooled liquid
                            ! correcting the liquid partial density
                            q_cons_vf( lp + cont_idx%beg - 1  )%sf( j, k, l ) = ( 1.0D0 - mixM ) * rM

                            ! correcting the vapor partial density
                            q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = mixM * rM

                            ! calling pT-equilibrium for subcooled liquid
                            CALL s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'SL', pk, pSSL, p_inf &
                            , pDumb, qv, q_cons_vf, rho, rhoe, TSSL )
                            
                            ! calculating Saturation temperature
                            CALL s_TSat(A, B, C, cv, D, g_min, lp, pSSL, p_inf, qv, qvp, TSatSL, TSSL, vp)
                            
                            ! Maybe I will have to change this for .GE. instead, since mass depletion can 
                            ! happen at pTg-eq. Make sure it works, first
                            IF ( TSOV .GT. TSatOV ) THEN

                                ! Assigning pressure
                                pS = pSOV

                                ! Assigning Temperature
                                TS = TSOV

                                ! correcting the liquid partial density
                                q_cons_vf( lp + cont_idx%beg - 1  )%sf( j, k, l ) = mixM * rM

                                ! correcting the vapor partial density
                                q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = ( 1.0D0 - mixM ) * rM

                            ELSE IF ( TSSL .LT. TSatSL ) THEN

                                ! Assigning pressure
                                pS = pSSL
                                
                                ! Assigning Temperature
                                TS = TSSL

                                ! correcting the liquid partial density
                                q_cons_vf( lp + cont_idx%beg - 1  )%sf( j, k, l ) = ( 1.0D0 - mixM ) * rM

                                ! correcting the vapor partial density
                                q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = mixM * rM

                            ELSE
                                
                                ! returning partial pressures to what they were from the homogeneous solver
                                ! liquid
                                q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) = m1
                                
                                ! vapor
                                q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) = m2

                                ! calling the pTg-equilibrium solver
                                CALL s_infinite_ptg_relaxation_k( A, B, C, cv, D, g_min, j, k, l, lp, pS &
                                , p_inf, p_infA, rhoe, qv, qvp, q_cons_vf, TS, vp )

                            END IF
                            
                        END IF
                        
                        ! Calculations AFTER the solvers
                        ! A priori, the following variables must be updated: alpha, alpha*rho, and alpha*rho*e now, given
                        ! alpha*rho remains constant in the pT-quilibrium, or the reacting masses have already been updated
                        ! in the pTg-equilibrium, (the remaining) alpha*rho does not need to be updated. I do need, however
                        ! to update alpha, and alpha*rho*e. So, in the end I update alpha, e, and then alpha*rho*e
                        ! entropy
                        sk( 1 : num_fluids ) = &
                        cv( 1 : num_fluids ) * DLOG( ( TS ** g_min( 1 : num_fluids ) ) &
                        / ( ( pS + p_inf( 1 : num_fluids ) ) ** ( g_min( 1 : num_fluids ) - 1.0D0 ) ) ) & 
                        + qvp( 1 : num_fluids )

                        ! enthalpy
                        hk( 1 : num_fluids ) = g_min( 1 : num_fluids ) * cv( 1 : num_fluids ) * TS &
                        + qv( 1 : num_fluids )

                        ! GIBBS-FREE ENERGY
                        gk( 1 : num_fluids ) = hk( 1 : num_fluids ) - TS * sk( 1 : num_fluids )

                        ! densities
                        rhok( 1 : num_fluids ) = ( pS + p_inf( 1 : num_fluids ) ) &
                        / ( ( g_min( 1 : num_fluids ) - 1 ) * cv( 1 : num_fluids ) * TS )

                        ! internal energy
                        ek( 1 : num_fluids ) = ( pS + g_min( 1 : num_fluids ) & 
                        * p_inf( 1 : num_fluids ) ) / ( pS + p_inf( 1 : num_fluids ) ) & 
                        * cv( 1 : num_fluids ) * TS + qv( 1 : num_fluids )

                        DO i = 1, num_fluids
                                                    
                            ! volume fractions
                            q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = &
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) / rhok( i ) 

                            ! alpha*rho*e
                            q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) = & 
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * ek( i )

                        END DO
                        
                        ! calculating the TOTAL VOLUME FRACTION, reacting mass, and TOTAL MASS, after the PC
                        TvFT = 0.0D0
                        rhoT = 0.0D0
                        DO i = 1, num_fluids
                            
                            ! Total volume Fraction Test
                            TvFT = TvFT + q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )
                            
                            ! Total Mass
                            rhoT = rhoT + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                        END DO

                        rMT = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) & 
                        + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

                        ! testing the total volume fraction
                        IF( DABS( TvF - TvFT ) .GT. 1.0D-2 ) THEN
                            
                            PRINT *, 'D total volume fractions AFTER PC: ', TvF - TvFT &
                            , 'j, k, l ', j, k, l

                            PRINT *, 'old', TvF

                            PRINT *, 'new', TvFT

                        END IF

                        ! testing the reacting mass
                        IF ( DABS( rM - rMT ) .GT. 1.0D-8 ) THEN

                            PRINT *, 'D Reacting Mass AFTER PC: ', rM - rMT &
                            , 'j, k, l ', j, k, l

                        END IF

                        ! testing the total mass
                        IF ( DABS( rho - rhoT ) .GT. 1.0D-6 ) THEN

                            PRINT *, 'D Total Mass AFTER PC: ', rho - rhoT &
                            , 'j, k, l ', j, k, l

                        END IF

                    END DO
                END DO
            END DO

        END SUBROUTINE s_infinite_relaxation_k ! ----------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! SUBROUTINES SUBROUTINES SUBROUTINES !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE s_infinite_pt_relaxation_k(cv, j, k, l, g_min, MDF, pk, pS, p_inf, p_infA, qv, q_cons_vf, rho, rhoe, TS )
            ! initializing variables
            TYPE( scalar_field ),  DIMENSION( sys_size ), INTENT( IN )    :: q_cons_vf         
            REAL( KIND( 0.0D0 ) ), DIMENSION( num_fluids ), INTENT( IN )  :: g_min, p_inf, pk, cv, qv
            REAL( KIND( 0.0D0 ) ), INTENT( OUT )                          :: pS, TS
            REAL( KIND( 0.0D0 ) ), DIMENSION( num_fluids ), INTENT( OUT ) :: p_infA
            REAL( KIND( 0.0D0 ) ), INTENT( IN )                           :: rho, rhoe
            INTEGER, INTENT( IN )                                         :: j, k, l
            INTEGER, DIMENSION( num_fluids )                              :: ig

            REAL( KIND( 0.0D0 ) )                                         :: gp, gpp, hp, pO, mCP, mQ
            CHARACTER( LEN = 2 )                                          :: MDF
            INTEGER                                                       :: i, ns

            ! auxiliary variables for the pT-equilibrium solver
            mCP = 0.0D0
            mQ = 0.0D0
            p_infA = p_inf

            ig = 0

            ! Performing tests before initializing the pT-equilibrium
            DO i = 1, num_fluids

                ! PRINT *, 'i', i

                ! PRINT *, q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                ! check if all alpha(i)*rho(i) are nonnegative. If so, abort
                IF ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN
                    
                    CALL s_tattletale( (/ 0.0D0, 0.0D0 /), RESHAPE( (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), (/ 2, 2 /) ) &
                    , j, (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), k, l, mQ, p_inf, pS, qv, (/ DABS( pS - pO ), DABS( pS - pO ) /) &
                    , rhoe, q_cons_vf, TS )
        
                    CALL s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                    . Please, check partial densities. Aborting!')
                
                ! check which indices I will ignore (no need to abort the solver in this case). Adjust this sgm_eps value for mixture cells
                ! ELSE IF ( ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .GE. 0.0D0 ) &
                !     .AND. ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LE. sgm_eps ) ) THEN

                !     ig( i ) = i

                !     ! this value is rather arbitrary, as I am interested in MINVAL( p_inf for the solver ). This way, I am ensuring
                !     ! this value will not be selected.
                !     p_infA( i ) = 2 * MAXVAL( p_inf )

                END IF
                
                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                * cv( i ) * g_min( i )
        
                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * qv( i )
        
            END DO
        
            ! Checking energy constraint
            IF ( ( rhoe - mQ - MINVAL( p_infA ) ) .LT. 0.0D0 ) THEN

                ! PRINT *, 'ERROR IN THE ENERGY CONSTRAINT'

                IF ( ( MDF .EQ. 'OV' ) .OR. ( MDF .EQ. 'SL' ) ) THEN

                    ! PRINT *, 'MDF', MDF

                    ! Assigning zero values for mass depletion cases
                    ! pressure 
                    pS = 0.0D0
                    
                    ! temperature
                    TS = 0.0D0

                    RETURN

                ELSE
            
                    PRINT *, 'MDF', MDF

                    CALL s_tattletale( (/ 0.0D0, 0.0D0 /), RESHAPE( (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), (/ 2, 2 /) ) &
                    , j, (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), k, l, mQ, p_inf, pS, qv, (/ DABS( pS - pO ), DABS( pS - pO ) /) &
                    , rhoe, q_cons_vf, TS )

                    CALL s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                    . Please, check energy constraint. Aborting!')

                END IF
        
            END IF
        
            ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to 
            ! iterate over the Newton's solver
            pO = 0.0D0
        
            ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(p_inf)
            ! and infinity, a solution should be able to be found.
            pS = MAX( MAXVAL( pk ), -1.0D0 * MINVAL( p_infA ) ) + 1.0D4
        
            ! Newton Solver for the pT-equilibrium
            ns = 0
            DO WHILE ( ( DABS( pS - pO ) .GT. palpha_eps ) .AND. ( DABS( ( pS - pO ) / pO ) .GT. palpha_eps ) &
                        .OR. ( ns .EQ. 0 ) )

                ! increasing counter
                ns = ns + 1
        
                ! updating old pressure 
                pO = pS
        
                ! updating functions used in the Newton's solver
                gpp = 0.0D0
                gp = 0.0D0
                hp = 0.0D0

                DO i = 1, num_fluids

                    ! given pS always change, I need ig( i ) and gp to be in here, as it dynamically updates.
                    ! not that I do not need to use p_infA here, but I will do it for consistency
                    IF ( i .NE. ig( i ) ) THEN

                        gp =   gp + ( g_min( i ) - 1.0D0 ) & 
                        * q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * cv( i ) &
                        * ( rhoe + pS - mQ ) / ( mCP * ( pS + p_infA( i ) ) )
            
                        gpp = gpp + ( g_min( i ) - 1.0D0 ) & 
                        * q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * cv( i ) &
                        * ( p_infA( i ) - rhoe + mQ ) / ( mCP * ( pS + p_infA( i ) ) ** 2 )

                    END IF

                END DO
        
                hp = 1.0D0 / ( rhoe + pS - mQ ) + 1.0D0 / ( pS + MINVAL( p_infA ) )
        
                ! updating common pressure for the newton solver
                pS = pO + ( ( 1.0D0 - gp ) / gpp ) / ( 1.0D0 - ( 1.0D0 - gp + DABS( 1.0D0 - gp ) ) &
                / ( 2.0D0 * gpp ) * hp )
        
                ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
                ! I am not an expert in numerical methods though)
                IF ( pS .LE. -1.0D0 * MINVAL( p_infA ) ) THEN

                    CALL s_tattletale( (/ 0.0D0, 0.0D0 /), RESHAPE( (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), (/ 2, 2 /) ) &
                    , j, (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /), k, l, mQ, p_infA, pS, qv, (/ pS - pO, pS + pO /) &
                    , rhoe, q_cons_vf, TS )
    
                    CALL s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                    . Please, check the pressure value. Aborting!')
                
                ELSEIF ( ISNAN( pS ) ) THEN

                    CALL s_mpi_abort( 'Solver for the pT-relaxation returned NaN values     &
                    (m_phase_change, s_infinite_pt_relaxation_k). Please, check the error.  &
                    Aborting!' )

                ! cheching is the maximum number of iterations was reached
                ELSEIF ( ns .GT. 1E4 ) THEN
        
                    CALL s_mpi_abort('Maximum number of iterations for the pT-relaxation solver reached &
                    (m_phase_change, s_infinite_pt_relaxation_k). Please, check the pressure value. Aborting!')

                END IF
        
            END DO
        
            ! common temperature
            TS = ( rhoe + pS - mQ ) / mCP
        
        END SUBROUTINE s_infinite_pt_relaxation_k ! -----------------------

        SUBROUTINE s_infinite_ptg_relaxation_k( A, B, C, cv, D, g_min, j, k, l, lp, pS, p_inf, p_infpT, rhoe, qv, qvp, q_cons_vf, TS, vp )

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)        :: q_cons_vf 
            REAL( KIND( 0.0D0 ) ), INTENT( INOUT )                        :: pS, TS
            REAL( KIND( 0.0D0 ) ), DIMENSION( num_fluids ), INTENT( IN )  :: g_min, p_inf, p_infpT, cv, qv, qvp
            REAL( KIND( 0.0D0 ) ), INTENT( IN )                           :: A, B, C, D, rhoe
            INTEGER, INTENT(IN)                                           :: j, k, l, lp, vp
            
            REAL( KIND( 0.0D0 ) ), DIMENSION( num_fluids )                :: p_infA
            REAL(KIND(0.0D0)), DIMENSION(2,2)                             :: Jac, InvJac, TJac
            REAL(KIND(0.0D0)), DIMENSION(2)                               :: R2D, DeltamP    
            REAL(KIND(0.0D0))                                             :: Om
            REAL(KIND(0.0D0))                                             :: mCP, mCPD, mQ, mQD, mCVGP, mCVGP2
            
            !< Generic loop iterators
            INTEGER :: i, nf, ns
        
            ! pTg-equilibrium solution procedure
            ! Newton Sover parameters
            ! counter
            ns = 0
        
            ! Relaxation factor
            Om = 1.0D-4
        
            ! Dummy guess to start the pTg-equilibrium problem.
            ! improve this initial condition
            R2D( 1 ) = 0.0D0
            R2D( 2 ) = 0.0D0
            DeltamP( 1 ) = 0.0D0
            DeltamP( 2 ) = 0.0D0
            mQ = 0.0D0

            p_infA = p_infpT

            IF ( ( ( pS .LT. 0.0D0 ) .AND. ( ( q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) ) .GT. ( ( rhoe                  & 
            - g_min( lp ) * p_inf( lp ) / ( g_min( lp ) - 1  ) ) / qv( i ) ) ) ) .OR.  &
            ( ( pS .GE. 0.0D0 ) .AND. ( pS .LT. 1.0D-1 ) ) ) THEN

                ! improve this initial condition
                pS = 1.0D4

            END IF

            ! Loop until the solution for F(X) is satisfied
            ! Check whether I need to use both absolute and relative values
            ! for the residual, and how to do it adequately.
            DO WHILE ( ( ( DSQRT( R2D( 1 ) ** 2 + R2D( 2 ) ** 2 ) .GT. ptgalpha_eps )                       &
                .AND.  ( ( DSQRT( R2D( 1 ) ** 2 + R2D( 2 ) ** 2 ) / rhoe ) .GT. ( ptgalpha_eps / 1E4 ) ) )  &
                .OR. ( ns .EQ. 0 ) )
                
                ! Updating counter for the iterative procedure
                ns = ns + 1

                ! Auxiliary variable to help in the calculation of the residue
                mCP = 0.0D0
                mCPD = 0.0D0
                mCVGP = 0.0D0
                mCVGP2 = 0.0D0
                mQ = 0.0D0
                mQD = 0.0D0
                
                ! Those must be updated through the iterations, as they either depend on 
                ! the partial masses for all fluids, or on the equilibrium pressure
                DO i = 1, num_fluids
        
                    ! sum of the total alpha*rho*cp of the system
                    mCP = mCP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                    * cv( i ) * g_min( i )

                    ! sum of the total alpha*rho*q of the system
                    mQ = mQ + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * qv( i )
        
                    ! These auxiliary variables now need to be updated, as the partial densities now 
                    ! vary at every iteration
                    IF ( ( i .NE. lp ) .AND. ( i .NE. vp ) ) THEN
        
                        mCVGP = mCVGP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                        * cv( i ) * ( g_min( i ) - 1 ) / ( pS + p_inf( i ) )
        
                        mCVGP2 = mCVGP2 + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                        * cv( i ) * ( g_min( i ) - 1 ) / ( ( pS + p_inf( i ) ) ** 2 )
        
                        mQD = mQD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * qv( i )
        
                        ! sum of the total alpha*rho*cp of the system
                        mCPD = mCPD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * cv( i ) &
                        * g_min( i )
        
                    END IF

                    ! IF ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .LT. 1.0D-10 ) THEN
    
                    !     p_infA( i ) = 2 * MAXVAL( p_inf )
    
                    ! ELSE

                    !     p_infA( i ) = p_inf( i )

                    ! END IF
        
                END DO

                ! Checking pressure and energy criteria for the (pT) solver to find a solution
                IF ( ( pS .LE. -1.0D0 * MINVAL( p_inf ) ) &
                .OR. ( ( rhoe - mQ - MINVAL( p_inf ) ) .LT. 0.0D0 ) ) THEN
                    
                    CALL s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_inf, pS, qv &
                                    , R2D, rhoe, q_cons_vf, TS )
        
                    CALL s_mpi_abort( 'Solver for the pTg-relaxation failed &
                    (m_phase_change, s_infinite_ptg_relaxation_k). Either pS is out of bounds. &
                    or the energy constraint has been violated. Aborting!' )
        
                END IF
                        
                ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
                CALL s_compute_jacobian_matrix( B, C, cv, D, rhoe, g_min, InvJac, j, Jac, k, l, &
                lp, mCPD, mCVGP, mCVGP2, mQD, p_inf, pS, qv, qvp, q_cons_vf, TJac, vp )
        
                ! calculating correction array for Newton's method
                DeltamP = - 1.0D0 * MATMUL( InvJac, R2D )

                ! checking if the correction in the mass/pressure will lead to negative values for those quantities
                ! If so, adjust the underrelaxation parameter Om
                IF ( q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) - Om * DeltamP( 1 ) .LT. 0 ) THEN
                    CALL s_mpi_abort('crapV')
                END IF

                IF ( q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) + Om * DeltamP( 1 ) .LT. 0 ) THEN
                    CALL s_mpi_abort('crapL')
                END IF

                ! updating two reacting 'masses'. Recall that the other 'masses' do not change during the phase change
                ! liquid
                q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) = &
                q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) + Om * DeltamP( 1 )
        
                ! gas
                q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = &
                q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) - Om * DeltamP( 1 )
        
                ! updating pressure
                pS = pS + Om * DeltamP( 2 ) 
        
                ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
                ! and (ii) the energy before and after the phase-change process.
                CALL s_compute_pTg_residue(A, B, C, cv, D, g_min, j, k, l, lp, mCPD, mCVGP, mQD &
                , p_inf, qv, qvp, q_cons_vf, pS, rhoe, R2D, vp)

                ! checking if the residue returned any NaN values
                IF ( ( ISNAN( R2D( 1 ) ) ) .OR. ( ISNAN( R2D( 2 ) ) ) ) THEN

                    CALL s_mpi_abort( 'Solver for the pTg-relaxation returned NaN values &
                    (m_phase_change, s_infinite_ptg_relaxation_k). Please, check the error. &
                    Aborting!' )

                ! checking if maximum number of iterations was reached
                ELSEIF ( ( ns .GT. 1E6 ) ) THEN

                    CALL s_mpi_abort( 'Maximum number of iterations reached for the (m_phase_change &
                    , s_infinite_ptg_relaxation_k). Aborting!' )

                END IF

            END DO

            ! common temperature
            TS = ( rhoe + pS - mQ ) / mCP

        END SUBROUTINE s_infinite_ptg_relaxation_k ! -----------------------

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial den/sity, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        SUBROUTINE s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            INTEGER, INTENT(IN)                                    :: j, k, l
            REAL(KIND(0d0))                                        :: sum_alpha
            !> @}
            INTEGER :: i           !< Generic loop iterators

            sum_alpha = 0.0D0
            
            DO i = 1, num_fluids

                ! I changed this line to .LE. so the phase cannot 'disappear'.
                ! Think of more appropriate conditions, later on
                IF ( ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LT. sgm_eps ) .OR. &
                        ( q_cons_vf( i + adv_idx%beg  - 1 )%sf( j, k, l ) .LT. sgm_eps ) ) THEN

                    q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) = sgm_eps
                    
                    q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )  = sgm_eps
                    
                    q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf(j,k,l) = 0.0D0

                END IF

                IF ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .GT. 1.0D0 ) THEN 
                
                    q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = 1.0D0

                END IF
                
                sum_alpha = sum_alpha + q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)

                END DO

                DO i = 1, num_fluids
                    
                q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = &                 
                q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) / sum_alpha

                END DO

        END SUBROUTINE s_mixture_volume_fraction_correction

        SUBROUTINE s_correct_partial_densities(g_min, MCT, mixM, pS, p_inf, q_cons_vf, rhoe, rM, j, k, l, lp, vp )
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)      :: q_cons_vf 
            REAL( KIND( 0.0D0 ) ), INTENT(OUT)                          :: MCT
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: mixM, rhoe, rM
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)    :: pS, g_min, p_inf
            INTEGER, INTENT(IN)                                         :: j, k, l, lp, vp
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids)                :: rho
            REAL( KIND( 0.0D0 ) )                                       :: TVF, mCP, mQ
            !> @}
            INTEGER :: i, AR           !< Generic loop iterators

            IF ( rM .LT. 0.0D0 ) THEN
                CALL s_mpi_abort( 'total reacting mass is negative on "s_correct_partial_densities". Aborting!')
            END IF

            ! Defining the correction in terms of an absolute value might not be the best practice.
            ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
            MCT = 2 * mixM

            ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
            IF ( q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN
                
                ! PRINT *, 'lp deficient'

                ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

                ! PRINT *, q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

                q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) = MCT * rM
                
                q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) = ( 1.0D0 - MCT ) * rM

            ELSE IF ( q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN

                ! PRINT *, 'vp deficient'

                ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

                ! PRINT *, q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

                q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) = ( 1.0D0 - MCT ) * rM
                
                q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) = MCT * rM

            END IF

        END SUBROUTINE s_correct_partial_densities

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        SUBROUTINE s_mixture_total_energy_correction(q_cons_vf, j, k, l )
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)  :: q_cons_vf 
            INTEGER, INTENT(IN)                                     :: j, k, l
            REAL(KIND(0d0))                                         :: rho, dyn_pres, E_We
            REAL(KIND(0d0))                                         :: gamma, pi_inf, pres_relax
            REAL(KIND(0d0))                                         :: qv
            REAL(KIND(0d0))                                         :: prT, rhoeT

            !> @}
            INTEGER :: i     !< Generic loop iterators
            CALL s_convert_to_mixture_variables( q_cons_vf, j, k, l, rho, &
                                                    gamma, pi_inf, qv )
            dyn_pres = 0.d0
            DO i = mom_idx%beg, mom_idx%end
                ! PRINT *, 'i from s_mixture_tec', i
                dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j,k,l) * & 
                q_cons_vf(i)%sf(j,k,l) / MAX(rho,sgm_eps)
            END DO
            E_We = 0.d0
            pres_relax = (q_cons_vf(E_idx)%sf(j,k,l) - dyn_pres - pi_inf - E_We)/gamma
            
            rhoeT = 0.0D0
            DO i = 1, num_fluids
                rhoeT = rhoeT + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
            END DO
            
            pRT = (rhoeT - pi_inf) / gamma

            ! PRINT *, 'P_RELAX', pres_relax
            ! PRINT *, 'pRT', pRT
            DO i = 1, num_fluids
                ! PRINT *, i
                ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = & 
                q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) * & 
                (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf) +  & 
                q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
                ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
            END DO
            ! ==================================================================
        END SUBROUTINE s_mixture_total_energy_correction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO PHASE PRESSURE FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_alpha1_prelax(p_k,alpha_k)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                                       ::  pstar, f_alpha1_prelax
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)    ::  p_k, alpha_k
            REAL(KIND(0d0))                                       ::  Z1, Z2, pI, C1, C2
            REAL(KIND(0d0))                                       ::  ap, bp, dp
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            Z1 = n1*(p_k(1)+pinf1)
            Z2 = n2*(p_k(2)+pinf2)
            pI = (Z2*p_k(1)+Z1*p_k(2))/(Z1+Z2)
            C1 = 2.d0*n1*pinf1+(n1-1.d0)*p_k(1)
            C2 = 2.d0*n2*pinf2+(n2-1.d0)*p_k(2)
            ap = 1.d0 + n2*alpha_k(1) + n1*alpha_k(2)
            bp = C1*alpha_k(2)+C2*alpha_k(1)-(n2+1.d0)*alpha_k(1)*p_k(1)-(n1+1.d0)*alpha_k(2)*p_k(2)
            dp = -(C2*alpha_k(1)*p_k(1) + C1*alpha_k(2)*p_k(2))
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_prelax = alpha_k(1)*((n1-1.d0)*pstar + 2.d0*p_k(1) + C1)/((n1+1.d0)*pstar+C1)
        END FUNCTION f_alpha1_prelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO PHASE PT FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_alpha1_ptrelax(rhoalpha1,rhoalpha2,E0)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                :: pstar, f_alpha1_ptrelax
            REAL(KIND(0d0)), INTENT(IN)    :: rhoalpha1, rhoalpha2, E0
            REAL(KIND(0d0))                ::  cv1, cv2, q1, q2
            REAL(KIND(0d0))                ::        ap, bp, dp
            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            ap = rhoalpha1*cv1 + rhoalpha2*cv2
            bp = q1*cv1*(n1-1.d0)*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*rhoalpha2*rhoalpha2 + &
                    rhoalpha1*cv1*(n1*pinf1+pinf2) + rhoalpha2*cv2*(n2*pinf2+pinf1) + &
                    rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)+q2*cv1*(n1-1.d0)) - &
                    E0*(cv1*(n1-1.d0)*rhoalpha1 + cv2*(n2-1.d0)*rhoalpha2)
            dp = q1*cv1*(n1-1.d0)*pinf2*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*pinf1*rhoalpha2*rhoalpha2 + &
                    pinf1*pinf2*(rhoalpha1*cv1*n1 + rhoalpha2*cv2*n2) + & 
                    rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)*pinf1 + q2*cv1*(n1-1.d0)*pinf2) - &
                    E0*(cv1*(n1-1.d0)*pinf2*rhoalpha1 + cv2*(n2-1.d0)*pinf1*rhoalpha2)
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_ptrelax = (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1)/&
                        (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1 + &
                        cv2*(n2-1.d0)*(pstar+pinf1)*rhoalpha2)

        END FUNCTION f_alpha1_ptrelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! SATURATION TEMPERATURE FUNCTIONS !!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_fdfTsat(fp,dfdp,pstar,Tstar)

            REAL(KIND(0d0)), INTENT(OUT)        :: fp, dfdp
            REAL(KIND(0d0)), INTENT(IN)         :: pstar, Tstar
            fp = gibbsA + gibbsB/Tstar + gibbsC*DLOG(Tstar) - & 
                    DLOG((pstar+pinf2)/(pstar+pinf1)**gibbsD)
            dfdp = -gibbsB/(Tstar*Tstar) + gibbsC/Tstar

        END SUBROUTINE s_compute_fdfTsat !-------------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_Tsat_bracket(TstarA,TstarB,pressure)
            !> @name In-subroutine variables: vapor and liquid material properties
            !> @{
            REAL(KIND(0d0)), INTENT(OUT)   :: TstarA, TstarB
            REAL(KIND(0d0)), INTENT(IN)    :: pressure
            REAL(KIND(0d0))                :: fA, fB, dfdp, factor
            ! Finding lower bound, getting the bracket 
            factor = 20.d0
            TstarA = TsatLv
            TstarB = TsatLv+factor
            CALL s_compute_fdfTsat(fA,dfdp,pressure,TstarA)
            CALL s_compute_fdfTsat(fB,dfdp,pressure,TstarB)
            !PRINT *, 'fA :: ',fA,', fB :: ',fB
            DO WHILE ( fA*fB .GT. 0.d0 )
                    IF (TstarA .GT. TsatHv) THEN
                            PRINT *, 'Tsat bracketing failed to find lower bound'
                            PRINT *, 'TstarA :: ',TstarA,', pressure :: ',pressure
                            PRINT *, 'fA :: ',fA,', fB :: ',fB
                            CALL s_mpi_abort()
                    END IF
                    fA = fB
                    TstarA = TstarB
                    TstarB = TstarA+factor
                    CALL s_compute_fdfTsat(fB,dfdp,pressure,TstarB)
                    !PRINT *, 'fB :: ',fB,', TstarB :: ',TstarB,', p :: ',pressure
                    !PRINT *,'fB :: ',fB,', TstarB :: ',TstarB
                    IF( ieee_is_nan(fB) ) THEN
                        fB = fA
                        TstarB = TstarA
                        factor = factor-10.d0
                    ELSE 
                        factor = 20.d0
                    END IF
            END DO
        END SUBROUTINE s_compute_Tsat_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_Tsat(pressure)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0)), INTENT(IN)    :: pressure
            REAL(KIND(0d0))                :: f_Tsat, Tstar
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, TstarL, TstarH, TsatA, TsatB
            INTEGER :: iter                !< Generic loop iterators
            CALL s_compute_Tsat_bracket(TsatA,TsatB,pressure)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_fdfTsat(fL,dfdp,pressure,TsatA)
            CALL s_compute_fdfTsat(fH,dfdp,pressure,TsatB)
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                TstarL  = TsatA; TstarH  = TsatB;
            ELSE
                TstarL  = TsatA; TstarH  = TsatB;
            END IF
            Tstar = 0.5d0*(TstarL+TstarH)
            delta_old = DABS(TstarH-TstarL)
            delta = delta_old
            CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((Tstar-TstarH)*dfdp-fp)*((Tstar-TstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                    delta_old = delta
                    delta = 0.5d0*(TstarH-TstarL)
                    Tstar = TstarL + delta
                    IF (delta .EQ. 0.d0) EXIT
                ELSE                    ! Newton step acceptable, take it
                    delta_old = delta
                    delta = fp/dfdp
                    Tstar = Tstar - delta
                    IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/Tstar) < pTsatnewton_eps) EXIT
                CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)           
                IF (fp < 0.d0) THEN     !Maintain the bracket on the root
                    TstarL = Tstar
                ELSE
                    TstarH = Tstar
                END IF
                IF (iter .EQ. newton_iter) THEN
                    PRINT *, 'Tsat : ',Tstar,', iter : ',iter
                    PRINT *, 'Tsat did not converge, stopping code'
                    CALL s_mpi_abort()
                END IF                  
            END DO
            f_Tsat = Tstar
        END FUNCTION f_Tsat !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO-PHASE PTG RELAXATION FUNCTIONS !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0)), INTENT(IN)    :: pstar, rho0, E0
            REAL(KIND(0d0)), INTENT(OUT)   :: fp, dfdp, Tstar
            REAL(KIND(0d0))                ::  cv1, cv2, q1, q2
            REAL(KIND(0d0))                ::        ap, bp, dp
            REAL(KIND(0d0))                ::  dadp, dbdp, dddp
            REAL(KIND(0d0))                ::              dTdp

            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            ap = rho0*cv1*cv2*((n2-1.d0)*(pstar+n1*pinf1)-(n1-1.d0)*(pstar+n2*pinf2))
            bp = E0*((n1-1.d0)*cv1*(pstar+pinf2) - (n2-1.d0)*cv2*(pstar+pinf1)) + &
                    rho0*((n2-1.d0)*cv2*q1*(pstar+pinf1) - (n1-1.d0)*cv1*q2*(pstar+pinf2)) + &
                    cv2*(pstar+pinf1)*(pstar+n2*pinf2) - cv1*(pstar+pinf2)*(pstar+n1*pinf1)
            dp = (q2-q1)*(pstar+pinf1)*(pstar+pinf2)
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            Tstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            ! PRINT *, 'ap: ', ap, 'bp: ', bp, 'dp:', dp
            ! PRINT *, 'sqrt(D)', DSQRT(bp*bp - 4.d0*ap*dp)
            ! PRINT *, 'Tstar: ', Tstar
            ! PRINT *, 'Tstar T : ', (-bp - DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            ! Calculating the derivatives wrt pressure of the coefficients
            dadp = rho0*cv1*cv2*((n2-1.d0)-(n1-1.d0))
            dbdp = E0*((n1-1.d0)*cv1 - (n2-1.d0)*cv2) + &
                    rho0*((n2-1.d0)*cv2*q1 - (n1-1.d0)*cv1*q2) + &
                    cv2*((pstar+pinf1)+(pstar+n2*pinf2)) - & 
                    cv1*((pstar+pinf2)+(pstar+n1*pinf1))
            dddp = (q2-q1)*((pstar+pinf1)+(pstar+pinf2))
            ! Derivative of the temperature wrt to pressure, needed for dfdp
            dTdp = (-dbdp + (0.5d0/DSQRT(bp*bp-4.d0*ap*dp))*(2.d0*bp*dbdp-&
                    4.d0*(ap*dddp+dp*dadp)))/(2.d0*ap) - (dadp/ap)*Tstar
            fp   = gibbsA + gibbsB/Tstar + gibbsC*DLOG(Tstar) + & 
                    gibbsD*DLOG(pstar+pinf1) - DLOG(pstar+pinf2)
            dfdp = -gibbsB/(Tstar*Tstar)*dTdp + gibbsC/Tstar*dTdp + & 
                    gibbsD/(pstar+pinf1) - 1.d0/(pstar+pinf2)

            ! PRINT *, Tstar
            ! PRINT *, dadp, dbdp, dddp, dTdp, fp, dfdp

        END SUBROUTINE s_compute_ptg_fdf !------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0)), INTENT(OUT)   :: pstarA, pstarB
            REAL(KIND(0d0)), INTENT(IN)    :: pstar, rho0, E0
            REAL(KIND(0d0))                :: fA, fB, dfdp, Tstar
            REAL(KIND(0d0))                :: pS, factor

            pstarA = 0.7d0*pstar
            pstarB = pstar
            ! PRINT *, 's_compute_ptg_fdf A'
            CALL s_compute_ptg_fdf(fA,dfdp,pstarA,Tstar,rho0,E0)
            ! PRINT *, 's_compute_ptg_fdf B'
            CALL s_compute_ptg_fdf(fB,dfdp,pstarB,Tstar,rho0,E0)
            factor = 1.05d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                ! PRINT *, 'fA', fA
                ! PRINT *, 'fB', fb
                    IF (pstarA .GT. 1.d13) THEN
                            PRINT *, 'ptg bracketing failed to find lower bound'
                            PRINT *, 'pstarA :: ',pstarA
                            CALL s_mpi_abort()
                    END IF
                    fA = fB
                    pstarA = pstarB
                    pstarB = pstarA*factor
                    CALL s_compute_ptg_fdf(fB,dfdp,pstarB,Tstar,rho0,E0)
                    IF( ieee_is_nan(fB) ) THEN
                        fB = fA
                        pstarB = pstarA
                        factor = factor*0.95d0
                    ELSE 
                        factor = 1.05d0
                    END IF
            END DO

        END SUBROUTINE s_compute_ptg_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_pTrelax(pstar,Tstar,rho0,E0)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0)), INTENT(INOUT) :: pstar
            REAL(KIND(0d0)), INTENT(OUT)   :: Tstar
            REAL(KIND(0d0)), INTENT(IN)    :: rho0, E0
            REAL(KIND(0d0))                ::  pstarA, pstarB, dTstar
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, pstarL, pstarH
            INTEGER :: iter                !< Generic loop iterators
            ! Computing the bracket of the root solution
            CALL s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)
            ! PRINT *, 'pstarA', pstarA, 'pstarB', pstarB 
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_ptg_fdf(fL,dfdp,pstarA,dTstar,rho0,E0)
            ! PRINT *, 'fL', fL, 'dfdpL', dfdp, 'TstarL', dTstar
            ! CALL s_compute_ptg_fdf(fH,dfdp,pstarB,dTstar,rho0,E0)
            ! PRINT *, 'fH', fH, 'dfdpH', dfdp, 'TstarH', dTstar
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                pstarL  = pstarA; pstarH  = pstarB;
            ELSE
                pstarL  = pstarB; pstarH  = pstarA;
            END IF
            pstar = 0.5d0*(pstarA+pstarB)
            delta_old = DABS(pstarB-pstarA)
            delta = delta_old
            CALL s_compute_ptg_fdf(fp,dfdp,pstar,dTstar,rho0,E0)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((pstar-pstarH)*dfdp-fp)*((pstar-pstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                    delta_old = delta
                    delta = 0.5d0*(pstarH-pstarL)
                    pstar = pstarL + delta
                    IF (delta .EQ. 0.d0) EXIT                    ! Change in root is negligible
                ELSE                                              ! Newton step acceptable, take it
                    delta_old = delta
                    delta = fp/dfdp
                    pstar = pstar - delta
                    IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/pstar) < ptgnewton_eps) EXIT           ! Stopping criteria
                ! Updating to next iteration
                CALL s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)
                IF (fp < 0.d0) THEN !Maintain the bracket on the root
                    pstarL = pstar
                ELSE
                    pstarH = pstar
                END IF
            END DO
        END SUBROUTINE s_compute_ptg_pTrelax !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!  K-PHASE P RELAXATION FUNCTIONS !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            ! Cell-average conservative variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0)), INTENT(OUT)                         :: fp, dfdp
            REAL(KIND(0d0)), INTENT(IN)                          :: pstar
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0))                                      :: numerator, denominator
            REAL(KIND(0d0))                                      :: drhodp
            INTEGER, INTENT(IN)                                  :: j, k, l
            INTEGER                                              :: i
            fp = -1.d0; dfdp = 0.d0;
            DO i = 1, num_fluids
                    numerator   = gamma_min(i)*(pstar+pres_inf(i))
                    denominator = numerator + pres_K_init(i)-pstar
                    rho_K_s(i)  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/&
                        MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps)*numerator/denominator
                    drhodp      = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / & 
                        MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) * & 
                        gamma_min(i)*(pres_K_init(i)+pres_inf(i)) / (denominator*denominator)
                    fp          = fp  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                    dfdp        = dfdp - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) * & 
                        drhodp / (rho_K_s(i)*rho_K_s(i))
            END DO

        END SUBROUTINE s_compute_pk_fdf !------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_pk_bracket(pstarA,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0)), INTENT(OUT)   :: pstarA, pstarB
            REAL(KIND(0d0))                :: fA, fB, dfdp
            REAL(KIND(0d0))                :: factor
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s
            INTEGER, INTENT(IN)            :: j,k,l
            !pstarA = 100.d0
            !pstarB = 5.d2
            pstarA = 1.d-15
            pstarB = 1.d1
            CALL s_compute_pk_fdf(fA,dfdp,pstarA,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            CALL s_compute_pk_fdf(fB,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            factor = 10.d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                    IF (pstarA .GT. 1.d13) THEN
                            PRINT *, 'P-K bracketing failed to find lower bound'
                            PRINT *, 'pstarA :: ',pstarA,', pstarB :: ',pstarB
                            PRINT *, 'alpha1 :: ',q_cons_vf(adv_idx%beg)%sf(j,k,l)
                            PRINT *, 'alpha2 :: ',q_cons_vf(adv_idx%beg)%sf(j,k,l)
                            CALL s_mpi_abort()
                    END IF
                    fA = fB
                    pstarA = pstarB
                    pstarB = pstarA*factor
                    CALL s_compute_pk_fdf(fB,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                    IF( ieee_is_nan(fB) ) THEN
                        fB = fA
                        pstarB = pstarA
                        factor = factor*0.5d0
                    ELSE 
                        factor = 10.d0
                    END IF
            END DO
        END SUBROUTINE s_compute_pk_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_p_relax_k(rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0))                                      :: pstar, pstarA, pstarB
            REAL(KIND(0d0))                                      :: delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                                      :: fL, fH, pstarL, pstarH
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s 
            INTEGER, INTENT(IN)                                  :: j,k,l
            INTEGER :: iter                !< Generic loop iterators
            ! Computing the bracket of the root solution
            CALL s_compute_pk_bracket(pstarA,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_pk_fdf(fL,dfdp,pstarA,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            CALL s_compute_pk_fdf(fH,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                pstarL  = pstarA; pstarH  = pstarB;
            ELSE
                pstarL  = pstarB; pstarH  = pstarA;
            END IF
            pstar = 0.5d0*(pstarA+pstarB)
            delta_old = DABS(pstarB-pstarA)
            delta = delta_old
            CALL s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((pstar-pstarH)*dfdp-fp)*((pstar-pstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                    delta_old = delta
                    delta = 0.5d0*(pstarH-pstarL)
                    pstar = pstarL + delta
                    IF (delta .EQ. 0.d0) EXIT                    ! Change in root is negligible
                ELSE                                              ! Newton step acceptable, take it
                    delta_old = delta
                    delta = fp/dfdp
                    pstar = pstar - delta
                    IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/pstar) < pknewton_eps) EXIT           ! Stopping criteria
                ! Updating to next iteration
                CALL s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                IF (fp < 0.d0) THEN !Maintain the bracket on the root
                    pstarL = pstar
                ELSE
                    pstarH = pstar
                END IF
            END DO
        END SUBROUTINE s_compute_p_relax_k !-------------------------------

        ! This SUBROUTINE IS USED TO CALCULATE THE 2X2 JACOBIAN AND ITS INVERSE FOR THE PROPOSED
        ! pTg-Equilibrium procedure
        SUBROUTINE s_compute_jacobian_matrix( B, C, cv, D, rhoe, g_min, InvJac, j, Jac, k, l, lp, mCPD, &
                                                mCVGP, mCVGP2, mQD, p_inf, pS, qv, qvp, q_cons_vf, TJac, vp )
        
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)        :: q_cons_vf
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                          :: B, C, D
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                          :: pS, rhoe, mCPD, mCVGP, mCVGP2, mQD
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)   :: g_min, p_inf, cv, qv, qvp
            REAL( KIND( 0.0D0 ) ), DIMENSION(2,2), INTENT(OUT)         :: Jac, InvJac, TJac
            INTEGER, INTENT(IN)                                        :: j, k, l, lp, vp
            REAL( KIND( 0.0D0 ) )                                      :: ml, mT, TS, dFdT, dTdm, dTdp
            INTEGER                                                    :: i
            
            ! mass of the reactant liquid
            ml = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )
            
            ! mass of the two participating fluids
            mT = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )         

            TS =  1 / ( mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )    &
            +           ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )    &
            -                  cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )  &
            +       mCVGP )

            dFdT =                                                              &
            - ( cv( lp ) * g_min( lp ) - cv( vp ) * g_min( vp ) ) * DLOG( TS )  &
            - ( qvp( lp ) - qvp( vp ) )                                         &
            + cv( lp ) * ( g_min( lp ) - 1 ) * DLOG( pS + p_inf( lp ) )         &
            - cv( vp ) * ( g_min( vp ) - 1 ) * DLOG( pS + p_inf( vp ) )

            dTdm = - ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )  &
            -          cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ) * TS ** 2

            dTdp =   ( mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ** 2      &
            +          ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) ** 2      &
            -                 cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ** 2 )    &
            +        mCVGP2 ) * TS ** 2

            ! F = (F1,F2) is the fuction whose roots we are looking for
            ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
            ! F1 = 0 is the gibbs free energy quality
            ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy           
            ! dF1dm
            Jac(1,1) = dFdT * dTdm

            Jac(1,2) = dFdT * dTdp + TS &
            * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) &
            -   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )

            ! dF2dm
            Jac(2,1) = ( qv( vp ) - qv( lp )                                                        &
            + ( cv( vp ) * g_min( vp ) - cv( lp ) * g_min( lp ) )                                   &
            / ( ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                        &
                        - cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                      & 
                + mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP )              &
            - ( ml * ( cv( vp ) * g_min( vp ) - cv( lp ) * g_min( lp ) )                            &
            -   mT * cv( vp ) * g_min( vp ) - mCPD )                                                &
            * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                               &
            -   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                             &
            / ( ( ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                      &
            -            cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                    &
            +     mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) ** 2 ) ) / 1
            
            Jac(2,2) = ( 1 + ( ml * ( cv( vp ) * g_min( vp ) - cv( lp ) * g_min( lp ) )             &
            -                mT *   cv( vp ) * g_min( vp ) - mCPD )                                 & 
            * ( ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )  ** 2                  & 
                        - cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )  ** 2 )                &
            +   mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )  ** 2 + mCVGP2 )       &
            / ( ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                        &   
                        - cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                      &
            +   mT *   cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) ** 2 ) / 1

            ! elements of J^{-1}
            InvJac(1,1) =  Jac(2,2)
            InvJac(1,2) = - 1.0D0 * Jac(1,2)
            InvJac(2,1) = - 1.0D0 * Jac(2,1)
            InvJac(2,2) = Jac(1,1)

            ! elements of J^{T}
            TJac(1,1) = Jac(1,1)
            TJac(1,2) = Jac(2,1)
            TJac(2,1) = Jac(1,2)
            TJac(2,2) = Jac(2,2)

            ! dividing by det(J)
            InvJac = InvJac / ( Jac(1,1) * Jac(2,2) - Jac(1,2) * Jac(2,1) )             

        END SUBROUTINE s_compute_jacobian_matrix

        SUBROUTINE s_compute_pTg_residue(A, B, C, cv, D, g_min, j, k, l, lp, mCPD, mCVGP, mQD &
                                            , p_inf, qv, qvp, q_cons_vf, pS, rhoe, R2D, vp)
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)         :: q_cons_vf
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: A, B, C, D
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: pS, rhoe, mCPD, mCVGP, mQD
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)    :: g_min, p_inf, cv, qv, qvp
            REAL( KIND( 0.0D0 ) ), DIMENSION(2), INTENT(OUT)            :: R2D
            INTEGER, INTENT(IN)                                         :: j, k, l, lp, vp
            REAL( KIND( 0.0D0 ) )                                       :: ml, mT, TS
            INTEGER                                                     :: i
            
            ! mass of the reactant liquid
            ml = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )
            
            ! mass of the two participating fluids
            mT = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )            

            TS =  1 / ( mT * cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )      &
            +           ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )    &
            -                  cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )  &
            +       mCVGP )

            ! Gibbs Free Energy Equality condition (DG)
            R2D( 1 ) = TS * ( ( cv( lp ) * g_min( lp ) - cv( vp ) * g_min( vp ) )   &
            * ( 1 - DLOG( TS ) ) - ( qvp( lp ) - qvp( vp ) )                        &
            + cv( lp ) * ( g_min( lp ) - 1 ) * DLOG( pS + p_inf( lp ) )             &
            - cv( vp ) * ( g_min( vp ) - 1 ) * DLOG( pS + p_inf( vp ) ) )           &
            + qv( lp ) - qv( vp )

            ! Constant Energy Process condition (DE)
            R2D( 2 ) = ( rhoe + pS                                                          &
            + ml * ( qv( vp ) - qv( lp ) ) - mT * qv( vp ) - mQD                            &
            + ( ml * ( g_min( vp ) * cv( vp ) - g_min( lp ) * cv( lp ) )                    &
            - mT * g_min( vp ) * cv( vp ) - mCPD )                                          &
            / ( ml * ( cv( lp ) * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                &
                        - cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )           &
            + mT * cv( vp ) * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) ) / 1

        END SUBROUTINE s_compute_pTg_residue

        ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
        SUBROUTINE s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_inf, pS, qv, R2D, rhoe, q_cons_vf, TS ) ! ----------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)         :: q_cons_vf 
            REAL( KIND( 0.0D0 ) ), DIMENSION(2,2), INTENT(IN)           :: Jac, InvJac
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)    :: p_inf, qv
            REAL( KIND( 0.0D0 ) ), DIMENSION(2), INTENT(IN)             :: R2D, DeltamP
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: pS, TS
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: rhoe, mQ
            INTEGER, INTENT(IN)                                         :: j, k, l
            REAL( KIND( 0.0D0 ) )                                       :: rho
            !< Generic loop iterator
            INTEGER                                                     :: i

            PRINT *, 'j, k, l', j, k, l

            PRINT *, 'rhoe', rhoe

            PRINT *, 'mQ', mQ

            PRINT *, 'Energy constrain', ( rhoe - mQ - MINVAL( p_inf ) )

            PRINT *, 'R2D', R2D

            PRINT *, 'l2(R2D)', DSQRT( R2D( 1 ) **2 + R2D( 2 ) ** 2 )

            PRINT *, 'DeltamP', DeltamP

            PRINT *, 'pS', pS

            PRINT *, '-min(p_inf)', - MINVAL( p_inf )

            PRINT *, 'TS', TS

            DO i = 1, num_fluids

                rho = rho + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )
                
                PRINT *, 'rho', rho

            END DO

            DO i = 1, num_fluids

                PRINT *, 'i', i

                PRINT *, 'alpha_i', q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                PRINT *, 'alpha_rho_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                PRINT *, 'mq_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                        * qv( i )

                PRINT *, 'rho', rho

                PRINT *, 'internal energies', q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l )

                PRINT *, 'Y_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) / rho

            END DO
            
            PRINT *, 'J', Jac, 'J-1', InvJac

        END SUBROUTINE s_tattletale

        ! Newton Solver for the finding the Saturation temperature TSat for a given saturation pressure
        SUBROUTINE s_TSat(A, B, C, cv, D, g_min, lp, pSat, p_inf, qv, qvp, TSat, TSIn, vp)

            REAL( KIND( 0.0D0 ) ), INTENT( OUT )                         :: TSat
            REAL( KIND( 0.0D0 ) ), DIMENSION( num_fluids ), INTENT( IN ) :: cv, g_min, p_inf, qv, qvp
            REAL( KIND( 0.0D0 ) ), INTENT( IN )                          :: pSat, TSIn, A, B, C, D
            INTEGER, INTENT( IN )                                        :: lp, vp
            REAL( KIND( 0.0D0 ) )                                        :: dFdT, FT

            ! Generic loop iterators
            INTEGER :: ns    

            IF ( ( pSat .EQ. 0.0D0 ) .AND. ( TSIn .EQ. 0.0D0 ) ) THEN
            
                ! assigning Saturation temperature
                TSat = 0.0D0

            ELSE

                ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to 
                ! iterate over the Newton's solver
                TSat = TSIn
                
                ! iteration counter
                ns = 0

                ! DO WHILE ( ( DABS( FT ) .GT. ptgalpha_eps .AND. DABS( FT / TSat ) .GT. ptgalpha_eps ) & 
                DO WHILE ( ( DABS( FT ) .GT. ptgalpha_eps ) & 
                        .OR. ( ns .EQ. 0 ) )
                    
                    ! increasing counter
                    ns = ns + 1

                    ! calculating residual
                    ! FT = A + B / TSat + C * DLOG( TSat ) & 
                    ! + D * DLOG( ( pSat + p_inf( lp ) ) ) - DLOG( pSat + p_inf( vp ) )

                    FT = TSat * ( ( cv( lp ) * g_min( lp ) - cv( vp ) * g_min( vp ) )   &
                    * ( 1 - DLOG( TSat ) ) - ( qvp( lp ) - qvp( vp ) )                  &
                    + cv( lp ) * ( g_min( lp ) - 1 ) * DLOG( pSat + p_inf( lp ) )       &
                    - cv( vp ) * ( g_min( vp ) - 1 ) * DLOG( pSat + p_inf( vp ) ) )     &
                    + qv( lp ) - qv( vp )

                    ! calculating the jacobian
                    ! dFdT = - B / ( TSat ** 2) + C / TSat

                    dFdT =                                                                  &
                    - ( cv( lp ) * g_min( lp ) - cv( vp ) * g_min( vp ) ) * DLOG( TSat )    &
                    - ( qvp( lp ) - qvp( vp ) )                                             &
                    + cv( lp ) * ( g_min( lp ) - 1 ) * DLOG( pSat + p_inf( lp ) )           &
                    - cv( vp ) * ( g_min( vp ) - 1 ) * DLOG( pSat + p_inf( vp ) )

                    ! updating saturation temperature
                    TSat = TSat - FT / dFdT

                    ! Checking if TSat returns a NaN
                    IF ( ISNAN( TSat ) ) THEN

                        CALL s_mpi_abort('TSat returned NaN values, when it should not (by assumption &
                        of first order transition - m_phase_change, s_TSat). Aborting!')
        
                    ! checking if the maximum number of iterations has been reached
                    ELSEIF ( ns .GT. 1E4 ) THEN
            
                            CALL s_mpi_abort('Maximum number of iterations reached for TSat &
                            (m_phase_change, s_TSat). Aborting!')
                    END IF

                END DO

            END IF

        END SUBROUTINE s_TSat

        SUBROUTINE s_finalize_relaxation_solver_module()

            s_relaxation_solver => NULL()

        END SUBROUTINE

#endif        

END MODULE m_phase_change