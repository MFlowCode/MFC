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

    implicit none

    private; public :: s_initialize_phasechange_module, &
 s_relaxation_solver, &
 s_relaxation_finite_solver, &
 s_finite_ptg_relaxation, &
 s_infinite_p_relaxation, &
 s_infinite_p_relaxation_k, &
 s_infinite_pt_relaxation, &
 s_infinite_relaxation_k, &
 s_infinite_ptg_relaxation, &
 s_finalize_relaxation_solver_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the finite relaxation solver
        !> @{
        subroutine s_abstract_relaxation_solver(q_cons_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        end subroutine

        !> @name Abstract subroutine for the finite relaxation solver
        !> @{
        subroutine s_abstract_relaxation_finite_solver(q_cons_vf, rhs_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
            type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        end subroutine

    end interface

    !> @name Parameters for the phase change part of the code
    !> @{
    integer, parameter :: newton_iter = 50        !< p_relaxk \alpha iter,                set to 25
    real(kind(0d0)), parameter :: pknewton_eps = 1.d-15    !< p_relaxk \alpha threshold,           set to 1E-15
    real(kind(0d0)), parameter :: pTsatnewton_eps = 1.d-10    !< Saturation temperature tol,          set to 1E-10
    real(kind(0d0)), parameter :: ptgnewton_eps = 1.d-8     !< Saturation p-T-mu tolerance,         set to 1.d-10
    real(kind(0d0)), parameter :: pres_crit = 22.09d6   !< Critical water pressure              set to 22.06d6
    real(kind(0d0)), parameter :: T_crit = 647.29d0  !< Critical water temperature           set to 648
    real(kind(0d0)), parameter :: TsatHv = 1000.d0   !< Saturation temperature threshold,    set to 900
    real(kind(0d0)), parameter :: TsatLv = 250.d0    !< Saturation temperature threshold,    set to 250
    integer, parameter         :: lp = 1    ! index for the liquid phase of the reacting fluid
    integer, parameter         :: vp = 2    ! index for the vapor phase of the reacting fluid

    !> @name Gibbs free energy phase change parameters
    !> @{b
    real(kind(0d0)) :: A, B, C, D
    real(kind(0.0d0)), dimension(2) :: g_minI, p_infI, cvI, qvI, qvpI
    !> @}

    procedure(s_abstract_relaxation_solver), &
        pointer :: s_relaxation_solver => null()

    procedure(s_abstract_relaxation_finite_solver), &
        pointer :: s_relaxation_finite_solver => null()

contains

    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_initialize_phasechange_module()

        ! renaming fluid properties JUST FOR EASINESS
        cvI(lp) = fluid_pp(lp)%cv
        cvI(vp) = fluid_pp(vp)%cv
        qvI(lp) = fluid_pp(lp)%qv
        qvI(vp) = fluid_pp(vp)%qv
        qvpI(lp) = fluid_pp(lp)%qvp
        qvpI(vp) = fluid_pp(vp)%qvp

        ! gamma_min
        g_minI(lp) = 1.0d0/fluid_pp(lp)%gamma + 1.0d0
        g_minI(vp) = 1.0d0/fluid_pp(vp)%gamma + 1.0d0

        ! p_inf
        p_infI(lp) = fluid_pp(lp)%pi_inf/(1.0d0 + fluid_pp(lp)%gamma)
        p_infI(vP) = fluid_pp(Vp)%pi_inf/(1.0d0 + fluid_pp(Vp)%gamma)

        g_minI(lp) = 1.d0/fluid_pp(1)%gamma + 1.d0
        p_infI(lp) = fluid_pp(1)%pi_inf/(1.d0 + fluid_pp(1)%gamma)

        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (g_minI(lp)*cvI(lp) - g_minI(vp)*cvI(vp) &
             + qvI(vp) - qvI(lp))/((g_minI(vp) - 1.0d0)*cvI(vp))

        B = (qvI(lp) - qvI(vp))/((g_minI(vp) - 1.0d0)*cvI(vp))

        C = (g_minI(vp)*cvI(vp) - g_minI(lp)*cvI(lp)) &
            /((g_minI(vp) - 1.0d0)*cvI(vp))

        D = ((g_minI(lp) - 1.0d0)*cvI(lp)) &
            /((g_minI(vp) - 1.0d0)*cvI(vp))

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem

        ! Opening and writing the header of the run-time information file
        if (relax_model == 0 .or. relax_model == 1) then
            s_relaxation_solver => s_infinite_p_relaxation
        elseif (relax_model == 2) then
            s_relaxation_solver => s_infinite_pt_relaxation
        elseif (relax_model == 3) then
            s_relaxation_solver => s_infinite_ptg_relaxation
        elseif (relax_model == 4) then
            s_relaxation_solver => s_infinite_p_relaxation_k
        elseif ((relax_model == 5) .or. (relax_model == 6)) then
            s_relaxation_solver => s_infinite_relaxation_k
        else
            print '(A)', 'relaxation solver was not set!'
            call s_mpi_abort()
        end if

        !IF (relax_model == 1) THEN
        !    s_relaxation_finite_solver => s_finite_ptg_relaxation
        !END IF

    end subroutine s_initialize_phasechange_module !-------------------------------

    !> The purpose of this procedure is to employ the inputted
        !!      cell-average conservative variables in order to compute
        !!      the cell-average RHS variables of the semidiscrete form
        !!      of the governing equations by utilizing the appropriate
        !!      Riemann solver.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param rhs_vf Cell-average RHS variables
    subroutine s_finite_ptg_relaxation(q_cons_vf, rhs_vf) ! -------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        real(kind(0d0)) :: sum_alpha, Tsat
        real(kind(0d0)), dimension(num_fluids) :: p_k, T_k, g_k, Z_k
        real(kind(0d0)) :: n_k, pinf_k
        real(kind(0d0)) :: rho_k, rhoeq_k, rhoe
        real(kind(0d0)) :: e_k, phi, psi
        real(kind(0d0)) :: f1, f2, f3, f4
        real(kind(0d0)) :: A1, B1, A2, B2
        real(kind(0d0)) :: rho_I, Q, kappa
        real(kind(0d0)) :: deltap, p_I, mu
        real(kind(0d0)) :: e_I, mdot, nu, theta
        real(kind(0d0)) :: mdotalpha, mdotrhoe

        integer :: i, j, k, l, r !< Generic loop iterators

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                    end if

                    if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps .or. &
                        q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) then

                        ! computing n_k, pinf_k, p_k, T_k, and g_k for finite relaxation
                        phi = 0.d0; psi = 0.d0; f1 = 0.d0; f2 = 0.d0; f3 = 0.d0; f4 = 0.d0; 
                        !rhoe = 0.d0;

                        do i = 1, num_fluids
                            n_k = 1.d0/fluid_pp(i)%gamma + 1.d0
                            pinf_k = fluid_pp(i)%pi_inf/(1.d0 + fluid_pp(i)%gamma)
                            rho_k = q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                                    /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            rhoeq_k = (q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                       - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                      /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            !rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            e_k = q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                  /q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)

                            p_k(i) = (rhoeq_k - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            T_k(i) = (rhoeq_k - fluid_pp(i)%pi_inf/(1.d0 + fluid_pp(i)%gamma)) &
                                     /(rho_k*fluid_pp(i)%cv)
                            g_k(i) = (n_k*fluid_pp(i)%cv - fluid_pp(i)%qvp)*T_k(i) &
                                     - fluid_pp(i)%cv*T_k(i)*log(T_k(i)**(n_k) &
                                                                 /((p_k(i) + pinf_k)**(n_k - 1.d0))) + fluid_pp(i)%qv
                            Z_k(i) = n_k*(p_k(i) + pinf_k); 
                            phi = phi + 1.d0/(q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%cv)
                            psi = psi + 1.d0/(fluid_pp(i)%gamma*q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l))

                            f1 = f1 + (p_k(i) + n_k*pinf_k)/q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            f2 = f2 + pinf_k/(q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%cv)
                            f3 = f3 + fluid_pp(i)%qv/(fluid_pp(i)%gamma &
                                                      *q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l))
                            f4 = f4 + (e_k - pinf_k/rho_k) &
                                 /(q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%cv)
                        end do
                        !TODO IMPROVE THIS APPROACH
                        Tsat = f_Tsat(p_k(1))
                        !PRINT *,'prelax =',p_k(1),'Trelax = ',T_k(1),', Tsat = ',Tsat

                        if (ieee_is_nan(p_k(1)) .or. p_k(1) < 0.d0) then
                            print *, 'code crashed'
                            call s_mpi_abort()
                        end if

                        kappa = f1/psi
                        rho_I = (phi*f1 - psi*f2)/(psi*f4 - phi*f3)
                        p_I = (Z_k(2)*p_k(1) + Z_k(1)*p_k(2))/(Z_k(1) + Z_k(2)); 
                        e_I = f4/phi + f2/(rho_I*phi)

                        !mu = 1.d8
                        mu = 0.d0; 
                        theta = 1.d8
                        nu = 0.d0
                        if (T_k(1) > Tsat) nu = 1d-3

                        deltap = mu*(p_k(1) - p_k(2))
                        Q = theta*(T_k(2) - T_k(1))
                        mdot = nu*(g_k(2) - g_k(1))
                        mdotalpha = mdot/rho_I
                        mdotrhoe = mdot*e_I

                        rhs_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(1 + adv_idx%beg - 1)%sf(j, k, l) + deltap + Q/kappa + mdotalpha
                        rhs_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(2 + adv_idx%beg - 1)%sf(j, k, l) - deltap - Q/kappa - mdotalpha
                        rhs_vf(1 + cont_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(1 + cont_idx%beg - 1)%sf(j, k, l) + mdot
                        rhs_vf(2 + cont_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(2 + cont_idx%beg - 1)%sf(j, k, l) - mdot
                        rhs_vf(1 + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(1 + internalEnergies_idx%beg - 1)%sf(j, k, l) - p_I*deltap + Q + mdotrhoe
                        rhs_vf(2 + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                            rhs_vf(2 + internalEnergies_idx%beg - 1)%sf(j, k, l) + p_I*deltap - Q - mdotrhoe
                    end if
                end do
            end do
        end do
    end subroutine s_finite_ptg_relaxation ! --------------------------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_p_relaxation(q_cons_vf) ! ----------------
        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0d0)) :: rhoeq_k
        real(kind(0d0)) :: a1
        real(kind(0d0)), dimension(num_fluids) :: p_k, alpha_k
        !> @}
        integer :: i, j, k, l           !< Generic loop iterators
        logical :: relax             !< Relaxation procedure determination variable
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! P RELAXATION ==================================
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) relax = .true.
                    if (relax) then
                        do i = 1, num_fluids
                            alpha_k(i) = q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            rhoeq_k = (q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                       - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                      /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            p_k(i) = (rhoeq_k - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                        end do
                        a1 = f_alpha1_prelax(p_k, alpha_k)
                        ! Cell update of the volume fraction
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = a1
                        q_cons_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = 1.d0 - a1
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                    end if
                end do
            end do
        end do

    end subroutine s_infinite_p_relaxation ! ----------------

    !> Description: The purpose of this procedure is to infinitely relax
        !!              the pressures from the internal-energy equations to a
        !!              unique pressure, from which the corresponding volume
        !!              fraction of each phase are recomputed. For conservation
        !!              purpose, this pressure is finally corrected using the
        !!              mixture-total-energy equation.
    subroutine s_infinite_p_relaxation_k(q_cons_vf) ! ----------------
        ! Relaxed pressure, initial partial pressures, function f(p) and its partial
        ! derivative df(p), isentropic partial density, sum of volume fractions,
        ! mixture density, dynamic pressure, surface energy, specific heat ratio
        ! function, liquid stiffness function (two variations of the last two
        ! ones), shear and volume Reynolds numbers and the Weber numbers
        ! Cell-average conservative variables
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s, pres_K_init
        real(kind(0d0)), dimension(num_fluids) :: gamma_min, pres_inf
        ! Generic loop iterators
        integer :: i, j, k, l
        ! Relaxation procedure determination variable
        logical :: relax
        do i = 1, num_fluids
            gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
            pres_inf(i) = fluid_pp(i)%pi_inf/(1d0 + fluid_pp(i)%gamma)
        end do
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    ! Pressures relaxation procedure ===================================
                    relax = .true.
                    do i = 1, num_fluids
                        if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) > (1.d0 - palpha_eps)) relax = .false.
                    end do
                    if (relax) then
                        ! Calculating the initial pressure
                        do i = 1, num_fluids
                            pres_K_init(i) = ((q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                               - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                              /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) &
                                              - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            if (pres_K_init(i) <= -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-0
                        end do
                        call s_compute_p_relax_k(rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
                        ! Cell update of the volume fraction
                        do i = 1, num_fluids
                            if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) &
                                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = &
                                q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/rho_K_s(i)
                        end do
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                end do
            end do
        end do
    end subroutine s_infinite_p_relaxation_k ! -----------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_pt_relaxation(q_cons_vf) ! ----------------
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)), dimension(num_fluids) :: p_k, alpha_k
        real(kind(0d0)) :: rhoeq_k, rhoe, a1
        real(kind(0d0)) :: rhoalpha1, rhoalpha2
        !> @}
        integer :: i, j, k, l        !< Generic loop iterators
        logical :: relax             !< Relaxation procedure determination variable
        !< Computing the constant saturation properties

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! P RELAXATION ==================================
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) relax = .true.
                    if (relax) then
                        do i = 1, num_fluids
                            alpha_k(i) = q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            rhoeq_k = (q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                       - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                      /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            p_k(i) = (rhoeq_k - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                        end do
                        a1 = f_alpha1_prelax(p_k, alpha_k)
                        ! Cell update of the volume fraction
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = a1
                        q_cons_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = 1.d0 - a1
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)

                    ! PT RELAXATION ==================================
                    rhoe = 0.d0
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) relax = .true.
                    if (relax) then
                        rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j, k, l)
                        rhoalpha2 = q_cons_vf(1 + cont_idx%beg)%sf(j, k, l)
                        do i = 1, num_fluids
                            rhoe = rhoe + q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)
                        end do
                        a1 = f_alpha1_ptrelax(rhoalpha1, rhoalpha2, rhoe)
                        ! Cell update of the volume fraction
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = a1
                        q_cons_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = 1.d0 - a1
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                end do
            end do
        end do
    end subroutine s_infinite_pt_relaxation ! -----------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_ptg_relaxation(q_cons_vf) ! ----------------
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)) :: pres_relax, Trelax
        real(kind(0d0)), dimension(num_fluids) :: p_k, alpha_k, Tk
        real(kind(0d0)) :: rhoalpha1, rhoalpha2
        real(kind(0d0)) :: rho, rhoe, rhoeq_k
        real(kind(0d0)) :: rho1, rho2
        real(kind(0d0)) :: a1, a2
        real(kind(0d0)) :: gamma, pi_inf, p_infk
        real(kind(0d0)) :: qv
        real(kind(0d0)) :: Tsat

        !> @}
        integer :: i, j, k, l        !< Generic loop iterators
        logical :: relax             !< Relaxation procedure determination variable
        !< Computing the constant saturation properties

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! P RELAXATION ========================================
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, rho, &
                                                        gamma, pi_inf, qv)

                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) relax = .true.
                    if (relax) then
                        do i = 1, num_fluids
                            alpha_k(i) = q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            rhoeq_k = (q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                       - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                      /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)
                            p_k(i) = (rhoeq_k - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                        end do
                        a1 = f_alpha1_prelax(p_k, alpha_k)
                        ! Cell update of the volume fraction
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = a1
                        q_cons_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = 1.d0 - a1
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                    ! PT RELAXATION ==========================================
                    rhoe = 0.d0
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, rho, &
                                                        gamma, pi_inf, qv)
                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > palpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - palpha_eps) relax = .true.
                    if (relax) then
                        rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j, k, l)
                        rhoalpha2 = q_cons_vf(1 + cont_idx%beg)%sf(j, k, l)
                        do i = 1, num_fluids
                            rhoe = rhoe + q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)
                        end do
                        a1 = f_alpha1_ptrelax(rhoalpha1, rhoalpha2, rhoe)
                        ! Cell update of the volume fraction
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) = a1
                        q_cons_vf(2 + adv_idx%beg - 1)%sf(j, k, l) = 1.d0 - a1
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                    ! CHECKING IF PTG RELAXATION IS NEEDED  =====================
                    rhoe = 0.d0
                    relax = .false.
                    if (mpp_lim) then
                        call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
                    end if
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l, rho, &
                                                        gamma, pi_inf, qv)
                    if ((q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) > ptgalpha_eps) .and. &
                        q_cons_vf(1 + adv_idx%beg - 1)%sf(j, k, l) < 1.d0 - ptgalpha_eps) relax = .true.
                    if (relax) then
                        do i = 1, num_fluids
                            rhoe = rhoe + q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)
                        end do
                        pres_relax = (rhoe - pi_inf)/gamma
                        Tsat = f_Tsat(pres_relax)
                        do i = 1, num_fluids
                            Tk(i) = ((q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                      - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv) &
                                     /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) &
                                     - fluid_pp(i)%pi_inf &
                                     /(1.d0 + fluid_pp(i)%gamma)) &
                                    /(q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%cv &
                                      /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l))
                        end do
                        if (Tk(1) < Tsat) relax = .false.
                    end if
                    ! PTG RELAXATION PROCEDURE ===========================
                    if (relax) then
                        call s_compute_ptg_pTrelax(pres_relax, Trelax, rho, rhoe)
                        p_infk = fluid_pp(1)%pi_inf/(1.d0 + fluid_pp(1)%gamma)
                        rho1 = (pres_relax + p_infk)*fluid_pp(1)%gamma/ &
                               (fluid_pp(1)%cv*Trelax)
                        p_infk = fluid_pp(2)%pi_inf/(1.d0 + fluid_pp(2)%gamma)
                        rho2 = (pres_relax + p_infk)*fluid_pp(2)%gamma/ &
                               (fluid_pp(2)%cv*Trelax)
                        ! Calculate vapor and liquid volume fractions
                        a1 = (rho - rho2)/(rho1 - rho2)
                        a2 = 1.d0 - a1
                        ! Cell update of the volume fraction
                        q_cons_vf(cont_idx%beg)%sf(j, k, l) = rho1*a1
                        q_cons_vf(1 + cont_idx%beg)%sf(j, k, l) = rho2*a2
                        q_cons_vf(adv_idx%beg)%sf(j, k, l) = a1
                        q_cons_vf(1 + adv_idx%beg)%sf(j, k, l) = a2
                    end if
                    call s_mixture_total_energy_correction(q_cons_vf, j, k, l)
                end do
            end do
        end do
    end subroutine s_infinite_ptg_relaxation ! -----------------------

    ! This subroutine is created to activate either the pT- or the pTg-equilibrium
    ! model, with changes such that mass depletion is taken into consideration
    subroutine s_infinite_relaxation_k(q_cons_vf) ! ----------------
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids) :: g_min, p_inf, p_infA, pDumb, cv, qv, qvp
        real(kind(0.0d0)), dimension(num_fluids) :: pk, sk, hk, gk
        real(kind(0.0d0)), dimension(num_fluids) :: ek, rhok
        real(kind(0.0d0)) :: pS, pSOV, pSSL
        real(kind(0.0d0)) :: TS, TSOV, TSatOV, TSatSL, TSSL
        real(kind(0.0d0)) :: rhoe, Ewe, mCP, mQ, TvF
        real(kind(0.0d0)) :: rho, rM, m1, m2, MCT, mixM, rhoT, rMT
        real(kind(0.0d0)) :: TvFT, rhos
        real(kind(0.0d0)) :: dynE, VFT

        !< Generic loop iterators
        integer :: i, j, k, l
        integer :: ns

        ! renaming fluid properties JUST FOR EASINESS
        cv(1:num_fluids) = fluid_pp(1:num_fluids)%cv
        qv(1:num_fluids) = fluid_pp(1:num_fluids)%qv
        qvp(1:num_fluids) = fluid_pp(1:num_fluids)%qvp

        ! gamma_min
        g_min(1:num_fluids) = 1.0d0/fluid_pp(1:num_fluids)%gamma + 1.0d0

        ! p_inf
        p_inf(1:num_fluids) = fluid_pp(1:num_fluids)%pi_inf &
                              /(1.0d0 + fluid_pp(1:num_fluids)%gamma)

        ! threshold for 'mixture cell'. If Y < mixM, the cell is not considered a mixture for pTg purposes
        mixM = 1.0d-08

        ! starting equilibrium solver
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0d0
                    rhos = 0.0d0
                    TvF = 0.0d0
                    do i = 1, num_fluids

                        ! pressure that comes from the homogeneous solver at that cell
                        pk(i) = (g_min(i) - 1) &
                                *(q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) &
                                  - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*qv(i)) &
                                /q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) &
                                - g_min(i)*p_inf(i)

                        ! Mixture density
                        rho = rho + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)

                        ! total entropy
                        rhos = rhos + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*sk(i)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)

                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    m1 = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l)

                    rM = m1 + m2

                    ! correcting negative volume and mass fraction values in case they happen
                    ! IF ( mpp_lim ) THEN

                    call s_correct_partial_densities(g_min, MCT, mixM, pK, p_inf, q_cons_vf, 0.0d0, rM, j, k, l)

                    ! END IF

                    ! kinetic energy so as to calculate the total internal energy as the total energy minus the
                    ! kinetic energy.rMT
                    dynE = 0.0d0
                    do i = mom_idx%beg, mom_idx%end

                        dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/rho

                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! at each of the cells. Note I calculate UE as TE - KE due to numerical reasons (stability at the dicontinuities)
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                    call s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'pT', pk, pS, p_inf &
                                                    , p_infA, qv, q_cons_vf, rho, rhoe, TS)

                    ! check if pTg-equilibrium is required
                    ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! given the outputs from the pT- and pTg-equilibrium solvers are just p, T (byproduct)
                    ! , and the partial masses (pTg- case)
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) &
                                                   > MixM*rM) .or. (q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) &
                                                                    > MixM*rM)) .and. (pS < pres_crit) .and. (TS < T_crit)) then

                        ! Checking if phase change can happen.
                        ! overheated vapor
                        ! correcting the liquid partial density
                        q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = mixM*rM

                        ! correcting the vapor partial density
                        q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! calling pT-equilibrium for overheated vapor
                        call s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'OV', pk, pSOV, p_inf &
                                                        , pDumb, qv, q_cons_vf, rho, rhoe, TSOV)

                        ! calculating Saturation temperature
                        call s_TSat(cv, g_min, pSOV, p_inf, qv, qvp, TSatOV, TSOV)

                        ! subcooled liquid
                        ! correcting the liquid partial density
                        q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! correcting the vapor partial density
                        q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = mixM*rM

                        ! calling pT-equilibrium for subcooled liquid
                        call s_infinite_pt_relaxation_k(cv, j, k, l, g_min, 'SL', pk, pSSL, p_inf &
                                                        , pDumb, qv, q_cons_vf, rho, rhoe, TSSL)

                        ! calculating Saturation temperature
                        call s_TSat(cv, g_min, pSSL, p_inf, qv, qvp, TSatSL, TSSL)

                        ! Maybe I will have to change this for .GE. instead, since mass depletion can
                        ! happen at pTg-eq. Make sure it works, first
                        if (TSOV > TSatOV) then

                            ! Assigning pressure
                            pS = pSOV

                            ! Assigning Temperature
                            TS = TSOV

                            ! correcting the liquid partial density
                            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = mixM*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        else if (TSSL < TSatSL) then

                            ! Assigning pressure
                            pS = pSSL

                            ! Assigning Temperature
                            TS = TSSL

                            ! correcting the liquid partial density
                            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = mixM*rM

                        else

                            ! returning partial pressures to what they were from the homogeneous solver
                            ! liquid
                            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = m1

                            ! vapor
                            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = m2

                            ! calling the pTg-equilibrium solver
                            call s_infinite_ptg_relaxation_k(cv, g_min, j, k, l, pS, p_inf, p_infA, rhoe, qv, qvp, q_cons_vf, TS)

                        end if

                    end if

                    ! Calculations AFTER the solvers
                    ! A priori, the following variables must be updated: alpha, alpha*rho, and alpha*rho*e now, given
                    ! alpha*rho remains constant in the pT-quilibrium, or the reacting masses have already been updated
                    ! in the pTg-equilibrium, (the remaining) alpha*rho does not need to be updated. I do need, however
                    ! to update alpha, and alpha*rho*e. So, in the end I update alpha, e, and then alpha*rho*e
                    ! entropy
                    sk(1:num_fluids) = &
                        cv(1:num_fluids)*DLOG((TS**g_min(1:num_fluids)) &
                                              /((pS + p_inf(1:num_fluids))**(g_min(1:num_fluids) - 1.0d0))) &
                        + qvp(1:num_fluids)

                    ! enthalpy
                    hk(1:num_fluids) = g_min(1:num_fluids)*cv(1:num_fluids)*TS &
                                       + qv(1:num_fluids)

                    ! GIBBS-FREE ENERGY
                    gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                    ! densities
                    rhok(1:num_fluids) = (pS + p_inf(1:num_fluids)) &
                                         /((g_min(1:num_fluids) - 1)*cv(1:num_fluids)*TS)

                    ! internal energy
                    ek(1:num_fluids) = (pS + g_min(1:num_fluids) &
                                        *p_inf(1:num_fluids))/(pS + p_inf(1:num_fluids)) &
                                       *cv(1:num_fluids)*TS + qv(1:num_fluids)

                    do i = 1, num_fluids

                        ! volume fractions
                        q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = &
                            q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/rhok(i)

                        ! alpha*rho*e
                        q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                            q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*ek(i)

                    end do

                    ! calculating the TOTAL VOLUME FRACTION, reacting mass, and TOTAL MASS, after the PC
                    TvFT = 0.0d0
                    rhoT = 0.0d0
                    do i = 1, num_fluids

                        ! Total volume Fraction Test
                        TvFT = TvFT + q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)

                        ! Total Mass
                        rhoT = rhoT + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)

                    end do

                    rMT = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) &
                          + q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l)

                    ! testing the total volume fraction
                    if (DABS(TvF - TvFT) > 1.0d-2) then

                        print *, 'D total volume fractions AFTER PC: ', TvF - TvFT &
                            , 'j, k, l ', j, k, l

                        print *, 'old', TvF

                        print *, 'new', TvFT

                    end if

                    ! testing the reacting mass
                    if (DABS(rM - rMT) > 1.0d-8) then

                        print *, 'D Reacting Mass AFTER PC: ', rM - rMT &
                            , 'j, k, l ', j, k, l

                    end if

                    ! testing the total mass
                    if (DABS(rho - rhoT) > 1.0d-6) then

                        print *, 'D Total Mass AFTER PC: ', rho - rhoT &
                            , 'j, k, l ', j, k, l

                    end if

                end do
            end do
        end do

    end subroutine s_infinite_relaxation_k ! ----------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! SUBROUTINES SUBROUTINES SUBROUTINES !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine s_infinite_pt_relaxation_k(cv, j, k, l, g_min, MDF, pk, pS, p_inf, p_infA, qv, q_cons_vf, rho, rhoe, TS)
        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: g_min, p_inf, pk, cv, qv
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infA
        real(kind(0.0d0)), intent(IN) :: rho, rhoe
        integer, intent(IN) :: j, k, l
        integer, dimension(num_fluids) :: ig

        real(kind(0.0d0)) :: gp, gpp, hp, pO, mCP, mQ
        character(LEN=2) :: MDF
        integer :: i, ns

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0
        mQ = 0.0d0
        p_infA = p_inf

        ig = 0

        ! Performing tests before initializing the pT-equilibrium
        do i = 1, num_fluids

            ! PRINT *, 'i', i

            ! PRINT *, q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

            ! check if all alpha(i)*rho(i) are nonnegative. If so, abort
            if (q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) < 0.0d0) then

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_inf, pS, qv, (/DABS(pS - pO), DABS(pS - pO)/) &
                                  , rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
&                    . Please, check partial densities. Aborting!')

                ! check which indices I will ignore (no need to abort the solver in this case). Adjust this sgm_eps value for mixture cells
                ! ELSE IF ( ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .GE. 0.0D0 ) &
                !     .AND. ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LE. sgm_eps ) ) THEN

                !     ig( i ) = i

                !     ! this value is rather arbitrary, as I am interested in MINVAL( p_inf for the solver ). This way, I am ensuring
                !     ! this value will not be selected.
                !     p_infA( i ) = 2 * MAXVAL( p_inf )

            end if

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                  *cv(i)*g_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*qv(i)

        end do

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infA)) < 0.0d0) then

            ! PRINT *, 'ERROR IN THE ENERGY CONSTRAINT'

            if ((MDF == 'OV') .or. (MDF == 'SL')) then

                ! PRINT *, 'MDF', MDF

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0d0

                ! temperature
                TS = 0.0d0

                return

            else

                print *, 'MDF', MDF

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_inf, pS, qv, (/DABS(pS - pO), DABS(pS - pO)/) &
                                  , rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
&                    . Please, check energy constraint. Aborting!')

            end if

        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(p_inf)
        ! and infinity, a solution should be able to be found.
        pS = max(maxval(pk), -1.0d0*minval(p_infA)) + 1.0d4

        ! Newton Solver for the pT-equilibrium
        ns = 0
        do while ((DABS(pS - pO) > palpha_eps) .and. (DABS((pS - pO)/pO) > palpha_eps) &
                  .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0d0
            gp = 0.0d0
            hp = 0.0d0

            do i = 1, num_fluids

                ! given pS always change, I need ig( i ) and gp to be in here, as it dynamically updates.
                ! not that I do not need to use p_infA here, but I will do it for consistency
                if (i /= ig(i)) then

                    gp = gp + (g_min(i) - 1.0d0) &
                         *q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*cv(i) &
                         *(rhoe + pS - mQ)/(mCP*(pS + p_infA(i)))

                    gpp = gpp + (g_min(i) - 1.0d0) &
                          *q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*cv(i) &
                          *(p_infA(i) - rhoe + mQ)/(mCP*(pS + p_infA(i))**2)

                end if

            end do

            hp = 1.0d0/(rhoe + pS - mQ) + 1.0d0/(pS + minval(p_infA))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0d0 - gp)/gpp)/(1.0d0 - (1.0d0 - gp + DABS(1.0d0 - gp)) &
                                          /(2.0d0*gpp)*hp)

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
            ! I am not an expert in numerical methods though)
            if (pS <= -1.0d0*minval(p_infA)) then

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infA, pS, qv, (/pS - pO, pS + pO/) &
                                  , rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
&                    . Please, check the pressure value. Aborting!')

            elseif (ieee_is_nan(pS)) then

                call s_mpi_abort('Solver for the pT-relaxation returned NaN values     &
&                    (m_phase_change, s_infinite_pt_relaxation_k). Please, check the error.  &
&                    Aborting!')

                ! cheching is the maximum number of iterations was reached
            elseif (ns > 1e4) then

                call s_mpi_abort('Maximum number of iterations for the pT-relaxation solver reached &
&                    (m_phase_change, s_infinite_pt_relaxation_k). Please, check the pressure value. Aborting!')

            end if

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_pt_relaxation_k ! -----------------------

    subroutine s_infinite_ptg_relaxation_k(cv, g_min, j, k, l, pS, p_inf, p_infpT, rhoe, qv, qvp, q_cons_vf, TS)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: g_min, p_inf, p_infpT, cv, qv, qvp
        real(kind(0.0d0)), intent(IN) :: rhoe
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infA
        real(kind(0.0d0)), dimension(2, 2) :: Jac, InvJac, TJac
        real(kind(0.0d0)), dimension(2) :: R2D, DeltamP
        real(kind(0.0d0)) :: Om
        real(kind(0.0d0)) :: mCP, mCPD, mQ, mQD, mCVGP, mCVGP2

        !< Generic loop iterators
        integer :: i, ns

        ! pTg-equilibrium solution procedure
        ! Newton Sover parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0d-4

        ! Dummy guess to start the pTg-equilibrium problem.
        ! improve this initial condition
        R2D(1) = 0.0d0
        R2D(2) = 0.0d0
        DeltamP(1) = 0.0d0
        DeltamP(2) = 0.0d0
        mQ = 0.0d0

        p_infA = p_infpT

        if (((pS < 0.0d0) .and. ((q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) &
                                  + q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l)) > ((rhoe &
                                                                              - g_min(lp)*p_inf(lp)/(g_min(lp) - 1))/qv(i)))) .or. &
            ((pS >= 0.0d0) .and. (pS < 1.0d-1))) then

            ! improve this initial condition
            pS = 1.0d4

        end if

        ! Loop until the solution for F(X) is satisfied
        ! Check whether I need to use both absolute and relative values
        ! for the residual, and how to do it adequately.
        do while (((DSQRT(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((DSQRT(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1e4))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variable to help in the calculation of the residue
            mCP = 0.0d0
            mCPD = 0.0d0
            mCVGP = 0.0d0
            mCVGP2 = 0.0d0
            mQ = 0.0d0
            mQD = 0.0d0

            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            do i = 1, num_fluids

                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                      *cv(i)*g_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*qv(i)

                ! These auxiliary variables now need to be updated, as the partial densities now
                ! vary at every iteration
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                            *cv(i)*(g_min(i) - 1)/(pS + p_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                             *cv(i)*(g_min(i) - 1)/((pS + p_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*qv(i)

                    ! sum of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*cv(i) &
                           *g_min(i)

                end if

                ! IF ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .LT. 1.0D-10 ) THEN

                !     p_infA( i ) = 2 * MAXVAL( p_inf )

                ! ELSE

                !     p_infA( i ) = p_inf( i )

                ! END IF

            end do

            ! Checking pressure and energy criteria for the (pT) solver to find a solution
            if ((pS <= -1.0d0*minval(p_inf)) &
                .or. ((rhoe - mQ - minval(p_inf)) < 0.0d0)) then

                call s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_inf, pS, qv &
                                  , R2D, rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pTg-relaxation failed &
&                    (m_phase_change, s_infinite_ptg_relaxation_k). Either pS is out of bounds. &
&                    or the energy constraint has been violated. Aborting!')

            end if

            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(cv, rhoe, g_min, InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2 &
            , mQD, p_inf, pS, qv, qvp, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0d0*matmul(InvJac, R2D)

            ! checking if the correction in the mass/pressure will lead to negative values for those quantities
            ! If so, adjust the underrelaxation parameter Om
            if (q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) - Om*DeltamP(1) < 0) then
                call s_mpi_abort('crapV')
            end if

            if (q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) + Om*DeltamP(1) < 0) then
                call s_mpi_abort('crapL')
            end if

            ! updating two reacting 'masses'. Recall that the other 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = &
                q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = &
                q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(cv, g_min, j, k, l, mCPD, mCVGP, mQD, p_inf, qv, qvp, q_cons_vf, pS, rhoe, R2D)

            ! checking if the residue returned any NaN values
            if ((ieee_is_nan(R2D(1))) .or. (ieee_is_nan(R2D(2)))) then

                call s_mpi_abort('Solver for the pTg-relaxation returned NaN values &
&                    (m_phase_change, s_infinite_ptg_relaxation_k). Please, check the error. &
&                    Aborting!')

                ! checking if maximum number of iterations was reached
            elseif ((ns > 1e6)) then

                call s_mpi_abort('Maximum number of iterations reached for the (m_phase_change &
&                    , s_infinite_ptg_relaxation_k). Aborting!')

            end if

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_ptg_relaxation_k ! -----------------------

    !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial den/sity, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
    !> @{
    subroutine s_mixture_volume_fraction_correction(q_cons_vf, j, k, l)
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: sum_alpha
        !> @}
        integer :: i           !< Generic loop iterators

        sum_alpha = 0.0d0

        do i = 1, num_fluids

            ! I changed this line to .LE. so the phase cannot 'disappear'.
            ! Think of more appropriate conditions, later on
            if ((q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) < sgm_eps) .or. &
                (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) < sgm_eps)) then

                q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) = sgm_eps

                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = sgm_eps

                q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = 0.0d0

            end if

            if (q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) > 1.0d0) then

                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = 1.0d0

            end if

            sum_alpha = sum_alpha + q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)

        end do

        do i = 1, num_fluids

            q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l) = &
                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)/sum_alpha

        end do

    end subroutine s_mixture_volume_fraction_correction

    subroutine s_correct_partial_densities(g_min, MCT, mixM, pS, p_inf, q_cons_vf, rhoe, rM, j, k, l)
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: MCT
        real(kind(0.0d0)), intent(IN) :: mixM, rhoe, rM
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: pS, g_min, p_inf
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: rho
        real(kind(0.0d0)) :: TVF, mCP, mQ
        !> @}
        integer :: i, AR           !< Generic loop iterators

        if (rM < 0.0d0) then
            call s_mpi_abort('total reacting mass is negative on "s_correct_partial_densities". Aborting!')
        end if

        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) < 0.0d0) then

            ! PRINT *, 'lp deficient'

            ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

            ! PRINT *, q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

        else if (q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) < 0.0d0) then

            ! PRINT *, 'vp deficient'

            ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

            ! PRINT *, q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )

            q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

            q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = MCT*rM

        end if

    end subroutine s_correct_partial_densities

    !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
    !> @{
    subroutine s_mixture_total_energy_correction(q_cons_vf, j, k, l)
        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: rho, dyn_pres, E_We
        real(kind(0d0)) :: gamma, pi_inf, pres_relax
        real(kind(0d0)) :: qv
        real(kind(0d0)) :: prT, rhoeT

        !> @}
        integer :: i     !< Generic loop iterators
        call s_convert_to_mixture_variables(q_cons_vf, j, k, l, rho, &
                                            gamma, pi_inf, qv)
        dyn_pres = 0.d0
        do i = mom_idx%beg, mom_idx%end
            ! PRINT *, 'i from s_mixture_tec', i
            dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                       q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do
        E_We = 0.d0
        pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - pi_inf - E_We)/gamma

        rhoeT = 0.0d0
        do i = 1, num_fluids
            rhoeT = rhoeT + q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)
        end do

        pRT = (rhoeT - pi_inf)/gamma

        ! PRINT *, 'P_RELAX', pres_relax
        ! PRINT *, 'pRT', pRT
        do i = 1, num_fluids
            ! PRINT *, i
            ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
            q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)* &
                (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf) + &
                q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv
            ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
        end do
        ! ==================================================================
    end subroutine s_mixture_total_energy_correction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO PHASE PRESSURE FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface
    function f_alpha1_prelax(p_k, alpha_k)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)) :: pstar, f_alpha1_prelax
        real(kind(0d0)), dimension(num_fluids), intent(IN) :: p_k, alpha_k
        real(kind(0d0)) :: Z1, Z2, pI, C1, C2
        real(kind(0d0)) :: ap, bp, dp
        ! Calculating coefficients, Eq. C.6, Pelanti 2014
        Z1 = g_minI(lp)*(p_k(1) + p_infI(lp))
        Z2 = g_minI(vp)*(p_k(2) + p_infI(vp))
        pI = (Z2*p_k(1) + Z1*p_k(2))/(Z1 + Z2)
        C1 = 2.d0*g_minI(lp)*p_infI(lp) + (g_minI(lp) - 1.d0)*p_k(1)
        C2 = 2.d0*g_minI(vp)*p_infI(vp) + (g_minI(vp) - 1.d0)*p_k(2)
        ap = 1.d0 + g_minI(vp)*alpha_k(1) + g_minI(lp)*alpha_k(2)
        bp = C1*alpha_k(2) + C2*alpha_k(1) - (g_minI(vp) + 1.d0)*alpha_k(1)*p_k(1) - (g_minI(lp) + 1.d0)*alpha_k(2)*p_k(2)
        dp = -(C2*alpha_k(1)*p_k(1) + C1*alpha_k(2)*p_k(2))
        ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
        pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
        f_alpha1_prelax = alpha_k(1)*((g_minI(lp) - 1.d0)*pstar + 2.d0*p_k(1) + C1)/((g_minI(lp) + 1.d0)*pstar + C1)
    end function f_alpha1_prelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO PHASE PT FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface
    function f_alpha1_ptrelax(rhoalpha1, rhoalpha2, E0)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)) :: pstar, f_alpha1_ptrelax
        real(kind(0d0)), intent(IN) :: rhoalpha1, rhoalpha2, E0
        real(kind(0d0)) :: cv1, cv2, q1, q2
        real(kind(0d0)) :: ap, bp, dp
        cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv; 
        cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv; 
        ! Calculating coefficients, Eq. C.6, Pelanti 2014
        ap = rhoalpha1*cv1 + rhoalpha2*cv2
        bp = q1*cv1*(g_minI(lp) - 1.d0)*rhoalpha1*rhoalpha1 + q2*cv2*(g_minI(vp) - 1.d0)*rhoalpha2*rhoalpha2 + &
             rhoalpha1*cv1*(g_minI(lp)*p_infI(lp) + p_infI(vp)) + rhoalpha2*cv2*(g_minI(vp)*p_infI(vp) + p_infI(lp)) + &
             rhoalpha1*rhoalpha2*(q1*cv2*(g_minI(vp) - 1.d0) + q2*cv1*(g_minI(lp) - 1.d0)) - &
             E0*(cv1*(g_minI(lp) - 1.d0)*rhoalpha1 + cv2*(g_minI(vp) - 1.d0)*rhoalpha2)
        dp = q1*cv1*(g_minI(lp) - 1.d0)*p_infI(vp)*rhoalpha1*rhoalpha1 + q2*cv2*(g_minI(vp) - 1.d0)*p_infI(lp)*rhoalpha2*rhoalpha2 + &
             p_infI(lp)*p_infI(vp)*(rhoalpha1*cv1*g_minI(lp) + rhoalpha2*cv2*g_minI(vp)) + &
             rhoalpha1*rhoalpha2*(q1*cv2*(g_minI(vp) - 1.d0)*p_infI(lp) + q2*cv1*(g_minI(lp) - 1.d0)*p_infI(vp)) - &
             E0*(cv1*(g_minI(lp) - 1.d0)*p_infI(vp)*rhoalpha1 + cv2*(g_minI(vp) - 1.d0)*p_infI(lp)*rhoalpha2)
        ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
        pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
        f_alpha1_ptrelax = (cv1*(g_minI(lp) - 1.d0)*(pstar + p_infI(vp))*rhoalpha1)/ &
                           (cv1*(g_minI(lp) - 1.d0)*(pstar + p_infI(vp))*rhoalpha1 + &
                            cv2*(g_minI(vp) - 1.d0)*(pstar + p_infI(lp))*rhoalpha2)

    end function f_alpha1_ptrelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! SATURATION TEMPERATURE FUNCTIONS !!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_fdfTsat(fp, dfdp, pstar, Tstar)

        real(kind(0d0)), intent(OUT) :: fp, dfdp
        real(kind(0d0)), intent(IN) :: pstar, Tstar
        fp = A + B/Tstar + C*DLOG(Tstar) - &
             DLOG((pstar + p_infI(vp))/(pstar + p_infI(lp))**D)
        dfdp = -B/(Tstar*Tstar) + C/Tstar

    end subroutine s_compute_fdfTsat !-------------------------------

    !>     The purpose of this subroutine is to determine the bracket of
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_Tsat_bracket(TstarA, TstarB, pressure)
        !> @name In-subroutine variables: vapor and liquid material properties
        !> @{
        real(kind(0d0)), intent(OUT) :: TstarA, TstarB
        real(kind(0d0)), intent(IN) :: pressure
        real(kind(0d0)) :: fA, fB, dfdp, factor
        ! Finding lower bound, getting the bracket
        factor = 20.d0
        TstarA = TsatLv
        TstarB = TsatLv + factor
        call s_compute_fdfTsat(fA, dfdp, pressure, TstarA)
        call s_compute_fdfTsat(fB, dfdp, pressure, TstarB)
        !PRINT *, 'fA :: ',fA,', fB :: ',fB
        do while (fA*fB > 0.d0)
            if (TstarA > TsatHv) then
                print *, 'Tsat bracketing failed to find lower bound'
                print *, 'TstarA :: ', TstarA, ', pressure :: ', pressure
                print *, 'fA :: ', fA, ', fB :: ', fB
                call s_mpi_abort()
            end if
            fA = fB
            TstarA = TstarB
            TstarB = TstarA + factor
            call s_compute_fdfTsat(fB, dfdp, pressure, TstarB)
            !PRINT *, 'fB :: ',fB,', TstarB :: ',TstarB,', p :: ',pressure
            !PRINT *,'fB :: ',fB,', TstarB :: ',TstarB
            if (ieee_is_nan(fB)) then
                fB = fA
                TstarB = TstarA
                factor = factor - 10.d0
            else
                factor = 20.d0
            end if
        end do
    end subroutine s_compute_Tsat_bracket

    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface
    function f_Tsat(pressure)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)), intent(IN) :: pressure
        real(kind(0d0)) :: f_Tsat, Tstar
        real(kind(0d0)) :: delta, delta_old, fp, dfdp
        real(kind(0d0)) :: fL, fH, TstarL, TstarH, TsatA, TsatB
        integer :: iter                !< Generic loop iterators
        call s_compute_Tsat_bracket(TsatA, TsatB, pressure)
        ! Computing f at lower and higher end of the bracket
        call s_compute_fdfTsat(fL, dfdp, pressure, TsatA)
        call s_compute_fdfTsat(fH, dfdp, pressure, TsatB)
        ! Establishing the direction of the descent to find zero
        if (fL < 0.d0) then
            TstarL = TsatA; TstarH = TsatB; 
        else
            TstarL = TsatA; TstarH = TsatB; 
        end if
        Tstar = 0.5d0*(TstarL + TstarH)
        delta_old = DABS(TstarH - TstarL)
        delta = delta_old
        call s_compute_fdfTsat(fp, dfdp, pressure, Tstar)
        ! Combining bisection and newton-raphson methods
        do iter = 0, newton_iter
            if ((((Tstar - TstarH)*dfdp - fp)*((Tstar - TstarL)*dfdp - fp) > 0.d0) & ! Bisect if Newton out of range,
                .or. (DABS(2.0*fp) > DABS(delta_old*dfdp))) then         ! or not decreasing fast enough.
                delta_old = delta
                delta = 0.5d0*(TstarH - TstarL)
                Tstar = TstarL + delta
                if (delta == 0.d0) exit
            else                    ! Newton step acceptable, take it
                delta_old = delta
                delta = fp/dfdp
                Tstar = Tstar - delta
                if (delta == 0.d0) exit
            end if
            if (DABS(delta/Tstar) < pTsatnewton_eps) exit
            call s_compute_fdfTsat(fp, dfdp, pressure, Tstar)
            if (fp < 0.d0) then     !Maintain the bracket on the root
                TstarL = Tstar
            else
                TstarH = Tstar
            end if
            if (iter == newton_iter) then
                print *, 'Tsat : ', Tstar, ', iter : ', iter
                print *, 'Tsat did not converge, stopping code'
                call s_mpi_abort()
            end if
        end do
        f_Tsat = Tstar
    end function f_Tsat !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! TWO-PHASE PTG RELAXATION FUNCTIONS !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_ptg_fdf(fp, dfdp, pstar, Tstar, rho0, E0)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)), intent(IN) :: pstar, rho0, E0
        real(kind(0d0)), intent(OUT) :: fp, dfdp, Tstar
        real(kind(0d0)) :: cv1, cv2, q1, q2
        real(kind(0d0)) :: ap, bp, dp
        real(kind(0d0)) :: dadp, dbdp, dddp
        real(kind(0d0)) :: dTdp

        cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv; 
        cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv; 
        ! Calculating coefficients, Eq. C.6, Pelanti 2014
        ap = rho0*cv1*cv2*((g_minI(vp) - 1.d0)*(pstar + g_minI(lp)*p_infI(lp)) - (g_minI(lp) - 1.d0)*(pstar + g_minI(vp)*p_infI(vp)))
        bp = E0*((g_minI(lp) - 1.d0)*cv1*(pstar + p_infI(vp)) - (g_minI(vp) - 1.d0)*cv2*(pstar + p_infI(lp))) + &
             rho0*((g_minI(vp) - 1.d0)*cv2*q1*(pstar + p_infI(lp)) - (g_minI(lp) - 1.d0)*cv1*q2*(pstar + p_infI(vp))) + &
             cv2*(pstar + p_infI(lp))*(pstar + g_minI(vp)*p_infI(vp)) - cv1*(pstar + p_infI(vp))*(pstar + g_minI(lp)*p_infI(lp))
        dp = (q2 - q1)*(pstar + p_infI(lp))*(pstar + p_infI(vp))
        ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
        Tstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
        ! PRINT *, 'ap: ', ap, 'bp: ', bp, 'dp:', dp
        ! PRINT *, 'sqrt(D)', DSQRT(bp*bp - 4.d0*ap*dp)
        ! PRINT *, 'Tstar: ', Tstar
        ! PRINT *, 'Tstar T : ', (-bp - DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
        ! Calculating the derivatives wrt pressure of the coefficients
        dadp = rho0*cv1*cv2*((g_minI(vp) - 1.d0) - (g_minI(lp) - 1.d0))
        dbdp = E0*((g_minI(lp) - 1.d0)*cv1 - (g_minI(vp) - 1.d0)*cv2) + &
               rho0*((g_minI(vp) - 1.d0)*cv2*q1 - (g_minI(lp) - 1.d0)*cv1*q2) + &
               cv2*((pstar + p_infI(lp)) + (pstar + g_minI(vp)*p_infI(vp))) - &
               cv1*((pstar + p_infI(vp)) + (pstar + g_minI(lp)*p_infI(lp)))
        dddp = (q2 - q1)*((pstar + p_infI(lp)) + (pstar + p_infI(vp)))
        ! Derivative of the temperature wrt to pressure, needed for dfdp
        dTdp = (-dbdp + (0.5d0/DSQRT(bp*bp - 4.d0*ap*dp))*(2.d0*bp*dbdp - &
                                                           4.d0*(ap*dddp + dp*dadp)))/(2.d0*ap) - (dadp/ap)*Tstar
        fp = A + B/Tstar + C*DLOG(Tstar) + &
             D*DLOG(pstar + p_infI(lp)) - DLOG(pstar + p_infI(vp))
        dfdp = -B/(Tstar*Tstar)*dTdp + C/Tstar*dTdp + &
               D/(pstar + p_infI(lp)) - 1.d0/(pstar + p_infI(vp))

        ! PRINT *, Tstar
        ! PRINT *, dadp, dbdp, dddp, dTdp, fp, dfdp

    end subroutine s_compute_ptg_fdf !------------------------

    !>     The purpose of this subroutine is to determine the bracket of
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_ptg_bracket(pstarA, pstarB, pstar, rho0, E0)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)), intent(OUT) :: pstarA, pstarB
        real(kind(0d0)), intent(IN) :: pstar, rho0, E0
        real(kind(0d0)) :: fA, fB, dfdp, Tstar
        real(kind(0d0)) :: pS, factor

        pstarA = 0.7d0*pstar
        pstarB = pstar
        ! PRINT *, 's_compute_ptg_fdf A'
        call s_compute_ptg_fdf(fA, dfdp, pstarA, Tstar, rho0, E0)
        ! PRINT *, 's_compute_ptg_fdf B'
        call s_compute_ptg_fdf(fB, dfdp, pstarB, Tstar, rho0, E0)
        factor = 1.05d0
        do while (fA*fB > 0.d0)
            ! PRINT *, 'fA', fA
            ! PRINT *, 'fB', fb
            if (pstarA > 1.d13) then
                print *, 'ptg bracketing failed to find lower bound'
                print *, 'pstarA :: ', pstarA
                call s_mpi_abort()
            end if
            fA = fB
            pstarA = pstarB
            pstarB = pstarA*factor
            call s_compute_ptg_fdf(fB, dfdp, pstarB, Tstar, rho0, E0)
            if (ieee_is_nan(fB)) then
                fB = fA
                pstarB = pstarA
                factor = factor*0.95d0
            else
                factor = 1.05d0
            end if
        end do

    end subroutine s_compute_ptg_bracket

    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_ptg_pTrelax(pstar, Tstar, rho0, E0)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        real(kind(0d0)), intent(INOUT) :: pstar
        real(kind(0d0)), intent(OUT) :: Tstar
        real(kind(0d0)), intent(IN) :: rho0, E0
        real(kind(0d0)) :: pstarA, pstarB, dTstar
        real(kind(0d0)) :: delta, delta_old, fp, dfdp
        real(kind(0d0)) :: fL, fH, pstarL, pstarH
        integer :: iter                !< Generic loop iterators
        ! Computing the bracket of the root solution
        call s_compute_ptg_bracket(pstarA, pstarB, pstar, rho0, E0)
        ! PRINT *, 'pstarA', pstarA, 'pstarB', pstarB
        ! Computing f at lower and higher end of the bracket
        call s_compute_ptg_fdf(fL, dfdp, pstarA, dTstar, rho0, E0)
        ! PRINT *, 'fL', fL, 'dfdpL', dfdp, 'TstarL', dTstar
        ! CALL s_compute_ptg_fdf(fH,dfdp,pstarB,dTstar,rho0,E0)
        ! PRINT *, 'fH', fH, 'dfdpH', dfdp, 'TstarH', dTstar
        ! Establishing the direction of the descent to find zero
        if (fL < 0.d0) then
            pstarL = pstarA; pstarH = pstarB; 
        else
            pstarL = pstarB; pstarH = pstarA; 
        end if
        pstar = 0.5d0*(pstarA + pstarB)
        delta_old = DABS(pstarB - pstarA)
        delta = delta_old
        call s_compute_ptg_fdf(fp, dfdp, pstar, dTstar, rho0, E0)
        ! Combining bisection and newton-raphson methods
        do iter = 0, newton_iter
            if ((((pstar - pstarH)*dfdp - fp)*((pstar - pstarL)*dfdp - fp) > 0.d0) & ! Bisect if Newton out of range,
                .or. (DABS(2.0*fp) > DABS(delta_old*dfdp))) then         ! or not decreasing fast enough.
                delta_old = delta
                delta = 0.5d0*(pstarH - pstarL)
                pstar = pstarL + delta
                if (delta == 0.d0) exit                    ! Change in root is negligible
            else                                              ! Newton step acceptable, take it
                delta_old = delta
                delta = fp/dfdp
                pstar = pstar - delta
                if (delta == 0.d0) exit
            end if
            if (DABS(delta/pstar) < ptgnewton_eps) exit           ! Stopping criteria
            ! Updating to next iteration
            call s_compute_ptg_fdf(fp, dfdp, pstar, Tstar, rho0, E0)
            if (fp < 0.d0) then !Maintain the bracket on the root
                pstarL = pstar
            else
                pstarH = pstar
            end if
        end do
    end subroutine s_compute_ptg_pTrelax !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!  K-PHASE P RELAXATION FUNCTIONS !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_pk_fdf(fp, dfdp, pstar, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        ! Cell-average conservative variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0d0)), intent(OUT) :: fp, dfdp
        real(kind(0d0)), intent(IN) :: pstar
        real(kind(0d0)), dimension(num_fluids), intent(OUT) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids), intent(IN) :: gamma_min, pres_inf, pres_K_init
        real(kind(0d0)) :: numerator, denominator
        real(kind(0d0)) :: drhodp
        integer, intent(IN) :: j, k, l
        integer :: i
        fp = -1.d0; dfdp = 0.d0; 
        do i = 1, num_fluids
            numerator = gamma_min(i)*(pstar + pres_inf(i))
            denominator = numerator + pres_K_init(i) - pstar
            rho_K_s(i) = q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/ &
                         max(q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l), sgm_eps)*numerator/denominator
            drhodp = q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/ &
                     max(q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l), sgm_eps)* &
                     gamma_min(i)*(pres_K_init(i) + pres_inf(i))/(denominator*denominator)
            fp = fp + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/rho_K_s(i)
            dfdp = dfdp - q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)* &
                   drhodp/(rho_K_s(i)*rho_K_s(i))
        end do

    end subroutine s_compute_pk_fdf !------------------------

    !>     The purpose of this subroutine is to determine the bracket of
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_pk_bracket(pstarA, pstarB, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0d0)), intent(OUT) :: pstarA, pstarB
        real(kind(0d0)) :: fA, fB, dfdp
        real(kind(0d0)) :: factor
        real(kind(0d0)), dimension(num_fluids), intent(IN) :: gamma_min, pres_inf, pres_K_init
        real(kind(0d0)), dimension(num_fluids), intent(OUT) :: rho_K_s
        integer, intent(IN) :: j, k, l
        !pstarA = 100.d0
        !pstarB = 5.d2
        pstarA = 1.d-15
        pstarB = 1.d1
        call s_compute_pk_fdf(fA, dfdp, pstarA, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        call s_compute_pk_fdf(fB, dfdp, pstarB, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        factor = 10.d0
        do while (fA*fB > 0.d0)
            if (pstarA > 1.d13) then
                print *, 'P-K bracketing failed to find lower bound'
                print *, 'pstarA :: ', pstarA, ', pstarB :: ', pstarB
                print *, 'alpha1 :: ', q_cons_vf(adv_idx%beg)%sf(j, k, l)
                print *, 'alpha2 :: ', q_cons_vf(adv_idx%beg)%sf(j, k, l)
                call s_mpi_abort()
            end if
            fA = fB
            pstarA = pstarB
            pstarB = pstarA*factor
            call s_compute_pk_fdf(fB, dfdp, pstarB, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
            if (ieee_is_nan(fB)) then
                fB = fA
                pstarB = pstarA
                factor = factor*0.5d0
            else
                factor = 10.d0
            end if
        end do
    end subroutine s_compute_pk_bracket

    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_compute_p_relax_k(rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
        !> @{
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0d0)) :: pstar, pstarA, pstarB
        real(kind(0d0)) :: delta, delta_old, fp, dfdp
        real(kind(0d0)) :: fL, fH, pstarL, pstarH
        real(kind(0d0)), dimension(num_fluids), intent(IN) :: gamma_min, pres_inf, pres_K_init
        real(kind(0d0)), dimension(num_fluids), intent(OUT) :: rho_K_s
        integer, intent(IN) :: j, k, l
        integer :: iter                !< Generic loop iterators
        ! Computing the bracket of the root solution
        call s_compute_pk_bracket(pstarA, pstarB, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        ! Computing f at lower and higher end of the bracket
        call s_compute_pk_fdf(fL, dfdp, pstarA, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        call s_compute_pk_fdf(fH, dfdp, pstarB, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        ! Establishing the direction of the descent to find zero
        if (fL < 0.d0) then
            pstarL = pstarA; pstarH = pstarB; 
        else
            pstarL = pstarB; pstarH = pstarA; 
        end if
        pstar = 0.5d0*(pstarA + pstarB)
        delta_old = DABS(pstarB - pstarA)
        delta = delta_old
        call s_compute_pk_fdf(fp, dfdp, pstar, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
        ! Combining bisection and newton-raphson methods
        do iter = 0, newton_iter
            if ((((pstar - pstarH)*dfdp - fp)*((pstar - pstarL)*dfdp - fp) > 0.d0) & ! Bisect if Newton out of range,
                .or. (DABS(2.0*fp) > DABS(delta_old*dfdp))) then         ! or not decreasing fast enough.
                delta_old = delta
                delta = 0.5d0*(pstarH - pstarL)
                pstar = pstarL + delta
                if (delta == 0.d0) exit                    ! Change in root is negligible
            else                                              ! Newton step acceptable, take it
                delta_old = delta
                delta = fp/dfdp
                pstar = pstar - delta
                if (delta == 0.d0) exit
            end if
            if (DABS(delta/pstar) < pknewton_eps) exit           ! Stopping criteria
            ! Updating to next iteration
            call s_compute_pk_fdf(fp, dfdp, pstar, rho_K_s, gamma_min, pres_inf, pres_K_init, q_cons_vf, j, k, l)
            if (fp < 0.d0) then !Maintain the bracket on the root
                pstarL = pstar
            else
                pstarH = pstar
            end if
        end do
    end subroutine s_compute_p_relax_k !-------------------------------

    ! This SUBROUTINE IS USED TO CALCULATE THE 2X2 JACOBIAN AND ITS INVERSE FOR THE PROPOSED
    ! pTg-Equilibrium procedure
    subroutine s_compute_jacobian_matrix(cv, rhoe, g_min, InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2 &
        , mQD, p_inf, pS, qv, qvp, q_cons_vf, TJac)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mCVGP2, mQD
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: g_min, p_inf, cv, qv, qvp
        real(kind(0.0d0)), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: ml, mT, TS, dFdT, dTdm, dTdp
        integer :: i

        ! mass of the reactant liquid
        ml = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) &
             + q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l)

        TS = 1/(mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) &
                + ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                      - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cv(lp)*g_min(lp) - cv(vp)*g_min(vp))*DLOG(TS) &
            - (qvp(lp) - qvp(vp)) &
            + cv(lp)*(g_min(lp) - 1)*DLOG(pS + p_inf(lp)) &
            - cv(vp)*(g_min(vp) - 1)*DLOG(pS + p_inf(vp))

        dTdm = -(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                 - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)))*TS**2

        dTdp = (mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))**2 &
                + ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp))**2 &
                      - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))**2) &
                + mCVGP2)*TS**2

        ! F = (F1,F2) is the fuction whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        Jac(1, 2) = dFdT*dTdp + TS &
                    *(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                      - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)))

        ! dF2dm
        Jac(2, 1) = (qv(vp) - qv(lp) &
                     + (cv(vp)*g_min(vp) - cv(lp)*g_min(lp)) &
                     /(ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                           - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                       + mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) + mCVGP) &
                     - (ml*(cv(vp)*g_min(vp) - cv(lp)*g_min(lp)) &
                        - mT*cv(vp)*g_min(vp) - mCPD) &
                     *(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                       - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                     /((ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                            - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                        + mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) + mCVGP)**2))/1

        Jac(2, 2) = (1 + (ml*(cv(vp)*g_min(vp) - cv(lp)*g_min(lp)) &
                          - mT*cv(vp)*g_min(vp) - mCPD) &
                     *(ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp))**2 &
                           - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))**2) &
                       + mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))**2 + mCVGP2) &
                     /(ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                           - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                       + mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) + mCVGP)**2)/1

        ! elements of J^{-1}
        InvJac(1, 1) = Jac(2, 2)
        InvJac(1, 2) = -1.0d0*Jac(1, 2)
        InvJac(2, 1) = -1.0d0*Jac(2, 1)
        InvJac(2, 2) = Jac(1, 1)

        ! elements of J^{T}
        TJac(1, 1) = Jac(1, 1)
        TJac(1, 2) = Jac(2, 1)
        TJac(2, 1) = Jac(1, 2)
        TJac(2, 2) = Jac(2, 2)

        ! dividing by det(J)
        InvJac = InvJac/(Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1))

    end subroutine s_compute_jacobian_matrix

    subroutine s_compute_pTg_residue(cv, g_min, j, k, l, mCPD, mCVGP, mQD, p_inf, qv, qvp, q_cons_vf, pS, rhoe, R2D)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: g_min, p_inf, cv, qv, qvp
        real(kind(0.0d0)), dimension(2), intent(OUT) :: R2D
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: ml, mT, TS
        integer :: i

        ! mass of the reactant liquid
        ml = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) &
             + q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l)

        TS = 1/(mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) &
                + ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                      - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cv(lp)*g_min(lp) - cv(vp)*g_min(vp)) &
                     *(1 - DLOG(TS)) - (qvp(lp) - qvp(vp)) &
                     + cv(lp)*(g_min(lp) - 1)*DLOG(pS + p_inf(lp)) &
                     - cv(vp)*(g_min(vp) - 1)*DLOG(pS + p_inf(vp))) &
                 + qv(lp) - qv(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + ml*(qv(vp) - qv(lp)) - mT*qv(vp) - mQD &
                  + (ml*(g_min(vp)*cv(vp) - g_min(lp)*cv(lp)) &
                     - mT*g_min(vp)*cv(vp) - mCPD) &
                  /(ml*(cv(lp)*(g_min(lp) - 1)/(pS + p_inf(lp)) &
                        - cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp))) &
                    + mT*cv(vp)*(g_min(vp) - 1)/(pS + p_inf(vp)) + mCVGP))/1

    end subroutine s_compute_pTg_residue

    ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
    subroutine s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_inf, pS, qv, R2D, rhoe, q_cons_vf, TS) ! ----------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), dimension(2, 2), intent(IN) :: Jac, InvJac
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_inf, qv
        real(kind(0.0d0)), dimension(2), intent(IN) :: R2D, DeltamP
        real(kind(0.0d0)), intent(IN) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe, mQ
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: rho
        !< Generic loop iterator
        integer :: i

        print *, 'j, k, l', j, k, l

        print *, 'rhoe', rhoe

        print *, 'mQ', mQ

        print *, 'Energy constrain', (rhoe - mQ - minval(p_inf))

        print *, 'R2D', R2D

        print *, 'l2(R2D)', DSQRT(R2D(1)**2 + R2D(2)**2)

        print *, 'DeltamP', DeltamP

        print *, 'pS', pS

        print *, '-min(p_inf)', -minval(p_inf)

        print *, 'TS', TS

        do i = 1, num_fluids

            rho = rho + q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)

            print *, 'rho', rho

        end do

        do i = 1, num_fluids

            print *, 'i', i

            print *, 'alpha_i', q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)

            print *, 'alpha_rho_i', q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)

            print *, 'mq_i', q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l) &
                *qv(i)

            print *, 'rho', rho

            print *, 'internal energies', q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l)

            print *, 'Y_i', q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)/rho

        end do

        print *, 'J', Jac, 'J-1', InvJac

    end subroutine s_tattletale

    ! Newton Solver for the finding the Saturation temperature TSat for a given saturation pressure
    subroutine s_TSat(cv, g_min, pSat, p_inf, qv, qvp, TSat, TSIn)

        real(kind(0.0d0)), intent(OUT) :: TSat
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: cv, g_min, p_inf, qv, qvp
        real(kind(0.0d0)), intent(IN) :: pSat, TSIn
        real(kind(0.0d0)) :: dFdT, FT

        ! Generic loop iterators
        integer :: ns

        if ((pSat == 0.0d0) .and. (TSIn == 0.0d0)) then

            ! assigning Saturation temperature
            TSat = 0.0d0

        else

            ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            TSat = TSIn

            ! iteration counter
            ns = 0

            ! DO WHILE ( ( DABS( FT ) .GT. ptgalpha_eps .AND. DABS( FT / TSat ) .GT. ptgalpha_eps ) &
            do while ((DABS(FT) > ptgalpha_eps) &
                      .or. (ns == 0))

                ! increasing counter
                ns = ns + 1

                ! calculating residual
                ! FT = A + B / TSat + C * DLOG( TSat ) &
                ! + D * DLOG( ( pSat + p_inf( lp ) ) ) - DLOG( pSat + p_inf( vp ) )

                FT = TSat*((cv(lp)*g_min(lp) - cv(vp)*g_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvp(lp) - qvp(vp)) &
                           + cv(lp)*(g_min(lp) - 1)*DLOG(pSat + p_inf(lp)) &
                           - cv(vp)*(g_min(vp) - 1)*DLOG(pSat + p_inf(vp))) &
                     + qv(lp) - qv(vp)

                ! calculating the jacobian
                ! dFdT = - B / ( TSat ** 2) + C / TSat

                dFdT = &
                    -(cv(lp)*g_min(lp) - cv(vp)*g_min(vp))*DLOG(TSat) &
                    - (qvp(lp) - qvp(vp)) &
                    + cv(lp)*(g_min(lp) - 1)*DLOG(pSat + p_inf(lp)) &
                    - cv(vp)*(g_min(vp) - 1)*DLOG(pSat + p_inf(vp))

                ! updating saturation temperature
                TSat = TSat - FT/dFdT

                ! Checking if TSat returns a NaN
                if (ieee_is_nan(TSat)) then

                    call s_mpi_abort('TSat returned NaN values, when it should not (by assumption &
&                        of first order transition - m_phase_change, s_TSat). Aborting!')

                    ! checking if the maximum number of iterations has been reached
                elseif (ns > 1e4) then

                    call s_mpi_abort('Maximum number of iterations reached for TSat &
&                            (m_phase_change, s_TSat). Aborting!')
                end if

            end do

        end if

    end subroutine s_TSat

    subroutine s_finalize_relaxation_solver_module()

        s_relaxation_solver => null()

    end subroutine

#endif

end module m_phase_change
