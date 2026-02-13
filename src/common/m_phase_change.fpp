!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.
#:include 'case.fpp'
#:include 'macros.fpp'

module m_phase_change

#ifndef MFC_POST_PROCESS

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; 
    public :: s_initialize_phasechange_module, &
              s_relaxation_solver, &
              s_infinite_relaxation_k, &
              s_finalize_relaxation_solver_module

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8_wp        !< max # of iterations
    real(wp), parameter :: pCr = 4.94e7_wp   !< Critical water pressure
    real(wp), parameter :: TCr = 385.05_wp + 273.15_wp  !< Critical water temperature
    real(wp), parameter :: mixM = 1.0e-8_wp !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(wp) :: A, B, C, D
    !> @}

    $:GPU_DECLARE(create='[A,B,C,D]')

contains

    !> This subroutine should dispatch to the correct relaxation solver based
        !! It replaces the procedure pointer, which CCE
        !! on.
    impure subroutine s_relaxation_solver(q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        ! is empty because in current master the procedure pointer
        ! never assigned
        @:ASSERT(.false., "s_relaxation_solver called but it currently does nothing")
    end subroutine s_relaxation_solver

    !>  The purpose of this subroutine is to initialize the phase change module
        !! the parameters needed for phase change and
        !! phase change module that will be used
        !! pTg-equilibrium)
    impure subroutine s_initialize_phasechange_module
        ! used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        D = ((gs_min(lp) - 1.0_wp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

    end subroutine s_initialize_phasechange_module

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !! fluids for g-equilibrium)
        !! considering mass depletion, depending on the incoming
        !!
        !! Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(wp) :: TS, TSOV, TSSL, TSatOV, TSatSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(wp) :: rhoe, dynE, rhos !< total internal energy, kinetic energy, and total entropy
        real(wp) :: rho, rM, m1, m2, MCT !< total density, total reacting mass, individual reacting masses
        real(wp) :: TvF !< total volume fraction

        !
        !
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok
        #:else
            real(wp), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok
        #:endif
        !

        !< Generic loop iterators
        integer :: i, j, k, l

        ! equilibrium solver
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i,j,k,l,p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok,pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF]')
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0_wp; TvF = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids

                        ! density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Volume Fraction
                        TvF = TvF + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                    end do

                    ! the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! the phase-change process.
                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! negative (reacting) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                    ! m1 and m2 AFTER correcting the partial densities. Note that these values must be stored for the phase
                    ! process that will happen a posteriori
                    m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = momxb, momxe

                        dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2/rho

                    end do

                    ! the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! each of the cells. The internal energy is calculated as the total energy minus the kinetic
                    ! to preserved its value at sharp interfaces
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                    ! this case, MFL cannot be either 0 or 1, so I chose it to be 2
                    call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, q_cons_vf, rhoe, TS)

                    ! if pTg-equilibrium is required
                    ! that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! the outputs from the pT- and pTg-equilibrium solvers are just p and one of the partial masses
                    ! case)
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) > mixM*rM) &
                                                  .and. (q_cons_vf(vp + contxb - 1)%sf(j, k, l) > mixM*rM)) &
                        .and. (pS < pCr) .and. (TS < TCr)) then

                        ! if phase change is needed, by checking whether the final solution is either subcoooled
                        ! or overheated vapor.

                        ! vapor case
                        ! the mass of liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! the total mass to vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        ! pT-equilibrium for overheated vapor, which is MFL = 0
                        call s_infinite_pt_relaxation_k(j, k, l, 0, pSOV, p_infOV, q_cons_vf, rhoe, TSOV)

                        ! Saturation temperature
                        call s_TSat(pSOV, TSatOV, TSOV)

                        ! liquid case
                        ! the total mass to liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        ! the mass of vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! pT-equilibrium for subcooled liquid, which is MFL = 1
                        call s_infinite_pt_relaxation_k(j, k, l, 1, pSSL, p_infSL, q_cons_vf, rhoe, TSSL)

                        ! Saturation temperature
                        call s_TSat(pSSL, TSatSL, TSSL)

                        ! the conditions for overheated vapor and subcooled liquide
                        if (TSOV > TSatOV) then

                            ! pressure
                            pS = pSOV

                            ! Temperature
                            TS = TSOV

                            ! the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                            ! the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        elseif (TSSL < TSatSL) then

                            ! pressure
                            pS = pSSL

                            ! Temperature
                            TS = TSSL

                            ! the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                            ! the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        else

                            ! partial pressures to what they were from the homogeneous solver
                            !
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m1

                            !
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m2

                            ! the pTg-equilibrium solver
                            call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

                        end if

                    end if

                    ! AFTER equilibrium

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        !
                        sk(i) = cvs(i)*log((TS**gs_min(i)) &
                                           /((pS + ps_inf(i))**(gs_min(i) - 1.0_wp))) + qvps(i)

                        !
                        hk(i) = gs_min(i)*cvs(i)*TS &
                                + qvs(i)

                        ! energy
                        gk(i) = hk(i) - TS*sk(i)

                        !
                        rhok(i) = (pS + ps_inf(i)) &
                                  /((gs_min(i) - 1)*cvs(i)*TS)

                        ! energy
                        ek(i) = (pS + gs_min(i) &
                                 *ps_inf(i))/(pS + ps_inf(i)) &
                                *cvs(i)*TS + qvs(i)
                    end do

                    ! volume fractions, internal energies, and total entropy
                    rhos = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids

                        ! fractions
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

                        !
                        if (model_eqns == 3) then
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)
                        end if

                        ! entropy
                        rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_infinite_relaxation_k

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !! generic loop iterator for x direction
        !! generic loop iterator for y direction
        !! generic loop iterator for z direction
        !! flag that tells whether the fluid is gas (0), liquid (1), or a mixture (2)
        !! equilibrium pressure at the interface
        !! stiffness for the participating fluids under pT-equilibrium
        !! sum of the reacting masses
        !! Cell-average conservative variables
        !! mixture energy
        !! equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, q_cons_vf, rhoe, TS)
        $:GPU_ROUTINE(function_name='s_infinite_pt_relaxation_k', &
            & parallelism='[seq]', cray_inline=True)

        ! variables
        integer, intent(in) :: j, k, l, MFL
        real(wp), intent(out) :: pS
        real(wp), dimension(1:), intent(out) :: p_infpT
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), intent(in) :: rhoe
        real(wp), intent(out) :: TS
        real(wp) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver
        real(wp) :: p_infpT_sum

        integer :: i, ns !< generic loop iterators

        ! variables for the pT-equilibrium solver
        mCP = 0.0_wp; mQ = 0.0_wp; p_infpT_sum = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            p_infpT(i) = ps_inf(i)
            p_infpT_sum = p_infpT_sum + abs(p_infpT(i))
        end do
        ! tests before initializing the pT-equilibrium
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids

            ! of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            if (num_fluids < 3) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = num_fluids + 1, 3
                    p_infpT(i) = p_infpT_sum
                end do
            end if
        #:endif

        ! energy constraint
        if ((rhoe - mQ - minval(p_infpT)) < 0.0_wp) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! zero values for mass depletion cases
                !
                pS = 0.0_wp

                !
                TS = 0.0_wp

                return
            end if

        end if

        ! initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! over the Newton's solver
        pO = 0.0_wp

        ! improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! infinity, a solution should be able to be found.
        pS = 1.0e4_wp

        ! Solver for the pT-equilibrium
        ns = 0
        ! this relative error metric. 1.e4_wp is just arbitrary
        do while ((abs(pS - pO) > palpha_eps) .and. (abs((pS - pO)/pO) > palpha_eps/1.e4_wp) .or. (ns == 0))

            ! counter
            ns = ns + 1

            ! old pressure
            pO = pS

            ! functions used in the Newton's solver
            gpp = 0.0_wp; gp = 0.0_wp; hp = 0.0_wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids

                gp = gp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                     *(rhoe + pS - mQ)/(mCP*(pS + p_infpT(i)))

                gpp = gpp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                      *(p_infpT(i) - rhoe + mQ)/(mCP*(pS + p_infpT(i))**2)

            end do

            hp = 1.0_wp/(rhoe + pS - mQ) + 1.0_wp/(pS + minval(p_infpT))

            ! common pressure for the newton solver
            pS = pO + ((1.0_wp - gp)/gpp)/(1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp)) &
                                           /(2.0_wp*gpp)*hp)
        end do

        ! temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_pt_relaxation_k

    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for N fluids under pT
        !! fluids under pTg-equilibrium. There is a final common p and T during relaxation
        !! generic loop iterator for x direction
        !! generic loop iterator for y direction
        !! generic loop iterator for z direction
        !! equilibrium pressure at the interface
        !! stiffness for the participating fluids under pT-equilibrium
        !! mixture energy
        !! Cell-average conservative variables
        !! equilibrium temperature at the interface
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)
        $:GPU_ROUTINE(function_name='s_infinite_ptg_relaxation_k', &
            & parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: j, k, l
        real(wp), intent(inout) :: pS
        real(wp), dimension(1:), intent(in) :: p_infpT
        real(wp), intent(in) :: rhoe
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout) :: TS
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        #:else
            real(wp), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        #:endif
        real(wp), dimension(2, 2) :: Jac, InvJac, TJac !< matrices for the Newton Solver
        real(wp), dimension(2) :: R2D, DeltamP !< residual and correction array
        real(wp) :: Om ! underrelaxation factor
        real(wp) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver
        real(wp) :: ml, mT, dFdT, dTdm, dTdp

        !< Generic loop iterators
        integer :: i, ns
        ! solution procedure
        ! Solver parameters
        !
        ns = 0

        ! factor
        Om = 1.0e-3_wp

        p_infpTg = p_infpT

        if (((pS < 0.0_wp) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                                   + q_cons_vf(vp + contxb - 1)%sf(j, k, l)) > ((rhoe &
                                                                                 - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp)))) .or. &
            ((pS >= 0.0_wp) .and. (pS < 1.0e-1_wp))) then

            ! this initial condition
            pS = 1.0e4_wp

        end if

        ! until the solution for F(X) is satisfied
        ! whether I need to use both absolute and relative values
        ! the residual, and how to do it adequately.
        ! guess to start the pTg-equilibrium problem.
        ! this initial condition
        R2D(1) = 0.0_wp; R2D(2) = 0.0_wp
        DeltamP(1) = 0.0_wp; DeltamP(2) = 0.0_wp
        do while (((sqrt(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((sqrt(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1.e6_wp))) &
                  .or. (ns == 0))

            ! counter for the iterative procedure
            ns = ns + 1

            ! variables to help in the calculation of the residue
            mCP = 0.0_wp; mCPD = 0.0_wp; mCVGP = 0.0_wp; mCVGP2 = 0.0_wp; mQ = 0.0_wp; mQD = 0.0_wp
            ! must be updated through the iterations, as they either depend on
            ! partial masses for all fluids, or on the equilibrium pressure
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids

                ! of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                      *cvs(i)*gs_min(i)

                ! of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                ! auxiliary variables now need to be updated, as the partial densities now
                ! at every iteration
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                            *cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                             *cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                    ! of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                           *gs_min(i)

                end if

            end do

            ! the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model

            ! of the reacting liquid
            ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

            ! of the two participating fluids
            mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                 + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                    + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                          - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mCVGP)

            dFdT = &
                -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TS) &
                - (qvps(lp) - qvps(vp)) &
                + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
                - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))

            dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

            dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                    + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                          - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                    + mCVGP2)*TS**2

            ! = (F1,F2) is the function whose roots we are looking for
            ! = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
            ! = 0 is the Gibbs free energy quality
            ! = 0 is the enforcement of the thermodynamic (total - kinectic) energy
            !
            Jac(1, 1) = dFdT*dTdm

            !
            Jac(1, 2) = dFdT*dTdp + TS &
                        *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                          - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))

            !
            Jac(2, 1) = (qvs(vp) - qvs(lp) &
                         + (cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                         /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                               - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                           + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP) &
                         - (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                            - mT*cvs(vp)*gs_min(vp) - mCPD) &
                         *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                         /((ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                                - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                            + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2))/1
            !
            Jac(2, 2) = (1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                              - mT*cvs(vp)*gs_min(vp) - mCPD) &
                         *(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                               - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                           + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                         /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                               - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                           + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

            ! elements of J^{-1}
            InvJac(1, 1) = Jac(2, 2)
            InvJac(1, 2) = -1.0_wp*Jac(1, 2)
            InvJac(2, 1) = -1.0_wp*Jac(2, 1)
            InvJac(2, 2) = Jac(1, 1)

            ! of J^{T}
            TJac(1, 1) = Jac(1, 1)
            TJac(1, 2) = Jac(2, 1)
            TJac(2, 1) = Jac(1, 2)
            TJac(2, 2) = Jac(2, 2)

            ! by det(J)
            InvJac = InvJac/(Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1))

            ! correction array for Newton's method
            DeltamP = -1.0_wp*(matmul(InvJac, R2D))

            ! two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            !
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            !
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! pressure
            pS = pS + Om*DeltamP(2)

            ! residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! (ii) the energy before and after the phase-change process.

            ! of the reacting liquid
            ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

            ! of the two participating fluids
            mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                 + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                    + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                          - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mCVGP)

            ! Free Energy Equality condition (DG)
            R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                         *(1 - log(TS)) - (qvps(lp) - qvps(vp)) &
                         + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
                         - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

            ! Energy Process condition (DE)
            R2D(2) = (rhoe + pS &
                      + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD &
                      + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                         - mT*gs_min(vp)*cvs(vp) - mCPD) &
                      /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                            - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                        + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP))/1

        end do

        ! temperature
        TS = (rhoe + pS - mQ)/mCP
    end subroutine s_infinite_ptg_relaxation_k

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !! sum is positive. Inert phases are not corrected at this moment
        !! partial density correction parameter
        !! Cell-average conservative variables
        !! sum of the reacting masses
        !! generic loop iterator for x direction
        !! generic loop iterator for y direction
        !! generic loop iterator for z direction
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)
        $:GPU_ROUTINE(function_name='s_correct_partial_densities', &
            & parallelism='[seq]', cray_inline=True)

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), intent(out) :: MCT
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout) :: rM
        integer, intent(in) :: j, k, l
        !> @}
        if (rM < 0.0_wp) then

            if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) >= -1.0_wp*mixM) .and. &
                (q_cons_vf(vp + contxb - 1)%sf(j, k, l) >= -1.0_wp*mixM)) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0.0_wp

                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0.0_wp

                rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            end if

        end if

        ! the correction in terms of an absolute value might not be the best practice.
        ! a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0_wp) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM

        elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0_wp) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if
    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine finds the Saturation temperature for a given
        !! through a newton solver
        !! Saturation Pressure
        !! Saturation Temperature
        !! equilibrium Temperature
    elemental subroutine s_TSat(pSat, TSat, TSIn)
        $:GPU_ROUTINE(function_name='s_TSat',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(in) :: pSat
        real(wp), intent(out) :: TSat
        real(wp), intent(in) :: TSIn

        real(wp) :: dFdT, FT, Om !< auxiliary variables

        ! loop iterators
        integer :: ns

        if ((f_approx_equal(pSat, 0.0_wp)) .and. (f_approx_equal(TSIn, 0.0_wp))) then

            ! Saturation temperature
            TSat = 0.0_wp

        else

            ! initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! over the Newton's solver
            TSat = TSIn

            ! counter
            ns = 0

            ! factor
            Om = 1.0e-3_wp
            do while ((abs(FT) > ptgalpha_eps) .or. (ns == 0))
                ! counter
                ns = ns + 1

                ! residual
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - log(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! the jacobian
                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))

                ! saturation temperature
                TSat = TSat - Om*FT/dFdT

            end do

        end if

    end subroutine s_TSat

    !>  This subroutine finalizes the phase change module
    impure subroutine s_finalize_relaxation_solver_module
    end subroutine s_finalize_relaxation_solver_module

#endif

end module m_phase_change
