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
 s_infinite_relaxation_k, &
 s_finalize_relaxation_solver_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the infinite relaxation solver
        !> @{
        subroutine s_abstract_relaxation_solver(q_cons_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        end subroutine

    end interface

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: newton_iter = 50        !< p_relaxk \alpha iter,                set to 25
    real(kind(0d0)), parameter :: pres_crit = 22.09d6   !< Critical water pressure              set to 22.06d6
    real(kind(0d0)), parameter :: T_crit = 647.29d0  !< Critical water temperature           set to 648
    integer, parameter         :: lp = 1    ! index for the liquid phase of the reacting fluid
    integer, parameter         :: vp = 2    ! index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{b
    real(kind(0d0)) :: A, B, C, D
    !> @}

    !$acc declare create(newton_iter,pres_crit,T_crit,lp,vp,A,B,C,D)

    real(kind(0d0)), allocatable, dimension(:) :: p_infA, pk, sk, hk, gk, ek, rhok

    !$acc declare create(p_infA, pk, sk, hk, gk, ek, rhok)

    procedure(s_abstract_relaxation_solver), &
        pointer :: s_relaxation_solver => null()

contains

    !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface
    subroutine s_initialize_phasechange_module()

        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvs(vp) - qvs(lp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        D = ((gs_min(lp) - 1.0d0)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        allocate (p_infA(1:num_fluids), pk(1:num_fluids), sk(1:num_fluids), &
        hk(1:num_fluids), gk(1:num_fluids), ek(1:num_fluids), rhok(1:num_fluids))

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        ! Opening and writing the header of the run-time information file
        if ((relax_model == 5) .or. (relax_model == 6)) then
            s_relaxation_solver => s_infinite_relaxation_k
        else
            call s_mpi_abort('relaxation solver was not set!')
        end if

    end subroutine s_initialize_phasechange_module !-------------------------------

    ! This subroutine is created to activate either the pT- or the pTg-equilibrium
    ! model, with changes such that mass depletion is taken into consideration
    subroutine s_infinite_relaxation_k(q_cons_vf) ! ----------------
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        ! real(kind(0.0d0)), dimension(num_fluids) :: p_infA
        ! real(kind(0.0d0)), dimension(num_fluids) :: pk, sk, hk, gk
        ! real(kind(0.0d0)), dimension(num_fluids) :: ek, rhok
        real(kind(0.0d0)) :: pS, pSOV, pSSL
        real(kind(0.0d0)) :: TS, TSOV, TSatOV, TSatSL, TSSL
        real(kind(0.0d0)) :: rhoe, mCP, mQ, TvF
        real(kind(0.0d0)) :: rho, rM, m1, m2, MCT, mixM, rhoT, rMT
        real(kind(0.0d0)) :: TvFT, rhos
        real(kind(0.0d0)) :: dynE, VFT

        !< Generic loop iterators
        integer :: i, j, k, l
        integer :: ns

        ! threshold for 'mixture cell'. If Y < mixM, the cell is not considered a mixture for pTg purposes
        mixM = 1.0d-08

!        !$acc declare create(p_infA, pk, sk, hk, gk, ek, rhok, pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, mCP, mQ, TvF, rho, rM, m1, m2, MCT, mixM, rhoT, rMT, TvFT, rhos, dynE, VFT)

        ! starting equilibrium solver
! !$acc parallel loop collapse(3) gang vector default(present)
        ! copy(p_infA, pk, sk, hk, gk, ek, rhok, pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, mCP, mQ, TvF, rho, rM, m1, m2, MCT, mixM, rhoT, rMT, TvFT, rhos, dynE, VFT)
        !$acc kernels
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0d0; TvF = 0.0d0
                    do i = 1, num_fluids

                        ! pressure that comes from the homogeneous solver at that cell
                        pk(i) = (gs_min(i) - 1)*(q_cons_vf(i + intxb - 1)%sf(j, k, l) &
                                  - q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)) &
                                /q_cons_vf(i + advxb - 1)%sf(j, k, l) &
                                - gs_min(i)*ps_inf(i)

                        ! Mixture density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    rM = m1 + m2

                    ! correcting negative volume and mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, mixM, q_cons_vf, rM, j, k, l)

                    ! kinetic energy so as to calculate the total internal energy as the total energy minus the
                    ! kinetic energy.
                    dynE = 0.0d0
                    do i = momxb, momxe

                        dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)**2/rho

                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! at each of the cells. Note I calculate UE as TE - KE due to numerical reasons (stability at the dicontinuities)
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                    ! for this case, MFL cannot be either 0 or 1
                    call s_infinite_pt_relaxation_k(j, k, l, 2, pk, pS &
                                                    , p_infA, q_cons_vf, rho, rhoe, TS)

                    ! check if pTg-equilibrium is required
                    ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! given the outputs from the pT- and pTg-equilibrium solvers are just p, T (byproduct)
                    ! , and the partial masses (pTg- case)
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                                                   > MixM*rM) .or. (q_cons_vf(vp + contxb - 1)%sf(j, k, l) &
                                                                    > MixM*rM)) .and. (pS < pres_crit) .and. (TS < T_crit)) then

                        ! Checking if phase change can happen.
                        ! overheated vapor
                        ! correcting the liquid partial density
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! correcting the vapor partial density
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                        call s_infinite_pt_relaxation_k(j, k, l, 0, pk, pSOV &
                                                        , p_infA, q_cons_vf, rho, rhoe, TSOV)

                        ! calculating Saturation temperature
                        call s_TSat(pSOV, TSatOV, TSOV)

                        ! subcooled liquid
                        ! correcting the liquid partial density
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! correcting the vapor partial density
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! calling pT-equilibrium for subcooled liquid, which is MFL = 1
                        call s_infinite_pt_relaxation_k(j, k, l, 1, pk, pSSL &
                                                        , p_infA, q_cons_vf, rho, rhoe, TSSL)

                        ! calculating Saturation temperature
                        call s_TSat(pSSL, TSatSL, TSSL)

                        ! Maybe I will have to change this for .GE. instead, since mass depletion can
                        ! happen at pTg-eq. Make sure it works, first
                        if (TSOV > TSatOV) then

                            ! Assigning pressure
                            pS = pSOV

                            ! Assigning Temperature
                            TS = TSOV

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        else if (TSSL < TSatSL) then

                            ! Assigning pressure
                            pS = pSSL

                            ! Assigning Temperature
                            TS = TSSL

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        else

                            ! returning partial pressures to what they were from the homogeneous solver
                            ! liquid
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m1

                            ! vapor
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m2

                            ! calling the pTg-equilibrium solver
                            call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infA, rhoe, q_cons_vf, TS)

                        end if

                    end if

                    ! Calculations AFTER the solvers
                    ! A priori, the following variables must be updated: alpha, alpha*rho, and alpha*rho*e now, given
                    ! alpha*rho remains constant in the pT-quilibrium, or the reacting masses have already been updated
                    ! in the pTg-equilibrium, (the remaining) alpha*rho does not need to be updated. I do need, however
                    ! to update alpha, and alpha*rho*e. So, in the end I update alpha, e, and then alpha*rho*e
                    ! entropy
                    sk(1:num_fluids) = cvs(1:num_fluids)*DLOG((TS**gs_min(1:num_fluids)) &
                        /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0d0))) + qvps(1:num_fluids)

                    ! enthalpy
                    hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                                       + qvs(1:num_fluids)

                    ! GIBBS-FREE ENERGY
                    gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                    ! densities
                    rhok(1:num_fluids) = (pS + ps_inf(1:num_fluids)) &
                                         /((gs_min(1:num_fluids) - 1)*cvs(1:num_fluids)*TS)

                    ! internal energy
                    ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                                        *ps_inf(1:num_fluids))/(pS + ps_inf(1:num_fluids)) &
                                       *cvs(1:num_fluids)*TS + qvs(1:num_fluids)

                    ! calculating the TOTAL VOLUME FRACTION, reacting mass, and TOTAL MASS, after the PC
                    TvFT = 0.0d0; rhoT = 0.0d0; rhos = 0.0d0
                    do i = 1, num_fluids

                        ! volume fractions
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

                        ! alpha*rho*e
                        q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)

                        ! Total volume Fraction Test
                        TvFT = TvFT + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        ! Total mixture density Test
                        rhoT = rhoT + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Total entropy
                        rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                    end do

                    rMT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)
#ifndef MFC_OpenACC
                    ! testing the total volume fraction
                    if (DABS(TvF - TvFT) > 1.0d-2) then

                        print *, 'D total volume fractions AFTER PC: ', TvF - TvFT &
                            , 'j, k, l ', j, k, l

                        print *, 'old', TvF, 'new', TvFT

                    end if

                    ! testing the reacting mass
                    if (DABS(rM - rMT) > 1.0d-8) then

                        print *, 'D Reacting Mass AFTER PC: ', rM - rMT &
                            , 'j, k, l ', j, k, l

                        print *, 'old', rM, 'new', rMT

                    end if

                    ! testing the total mass
                    if (DABS(rho - rhoT) > 1.0d-6) then

                        print *, 'D Total Mass AFTER PC: ', rho - rhoT &
                            , 'j, k, l ', j, k, l

                        print *, 'old', rho, 'new', rhoT

                    end if
#endif
                end do
            end do
        end do
        !$acc end kernels

    end subroutine s_infinite_relaxation_k ! ----------------

    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pk, pS, p_infA, q_cons_vf, rho, rhoe, TS)
!$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: pk 
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infA
        real(kind(0.0d0)), intent(IN) :: rho, rhoe
        integer, intent(IN) :: j, k, l, MFL
        integer, dimension(num_fluids) :: ig

        real(kind(0.0d0)) :: gp, gpp, hp, pO, mCP, mQ
        integer :: i, ns

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0; mQ = 0.0d0; p_infA = ps_inf

        ig = 0

        ! Performing tests before initializing the pT-equilibrium
        do i = 1, num_fluids

            ! check if all alpha(i)*rho(i) are nonnegative. If so, abort
#ifndef MFC_OpenACC
            if (q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0.0d0) then

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infA, pS, (/DABS(pS - pO), DABS(pS - pO)/) &
                                  , rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                &                    . Please, check partial densities. Aborting!')

                ! check which indices I will ignore (no need to abort the solver in this case). Adjust this sgm_eps value for mixture cells
                ! ELSE IF ( ( q_cons_vf( i + contxb - 1 )%sf( j, k, l ) .GE. 0.0D0 ) &
                !     .AND. ( q_cons_vf( i + contxb - 1 )%sf( j, k, l ) .LE. sgm_eps ) ) THEN

                !     ig( i ) = i

                !     ! this value is rather arbitrary, as I am interested in MINVAL( ps_inf for the solver ). This way, I am ensuring
                !     ! this value will not be selected.
                !     p_infA( i ) = 2 * MAXVAL( ps_inf )

            end if
#endif
            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                  *cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infA)) < 0.0d0) then

            ! PRINT *, 'ERROR IN THE ENERGY CONSTRAINT'

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0d0

                ! temperature
                TS = 0.0d0

                return
#ifndef MFC_OpenACC
            else

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infA, pS, (/DABS(pS - pO), DABS(pS - pO)/) &
                                  , rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
&                    . Please, check energy constraint. Aborting!')

#endif
            end if

        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
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
                ! if (i /= ig(i)) then

                    gp = gp + (gs_min(i) - 1.0d0) &
                         *q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                         *(rhoe + pS - mQ)/(mCP*(pS + p_infA(i)))

                    gpp = gpp + (gs_min(i) - 1.0d0) &
                          *q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                          *(p_infA(i) - rhoe + mQ)/(mCP*(pS + p_infA(i))**2)

                ! end if

            end do

            hp = 1.0d0/(rhoe + pS - mQ) + 1.0d0/(pS + minval(p_infA))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0d0 - gp)/gpp)/(1.0d0 - (1.0d0 - gp + DABS(1.0d0 - gp)) &
                                          /(2.0d0*gpp)*hp)

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
#ifndef MFC_OpenACC                                          
            if (pS <= -1.0d0*minval(p_infA)) then

                call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                  , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infA, pS, (/pS - pO, pS + pO/) &
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
#endif

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_pt_relaxation_k ! -----------------------

    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)
!$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
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
        R2D(1) = 0.0d0; R2D(2) = 0.0d0
        DeltamP(1) = 0.0d0; DeltamP(2) = 0.0d0
        mQ = 0.0d0

        p_infA = p_infpT

        if (((pS < 0.0d0) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                                  + q_cons_vf(vp + contxb - 1)%sf(j, k, l)) > ((rhoe &
                                                                              - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp)))) .or. &
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

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0d0; mCPD = 0.0d0; mCVGP = 0.0d0; mCVGP2 = 0.0d0; mQ = 0.0d0; mQD = 0.0d0

            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            do i = 1, num_fluids

                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                      *cvs(i)*gs_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                ! These auxiliary variables now need to be updated, as the partial densities now
                ! vary at every iteration
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                            *cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                             *cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                    ! sum of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                           *gs_min(i)

                end if

                ! IF ( q_cons_vf( i + advxb - 1 )%sf( j, k, l ) .LT. 1.0D-10 ) THEN

                !     p_infA( i ) = 2 * MAXVAL( ps_inf )

                ! ELSE

                !     p_infA( i ) = ps_inf( i )

                ! END IF

            end do

            ! Checking pressure and energy criteria for the (pT) solver to find a solution
#ifndef MFC_OpenACC            
            if ((pS <= -1.0d0*minval(ps_inf)) &
                .or. ((rhoe - mQ - minval(ps_inf)) < 0.0d0)) then

                call s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_infA, pS &
                                  , R2D, rhoe, q_cons_vf, TS)

                call s_mpi_abort('Solver for the pTg-relaxation failed &
&                    (m_phase_change, s_infinite_ptg_relaxation_k). Either pS is out of bounds. &
&                    or the energy constraint has been violated. Aborting!')

            end if
#endif
            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(rhoe, InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2 &
            , mQD, pS, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0d0*matmul(InvJac, R2D)

            ! checking if the correction in the mass/pressure will lead to negative values for those quantities
            ! If so, adjust the underrelaxation parameter Om

#ifndef MFC_OpenACC
            if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1) < 0) then
                call s_mpi_abort('Liquid phase will become negative. Aborting (for the moment).')
            else if (q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1) < 0) then
                call s_mpi_abort('Vapor phase will become negative. Aborting (for the moment).')
            end if
#endif
            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = &
                q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = &
                q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

            ! checking if the residue returned any NaN values
#ifndef MFC_OpenACC            
            if ((ieee_is_nan(R2D(1))) .or. (ieee_is_nan(R2D(2)))) then

                call s_mpi_abort('Solver for the pTg-relaxation returned NaN values &
&                    (m_phase_change, s_infinite_ptg_relaxation_k). Please, check the error. &
&                    Aborting!')

                ! checking if maximum number of iterations was reached
            elseif ((ns > 1e6)) then

                call s_mpi_abort('Maximum number of iterations reached for the (m_phase_change &
&                    , s_infinite_ptg_relaxation_k). Aborting!')

            end if
#endif

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_ptg_relaxation_k ! -----------------------

    subroutine s_correct_partial_densities(MCT, mixM, q_cons_vf, rM, j, k, l)
!$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: MCT
        real(kind(0.0d0)), intent(IN) :: mixM, rM
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: rho
        real(kind(0.0d0)) :: TVF, mCP, mQ
        !> @}
        integer :: i, AR           !< Generic loop iterators

#ifndef MFC_OpenACC        
        if (rM < 0.0d0) then
            call s_mpi_abort('total reacting mass is negative on "s_correct_partial_densities". Aborting!')
        end if
#endif
        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

        else if (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if

    end subroutine s_correct_partial_densities

    ! This SUBROUTINE IS USED TO CALCULATE THE 2X2 JACOBIAN AND ITS INVERSE FOR THE PROPOSED
    ! pTg-Equilibrium procedure
    subroutine s_compute_jacobian_matrix(rhoe, InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2 &
        , mQD, pS, q_cons_vf, TJac)
!$acc routine seq
        
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mCVGP2, mQD
        real(kind(0.0d0)), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: ml, mT, TS, dFdT, dTdm, dTdp
        integer :: i

        ! mass of the reactant liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TS) &
            - (qvps(lp) - qvps(vp)) &
            + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
            - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))

        dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

        dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                + mCVGP2)*TS**2

        ! F = (F1,F2) is the fuction whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        Jac(1, 2) = dFdT*dTdp + TS &
                    *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))

        ! dF2dm
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

    subroutine s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)
!$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD
        real(kind(0.0d0)), dimension(2), intent(OUT) :: R2D
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: ml, mT, TS
        integer :: i

        ! mass of the reactant liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                     *(1 - DLOG(TS)) - (qvps(lp) - qvps(vp)) &
                     + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))) &
                 + qvs(lp) - qvs(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD &
                  + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                     - mT*gs_min(vp)*cvs(vp) - mCPD) &
                  /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                        - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP))/1

    end subroutine s_compute_pTg_residue

    ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
    subroutine s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_infA, pS, R2D, rhoe, q_cons_vf, TS) ! ----------------
        
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), dimension(2, 2), intent(IN) :: Jac, InvJac
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infA
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

        print *, 'Energy constrain', (rhoe - mQ - minval(p_infA))

        print *, 'R2D', R2D

        print *, 'l2(R2D)', DSQRT(R2D(1)**2 + R2D(2)**2)

        print *, 'DeltamP', DeltamP

        print *, 'pS', pS

        print *, '-min(ps_inf)', -minval(p_infA)

        print *, 'TS', TS

        do i = 1, num_fluids

            rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

        end do

        print *, 'rho', rho

        do i = 1, num_fluids

            print *, 'i', i

            print *, 'alpha_i', q_cons_vf(i + advxb - 1)%sf(j, k, l)

            print *, 'alpha_rho_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)

            print *, 'mq_i', q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                *qvs(i)

            print *, 'internal energies', q_cons_vf(i + intxb - 1)%sf(j, k, l)

            print *, 'Y_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho

        end do

        print *, 'J', Jac, 'J-1', InvJac

    end subroutine s_tattletale

    ! Newton Solver for the finding the Saturation temperature TSat for a given saturation pressure
    subroutine s_TSat(pSat, TSat, TSIn)
!$acc routine seq

        real(kind(0.0d0)), intent(OUT) :: TSat
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
                ! + D * DLOG( ( pSat + ps_inf( lp ) ) ) - DLOG( pSat + ps_inf( vp ) )

                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                ! dFdT = - B / ( TSat ** 2) + C / TSat

                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - FT/dFdT
                
#ifndef MFC_OpenACC
                ! Checking if TSat returns a NaN
                if (ieee_is_nan(TSat)) then

                    call s_mpi_abort('TSat returned NaN values, when it should not (by assumption &
&                        of first order transition - m_phase_change, s_TSat). Aborting!')

                    ! checking if the maximum number of iterations has been reached
                elseif (ns > 1e4) then

                    call s_mpi_abort('Maximum number of iterations reached for TSat &
&                            (m_phase_change, s_TSat). Aborting!')
                end if

#endif
            end do

        end if

    end subroutine s_TSat

    subroutine s_finalize_relaxation_solver_module()

        s_relaxation_solver => null()

        deallocate (p_infA, pk, sk, hk, gk, ek, rhok)
        
    end subroutine

#endif

end module m_phase_change