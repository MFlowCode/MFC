!>
!! @file m_phase_change.fpp
!! @brief Contains module m_phasechange

#:include 'macros.fpp'

!> @brief This module is used to relax the model equations (6-eqn model)
!> towards pressure and temperature (6-eqn to 4-eqn), and (if wanted) Gibbs free
!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.
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
        !> @}

    end interface
    !> @}

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8        !< max # of iterations
    real(kind(0d0)), parameter :: pCr = 4.94d7   !< Critical water pressure
    real(kind(0d0)), parameter :: TCr = 385.05 + 273.15  !< Critical water temperature
    real(kind(0d0)), parameter :: mixM = 1.0d-8 !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(kind(0d0)) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)

    procedure(s_abstract_relaxation_solver), pointer :: s_relaxation_solver => null()

contains

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    subroutine s_initialize_phasechange_module()
        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        D = ((gs_min(lp) - 1.0d0)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

    end subroutine s_initialize_phasechange_module !-------------------------------

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(kind(0.0d0)) :: TS, TSOV, TSSL, TSatOV, TSatSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(kind(0.0d0)) :: rhoe, dynE, rhos !< total internal energy, kinetic energy, and total entropy
        real(kind(0.0d0)) :: rho, rM, m1, m2, MCT !< total density, total reacting mass, individual reacting masses
        real(kind(0.0d0)) :: TvF !< total volume fraction

        !$acc declare create(pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)

        real(kind(0d0)), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok

        !< Generic loop iterators
        integer :: i, j, k, l

        !$acc declare create(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok)

        ! starting equilibrium solver
        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok,pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0d0; TvF = 0.0d0
                    !$acc loop seq
                    do i = 1, num_fluids

                        ! Mixture density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! correcting negative (recating) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                    ! fixing m1 and m2 AFTER correcting the partial densities. Note that these values must be stored for the phase
                    ! change process that will happen a posteriori
                    m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0d0
                    !$acc loop seq
                    do i = momxb, momxe

                        dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)**2/rho

                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! at each of the cells. The internal energy is calculated as the total energy minus the kinetic
                    ! energy to preserved its value at sharp interfaces
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                    ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                    call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

                    ! check if pTg-equilibrium is required
                    ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! given the outputs from the pT- and pTg-equilibrium solvers are just p and one of the partial masses
                    ! (pTg- case)
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) > mixM*rM) &
                                                  .and. (q_cons_vf(vp + contxb - 1)%sf(j, k, l) > mixM*rM)) &
                        .and. (pS < pCr) .and. (TS < TCr)) then

                        ! Checking if phase change is needed, by checking whether the final solution is either subcoooled
                        ! liquid or overheated vapor.

                        ! overheated vapor case
                        ! depleting the mass of liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! tranferring the total mass to vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                        call s_infinite_pt_relaxation_k(j, k, l, 0, pSOV, p_infOV, rM, q_cons_vf, rhoe, TSOV)

                        ! calculating Saturation temperature
                        call s_TSat(pSOV, TSatOV, TSOV)

                        ! subcooled liquid case
                        ! tranferring the total mass to liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! depleting the mass of vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! calling pT-equilibrium for subcooled liquid, which is MFL = 1
                        call s_infinite_pt_relaxation_k(j, k, l, 1, pSSL, p_infSL, rM, q_cons_vf, rhoe, TSSL)

                        ! calculating Saturation temperature
                        call s_TSat(pSSL, TSatSL, TSSL)

                        ! checking the conditions for overheated vapor and subcooled liquide
                        if (TSOV > TSatOV) then

                            ! Assigning pressure
                            pS = pSOV

                            ! Assigning Temperature
                            TS = TSOV

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        elseif (TSSL < TSatSL) then

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
                            call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

                        end if

                    end if

                    ! Calculations AFTER equilibrium

                    ! entropy
                    sk(1:num_fluids) = cvs(1:num_fluids)*DLOG((TS**gs_min(1:num_fluids)) &
                                                              /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0d0))) + qvps(1:num_fluids)

                    ! enthalpy
                    hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                                       + qvs(1:num_fluids)

                    ! Gibbs-free energy
                    gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                    ! densities
                    rhok(1:num_fluids) = (pS + ps_inf(1:num_fluids)) &
                                         /((gs_min(1:num_fluids) - 1)*cvs(1:num_fluids)*TS)

                    ! internal energy
                    ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                                        *ps_inf(1:num_fluids))/(pS + ps_inf(1:num_fluids)) &
                                       *cvs(1:num_fluids)*TS + qvs(1:num_fluids)

                    ! calculating volume fractions, internal energies, and total entropy
                    rhos = 0.0d0
                    !$acc loop seq
                    do i = 1, num_fluids

                        ! volume fractions
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

                        ! alpha*rho*e
                        q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)

                        ! Total entropy
                        rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                    end do
                end do
            end do
        end do

    end subroutine s_infinite_relaxation_k ! ----------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infpT
        real(kind(0.0d0)), intent(IN) :: rM, rhoe
        integer, intent(IN) :: j, k, l, MFL
        real(kind(0.0d0)), dimension(num_fluids) :: pk !< individual initial pressures
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium

        real(kind(0.0d0)) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0; mQ = 0.0d0; p_infpT = ps_inf; 
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infpT)) < 0.0d0) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0d0

                ! temperature
                TS = 0.0d0

                return
            end if

        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0d4

        ! Newton Solver for the pT-equilibrium
        ns = 0
        ! change this relative error metric. 1E4 is just arbitrary
        do while ((DABS(pS - pO) > palpha_eps) .and. (DABS((pS - pO)/pO) > palpha_eps/1e4) .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0d0; gp = 0.0d0; hp = 0.0d0
            !$acc loop seq
            do i = 1, num_fluids

                gp = gp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                     *(rhoe + pS - mQ)/(mCP*(pS + p_infpT(i)))

                gpp = gpp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                      *(p_infpT(i) - rhoe + mQ)/(mCP*(pS + p_infpT(i))**2)

            end do

            hp = 1.0d0/(rhoe + pS - mQ) + 1.0d0/(pS + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0d0 - gp)/gpp)/(1.0d0 - (1.0d0 - gp + DABS(1.0d0 - gp)) &
                                          /(2.0d0*gpp)*hp)
        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_pt_relaxation_k ! -----------------------

    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for N fluids under pT
        !!      and 2 fluids under pTg-equilibrium. There is a final common p and T during relaxation
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(kind(0.0d0)), dimension(2, 2) :: Jac, InvJac, TJac !< matrices for the Newton Solver
        real(kind(0.0d0)), dimension(2) :: R2D, DeltamP !< residual and correction array
        real(kind(0.0d0)) :: Om ! underrelaxation factor
        real(kind(0.0d0)) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure
        ! Newton Solver parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0d-3

        p_infpTg = p_infpT

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
        ! Dummy guess to start the pTg-equilibrium problem.
        ! improve this initial condition
        R2D(1) = 0.0d0; R2D(2) = 0.0d0
        DeltamP(1) = 0.0d0; DeltamP(2) = 0.0d0
        do while (((DSQRT(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((DSQRT(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1d6))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0d0; mCPD = 0.0d0; mCVGP = 0.0d0; mCVGP2 = 0.0d0; mQ = 0.0d0; mQD = 0.0d0
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            !$acc loop seq
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

            end do

            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0d0*matmul(InvJac, R2D)

            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP
    end subroutine s_infinite_ptg_relaxation_k ! -----------------------

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param MCT partial density correction parameter
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)
        !$acc routine seq

        !> @name variables for the correction of the reacting partial densities
        !> @{
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(INOUT) :: rM
        real(kind(0.0d0)), intent(OUT) :: MCT
        integer, intent(IN) :: j, k, l
        !> @}
        if (rM < 0.0d0) then

            if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM) .and. &
                (q_cons_vf(vp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM)) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0.0d0

                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0.0d0

                rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            end if

        end if

        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

        elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if
    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine calculates the 2 x 2 Jacobian and, its inverse and transpose
        !!      to be used in the pTg-equilibirium procedure
        !!  @param InvJac Inverse of the Jacobian Matrix
        !!  @param j generic loop iterator for x direction
        !!  @param Jac Jacobian Matrix
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mCVGP2 auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press^2
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TJac Transpose of the Jacobian Matrix
    subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, mCPD, mCVGP, mCVGP2
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        real(kind(0.0d0)) :: ml, mT, TS, dFdT, dTdm, dTdp ! mass of the reacting fluid, total reacting mass, and auxiliary variables

        ! mass of the reacting liquid
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

        ! F = (F1,F2) is the function whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the Gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        ! dF1dp
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
        ! dF2dp
        Jac(2, 2) = (1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                          - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

        ! intermediate elements of J^{-1}
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

    !>  This auxiliary subroutine computes the residue of the pTg-equilibrium procedure
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mQD sum of the total alpha*rho*qv
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param pS equilibrium pressure at the interface
        !!  @param rhoe mixture energy
        !!  @param R2D (2D) residue array
    subroutine s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2), intent(OUT) :: R2D
        real(kind(0.0d0)) :: ml, mT, TS !< mass of the reacting liquid, total reacting mass, equilibrium temperature

        ! mass of the reacting liquid
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

    !>  This auxiliary subroutine finds the Saturation temperature for a given
        !!      saturation pressure through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param TSIn equilibrium Temperature
    subroutine s_TSat(pSat, TSat, TSIn)
        !$acc routine seq

        real(kind(0.0d0)), intent(OUT) :: TSat
        real(kind(0.0d0)), intent(IN) :: pSat, TSIn
        real(kind(0.0d0)) :: dFdT, FT, Om !< auxiliary variables

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

            ! underrelaxation factor
            Om = 1.0d-3
            do while ((DABS(FT) > ptgalpha_eps) .or. (ns == 0))
                ! increasing counter
                ns = ns + 1

                ! calculating residual
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

            end do

        end if

    end subroutine s_TSat

    !>  This subroutine finalizes the phase change module
    subroutine s_finalize_relaxation_solver_module()
    end subroutine

#endif

end module m_phase_change
