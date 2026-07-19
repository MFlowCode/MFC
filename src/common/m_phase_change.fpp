!>
!! @file
!! @brief Contains module m_phase_change

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Phase transition relaxation solvers for liquid-vapor flows with cavitation and boiling
module m_phase_change

#ifndef MFC_POST_PROCESS
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use ieee_arithmetic
    use m_helper_basic
    use m_constants, only: model_eqns_6eq

    implicit none

    private
    public :: s_initialize_phasechange_module, s_relaxation_solver, s_infinite_relaxation_k, s_finalize_relaxation_solver_module, &
        & pc_iter_count

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter  :: max_iter = 100000            !< max Newton iterations before accepting the last iterate
    real(wp), parameter :: pCr = 4.94e7_wp              !< Critical pressure of water [Pa]
    real(wp), parameter :: TCr = 385.05_wp + 273.15_wp  !< Critical temperature of water [K]
    integer, parameter  :: ptg_ls_max = 30              !< max backtracking-line-search halvings in the pTg solver
    real(wp), parameter :: mixM = 1.0e-8_wp             !< Mixture mass fraction threshold for triggering phase change
    integer, parameter  :: lp = 1                       !< index for the liquid phase of the reacting fluid
    integer, parameter  :: vp = 2                       !< index for the vapor phase of the reacting fluid
    !> @}

    !> Per-cell Newton-iteration count for the current step; allocated only when relax .and. (load_weight_wrt .or.
    !! sfc_partition_wrt).
    real(stp), allocatable :: pc_iter_count(:,:,:)
    $:GPU_DECLARE(create='[pc_iter_count]')

contains

    !> Dispatch to the correct relaxation solver. Replaces the procedure pointer, which CCE is breaking on.
    impure subroutine s_relaxation_solver(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        ! This is empty because in current master the procedure pointer was never assigned

        @:ASSERT(.false., "s_relaxation_solver called but it currently does nothing")

    end subroutine s_relaxation_solver

    !> Initialize the phase change module (no module-level state to set up; the pT/pTg relaxation solvers are self-contained)
    impure subroutine s_initialize_phasechange_module

#ifdef MFC_SIMULATION
        ! the load-weight field is computed for load_weight_wrt AND for sfc_partition_wrt (which
        ! calls s_compute_load_weight too): allocate and populate under the same condition, else
        ! the sfc-only path reads unallocated (or allocated-but-never-written) iteration counts
        if (relax .and. (load_weight_wrt .or. sfc_partition_wrt)) then
            @:ALLOCATE(pc_iter_count(0:m, 0:n, 0:p))
        end if
#endif

    end subroutine s_initialize_phasechange_module

    !> Apply pT- or pTg-equilibrium relaxation with mass depletion based on the incoming state conditions.
    subroutine s_infinite_relaxation_k(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: pS                    !< equilibrium pressure
        real(wp) :: TS                    !< equilibrium temperature
        real(wp) :: rhoe, dynE, rhos      !< total internal energy, kinetic energy, and total entropy
        real(wp) :: rho, rM, m1, m2, MCT  !< total density, total reacting mass, individual reacting masses
        real(wp) :: TvF                   !< total volume fraction
        ! $:GPU_DECLARE(create='[pS,TS,rhoe,dynE,rhos,rho,rM,m1,m2,MCT,TvF]')

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: p_infpT, sk, hk, gk, ek, rhok
        #:else
            real(wp), dimension(num_fluids) :: p_infpT, sk, hk, gk, ek, rhok
        #:endif
        ! $:GPU_DECLARE(create='[p_infpT,sk,hk,gk,ek,rhok]')

        !> Generic loop iterators
        integer :: i, j, k, l
        integer :: ns_pc, ns_tmp   !< per-cell Newton-iteration accumulators for load-weight diagnostic
        logical :: count_pc_iters  !< host-evaluated kernel guard (a device copy of the wrt flags is never synced)

#ifdef _CRAYFTN
#ifdef MFC_OpenACC
        ! CCE 19 IPA workaround: prevent bring_routine_resident SIGSEGV DIR$ NOINLINE s_infinite_pt_relaxation_k DIR$ NOINLINE
        ! s_infinite_ptg_relaxation_k DIR$ NOINLINE s_correct_partial_densities
#endif
#endif

        ! starting equilibrium solver

#ifdef MFC_SIMULATION
        count_pc_iters = load_weight_wrt .or. sfc_partition_wrt
#else
        count_pc_iters = .false.
#endif

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, p_infpT, sk, hk, gk, ek, rhok, pS, TS, rhoe, dynE, rhos, rho, rM, &
                            & m1, m2, MCT, TvF, ns_pc, ns_tmp]', copyin='[count_pc_iters]')
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    rho = 0.0_wp; TvF = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! Mixture density
                        rho = rho + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) + q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)

                    ! correcting negative (reacting) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                    ! fixing m1 and m2 AFTER correcting the partial densities. Note that these values must be stored for the phase
                    ! change process that will happen a posteriori
                    m1 = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2/rho
                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures at each
                    ! of the cells. The internal energy is calculated as the total energy minus the kinetic energy to preserved its
                    ! value at sharp interfaces
                    rhoe = q_cons_vf(eqn_idx%E)%sf(j, k, l) - dynE

                    ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium for this
                    ! case, MFL cannot be either 0 or 1, so I chose it to be 2
                    call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, q_cons_vf, rhoe, TS, ns_pc)

                    ! Check if pTg-equilibrium needed; only partial densities require updating
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, &
                        & l) > mixM*rM) .and. (q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, &
                        & l) > mixM*rM)) .and. (pS < pCr) .and. (TS < TCr)) then
                        ! Solve pTg-equilibrium directly on the actual reacting masses. The Newton solver projects
                        ! the liquid mass onto [0, mT], so it recovers the single-phase limits itself (ml -> 0 for
                        ! all-vapor, ml -> mT for all-liquid). The former overheated-vapor / subcooled-liquid pT
                        ! shortcuts were removed: their pT states differ O(1) from the pTg equilibrium, so the
                        ! sub-ULP shortcut/pTg branch decision flipped across backends (CPU vs GPU) near a phase
                        ! boundary and destroyed cross-backend reproducibility.
                        q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = m1
                        q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = m2

                        call s_infinite_ptg_relaxation_k(j, k, l, pS, rhoe, q_cons_vf, TS, ns_tmp)
                        ns_pc = ns_pc + ns_tmp
                    end if

                    ! Calculations AFTER equilibrium

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! entropy
                        sk(i) = cvs(i)*log((TS**gs_min(i))/((pS + ps_inf(i))**(gs_min(i) - 1.0_wp))) + qvps(i)

                        ! enthalpy
                        hk(i) = gs_min(i)*cvs(i)*TS + qvs(i)

                        ! Gibbs-free energy
                        gk(i) = hk(i) - TS*sk(i)

                        ! densities
                        rhok(i) = (pS + ps_inf(i))/((gs_min(i) - 1)*cvs(i)*TS)

                        ! internal energy
                        ek(i) = (pS + gs_min(i)*ps_inf(i))/(pS + ps_inf(i))*cvs(i)*TS + qvs(i)
                    end do

                    ! calculating volume fractions, internal energies, and total entropy
                    rhos = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! volume fractions
                        q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rhok(i)

                        ! alpha*rho*e
                        if (model_eqns == model_eqns_6eq) then
                            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                      & l)*ek(i)
                        end if

                        ! Total entropy
                        rhos = rhos + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*sk(i)
                    end do

#ifdef MFC_SIMULATION
                    ! Accumulate Newton iteration count for the load-weight diagnostic (matches the
                    ! allocation condition; no-op when neither writer is enabled).
                    if (count_pc_iters) pc_iter_count(j, k, l) = real(ns_pc, stp)
#endif
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_infinite_relaxation_k

    !> Apply pT-equilibrium relaxation for N fluids
    !! @param MFL flag: 0=gas, 1=liquid, 2=mixture
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, q_cons_vf, rhoe, TS, ns_out)

        $:GPU_ROUTINE(function_name='s_infinite_pt_relaxation_k', parallelism='[seq]', cray_noinline=True)

        ! initializing variables
        integer, intent(in)                                 :: j, k, l, MFL
        real(wp), intent(out)                               :: pS
        real(wp), dimension(1:), intent(out)                :: p_infpT
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), intent(in)                                :: rhoe
        real(wp), intent(out)                               :: TS
        integer, intent(out)                                :: ns_out                    !< Newton iteration count for this call
        real(wp)                                            :: gp, gpp, hp, pO, mCP, mQ  !< variables for the Newton Solver
        real(wp)                                            :: p_infpT_sum
        integer                                             :: i, ns                     !< generic loop iterators
        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0_wp; mQ = 0.0_wp; p_infpT_sum = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            p_infpT(i) = ps_inf(i)
            p_infpT_sum = p_infpT_sum + abs(p_infpT(i))
        end do
        ! Performing tests before initializing the pT-equilibrium
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
        end do

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            if (num_fluids < 3) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = num_fluids + 1, 3
                    p_infpT(i) = p_infpT_sum
                end do
            end if
        #:endif

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infpT)) < 0.0_wp) then
            if ((MFL == 0) .or. (MFL == 1)) then
                ! Assigning zero values for mass depletion cases pressure
                pS = 0.0_wp

                ! temperature
                TS = 0.0_wp

                ns_out = 0
                return
            end if
        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to iterate over
        ! the Newton's solver
        pO = 0.0_wp

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf) and infinity, a solution
        ! should be able to be found.
        pS = 1.0e4_wp

        ! Newton Solver for the pT-equilibrium
        ns = 0
        ! change this relative error metric. 1.e4_wp is just arbitrary
        ! Relative criterion written in multiply form to avoid dividing by pO (pO = 0 on the first pass).
        do while ((abs(pS - pO) > palpha_eps) .and. (abs(pS - pO) > (palpha_eps/1.e4_wp)*abs(pO)) .or. (ns == 0))
            ! increasing counter
            ns = ns + 1
            ! guard against non-convergence: accept the last iterate rather than looping forever
            if (ns >= max_iter) exit

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0_wp; gp = 0.0_wp; hp = 0.0_wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                gp = gp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                           & l)*cvs(i)*(rhoe + pS - mQ)/(mCP*(pS + p_infpT(i)))

                gpp = gpp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                             & l)*cvs(i)*(p_infpT(i) - rhoe + mQ)/(mCP*(pS + p_infpT(i))**2)
            end do

            hp = 1.0_wp/(rhoe + pS - mQ) + 1.0_wp/(pS + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0_wp - gp)/gpp)/(1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp))/(2.0_wp*gpp)*hp)
        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP
        ns_out = ns

    end subroutine s_infinite_pt_relaxation_k

    !> Evaluate the pTg-equilibrium residual R2D and temperature TS at a trial state (ml, pS) WITHOUT mutating q_cons_vf, so the
    !! Newton driver can line-search. The total reacting mass mT is conserved, so the reacting masses are (ml, mT - ml) and only the
    !! inert fluids are read from q_cons_vf. Also returns the mixture sums the Jacobian and the final temperature need.
    subroutine s_compute_ptg_residual(ml, mT, pS, j, k, l, q_cons_vf, rhoe, R2D, TS, mCP, mQ, mCVGP, mCVGP2, mCPD)

        $:GPU_ROUTINE(function_name='s_compute_ptg_residual', parallelism='[seq]')

        real(wp), intent(in)                                :: ml, mT, pS, rhoe
        integer, intent(in)                                 :: j, k, l
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(2), intent(out)                 :: R2D
        real(wp), intent(out)                               :: TS, mCP, mQ, mCVGP, mCVGP2, mCPD
        real(wp)                                            :: mQD
        integer                                             :: i

        ! reacting fluids contribute via (ml, mT - ml); inert fluids are summed from q_cons_vf
        mCP = ml*cvs(lp)*gs_min(lp) + (mT - ml)*cvs(vp)*gs_min(vp)
        mQ = ml*qvs(lp) + (mT - ml)*qvs(vp)
        mCVGP = 0.0_wp; mCVGP2 = 0.0_wp; mCPD = 0.0_wp; mQD = 0.0_wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if ((i /= lp) .and. (i /= vp)) then
                mCP = mCP + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*cvs(i)*gs_min(i)
                mQ = mQ + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
                mCVGP = mCVGP + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))
                mCVGP2 = mCVGP2 + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)
                mQD = mQD + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
                mCPD = mCPD + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*cvs(i)*gs_min(i)
            end if
        end do

        TS = 1.0_wp/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) - cvs(vp) &
                     & *(gs_min(vp) - 1)/(pS + ps_inf(vp))) + mCVGP)

        ! (i) Gibbs free-energy equality
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*(1 - log(TS)) - (qvps(lp) - qvps(vp)) + cvs(lp)*(gs_min(lp) - 1) &
            & *log(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))) + qvs(lp) - qvs(vp)

        ! (ii) constant-energy condition
        R2D(2) = rhoe + pS + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            & - mT*gs_min(vp)*cvs(vp) - mCPD)/(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS &
            & + ps_inf(vp))) + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)

    end subroutine s_compute_ptg_residual

    !> Apply pTg-equilibrium relaxation: a damped (backtracking line search) Newton solve for the reacting liquid mass ml and
    !! pressure pS enforcing Gibbs equality and energy conservation, converging on the residual norm (absolute ptgalpha_eps, or the
    !! rhoe-relative branch). Every step is projected onto the physical bounds 0 <= ml <= mT, pS > pmin. This converges in a handful
    !! of iterations with a bounded, uniform count (no GPU warp divergence), unlike the former fixed 1e-3 underrelaxation that
    !! stalled far from the root.
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, rhoe, q_cons_vf, TS, ns_out)

        $:GPU_ROUTINE(function_name='s_infinite_ptg_relaxation_k', parallelism='[seq]', cray_noinline=True)

        integer, intent(in)                                    :: j, k, l
        real(wp), intent(inout)                                :: pS
        real(wp), intent(in)                                   :: rhoe
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout)                                :: TS
        integer, intent(out)                                   :: ns_out  !< Newton iteration count for this call
        real(wp), dimension(2, 2)                              :: Jac, InvJac
        real(wp), dimension(2)                                 :: R2D, R2D_try, DeltamP
        real(wp)                                               :: mCP, mCPD, mCVGP, mCVGP2, mQ
        real(wp)                                               :: ml, ml_try, mT, pS_try, pmin, lambda, resnorm, resnorm_try
        real(wp)                                               :: dFdT, dTdm, dTdp, detJ
        integer                                                :: ns, ls

        ! total reacting mass is conserved; the liquid mass ml is the primary unknown, vapor mass = mT - ml
        mT = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) + q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)
        ml = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l)

        ! recover a physical pressure guess when the incoming pS is non-physical
        if (((pS < 0.0_wp) .and. ((q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) + q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, &
            & k, &
            & l)) > ((rhoe - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp)))) .or. ((pS >= 0.0_wp) .and. (pS < 1.0e-1_wp))) then
            pS = 1.0e4_wp
        end if

        ! pressure floor (stiffened gas requires pS + ps_inf > 0 for both phases)
        pmin = -min(ps_inf(lp), ps_inf(vp)) + 1.0_wp

        call s_compute_ptg_residual(ml, mT, pS, j, k, l, q_cons_vf, rhoe, R2D, TS, mCP, mQ, mCVGP, mCVGP2, mCPD)
        resnorm = sqrt(R2D(1)**2 + R2D(2)**2)

        do ns = 1, max_iter
            ! converged on the absolute residual, or on the rhoe-relative residual (multiply form, rhoe > 0)
            if ((resnorm <= ptgalpha_eps) .or. (resnorm <= (ptgalpha_eps/1.e6_wp)*rhoe)) exit

            ! 2x2 Jacobian of (Gibbs equality, energy) with respect to (ml, pS) at the current state
            dFdT = -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TS) - (qvps(lp) - qvps(vp)) + cvs(lp)*(gs_min(lp) - 1)*log(pS &
                     & + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))
            dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2
            dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 - cvs(vp) &
                    & *(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) + mCVGP2)*TS**2

            Jac(1, 1) = dFdT*dTdm
            Jac(1, 2) = dFdT*dTdp + TS*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))
            Jac(2, &
                & 1) = qvs(vp) - qvs(lp) + (cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp))/(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS &
                & + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                & + mCVGP) - (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) - mT*cvs(vp)*gs_min(vp) - mCPD)*(cvs(lp)*(gs_min(lp) &
                & - 1)/(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))/((ml*(cvs(lp)*(gs_min(lp) - 1)/(pS &
                & + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                & + mCVGP)**2)
            Jac(2, &
                & 2) = 1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) - mT*cvs(vp)*gs_min(vp) - mCPD)*(ml*(cvs(lp)*(gs_min(lp) &
                & - 1)/(pS + ps_inf(lp))**2 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) + mT*cvs(vp)*(gs_min(vp) - 1)/(pS &
                & + ps_inf(vp))**2 + mCVGP2)/(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) - cvs(vp)*(gs_min(vp) - 1)/(pS &
                & + ps_inf(vp))) + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2

            detJ = Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1)
            ! singular Jacobian: no usable Newton direction, accept the current (best) state
            if (detJ == 0.0_wp) exit

            InvJac(1, 1) = Jac(2, 2)/detJ
            InvJac(1, 2) = -Jac(1, 2)/detJ
            InvJac(2, 1) = -Jac(2, 1)/detJ
            InvJac(2, 2) = Jac(1, 1)/detJ

            DeltamP(1) = -(InvJac(1, 1)*R2D(1) + InvJac(1, 2)*R2D(2))
            DeltamP(2) = -(InvJac(2, 1)*R2D(1) + InvJac(2, 2)*R2D(2))

            ! backtracking line search: halve the step until the residual decreases, keeping the state
            ! physical (0 <= ml <= mT, pS above the stiffened-gas floor)
            lambda = 1.0_wp
            do ls = 1, ptg_ls_max
                ml_try = min(max(ml + lambda*DeltamP(1), 0.0_wp), mT)
                pS_try = max(pS + lambda*DeltamP(2), pmin)
                call s_compute_ptg_residual(ml_try, mT, pS_try, j, k, l, q_cons_vf, rhoe, R2D_try, TS, mCP, mQ, mCVGP, mCVGP2, mCPD)
                resnorm_try = sqrt(R2D_try(1)**2 + R2D_try(2)**2)
                if ((resnorm_try < resnorm) .or. (ls == ptg_ls_max)) exit
                lambda = 0.5_wp*lambda
            end do

            ! accept the trial state (TS, mCP, mQ, mCVGP, mCVGP2, mCPD already set to it by the last call)
            ml = ml_try; pS = pS_try; R2D = R2D_try; resnorm = resnorm_try
        end do

        ! commit the reacting masses (mT conserved) and set the common temperature
        q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = ml
        q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = mT - ml

        TS = (rhoe + pS - mQ)/mCP
        ns_out = ns

    end subroutine s_infinite_ptg_relaxation_k

    !> Correct the partial densities of the reacting fluids in case one of them is negative but their sum is positive. Inert phases
    !! are not corrected at this moment
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

        $:GPU_ROUTINE(function_name='s_correct_partial_densities', parallelism='[seq]', cray_noinline=True)

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), intent(out)                                  :: MCT
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout)                                :: rM
        integer, intent(in)                                    :: j, k, l
        !> @}
        if (rM < 0.0_wp) then
            if ((q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, &
                & l) >= -1.0_wp*mixM) .and. (q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) >= -1.0_wp*mixM)) then
                q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = 0.0_wp

                q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = 0.0_wp

                rM = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) + q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)
            end if
        end if

        ! TODO: Consider partitioning partial densities instead of absolute-value correction
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) < 0.0_wp) then
            q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM
        else if (q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) < 0.0_wp) then
            q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM

            q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = MCT*rM
        end if

    end subroutine s_correct_partial_densities

    !> Finalize the phase change module
    impure subroutine s_finalize_relaxation_solver_module

        if (allocated(pc_iter_count)) then
            @:DEALLOCATE(pc_iter_count)
        end if

    end subroutine s_finalize_relaxation_solver_module
#endif
end module m_phase_change
