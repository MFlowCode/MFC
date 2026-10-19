!>
!! @file
!! @brief Contains module @ref m_bubbles_ee "m_bubbles_EE"

#:include 'macros.fpp'
#:set MFC_BUBBLE_VECTOR_LENGTH = MFC_BUBBLE_VECTOR_LENGTH if defined('MFC_BUBBLE_VECTOR_LENGTH') else 128
#:include 'case.fpp'

!> @brief Computes ensemble-averaged (Euler--Euler) bubble source terms for radius, velocity, pressure, and mass transfer
module m_bubbles_EE

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_eos
    use m_bubbles
#ifdef MFC_FUSED_BUBBLE_SOURCE
    use m_bubbles_EE_openacc, only: s_compute_bubble_EE_source_rp_one
    use m_helper_basic, only: f_is_default
    use m_constants, only: bubble_model_keller_miksis, bubble_model_rayleigh_plesset, pi, small_alf
#endif

    implicit none

    real(wp), allocatable, dimension(:,:,:)   :: bub_adv_src
    real(wp), allocatable, dimension(:,:,:,:) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    $:GPU_DECLARE(create='[bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src]')

    type(scalar_field) :: divu  !< matrix for div(u)
    $:GPU_DECLARE(create='[divu]')

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    $:GPU_DECLARE(create='[rs, vs, ms, ps]')

contains

    logical function f_use_fused_bubble_source()

#ifdef MFC_FUSED_BUBBLE_SOURCE
        f_use_fused_bubble_source = bubbles_euler .and. .not. qbmm .and. .not. adap_dt .and. .not. adv_n &
            & .and. .not. bubbles_lagrange .and. polytropic .and. (bubble_model == bubble_model_keller_miksis &
            & .or. bubble_model == bubble_model_rayleigh_plesset) .and. model_eqns == model_eqns_5eq .and. num_fluids == 1 &
            & .and. .not. mpp_lim .and. nb >= 1
#else
        f_use_fused_bubble_source = .false.
#endif

    end function f_use_fused_bubble_source

    !> Initialize the Euler-Euler bubble module
    impure subroutine s_initialize_bubbles_EE_module

        integer :: l

        @:ALLOCATE(rs(1:nb))
        @:ALLOCATE(vs(1:nb))
        @:ALLOCATE(ps(1:nb))
        @:ALLOCATE(ms(1:nb))

        do l = 1, nb
            rs(l) = qbmm_idx%rs(l)
            vs(l) = qbmm_idx%vs(l)
            if (.not. polytropic) then
                ps(l) = qbmm_idx%ps(l)
                ms(l) = qbmm_idx%ms(l)
            else
                ps(l) = rs(l)
                ms(l) = rs(l)
            end if
        end do

        $:GPU_UPDATE(device='[rs, vs]')
        $:GPU_UPDATE(device='[ps, ms]')

        @:ALLOCATE(divu%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(divu)

#ifdef MFC_FUSED_RHS
        @:ALLOCATE(bub_adv_src(state_x_storage_beg:state_x_storage_end, 0:n, 0:p))
        @:ALLOCATE(bub_r_src(state_x_storage_beg:state_x_storage_end, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_v_src(state_x_storage_beg:state_x_storage_end, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_p_src(state_x_storage_beg:state_x_storage_end, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_m_src(state_x_storage_beg:state_x_storage_end, 0:n, 0:p, 1:nb))
#else
        @:ALLOCATE(bub_adv_src(0:m, 0:n, 0:p))
        @:ALLOCATE(bub_r_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_v_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_p_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_m_src(0:m, 0:n, 0:p, 1:nb))
#endif

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

    end subroutine s_initialize_bubbles_EE_module

    !> Compute the bubble volume fraction alpha from the bubble number density
    subroutine s_comp_alpha_from_n(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp)                                               :: nR3bar
        integer                                                :: i, j, k, l

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, nR3bar]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    nR3bar = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        nR3bar = nR3bar + weight(i)*(q_cons_vf(rs(i))%sf(j, k, l))**3._wp
                    end do
                    q_cons_vf(eqn_idx%alf)%sf(j, k, l) = (4._wp*pi*nR3bar)/(3._wp*q_cons_vf(eqn_idx%n)%sf(j, k, l)**2._wp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_comp_alpha_from_n

    !> Compute the right-hand side for Euler-Euler bubble transport
    subroutine s_compute_bubbles_EE_rhs(idir, q_prim_vf, divu_in)

        integer, intent(in)                                 :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(inout)                   :: divu_in  !< matrix for div(u)
        integer                                             :: j, k, l

        if (idir == 1) then
            if (.not. qbmm) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu_in%sf(j, k, l) = 0._wp
                            divu_in%sf(j, k, l) = 5.e-1_wp/dx(j)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j + 1, k, &
                                       & l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j - 1, k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (idir == 2) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dy(k)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, &
                                   & k + 1, l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k - 1, l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (idir == 3) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dz(l)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, &
                                   & l + 1) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubbles_EE_rhs

    !> Compute the Euler-Euler bubble source terms
    impure subroutine s_compute_bubble_EE_source(q_cons_vf, q_prim_vf, rhs_vf, divu_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(in) :: divu_in  !< matrix for div(u)
        real(wp) :: rddot
        real(wp) :: pb_local, mv_local, vflux, pbdot
        real(wp) :: n_tait, B_tait, qv_bub
        real(wp) :: chi_vw_l, k_mw_l, rho_mw_l     !< Per-thread bubble-wall scratch (avoid module-scalar race)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: Rtmp, Vtmp
            real(wp), dimension(3) :: myalpha, myalpha_rho
        #:else
            real(wp), dimension(nb)         :: Rtmp, Vtmp
            real(wp), dimension(num_fluids) :: myalpha, myalpha_rho
        #:endif
        real(wp) :: myR, myV, alf, myP, myRho, R2Vav, R3
        real(wp) :: nbub  !< Bubble number density
#ifdef MFC_FUSED_BUBBLE_SOURCE
        real(wp) :: c_liquid
        real(wp) :: pR0, pR, pV, pP, pRho, pC
        real(wp) :: pCa, pWeb, pRe, pGam, pEu, pRatio, pCpbw, pCstar, pMach, pNumerator
        logical  :: keller_miksis, web_active, re_active
#endif
        integer  :: i, j, k, l, q, ii               !< Loop variables
        integer  :: adap_dt_stop_sum, adap_dt_stop  !< Fail-safe exit if max iteration count reached
        integer  :: dmBub_id                        !< Dummy variables for unified subgrid bubble subroutines
        real(wp) :: dmMass_v, dmMass_n, dmBeta_c, dmBeta_t, dmCson
#ifdef MFC_FUSED_BUBBLE_SOURCE
        #:if MFC_CASE_OPTIMIZATION
            if (f_use_fused_bubble_source() .and. bubble_model == bubble_model_rayleigh_plesset .and. nb == 1) then
                call s_compute_bubble_EE_source_rp_one(m, n, p, lbound(q_prim_vf(rs(1))%sf), shape(q_prim_vf(rs(1))%sf), &
                                                       & q_cons_vf(1)%sf, q_cons_vf(vs(1))%sf, q_prim_vf(rs(1))%sf, &
                                                       & q_prim_vf(vs(1))%sf, q_prim_vf(eqn_idx%alf)%sf, q_prim_vf(eqn_idx%E)%sf, &
                                                       & rhs_vf(rs(1))%sf, rhs_vf(vs(1))%sf, rhs_vf(eqn_idx%alf)%sf, weight(1), &
                                                       & R0(1), Ca, Web, Re_inv, gam, Eu, pi, small_alf, .not. f_is_default(Web), &
                                                       & .not. f_is_default(Re_inv))
                return
            end if
        #:endif

        if (f_use_fused_bubble_source()) then
            keller_miksis = bubble_model == bubble_model_keller_miksis
            web_active = .not. f_is_default(Web)
            re_active = .not. f_is_default(Re_inv)
            B_tait = f_isentrope_pressure(pi_infs(1)/pi_fac, gammas(1))
            n_tait = f_isentrope_exponent(gammas(1))
            pCa = real(Ca, wp)
            pWeb = real(Web, wp)
            pRe = real(Re_inv, wp)
            pGam = real(gam, wp)
            pEu = real(Eu, wp)
            $:GPU_PARALLEL_LOOP(private='[j, k, l, myR, myV, alf, myP, myRho, R2Vav, R3, nbub, rddot, c_liquid, pR0, pR, pV, pP, &
                                & pRho, pC, pRatio, pCpbw, pCstar, pMach, pNumerator]', firstprivate='[keller_miksis, web_active, &
                                & re_active, n_tait, B_tait, pCa, pWeb, pRe, pGam, pEu]', collapse=3, extraAccArgs='vector_length(%d)' % MFC_BUBBLE_VECTOR_LENGTH)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                        myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                        myRho = q_cons_vf(1)%sf(j, k, l)
                        R3 = 0._wp
                        R2Vav = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            myR = q_prim_vf(rs(q))%sf(j, k, l)
                            myV = q_prim_vf(vs(q))%sf(j, k, l)
                            R3 = R3 + weight(q)*myR**3._wp
                            R2Vav = R2Vav + weight(q)*myR**2._wp*myV
                        end do

                        nbub = (3._wp/(4._wp*pi))*alf/R3
                        if (alf >= small_alf) then
                            R2Vav = 4._wp*pi*nbub*R2Vav
                            rhs_vf(eqn_idx%alf)%sf(j, k, l) = rhs_vf(eqn_idx%alf)%sf(j, k, l) + R2Vav

                            pP = real(myP, wp)
                            pRho = real(myRho, wp)
                            if (keller_miksis) then
                                c_liquid = sqrt(n_tait*(myP + B_tait)/(myRho*(1._wp - alf)))
                                pC = real(c_liquid, wp)
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, nb
                                myR = q_prim_vf(rs(q))%sf(j, k, l)
                                myV = q_prim_vf(vs(q))%sf(j, k, l)
                                pR0 = real(R0(q), wp)
                                pR = real(myR, wp)
                                pV = real(myV, wp)
                                pRatio = (pR0/pR)**(3._wp*pGam)
                                pCpbw = pCa*pRatio - pCa + pEu
                                if (web_active) pCpbw = pCpbw + (2._wp/(pWeb*pR0))*pRatio
                                if (web_active) pCpbw = pCpbw - 2._wp/(pR*pWeb)
                                if (re_active) pCpbw = pCpbw - 4._wp*pRe*pV/pR

                                if (keller_miksis) then
                                    pCstar = -3._wp*pGam*pCa*pRatio*pV/pR
                                    if (web_active) pCstar = pCstar - 3._wp*pGam*(2._wp/(pWeb*pR0))*pRatio*pV/pR
                                    if (web_active) pCstar = pCstar + (2._wp/pWeb)*pV/(pR**2._wp)
                                    if (re_active) pCstar = pCstar + 4._wp*pRe*(pV/pR)**2._wp

                                    pMach = pV/pC
                                    pNumerator = 1.5_wp*pV**2._wp*(pMach/3._wp - 1._wp) + (1._wp + pMach)*(pCpbw - pP)/pRho &
                                                                   & + pCstar*pR/(pRho*pC)
                                    if (re_active) then
                                        rddot = real(pNumerator/(pR*(1._wp - pMach) + 4._wp*pRe/(pRho*pC)), wp)
                                    else
                                        rddot = real(pNumerator/(pR*(1._wp - pMach)), wp)
                                    end if
                                else
                                    rddot = real((-1.5_wp*pV**2._wp + (pCpbw - pP)/pRho)/pR, wp)
                                end if
                                rhs_vf(vs(q))%sf(j, k, l) = rhs_vf(vs(q))%sf(j, k, l) + nbub*rddot
                                rhs_vf(rs(q))%sf(j, k, l) = rhs_vf(rs(q))%sf(j, k, l) + q_cons_vf(vs(q))%sf(j, k, l)
                            end do
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if
#endif

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    bub_adv_src(j, k, l) = 0._wp

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        bub_r_src(j, k, l, q) = 0._wp
                        bub_v_src(j, k, l, q) = 0._wp
                        bub_p_src(j, k, l, q) = 0._wp
                        bub_m_src(j, k, l, q) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        #:if USING_AMD
            ! LLVM Flang miscompiles the adap_dt-specialized kernel pair below, so AMD builds use one runtime-branched kernel.
            adap_dt_stop_sum = 0
            $:GPU_PARALLEL_LOOP(private='[j, k, l, Rtmp, Vtmp, myalpha_rho, myalpha, myR, myV, alf, myP, myRho, R2Vav, R3, nbub, &
                                & pb_local, mv_local, vflux, pbdot, rddot, n_tait, B_tait, qv_bub, adap_dt_stop, chi_vw_l, &
                                & k_mw_l, rho_mw_l]', collapse=3, copy='[adap_dt_stop_sum]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            Rtmp(q) = q_prim_vf(rs(q))%sf(j, k, l)
                            Vtmp(q) = q_prim_vf(vs(q))%sf(j, k, l)
                        end do

                        if (adv_n) then
                            nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                        else
                            R3 = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, nb
                                R3 = R3 + weight(q)*Rtmp(q)**3._wp
                            end do

                            nbub = (3._wp/(4._wp*pi))*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/R3
                        end if

                        if (.not. adap_dt) then
                            R2Vav = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, nb
                                R2Vav = R2Vav + weight(q)*Rtmp(q)**2._wp*Vtmp(q)
                            end do

                            bub_adv_src(j, k, l) = 4._wp*pi*nbub*R2Vav
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            $:GPU_LOOP(parallelism='[seq]')
                            do ii = 1, num_fluids
                                myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                                myalpha(ii) = q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)
                            end do

                            call s_compute_mixture_coefficients(myalpha_rho, myalpha, myRho, n_tait, B_tait, qv_bub)
                            B_tait = B_tait/pi_fac

                            ! n_tait still holds the mixture gamma here, which the reference pressure needs
                            B_tait = f_isentrope_pressure(B_tait, n_tait)
                            n_tait = f_isentrope_exponent(n_tait)

                            myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                            alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                            myR = q_prim_vf(rs(q))%sf(j, k, l)
                            myV = q_prim_vf(vs(q))%sf(j, k, l)

                            if (alf < small_alf) then
                                bub_adv_src(j, k, l) = 0._wp
                                bub_r_src(j, k, l, q) = 0._wp
                                bub_v_src(j, k, l, q) = 0._wp
                                if (.not. polytropic) then
                                    bub_p_src(j, k, l, q) = 0._wp
                                    bub_m_src(j, k, l, q) = 0._wp
                                end if
                            else
                                if (.not. polytropic) then
                                    pb_local = q_prim_vf(ps(q))%sf(j, k, l)
                                    mv_local = q_prim_vf(ms(q))%sf(j, k, l)
                                    call s_bwproperty(pb_local, q, chi_vw_l, k_mw_l, rho_mw_l)
                                    call s_vflux(myR, myV, pb_local, mv_local, q, vflux, fchi_vw=chi_vw_l, frho_mw=rho_mw_l)
                                    pbdot = f_bpres_dot(vflux, myR, myV, pb_local, mv_local, q, fk_mw=k_mw_l)
                                    bub_p_src(j, k, l, q) = nbub*pbdot
                                    bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                                else
                                    pb_local = 0._wp; mv_local = 0._wp; vflux = 0._wp; pbdot = 0._wp
                                end if

                                adap_dt_stop = 0

                                ! Adaptive time stepping
                                if (adap_dt) then
                                    adap_dt_stop = f_advance_step(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, &
                                                                  & B_tait, bub_adv_src(j, k, l), real(divu_in%sf(j, k, l), &
                                                                  & kind=wp), dmBub_id, dmMass_v, dmMass_n, dmBeta_c, dmBeta_t, &
                                                                  & dmCson)

                                    q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                                    q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV
                                else
                                    rddot = f_rddot(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, &
                                                    & bub_adv_src(j, k, l), real(divu_in%sf(j, k, l), kind=wp), dmCson)
                                    bub_v_src(j, k, l, q) = nbub*rddot
                                    bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                                end if

                                $:GPU_ATOMIC(atomic='update')
                                adap_dt_stop_sum = adap_dt_stop_sum + adap_dt_stop
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        #:else
            adap_dt_stop_sum = 0
            #:for ADAPTIVE in [True, False]
                #:if ADAPTIVE
                    if (adap_dt) then
                    #:else
                    else
                    #:endif
                    $:GPU_PARALLEL_LOOP(private='[j, k, l, Rtmp, Vtmp, myalpha_rho, myalpha, myR, myV, alf, myP, myRho, R2Vav, &
                                        & R3, nbub, pb_local, mv_local, vflux, pbdot, rddot, n_tait, B_tait, qv_bub, &
                                        & adap_dt_stop, chi_vw_l, k_mw_l, rho_mw_l]', collapse=3, copy='[adap_dt_stop_sum]' if ADAPTIVE else None)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, nb
                                    Rtmp(q) = q_prim_vf(rs(q))%sf(j, k, l)
                                    Vtmp(q) = q_prim_vf(vs(q))%sf(j, k, l)
                                end do

                                if (adv_n) then
                                    nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                                else
                                    R3 = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, nb
                                        R3 = R3 + weight(q)*Rtmp(q)**3._wp
                                    end do

                                    nbub = (3._wp/(4._wp*pi))*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/R3
                                end if

                                if (${' .false.' if ADAPTIVE else '.true.'}$) then
                                    R2Vav = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, nb
                                        R2Vav = R2Vav + weight(q)*Rtmp(q)**2._wp*Vtmp(q)
                                    end do

                                    bub_adv_src(j, k, l) = 4._wp*pi*nbub*R2Vav
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, nb
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do ii = 1, num_fluids
                                        myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                                        myalpha(ii) = q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)
                                    end do

                                    call s_compute_mixture_coefficients(myalpha_rho, myalpha, myRho, n_tait, B_tait, qv_bub)
                                    B_tait = B_tait/pi_fac

                                    ! n_tait still holds the mixture gamma here, which the reference pressure needs
                                    B_tait = f_isentrope_pressure(B_tait, n_tait)
                                    n_tait = f_isentrope_exponent(n_tait)

                                    myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                                    alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                                    myR = q_prim_vf(rs(q))%sf(j, k, l)
                                    myV = q_prim_vf(vs(q))%sf(j, k, l)

                                    if (alf < small_alf) then
                                        bub_adv_src(j, k, l) = 0._wp
                                        bub_r_src(j, k, l, q) = 0._wp
                                        bub_v_src(j, k, l, q) = 0._wp
                                        if (.not. polytropic) then
                                            bub_p_src(j, k, l, q) = 0._wp
                                            bub_m_src(j, k, l, q) = 0._wp
                                        end if
                                    else
                                        if (.not. polytropic) then
                                            pb_local = q_prim_vf(ps(q))%sf(j, k, l)
                                            mv_local = q_prim_vf(ms(q))%sf(j, k, l)
                                            call s_bwproperty(pb_local, q, chi_vw_l, k_mw_l, rho_mw_l)
                                            call s_vflux(myR, myV, pb_local, mv_local, q, vflux, fchi_vw=chi_vw_l, frho_mw=rho_mw_l)
                                            pbdot = f_bpres_dot(vflux, myR, myV, pb_local, mv_local, q, fk_mw=k_mw_l)
                                            bub_p_src(j, k, l, q) = nbub*pbdot
                                            bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                                        else
                                            pb_local = 0._wp; mv_local = 0._wp; vflux = 0._wp; pbdot = 0._wp
                                        end if

                                        adap_dt_stop = 0

                                        ! Adaptive time stepping
                                        if (${' .true.' if ADAPTIVE else '.false.'}$) then
                                            adap_dt_stop = f_advance_step(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, &
                                                                          & n_tait, B_tait, bub_adv_src(j, k, l), &
                                                                          & real(divu_in%sf(j, k, l), kind=wp), dmBub_id, &
                                                                          & dmMass_v, dmMass_n, dmBeta_c, dmBeta_t, dmCson)

                                            q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                                            q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV
                                        else
                                            rddot = f_rddot(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, &
                                                            & bub_adv_src(j, k, l), real(divu_in%sf(j, k, l), kind=wp), dmCson)
                                            bub_v_src(j, k, l, q) = nbub*rddot
                                            bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                                        end if

                                        #:if ADAPTIVE
                                            $:GPU_ATOMIC(atomic='update')
                                            adap_dt_stop_sum = adap_dt_stop_sum + adap_dt_stop
                                        #:endif
                                    end if
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endfor
            end if
        #:endif

        if (adap_dt .and. adap_dt_stop_sum > 0) call s_mpi_abort("Adaptive time stepping failed to converge.")

        if (.not. adap_dt) then
            $:GPU_PARALLEL_LOOP(private='[i, k, l, q]', collapse=3)
            do l = 0, p
                do q = 0, n
                    do i = 0, m
                        rhs_vf(eqn_idx%alf)%sf(i, q, l) = rhs_vf(eqn_idx%alf)%sf(i, q, l) + bub_adv_src(i, q, l)
                        if (num_fluids > 1) rhs_vf(eqn_idx%adv%beg)%sf(i, q, l) = rhs_vf(eqn_idx%adv%beg)%sf(i, q, &
                            & l) - bub_adv_src(i, q, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do k = 1, nb
                            rhs_vf(rs(k))%sf(i, q, l) = rhs_vf(rs(k))%sf(i, q, l) + bub_r_src(i, q, l, k)
                            rhs_vf(vs(k))%sf(i, q, l) = rhs_vf(vs(k))%sf(i, q, l) + bub_v_src(i, q, l, k)
                            if (polytropic .neqv. .true.) then
                                rhs_vf(ps(k))%sf(i, q, l) = rhs_vf(ps(k))%sf(i, q, l) + bub_p_src(i, q, l, k)
                                rhs_vf(ms(k))%sf(i, q, l) = rhs_vf(ms(k))%sf(i, q, l) + bub_m_src(i, q, l, k)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubble_EE_source

    !> Evaluate the Euler-Euler pressure correction at cell centers and store it in ptil. This mirrors the modified mixture pressure
    !! that the HLLC solver applies to the reconstructed face states in m_riemann_solver_hllc.fpp, following Ando (2010), Eq. (2.5).
    !! It is a diagnostic: ptil is read only by the probe output path and never feeds back into the solution.
    impure subroutine s_compute_ptilde(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        real(wp)                                            :: R3bar, R3V2bar, PbwR3bar
        real(wp)                                            :: myRho, myP, alf, myR, myV, myPb
        integer                                             :: j, k, l, q, ii

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q, ii, R3bar, R3V2bar, PbwR3bar, myRho, myP, alf, myR, myV, myPb]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                    myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)

                    myRho = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do ii = 1, num_fluids
                        myRho = myRho + q_cons_vf(ii)%sf(j, k, l)
                    end do

                    if (qbmm) then
                        R3bar = mom_sp(1)%sf(j, k, l)
                        R3V2bar = mom_sp(3)%sf(j, k, l)
                        PbwR3bar = mom_sp(4)%sf(j, k, l)
                    else
                        R3bar = 0._wp
                        R3V2bar = 0._wp
                        PbwR3bar = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            myR = q_prim_vf(rs(q))%sf(j, k, l)
                            myV = q_prim_vf(vs(q))%sf(j, k, l)
                            if (polytropic) then
                                myPb = 0._wp
                            else
                                myPb = q_prim_vf(ps(q))%sf(j, k, l)
                            end if

                            R3bar = R3bar + weight(q)*(myR**3._wp)
                            R3V2bar = R3V2bar + weight(q)*(myR**3._wp)*(myV**2._wp)
                            PbwR3bar = PbwR3bar + weight(q)*f_cpbw_KM(R0(q), myR, myV, myPb)*(myR**3._wp)
                        end do
                    end if

                    ! Same degenerate-state guard the HLLC solver uses
                    if (alf < small_alf .or. R3bar < small_alf) then
                        ptil(j, k, l) = alf*myP
                    else
                        ptil(j, k, l) = alf*(myP - PbwR3bar/R3bar - myRho*R3V2bar/R3bar)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_ptilde

end module m_bubbles_EE
