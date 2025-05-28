!>
!! @file m_bubbles_EE.f90
!! @brief Contains module m_bubbles_EE

#:include 'macros.fpp'

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
module m_bubbles_EE

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< General bubble dynamics procedures

    implicit none

    real(wp), allocatable, dimension(:, :, :) :: bub_adv_src
    real(wp), allocatable, dimension(:, :, :, :) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    !$acc declare create(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

    type(scalar_field) :: divu !< matrix for div(u)
    !$acc declare create(divu)

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    !$acc declare create(rs, vs, ms, ps)

contains

    impure subroutine s_initialize_bubbles_EE_module

        integer :: l

        @:ALLOCATE(rs(1:nb))
        @:ALLOCATE(vs(1:nb))
        if (.not. polytropic) then
            @:ALLOCATE(ps(1:nb))
            @:ALLOCATE(ms(1:nb))
        end if

        do l = 1, nb
            rs(l) = bub_idx%rs(l)
            vs(l) = bub_idx%vs(l)
            if (.not. polytropic) then
                ps(l) = bub_idx%ps(l)
                ms(l) = bub_idx%ms(l)
            end if
        end do

        !$acc update device(rs, vs)
        if (.not. polytropic) then
            !$acc update device(ps, ms)
        end if

        @:ALLOCATE(divu%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(divu)

        @:ALLOCATE(bub_adv_src(0:m, 0:n, 0:p))
        @:ALLOCATE(bub_r_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_v_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_p_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_m_src(0:m, 0:n, 0:p, 1:nb))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

    end subroutine s_initialize_bubbles_EE_module

    ! Compute the bubble volume fraction alpha from the bubble number density n
        !! @param q_cons_vf is the conservative variable
    pure subroutine s_comp_alpha_from_n(q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: nR3bar
        integer(wp) :: i, j, k, l

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    nR3bar = 0._wp
                    !$acc loop seq
                    do i = 1, nb
                        nR3bar = nR3bar + weight(i)*(q_cons_vf(rs(i))%sf(j, k, l))**3._wp
                    end do
                    q_cons_vf(alf_idx)%sf(j, k, l) = (4._wp*pi*nR3bar)/(3._wp*q_cons_vf(n_idx)%sf(j, k, l)**2._wp)
                end do
            end do
        end do

    end subroutine s_comp_alpha_from_n

    pure subroutine s_compute_bubbles_EE_rhs(idir, q_prim_vf, divu)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(inout) :: divu !< matrix for div(u)

        integer :: j, k, l

        if (idir == 1) then

            if (.not. qbmm) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu%sf(j, k, l) = 0._wp
                            divu%sf(j, k, l) = &
                                5e-1_wp/dx(j)*(q_prim_vf(contxe + idir)%sf(j + 1, k, l) - &
                                               q_prim_vf(contxe + idir)%sf(j - 1, k, l))

                        end do
                    end do
                end do
            end if

        elseif (idir == 2) then

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu%sf(j, k, l) = divu%sf(j, k, l) + &
                                           5e-1_wp/dy(k)*(q_prim_vf(contxe + idir)%sf(j, k + 1, l) - &
                                                          q_prim_vf(contxe + idir)%sf(j, k - 1, l))

                    end do
                end do
            end do

        elseif (idir == 3) then

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu%sf(j, k, l) = divu%sf(j, k, l) + &
                                           5e-1_wp/dz(l)*(q_prim_vf(contxe + idir)%sf(j, k, l + 1) - &
                                                          q_prim_vf(contxe + idir)%sf(j, k, l - 1))

                    end do
                end do
            end do

        end if

    end subroutine s_compute_bubbles_EE_rhs

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the bubble modeling
        !!  @param q_prim_vf Primitive variables
        !!  @param q_cons_vf Conservative variables
    impure subroutine s_compute_bubble_EE_source(q_cons_vf, q_prim_vf, t_step, rhs_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(wp) :: rddot
        real(wp) :: pb, mv, vflux, pbdot
        real(wp) :: n_tait, B_tait
        real(wp), dimension(nb) :: Rtmp, Vtmp
        real(wp) :: myR, myV, alf, myP, myRho, R2Vav, R3
        real(wp), dimension(num_fluids) :: myalpha, myalpha_rho
        real(wp) :: nbub !< Bubble number density

        integer :: i, j, k, l, q, ii !< Loop variables

        integer :: adap_dt_stop_max, adap_dt_stop !< Fail-safe exit if max iteration count reached
        integer :: dmBub_id !< Dummy variables for unified subgrid bubble subroutines
        real(wp) :: dmMass_v, dmMass_n, dmBeta_c, dmBeta_t, dmCson

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    bub_adv_src(j, k, l) = 0._wp

                    !$acc loop seq
                    do q = 1, nb
                        bub_r_src(j, k, l, q) = 0._wp
                        bub_v_src(j, k, l, q) = 0._wp
                        bub_p_src(j, k, l, q) = 0._wp
                        bub_m_src(j, k, l, q) = 0._wp
                    end do
                end do
            end do
        end do

        adap_dt_stop_max = 0
        !$acc parallel loop collapse(3) gang vector default(present) private(Rtmp, Vtmp, myalpha_rho, myalpha)  &
        !$acc reduction(MAX:adap_dt_stop_max) copy(adap_dt_stop_max)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    if (adv_n) then
                        nbub = q_prim_vf(n_idx)%sf(j, k, l)
                    else
                        !$acc loop seq
                        do q = 1, nb
                            Rtmp(q) = q_prim_vf(rs(q))%sf(j, k, l)
                            Vtmp(q) = q_prim_vf(vs(q))%sf(j, k, l)
                        end do

                        R3 = 0._wp

                        !$acc loop seq
                        do q = 1, nb
                            R3 = R3 + weight(q)*Rtmp(q)**3._wp
                        end do

                        nbub = (3._wp/(4._wp*pi))*q_prim_vf(alf_idx)%sf(j, k, l)/R3
                    end if

                    if (.not. adap_dt) then
                        R2Vav = 0._wp

                        !$acc loop seq
                        do q = 1, nb
                            R2Vav = R2Vav + weight(q)*Rtmp(q)**2._wp*Vtmp(q)
                        end do

                        bub_adv_src(j, k, l) = 4._wp*pi*nbub*R2Vav
                    end if

                    !$acc loop seq
                    do q = 1, nb

                        !$acc loop seq
                        do ii = 1, num_fluids
                            myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                            myalpha(ii) = q_cons_vf(advxb + ii - 1)%sf(j, k, l)
                        end do

                        myRho = 0._wp
                        n_tait = 0._wp
                        B_tait = 0._wp

                        if (mpp_lim .and. (num_fluids > 2)) then
                            !$acc loop seq
                            do ii = 1, num_fluids
                                myRho = myRho + myalpha_rho(ii)
                                n_tait = n_tait + myalpha(ii)*gammas(ii)
                                B_tait = B_tait + myalpha(ii)*pi_infs(ii)
                            end do
                        else if (num_fluids > 2) then
                            !$acc loop seq
                            do ii = 1, num_fluids - 1
                                myRho = myRho + myalpha_rho(ii)
                                n_tait = n_tait + myalpha(ii)*gammas(ii)
                                B_tait = B_tait + myalpha(ii)*pi_infs(ii)
                            end do
                        else
                            myRho = myalpha_rho(1)
                            n_tait = gammas(1)
                            B_tait = pi_infs(1)/pi_fac
                        end if

                        n_tait = 1._wp/n_tait + 1._wp !make this the usual little 'gamma'
                        B_tait = B_tait*(n_tait - 1)/n_tait ! make this the usual pi_inf

                        myRho = q_prim_vf(1)%sf(j, k, l)
                        myP = q_prim_vf(E_idx)%sf(j, k, l)
                        alf = q_prim_vf(alf_idx)%sf(j, k, l)
                        myR = q_prim_vf(rs(q))%sf(j, k, l)
                        myV = q_prim_vf(vs(q))%sf(j, k, l)

                        if (.not. polytropic) then
                            pb = q_prim_vf(ps(q))%sf(j, k, l)
                            mv = q_prim_vf(ms(q))%sf(j, k, l)
                            call s_bwproperty(pb, q, chi_vw, k_mw, rho_mw)
                            call s_vflux(myR, myV, pb, mv, q, vflux)
                            pbdot = f_bpres_dot(vflux, myR, myV, pb, mv, q)

                            bub_p_src(j, k, l, q) = nbub*pbdot
                            bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                        else
                            pb = 0._wp; mv = 0._wp; vflux = 0._wp; pbdot = 0._wp
                        end if

                        ! Adaptive time stepping
                        adap_dt_stop = 0

                        if (adap_dt) then

                            call s_advance_step(myRho, myP, myR, myV, R0(q), &
                                                pb, pbdot, alf, n_tait, B_tait, &
                                                bub_adv_src(j, k, l), divu%sf(j, k, l), &
                                                dmBub_id, dmMass_v, dmMass_n, dmBeta_c, &
                                                dmBeta_t, dmCson, adap_dt_stop)

                            q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                            q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV

                        else
                            rddot = f_rddot(myRho, myP, myR, myV, R0(q), &
                                            pb, pbdot, alf, n_tait, B_tait, &
                                            bub_adv_src(j, k, l), divu%sf(j, k, l), &
                                            dmCson)
                            bub_v_src(j, k, l, q) = nbub*rddot
                            bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                        end if

                        adap_dt_stop_max = max(adap_dt_stop_max, adap_dt_stop)

                        if (alf < 1.e-11_wp) then
                            bub_adv_src(j, k, l) = 0._wp
                            bub_r_src(j, k, l, q) = 0._wp
                            bub_v_src(j, k, l, q) = 0._wp
                            if (.not. polytropic) then
                                bub_p_src(j, k, l, q) = 0._wp
                                bub_m_src(j, k, l, q) = 0._wp
                            end if
                        end if
                    end do
                end do
            end do
        end do

        if (adap_dt .and. adap_dt_stop_max > 0) call s_mpi_abort("Adaptive time stepping failed to converge.")

        if (.not. adap_dt) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do q = 0, n
                    do i = 0, m
                        rhs_vf(alf_idx)%sf(i, q, l) = rhs_vf(alf_idx)%sf(i, q, l) + bub_adv_src(i, q, l)
                        if (num_fluids > 1) rhs_vf(advxb)%sf(i, q, l) = &
                            rhs_vf(advxb)%sf(i, q, l) - bub_adv_src(i, q, l)
                        !$acc loop seq
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
        end if
    end subroutine s_compute_bubble_EE_source

end module m_bubbles_EE
