!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles

#:include 'macros.fpp'

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
module m_bubbles

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    implicit none

    real(kind(0._wp)) :: chi_vw  !< Bubble wall properties (Ando 2010)
    real(kind(0._wp)) :: k_mw    !< Bubble wall properties (Ando 2010)
    real(kind(0._wp)) :: rho_mw  !< Bubble wall properties (Ando 2010)
!$acc declare create(chi_vw, k_mw, rho_mw)

#ifdef CRAY_ACC_WAR
    !> @name Bubble dynamic source terms
    !> @{

    @:CRAY_DECLARE_GLOBAL(real(wp), dimension(:, :, :), bub_adv_src)
    !$acc declare link(bub_adv_src)

    @:CRAY_DECLARE_GLOBAL(real(wp), dimension(:, :, :, :), bub_r_src, bub_v_src, bub_p_src, bub_m_src)
    !$acc declare link(bub_r_src, bub_v_src, bub_p_src, bub_m_src)

    type(scalar_field) :: divu !< matrix for div(u)
    !$acc declare create(divu)

    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), rs, vs, ms, ps)
    !$acc declare link(rs, vs, ms, ps)
#else
    real(wp), allocatable, dimension(:, :, :) :: bub_adv_src
    real(wp), allocatable, dimension(:, :, :, :) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    !$acc declare create(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

    type(scalar_field) :: divu !< matrix for div(u)
    !$acc declare create(divu)

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    !$acc declare create(rs, vs, ms, ps)
#endif

contains

    subroutine s_initialize_bubbles_module

        integer :: i, j, k, l, q
        type(int_bounds_info) :: ix, iy, iz

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        @:ALLOCATE_GLOBAL(rs(1:nb))
        @:ALLOCATE_GLOBAL(vs(1:nb))
        if (.not. polytropic) then
            @:ALLOCATE_GLOBAL(ps(1:nb))
            @:ALLOCATE_GLOBAL(ms(1:nb))
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

        @:ALLOCATE(divu%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        @:ACC_SETUP_SFs(divu)

        @:ALLOCATE_GLOBAL(bub_adv_src(0:m, 0:n, 0:p))
        @:ALLOCATE_GLOBAL(bub_r_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE_GLOBAL(bub_v_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE_GLOBAL(bub_p_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE_GLOBAL(bub_m_src(0:m, 0:n, 0:p, 1:nb))

    end subroutine s_initialize_bubbles_module

    ! Compute the bubble volume fraction alpha from the bubble number density n
        !! @param q_cons_vf is the conservative variable
    subroutine s_comp_alpha_from_n(q_cons_vf)
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

    subroutine s_compute_bubbles_rhs(idir, q_prim_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: i, j, k, l, q

        if (idir == 1) then

            if (.not. qbmm) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu%sf(j, k, l) = 0._wp
                            divu%sf(j, k, l) = &
                                5d-1/dx(j)*(q_prim_vf(contxe + idir)%sf(j + 1, k, l) - &
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
                                           5d-1/dy(k)*(q_prim_vf(contxe + idir)%sf(j, k + 1, l) - &
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
                                           5d-1/dz(l)*(q_prim_vf(contxe + idir)%sf(j, k, l + 1) - &
                                                       q_prim_vf(contxe + idir)%sf(j, k, l - 1))

                    end do
                end do
            end do

        end if

    end subroutine s_compute_bubbles_rhs

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the bubble modeling
        !!  @param q_prim_vf Primitive variables
        !!  @param q_cons_vf Conservative variables
    subroutine s_compute_bubble_source(q_cons_vf, q_prim_vf, t_step, rhs_vf)
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
        real(wp) :: start, finish

        real(wp) :: nbub !< Bubble number density

        real(wp), dimension(2) :: Re !< Reynolds number

        integer :: i, j, k, l, q, ii !< Loop variables
        integer :: ndirs  !< Number of coordinate directions

        real(wp) :: err1, err2, err3, err4, err5 !< Error estimates for adaptive time stepping
        real(wp) :: t_new !< Updated time step size
        real(wp) :: h !< Time step size
        real(wp), dimension(4) :: myR_tmp1, myV_tmp1, myR_tmp2, myV_tmp2 !< Bubble radius, radial velocity, and radial acceleration for the inner loop

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

        !$acc parallel loop collapse(3) gang vector default(present) private(Rtmp, Vtmp, myalpha_rho, myalpha, myR_tmp1, myV_tmp1, myR_tmp2, myV_tmp2)
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
                            call s_bwproperty(pb, q)
                            vflux = f_vflux(myR, myV, mv, q)
                            pbdot = f_bpres_dot(vflux, myR, myV, pb, mv, q)

                            bub_p_src(j, k, l, q) = nbub*pbdot
                            bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                        else
                            pb = 0._wp; mv = 0._wp; vflux = 0._wp; pbdot = 0._wp
                        end if

                        ! Adaptive time stepping
                        if (adap_dt) then
                            ! Determine the starting time step
                            call s_initialize_adap_dt(myRho, myP, myR, myV, R0(q), &
                                                      pb, pbdot, alf, n_tait, B_tait, &
                                                      bub_adv_src(j, k, l), divu%sf(j, k, l), h)

                            ! Advancing one step
                            t_new = 0._wp
                            do while (.true.)
                                if (t_new + h > 0.5_wp*dt) then
                                    h = 0.5_wp*dt - t_new
                                end if

                                ! Advancing one sub-step
                                do while (.true.)
                                    ! Advance one sub-step
                                    call s_advance_substep(myRho, myP, myR, myV, R0(q), &
                                                           pb, pbdot, alf, n_tait, B_tait, &
                                                           bub_adv_src(j, k, l), divu%sf(j, k, l), h, &
                                                           myR_tmp1, myV_tmp1, err1)

                                    ! Advance one sub-step by advancing two half steps
                                    call s_advance_substep(myRho, myP, myR, myV, R0(q), &
                                                           pb, pbdot, alf, n_tait, B_tait, &
                                                           bub_adv_src(j, k, l), divu%sf(j, k, l), 0.5_wp*h, &
                                                           myR_tmp2, myV_tmp2, err2)

                                    call s_advance_substep(myRho, myP, myR_tmp2(4), myV_tmp2(4), R0(q), &
                                                           pb, pbdot, alf, n_tait, B_tait, &
                                                           bub_adv_src(j, k, l), divu%sf(j, k, l), 0.5_wp*h, &
                                                           myR_tmp2, myV_tmp2, err3)

                                    err4 = abs((myR_tmp1(4) - myR_tmp2(4))/myR_tmp1(4))
                                    err5 = abs((myV_tmp1(4) - myV_tmp2(4))/myV_tmp1(4))
                                    if (abs(myV_tmp1(4)) < 1e-12) err5 = 0._wp

                                    ! Determine acceptance/rejection and update step size
                                    !   Rule 1: err1, err2, err3 < tol
                                    !   Rule 2: myR_tmp1(4) > 0._wp
                                    !   Rule 3: abs((myR_tmp1(4) - myR_tmp2(4))/myR) < tol
                                    !   Rule 4: abs((myV_tmp1(4) - myV_tmp2(4))/myV) < tol
                                    if ((err1 <= 1d-4) .and. (err2 <= 1d-4) .and. (err3 <= 1d-4) &
                                        .and. (err4 < 1d-4) .and. (err5 < 1d-4) &
                                        .and. myR_tmp1(4) > 0._wp) then

                                        ! Accepted. Finalize the sub-step
                                        t_new = t_new + h

                                        ! Update R and V
                                        myR = myR_tmp1(4)
                                        myV = myV_tmp1(4)

                                        ! Update step size for the next sub-step
                                        h = h*min(2._wp, max(0.5_wp, (1d-4/err1)**(1._wp/3._wp)))

                                        exit
                                    else
                                        ! Rejected. Update step size for the next try on sub-step
                                        if (err2 <= 1d-4) then
                                            h = 0.5_wp*h
                                        else
                                            h = 0.25_wp*h
                                        end if

                                    end if
                                end do

                                ! Exit the loop if the final time reached dt
                                if (t_new == 0.5_wp*dt) exit

                            end do

                            q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                            q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV

                        else
                            rddot = f_rddot(myRho, myP, myR, myV, R0(q), &
                                            pb, pbdot, alf, n_tait, B_tait, &
                                            bub_adv_src(j, k, l), divu%sf(j, k, l))
                            bub_v_src(j, k, l, q) = nbub*rddot
                            bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                        end if

                        if (alf < 1.d-11) then
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
    end subroutine s_compute_bubble_source

    !> Choose the initial time step size for the adaptive time stepping routine
        !!  (See Heirer, E. Hairer S.P.Nørsett G. Wanner, Solving Ordinary
        !!  Differential Equations I, Chapter II.4)
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param h Time step size
    subroutine s_initialize_adap_dt(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                    fntait, fBtait, f_bub_adv_src, f_divu, h)
        !$acc routine seq
        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(wp), intent(out) :: h

        real(wp) :: h0, h1, h_min !< Time step size
        real(wp) :: d0, d1, d2 !< norms
        real(wp), dimension(2) :: myR_tmp, myV_tmp, myA_tmp !< Bubble radius, radial velocity, and radial acceleration

        ! Determine the starting time step
        ! Evaluate f(x0,y0)
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Compute d0 = ||y0|| and d1 = ||f(x0,y0)||
        d0 = DSQRT((myR_tmp(1)**2._wp + myV_tmp(1)**2._wp)/2._wp)
        d1 = DSQRT((myV_tmp(1)**2._wp + myA_tmp(1)**2._wp)/2._wp)
        if (d0 < 1d-5 .or. d1 < 1d-5) then
            h0 = 1d-6
        else
            h0 = 1d-2*(d0/d1)
        end if

        ! Evaluate f(x0+h0,y0+h0*f(x0,y0))
        myR_tmp(2) = myR_tmp(1) + h0*myV_tmp(1)
        myV_tmp(2) = myV_tmp(1) + h0*myA_tmp(1)
        myA_tmp(2) = f_rddot(fRho, fP, myR_tmp(2), myV_tmp(2), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Compute d2 = ||f(x0+h0,y0+h0*f(x0,y0))-f(x0,y0)||/h0
        d2 = DSQRT(((myV_tmp(2) - myV_tmp(1))**2._wp + (myA_tmp(2) - myA_tmp(1))**2._wp)/2._wp)/h0

        ! Set h1 = (0.01/max(d1,d2))^{1/(p+1)}
        !      if max(d1,d2) < 1e-15, h1 = max(1e-6, h0*1e-3)
        if (max(d1, d2) < 1d-15) then
            h1 = max(1d-6, h0*1d-3)
        else
            h1 = (1d-2/max(d1, d2))**(1._wp/3._wp)
        end if

        ! Set h = min(100*h0,h1)
        h = min(100._wp*h0, h1)

    end subroutine s_initialize_adap_dt

    !>  Integrate bubble variables over the given time step size, h
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param h Time step size
        !!  @param myR_tmp Bubble radius at each stage
        !!  @param myV_tmp Bubble radial velocity at each stage
        !!  @param err Estimated error
    subroutine s_advance_substep(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                 fntait, fBtait, f_bub_adv_src, f_divu, h, &
                                 myR_tmp, myV_tmp, err)
        !$acc routine seq
        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu, h
        real(wp), dimension(4), intent(OUT) :: myR_tmp, myV_tmp
        real(wp), dimension(4) :: myA_tmp
        real(wp), intent(OUT) :: err
        real(wp) :: err_R, err_V

        ! Stage 0
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Stage 1
        myR_tmp(2) = myR_tmp(1) + h*myV_tmp(1)
        myV_tmp(2) = myV_tmp(1) + h*myA_tmp(1)
        myA_tmp(2) = f_rddot(fRho, fP, myR_tmp(2), myV_tmp(2), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Stage 2
        myR_tmp(3) = myR_tmp(1) + (h/4._wp)*(myV_tmp(1) + myV_tmp(2))
        myV_tmp(3) = myV_tmp(1) + (h/4._wp)*(myA_tmp(1) + myA_tmp(2))
        myA_tmp(3) = f_rddot(fRho, fP, myR_tmp(3), myV_tmp(3), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Stage 3
        myR_tmp(4) = myR_tmp(1) + (h/6._wp)*(myV_tmp(1) + myV_tmp(2) + 4._wp*myV_tmp(3))
        myV_tmp(4) = myV_tmp(1) + (h/6._wp)*(myA_tmp(1) + myA_tmp(2) + 4._wp*myA_tmp(3))
        myA_tmp(4) = f_rddot(fRho, fP, myR_tmp(4), myV_tmp(4), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Estimate error
        err_R = (-5._wp*h/24._wp)*(myV_tmp(2) + myV_tmp(3) - 2._wp*myV_tmp(4)) &
                /max(abs(myR_tmp(1)), abs(myR_tmp(4)))
        err_V = (-5._wp*h/24._wp)*(myA_tmp(2) + myA_tmp(3) - 2._wp*myA_tmp(4)) &
                /max(abs(myV_tmp(1)), abs(myV_tmp(4)))
        err = DSQRT((err_R**2._wp + err_V**2._wp)/2._wp)

    end subroutine s_advance_substep

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw(fR0, fR, fV, fpb)
        !$acc routine seq
        real(wp), intent(in) :: fR0, fR, fV, fpb

        real(wp) :: f_cpbw

        if (polytropic) then
            f_cpbw = (Ca + 2._wp/Web/fR0)*((fR0/fR)**(3._wp*gam)) - Ca - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        else
            f_cpbw = fpb - 1._wp - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        end if

    end function f_cpbw

    !>  Function that computes the bubble enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    function f_H(fCpbw, fCpinf, fntait, fBtait)
        !$acc routine seq
        real(wp), intent(in) :: fCpbw, fCpinf, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_H

        tmp1 = (fntait - 1._wp)/fntait
        tmp2 = (fCpbw/(1._wp + fBtait) + 1._wp)**tmp1
        tmp3 = (fCpinf/(1._wp + fBtait) + 1._wp)**tmp1

        f_H = (tmp2 - tmp3)*fntait*(1._wp + fBtait)/(fntait - 1._wp)

    end function f_H

    !> Function that computes the sound speed for the bubble
        !! @param fCpinf Driving bubble pressure
        !! @param fntait Tait EOS parameter
        !! @param fBtait Tait EOS parameter
        !! @param fH Bubble enthalpy
    function f_cgas(fCpinf, fntait, fBtait, fH)
        !$acc routine seq
        real(wp), intent(in) :: fCpinf, fntait, fBtait, fH

        real(wp) :: tmp
        real(wp) :: f_cgas

        ! get sound speed for Gilmore equations "C" -> c_gas
        tmp = (fCpinf/(1._wp + fBtait) + 1._wp)**((fntait - 1._wp)/fntait)
        tmp = fntait*(1._wp + fBtait)*tmp

        f_cgas = dsqrt(tmp + (fntait - 1._wp)*fH)

    end function f_cgas

    !>  Function that computes the time derivative of the driving pressure
        !!  @param fRho Local liquid density
        !!  @param fP Local pressure
        !!  @param falf Local void fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param advsrc Advection equation source term
        !!  @param divu Divergence of velocity
    function f_cpinfdot(fRho, fP, falf, fntait, fBtait, advsrc, divu)
        !$acc routine seq
        real(wp), intent(in) :: fRho, fP, falf, fntait, fBtait, advsrc, divu

        real(wp) :: c2_liquid
        real(wp) :: f_cpinfdot

        ! get sound speed squared for liquid (only needed for pbdot)
        ! c_l^2 = gam (p+B) / (rho*(1-alf))
        if (mpp_lim) then
            c2_liquid = fntait*(fP + fBtait)/fRho
        else
            c2_liquid = fntait*(fP + fBtait)/(fRho*(1._wp - falf))
        end if

        ! \dot{Cp_inf} = rho sound^2 (alf_src - divu)
        f_cpinfdot = fRho*c2_liquid*(advsrc - divu)

    end function f_cpinfdot

    !>  Function that computes the time derivative of the enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fCpinf_dot Time derivative of the driving pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpbdot Time derivative of the internal bubble pressure
    function f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
        !$acc routine seq
        real(wp), intent(in) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
        real(wp), intent(in) :: fR, fV, fR0, fpbdot

        real(wp) :: tmp1, tmp2
        real(wp) :: f_Hdot

        if (polytropic) then
            tmp1 = (fR0/fR)**(3._wp*gam)
            tmp1 = -3._wp*gam*(Ca + 2._wp/Web/fR0)*tmp1*fV/fR
        else
            tmp1 = fpbdot
        end if
        tmp2 = (2._wp/Web + 4._wp*Re_inv*fV)*fV/(fR**2._wp)

        f_Hdot = &
            (fCpbw/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*(tmp1 + tmp2) &
            - (fCpinf/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*fCpinf_dot

        ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
        !f_Hdot = ((fCpbw/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*(-3._wp)*gam * &
        !            ( (fR0/fR)**(3._wp*gam ))*(fV/fR)

        ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
        !f_Hdot = f_Hdot - ((fCpinf/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*fCpinf_dot

    end function f_Hdot

    !> Function that computes the bubble radial acceleration based on bubble models
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
    function f_rddot(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, fntait, fBtait, f_bub_adv_src, f_divu)
        !$acc routine seq
        real(wp), intent(in) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu

        real(wp) :: fCpbw, fCpinf, fCpinf_dot, fH, fHdot, c_gas, c_liquid
        real(wp) :: f_rddot

        if (bubble_model == 1) then
            ! Gilmore bubbles
            fCpinf = fP - pref
            fCpbw = f_cpbw(fR0, fR, fV, fpb)
            fH = f_H(fCpbw, fCpinf, fntait, fBtait)
            c_gas = f_cgas(fCpinf, fntait, fBtait, fH)
            fCpinf_dot = f_cpinfdot(fRho, fP, alf, fntait, fBtait, f_bub_adv_src, f_divu)
            fHdot = f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
            f_rddot = f_rddot_G(fCpbw, fR, fV, fH, fHdot, c_gas, fntait, fBtait)
        else if (bubble_model == 2) then
            ! Keller-Miksis bubbles
            fCpinf = fP
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            c_liquid = dsqrt(fntait*(fP + fBtait)/(fRho*(1._wp - alf)))
            f_rddot = f_rddot_KM(fpbdot, fCpinf, fCpbw, fRho, fR, fV, fR0, c_liquid)
        else if (bubble_model == 3) then
            ! Rayleigh-Plesset bubbles
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            f_rddot = f_rddot_RP(fP, fRho, fR, fV, fR0, fCpbw)
        end if

    end function f_rddot

    !>  Function that computes the bubble radial acceleration for Rayleigh-Plesset bubbles
        !!  @param fCp Driving pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fCpbw Boundary wall pressure
    function f_rddot_RP(fCp, fRho, fR, fV, fR0, fCpbw)
        !$acc routine seq
        real(wp), intent(in) :: fCp, fRho, fR, fV, fR0, fCpbw

        real(wp) :: f_rddot_RP

            !! rddot = (1/r) (  -3/2 rdot^2 + ((r0/r)^3\gamma - Cp)/rho )
            !! rddot = (1/r) (  -3/2 rdot^2 + (tmp1 - Cp)/rho )
            !! rddot = (1/r) (  tmp2 )

        f_rddot_RP = (-1.5_wp*(fV**2._wp) + (fCpbw - fCp)/fRho)/fR

    end function f_rddot_RP

    !>  Function that computes the bubble radial acceleration
        !!  @param fCpbw Bubble wall pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fH Current enthalpy
        !!  @param fHdot Current time derivative of the enthalpy
        !!  @param fcgas Current gas sound speed
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    function f_rddot_G(fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait)
        !$acc routine seq
        real(wp), intent(in) :: fCpbw, fR, fV, fH, fHdot
        real(wp), intent(in) :: fcgas, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_rddot_G

        tmp1 = fV/fcgas
        tmp2 = 1._wp + 4._wp*Re_inv/fcgas/fR*(fCpbw/(1._wp + fBtait) + 1._wp) &
               **(-1._wp/fntait)
        tmp3 = 1.5_wp*fV**2._wp*(tmp1/3._wp - 1._wp) + fH*(1._wp + tmp1) &
               + fR*fHdot*(1._wp - tmp1)/fcgas

        f_rddot_G = tmp3/(fR*(1._wp - tmp1)*tmp2)

    end function f_rddot_G

    !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw_KM(fR0, fR, fV, fpb)
        !$acc routine seq
        real(wp), intent(in) :: fR0, fR, fV, fpb

        real(wp) :: f_cpbw_KM

        if (polytropic) then
            f_cpbw_KM = Ca*((fR0/fR)**(3._wp*gam)) - Ca + 1._wp
            if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM + &
                                                     (2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))
        else
            f_cpbw_KM = fpb
        end if

        if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM - 2._wp/(fR*Web)
        if (.not. f_is_default(Re_inv)) f_cpbw_KM = f_cpbw_KM - 4._wp*Re_inv*fV/fR

    end function f_cpbw_KM

    !>  Function that computes the bubble radial acceleration for Keller--Miksis bubbles
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param fCp Driving pressure
        !!  @param fCpbw Bubble wall pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fC Current sound speed
    function f_rddot_KM(fpbdot, fCp, fCpbw, fRho, fR, fV, fR0, fC)
        !$acc routine seq
        real(wp), intent(in) :: fpbdot, fCp, fCpbw
        real(wp), intent(in) :: fRho, fR, fV, fR0, fC

        real(wp) :: tmp1, tmp2, cdot_star
        real(wp) :: f_rddot_KM

        if (polytropic) then
            cdot_star = -3._wp*gam*Ca*((fR0/fR)**(3._wp*gam))*fV/fR
            if (.not. f_is_default(Web)) cdot_star = cdot_star - &
                                                     3._wp*gam*(2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))*fV/fR
        else
            cdot_star = fpbdot
        end if

        if (.not. f_is_default(Web)) cdot_star = cdot_star + (2._wp/Web)*fV/(fR**2._wp)
        if (.not. f_is_default(Re_inv)) cdot_star = cdot_star + 4._wp*Re_inv*((fV/fR)**2._wp)

        tmp1 = fV/fC
        tmp2 = 1.5_wp*(fV**2._wp)*(tmp1/3._wp - 1._wp) + &
               (1._wp + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (f_is_default(Re_inv)) then
            f_rddot_KM = tmp2/(fR*(1._wp - tmp1))
        else
            f_rddot_KM = tmp2/(fR*(1._wp - tmp1) + 4._wp*Re_inv/(fRho*fC))
        end if

    end function f_rddot_KM

    !>  Subroutine that computes bubble wall properties for vapor bubbles
        !!  @param pb Internal bubble pressure
        !!  @param iR0 Current bubble size index
    subroutine s_bwproperty(pb, iR0)
        !$acc routine seq
        real(kind(0._wp)), intent(in) :: pb
        integer, intent(in) :: iR0

        real(kind(0._wp)) :: x_vw

        ! mass fraction of vapor
        chi_vw = 1._wp/(1._wp + R_v/R_n*(pb/pv - 1._wp))
        ! mole fraction of vapor & thermal conductivity of gas mixture
        x_vw = M_n*chi_vw/(M_v + (M_n - M_v)*chi_vw)
        k_mw = x_vw*k_v(iR0)/(x_vw + (1._wp - x_vw)*phi_vn) &
               + (1._wp - x_vw)*k_n(iR0)/(x_vw*phi_nv + 1._wp - x_vw)
        ! gas mixture density
        rho_mw = pv/(chi_vw*R_v*Tw)

    end subroutine s_bwproperty

    !>  Function that computes the vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
    function f_vflux(fR, fV, fmass_v, iR0)
        !$acc routine seq
        real(kind(0._wp)), intent(in) :: fR
        real(kind(0._wp)), intent(in) :: fV
        real(kind(0._wp)), intent(in) :: fmass_v
        integer, intent(in) :: iR0

        real(kind(0._wp)) :: chi_bar
        real(kind(0._wp)) :: grad_chi
        real(kind(0._wp)) :: f_vflux

        if (thermal == 3) then !transfer
            ! constant transfer model
            chi_bar = fmass_v/(fmass_v + mass_n0(iR0))
            grad_chi = -Re_trans_c(iR0)*(chi_bar - chi_vw)
            f_vflux = rho_mw*grad_chi/Pe_c/(1._wp - chi_vw)/fR
        else
            ! polytropic
            f_vflux = pv*fV/(R_v*Tw)
        end if

    end function f_vflux

    !>  Function that computes the time derivative of
        !!  the internal bubble pressure
        !!  @param fvflux Vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Current internal bubble pressure
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
    function f_bpres_dot(fvflux, fR, fV, fpb, fmass_v, iR0)
        !$acc routine seq
        real(kind(0._wp)), intent(in) :: fvflux
        real(kind(0._wp)), intent(in) :: fR
        real(kind(0._wp)), intent(in) :: fV
        real(kind(0._wp)), intent(in) :: fpb
        real(kind(0._wp)), intent(in) :: fmass_v
        integer, intent(in) :: iR0

        real(kind(0._wp)) :: T_bar
        real(kind(0._wp)) :: grad_T
        real(kind(0._wp)) :: tmp1, tmp2
        real(kind(0._wp)) :: f_bpres_dot

        if (thermal == 3) then
            T_bar = Tw*(fpb/pb0(iR0))*(fR/R0(iR0))**3 &
                    *(mass_n0(iR0) + mass_v0(iR0))/(mass_n0(iR0) + fmass_v)
            grad_T = -Re_trans_T(iR0)*(T_bar - Tw)
            f_bpres_dot = 3._wp*gamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                         + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR)/fR
        else
            f_bpres_dot = -3._wp*gamma_m*fV/fR*(fpb - pv)
        end if

    end function f_bpres_dot

end module m_bubbles
