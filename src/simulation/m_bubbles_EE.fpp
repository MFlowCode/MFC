!>
!! @file m_bubbles_EE.f90
!! @brief Contains module m_bubbles_EE

#:include 'macros.fpp'

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
module m_bubbles_EE

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< General bubble dynamics procedures

    ! ==========================================================================

    implicit none

    real(kind(0d0)), allocatable, dimension(:, :, :) :: bub_adv_src
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    !$acc declare create(bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

    type(scalar_field) :: divu !< matrix for div(u)
    !$acc declare create(divu)

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    !$acc declare create(rs, vs, ms, ps)

contains

    subroutine s_initialize_bubbles_EE_module

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

    end subroutine s_initialize_bubbles_EE_module

    ! Compute the bubble volume fraction alpha from the bubble number density n
        !! @param q_cons_vf is the conservative variable
    subroutine s_comp_alpha_from_n(q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(kind(0d0)) :: nR3bar
        integer(kind(0d0)) :: i, j, k, l

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    nR3bar = 0d0
                    !$acc loop seq
                    do i = 1, nb
                        nR3bar = nR3bar + weight(i)*(q_cons_vf(rs(i))%sf(j, k, l))**3d0
                    end do
                    q_cons_vf(alf_idx)%sf(j, k, l) = (4d0*pi*nR3bar)/(3d0*q_cons_vf(n_idx)%sf(j, k, l)**2d0)
                end do
            end do
        end do

    end subroutine s_comp_alpha_from_n

    subroutine s_compute_bubbles_EE_rhs(idir, q_prim_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: j, k, l

        if (idir == 1) then

            if (.not. qbmm) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu%sf(j, k, l) = 0d0
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

    end subroutine s_compute_bubbles_EE_rhs

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the bubble modeling
        !!  @param q_prim_vf Primitive variables
        !!  @param q_cons_vf Conservative variables
    subroutine s_compute_bubble_EE_source(q_cons_vf, q_prim_vf, t_step, rhs_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rddot
        real(kind(0d0)) :: pb, mv, vflux, pbdot
        real(kind(0d0)) :: n_tait, B_tait
        real(kind(0d0)), dimension(nb) :: Rtmp, Vtmp
        real(kind(0d0)) :: myR, myV, alf, myP, myRho, R2Vav, R3
        real(kind(0d0)), dimension(num_fluids) :: myalpha, myalpha_rho
        real(kind(0d0)) :: nbub !< Bubble number density
        integer :: i, j, k, l, q, ii !< Loop variables

        real(kind(0d0)) :: err1, err2, err3, err4, err5 !< Error estimates for adaptive time stepping
        real(kind(0d0)) :: t_new !< Updated time step size
        real(kind(0d0)) :: h !< Time step size
        real(kind(0d0)), dimension(4) :: myR_tmp1, myV_tmp1, myR_tmp2, myV_tmp2 !< Bubble radius, radial velocity, and radial acceleration for the inner loop

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    bub_adv_src(j, k, l) = 0d0

                    !$acc loop seq
                    do q = 1, nb
                        bub_r_src(j, k, l, q) = 0d0
                        bub_v_src(j, k, l, q) = 0d0
                        bub_p_src(j, k, l, q) = 0d0
                        bub_m_src(j, k, l, q) = 0d0
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

                        R3 = 0d0

                        !$acc loop seq
                        do q = 1, nb
                            R3 = R3 + weight(q)*Rtmp(q)**3.d0
                        end do

                        nbub = (3.d0/(4.d0*pi))*q_prim_vf(alf_idx)%sf(j, k, l)/R3
                    end if

                    if (.not. adap_dt) then
                        R2Vav = 0d0

                        !$acc loop seq
                        do q = 1, nb
                            R2Vav = R2Vav + weight(q)*Rtmp(q)**2.d0*Vtmp(q)
                        end do

                        bub_adv_src(j, k, l) = 4.d0*pi*nbub*R2Vav
                    end if

                    !$acc loop seq
                    do q = 1, nb

                        !$acc loop seq
                        do ii = 1, num_fluids
                            myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                            myalpha(ii) = q_cons_vf(advxb + ii - 1)%sf(j, k, l)
                        end do

                        myRho = 0d0
                        n_tait = 0d0
                        B_tait = 0d0

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

                        n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'
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
                            bub_m_src(j, k, l, q) = nbub*vflux*4.d0*pi*(myR**2.d0)
                        else
                            pb = 0d0; mv = 0d0; vflux = 0d0; pbdot = 0d0
                        end if

                        ! Adaptive time stepping
                        if (adap_dt) then
                            ! Determine the starting time step
                            call s_initialize_adap_dt(myRho, myP, myR, myV, R0(q), &
                                                      pb, pbdot, alf, n_tait, B_tait, &
                                                      bub_adv_src(j, k, l), divu%sf(j, k, l), h)

                            ! Advancing one step
                            t_new = 0d0
                            do while (.true.)
                                if (t_new + h > 0.5d0*dt) then
                                    h = 0.5d0*dt - t_new
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
                                                           bub_adv_src(j, k, l), divu%sf(j, k, l), 0.5d0*h, &
                                                           myR_tmp2, myV_tmp2, err2)

                                    call s_advance_substep(myRho, myP, myR_tmp2(4), myV_tmp2(4), R0(q), &
                                                           pb, pbdot, alf, n_tait, B_tait, &
                                                           bub_adv_src(j, k, l), divu%sf(j, k, l), 0.5d0*h, &
                                                           myR_tmp2, myV_tmp2, err3)

                                    err4 = abs((myR_tmp1(4) - myR_tmp2(4))/myR_tmp1(4))
                                    err5 = abs((myV_tmp1(4) - myV_tmp2(4))/myV_tmp1(4))
                                    if (abs(myV_tmp1(4)) < 1e-12) err5 = 0d0

                                    ! Determine acceptance/rejection and update step size
                                    !   Rule 1: err1, err2, err3 < tol
                                    !   Rule 2: myR_tmp1(4) > 0d0
                                    !   Rule 3: abs((myR_tmp1(4) - myR_tmp2(4))/myR) < tol
                                    !   Rule 4: abs((myV_tmp1(4) - myV_tmp2(4))/myV) < tol
                                    if ((err1 <= 1d-4) .and. (err2 <= 1d-4) .and. (err3 <= 1d-4) &
                                        .and. (err4 < 1d-4) .and. (err5 < 1d-4) &
                                        .and. myR_tmp1(4) > 0d0) then

                                        ! Accepted. Finalize the sub-step
                                        t_new = t_new + h

                                        ! Update R and V
                                        myR = myR_tmp1(4)
                                        myV = myV_tmp1(4)

                                        ! Update step size for the next sub-step
                                        h = h*min(2d0, max(0.5d0, (1d-4/err1)**(1d0/3d0)))

                                        exit
                                    else
                                        ! Rejected. Update step size for the next try on sub-step
                                        if (err2 <= 1d-4) then
                                            h = 0.5d0*h
                                        else
                                            h = 0.25d0*h
                                        end if

                                    end if
                                end do

                                ! Exit the loop if the final time reached dt
                                if (t_new == 0.5d0*dt) exit

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
                            bub_adv_src(j, k, l) = 0d0
                            bub_r_src(j, k, l, q) = 0d0
                            bub_v_src(j, k, l, q) = 0d0
                            if (.not. polytropic) then
                                bub_p_src(j, k, l, q) = 0d0
                                bub_m_src(j, k, l, q) = 0d0
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
    end subroutine s_compute_bubble_EE_source

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
        real(kind(0d0)), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(kind(0d0)), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(kind(0d0)), intent(out) :: h

        real(kind(0d0)) :: h0, h1 !< Time step size
        real(kind(0d0)) :: d0, d1, d2 !< norms
        real(kind(0d0)), dimension(2) :: myR_tmp, myV_tmp, myA_tmp !< Bubble radius, radial velocity, and radial acceleration

        ! Determine the starting time step
        ! Evaluate f(x0,y0)
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Compute d0 = ||y0|| and d1 = ||f(x0,y0)||
        d0 = DSQRT((myR_tmp(1)**2d0 + myV_tmp(1)**2d0)/2d0)
        d1 = DSQRT((myV_tmp(1)**2d0 + myA_tmp(1)**2d0)/2d0)
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
        d2 = DSQRT(((myV_tmp(2) - myV_tmp(1))**2d0 + (myA_tmp(2) - myA_tmp(1))**2d0)/2d0)/h0

        ! Set h1 = (0.01/max(d1,d2))^{1/(p+1)}
        !      if max(d1,d2) < 1e-15, h1 = max(1e-6, h0*1e-3)
        if (max(d1, d2) < 1d-15) then
            h1 = max(1d-6, h0*1d-3)
        else
            h1 = (1d-2/max(d1, d2))**(1d0/3d0)
        end if

        ! Set h = min(100*h0,h1)
        h = min(100d0*h0, h1)

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
        real(kind(0d0)), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(kind(0d0)), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu, h
        real(kind(0d0)), dimension(4), intent(OUT) :: myR_tmp, myV_tmp
        real(kind(0d0)), dimension(4) :: myA_tmp
        real(kind(0d0)), intent(OUT) :: err
        real(kind(0d0)) :: err_R, err_V

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
        myR_tmp(3) = myR_tmp(1) + (h/4d0)*(myV_tmp(1) + myV_tmp(2))
        myV_tmp(3) = myV_tmp(1) + (h/4d0)*(myA_tmp(1) + myA_tmp(2))
        myA_tmp(3) = f_rddot(fRho, fP, myR_tmp(3), myV_tmp(3), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Stage 3
        myR_tmp(4) = myR_tmp(1) + (h/6d0)*(myV_tmp(1) + myV_tmp(2) + 4d0*myV_tmp(3))
        myV_tmp(4) = myV_tmp(1) + (h/6d0)*(myA_tmp(1) + myA_tmp(2) + 4d0*myA_tmp(3))
        myA_tmp(4) = f_rddot(fRho, fP, myR_tmp(4), myV_tmp(4), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu)

        ! Estimate error
        err_R = (-5d0*h/24d0)*(myV_tmp(2) + myV_tmp(3) - 2d0*myV_tmp(4)) &
                /max(abs(myR_tmp(1)), abs(myR_tmp(4)))
        err_V = (-5d0*h/24d0)*(myA_tmp(2) + myA_tmp(3) - 2d0*myA_tmp(4)) &
                /max(abs(myV_tmp(1)), abs(myV_tmp(4)))
        err = DSQRT((err_R**2d0 + err_V**2d0)/2d0)

    end subroutine s_advance_substep

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
        real(kind(0d0)), intent(in) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(kind(0d0)), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu

        real(kind(0d0)) :: fCpbw, fCpinf, fCpinf_dot, fH, fHdot, c_gas, c_liquid
        real(kind(0d0)) :: f_rddot

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
            c_liquid = dsqrt(fntait*(fP + fBtait)/(fRho*(1.d0 - alf)))
            f_rddot = f_rddot_KM(fpbdot, fCpinf, fCpbw, fRho, fR, fV, fR0, c_liquid)
        else if (bubble_model == 3) then
            ! Rayleigh-Plesset bubbles
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            f_rddot = f_rddot_RP(fP, fRho, fR, fV, fR0, fCpbw)
        end if

    end function f_rddot

end module m_bubbles_EE
