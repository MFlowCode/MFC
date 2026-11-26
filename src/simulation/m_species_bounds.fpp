!>
!! @file
!! @brief Contains module m_species_bounds

#:include 'macros.fpp'

!> @brief Keeps reconstructed species mass fractions on the simplex (Y >= 0, sum Y = 1). Each cell's face values are first shifted
!! to sum to 1, then scaled toward the cell average by one theta shared by all species and both faces (Zhang & Shu, JCP 2010). With
!! face Y on the simplex the HLLC species fluxes are Y_upwind times the mass flux (Larrouturou, JCP 1991), so they sum to it and
!! keep cell averages non-negative under the Zhang-Shu CFL.
module m_species_bounds

    use m_derived_types
    use m_global_parameters
    use m_mpi_common, only: s_mpi_allreduce_sum
    use m_constants, only: T_surface_max, verysmall
    use m_thermochem, only: gas_constant, get_mixture_molecular_weight
    use m_ibm, only: ib_markers, ib_markers_fine

    implicit none

    private
    public :: s_bound_species_faces, s_bound_face_states, s_bound_face_states_zs, s_clean_species, s_report_species_cleanup

    !> face_bound = 2 caps a face's (p + p_inf)/rho at this multiple of the largest over the cell and its neighbours (sound speed at
    !! ~32x): high enough never to act on a resolved smooth field, low enough to rule out a floor-density face at finite pressure
    real(wp), parameter :: face_T_cap = 1000._wp
    !> Relative positivity floor of the face bounds: verysmall, raised to 100 ulp where wp cannot resolve it (single precision); a
    !! floor below roundoff lets the bounded face land on the wrong side of zero
    real(wp), parameter :: pos_floor = max(verysmall, 100._wp*epsilon(1._wp))
    !> Violations below this fraction of rho are roundoff: cleaned but not counted
    real(wp), parameter :: clean_report_tol = 1.e-12_wp
    !> Fluid-cell cleanings, and sum of |rho_new - rho_old|, since the last report, per AMR level (0 = the coarse grid)
    real(wp), allocatable :: clean_cells(:), clean_mass(:)

contains

    !> Bound the species face values vL/vR of direction id, over the range the reconstruction filled.
    subroutine s_bound_species_faces(q_prim_vf, vL, vR, id)

        type(scalar_field), dimension(sys_size), intent(in)                                    :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL, vR
        integer, intent(in)                                                                    :: id
        type(int_bounds_info), dimension(3)                                                    :: b
        real(wp)                                                                               :: w, theta, sL, sR, ybar, yL, yR, xi
        integer                                                                                :: i, j, k, l, polyn

        ! Same range as s_reconstruct_cell_boundary_values (m_rhs.fpp)

        polyn = merge(weno_polyn, muscl_polyn, recon_type == recon_type_weno)
        b = idwbuff
        b(id)%beg = b(id)%beg + polyn; b(id)%end = b(id)%end - polyn
        w = merge(f_lobatto_weight(polyn), 0.5_wp, recon_type == recon_type_weno)  ! MUSCL is linear: no interior point

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, theta, sL, sR, ybar, yL, yR, xi]', copyin='[b, w]')
        do l = b(3)%beg, b(3)%end
            do k = b(2)%beg, b(2)%end
                do j = b(1)%beg, b(1)%end
                    sL = 0._wp; sR = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        sL = sL + vL(j, k, l, i); sR = sR + vR(j, k, l, i)
                    end do
                    theta = 1._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        ybar = q_prim_vf(i)%sf(j, k, l)
                        yL = vL(j, k, l, i) + ybar*(1._wp - sL)
                        yR = vR(j, k, l, i) + ybar*(1._wp - sR)
                        vL(j, k, l, i) = yL; vR(j, k, l, i) = yR
                        theta = min(theta, f_theta(ybar, yL), f_theta(ybar, yR))
                        if (w < 0.5_wp) then
                            xi = (ybar - w*(yL + yR))/(1._wp - 2._wp*w)
                            theta = min(theta, f_theta(ybar, xi))
                        end if
                    end do
                    if (theta < 1._wp) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            ybar = q_prim_vf(i)%sf(j, k, l)
                            vL(j, k, l, i) = ybar + theta*(vL(j, k, l, i) - ybar)
                            vR(j, k, l, i) = ybar + theta*(vR(j, k, l, i) - ybar)
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_bound_species_faces

    !> Scales each cell's primitive face states vL/vR of direction id toward its cell state by one theta (Zhang & Shu, JCP 2010) so
    !! that face rho and p + p_inf stay positive and, with chemistry, the face gas constant stays positive and the face T <=
    !! T_surface_max (past which cp can go negative and the sound speed NaN). Inactive where no face state leaves those bounds.
    subroutine s_bound_face_states(q_prim_vf, vL, vR, id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL, vR
        integer, intent(in) :: id
        type(int_bounds_info), dimension(3) :: b
        real(wp), dimension(${NUM_SPECIES}$) :: Ys
        real(wp), dimension(3) :: vb, wL, wR  !< (rho, p + p_inf, R_mix) of the cell and its two faces
        real(wp) :: theta, pinf, mw
        integer :: i, j, k, l, polyn

        ! Same range as s_reconstruct_cell_boundary_values (m_rhs.fpp)

        polyn = merge(weno_polyn, muscl_polyn, recon_type == recon_type_weno)
        b = idwbuff
        b(id)%beg = b(id)%beg + polyn; b(id)%end = b(id)%end - polyn
        pinf = maxval(pi_infs(1:num_fluids)/(gammas(1:num_fluids) + 1._wp))  ! stiffest fluid: liquid tension is not limited

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, Ys, vb, wL, wR, theta, mw]', copyin='[b, pinf]')
        do l = b(3)%beg, b(3)%end
            do k = b(2)%beg, b(2)%end
                do j = b(1)%beg, b(1)%end
                    #:for V, Q in [('vb', 'q_prim_vf({})%sf(j, k, l)'), ('wL', 'vL(j, k, l, {})'), ('wR', 'vR(j, k, l, {})')]
                        ${V}$(1) = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%cont%beg, eqn_idx%cont%end
                            ${V}$(1) = ${V}$(1) + ${Q.format('i')}$
                        end do
                        ${V}$(2) = ${Q.format('eqn_idx%E')}$ + pinf
                        if (chemistry) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_species
                                Ys(i) = ${Q.format('eqn_idx%species%beg + i - 1')}$
                            end do
                            call get_mixture_molecular_weight(Ys, mw)
                            ${V}$(3) = gas_constant/mw
                        end if
                    #:endfor
                    ! rho, p and R_mix are linear in theta, so positivity is closed form and T <= T_max quadratic
                    theta = min(1._wp, f_theta_floor(vb(1), wL(1)), f_theta_floor(vb(1), wR(1)), f_theta_floor(vb(2), wL(2)), &
                                & f_theta_floor(vb(2), wR(2)))
                    if (chemistry) then
                        theta = min(theta, f_theta_floor(vb(3), wL(3)), f_theta_floor(vb(3), wR(3)))
                        theta = min(f_theta_T(theta, vb, wL), f_theta_T(theta, vb, wR))
                    end if
                    if (theta < 1._wp) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, sys_size
                            vL(j, k, l, i) = q_prim_vf(i)%sf(j, k, l) + theta*(vL(j, k, l, i) - q_prim_vf(i)%sf(j, k, l))
                            vR(j, k, l, i) = q_prim_vf(i)%sf(j, k, l) + theta*(vR(j, k, l, i) - q_prim_vf(i)%sf(j, k, l))
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_bound_face_states

    !> face_bound = 2, for stiffened-gas 5- and 6-equation states with one or more fluids (Zhang & Shu, JCP 2010). Each cell's
    !! primitive face states vL/vR of direction id are scaled toward its cell state by the largest theta that keeps, on both faces,
    !! rho and p + p_inf positive and (p + p_inf)/rho at most face_T_cap times its largest value over the cell and its two
    !! neighbours along id. The last bound limits the face sound speed: without it, a face scaled down to a density floor at finite
    !! pressure sets an unbounded wave speed. theta is then halved until the conservative interior point (U_bar - w(U_L + U_R))/(1 -
    !! 2w), with w the Gauss-Lobatto end weight, has positive density and pressure. Volume fractions get their own theta into [0,
    !! 1]. The limiter is inactive wherever every bound already holds. Cell averages then stay positive for dt sum_d max(|u_d| +
    !! c)/dx_d <= w.
    subroutine s_bound_face_states_zs(q_prim_vf, vL, vR, id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL, vR
        integer, intent(in) :: id
        type(int_bounds_info), dimension(3) :: b
        real(wp), dimension(3) :: mx  !< interior-point momentum
        real(wp) :: w, pinf, theta, theta_a, rb, pb, tc, rf, pf, ab, s, r, rx, ex, gx, pix, qvx, gam, pim, qv, e, u, a, gb
        integer :: i, j, k, l, o, it, dj, dk, dl, polyn
        logical :: ok

        ! Range of s_reconstruct_cell_boundary_values (m_rhs.fpp), narrowed by at least 1 for the +-1 neighbours along id

        polyn = merge(weno_polyn, muscl_polyn, recon_type == recon_type_weno)
        b = idwbuff
        b(id)%beg = b(id)%beg + max(polyn, 1); b(id)%end = b(id)%end - max(polyn, 1)
        w = merge(f_lobatto_weight(polyn), 0.5_wp, recon_type == recon_type_weno)  ! MUSCL is linear: no interior point
        pinf = maxval(pi_infs(1:num_fluids)/(gammas(1:num_fluids) + 1._wp))  ! stiffest fluid: liquid tension is not limited
        dj = merge(1, 0, id == 1); dk = merge(1, 0, id == 2); dl = merge(1, 0, id == 3)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, o, it, mx, theta, theta_a, rb, pb, tc, rf, pf, ab, s, r, rx, ex, &
                            & gx, pix, qvx, gam, pim, qv, e, u, a, gb, ok]', copyin='[b, w, pinf, dj, dk, dl]')
        do l = b(3)%beg, b(3)%end
            do k = b(2)%beg, b(2)%end
                do j = b(1)%beg, b(1)%end
                    tc = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do o = -1, 1
                        r = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%cont%beg, eqn_idx%cont%end
                            r = r + q_prim_vf(i)%sf(j + o*dj, k + o*dk, l + o*dl)
                        end do
                        tc = max(tc, (q_prim_vf(eqn_idx%E)%sf(j + o*dj, k + o*dk, l + o*dl) + pinf)/max(r, sgm_eps))
                    end do
                    tc = face_T_cap*tc

                    rb = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%cont%beg, eqn_idx%cont%end
                        rb = rb + q_prim_vf(i)%sf(j, k, l)
                    end do
                    pb = q_prim_vf(eqn_idx%E)%sf(j, k, l) + pinf
                    theta = 1._wp
                    #:for V in ['vL', 'vR']
                        rf = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%cont%beg, eqn_idx%cont%end
                            rf = rf + ${V}$(j, k, l, i)
                        end do
                        pf = ${V}$(j, k, l, eqn_idx%E) + pinf
                        theta = min(theta, f_theta_floor(rb, rf), f_theta_floor(pb, pf), f_theta_floor(tc*rb - pb, tc*rf - pf))
                    #:endfor

                    theta_a = 1._wp
                    if (num_fluids > 1) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%beg + num_fluids - 1
                            ab = q_prim_vf(i)%sf(j, k, l)
                            theta_a = min(theta_a, f_theta_floor(ab, vL(j, k, l, i)), f_theta_floor(1._wp - ab, 1._wp - vL(j, k, &
                                          & l, i)), f_theta_floor(ab, vR(j, k, l, i)), f_theta_floor(1._wp - ab, 1._wp - vR(j, k, &
                                          & l, i)))
                        end do
                    end if

                    if (w < 0.5_wp) then
                        ! Interior point from the reduced conservative state (rho, m, E, Gamma, Pi_inf, q_v) of the cell and faces
                        gb = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            gb = gb + q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l)*gammas(i)
                        end do
                        gb = gb*pb  ! Gamma (p + p_inf)
                        ok = .false.
                        $:GPU_LOOP(parallelism='[seq]')
                        do it = 1, 30
                            ! xi*(1 - 2w) = U_bar - w(U_L + U_R); rho and g = E - |m|^2/(2 rho) - Pi_inf - q_v + Gamma p_inf are
                            ! homogeneous of degree one, so the positive factor 1 - 2w is kept on the bounds instead
                            rx = 0._wp; ex = 0._wp; gx = 0._wp; pix = 0._wp; qvx = 0._wp; mx = 0._wp
                            #:for V in ['', 'vL', 'vR']
                                #:set C = '1._wp' if V == '' else '(-w)'
                                #:if V == ''
                                    #:set P = 'q_prim_vf({0})%sf(j, k, l)'
                                #:else
                                    #:set P = '(q_prim_vf({0})%sf(j, k, l) + {1}*(' + V &
                                        & + '(j, k, l, {0}) - q_prim_vf({0})%sf(j, k, l)))'
                                #:endif
                                r = 0._wp; qv = 0._wp; gam = 0._wp; pim = 0._wp; e = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    a = ${P.format('eqn_idx%cont%beg + i - 1', 'theta')}$
                                    r = r + a; qv = qv + a*qvs(i)
                                    a = ${P.format('eqn_idx%adv%beg + i - 1', 'theta*theta_a')}$
                                    gam = gam + a*gammas(i); pim = pim + a*pi_infs(i)
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    u = ${P.format('eqn_idx%mom%beg + i - 1', 'theta')}$
                                    mx(i) = mx(i) + ${C}$*r*u
                                    e = e + 0.5_wp*r*u*u
                                end do
                                e = e + gam*${P.format('eqn_idx%E', 'theta')}$ + pim + qv
                                rx = rx + ${C}$*r; ex = ex + ${C}$*e; gx = gx + ${C}$*gam; pix = pix + ${C}$*pim; qvx = qvx &
                                    & + ${C}$*qv
                            #:endfor
                            if (rx >= pos_floor*(1._wp - 2._wp*w)*rb) then
                                ok = ex - 0.5_wp*(mx(1)**2 + mx(2)**2 + mx(3)**2)/rx - pix - qvx + gx*pinf >= pos_floor*(1._wp &
                                                  & - 2._wp*w)*gb
                            end if
                            if (ok) exit
                            theta = 0.5_wp*theta
                        end do
                        if (.not. ok) theta = 0._wp
                    end if

                    if (theta < 1._wp .or. theta_a < 1._wp) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, sys_size
                            s = theta
                            if (i >= eqn_idx%adv%beg .and. i < eqn_idx%adv%beg + num_fluids) s = theta*theta_a
                            vL(j, k, l, i) = q_prim_vf(i)%sf(j, k, l) + s*(vL(j, k, l, i) - q_prim_vf(i)%sf(j, k, l))
                            vR(j, k, l, i) = q_prim_vf(i)%sf(j, k, l) + s*(vR(j, k, l, i) - q_prim_vf(i)%sf(j, k, l))
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_bound_face_states_zs

    !> Largest theta in [0, 1] keeping vbar + theta*(v - vbar) >= pos_floor*vbar (0 if vbar <= 0).
    pure function f_theta_floor(vbar, v) result(theta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: vbar, v
        real(wp)             :: theta

        theta = 1._wp
        if (v < pos_floor*vbar) theta = max(0._wp, (1._wp - pos_floor)*vbar)/max(vbar - v, sgm_eps)

    end function f_theta_floor

    !> Largest theta' <= theta keeping the face (rho, p + p_inf, R_mix) scaled toward vb by theta' at T <= T_surface_max (0 if vb is
    !! above it; p_inf is zero for the gases that carry chemistry). g = T_max rho R_mix - p is quadratic in theta'; where g(theta) <
    !! 0 its first root is the stable smaller-root form.
    pure function f_theta_T(theta, vb, vf) result(t)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)               :: theta
        real(wp), dimension(3), intent(in) :: vb, vf
        real(wp), dimension(3)             :: d
        real(wp)                           :: a, b, c, t

        d = vf - vb
        a = T_surface_max*d(1)*d(3)
        b = T_surface_max*(vb(1)*d(3) + vb(3)*d(1)) - d(2)
        c = T_surface_max*vb(1)*vb(3) - vb(2)
        t = theta
        if ((a*theta + b)*theta + c < 0._wp) t = max(0._wp, 2._wp*c/max(sqrt(max(b*b - 4._wp*a*c, 0._wp)) - b, sgm_eps))

    end function f_theta_T

    !> Fallback (as PeleC's clean_massfrac): clip each rho*Y_k to [0, rho], set rho = sum rho*Y_k, and scale momentum and energy by
    !! rho_new/rho_old so velocity and specific energy are kept. A cell whose clipped sum is not positive is left untouched. Every
    !! cell is cleaned, but only fluid cells (not IB cells) off by > clean_report_tol*rho, or left untouched, are counted as cell
    !! cleanings (per cell per RK stage) for s_report_species_cleanup, under AMR level lvl (0 = the coarse grid). A fine block's
    !! markers are in ib_markers_fine.
    subroutine s_clean_species(q_cons_vf, lvl)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: lvl
        real(wp)                                               :: rho_old, rho_new, f, n_cells, d_mass, viol, rhoY
        logical                                                :: is_fluid
        integer                                                :: i, j, k, l

        n_cells = 0._wp; d_mass = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho_old, rho_new, f, viol, rhoY, is_fluid]', reduction='[[n_cells, &
                            & d_mass]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho_old = real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp)
                    viol = 0._wp
                    rho_new = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        rhoY = real(q_cons_vf(i)%sf(j, k, l), wp)
                        viol = max(viol, -rhoY, rhoY - rho_old)
                        rho_new = rho_new + min(max(rhoY, 0._wp), rho_old)
                    end do
                    if (viol > 0._wp .and. rho_old > 0._wp) then
                        is_fluid = .true.
                        if (ib .and. lvl == 0) is_fluid = (ib_markers%sf(j, k, l) == 0)
                        if (ib .and. lvl > 0) is_fluid = (ib_markers_fine%sf(j, k, l) == 0)
                        if (rho_new > 0._wp) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                rhoY = min(max(real(q_cons_vf(i)%sf(j, k, l), wp), 0._wp), rho_old)
                                q_cons_vf(i)%sf(j, k, l) = real(rhoY, stp)
                            end do
                            f = rho_new/rho_old
                            q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l) = real(rho_new, stp)
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%mom%beg, eqn_idx%E
                                q_cons_vf(i)%sf(j, k, l) = real(f*real(q_cons_vf(i)%sf(j, k, l), wp), stp)
                            end do
                            if (is_fluid) then
                                if (viol > clean_report_tol*rho_old) n_cells = n_cells + 1._wp
                                d_mass = d_mass + abs(rho_new - rho_old)
                            end if
                        else if (is_fluid) then
                            n_cells = n_cells + 1._wp
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (.not. allocated(clean_cells)) then
            allocate (clean_cells(0:amr_max_level), clean_mass(0:amr_max_level))
            clean_cells = 0._wp; clean_mass = 0._wp
        end if
        clean_cells(lvl) = clean_cells(lvl) + n_cells; clean_mass(lvl) = clean_mass(lvl) + d_mass

    end subroutine s_clean_species

    !> Print, once per time step and only if nonzero, how many cell cleanings s_clean_species made across all ranks, per AMR level.
    !! Coarse cells under a fine block are counted on level 0 too.
    impure subroutine s_report_species_cleanup(t_step)

        integer, intent(in) :: t_step
        real(wp)            :: cells_glb, mass_glb
        integer             :: lvl
        character(len=16)   :: tag

        if (.not. allocated(clean_cells)) return
        do lvl = 0, amr_max_level
            call s_mpi_allreduce_sum(clean_cells(lvl), cells_glb)
            call s_mpi_allreduce_sum(clean_mass(lvl), mass_glb)
            tag = ''
            if (lvl > 0) write (tag, '(A,I0)') ', AMR level ', lvl
            if (proc_rank == 0 .and. cells_glb > 0._wp) print '(A,I0,A,I0,A,ES10.3)', 'Species cleanup at step ', t_step, &
                & trim(tag) // ': ', nint(cells_glb), ' cell cleanings, sum |d rho| = ', mass_glb
        end do
        clean_cells = 0._wp; clean_mass = 0._wp

    end subroutine s_report_species_cleanup

    !> Largest theta in [0, 1] keeping ybar + theta*(y - ybar) >= 0, given ybar >= 0.
    pure function f_theta(ybar, y) result(theta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: ybar, y
        real(wp)             :: theta

        theta = 1._wp
        if (y < 0._wp) theta = max(0._wp, ybar)/max(ybar - y, sgm_eps)

    end function f_theta

    !> Endpoint weight of the Gauss-Lobatto rule exact for the reconstruction's degree 2*polyn: 1/(N(N-1)), N = polyn + 2.
    pure function f_lobatto_weight(polyn) result(w)

        integer, intent(in) :: polyn
        real(wp)            :: w

        w = 1._wp/real((polyn + 2)*(polyn + 1), wp)

    end function f_lobatto_weight

end module m_species_bounds
