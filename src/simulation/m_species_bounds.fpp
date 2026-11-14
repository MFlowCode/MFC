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
    use m_ibm, only: ib_markers

    implicit none

    private
    public :: s_bound_species_faces, s_bound_face_states, s_clean_species, s_report_species_cleanup

    !> Violations below this fraction of rho are roundoff: cleaned but not counted
    real(wp), parameter :: clean_report_tol = 1.e-12_wp
    !> Fluid-cell cleanings, and sum of |rho_new - rho_old|, since the last report
    real(wp) :: clean_cells = 0._wp, clean_mass = 0._wp

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

    !> Largest theta in [0, 1] keeping vbar + theta*(v - vbar) >= verysmall*vbar (0 if vbar <= 0).
    pure function f_theta_floor(vbar, v) result(theta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: vbar, v
        real(wp)             :: theta

        theta = 1._wp
        if (v < verysmall*vbar) theta = max(0._wp, (1._wp - verysmall)*vbar)/max(vbar - v, sgm_eps)

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
    !! cleanings (per cell per RK stage) for s_report_species_cleanup. tally = F (an AMR fine block, whose markers ib_markers does
    !! not hold) cleans without counting.
    subroutine s_clean_species(q_cons_vf, tally)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        logical, intent(in), optional                          :: tally
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
                        if (ib) is_fluid = (ib_markers%sf(j, k, l) == 0)
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
        if (present(tally)) then
            if (.not. tally) return
        end if
        clean_cells = clean_cells + n_cells; clean_mass = clean_mass + d_mass

    end subroutine s_clean_species

    !> Print, once per time step and only if nonzero, how many cell cleanings s_clean_species made across all ranks.
    impure subroutine s_report_species_cleanup(t_step)

        integer, intent(in) :: t_step
        real(wp)            :: cells_glb, mass_glb

        call s_mpi_allreduce_sum(clean_cells, cells_glb)
        call s_mpi_allreduce_sum(clean_mass, mass_glb)
        if (proc_rank == 0 .and. cells_glb > 0._wp) print '(A,I0,A,I0,A,ES10.3)', 'Species cleanup at step ', t_step, ': ', &
            & nint(cells_glb), ' cell cleanings, sum |d rho| = ', mass_glb
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
