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

    implicit none

    private
    public :: s_bound_species_faces, s_clean_species, s_report_species_cleanup

    !> Violations below this fraction of rho are roundoff: cleaned but not counted
    real(wp), parameter :: clean_report_tol = 1.e-12_wp
    real(wp) :: clean_cells = 0._wp, clean_mass = 0._wp  !< Cells cleaned, and sum of |rho_new - rho_old|, since the last report

contains

    !> Bound the species face values vL/vR of direction id, over the range the reconstruction filled.
    subroutine s_bound_species_faces(q_prim_vf, vL, vR, id)

        type(scalar_field), dimension(sys_size), intent(in)                                    :: q_prim_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL, vR
        integer, intent(in)                                                                    :: id
        type(int_bounds_info), dimension(3)                                                    :: b
        real(wp)                                                                               :: w, theta, sL, sR, ybar, yL, yR, xi
        integer                                                                                :: i, j, k, l, polyn

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

    !> Fallback (as PeleC's clean_massfrac): clip each rho*Y_k to [0, rho], set rho = sum rho*Y_k, and scale momentum and energy by
    !! rho_new/rho_old so velocity and specific energy are kept. Counts for s_report_species_cleanup only cells off by >
    !! clean_report_tol*rho.
    subroutine s_clean_species(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp)                                               :: rho_old, rho_new, f, n_cells, d_mass, viol, rhoY
        integer                                                :: i, j, k, l

        n_cells = 0._wp; d_mass = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho_old, rho_new, f, viol, rhoY]', reduction='[[n_cells, d_mass]]', &
                            & reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho_old = real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp)
                    viol = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        rhoY = real(q_cons_vf(i)%sf(j, k, l), wp)
                        viol = max(viol, -rhoY, rhoY - rho_old)
                    end do
                    if (viol > 0._wp .and. rho_old > 0._wp) then
                        rho_new = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            rhoY = min(max(real(q_cons_vf(i)%sf(j, k, l), wp), 0._wp), rho_old)
                            q_cons_vf(i)%sf(j, k, l) = real(rhoY, stp)
                            rho_new = rho_new + rhoY
                        end do
                        f = rho_new/rho_old
                        q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l) = real(rho_new, stp)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%E
                            q_cons_vf(i)%sf(j, k, l) = real(f*real(q_cons_vf(i)%sf(j, k, l), wp), stp)
                        end do
                        if (viol > clean_report_tol*rho_old) n_cells = n_cells + 1._wp
                        d_mass = d_mass + abs(rho_new - rho_old)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        clean_cells = clean_cells + n_cells; clean_mass = clean_mass + d_mass

    end subroutine s_clean_species

    !> Print, once per time step and only if nonzero, how many cells s_clean_species changed across all ranks.
    impure subroutine s_report_species_cleanup(t_step)

        integer, intent(in) :: t_step
        real(wp)            :: cells_glb, mass_glb

        call s_mpi_allreduce_sum(clean_cells, cells_glb)
        call s_mpi_allreduce_sum(clean_mass, mass_glb)
        if (proc_rank == 0 .and. cells_glb > 0._wp) print '(A,I0,A,I0,A,ES10.3)', 'Species cleanup at step ', t_step, ': ', &
            & nint(cells_glb), ' cells, sum |d rho| = ', mass_glb
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
