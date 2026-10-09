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

    implicit none

    private
    public :: s_bound_species_faces

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
