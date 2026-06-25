!>
!! @file m_acoustic_substep.fpp
!! @brief Contains module m_acoustic_substep

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Cheap forward-backward acoustic substep kernel for the split-explicit low-Mach time integrator (acoustic_substepping).
!!
!! In split mode the HLLC solver returns only the slow (advective) flux:
!! contact-upwinded momentum advection and the volume-fraction flux; the mass
!! and total-energy face fluxes are zeroed there. This module advances the
!! conserved state by `n_micro` cheap second-order forward-backward microsteps
!! that transport mass and total energy and perform the acoustic (pressure)
!! work, with the slow forcing held frozen.
!!
!! Per microstep (5-equation model, total-energy formulation, step size dtau):
!!   1. narrow halo exchange of the conserved field
!!   2. forward : advance alpha_k*rho_k and rho*E by the second-order centered
!!                divergence of (alpha_k*rho_k * u) and ((rho*E + p) * u), using
!!                the current (frozen-this-microstep) velocity u, plus
!!                dtau * rhs_slow
!!   3. EOS     : recompute the cell pressure p from the updated conserved state
!!   4. backward: advance rho*u by  -dtau*grad(p_new) + dtau*rhs_slow_mom
!!                + dtau*acoustic_div_damp*grad(div(rho*u))   (divergence damping)
!!
!! The mixture coefficients (gamma, pi_inf) come from the FROZEN volume
!! fractions (the alpha_k are pure transport / slow forcing and are not advanced
!! here), so the EOS pressure is cell-local and cheap.
module m_acoustic_substep

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_boundary_common

    implicit none

    private
    public :: s_initialize_acoustic_substep_module, s_acoustic_substep, s_finalize_acoustic_substep_module

    !> Cell pressure scratch (working precision), spanning the buffered domain.
    real(wp), allocatable, dimension(:,:,:) :: p_sf
    $:GPU_DECLARE(create='[p_sf]')

contains

    !> Allocate scratch storage for the acoustic substep kernel.
    impure subroutine s_initialize_acoustic_substep_module

        @:ALLOCATE(p_sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

    end subroutine s_initialize_acoustic_substep_module

    !> Advance the conserved field in place by `n_micro` forward-backward acoustic microsteps.
    !! @param q_cons_vf    Conserved field (advanced in place).
    !! @param rhs_slow_vf  Frozen slow (advective) forcing, added dtau-scaled each microstep.
    !! @param bc_type      Boundary-condition type fields for the halo exchange.
    !! @param dtau         Microstep size (= dt_stage / n_substeps_stage).
    !! @param n_micro      Number of microsteps to take.
    impure subroutine s_acoustic_substep(q_cons_vf, rhs_slow_vf, bc_type, dtau, n_micro)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: rhs_slow_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp), intent(in) :: dtau
        integer, intent(in) :: n_micro
        real(wp) :: rho, gamma_K, pi_inf_K, qv_K                 !< Mixture EOS coefficients (frozen alpha_k)
        real(wp), dimension(num_fluids) :: alpha_rho_K, alpha_K  !< Per-fluid partial densities / volume fractions
        real(wp), dimension(2) :: Re_K                           !< Reynolds numbers (unused here; required by mixture routine)
        real(wp) :: dyn_p                                        !< Dynamic pressure 0.5*rho*|u|^2
        real(wp) :: T                                            !< Temperature (unused for the 5-eq stiffened-gas EOS branch)
        real(wp), dimension(1:num_species) :: rhoYks             !< Species partial densities (unused for non-reacting)
        real(wp) :: u_x, u_y, u_z                                !< Cell-center velocity components
        real(wp) :: div                                          !< Centered flux divergence accumulator
        real(wp) :: dpdx                                         !< Centered pressure-gradient component
        real(wp) :: damp                                         !< Divergence-damping (grad of div(rho*u)) accumulator
        real(wp) :: inv_2dx, inv_2dy, inv_2dz                    !< 1/(2*dx) etc. for centered first differences
        integer :: q                                             !< Microstep iterator
        integer :: i, j, k, l                                    !< Generic / spatial iterators

        T = 0._wp
        rhoYks = 0._wp
        Re_K = 0._wp

        do q = 1, n_micro
            ! 1. Narrow halo exchange of the conserved field (alpha_k*rho_k, rho*u, rho*E
            !    are all carried in q_cons_vf; the centered stencil needs one ghost cell).
            call s_populate_variables_buffers(bc_type, q_cons_vf)

            ! 2a. Cell pressure from the CURRENT conserved state (used by the energy flux).
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, alpha_rho_K, alpha_K, Re_K, rho, gamma_K, pi_inf_K, qv_K, &
                                & dyn_p, T, rhoYks]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        call s_compute_species_fraction(q_cons_vf, j, k, l, alpha_rho_K, alpha_K)
                        call s_convert_species_to_mixture_variables_acc(rho, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)
                        rho = max(rho, sgm_eps)

                        dyn_p = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/rho
                        end do

                        T = 0._wp
                        call s_compute_pressure(q_cons_vf(eqn_idx%E)%sf(j, k, l), q_cons_vf(eqn_idx%alf)%sf(j, k, l), dyn_p, &
                                                & pi_inf_K, gamma_K, rho, qv_K, rhoYks, p_sf(j, k, l), T)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! p_sf is filled over the full buffered range (idwbuff) from the already-exchanged
            ! q_cons_vf, so the centered energy flux at interior edges sees a consistent p in
            ! the one ghost cell it touches. No separate pressure halo exchange is needed.

            ! 2b. FORWARD: advance alpha_k*rho_k and rho*E by the centered divergence of
            !     (alpha_k*rho_k * u) and ((rho*E + p) * u), plus dtau*rhs_slow. Velocity u
            !     is held frozen (computed from the current momentum) for this forward sweep.
            inv_2dx = 0._wp; inv_2dy = 0._wp; inv_2dz = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, u_x, u_y, u_z, div, inv_2dx, inv_2dy, inv_2dz, rho]')
            do l = idwint(3)%beg, idwint(3)%end
                do k = idwint(2)%beg, idwint(2)%end
                    do j = idwint(1)%beg, idwint(1)%end
                        inv_2dx = 1._wp/(2._wp*dx(j))

                        ! Mass equations: d/dt(alpha_k*rho_k) = -div(alpha_k*rho_k * u) + slow
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%cont%beg, eqn_idx%cont%end
                            div = inv_2dx*(f_velx(q_cons_vf, j + 1, k, l)*q_cons_vf(i)%sf(j + 1, k, l) - f_velx(q_cons_vf, j - 1, &
                                           & k, l)*q_cons_vf(i)%sf(j - 1, k, l))
                            if (n > 0) then
                                inv_2dy = 1._wp/(2._wp*dy(k))
                                div = div + inv_2dy*(f_vely(q_cons_vf, j, k + 1, l)*q_cons_vf(i)%sf(j, k + 1, &
                                                     & l) - f_vely(q_cons_vf, j, k - 1, l)*q_cons_vf(i)%sf(j, k - 1, l))
                            end if
                            if (p > 0) then
                                inv_2dz = 1._wp/(2._wp*dz(l))
                                div = div + inv_2dz*(f_velz(q_cons_vf, j, k, l + 1)*q_cons_vf(i)%sf(j, k, &
                                                     & l + 1) - f_velz(q_cons_vf, j, k, l - 1)*q_cons_vf(i)%sf(j, k, l - 1))
                            end if
                            q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + dtau*(-div + rhs_slow_vf(i)%sf(j, k, l))
                        end do

                        ! Total-energy equation: d/dt(rho*E) = -div((rho*E + p) * u) + slow
                        div = inv_2dx*((q_cons_vf(eqn_idx%E)%sf(j + 1, k, l) + p_sf(j + 1, k, l))*f_velx(q_cons_vf, j + 1, k, &
                                       & l) - (q_cons_vf(eqn_idx%E)%sf(j - 1, k, l) + p_sf(j - 1, k, l))*f_velx(q_cons_vf, j - 1, &
                                       & k, l))
                        if (n > 0) then
                            inv_2dy = 1._wp/(2._wp*dy(k))
                            div = div + inv_2dy*((q_cons_vf(eqn_idx%E)%sf(j, k + 1, l) + p_sf(j, k + 1, l))*f_vely(q_cons_vf, j, &
                                                 & k + 1, l) - (q_cons_vf(eqn_idx%E)%sf(j, k - 1, l) + p_sf(j, k - 1, &
                                                 & l))*f_vely(q_cons_vf, j, k - 1, l))
                        end if
                        if (p > 0) then
                            inv_2dz = 1._wp/(2._wp*dz(l))
                            div = div + inv_2dz*((q_cons_vf(eqn_idx%E)%sf(j, k, l + 1) + p_sf(j, k, l + 1))*f_velz(q_cons_vf, j, &
                                                 & k, l + 1) - (q_cons_vf(eqn_idx%E)%sf(j, k, l - 1) + p_sf(j, k, &
                                                 & l - 1))*f_velz(q_cons_vf, j, k, l - 1))
                        end if
                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, &
                                  & l) + dtau*(-div + rhs_slow_vf(eqn_idx%E)%sf(j, k, l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! 3. Halo exchange + EOS on the UPDATED state, so the backward pressure gradient
            !    sees the freshly advanced mass and energy.
            call s_populate_variables_buffers(bc_type, q_cons_vf)

            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, alpha_rho_K, alpha_K, Re_K, rho, gamma_K, pi_inf_K, qv_K, &
                                & dyn_p, T, rhoYks]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        call s_compute_species_fraction(q_cons_vf, j, k, l, alpha_rho_K, alpha_K)
                        call s_convert_species_to_mixture_variables_acc(rho, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)
                        rho = max(rho, sgm_eps)

                        dyn_p = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/rho
                        end do

                        T = 0._wp
                        call s_compute_pressure(q_cons_vf(eqn_idx%E)%sf(j, k, l), q_cons_vf(eqn_idx%alf)%sf(j, k, l), dyn_p, &
                                                & pi_inf_K, gamma_K, rho, qv_K, rhoYks, p_sf(j, k, l), T)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! 4. BACKWARD: advance rho*u by  -dtau*grad(p_new) + dtau*rhs_slow_mom
            !    + dtau*acoustic_div_damp*grad(div(rho*u))   (divergence damping).
            !    The momentum needs its halo current (div(rho*u) uses neighbours) - q_cons_vf
            !    was just exchanged above and momentum was untouched by the forward sweep.
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, dpdx, damp, inv_2dx, inv_2dy, inv_2dz]')
            do l = idwint(3)%beg, idwint(3)%end
                do k = idwint(2)%beg, idwint(2)%end
                    do j = idwint(1)%beg, idwint(1)%end
                        inv_2dx = 1._wp/(2._wp*dx(j))

                        ! x-momentum
                        dpdx = inv_2dx*(p_sf(j + 1, k, l) - p_sf(j - 1, k, l))
                        damp = (q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, l) - 2._wp*q_cons_vf(eqn_idx%mom%beg)%sf(j, k, &
                                & l) + q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, l))/(dx(j)*dx(j))
                        q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%beg)%sf(j, k, &
                                  & l) + dtau*(-dpdx + rhs_slow_vf(eqn_idx%mom%beg)%sf(j, k, l) + acoustic_div_damp*damp)

                        ! y-momentum
                        if (n > 0) then
                            inv_2dy = 1._wp/(2._wp*dy(k))
                            dpdx = inv_2dy*(p_sf(j, k + 1, l) - p_sf(j, k - 1, l))
                            damp = (q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, l) - 2._wp*q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, &
                                    & k, l) + q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, l))/(dy(k)*dy(k))
                            q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                      & l) + dtau*(-dpdx + rhs_slow_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) + acoustic_div_damp*damp)
                        end if

                        ! z-momentum
                        if (p > 0) then
                            inv_2dz = 1._wp/(2._wp*dz(l))
                            dpdx = inv_2dz*(p_sf(j, k, l + 1) - p_sf(j, k, l - 1))
                            damp = (q_cons_vf(eqn_idx%mom%end)%sf(j, k, l + 1) - 2._wp*q_cons_vf(eqn_idx%mom%end)%sf(j, k, &
                                    & l) + q_cons_vf(eqn_idx%mom%end)%sf(j, k, l - 1))/(dz(l)*dz(l))
                            q_cons_vf(eqn_idx%mom%end)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%end)%sf(j, k, &
                                      & l) + dtau*(-dpdx + rhs_slow_vf(eqn_idx%mom%end)%sf(j, k, l) + acoustic_div_damp*damp)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_acoustic_substep

    !> Cell-center x-velocity = rho*u_x / rho, with rho the mixture density.
    pure function f_velx(q_cons_vf, j, k, l) result(u)

        $:GPU_ROUTINE(function_name='f_velx', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        real(wp)                                            :: u, rho
        integer                                             :: i
        rho = 0._wp
        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            rho = rho + q_cons_vf(i)%sf(j, k, l)
        end do
        u = q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l)/max(rho, sgm_eps)

    end function f_velx

    !> Cell-center y-velocity.
    pure function f_vely(q_cons_vf, j, k, l) result(u)

        $:GPU_ROUTINE(function_name='f_vely', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        real(wp)                                            :: u, rho
        integer                                             :: i
        rho = 0._wp
        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            rho = rho + q_cons_vf(i)%sf(j, k, l)
        end do
        u = q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/max(rho, sgm_eps)

    end function f_vely

    !> Cell-center z-velocity.
    pure function f_velz(q_cons_vf, j, k, l) result(u)

        $:GPU_ROUTINE(function_name='f_velz', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        real(wp)                                            :: u, rho
        integer                                             :: i
        rho = 0._wp
        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            rho = rho + q_cons_vf(i)%sf(j, k, l)
        end do
        u = q_cons_vf(eqn_idx%mom%end)%sf(j, k, l)/max(rho, sgm_eps)

    end function f_velz

    !> Deallocate scratch storage.
    impure subroutine s_finalize_acoustic_substep_module

        @:DEALLOCATE(p_sf)

    end subroutine s_finalize_acoustic_substep_module

end module m_acoustic_substep
