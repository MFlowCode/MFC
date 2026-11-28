!>
!! @file
!! @brief Contains module m_ib_surface_chem

#:include 'macros.fpp'

!> @brief Reacting-surface state of an immersed boundary: the Newton solve for the surface composition and temperature, and the
!! bounded ghost-state blend. Takes no ghost-point arguments; m_ibm couples it to the ghost points.
module m_ib_surface_chem

    use m_precision_select, only: wp
    use m_constants, only: T_surface_min, T_surface_max
    use m_thermochem, only: num_species, gas_constant, get_mixture_molecular_weight, get_mixture_thermal_conductivity_mixavg, &
        & get_species_mass_diffusivities_mixavg
    use m_surface_thermochem, only: get_surface_net_production_rates, get_surface_reaction_heat_flux

    implicit none

    private; public :: s_solve_surface, s_blend_ghost_temperature, s_blend_ghost_species

contains

    !> Species flux residual at a heterogeneous reacting surface.
    subroutine s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_s, d
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        real(wp), intent(out) :: R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: mw_IP, mw_s, rho_s, sum_BG
        real(wp)              :: Xs_IP(num_species), Xs_s(num_species)
        real(wp)              :: D_s(num_species), B_s(num_species), G_s(num_species)
        integer               :: k

        call get_mixture_molecular_weight(Ys_IP, mw_IP)
        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        do k = 1, num_species
            Xs_IP(k) = Ys_IP(k)*mw_IP/W_species(k)
            Xs_s(k) = Ys_s(k)*mw_s/W_species(k)
        end do

        call get_species_mass_diffusivities_mixavg(pres, T_s, Ys_s, D_s)
        call get_surface_net_production_rates(rho_s, T_s, Ys_s, omega_s)

        mdot_s = 0._wp
        do k = 1, num_species
            mdot_s = mdot_s + W_species(k)*omega_s(k)
        end do

        sum_BG = 0._wp
        do k = 1, num_species
            B_s(k) = rho_s*D_s(k)*W_species(k)/mw_s
            G_s(k) = (Xs_IP(k) - Xs_s(k))/d
            sum_BG = sum_BG + B_s(k)*G_s(k)
        end do

        do k = 1, num_species
            R_species(k) = -B_s(k)*G_s(k) + Ys_s(k)*(sum_BG + mdot_s) - W_species(k)*omega_s(k)
        end do

    end subroutine s_surface_species_residual

    !> Surface energy residual: gas-side conduction balances heterogeneous reaction heat. Radiation and solid-side conduction are
    !! omitted.
    subroutine s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d
        real(wp), intent(in)  :: Ys_s(num_species)
        real(wp), intent(out) :: R_energy
        real(wp)              :: mw_s, rho_s, k_s, q_rxn

        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        call get_mixture_thermal_conductivity_mixavg(T_s, Ys_s, k_s)
        call get_surface_reaction_heat_flux(rho_s, T_s, Ys_s, q_rxn)

        R_energy = k_s*(T_s - T_IP)/d - q_rxn

    end subroutine s_surface_energy_residual

    !> Assemble the Newton residual for Ns species, with temperature appended only when it is solved. k_bath selects which species
    !! balance the sum constraint displaces; see s_solve_surface for why it must be the most abundant one.
    subroutine s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                  & energy_scale, R, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d, flux_scale, energy_scale
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        integer, intent(in)   :: k_bath
        logical, intent(in)   :: solve_temperature
        real(wp), intent(out) :: R(num_species + 1), R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: R_energy
        integer               :: k

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        R = 0._wp
        do k = 1, num_species
            R(k) = R_species(k)/flux_scale
        end do
        R(k_bath) = sum(Ys_s) - 1._wp

        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            R(num_species + 1) = R_energy/energy_scale
        end if

    end subroutine s_surface_residual

    !> Ghost state for an immersed surface, from the image-point state and the surface state.
    !!
    !! The natural closure is the linear mirror phi_g = 2*phi_s - phi_IP, which is what makes the midpoint of the ghost/image
    !! pair reproduce the surface value. Its -1 coefficient on phi_IP is unconditional, though, so it drives strictly positive
    !! quantities negative whenever the surface value sits far below the image-point value: a cold wall in hot gas gives a
    !! negative ghost temperature and hence a negative ghost density, and a species the surface consumes faster than half the
    !! free-stream value gives a negative mass fraction. Rather than clamp the result -- which would break sum(Y) = 1 and hide
    !! the excursion -- the ghost state is the largest convex blend of the mirror back toward the surface value that stays
    !! physical:
    !!
    !!     phi_g(theta) = phi_s + theta*(phi_s - phi_IP),   theta in [0, 1]
    !!
    !! theta = 1 is the full second-order mirror; theta = 0 is the first-order Dirichlet ghost phi_g = phi_s, which is
    !! Gibou et al. (JCP 176:205, 2002) Eq. 17 and is what they likewise fall back to where the linear form is ill-behaved.
    !! Each species and the temperature get their own theta (s_blend_ghost_species keeps sum(Y) = 1 through the bath species).
    !! Sharing one would let a trace radical the surface consumes to ~1e-9 (O and OH at a burning carbon wall) drag every other
    !! mirror down with it, turning a 2100 K ghost into 1205 K over a species whose own excursion was 1e-9. The cost of limiting
    !! is small: the mirror already delivers only first-order wall flux (Ezra et al., Int. J. Heat Mass Transfer, 2025),
    !! and a boundary closure one order below the interior scheme retains the interior convergence rate (Gustafsson,
    !! Math. Comp. 29:396, 1975).
    !!
    !! This bounds the ghost state; it does not make the boundary strictly conservative. Baskaya et al. (Computers & Fluids
    !! 270:106134, 2024) trace the same cold-wall, large-gradient regime in ablation to ghost-cell mass conservation error
    !! that surfaces as spurious blowing, and resolve it only by moving to a flux-based cut-cell boundary -- a different
    !! discretization from this one, not a tuning of it.
    subroutine s_blend_ghost_temperature(T_IP, T_s, T_lo, T_hi, theta_max, T_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: T_IP, T_s, T_lo, T_hi
        real(wp), intent(in)  :: theta_max  !< Largest mirror weight the caller allows (f_face_theta), at most 1
        real(wp), intent(out) :: T_g
        ! Stop short of the bound so roundoff cannot cross it.
        real(wp), parameter :: blend_safety = 0.9_wp
        real(wp)            :: theta_T

        ! Largest theta <= theta_max keeping T_g in [T_lo, T_hi]. The bounds are compared before dividing, so +-huge() cannot
        ! overflow.
        theta_T = theta_max
        if (T_s < T_IP .and. blend_safety*(T_s - T_lo) < theta_T*(T_IP - T_s)) then
            theta_T = blend_safety*(T_s - T_lo)/(T_IP - T_s)
        else if (T_s > T_IP .and. blend_safety*(T_hi - T_s) < theta_T*(T_s - T_IP)) then
            theta_T = blend_safety*(T_hi - T_s)/(T_s - T_IP)
        end if
        theta_T = max(theta_T, 0._wp)

        T_g = T_s + theta_T*(T_s - T_IP)

    end subroutine s_blend_ghost_temperature

    !> Species half of the ghost-state blend: each species takes its own positivity-limited theta and the bath species (the most
    !! abundant at the surface) absorbs the difference, so sum(Y_g) = sum(Y_s). One theta shared by every species collapsed to
    !! ~0.006 at a burning carbon wall, where the surface consumes the trace radicals (O, OH), and made every species ghost first
    !! order. The shared theta is the fallback if the bath species would go negative.
    subroutine s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species)
        real(wp), intent(out) :: Ys_g(num_species)
        real(wp), parameter   :: blend_safety = 0.9_wp
        real(wp)              :: theta_Y
        integer               :: k, k_bath

        call s_pick_bath_species(Ys_s, k_bath)
        Ys_g(k_bath) = Ys_s(k_bath)
        do k = 1, num_species
            if (k == k_bath) cycle
            theta_Y = 1._wp
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(1._wp, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
            Ys_g(k_bath) = Ys_g(k_bath) + Ys_s(k) - Ys_g(k)
        end do
        if (Ys_g(k_bath) >= 0._wp) return

        theta_Y = 1._wp
        do k = 1, num_species
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(theta_Y, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
        end do
        do k = 1, num_species
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
        end do

    end subroutine s_blend_ghost_species

    !> Index of the most abundant species, i.e. the balance the sum constraint displaces.
    !!
    !! Closing sum(Y) = 1 by dropping one species' flux balance is standard (Surface CHEMKIN), but the dropped species absorbs
    !! the roundoff of every other balance, so it must be the bath gas. Dropping a fixed index instead -- the last species in
    !! the mechanism -- lands that error on whatever the mechanism happens to list last, which for the reduced GRI mechanism
    !! shipped with the reacting-surface example is H2O2, a trace radical whose own surface balance then goes unenforced.
    !! Cantera's solveSP re-scans for the largest species each iteration (evalSurfLarge) for this reason.
    subroutine s_pick_bath_species(Ys_s, k_bath)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: Ys_s(num_species)
        integer, intent(out) :: k_bath
        integer              :: k

        k_bath = 1
        do k = 2, num_species
            if (Ys_s(k) > Ys_s(k_bath)) k_bath = k
        end do

    end subroutine s_pick_bath_species

    !> Newton solve for a reacting surface. thermal_bc=0: zero-normal-gradient T; 1: prescribed T; 2: energy balance.
    subroutine s_solve_surface(pres, T_IP, T_wall, d, Ys_IP, W_species, thermal_bc, Ys_s, T_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_wall, d
        real(wp), intent(in)  :: Ys_IP(num_species), W_species(num_species)
        integer, intent(in)   :: thermal_bc
        real(wp), intent(out) :: Ys_s(num_species), T_s, mdot_s
        logical, intent(out)  :: converged
        integer, parameter    :: max_iter = 30, max_backtrack = 20
        real(wp), parameter   :: fd_eps_Y = 1.e-7_wp, fd_eps_T = 1.e-6_wp
        ! The Jacobian is a forward difference, so the Newton update carries O(sqrt(epsilon)) noise and the iteration turns
        ! linear near the root. The residual itself is evaluated exactly, so a tighter tolerance is not unreachable, but it is
        ! reached by stalling rather than converging -- which the caller then reports as a failed surface. Cantera's solveSP
        ! runs the same 1e-7 forward-difference step against a 1e-4 relative tolerance; 1e-6 here stays well inside that.
        real(wp), parameter :: tol = 1.e-6_wp
        real(wp), parameter :: Y_tol = 100._wp*epsilon(1._wp)
        real(wp)            :: A(num_species + 1, num_species + 1), rhs(num_species + 1), delta(num_species + 1)
        real(wp)            :: R(num_species + 1), R_pert(num_species + 1), R_trial(num_species + 1)
        real(wp)            :: R_species(num_species), R_species_pert(num_species), R_species_trial(num_species)
        real(wp)            :: omega_s(num_species), omega_pert(num_species), omega_trial(num_species)
        real(wp)            :: Ys_pert(num_species), Ys_trial(num_species)
        real(wp)            :: mdot_pert, mdot_trial, T_pert, T_trial
        real(wp)            :: flux_scale, energy_scale, R_energy, dx, lambda, norm_R, norm_trial
        logical             :: solve_temperature, linear_success, accepted
        integer             :: nsolve, iter, j, iback, k_bath, k_bath_prev

        converged = .false.
        Ys_s = Ys_IP

        select case (thermal_bc)
        case (0)
            T_s = T_IP
            solve_temperature = .false.
        case (1)
            T_s = T_wall
            solve_temperature = .false.
        case (2)
            T_s = min(max(T_IP, T_surface_min), T_surface_max)
            solve_temperature = .true.
        case default
            T_s = T_IP
            omega_s = 0._wp
            mdot_s = 0._wp
            return
        end select

        nsolve = num_species + merge(1, 0, solve_temperature)

        call s_pick_bath_species(Ys_s, k_bath)

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)
        flux_scale = max(maxval(abs(R_species)), 1.e-12_wp)
        energy_scale = 1._wp
        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            energy_scale = max(abs(R_energy), 1._wp)
        end if

        call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, energy_scale, &
                                & R, R_species, omega_s, mdot_s)
        norm_R = maxval(abs(R(1:nsolve)))
        if (norm_R < tol) then
            converged = .true.
            return
        end if

        do iter = 1, max_iter
            do j = 1, num_species
                Ys_pert = Ys_s
                T_pert = T_s
                dx = fd_eps_Y*max(abs(Ys_s(j)), 1._wp)
                if (Ys_s(j) + dx > 1._wp) dx = -dx
                Ys_pert(j) = Ys_pert(j) + dx

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,j) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end do

            if (solve_temperature) then
                Ys_pert = Ys_s
                dx = fd_eps_T*max(abs(T_s), 1._wp)
                T_pert = T_s + dx
                if (T_pert > T_surface_max) then
                    dx = -dx
                    T_pert = T_s + dx
                end if

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,nsolve) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end if

            rhs(1:nsolve) = -R(1:nsolve)
            call s_solve_surface_linear_system(A, rhs, delta, nsolve, linear_success)
            if (.not. linear_success) return

            lambda = 1._wp
            accepted = .false.
            do iback = 1, max_backtrack
                Ys_trial = Ys_s + lambda*delta(1:num_species)
                T_trial = T_s
                if (solve_temperature) T_trial = T_s + lambda*delta(nsolve)

                if (minval(Ys_trial) < -Y_tol .or. maxval(Ys_trial) > 1._wp + Y_tol .or. T_trial < T_surface_min &
                    & .or. T_trial > T_surface_max) then
                    lambda = 0.5_wp*lambda
                    cycle
                end if

                where (Ys_trial < 0._wp) Ys_trial = 0._wp
                where (Ys_trial > 1._wp) Ys_trial = 1._wp

                call s_surface_residual(pres, T_IP, T_trial, d, Ys_IP, Ys_trial, W_species, k_bath, solve_temperature, &
                                        & flux_scale, energy_scale, R_trial, R_species_trial, omega_trial, mdot_trial)
                norm_trial = maxval(abs(R_trial(1:nsolve)))
                if (norm_trial < norm_R) then
                    accepted = .true.
                    exit
                end if
                lambda = 0.5_wp*lambda
            end do

            if (.not. accepted) return

            Ys_s = Ys_trial
            T_s = T_trial
            R = R_trial
            R_species = R_species_trial
            omega_s = omega_trial
            mdot_s = mdot_trial
            norm_R = norm_trial

            ! The composition has moved, so the species the sum constraint displaces may no longer be the most abundant one.
            ! Re-closing changes what R means, so the residual is rebuilt on the rare steps where the choice actually changes.
            k_bath_prev = k_bath
            call s_pick_bath_species(Ys_s, k_bath)
            if (k_bath /= k_bath_prev) then
                call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R, R_species, omega_s, mdot_s)
                norm_R = maxval(abs(R(1:nsolve)))
            end if

            if (norm_R < tol) then
                converged = .true.
                return
            end if
        end do

    end subroutine s_solve_surface

    !> Small dense linear solve with partial pivoting for the local surface Newton system.
    subroutine s_solve_surface_linear_system(A, b, x, nsolve, success)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)     :: nsolve
        real(wp), intent(inout) :: A(num_species + 1, num_species + 1)
        real(wp), intent(inout) :: b(num_species + 1)
        real(wp), intent(out)   :: x(num_species + 1)
        logical, intent(out)    :: success
        real(wp)                :: factor, pivot_value, tmp, row_tmp(num_species + 1)
        integer                 :: i, j, k, pivot

        success = .true.
        x = 0._wp

        do k = 1, nsolve - 1
            pivot = k
            pivot_value = abs(A(k, k))
            do i = k + 1, nsolve
                if (abs(A(i, k)) > pivot_value) then
                    pivot = i
                    pivot_value = abs(A(i, k))
                end if
            end do
            if (pivot_value <= epsilon(1._wp)) then
                success = .false.
                return
            end if

            if (pivot /= k) then
                row_tmp(1:nsolve) = A(k,1:nsolve)
                A(k,1:nsolve) = A(pivot,1:nsolve)
                A(pivot,1:nsolve) = row_tmp(1:nsolve)
                tmp = b(k)
                b(k) = b(pivot)
                b(pivot) = tmp
            end if

            do i = k + 1, nsolve
                factor = A(i, k)/A(k, k)
                A(i, k) = 0._wp
                do j = k + 1, nsolve
                    A(i, j) = A(i, j) - factor*A(k, j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(A(nsolve, nsolve)) <= epsilon(1._wp)) then
            success = .false.
            return
        end if

        x(nsolve) = b(nsolve)/A(nsolve, nsolve)
        do i = nsolve - 1, 1, -1
            tmp = b(i)
            do j = i + 1, nsolve
                tmp = tmp - A(i, j)*x(j)
            end do
            if (abs(A(i, i)) <= epsilon(1._wp)) then
                success = .false.
                return
            end if
            x(i) = tmp/A(i, i)
        end do

    end subroutine s_solve_surface_linear_system

end module m_ib_surface_chem
