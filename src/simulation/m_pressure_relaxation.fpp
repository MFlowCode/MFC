!>
!! @file m_pressure_relaxation.fpp
!! @brief Contains module m_pressure_relaxation

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The module contains the subroutines used to perform pressure relaxation
!!        for multi-component flows using the 6-equation model. This includes
!!        volume fraction correction, Newton-Raphson pressure equilibration, and
!!        internal energy correction to maintain thermodynamic consistency.
module m_pressure_relaxation

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters

    implicit none

    private; public :: s_pressure_relaxation_procedure, &
 s_initialize_pressure_relaxation_module, &
 s_finalize_pressure_relaxation_module

    real(wp), allocatable, dimension(:) :: gamma_min, pres_inf
    $:GPU_DECLARE(create='[gamma_min, pres_inf]')

    real(wp), allocatable, dimension(:, :) :: Res_pr
    $:GPU_DECLARE(create='[Res_pr]')

contains

    !> Initialize the pressure relaxation module
    impure subroutine s_initialize_pressure_relaxation_module

        integer :: i, j

        @:ALLOCATE(gamma_min(1:num_fluids), pres_inf(1:num_fluids))

        do i = 1, num_fluids
            gamma_min(i) = 1._wp/fluid_pp(i)%gamma + 1._wp
            pres_inf(i) = fluid_pp(i)%pi_inf/(1._wp + fluid_pp(i)%gamma)
        end do
        $:GPU_UPDATE(device='[gamma_min, pres_inf]')

        if (viscous) then
            @:ALLOCATE(Res_pr(1:2, 1:Re_size_max))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_pr(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res_pr, Re_idx, Re_size]')
        end if

    end subroutine s_initialize_pressure_relaxation_module

    !> Finalize the pressure relaxation module
    impure subroutine s_finalize_pressure_relaxation_module

        @:DEALLOCATE(gamma_min, pres_inf)
        if (viscous) then
            @:DEALLOCATE(Res_pr)
        end if

    end subroutine s_finalize_pressure_relaxation_module

    !> The main pressure relaxation procedure
    !! @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_procedure(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer :: j, k, l

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_relax_cell_pressure(q_cons_vf, j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pressure_relaxation_procedure

    !> Process pressure relaxation for a single cell
    subroutine s_relax_cell_pressure(q_cons_vf, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in) :: j, k, l

        ! Volume fraction correction
        if (mpp_lim) call s_correct_volume_fractions(q_cons_vf, j, k, l)

        ! Pressure equilibration
        if (s_needs_pressure_relaxation(q_cons_vf, j, k, l)) then
            call s_equilibrate_pressure(q_cons_vf, j, k, l)
        end if

        ! Internal energy correction
        call s_correct_internal_energies(q_cons_vf, j, k, l)

    end subroutine s_relax_cell_pressure

    !> Check if pressure relaxation is needed for this cell
    logical function s_needs_pressure_relaxation(q_cons_vf, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in) :: j, k, l
        integer :: i

        s_needs_pressure_relaxation = .true.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > (1._wp - sgm_eps)) then
                s_needs_pressure_relaxation = .false.
            end if
        end do

    end function s_needs_pressure_relaxation

    !> Correct volume fractions to physical bounds
    subroutine s_correct_volume_fractions(q_cons_vf, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in) :: j, k, l
        real(wp) :: sum_alpha
        integer :: i

        sum_alpha = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0._wp) .or. &
                (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0._wp)) then
                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0._wp
            end if
            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > 1._wp) &
                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1._wp
            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
        end do

    end subroutine s_correct_volume_fractions

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in) :: j, k, l

        real(wp) :: pres_relax, f_pres, df_pres
        real(wp), dimension(num_fluids) :: pres_K_init, rho_K_s
        integer, parameter :: MAX_ITER = 50
        real(wp), parameter :: TOLERANCE = 1.e-10_wp
        integer :: iter, i

        ! Initialize pressures
        pres_relax = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                pres_K_init(i) = (q_cons_vf(i + intxb - 1)%sf(j, k, l)/ &
                                  q_cons_vf(i + advxb - 1)%sf(j, k, l) - pi_infs(i))/gammas(i)
                if (pres_K_init(i) <= -(1._wp - 1.e-8_wp)*pres_inf(i) + 1.e-8_wp) &
                    pres_K_init(i) = -(1._wp - 1.e-8_wp)*pres_inf(i) + 1.e-8_wp
            else
                pres_K_init(i) = 0._wp
            end if
            pres_relax = pres_relax + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pres_K_init(i)
        end do

        ! Newton-Raphson iteration
        f_pres = 1.e-9_wp
        df_pres = 1.e9_wp
        $:GPU_LOOP(parallelism='[seq]')
        do iter = 0, MAX_ITER - 1
            if (abs(f_pres) > TOLERANCE) then
                pres_relax = pres_relax - f_pres/df_pres

                ! Enforce pressure bounds
                do i = 1, num_fluids
                    if (pres_relax <= -(1._wp - 1.e-8_wp)*pres_inf(i) + 1.e-8_wp) &
                        pres_relax = -(1._wp - 1.e-8_wp)*pres_inf(i) + 1._wp
                end do

                ! Newton-Raphson step
                f_pres = -1._wp
                df_pres = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                        rho_K_s(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/ &
                                     max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                     *((pres_relax + pres_inf(i))/(pres_K_init(i) + &
                                                                   pres_inf(i)))**(1._wp/gamma_min(i))
                        f_pres = f_pres + q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho_K_s(i)
                        df_pres = df_pres - q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                  /(gamma_min(i)*rho_K_s(i)*(pres_relax + pres_inf(i)))
                    end if
                end do
            end if
        end do

        ! Update volume fractions
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) &
                q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho_K_s(i)
        end do

    end subroutine s_equilibrate_pressure

    !> Correct internal energies using equilibrated pressure
    subroutine s_correct_internal_energies(q_cons_vf, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in) :: j, k, l

        real(wp), dimension(num_fluids) :: alpha_rho, alpha
        real(wp) :: rho, dyn_pres, gamma, pi_inf, pres_relax, sum_alpha
        real(wp), dimension(2) :: Re
        integer :: i, q

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
            alpha(i) = q_cons_vf(E_idx + i)%sf(j, k, l)
        end do

        ! Compute mixture properties (combined bubble and standard logic)
        rho = 0._wp
        gamma = 0._wp
        pi_inf = 0._wp

        if (bubbles_euler) then
            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    rho = rho + alpha_rho(i)
                    gamma = gamma + alpha(i)*gammas(i)
                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                end do
            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids - 1
                    rho = rho + alpha_rho(i)
                    gamma = gamma + alpha(i)*gammas(i)
                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                end do
            else
                rho = alpha_rho(1)
                gamma = gammas(1)
                pi_inf = pi_infs(1)
            end if
        else
            sum_alpha = 0._wp
            if (mpp_lim) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    alpha_rho(i) = max(0._wp, alpha_rho(i))
                    alpha(i) = min(max(0._wp, alpha(i)), 1._wp)
                    sum_alpha = sum_alpha + alpha(i)
                end do
                alpha = alpha/max(sum_alpha, sgm_eps)
            end if

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                rho = rho + alpha_rho(i)
                gamma = gamma + alpha(i)*gammas(i)
                pi_inf = pi_inf + alpha(i)*pi_infs(i)
            end do

            if (viscous) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, 2
                    Re(i) = dflt_real
                    if (Re_size(i) > 0) Re(i) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, Re_size(i)
                        Re(i) = alpha(Re_idx(i, q))/Res_pr(i, q) + Re(i)
                    end do
                    Re(i) = 1._wp/max(Re(i), sgm_eps)
                end do
            end if
        end if

        ! Compute dynamic pressure and update internal energies
        dyn_pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = momxb, momxe
            dyn_pres = dyn_pres + 5.e-1_wp*q_cons_vf(i)%sf(j, k, l)* &
                       q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do

        pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - pi_inf)/gamma

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                q_cons_vf(i + advxb - 1)%sf(j, k, l)*(gammas(i)*pres_relax + pi_infs(i))
        end do

    end subroutine s_correct_internal_energies

end module m_pressure_relaxation
