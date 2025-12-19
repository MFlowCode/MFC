!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive ones and vice versa. In
!!              addition, the module also contains the subroutines used to obtain
!!              the mixture variables and the subroutines used to compute pressure.
module m_variables_conversion

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    use m_thermochem, only: &
        num_species, get_temperature, get_pressure, gas_constant, &
        get_mixture_molecular_weight, get_mixture_energy_mass

    implicit none

    private; 
    public :: s_initialize_variables_conversion_module, &
              s_initialize_pb, &
              s_initialize_mv, &
              s_convert_to_mixture_variables, &
              s_convert_mixture_to_mixture_variables, &
              s_convert_species_to_mixture_variables, &
              s_convert_species_to_mixture_variables_acc, &
              s_convert_conservative_to_primitive_variables, &
              s_convert_primitive_to_conservative_variables, &
              s_convert_primitive_to_flux_variables, &
              s_compute_pressure, &
              s_compute_species_fraction, &
#ifndef MFC_PRE_PROCESS
              s_compute_speed_of_sound, &
              s_compute_fast_magnetosonic_speed, &
#endif
              s_finalize_variables_conversion_module

    !! In simulation, gammas, pi_infs, and qvs are already declared in m_global_variables
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    $:GPU_DECLARE(create='[gammas,gs_min,pi_infs,ps_inf,cvs,qvs,qvps]')
#endif

    real(wp), allocatable, dimension(:) :: Gs_vc
    integer, allocatable, dimension(:) :: bubrs_vc
    real(wp), allocatable, dimension(:, :) :: Res_vc
    $:GPU_DECLARE(create='[bubrs_vc,Gs_vc,Res_vc]')

    integer :: is1b, is2b, is3b, is1e, is2e, is3e
    $:GPU_DECLARE(create='[is1b,is2b,is3b,is1e,is2e,is3e]')

    real(wp), allocatable, dimension(:, :, :), public :: rho_sf !< Scalar density function
    real(wp), allocatable, dimension(:, :, :), public :: gamma_sf !< Scalar sp. heat ratio function
    real(wp), allocatable, dimension(:, :, :), public :: pi_inf_sf !< Scalar liquid stiffness function
    real(wp), allocatable, dimension(:, :, :), public :: qv_sf !< Scalar liquid energy reference function

contains

    !> Dispatch to the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines.
        !!      Replaces a procedure pointer.
        !!  @param q_vf Conservative or primitive variables
        !!  @param i First-coordinate cell index
        !!  @param j First-coordinate cell index
        !!  @param k First-coordinate cell index
        !!  @param rho Density
        !!  @param gamma Specific heat ratio function
        !!  @param pi_inf Liquid stiffness function
        !!  @param qv Fluid reference energy
    subroutine s_convert_to_mixture_variables(q_vf, i, j, k, &
                                              rho, gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in) :: i, j, k
        real(wp), intent(out), target :: rho, gamma, pi_inf, qv
        real(wp), optional, dimension(2), intent(out) :: Re_K
        real(wp), optional, intent(out) :: G_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            call s_convert_mixture_to_mixture_variables(q_vf, i, j, k, &
                                                        rho, gamma, pi_inf, qv)

        else  ! Volume fraction model
            call s_convert_species_to_mixture_variables(q_vf, i, j, k, &
                                                        rho, gamma, pi_inf, qv, Re_K, G_K, G)
        end if

    end subroutine s_convert_to_mixture_variables

    !>  This procedure conditionally calculates the appropriate pressure
        !! @param energy Energy
        !! @param alf Void Fraction
        !! @param dyn_p Dynamic Pressure
        !! @param pi_inf Liquid Stiffness
        !! @param gamma Specific Heat Ratio
        !! @param rho Density
        !! @param qv fluid reference energy
        !! @param pres Pressure to calculate
        !! @param stress Shear Stress
        !! @param mom Momentum
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, rhoYks, pres, T, stress, mom, G, pres_mag)
        $:GPU_ROUTINE(function_name='s_compute_pressure',parallelism='[seq]', &
            & cray_inline=True)

        real(stp), intent(in) :: energy, alf
        real(wp), intent(in) :: dyn_p
        real(wp), intent(in) :: pi_inf, gamma, rho, qv
        real(wp), intent(out) :: pres
        real(wp), intent(inout) :: T
        real(stp), intent(in), optional :: stress, mom
        real(wp), intent(in), optional :: G, pres_mag

        ! Chemistry
        real(wp), dimension(1:num_species), intent(in) :: rhoYks
        real(wp) :: E_e
        real(wp) :: e_Per_Kg, Pdyn_Per_Kg
        real(wp) :: T_guess
        real(wp), dimension(1:num_species) :: Y_rs

        integer :: s !< Generic loop iterator

        #:if not chemistry
            ! Depending on model_eqns and bubbles_euler, the appropriate procedure
            ! for computing pressure is targeted by the procedure pointer

            if (mhd) then
                pres = (energy - dyn_p - pi_inf - qv - pres_mag)/gamma
            elseif ((model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                pres = (energy - dyn_p - pi_inf - qv)/gamma
            else if ((model_eqns /= 4) .and. bubbles_euler) then
                pres = ((energy - dyn_p)/(1._wp - alf) - pi_inf - qv)/gamma
            else
                pres = (pref + pi_inf)* &
                       (energy/ &
                        (rhoref*(1 - alf)) &
                        )**(1/gamma + 1) - pi_inf
            end if

            if (hypoelasticity .and. present(G)) then
                ! calculate elastic contribution to Energy
                E_e = 0._wp
                do s = stress_idx%beg, stress_idx%end
                    if (G > 0) then
                        E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        ! Double for shear stresses
                        if (any(s == shear_indices)) then
                            E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        end if
                    end if
                end do

                pres = ( &
                       energy - &
                       0.5_wp*(mom**2._wp)/rho - &
                       pi_inf - qv - E_e &
                       )/gamma

            end if

        #:else

            Y_rs(:) = rhoYks(:)/rho
            e_Per_Kg = energy/rho
            Pdyn_Per_Kg = dyn_p/rho

            T_guess = T

            call get_temperature(e_Per_Kg - Pdyn_Per_Kg, T_guess, Y_rs, .true., T)
            call get_pressure(rho, T, Y_rs, pres)

        #:endif

    end subroutine s_compute_pressure

    !>  This subroutine is designed for the gamma/pi_inf model
        !!      and provided a set of either conservative or primitive
        !!      variables, transfers the density, specific heat ratio
        !!      function and the liquid stiffness function from q_vf to
        !!      rho, gamma and pi_inf.
        !! @param q_vf conservative or primitive variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        !! @param rho density
        !! @param gamma  specific heat ratio function
        !! @param pi_inf liquid stiffness
        !! @param qv fluid reference energy
    subroutine s_convert_mixture_to_mixture_variables(q_vf, i, j, k, &
                                                      rho, gamma, pi_inf, qv)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in) :: i, j, k

        real(wp), intent(out), target :: rho
        real(wp), intent(out), target :: gamma
        real(wp), intent(out), target :: pi_inf
        real(wp), intent(out), target :: qv

        ! Transferring the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        rho = q_vf(1)%sf(i, j, k)
        gamma = q_vf(gamma_idx)%sf(i, j, k)
        pi_inf = q_vf(pi_inf_idx)%sf(i, j, k)
        qv = 0._wp ! keep this value nill for now. For future adjustment

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(i, j, k) = rho
        gamma_sf(i, j, k) = gamma
        pi_inf_sf(i, j, k) = pi_inf
        qv_sf(i, j, k) = qv
#endif

    end subroutine s_convert_mixture_to_mixture_variables

    !>  This subroutine is designed for the volume fraction model
        !!              and provided a set of either conservative or primitive
        !!              variables, computes the density, the specific heat ratio
        !!              function and the liquid stiffness function from q_vf and
        !!              stores the results into rho, gamma and pi_inf.
        !! @param q_vf primitive variables
        !! @param k Cell index
        !! @param l Cell index
        !! @param r Cell index
        !! @param rho density
        !! @param gamma specific heat ratio
        !! @param pi_inf liquid stiffness
        !! @param qv fluid reference energy
    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, r, rho, &
                                                      gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf

        integer, intent(in) :: k, l, r

        real(wp), intent(out), target :: rho
        real(wp), intent(out), target :: gamma
        real(wp), intent(out), target :: pi_inf
        real(wp), intent(out), target :: qv

        real(wp), optional, dimension(2), intent(out) :: Re_K
        real(wp), optional, intent(out) :: G_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G
        real(wp), dimension(num_fluids) :: alpha_rho_K, alpha_K !<

        integer :: i, j !< Generic loop iterator

        ! Computing the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        call s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        ! Calculating the density, the specific heat ratio function, the
        ! liquid stiffness function, and the energy reference function,
        ! respectively, from the species analogs
        if (num_fluids == 1 .and. bubbles_euler) then
            rho = alpha_rho_K(1)
            gamma = gammas(1)
            pi_inf = pi_infs(1)
            qv = qvs(1)
        else
            rho = 0._wp; gamma = 0._wp; pi_inf = 0._wp; qv = 0._wp
            do i = 1, num_fluids
                rho = rho + alpha_rho_K(i)
                gamma = gamma + alpha_K(i)*gammas(i)
                pi_inf = pi_inf + alpha_K(i)*pi_infs(i)
                qv = qv + alpha_rho_K(i)*qvs(i)
            end do
        end if

#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) &
                              + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if
#endif

        if (present(G_K)) then
            G_K = 0._wp
            do i = 1, num_fluids
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0._wp, G_K)
        end if

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(k, l, r) = rho
        gamma_sf(k, l, r) = gamma
        pi_inf_sf(k, l, r) = pi_inf
        qv_sf(k, l, r) = qv
#endif

    end subroutine s_convert_species_to_mixture_variables

    subroutine s_convert_species_to_mixture_variables_acc(rho_K, &
                                                          gamma_K, pi_inf_K, qv_K, &
                                                          alpha_K, alpha_rho_K, Re_K, &
                                                          G_K, G)
        $:GPU_ROUTINE(function_name='s_convert_species_to_mixture_variables_acc', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K
        real(wp), dimension(num_fluids), intent(inout) :: alpha_rho_K, alpha_K !<
        real(wp), dimension(2), intent(out) :: Re_K
        real(wp), optional, intent(out) :: G_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G

        integer :: i, j !< Generic loop iterators

#ifdef MFC_SIMULATION
        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        if (num_fluids == 1 .and. bubbles_euler) then
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
            qv_K = qvs(1)
        else
            if (mpp_lim) then
                do i = 1, num_fluids
                    alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                    alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                end do
                alpha_K = alpha_K/max(sum(alpha_K), sgm_eps)
            end if
            rho_K = 0._wp; gamma_K = 0._wp; pi_inf_K = 0._wp; qv_K = 0._wp
            do i = 1, num_fluids
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
                qv_K = qv_K + alpha_rho_K(i)*qvs(i)
            end do
        end if

        if (present(G_K)) then
            G_K = 0._wp
            do i = 1, num_fluids
                !TODO: change to use Gs_vc directly here?
                !TODO: Make this changes as well for GPUs
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0._wp, G_K)
        end if

        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res_vc(i, j) &
                              + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if
#endif

    end subroutine s_convert_species_to_mixture_variables_acc

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_variables_conversion_module

        integer :: i, j

        $:GPU_ENTER_DATA(copyin='[is1b,is1e,is2b,is2e,is3b,is3e]')

        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(gs_min (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(ps_inf(1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps    (1:num_fluids))
        @:ALLOCATE(Gs_vc     (1:num_fluids))

        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            gs_min(i) = 1.0_wp/gammas(i) + 1.0_wp
            pi_infs(i) = fluid_pp(i)%pi_inf
            Gs_vc(i) = fluid_pp(i)%G
            ps_inf(i) = pi_infs(i)/(1.0_wp + gammas(i))
            cvs(i) = fluid_pp(i)%cv
            qvs(i) = fluid_pp(i)%qv
            qvps(i) = fluid_pp(i)%qvp
        end do
        $:GPU_UPDATE(device='[gammas,gs_min,pi_infs,ps_inf,cvs,qvs,qvps,Gs_vc]')

#ifdef MFC_SIMULATION

        if (viscous) then
            @:ALLOCATE(Res_vc(1:2, 1:Re_size_max))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_vc(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do

            $:GPU_UPDATE(device='[Res_vc,Re_idx,Re_size]')
        end if
#endif

        if (bubbles_euler) then
            @:ALLOCATE(bubrs_vc(1:nb))
            do i = 1, nb
                bubrs_vc(i) = bub_idx%rs(i)
            end do
            $:GPU_UPDATE(device='[bubrs_vc]')
        end if

#ifdef MFC_POST_PROCESS
        ! Allocating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then

            ! Simulation is 3D
            if (p > 0) then

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 -buff_size:p + buff_size))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   -buff_size:p + buff_size))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    -buff_size:p + buff_size))
                allocate (qv_sf(-buff_size:m + buff_size, &
                                -buff_size:n + buff_size, &
                                -buff_size:p + buff_size))

                ! Simulation is 2D
            else

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 0:0))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    0:0))
                allocate (qv_sf(-buff_size:m + buff_size, &
                                -buff_size:n + buff_size, &
                                0:0))
            end if

            ! Simulation is 1D
        else

            allocate (rho_sf(-buff_size:m + buff_size, &
                             0:0, &
                             0:0))
            allocate (gamma_sf(-buff_size:m + buff_size, &
                               0:0, &
                               0:0))
            allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                0:0, &
                                0:0))
            allocate (qv_sf(-buff_size:m + buff_size, &
                            0:0, &
                            0:0))

        end if
#endif

    end subroutine s_initialize_variables_conversion_module

    !Initialize mv at the quadrature nodes based on the initialized moments and sigma
    subroutine s_initialize_mv(qK_cons_vf, mv)

        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf

        real(stp), dimension(idwint(1)%beg:, idwint(2)%beg:, idwint(3)%beg:, 1:, 1:), intent(inout) :: mv

        integer :: i, j, k, l
        real(wp) :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end

                    nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        mv(j, k, l, 1, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 2, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 3, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 4, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                    end do

                end do
            end do
        end do

    end subroutine s_initialize_mv

    !Initialize pb at the quadrature nodes using isothermal relations (Preston model)
    subroutine s_initialize_pb(qK_cons_vf, mv, pb)
        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf

        real(stp), dimension(idwint(1)%beg:, idwint(2)%beg:, idwint(3)%beg:, 1:, 1:), intent(in) :: mv
        real(stp), dimension(idwint(1)%beg:, idwint(2)%beg:, idwint(3)%beg:, 1:, 1:), intent(inout) :: pb

        integer :: i, j, k, l
        real(wp) :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end

                    nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        !PRESTON (ISOTHERMAL)
                        pb(j, k, l, 1, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 1, i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 2, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 2, i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 3, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 3, i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 4, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 4, i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_pb

    !> The following procedure handles the conversion between
        !!      the conservative variables and the primitive variables.
        !! @param qK_cons_vf Conservative variables
        !! @param qK_prim_vf Primitive variables
        !! @param gm_alphaK_vf Gradient magnitude of the volume fraction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, &
                                                             q_T_sf, &
                                                             qK_prim_vf, &
                                                             ibounds)

        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: qK_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: ibounds

        real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(wp), dimension(2) :: Re_K
        real(wp) :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K
        real(wp), dimension(nb) :: nRtmp

        real(wp) :: rhoYks(1:num_species)

        real(wp) :: vftmp, nbub_sc

        real(wp) :: G_K

        real(wp) :: pres

        integer :: i, j, k, l !< Generic loop iterators

        real(wp) :: T
        real(wp) :: pres_mag

        real(wp) :: Ga ! Lorentz factor (gamma in relativity)
        real(wp) :: B2 ! Magnetic field magnitude squared
        real(wp) :: B(3) ! Magnetic field components
        real(wp) :: m2 ! Relativistic momentum magnitude squared
        real(wp) :: S ! Dot product of the magnetic field and the relativistic momentum
        real(wp) :: W, dW ! W := rho*v*Ga**2; f = f(W) in Newton-Raphson
        real(wp) :: E, D ! Prim/Cons variables within Newton-Raphson iteration
        real(wp) :: f, dGa_dW, dp_dW, df_dW ! Functions within Newton-Raphson iteration
        integer :: iter ! Newton-Raphson iteration counter

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_K, alpha_rho_K, Re_K, nRtmp, rho_K, gamma_K, pi_inf_K,qv_K, dyn_pres_K, rhoYks, B, pres, vftmp, nbub_sc, G_K, T, pres_mag, Ga, B2, m2, S, W, dW, E, D, f, dGa_dW, dp_dW, df_dW, iter ]')
        do l = ibounds(3)%beg, ibounds(3)%end
            do k = ibounds(2)%beg, ibounds(2)%end
                do j = ibounds(1)%beg, ibounds(1)%end
                    dyn_pres_K = 0._wp

                    call s_compute_species_fraction(qK_cons_vf, j, k, l, alpha_rho_K, alpha_K)

                    if (model_eqns /= 4) then
#ifdef MFC_SIMULATION
                        ! If in simulation, use acc mixture subroutines
                        if (elasticity) then
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, &
                                                                            alpha_rho_K, Re_K, G_K, Gs_vc)
                        else
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                            alpha_K, alpha_rho_K, Re_K)
                        end if
#else
                        ! If pre-processing, use non acc mixture subroutines
                        if (elasticity) then
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, &
                                                                rho_K, gamma_K, pi_inf_K, qv_K, Re_K, G_K, fluid_pp(:)%G)
                        else
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, &
                                                                rho_K, gamma_K, pi_inf_K, qv_K)
                        end if
#endif
                    end if

                    if (relativity) then
                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = qK_cons_vf(B_idx%beg)%sf(j, k, l)
                            B(3) = qK_cons_vf(B_idx%beg + 1)%sf(j, k, l)
                        else
                            B(1) = qK_cons_vf(B_idx%beg)%sf(j, k, l)
                            B(2) = qK_cons_vf(B_idx%beg + 1)%sf(j, k, l)
                            B(3) = qK_cons_vf(B_idx%beg + 2)%sf(j, k, l)
                        end if
                        B2 = B(1)**2 + B(2)**2 + B(3)**2

                        m2 = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = momxb, momxe
                            m2 = m2 + qK_cons_vf(i)%sf(j, k, l)**2
                        end do

                        S = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            S = S + qK_cons_vf(momxb + i - 1)%sf(j, k, l)*B(i)
                        end do

                        E = qK_cons_vf(E_idx)%sf(j, k, l)

                        D = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            D = D + qK_cons_vf(i)%sf(j, k, l)
                        end do

                        ! Newton-Raphson
                        W = E + D
                        $:GPU_LOOP(parallelism='[seq]')
                        do iter = 1, relativity_cons_to_prim_max_iter
                            Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                            pres = (W - D*Ga)/((gamma_K + 1)*Ga**2) ! Thermal pressure from EOS
                            f = W - pres + (1 - 1/(2*Ga**2))*B2 - S**2/(2*W**2) - E - D

                            ! The first equation below corrects a typo in (Mignone & Bodo, 2006)
                            ! m2*W**2 → 2*m2*W**2, which would cancel with the 2* in other terms
                            ! This corrected version is not used as the second equation empirically converges faster.
                            ! First equation is kept for further investigation.
                            ! dGa_dW = -Ga**3 * ( S**2*(3*W**2+3*W*B2+B2**2) + m2*W**2 ) / (W**3 * (W+B2)**3) ! first (corrected)
                            dGa_dW = -Ga**3*(2*S**2*(3*W**2 + 3*W*B2 + B2**2) + m2*W**2)/(2*W**3*(W + B2)**3) ! second (in paper)

                            dp_dW = (Ga*(1 + D*dGa_dW) - 2*W*dGa_dW)/((gamma_K + 1)*Ga**3)
                            df_dW = 1 - dp_dW + (B2/Ga**3)*dGa_dW + S**2/W**3

                            dW = -f/df_dW
                            W = W + dW
                            if (abs(dW) < 1.e-12_wp*W) exit
                        end do

                        ! Recalculate pressure using converged W
                        Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                        qK_prim_vf(E_idx)%sf(j, k, l) = (W - D*Ga)/((gamma_K + 1)*Ga**2)

                        ! Recover the other primitive variables
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            qK_prim_vf(momxb + i - 1)%sf(j, k, l) = (qK_cons_vf(momxb + i - 1)%sf(j, k, l) + (S/W)*B(i))/(W + B2)
                        end do
                        qK_prim_vf(1)%sf(j, k, l) = D/Ga ! Hard-coded for single-component for now

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = B_idx%beg, B_idx%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do

                        cycle ! skip all the non-relativistic conversions below
                    end if

                    if (chemistry) then
                        rho_K = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = chemxb, chemxe
                            rho_K = rho_K + max(0._wp, qK_cons_vf(i)%sf(j, k, l))
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            qK_prim_vf(i)%sf(j, k, l) = rho_K
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = chemxb, chemxe
                            qK_prim_vf(i)%sf(j, k, l) = max(0._wp, qK_cons_vf(i)%sf(j, k, l)/rho_K)
                        end do
                    else
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

#ifdef MFC_SIMULATION
                    rho_K = max(rho_K, sgm_eps)
#endif

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = momxb, momxe
                        if (model_eqns /= 4) then
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                        /rho_K
                            dyn_pres_K = dyn_pres_K + 5.e-1_wp*qK_cons_vf(i)%sf(j, k, l) &
                                         *qK_prim_vf(i)%sf(j, k, l)
                        else
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                        /qK_cons_vf(1)%sf(j, k, l)
                        end if
                    end do

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            rhoYks(i) = qK_cons_vf(chemxb + i - 1)%sf(j, k, l)
                        end do

                        T = q_T_sf%sf(j, k, l)
                    end if

                    if (mhd) then
                        if (n == 0) then
                            pres_mag = 0.5_wp*(Bx0**2 + qK_cons_vf(B_idx%beg)%sf(j, k, l)**2 + qK_cons_vf(B_idx%beg + 1)%sf(j, k, l)**2)
                        else
                            pres_mag = 0.5_wp*(qK_cons_vf(B_idx%beg)%sf(j, k, l)**2 + qK_cons_vf(B_idx%beg + 1)%sf(j, k, l)**2 + qK_cons_vf(B_idx%beg + 2)%sf(j, k, l)**2)
                        end if
                    else
                        pres_mag = 0._wp
                    end if

                    call s_compute_pressure(qK_cons_vf(E_idx)%sf(j, k, l), &
                                            qK_cons_vf(alf_idx)%sf(j, k, l), &
                                            dyn_pres_K, pi_inf_K, gamma_K, rho_K, &
                                            qv_K, rhoYks, pres, T, pres_mag=pres_mag)

                    qK_prim_vf(E_idx)%sf(j, k, l) = pres

                    if (chemistry) then
                        q_T_sf%sf(j, k, l) = T
                    end if

                    if (bubbles_euler) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, nb
                            nRtmp(i) = qK_cons_vf(bubrs_vc(i))%sf(j, k, l)
                        end do

                        vftmp = qK_cons_vf(alf_idx)%sf(j, k, l)

                        if (qbmm) then
                            !Get nb (constant across all R0 bins)
                            nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                            !Convert cons to prim
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = bubxb, bubxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                            !Need to keep track of nb in the primitive variable list (converted back to true value before output)
#ifdef MFC_SIMULATION
                            qK_prim_vf(bubxb)%sf(j, k, l) = qK_cons_vf(bubxb)%sf(j, k, l)
#endif

                        else
                            if (adv_n) then
                                qK_prim_vf(n_idx)%sf(j, k, l) = qK_cons_vf(n_idx)%sf(j, k, l)
                                nbub_sc = qK_prim_vf(n_idx)%sf(j, k, l)
                            else
                                call s_comp_n_from_cons(vftmp, nRtmp, nbub_sc, weight)
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = bubxb, bubxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                        end if
                    end if

                    if (mhd) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = B_idx%beg, B_idx%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (elasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = strxb, strxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (hypoelasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = strxb, strxe
                            ! subtracting elastic contribution for pressure calculation
                            if (G_K > verysmall) then
                                if (cont_damage) G_K = G_K*max((1._wp - qK_cons_vf(damage_idx)%sf(j, k, l)), 0._wp)
                                qK_prim_vf(E_idx)%sf(j, k, l) = qK_prim_vf(E_idx)%sf(j, k, l) - &
                                                                ((qK_prim_vf(i)%sf(j, k, l)**2._wp)/(4._wp*G_K))/gamma_K
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    qK_prim_vf(E_idx)%sf(j, k, l) = qK_prim_vf(E_idx)%sf(j, k, l) - &
                                                                    ((qK_prim_vf(i)%sf(j, k, l)**2._wp)/(4._wp*G_K))/gamma_K
                                end if
                            end if
                        end do
                    end if

                    if (hyperelasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = xibeg, xiend
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (.not. igr .or. num_fluids > 1) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = advxb, advxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (surface_tension) then
                        qK_prim_vf(c_idx)%sf(j, k, l) = qK_cons_vf(c_idx)%sf(j, k, l)
                    end if

                    if (cont_damage) qK_prim_vf(damage_idx)%sf(j, k, l) = qK_cons_vf(damage_idx)%sf(j, k, l)

#ifdef MFC_POST_PROCESS
                    if (bubbles_lagrange) qK_prim_vf(beta_idx)%sf(j, k, l) = qK_cons_vf(beta_idx)%sf(j, k, l)
#endif

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_convert_conservative_to_primitive_variables

    !>  The following procedure handles the conversion between
        !!      the primitive variables and the conservative variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Conservative variables
        !!  @param gm_alphaK_vf Gradient magnitude of the volume fractions
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    impure subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                                    q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function
        ! and dynamic pressure, as defined in the incompressible flow sense,
        ! respectively
        real(wp) :: rho
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: qv
        real(wp) :: dyn_pres
        real(wp) :: nbub, R3tmp
        real(wp), dimension(nb) :: Rtmp
        real(wp) :: G
        real(wp), dimension(2) :: Re_K

        integer :: i, j, k, l !< Generic loop iterators

        real(wp), dimension(num_species) :: Ys
        real(wp) :: e_mix, mix_mol_weight, T
        real(wp) :: pres_mag

        real(wp) :: Ga ! Lorentz factor (gamma in relativity)
        real(wp) :: h ! relativistic enthalpy
        real(wp) :: v2 ! Square of the velocity magnitude
        real(wp) :: B2 ! Square of the magnetic field magnitude
        real(wp) :: vdotB ! Dot product of the velocity and magnetic field vectors
        real(wp) :: B(3) ! Magnetic field components

        pres_mag = 0._wp

        G = 0._wp

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Obtaining the density, specific heat ratio function
                    ! and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                                                        rho, gamma, pi_inf, qv, Re_K, G, fluid_pp(:)%G)

                    if (.not. igr .or. num_fluids > 1) then
                        ! Transferring the advection equation(s) variable(s)
                        do i = adv_idx%beg, adv_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (relativity) then

                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = q_prim_vf(B_idx%beg)%sf(j, k, l)
                            B(3) = q_prim_vf(B_idx%beg + 1)%sf(j, k, l)
                        else
                            B(1) = q_prim_vf(B_idx%beg)%sf(j, k, l)
                            B(2) = q_prim_vf(B_idx%beg + 1)%sf(j, k, l)
                            B(3) = q_prim_vf(B_idx%beg + 2)%sf(j, k, l)
                        end if

                        v2 = 0._wp
                        do i = momxb, momxe
                            v2 = v2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (v2 >= 1._wp) call s_mpi_abort('Error: v squared > 1 in s_convert_primitive_to_conservative_variables')

                        Ga = 1._wp/sqrt(1._wp - v2)

                        h = 1._wp + (gamma + 1)*q_prim_vf(E_idx)%sf(j, k, l)/rho ! Assume perfect gas for now

                        B2 = 0._wp
                        do i = B_idx%beg, B_idx%end
                            B2 = B2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (n == 0) B2 = B2 + Bx0**2

                        vdotB = 0._wp
                        do i = 1, 3
                            vdotB = vdotB + q_prim_vf(momxb + i - 1)%sf(j, k, l)*B(i)
                        end do

                        do i = 1, contxe
                            q_cons_vf(i)%sf(j, k, l) = Ga*q_prim_vf(i)%sf(j, k, l)
                        end do

                        do i = momxb, momxe
                            q_cons_vf(i)%sf(j, k, l) = (rho*h*Ga**2 + B2)*q_prim_vf(i)%sf(j, k, l) &
                                                       - vdotB*B(i - momxb + 1)
                        end do

                        q_cons_vf(E_idx)%sf(j, k, l) = rho*h*Ga**2 - q_prim_vf(E_idx)%sf(j, k, l) &
                                                       + 0.5_wp*(B2 + v2*B2 - vdotB**2)
                        ! Remove rest energy
                        do i = 1, contxe
                            q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) - q_cons_vf(i)%sf(j, k, l)
                        end do

                        do i = B_idx%beg, B_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do

                        cycle ! skip all the non-relativistic conversions below

                    end if

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, contxe
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed
                    ! iteratively by cycling through the velocity equations
                    dyn_pres = 0._wp

                    ! Computing momenta and dynamic pressure from velocity
                    do i = momxb, momxe
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)* &
                                   q_prim_vf(i)%sf(j, k, l)/2._wp
                    end do

                    if (chemistry) then
                        do i = chemxb, chemxe
                            Ys(i - chemxb + 1) = q_prim_vf(i)%sf(j, k, l)
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do

                        call get_mixture_molecular_weight(Ys, mix_mol_weight)
                        T = q_prim_vf(E_idx)%sf(j, k, l)*mix_mol_weight/(gas_constant*rho)
                        call get_mixture_energy_mass(T, Ys, e_mix)

                        q_cons_vf(E_idx)%sf(j, k, l) = &
                            dyn_pres + rho*e_mix
                    else
                        ! Computing the energy from the pressure
                        if (mhd) then
                            if (n == 0) then
                                pres_mag = 0.5_wp*(Bx0**2 + q_prim_vf(B_idx%beg)%sf(j, k, l)**2 + q_prim_vf(B_idx%beg + 1)%sf(j, k, l)**2)
                            else
                                pres_mag = 0.5_wp*(q_prim_vf(B_idx%beg)%sf(j, k, l)**2 + q_prim_vf(B_idx%beg + 1)%sf(j, k, l)**2 + q_prim_vf(B_idx%beg + 2)%sf(j, k, l)**2)
                            end if
                            q_cons_vf(E_idx)%sf(j, k, l) = &
                                gamma*q_prim_vf(E_idx)%sf(j, k, l) + dyn_pres + pres_mag &
                                + pi_inf + qv
                        elseif ((model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                            ! E = Gamma*P + \rho u u /2 + \pi_inf + (\alpha\rho qv)
                            q_cons_vf(E_idx)%sf(j, k, l) = &
                                gamma*q_prim_vf(E_idx)%sf(j, k, l) + dyn_pres + pi_inf + qv
                        else if ((model_eqns /= 4) .and. (bubbles_euler)) then
                            ! \tilde{E} = dyn_pres + (1-\alf)(\Gamma p_l + \Pi_inf)
                            q_cons_vf(E_idx)%sf(j, k, l) = dyn_pres + &
                                                           (1._wp - q_prim_vf(alf_idx)%sf(j, k, l))* &
                                                           (gamma*q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)
                        else
                            !Tait EOS, no conserved energy variable
                            q_cons_vf(E_idx)%sf(j, k, l) = 0._wp
                        end if
                    end if

                    ! Computing the internal energies from the pressure and continuities
                    if (model_eqns == 3) then
                        do i = 1, num_fluids
                            ! internal energy calculation for each of the fluids
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)* &
                                                                   (gammas(i)*q_prim_vf(E_idx)%sf(j, k, l) + pi_infs(i)) + &
                                                                   q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)
                        end do
                    end if

                    if (bubbles_euler) then
                        ! From prim: Compute nbub = (3/4pi) * \alpha / \bar{R^3}
                        do i = 1, nb
                            Rtmp(i) = q_prim_vf(bub_idx%rs(i))%sf(j, k, l)
                        end do

                        if (.not. qbmm) then
                            if (adv_n) then
                                q_cons_vf(n_idx)%sf(j, k, l) = q_prim_vf(n_idx)%sf(j, k, l)
                                nbub = q_prim_vf(n_idx)%sf(j, k, l)
                            else
                                call s_comp_n_from_prim(real(q_prim_vf(alf_idx)%sf(j, k, l), kind=wp), Rtmp, nbub, weight)
                            end if
                        else
                            !Initialize R3 averaging over R0 and R directions
                            R3tmp = 0._wp
                            do i = 1, nb
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) + sigR)**3._wp
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) - sigR)**3._wp
                            end do
                            !Initialize nb
                            nbub = 3._wp*q_prim_vf(alf_idx)%sf(j, k, l)/(4._wp*pi*R3tmp)
                        end if

                        if (j == 0 .and. k == 0 .and. l == 0) print *, 'In convert, nbub:', nbub

                        do i = bub_idx%beg, bub_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)*nbub
                        end do
                    end if

                    if (mhd) then
                        do i = B_idx%beg, B_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (elasticity) then
                        ! adding the elastic contribution
                        ! Multiply \tau to \rho \tau
                        do i = strxb, strxe
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        do i = strxb, strxe
                            ! adding elastic contribution
                            if (G > verysmall) then
                                if (cont_damage) G = G*max((1._wp - q_prim_vf(damage_idx)%sf(j, k, l)), 0._wp)

                                q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                               (q_prim_vf(i)%sf(j, k, l)**2._wp)/(4._wp*G)
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                                   (q_prim_vf(i)%sf(j, k, l)**2._wp)/(4._wp*G)
                                end if
                            end if
                        end do
                    end if

                    ! using \rho xi as the conservative formulation stated in Kamrin et al. JFM 2022
                    if (hyperelasticity) then
                        ! Multiply \xi to \rho \xi
                        do i = xibeg, xiend
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (surface_tension) then
                        q_cons_vf(c_idx)%sf(j, k, l) = q_prim_vf(c_idx)%sf(j, k, l)
                    end if

                    if (cont_damage) q_cons_vf(damage_idx)%sf(j, k, l) = q_prim_vf(damage_idx)%sf(j, k, l)

                end do
            end do
        end do
#else
        if (proc_rank == 0) then
            call s_mpi_abort('Conversion from primitive to '// &
                             'conservative variables not '// &
                             'implemented. Exiting.')
        end if
#endif
    end subroutine s_convert_primitive_to_conservative_variables

    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, &
                                                     FK_vf, &
                                                     FK_src_vf, &
                                                     is1, is2, is3, s2b, s3b)

        integer, intent(in) :: s2b, s3b
        real(wp), dimension(0:, s2b:, s3b:, 1:), intent(in) :: qK_prim_vf
        real(wp), dimension(0:, s2b:, s3b:, 1:), intent(inout) :: FK_vf
        real(wp), dimension(0:, s2b:, s3b:, advxb:), intent(inout) :: FK_src_vf

        type(int_bounds_info), intent(in) :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(wp), dimension(num_fluids) :: alpha_rho_K
        real(wp), dimension(num_fluids) :: alpha_K
        real(wp) :: rho_K
        real(wp), dimension(num_vels) :: vel_K
        real(wp) :: vel_K_sum
        real(wp) :: pres_K
        real(wp) :: E_K
        real(wp) :: gamma_K
        real(wp) :: pi_inf_K
        real(wp) :: qv_K
        real(wp), dimension(2) :: Re_K
        real(wp) :: G_K
        real(wp), dimension(num_species) :: Y_K
        real(wp) :: T_K, mix_mol_weight, R_gas

        integer :: i, j, k, l !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        $:GPU_UPDATE(device='[is1b,is2b,is3b,is1e,is2e,is3e]')

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity
#ifdef MFC_SIMULATION
        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_K, vel_K, alpha_K, Re_K, Y_K, rho_K, vel_K_sum, pres_K, E_K, gamma_K, pi_inf_K, qv_K, G_K, T_K, mix_mol_weight, R_gas]')
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, contxe
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advxb, advxe
                        alpha_K(i - E_idx) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K(i) = qK_prim_vf(j, k, l, contxe + i)
                    end do

                    vel_K_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K_sum = vel_K_sum + vel_K(i)**2._wp
                    end do

                    pres_K = qK_prim_vf(j, k, l, E_idx)
                    if (elasticity) then
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                        alpha_K, alpha_rho_K, Re_K, &
                                                                        G_K, Gs_vc)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                        alpha_K, alpha_rho_K, Re_K)
                    end if

                    ! Computing the energy from the pressure

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = chemxb, chemxe
                            Y_K(i - chemxb + 1) = qK_prim_vf(j, k, l, i)
                        end do
                        !Computing the energy from the internal energy of the mixture
                        call get_mixture_molecular_weight(Y_k, mix_mol_weight)
                        R_gas = gas_constant/mix_mol_weight
                        T_K = pres_K/rho_K/R_gas
                        call get_mixture_energy_mass(T_K, Y_K, E_K)
                        E_K = rho_K*E_K + 5.e-1_wp*rho_K*vel_K_sum
                    else
                        ! Computing the energy from the pressure
                        E_K = gamma_K*pres_K + pi_inf_K &
                              + 5.e-1_wp*rho_K*vel_K_sum + qv_K
                    end if

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, contxe
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        FK_vf(j, k, l, contxe + dir_idx(i)) = &
                            rho_K*vel_K(dir_idx(1)) &
                            *vel_K(dir_idx(i)) &
                            + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, E_idx) = vel_K(dir_idx(1))*(E_K + pres_K)

                    ! Species advection Flux, \rho*u*Y
                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            FK_vf(j, k, l, i - 1 + chemxb) = vel_K(dir_idx(1))*(rho_K*Y_K(i))
                        end do
                    end if

                    if (riemann_solver == 1 .or. riemann_solver == 4) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = 0._wp
                            FK_src_vf(j, k, l, i) = alpha_K(i - E_idx)
                        end do

                    else
                        ! Could be bubbles_euler!
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = vel_K(dir_idx(1))*alpha_K(i - E_idx)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = advxb, advxe
                            FK_src_vf(j, k, l, i) = vel_K(dir_idx(1))
                        end do

                    end if

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#endif
    end subroutine s_convert_primitive_to_flux_variables

    !>  This subroutine computes partial densities and volume fractions
    subroutine s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)
        $:GPU_ROUTINE(function_name='s_compute_species_fraction', &
            & parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in) :: k, l, r
        real(wp), dimension(num_fluids), intent(out) :: alpha_rho_K, alpha_K
        integer :: i

        if (num_fluids == 1) then
            alpha_rho_K(1) = q_vf(contxb)%sf(k, l, r)
            if (igr .or. bubbles_euler) then
                alpha_K(1) = 1._wp
            else
                alpha_K(1) = q_vf(advxb)%sf(k, l, r)
            end if
        else
            if (igr) then
                do i = 1, num_fluids - 1
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(advxb + i - 1)%sf(k, l, r)
                end do
                alpha_rho_K(num_fluids) = q_vf(num_fluids)%sf(k, l, r)
                alpha_K(num_fluids) = 1._wp - sum(alpha_K(1:num_fluids - 1))
            else
                do i = 1, num_fluids
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(advxb + i - 1)%sf(k, l, r)
                end do
            end if
        end if

        if (mpp_lim) then
            do i = 1, num_fluids
                alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
            end do
            alpha_K = alpha_K/max(sum(alpha_K), 1.e-16_wp)
        end if

        if (num_fluids == 1 .and. bubbles_euler) alpha_K(1) = q_vf(advxb)%sf(k, l, r)

    end subroutine s_compute_species_fraction

    impure subroutine s_finalize_variables_conversion_module()

        ! Deallocating the density, the specific heat ratio function and the
        ! liquid stiffness function
#ifdef MFC_POST_PROCESS
        deallocate (rho_sf, gamma_sf, pi_inf_sf, qv_sf)
#endif

#ifdef MFC_SIMULATION
        @:DEALLOCATE(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs_vc)
        if (bubbles_euler) then
            @:DEALLOCATE(bubrs_vc)
        end if
#else
        @:DEALLOCATE(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs_vc)
        if (bubbles_euler) then
            @:DEALLOCATE(bubrs_vc)
        end if
#endif

    end subroutine s_finalize_variables_conversion_module

#ifndef MFC_PRE_PROCESS
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c_c, c, qv)
        $:GPU_ROUTINE(function_name='s_compute_speed_of_sound', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: pres
        real(wp), intent(in) :: rho, gamma, pi_inf, qv
        real(wp), intent(in) :: H
        real(wp), dimension(num_fluids), intent(in) :: adv
        real(wp), intent(in) :: vel_sum
        real(wp), intent(in) :: c_c
        real(wp), intent(out) :: c

        real(wp) :: blkmod1, blkmod2

        integer :: q

        if (chemistry) then
            if (avg_state == 1 .and. abs(c_c) > verysmall) then
                c = sqrt(c_c - (gamma - 1.0_wp)*(vel_sum - H))
            else
                c = sqrt((1.0_wp + 1.0_wp/gamma)*pres/rho)
            end if
        elseif (relativity) then
            ! Only supports perfect gas for now
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/H)
        else
            if (alt_soundspeed) then
                blkmod1 = ((gammas(1) + 1._wp)*pres + &
                           pi_infs(1))/gammas(1)
                blkmod2 = ((gammas(2) + 1._wp)*pres + &
                           pi_infs(2))/gammas(2)
                c = (1._wp/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
            elseif (model_eqns == 3) then
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    c = c + adv(q)*gs_min(q)* &
                        (pres + pi_infs(q)/(gammas(q) + 1._wp))
                end do
                c = c/rho
            elseif (((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles_euler))) then
                ! Sound speed for bubble mixture to order O(\alpha)

                if (mpp_lim .and. (num_fluids > 1)) then
                    c = (1._wp/gamma + 1._wp)* &
                        (pres + pi_inf/(gamma + 1._wp))/rho
                else
                    c = &
                        (1._wp/gamma + 1._wp)* &
                        (pres + pi_inf/(gamma + 1._wp))/ &
                        (rho*(1._wp - adv(num_fluids)))
                end if
            else
                c = (H - 5.e-1*vel_sum - qv/rho)/gamma
            end if

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if
    end subroutine s_compute_speed_of_sound
#endif

#ifndef MFC_PRE_PROCESS
    subroutine s_compute_fast_magnetosonic_speed(rho, c, B, norm, c_fast, h)
        $:GPU_ROUTINE(function_name='s_compute_fast_magnetosonic_speed', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: B(3), rho, c
        real(wp), intent(in) :: h ! only used for relativity
        real(wp), intent(out) :: c_fast
        integer, intent(in) :: norm

        real(wp) :: B2, term, disc

        B2 = sum(B**2)

        if (.not. relativity) then
            term = c**2 + B2/rho
            disc = term**2 - 4*c**2*(B(norm)**2/rho)
        else
            ! Note: this is approximation for the non-relatisitic limit; accurate solution requires solving a quartic equation
            term = (c**2*(B(norm)**2 + rho*h) + B2)/(rho*h + B2)
            disc = term**2 - 4*c**2*B(norm)**2/(rho*h + B2)
        end if

#ifdef DEBUG
        if (disc < 0._wp) then
            print *, 'rho, c, Bx, By, Bz, h, term, disc:', rho, c, B(1), B(2), B(3), h, term, disc
            call s_mpi_abort('Error: negative discriminant in s_compute_fast_magnetosonic_speed')
        end if
#endif

        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))

    end subroutine s_compute_fast_magnetosonic_speed
#endif

end module m_variables_conversion
