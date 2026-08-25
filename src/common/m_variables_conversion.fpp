!>
!! @file
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Conservative-to-primitive variable conversion, mixture property evaluation, and pressure computation
module m_variables_conversion

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper
    use m_constants, only: riemann_solver_hll, riemann_solver_hlld, model_eqns_gamma_law, model_eqns_5eq, model_eqns_6eq, &
        & avg_state_roe
    use m_thermochem, only: num_species, get_temperature, get_pressure, gas_constant, get_mixture_molecular_weight, &
        & get_mixture_energy_mass

    implicit none

    private
    public :: s_initialize_variables_conversion_module, s_initialize_pb, s_initialize_mv, s_convert_to_mixture_variables, &
        & s_convert_mixture_to_mixture_variables, s_convert_species_to_mixture_variables, &
        & s_convert_species_to_mixture_variables_kernel, s_convert_conservative_to_primitive_variables, &
        & s_convert_primitive_to_conservative_variables, s_convert_primitive_to_flux_variables, s_compute_pressure, &
        & s_compute_species_fraction, s_accumulate_mixture_properties, s_compute_energy, s_compute_speed_of_sound, &
        & s_compute_speed_of_sound_avg, &
        & s_compute_fast_magnetosonic_speed, s_finalize_variables_conversion_module, gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps

    real(wp), allocatable, dimension(:)   :: Gs_vc
    integer, allocatable, dimension(:)    :: bubrs_vc
    real(wp), allocatable, dimension(:,:) :: Res_vc
    $:GPU_DECLARE(create='[bubrs_vc, Gs_vc, Res_vc]')

    integer :: is1b, is2b, is3b, is1e, is2e, is3e
    $:GPU_DECLARE(create='[is1b, is2b, is3b, is1e, is2e, is3e]')

    logical :: enforce_density_floor_vc = .false.
    logical :: preserve_qbmm_number_vc = .false.
    integer :: lagrange_beta_index_vc = 0
    $:GPU_DECLARE(create='[enforce_density_floor_vc, preserve_qbmm_number_vc, lagrange_beta_index_vc]')

    real(wp), allocatable, dimension(:,:,:), public :: rho_sf     !< Scalar density function
    real(wp), allocatable, dimension(:,:,:), public :: gamma_sf   !< Scalar sp. heat ratio function
    real(wp), allocatable, dimension(:,:,:), public :: pi_inf_sf  !< Scalar liquid stiffness function

contains

    !> Dispatch to the s_convert_mixture_to_mixture_variables and s_convert_species_to_mixture_variables subroutines. Replaces a
    !! procedure pointer.
    subroutine s_convert_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in)   :: q_vf
        integer, intent(in)                                   :: i, j, k
        real(wp), intent(out), target                         :: rho, gamma, pi_inf, qv
        real(wp), optional, dimension(2), intent(out)         :: Re_K
        real(wp), optional, intent(out)                       :: G_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G

        if (model_eqns == model_eqns_gamma_law) then  ! Gamma/pi_inf model
            call s_convert_mixture_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv)
        else  ! Volume fraction model
            call s_convert_species_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K, G_K, G)
        end if

    end subroutine s_convert_to_mixture_variables

    !> Compute the pressure from the appropriate equation of state
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, rhoYks, pres, T, stress, mom, G, pres_mag)

        $:GPU_ROUTINE(function_name='s_compute_pressure',parallelism='[seq]', cray_noinline=True)

        real(stp), intent(in)           :: energy, alf
        real(wp), intent(in)            :: dyn_p
        real(wp), intent(in)            :: pi_inf, gamma, rho, qv
        real(wp), intent(out)           :: pres
        real(wp), intent(inout)         :: T
        real(stp), intent(in), optional :: stress, mom
        real(wp), intent(in), optional  :: G, pres_mag

        ! Chemistry
        real(wp), dimension(1:num_species), intent(in) :: rhoYks
        real(wp), dimension(1:num_species)             :: Y_rs
        real(wp)                                       :: E_e
        real(wp)                                       :: e_Per_Kg, Pdyn_Per_Kg
        real(wp)                                       :: T_guess
        integer                                        :: s  !< Generic loop iterator
        #:if not chemistry
            ! Depending on model_eqns and bubbles_euler, the appropriate procedure for computing pressure is targeted by the
            ! procedure pointer

            if (mhd) then
                ! MHD pressure: subtract magnetic pressure from total energy
                pres = (energy - dyn_p - pi_inf - qv - pres_mag)/gamma
            else if (bubbles_euler .neqv. .true.) then
                ! Gamma/pi_inf model or five-equation model (Allaire et al. JCP 2002): p from mixture EOS
                pres = (energy - dyn_p - pi_inf - qv)/gamma
            else
                ! Bubble-augmented pressure with void fraction correction
                pres = ((energy - dyn_p)/(1._wp - alf) - pi_inf - qv)/gamma
            end if

            if (hypoelasticity .and. present(G)) then
                ! Subtract elastic strain energy before computing pressure (hypoelastic model)
                E_e = 0._wp
                do s = eqn_idx%stress%beg, eqn_idx%stress%end
                    if (G > 0) then
                        E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        ! Double for shear stresses
                        if (any(s == shear_indices)) then
                            E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        end if
                    end if
                end do

                pres = (energy - 0.5_wp*(mom**2._wp)/rho - pi_inf - qv - E_e)/gamma
            end if
        #:else
            ! Reacting mixture pressure from temperature and species
            Y_rs(:) = rhoYks(:)/rho
            e_Per_Kg = energy/rho
            Pdyn_Per_Kg = dyn_p/rho

            T_guess = T

            call get_temperature(e_Per_Kg - Pdyn_Per_Kg, T_guess, Y_rs, .true., T)
            call get_pressure(rho, T, Y_rs, pres)
        #:endif

    end subroutine s_compute_pressure

    !> Convert mixture variables to density, gamma, pi_inf, and qv for the gamma/pi_inf model. Given conservative or primitive
    !! variables, transfers the density, specific heat ratio function and the liquid stiffness function from q_vf to rho, gamma and
    !! pi_inf.
    subroutine s_convert_mixture_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: i, j, k
        real(wp), intent(out), target                       :: rho
        real(wp), intent(out), target                       :: gamma
        real(wp), intent(out), target                       :: pi_inf
        real(wp), intent(out), target                       :: qv

        ! Transferring the density, the specific heat ratio function and the liquid stiffness function, respectively

        rho = q_vf(1)%sf(i, j, k)
        gamma = q_vf(eqn_idx%gamma)%sf(i, j, k)
        pi_inf = q_vf(eqn_idx%pi_inf)%sf(i, j, k)
        qv = 0._wp  ! keep this value nil for now. For future adjustment

        ! Store derived mixture fields when requested during module initialization.
        if (allocated(rho_sf)) then
            rho_sf(i, j, k) = rho
            gamma_sf(i, j, k) = gamma
            pi_inf_sf(i, j, k) = pi_inf
        end if

    end subroutine s_convert_mixture_to_mixture_variables

    !> Convert species volume fractions and partial densities to mixture density, gamma, pi_inf, and qv. Given conservative or
    !! primitive variables, computes the density, the specific heat ratio function and the liquid stiffness function from q_vf and
    !! stores the results into rho, gamma and pi_inf.
    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, r, rho, gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in)   :: q_vf
        integer, intent(in)                                   :: k, l, r
        real(wp), intent(out), target                         :: rho
        real(wp), intent(out), target                         :: gamma
        real(wp), intent(out), target                         :: pi_inf
        real(wp), intent(out), target                         :: qv
        real(wp), optional, dimension(2), intent(out)         :: Re_K
        real(wp), optional, intent(out)                       :: G_K
        real(wp), dimension(num_fluids)                       :: alpha_rho_K, alpha_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G
        integer                                               :: i, j  !< Generic loop iterator
        ! Computing the density, the specific heat ratio function and the liquid stiffness function, respectively

        call s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        ! Use the same scalar kernel on host and device so mixture semantics do not depend on the executable or accelerator backend.
        ! Absent optional dummies forward as absent, so the optional arguments need no dispatch here.
        call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha_K, alpha_rho_K, Re_K, G_K, G)

        ! Store derived mixture fields when requested during module initialization.
        if (allocated(rho_sf)) then
            rho_sf(k, l, r) = rho
            gamma_sf(k, l, r) = gamma
            pi_inf_sf(k, l, r) = pi_inf
        end if

    end subroutine s_convert_species_to_mixture_variables

    !> Host- and device-callable conversion kernel for species and mixture variables.
    subroutine s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K, G_K, G)

        $:GPU_ROUTINE(function_name='s_convert_species_to_mixture_variables_kernel', parallelism='[seq]', cray_noinline=True)

        real(wp), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(inout)        :: alpha_rho_K, alpha_K
            real(wp), optional, dimension(3), intent(in) :: G
        #:else
            real(wp), dimension(num_fluids), intent(inout)        :: alpha_rho_K, alpha_K
            real(wp), optional, dimension(num_fluids), intent(in) :: G
        #:endif
        real(wp), optional, dimension(2), intent(out) :: Re_K
        real(wp), optional, intent(out)               :: G_K
        real(wp)                                      :: alpha_K_sum
        integer                                       :: i, j  !< Generic loop iterators

        rho_K = 0._wp
        gamma_K = 0._wp
        pi_inf_K = 0._wp
        qv_K = 0._wp
        if (present(Re_K)) Re_K = dflt_real
        if (present(G_K)) G_K = 0._wp

        ! Constrain partial densities and volume fractions within physical bounds
        if (num_fluids == 1 .and. bubbles_euler) then
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
            qv_K = qvs(1)
        else
            if (mpp_lim) then
                alpha_K_sum = 0._wp
                do i = 1, num_fluids
                    alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                    alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                    alpha_K_sum = alpha_K_sum + alpha_K(i)
                end do
                alpha_K = alpha_K/max(alpha_K_sum, sgm_eps)
            end if
            call s_accumulate_mixture_properties(num_fluids, alpha_rho_K, alpha_K, rho_K, gamma_K, pi_inf_K, qv_K)
        end if

        if (present(G_K)) then
            G_K = 0._wp
            do i = 1, num_fluids
                ! TODO: change to use Gs_vc directly here? TODO: Make this change as well for GPUs
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0._wp, G_K)
        end if

        if (viscous .and. present(Re_K)) then
            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res_vc(i, j) + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if

    end subroutine s_convert_species_to_mixture_variables_kernel

    !> Initialize the variables conversion module.
    impure subroutine s_initialize_variables_conversion_module(store_mixture_fields, enforce_density_floor, preserve_qbmm_number, &
        & lagrange_beta_index)

        integer                       :: i, j
        logical, optional, intent(in) :: store_mixture_fields
        logical, optional, intent(in) :: enforce_density_floor, preserve_qbmm_number
        integer, optional, intent(in) :: lagrange_beta_index
        logical                       :: allocate_mixture_fields

        allocate_mixture_fields = .false.
        if (present(store_mixture_fields)) allocate_mixture_fields = store_mixture_fields
        enforce_density_floor_vc = .false.
        if (present(enforce_density_floor)) enforce_density_floor_vc = enforce_density_floor
        preserve_qbmm_number_vc = .false.
        if (present(preserve_qbmm_number)) preserve_qbmm_number_vc = preserve_qbmm_number
        lagrange_beta_index_vc = 0
        if (present(lagrange_beta_index)) lagrange_beta_index_vc = lagrange_beta_index

        $:GPU_ENTER_DATA(copyin='[is1b, is1e, is2b, is2e, is3b, is3e]')
        $:GPU_UPDATE(device='[enforce_density_floor_vc, preserve_qbmm_number_vc, lagrange_beta_index_vc]')

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
        $:GPU_UPDATE(device='[gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs_vc]')

        @:ALLOCATE(Res_vc(1:2, 1:max(1, Re_size_max)))
        Res_vc = dflt_real
        if (allocated(Re_idx)) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_vc(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Re_idx]')
        end if
        $:GPU_UPDATE(device='[Res_vc, Re_size]')

        if (bubbles_euler) then
            @:ALLOCATE(bubrs_vc(1:nb))
            do i = 1, nb
                bubrs_vc(i) = qbmm_idx%rs(i)
            end do
            $:GPU_UPDATE(device='[bubrs_vc]')
        end if

        if (allocate_mixture_fields) then
            ! Allocate derived mixture fields over the available grid storage.
            if (n > 0) then
                if (p > 0) then
                    allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                    allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                    allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                else
                    allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                    allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                    allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                end if
            else
                allocate (rho_sf(-buff_size:m + buff_size,0:0,0:0))
                allocate (gamma_sf(-buff_size:m + buff_size,0:0,0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size,0:0,0:0))
            end if
        end if

    end subroutine s_initialize_variables_conversion_module

    !> Initialize bubble mass-vapor values at quadrature nodes from the conserved moment statistics.
    subroutine s_initialize_mv(qK_cons_vf, mv)

        type(scalar_field), dimension(sys_size), intent(in)                                     :: qK_cons_vf
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(inout) :: mv
        integer                                                                                 :: i, j, k, l
        real(wp)                                                                                :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end
                    nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(eqn_idx%bub%beg + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(eqn_idx%bub%beg + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        mv(j, k, l, 1, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 2, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 3, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 4, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_mv

    !> Initialize bubble internal pressures at quadrature nodes using isothermal relations from the Preston model.
    subroutine s_initialize_pb(qK_cons_vf, mv, pb)

        type(scalar_field), dimension(sys_size), intent(in)                                     :: qK_cons_vf
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(in)    :: mv
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(inout) :: pb
        integer                                                                                 :: i, j, k, l
        real(wp)                                                                                :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end
                    nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(eqn_idx%bub%beg + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(eqn_idx%bub%beg + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        ! PRESTON (ISOTHERMAL)
                        pb(j, k, l, 1, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 1, &
                           & i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 2, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 2, &
                           & i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 3, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 3, &
                           & i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 4, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 4, &
                           & i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_pb

    !> Convert conserved variables (rho*alpha, rho*u, E, alpha) to primitives (rho, u, p, alpha). Conversion depends on model_eqns:
    !! each model has different variable sets and EOS.
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, q_T_sf, qK_prim_vf, ibounds)

        use m_global_parameters_common, only: shear_indices  ! Performance fix with AMDFlang

        type(scalar_field), dimension(sys_size), intent(in)    :: qK_cons_vf
        type(scalar_field), intent(inout)                      :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: qK_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in)      :: ibounds

        #:if USING_AMD and not MFC_CASE_OPTIMIZATION
            real(wp), dimension(3) :: alpha_K, alpha_rho_K
            real(wp), dimension(3) :: nRtmp
            real(wp)               :: rhoYks(1:10)
        #:else
            real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
            real(wp), dimension(nb)         :: nRtmp
            real(wp)                        :: rhoYks(1:num_species)
        #:endif
        real(wp), dimension(2) :: Re_K
        real(wp)               :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K
        real(wp)               :: vftmp, nbub_sc
        real(wp)               :: G_K
        real(wp)               :: pres
        integer                :: i, j, k, l               !< Generic loop iterators
        real(wp)               :: T
        real(wp)               :: pres_mag
        real(wp)               :: Ga                       !< Lorentz factor (gamma in relativity)
        real(wp)               :: B2                       !< Magnetic field magnitude squared
        real(wp)               :: B(3)                     !< Magnetic field components
        real(wp)               :: m2                       !< Relativistic momentum magnitude squared
        real(wp)               :: S                        !< Dot product of the magnetic field and the relativistic momentum
        real(wp)               :: W, dW                    !< W := rho*v*Ga**2; f = f(W) in Newton-Raphson
        real(wp)               :: E, D                     !< Prim/Cons variables within Newton-Raphson iteration
        real(wp)               :: f, dGa_dW, dp_dW, df_dW  !< Functions within Newton-Raphson iteration
        integer                :: iter                     !< Newton-Raphson iteration counter

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_K, alpha_rho_K, Re_K, nRtmp, rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, &
                            & rhoYks, B, pres, vftmp, nbub_sc, G_K, T, pres_mag, Ga, B2, m2, S, W, dW, E, D, f, dGa_dW, dp_dW, &
                            & df_dW, iter]')
        do l = ibounds(3)%beg, ibounds(3)%end
            do k = ibounds(2)%beg, ibounds(2)%end
                do j = ibounds(1)%beg, ibounds(1)%end
                    dyn_pres_K = 0._wp

                    call s_compute_species_fraction(qK_cons_vf, j, k, l, alpha_rho_K, alpha_K)

#ifdef MFC_GPU
                    ! Device regions call the device-compiled scalar kernel directly.
                    if (hypoelasticity) then
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K, G_K, Gs_vc)
                    else
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K)
                    end if
#else
                    ! Host execution uses the wrapper, which also stores requested diagnostics.
                    if (hypoelasticity) then
                        call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, rho_K, gamma_K, pi_inf_K, qv_K, Re_K, G_K, &
                                                            & fluid_pp(:)%G)
                    else
                        call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, rho_K, gamma_K, pi_inf_K, qv_K)
                    end if
#endif

                    ! Relativistic MHD primitive variable recovery, Mignone & Bodo A&A (2006)
                    if (relativity) then
                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(3) = qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                        else
                            B(1) = qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(2) = qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                            B(3) = qK_cons_vf(eqn_idx%B%beg + 2)%sf(j, k, l)
                        end if
                        B2 = B(1)**2 + B(2)**2 + B(3)**2

                        m2 = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            m2 = m2 + qK_cons_vf(i)%sf(j, k, l)**2
                        end do

                        S = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            S = S + qK_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)*B(i)
                        end do

                        E = qK_cons_vf(eqn_idx%E)%sf(j, k, l)

                        D = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            D = D + qK_cons_vf(i)%sf(j, k, l)
                        end do

                        ! Newton-Raphson
                        W = E + D
                        $:GPU_LOOP(parallelism='[seq]')
                        do iter = 1, relativity_cons_to_prim_max_iter
                            ! Lorentz factor from total enthalpy and magnetic field
                            Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                            ! Thermal pressure from EOS
                            pres = (W - D*Ga)/((gamma_K + 1)*Ga**2)
                            f = W - pres + (1 - 1/(2*Ga**2))*B2 - S**2/(2*W**2) - E - D

                            ! The first equation below corrects a typo in (Mignone & Bodo, 2006) m2*W**2 -> 2*m2*W**2, which would
                            ! cancel with the 2* in other terms This corrected version is not used as the second equation
                            ! empirically converges faster. First equation is kept for further investigation. dGa_dW = -Ga**3 * (
                            ! S**2*(3*W**2+3*W*B2+B2**2) + m2*W**2 ) / (W**3 * (W+B2)**3) ! first (corrected)
                            dGa_dW = -Ga**3*(2*S**2*(3*W**2 + 3*W*B2 + B2**2) + m2*W**2)/(2*W**3*(W + B2)**3)  ! second (in paper)

                            dp_dW = (Ga*(1 + D*dGa_dW) - 2*W*dGa_dW)/((gamma_K + 1)*Ga**3)
                            df_dW = 1 - dp_dW + (B2/Ga**3)*dGa_dW + S**2/W**3

                            dW = -f/df_dW
                            W = W + dW
                            if (abs(dW) < 1.e-12_wp*W) exit  ! Relative convergence criterion
                        end do

                        ! Recalculate pressure using converged W
                        Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                        qK_prim_vf(eqn_idx%E)%sf(j, k, l) = (W - D*Ga)/((gamma_K + 1)*Ga**2)

                        ! Recover the other primitive variables
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            qK_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = (qK_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                       & l) + (S/W)*B(i))/(W + B2)
                        end do
                        qK_prim_vf(1)%sf(j, k, l) = D/Ga  ! Hard-coded for single-component for now

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do

                        cycle  ! skip all the non-relativistic conversions below
                    end if

                    if (chemistry) then
                        ! Reacting flow: recover density from species partial densities, compute mass fractions Y_k = rhoY_k / rho
                        rho_K = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            rho_K = rho_K + max(0._wp, qK_cons_vf(i)%sf(j, k, l))
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            qK_prim_vf(i)%sf(j, k, l) = rho_K
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            qK_prim_vf(i)%sf(j, k, l) = max(0._wp, qK_cons_vf(i)%sf(j, k, l)/rho_K)
                        end do
                    else
                        ! Non-reacting: partial densities are directly primitive (alpha_i * rho_i)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (enforce_density_floor_vc) rho_K = max(rho_K, sgm_eps)

                    ! Recover velocity from momentum: u = rho*u / rho, and accumulate dynamic pressure 0.5*rho*|u|^2
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        dyn_pres_K = dyn_pres_K + 5.e-1_wp*qK_cons_vf(i)%sf(j, k, l)*qK_prim_vf(i)%sf(j, k, l)
                    end do

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            rhoYks(i) = qK_cons_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l)
                        end do

                        T = q_T_sf%sf(j, k, l)
                    end if

                    if (mhd) then
                        if (n == 0) then
                            pres_mag = 0.5_wp*(Bx0**2 + qK_cons_vf(eqn_idx%B%beg)%sf(j, k, &
                                               & l)**2 + qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)**2)
                        else
                            pres_mag = 0.5_wp*(qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, &
                                               & l)**2 + qK_cons_vf(eqn_idx%B%beg + 2)%sf(j, k, l)**2)
                        end if
                    else
                        pres_mag = 0._wp
                    end if

                    call s_compute_pressure(qK_cons_vf(eqn_idx%E)%sf(j, k, l), qK_cons_vf(eqn_idx%alf)%sf(j, k, l), dyn_pres_K, &
                                            & pi_inf_K, gamma_K, rho_K, qv_K, rhoYks, pres, T, pres_mag=pres_mag)

                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = pres

                    if (chemistry) then
                        q_T_sf%sf(j, k, l) = T
                    end if

                    if (bubbles_euler) then
                        ! Recover bubble primitive variables: divide conserved moments by bubble number density
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, nb
                            nRtmp(i) = qK_cons_vf(bubrs_vc(i))%sf(j, k, l)
                        end do

                        vftmp = qK_cons_vf(eqn_idx%alf)%sf(j, k, l)

                        if (qbmm) then
                            ! Get nb (constant across all R0 bins)
                            nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                            ! Convert cons to prim
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                            ! Need to keep track of nb in the primitive variable list (converted back to true value before output)
                            if (preserve_qbmm_number_vc) then
                                qK_prim_vf(eqn_idx%bub%beg)%sf(j, k, l) = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)
                            end if
                        else
                            if (adv_n) then
                                qK_prim_vf(eqn_idx%n)%sf(j, k, l) = qK_cons_vf(eqn_idx%n)%sf(j, k, l)
                                nbub_sc = qK_prim_vf(eqn_idx%n)%sf(j, k, l)
                            else
                                call s_comp_n_from_cons(vftmp, nRtmp, nbub_sc, weight)
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                        end if
                    end if

                    if (mhd) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (hypoelasticity) then
                        if (cont_damage) G_K = G_K*max((1._wp - qK_cons_vf(eqn_idx%damage)%sf(j, k, l)), 0._wp)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            ! Elastic energy subtraction (guard skips when G near zero from alpha undershoot)
                            if (G_K > verysmall) then
                                qK_prim_vf(eqn_idx%E)%sf(j, k, l) = qK_prim_vf(eqn_idx%E)%sf(j, k, l) - ((qK_prim_vf(i)%sf(j, k, &
                                           & l)**2._wp)/max(4._wp*G_K, verysmall))/gamma_K
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = qK_prim_vf(eqn_idx%E)%sf(j, k, l) - ((qK_prim_vf(i)%sf(j, &
                                               & k, l)**2._wp)/max(4._wp*G_K, verysmall))/gamma_K
                                end if
                            end if
                        end do
                    end if

                    if (.not. igr .or. num_fluids > 1) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (surface_tension) then
                        qK_prim_vf(eqn_idx%c)%sf(j, k, l) = qK_cons_vf(eqn_idx%c)%sf(j, k, l)
                    end if

                    if (cont_damage) qK_prim_vf(eqn_idx%damage)%sf(j, k, l) = qK_cons_vf(eqn_idx%damage)%sf(j, k, l)

                    if (hyper_cleaning) qK_prim_vf(eqn_idx%psi)%sf(j, k, l) = qK_cons_vf(eqn_idx%psi)%sf(j, k, l)
                    if (bubbles_lagrange .and. lagrange_beta_index_vc > 0) then
                        qK_prim_vf(lagrange_beta_index_vc)%sf(j, k, l) = qK_cons_vf(lagrange_beta_index_vc)%sf(j, k, l)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_convert_conservative_to_primitive_variables

    !> Convert primitives (rho, u, p, alpha) to conserved variables (rho*alpha, rho*u, E, alpha).
    impure subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function and dynamic pressure, as defined in the incompressible
        ! flow sense, respectively
        real(wp)                         :: rho
        real(wp)                         :: gamma
        real(wp)                         :: pi_inf
        real(wp)                         :: qv
        real(wp)                         :: dyn_pres
        real(wp)                         :: nbub, R3tmp
        real(wp), dimension(nb)          :: Rtmp
        real(wp)                         :: G
        real(wp), dimension(2)           :: Re_K
        integer                          :: i, j, k, l  !< Generic loop iterators
        real(wp), dimension(num_species) :: Ys
        real(wp)                         :: e_mix, mix_mol_weight, T
        real(wp)                         :: pres_mag
        real(wp)                         :: Ga          !< Lorentz factor (gamma in relativity)
        real(wp)                         :: h           !< relativistic enthalpy
        real(wp)                         :: v2          !< Square of the velocity magnitude
        real(wp)                         :: B2          !< Square of the magnetic field magnitude
        real(wp)                         :: vdotB       !< Dot product of the velocity and magnetic field vectors
        real(wp)                         :: B(3)        !< Magnetic field components

        pres_mag = 0._wp

        G = 0._wp

        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    ! Obtaining the density, specific heat ratio function and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, rho, gamma, pi_inf, qv, Re_K, G, fluid_pp(:)%G)

                    if (.not. igr .or. num_fluids > 1) then
                        ! Transferring the advection equation(s) variable(s)
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (relativity) then
                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(3) = q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                        else
                            B(1) = q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(2) = q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                            B(3) = q_prim_vf(eqn_idx%B%beg + 2)%sf(j, k, l)
                        end if

                        v2 = 0._wp
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            v2 = v2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (v2 >= 1._wp) call s_mpi_abort('Error: v squared > 1 in s_convert_primitive_to_conservative_variables')

                        Ga = 1._wp/sqrt(1._wp - v2)

                        h = 1._wp + (gamma + 1)*q_prim_vf(eqn_idx%E)%sf(j, k, l)/rho  ! Assume perfect gas for now

                        B2 = 0._wp
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            B2 = B2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (n == 0) B2 = B2 + Bx0**2

                        vdotB = 0._wp
                        do i = 1, 3
                            vdotB = vdotB + q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)*B(i)
                        end do

                        do i = 1, eqn_idx%cont%end
                            q_cons_vf(i)%sf(j, k, l) = Ga*q_prim_vf(i)%sf(j, k, l)
                        end do

                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            q_cons_vf(i)%sf(j, k, l) = (rho*h*Ga**2 + B2)*q_prim_vf(i)%sf(j, k, &
                                      & l) - vdotB*B(i - eqn_idx%mom%beg + 1)
                        end do

                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = rho*h*Ga**2 - q_prim_vf(eqn_idx%E)%sf(j, k, &
                                  & l) + 0.5_wp*(B2 + v2*B2 - vdotB**2)
                        ! Remove rest energy
                        do i = 1, eqn_idx%cont%end
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) - q_cons_vf(i)%sf(j, k, l)
                        end do

                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do

                        cycle  ! skip all the non-relativistic conversions below
                    end if

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, eqn_idx%cont%end
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed iteratively by cycling through the velocity equations
                    dyn_pres = 0._wp

                    ! Computing momenta and dynamic pressure from velocity
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)*q_prim_vf(i)%sf(j, k, l)/2._wp
                    end do

                    if (chemistry) then
                        ! Reacting mixture: compute conserved energy from species mass fractions and temperature
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            Ys(i - eqn_idx%species%beg + 1) = q_prim_vf(i)%sf(j, k, l)
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do

                        call get_mixture_molecular_weight(Ys, mix_mol_weight)
                        T = q_prim_vf(eqn_idx%E)%sf(j, k, l)*mix_mol_weight/(gas_constant*rho)
                        call get_mixture_energy_mass(T, Ys, e_mix)

                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + rho*e_mix
                    else
                        ! Computing the energy from the pressure
                        if (mhd) then
                            if (n == 0) then
                                pres_mag = 0.5_wp*(Bx0**2 + q_prim_vf(eqn_idx%B%beg)%sf(j, k, &
                                                   & l)**2 + q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)**2)
                            else
                                pres_mag = 0.5_wp*(q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + q_prim_vf(eqn_idx%B%beg + 1)%sf(j, &
                                                   & k, l)**2 + q_prim_vf(eqn_idx%B%beg + 2)%sf(j, k, l)**2)
                            end if
                            ! MHD energy includes magnetic pressure contribution
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = gamma*q_prim_vf(eqn_idx%E)%sf(j, k, &
                                      & l) + dyn_pres + pres_mag + pi_inf + qv
                        else if (bubbles_euler .neqv. .true.) then
                            ! Five-equation model (Allaire et al. JCP 2002): E = Gamma*p + 0.5*rho*|u|^2 + pi_inf + qv
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = gamma*q_prim_vf(eqn_idx%E)%sf(j, k, l) + dyn_pres + pi_inf + qv
                        else
                            ! Bubble-augmented energy with void fraction correction
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + (1._wp - q_prim_vf(eqn_idx%alf)%sf(j, k, &
                                      & l))*(gamma*q_prim_vf(eqn_idx%E)%sf(j, k, l) + pi_inf)
                        end if
                    end if

                    ! Six-equation model (Saurel et al. JCP 2009): compute per-phase internal energies
                    if (model_eqns == model_eqns_6eq) then
                        do i = 1, num_fluids
                            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                                      & l)*(gammas(i)*q_prim_vf(eqn_idx%E)%sf(j, k, &
                                      & l) + pi_infs(i)) + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
                        end do
                    end if

                    if (bubbles_euler) then
                        ! From prim: Compute nbub = (3/4pi) * \alpha / \bar{R^3}
                        do i = 1, nb
                            Rtmp(i) = q_prim_vf(qbmm_idx%rs(i))%sf(j, k, l)
                        end do

                        if (.not. qbmm) then
                            if (adv_n) then
                                q_cons_vf(eqn_idx%n)%sf(j, k, l) = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                                nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                            else
                                call s_comp_n_from_prim(real(q_prim_vf(eqn_idx%alf)%sf(j, k, l), kind=wp), Rtmp, nbub, weight)
                            end if
                        else
                            ! Initialize R3 averaging over R0 and R directions
                            R3tmp = 0._wp
                            do i = 1, nb
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) + sigR)**3._wp
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) - sigR)**3._wp
                            end do
                            ! Initialize nb
                            nbub = 3._wp*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/(4._wp*pi*R3tmp)
                        end if

                        do i = eqn_idx%bub%beg, eqn_idx%bub%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)*nbub
                        end do
                    end if

                    if (mhd) then
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        ! adding the elastic contribution Multiply \tau to \rho \tau
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        if (cont_damage) G = G*max((1._wp - q_prim_vf(eqn_idx%damage)%sf(j, k, l)), 0._wp)
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            ! Elastic energy addition (guard skips when G near zero from alpha undershoot)
                            if (G > verysmall) then
                                q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) + (q_prim_vf(i)%sf(j, k, &
                                          & l)**2._wp)/max(4._wp*G, verysmall)
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) + (q_prim_vf(i)%sf(j, k, &
                                              & l)**2._wp)/max(4._wp*G, verysmall)
                                end if
                            end if
                        end do
                    end if

                    if (surface_tension) then
                        q_cons_vf(eqn_idx%c)%sf(j, k, l) = q_prim_vf(eqn_idx%c)%sf(j, k, l)
                    end if

                    if (cont_damage) q_cons_vf(eqn_idx%damage)%sf(j, k, l) = q_prim_vf(eqn_idx%damage)%sf(j, k, l)

                    if (hyper_cleaning) q_cons_vf(eqn_idx%psi)%sf(j, k, l) = q_prim_vf(eqn_idx%psi)%sf(j, k, l)
                end do
            end do
        end do

    end subroutine s_convert_primitive_to_conservative_variables

    !> Convert primitive variables to Eulerian flux variables.
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, FK_vf, FK_src_vf, is1, is2, is3, s2b, s3b, dir_idx_in, &
        & dir_flg_in, hll_u_interface_in)

        integer, intent(in) :: s2b, s3b
        !> Working-direction mapping, passed explicitly: it is simulation state (m_global_parameters), and use-associating it into
        !! this common kernel spills registers on AMD OpenMP offload.
        integer, dimension(3), intent(in)                                                       :: dir_idx_in
        real(wp), dimension(3), intent(in)                                                      :: dir_flg_in
        logical, intent(in)                                                                     :: hll_u_interface_in
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in)                  :: qK_prim_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout)               :: FK_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,eqn_idx%adv%beg:), intent(inout) :: FK_src_vf
        type(int_bounds_info), intent(in)                                                       :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection variables, the specific heat ratio and liquid stiffness
        ! functions, the shear and volume Reynolds numbers and the Weber numbers

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_K
            real(wp), dimension(3)  :: alpha_K
            real(wp), dimension(3)  :: vel_K
            real(wp), dimension(10) :: Y_K
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_K
            real(wp), dimension(num_fluids)  :: alpha_K
            real(wp), dimension(num_vels)    :: vel_K
            real(wp), dimension(num_species) :: Y_K
        #:endif
        real(wp)               :: rho_K
        real(wp)               :: vel_K_sum
        real(wp)               :: pres_K
        real(wp)               :: E_K
        real(wp)               :: gamma_K
        real(wp)               :: pi_inf_K
        real(wp)               :: qv_K
        real(wp), dimension(2) :: Re_K
        real(wp)               :: G_K
        real(wp)               :: blkmod1_K, blkmod2_K, K_K
        real(wp)               :: T_K, mix_mol_weight, R_gas
        integer                :: i, j, k, l  !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        $:GPU_UPDATE(device='[is1b, is2b, is3b, is1e, is2e, is3e]')

        ! Computing the flux variables from the primitive variables, without accounting for the contribution of either viscosity or
        ! capillarity
        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_K, vel_K, alpha_K, Re_K, Y_K, rho_K, vel_K_sum, pres_K, E_K, gamma_K, &
                            & pi_inf_K, qv_K, G_K, blkmod1_K, blkmod2_K, K_K, T_K, mix_mol_weight, R_gas]', copyin='[dir_idx_in, &
                            & dir_flg_in, hll_u_interface_in]')
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, eqn_idx%cont%end
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                        alpha_K(i - eqn_idx%E) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K(i) = qK_prim_vf(j, k, l, eqn_idx%cont%end + i)
                    end do

                    vel_K_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K_sum = vel_K_sum + vel_K(i)**2._wp
                    end do

                    pres_K = qK_prim_vf(j, k, l, eqn_idx%E)
                    if (hypoelasticity) then
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K, G_K, Gs_vc)
                    else
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K)
                    end if

                    ! Computing the energy from the pressure

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            Y_K(i - eqn_idx%species%beg + 1) = qK_prim_vf(j, k, l, i)
                        end do
                        ! Computing the energy from the internal energy of the mixture
                        call get_mixture_molecular_weight(Y_k, mix_mol_weight)
                        R_gas = gas_constant/mix_mol_weight
                        T_K = pres_K/rho_K/R_gas
                        call get_mixture_energy_mass(T_K, Y_K, E_K)
                        E_K = rho_K*E_K + 5.e-1_wp*rho_K*vel_K_sum
                    else
                        ! Computing the energy from the pressure
                        E_K = gamma_K*pres_K + pi_inf_K + 5.e-1_wp*rho_K*vel_K_sum + qv_K
                    end if

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, eqn_idx%cont%end
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx_in(1))
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        FK_vf(j, k, l, &
                              & eqn_idx%cont%end + dir_idx_in(i)) = rho_K*vel_K(dir_idx_in(1))*vel_K(dir_idx_in(i)) &
                              & + pres_K*dir_flg_in(dir_idx_in(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, eqn_idx%E) = vel_K(dir_idx_in(1))*(E_K + pres_K)

                    ! Species advection Flux, \rho*u*Y
                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            FK_vf(j, k, l, i - 1 + eqn_idx%species%beg) = vel_K(dir_idx_in(1))*(rho_K*Y_K(i))
                        end do
                    end if

                    ! Match the volume-fraction flux representation exported by the Riemann solver. HLL Method 1: zero alpha
                    ! flux plus per-fluid interface-alpha source traces. Hypoelastic HLLD folds every non-conservative term
                    ! into its augmented flux (adv_src_mode_none), so its source trace is zero; for this cell-local conversion
                    ! the fold collapses exactly to -/+ K*u_n on the two volume-fraction rows (K = 0 without alt_soundspeed),
                    ! with the same two-fluid longitudinal-modulus K as the HLLD kernel (num_fluids = 2 is checker-enforced).
                    ! MHD HLLD keeps the per-fluid-trace representation it has always used. HLL Method 2, HLLC, and LF use the
                    ! shared-velocity representation below.
                    if (riemann_solver == riemann_solver_hlld) then
                        if (hypoelasticity) then
                            K_K = 0._wp
                            ! The fluid-2 subscripts must not be compiled when case optimization
                            ! bakes num_fluids = 1 (amdflang rejects them at compile time); the
                            ! checker prohibits hypoelastic HLLD there, so the block is dead code.
                            #:if not MFC_CASE_OPTIMIZATION or num_fluids > 1
                                if (alt_soundspeed) then
                                    blkmod1_K = ((gammas(1) + 1._wp)*pres_K + pi_infs(1))/gammas(1) + (4._wp/3._wp)*Gs_vc(1)
                                    blkmod2_K = ((gammas(2) + 1._wp)*pres_K + pi_infs(2))/gammas(2) + (4._wp/3._wp)*Gs_vc(2)
                                    K_K = alpha_K(1)*alpha_K(2)*(blkmod2_K - blkmod1_K)/(alpha_K(1)*blkmod2_K + alpha_K(2) &
                                                  & *blkmod1_K + verysmall)
                                end if
                            #:endif
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                FK_vf(j, k, l, i) = 0._wp
                                FK_src_vf(j, k, l, i) = 0._wp
                            end do
                            FK_vf(j, k, l, eqn_idx%adv%beg) = -K_K*vel_K(dir_idx_in(1))
                            FK_vf(j, k, l, eqn_idx%adv%end) = K_K*vel_K(dir_idx_in(1))
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                FK_vf(j, k, l, i) = 0._wp
                                FK_src_vf(j, k, l, i) = alpha_K(i - eqn_idx%E)
                            end do
                        end if
                    else if (riemann_solver == riemann_solver_hll .and. .not. hll_u_interface_in) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_vf(j, k, l, i) = 0._wp
                            FK_src_vf(j, k, l, i) = alpha_K(i - eqn_idx%E)
                        end do
                    else
                        ! Could be bubbles_euler!
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_vf(j, k, l, i) = vel_K(dir_idx_in(1))*alpha_K(i - eqn_idx%E)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_src_vf(j, k, l, i) = vel_K(dir_idx_in(1))
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_convert_primitive_to_flux_variables

    !> Compute partial densities and volume fractions
    subroutine s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        $:GPU_ROUTINE(function_name='s_compute_species_fraction', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: k, l, r
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(out) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(out) :: alpha_rho_K, alpha_K
        #:endif
        integer  :: i
        real(wp) :: alpha_K_sum

        if (num_fluids == 1) then
            alpha_rho_K(1) = q_vf(eqn_idx%cont%beg)%sf(k, l, r)
            if (igr .or. bubbles_euler) then
                alpha_K(1) = 1._wp
            else
                alpha_K(1) = q_vf(eqn_idx%adv%beg)%sf(k, l, r)
            end if
        else
            if (igr) then
                do i = 1, num_fluids - 1
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(eqn_idx%adv%beg + i - 1)%sf(k, l, r)
                end do
                alpha_rho_K(num_fluids) = q_vf(num_fluids)%sf(k, l, r)
                alpha_K(num_fluids) = 1._wp - sum(alpha_K(1:num_fluids - 1))
            else
                do i = 1, num_fluids
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(eqn_idx%adv%beg + i - 1)%sf(k, l, r)
                end do
            end if
        end if

        if (mpp_lim) then
            alpha_K_sum = 0._wp
            do i = 1, num_fluids
                alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                alpha_K_sum = alpha_K_sum + alpha_K(i)
            end do
            alpha_K = alpha_K/max(alpha_K_sum, 1.e-16_wp)
        end if

        if (num_fluids == 1 .and. bubbles_euler) alpha_K(1) = q_vf(eqn_idx%adv%beg)%sf(k, l, r)

    end subroutine s_compute_species_fraction

    !> Deallocate fluid property arrays and post-processing fields allocated during module initialization.
    impure subroutine s_finalize_variables_conversion_module()

        if (allocated(rho_sf)) deallocate (rho_sf, gamma_sf, pi_inf_sf)

        @:DEALLOCATE(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs_vc)
        if (allocated(bubrs_vc)) then
            @:DEALLOCATE(bubrs_vc)
        end if
        if (allocated(Res_vc)) then
            @:DEALLOCATE(Res_vc)
        end if

    end subroutine s_finalize_variables_conversion_module

    !> Accumulate stiffened-gas mixture coefficients over the first nf fluids.
    !!
    !! This is the only implementation of the stiffened-gas mixture rule. It is deliberately not merged
    !! with s_convert_species_to_mixture_variables_kernel, which additionally clips and renormalises
    !! alpha_K in place under mpp_lim, special-cases num_fluids == 1 with bubbles_euler, and optionally
    !! returns Re_K and G_K. Merging would need a clipping flag and optional dummies on a [seq] device
    !! routine, which is not portable across the offload backends.
    !!
    !! nf is not always num_fluids: the bubbles path in m_riemann_solver_hllc passes num_fluids - 1 to
    !! exclude the gas phase, and one site passes limited volume fractions rather than the raw ones.
    subroutine s_accumulate_mixture_properties(nf, alpha_rho_K, alpha_K, rho_K, gamma_K, pi_inf_K, qv_K)

        $:GPU_ROUTINE(function_name='s_accumulate_mixture_properties', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                 :: nf  !< Number of fluids to accumulate over
        real(wp), dimension(nf), intent(in) :: alpha_rho_K, alpha_K
        real(wp), intent(out)               :: rho_K, gamma_K, pi_inf_K, qv_K
        integer                             :: i   !< Loop iterator over fluids

        rho_K = 0._wp
        gamma_K = 0._wp
        pi_inf_K = 0._wp
        qv_K = 0._wp

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, nf
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            qv_K = qv_K + alpha_rho_K(i)*qvs(i)
        end do

    end subroutine s_accumulate_mixture_properties

    !> Total energy per unit volume of a stiffened-gas state.
    !!
    !! Thermodynamic contributions only. Magnetic energy (pres_mag) and elastic energy are added by the
    !! caller, because they are not equation-of-state terms: folding them in here would make the
    !! operator impossible to reuse for a second equation of state.
    !!
    !! The chemistry and relativistic branches of the Riemann solvers do not use this relation at all -
    !! chemistry builds E from the mixture internal energy and the relativistic form is unrelated - so
    !! those sites are deliberately left open-coded.
    subroutine s_compute_energy(pres, rho, gamma, pi_inf, qv, vel_sum, E)

        $:GPU_ROUTINE(function_name='s_compute_energy', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: pres, rho, gamma, pi_inf, qv, vel_sum
        real(wp), intent(out) :: E

        E = gamma*pres + pi_inf + 5.e-1_wp*rho*vel_sum + qv

    end subroutine s_compute_energy

    !> Compute the speed of sound of a thermodynamic state.
    !!
    !! The enthalpy is not an argument. For a real state it is fixed by the arguments that are, and it
    !! cancels out of the stiffened-gas branch entirely: substituting
    !! H = ((Gamma + 1)p + Pi + qv)/rho + |u|^2/2 into c^2 = (H - |u|^2/2 - qv/rho)/Gamma leaves
    !! c^2 = ((Gamma + 1)p + Pi)/(Gamma rho), free of H, |u|^2 and qv alike. Asking callers to supply
    !! those three is what let #1707 happen - five sites open-coded H and three of them dropped the qv
    !! this routine then subtracted - and here there is nothing to open-code, so it cannot recur.
    !!
    !! An average of two states is not a state: its enthalpy is a free input rather than a function of
    !! its pressure and density. Those callers use s_compute_speed_of_sound_avg.
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres, rho, gamma, pi_inf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(out) :: c
        real(wp)              :: blkmod1, blkmod2
        integer               :: q

        if (chemistry) then  ! Reacting mixture sound speed
            c = sqrt((1.0_wp + 1.0_wp/gamma)*pres/rho)
        else if (relativity) then  ! Relativistic sound speed, whose enthalpy is 1 + (Gamma + 1)p/rho
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/(1._wp + (gamma + 1._wp)*pres/rho))
        else
            if (alt_soundspeed) then  ! Wood's mixture sound speed via bulk moduli
                blkmod1 = ((gammas(1) + 1._wp)*pres + pi_infs(1))/gammas(1)
                blkmod2 = ((gammas(2) + 1._wp)*pres + pi_infs(2))/gammas(2)
                c = (1._wp/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
            else if (model_eqns == model_eqns_6eq) then  ! Six-equation model sound speed
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    c = c + adv(q)*gs_min(q)*(pres + pi_infs(q)/(gammas(q) + 1._wp))
                end do
                c = c/rho
            else if (model_eqns == model_eqns_5eq .and. bubbles_euler) then
                ! Sound speed for bubble mixture to order O(\alpha)

                if (mpp_lim .and. (num_fluids > 1)) then
                    c = (1._wp/gamma + 1._wp)*(pres + pi_inf/(gamma + 1._wp))/rho
                else
                    c = (1._wp/gamma + 1._wp)*(pres + pi_inf/(gamma + 1._wp))/(rho*(1._wp - adv(num_fluids)))
                end if
            else  ! Stiffened-gas mixture, with H, |u|^2 and qv cancelled out
                c = ((gamma + 1._wp)*pres + pi_inf)/(gamma*rho)
            end if

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if

    end subroutine s_compute_speed_of_sound

    !> Compute the speed of sound of an interface-averaged state.
    !!
    !! A Roe or arithmetic average of two states is not itself a state: the enthalpy it carries is an
    !! average of two enthalpies, not the one its own pressure and density imply, so the caller has to
    !! supply it - and with it |u|^2 and qv, which no longer cancel against it.
    !!
    !! Only the three branches below read an enthalpy. Every other equation of state is a function of
    !! the state alone, for which an average behaves exactly like a state, so those are deferred to
    !! s_compute_speed_of_sound rather than written out a second time. Keep the condition below in
    !! step with the branch list there.
    subroutine s_compute_speed_of_sound_avg(pres, rho, gamma, pi_inf, qv, vel_sum, H, c_c, adv, c)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres, rho, gamma, pi_inf, qv, vel_sum, H, c_c
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(out) :: c

        if (chemistry) then  ! Reacting mixture sound speed
            if (avg_state == avg_state_roe .and. abs(c_c) > verysmall) then
                c = sqrt(c_c - (gamma - 1.0_wp)*(vel_sum - H))
            else
                call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c)
            end if
        else if (relativity) then  ! Relativistic sound speed
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/H)
        else if (alt_soundspeed .or. model_eqns == model_eqns_6eq .or. (model_eqns == model_eqns_5eq .and. bubbles_euler)) then
            call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c)
        else  ! Stiffened-gas mixture, the one branch where the averaged enthalpy survives
            c = (H - 5.e-1*vel_sum - qv/rho)/gamma

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if

    end subroutine s_compute_speed_of_sound_avg

    !> Compute the fast magnetosonic wave speed from the sound speed, density, and magnetic field components.
    subroutine s_compute_fast_magnetosonic_speed(rho, c, B, norm, c_fast, h)

        $:GPU_ROUTINE(function_name='s_compute_fast_magnetosonic_speed', parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: B(3), rho, c
        real(wp), intent(in)  :: h  !< only used for relativity
        real(wp), intent(out) :: c_fast
        integer, intent(in)   :: norm
        real(wp)              :: B2, term, disc

        B2 = sum(B**2)

        if (.not. relativity) then
            term = c**2 + B2/rho
            disc = term**2 - 4*c**2*(B(norm)**2/rho)
        else
            ! Note: this is approximation for the non-relatisitic limit; accurate solution requires solving a quartic equation
            term = (c**2*(B(norm)**2 + rho*h) + B2)/(rho*h + B2)
            disc = term**2 - 4*c**2*B(norm)**2/(rho*h + B2)
        end if

#ifdef MFC_DEBUG
        if (disc < 0._wp) then
            print *, 'rho, c, Bx, By, Bz, h, term, disc:', rho, c, B(1), B(2), B(3), h, term, disc
            ! s_mpi_abort is a host routine and cannot be called from device code
            ! (this is a GPU routine); on GPU builds, emit the diagnostic print only.
#ifndef MFC_GPU
            call s_mpi_abort('Error: negative discriminant in s_compute_fast_magnetosonic_speed')
#endif
        end if
#endif

        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))

    end subroutine s_compute_fast_magnetosonic_speed

end module m_variables_conversion
