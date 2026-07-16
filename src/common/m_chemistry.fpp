!>
!! @file
!! @brief Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Multi-species chemistry interface for thermodynamic properties, reaction rates, and transport coefficients
module m_chemistry

    use m_thermochem, only: num_species, molecular_weights, get_temperature, get_net_production_rates, &
        & get_creation_destruction_rates, get_mole_fractions, get_species_binary_mass_diffusivities, &
        & get_species_mass_diffusivities_mixavg, gas_constant, get_mixture_molecular_weight, get_mixture_energy_mass, &
        & get_mixture_thermal_conductivity_mixavg, get_species_enthalpies_rt, get_mixture_viscosity_mixavg, &
        & get_mixture_specific_heat_cp_mass, get_mixture_enthalpy_mass

    use m_global_parameters

    implicit none

    type(int_bounds_info) :: isc1, isc2, isc3
    $:GPU_DECLARE(create='[isc1, isc2, isc3]')
    integer, dimension(3) :: offsets
    $:GPU_DECLARE(create='[offsets]')

contains

    !> Compute mixture viscosities for left and right states and invert them for use as reciprocal Reynolds numbers.
    subroutine compute_viscosity_and_inversion(T_L, Ys_L, T_R, Ys_R, Re_L, Re_R)

        $:GPU_ROUTINE(function_name='compute_viscosity_and_inversion',parallelism='[seq]', cray_inline=True)

        real(wp), intent(inout)                         :: T_L, T_R, Re_L, Re_R
        real(wp), dimension(num_species), intent(inout) :: Ys_R, Ys_L

        call get_mixture_viscosity_mixavg(T_L, Ys_L, Re_L)
        call get_mixture_viscosity_mixavg(T_R, Ys_R, Re_R)
        ! Convert dynamic viscosity to inverse (MFC stores 1/mu for Reynolds number convention)
        Re_L = 1.0_wp/Re_L
        Re_R = 1.0_wp/Re_R

    end subroutine compute_viscosity_and_inversion

    !> Initialize the temperature field from conservative variables by inverting the energy equation.
    subroutine s_compute_q_T_sf(q_T_sf, q_cons_vf, bounds)

        ! Initialize the temperature field at the start of the simulation to reasonable values. Temperature is computed the regular
        ! way using the conservative variables.

        type(scalar_field), intent(inout)                   :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in)   :: bounds
        integer                                             :: x, y, z, eqn
        real(wp)                                            :: energy, T_in
        real(wp), dimension(num_species)                    :: Ys

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    do eqn = eqn_idx%species%beg, eqn_idx%species%end
                        Ys(eqn - eqn_idx%species%beg + 1) = q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z)
                    end do

                    ! e = E - 1/2*|u|^2 cons. eqn_idx%E = \rho E cons. eqn_idx%cont%beg = \rho (1-fluid model) cons. eqn_idx%mom%beg
                    ! + i = \rho u_i
                    energy = q_cons_vf(eqn_idx%E)%sf(x, y, z)/q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z)
                    do eqn = eqn_idx%mom%beg, eqn_idx%mom%end
                        energy = energy - 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z))**2._wp
                    end do

                    T_in = real(q_T_sf%sf(x, y, z), kind=wp)
                    call get_temperature(energy, dflt_T_guess, Ys, .true., T_in)
                    q_T_sf%sf(x, y, z) = T_in
                end do
            end do
        end do

    end subroutine s_compute_q_T_sf

    !> Compute the temperature field from primitive variables using the ideal gas law and mixture molecular weight.
    subroutine s_compute_T_from_primitives(q_T_sf, q_prim_vf, bounds)

        type(scalar_field), intent(inout)                   :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in)   :: bounds
        integer                                             :: x, y, z, i
        real(wp), dimension(num_species)                    :: Ys
        real(wp)                                            :: mix_mol_weight

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        Ys(i - eqn_idx%species%beg + 1) = q_prim_vf(i)%sf(x, y, z)
                    end do

                    call get_mixture_molecular_weight(Ys, mix_mol_weight)
                    q_T_sf%sf(x, y, z) = q_prim_vf(eqn_idx%E)%sf(x, y, z)*mix_mol_weight/(gas_constant*q_prim_vf(1)%sf(x, y, z))
                end do
            end do
        end do

    end subroutine s_compute_T_from_primitives

    !> Add chemical reaction source terms to the species transport RHS using net production rates.
    subroutine s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp, q_T_sf, q_prim_qp, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(inout)                      :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_qp, q_prim_qp
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: x, y, z
        integer                                                :: eqn
        real(wp)                                               :: T
        real(wp)                                               :: rho, omega_m

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10) :: Ys
            real(wp), dimension(10) :: omega
        #:else
            real(wp), dimension(num_species) :: Ys
            real(wp), dimension(num_species) :: omega
        #:endif

        $:GPU_PARALLEL_LOOP(collapse=3, private='[Ys, omega, eqn, T, rho, omega_m]', copyin='[bounds]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = eqn_idx%species%beg, eqn_idx%species%end
                        Ys(eqn - eqn_idx%species%beg + 1) = q_prim_qp(eqn)%sf(x, y, z)
                    end do

                    rho = q_cons_qp(eqn_idx%cont%end)%sf(x, y, z)
                    T = q_T_sf%sf(x, y, z)

                    call get_net_production_rates(rho, T, Ys, omega)

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = eqn_idx%species%beg, eqn_idx%species%end
                        omega_m = molecular_weights(eqn - eqn_idx%species%beg + 1)*omega(eqn - eqn_idx%species%beg + 1)
                        rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_chemistry_reaction_flux

    !> Operator-split integration of the reaction source with an alpha-QSS (Mott quasi-steady-state) reactor. Called after the flow
    !! update: each cell's constant-(rho, e) reactor is advanced over dtime with chem_params%reaction_substeps predictor-corrector
    !! sub-steps, updating the species partial densities and temperature in place. Mixture density, momentum, and total energy are
    !! unchanged (reactions convert chemical to thermal energy at fixed internal energy). The alpha-QSS update treats each species'
    !! destruction as a pseudo-first-order loss, so it is stable for stiff ignition where an explicit source would overshoot and
    !! diverge, and relaxes to the correct chemical equilibrium rather than over-heating.
    subroutine s_chemistry_reaction_substep(q_cons_vf, q_T_sf, dtime, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), intent(inout)                      :: q_T_sf
        real(wp), intent(in)                                   :: dtime
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: x, y, z, eqn, s, nsub
        real(wp)                                               :: rho, energy, T, T_new, dt_sub, Ysum
        real(wp)                                               :: r, r2, wr, loss_i, prod_p, loss_p, Lbar, pbar
        real(wp)                                               :: stiff_max, cell_stiff
        real(wp), parameter                                    :: y_floor = 1.e-16_wp
        real(wp), parameter                                    :: stiff_target = 0.5_wp

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10) :: Ys, cdot, ddot, y0, prod0, Lloss, alp
        #:else
            real(wp), dimension(num_species) :: Ys, cdot, ddot, y0, prod0, Lloss, alp
        #:endif

        if (chem_params%adap_substeps) then
            ! Pass 1: per-rank local stiffness probe -> adapt nsub for this step, no MPI. Each rank
            ! sizes its own work from the largest fractional net species change any of its cells sees.
            stiff_max = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[Ys, cdot, eqn, rho, T, wr, cell_stiff]', reduction='[[stiff_max]]', &
                                & reductionOp='[MAX]', copyin='[bounds, dtime]')
            do z = bounds(3)%beg, bounds(3)%end
                do y = bounds(2)%beg, bounds(2)%end
                    do x = bounds(1)%beg, bounds(1)%end
                        rho = q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z)
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = eqn_idx%species%beg, eqn_idx%species%end
                            Ys(eqn - eqn_idx%species%beg + 1) = q_cons_vf(eqn)%sf(x, y, z)/rho
                        end do
                        T = q_T_sf%sf(x, y, z)
                        ! Net rate (creation - destruction) on purpose: nsub sizes the accuracy of the
                        ! composition trajectory, which is set by how fast Ys actually moves, not by the
                        ! raw forward/reverse magnitudes. In fast partial equilibrium the net is ~0 and the
                        ! composition is static, so the floor is adequate; the alpha-QSS update is itself
                        ! unconditionally stable, so an under-sized nsub loses accuracy, never stability.
                        call get_net_production_rates(rho, T, Ys, cdot)  ! net omega in cdot
                        cell_stiff = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = 1, num_species
                            wr = molecular_weights(eqn)/rho
                            cell_stiff = max(cell_stiff, dtime*abs(wr*cdot(eqn))/max(Ys(eqn), y_floor))
                        end do
                        stiff_max = max(stiff_max, cell_stiff)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            nsub = ceiling(max(real(chem_params%reaction_substeps, wp), min(real(chem_params%reaction_substeps_max, wp), &
                           & stiff_max/stiff_target)))
        else
            nsub = chem_params%reaction_substeps
        end if
        dt_sub = dtime/real(nsub, wp)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[Ys, cdot, ddot, y0, prod0, Lloss, alp, eqn, s, rho, energy, T, T_new, Ysum, r, &
                            & r2, wr, loss_i, prod_p, loss_p, Lbar, pbar]', copyin='[bounds, dt_sub, nsub]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    rho = q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z)

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = eqn_idx%species%beg, eqn_idx%species%end
                        Ys(eqn - eqn_idx%species%beg + 1) = q_cons_vf(eqn)%sf(x, y, z)/rho
                    end do

                    ! internal energy per mass, held fixed through the reactor sub-steps
                    energy = q_cons_vf(eqn_idx%E)%sf(x, y, z)/rho
                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = eqn_idx%mom%beg, eqn_idx%mom%end
                        energy = energy - 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/rho)**2
                    end do

                    T = q_T_sf%sf(x, y, z)

                    do s = 1, nsub
                        ! predictor: rates at the start of the sub-step (one fused pass fills both)
                        call get_creation_destruction_rates(rho, T, Ys, cdot, ddot)
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = 1, num_species
                            y0(eqn) = Ys(eqn)
                            wr = molecular_weights(eqn)/rho
                            prod0(eqn) = wr*cdot(eqn)  ! mass-fraction production
                            loss_i = wr*ddot(eqn)  ! mass-fraction loss
                            Lloss(eqn) = loss_i/max(Ys(eqn), y_floor)  ! pseudo-first-order loss rate
                            r = dt_sub*Lloss(eqn); r2 = r*r
                            alp(eqn) = (180._wp + 60._wp*r + 11._wp*r2 + r2*r)/(360._wp + 60._wp*r + 12._wp*r2 + r2*r)
                            Ys(eqn) = y0(eqn) + dt_sub*(prod0(eqn) - loss_i)/(1._wp + alp(eqn)*dt_sub*Lloss(eqn))
                            if (Ys(eqn) < 0._wp) Ys(eqn) = 0._wp
                        end do
                        ! corrector: re-evaluate rates at the predicted state (T is the Newton guess; T_new is intent(out))
                        call get_temperature(energy, T, Ys, .true., T_new)
                        call get_creation_destruction_rates(rho, T_new, Ys, cdot, ddot)
                        Ysum = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = 1, num_species
                            wr = molecular_weights(eqn)/rho
                            prod_p = wr*cdot(eqn)
                            loss_p = wr*ddot(eqn)
                            Lbar = 0.5_wp*(Lloss(eqn) + loss_p/max(Ys(eqn), y_floor))
                            pbar = alp(eqn)*prod_p + (1._wp - alp(eqn))*prod0(eqn)
                            Ys(eqn) = y0(eqn) + dt_sub*(pbar - Lbar*y0(eqn))/(1._wp + alp(eqn)*dt_sub*Lbar)
                            if (Ys(eqn) < 0._wp) Ys(eqn) = 0._wp
                            Ysum = Ysum + Ys(eqn)
                        end do
                        if (Ysum > y_floor) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = 1, num_species
                                Ys(eqn) = Ys(eqn)/Ysum
                            end do
                        else
                            ! Degenerate corrector (every species clipped to zero): fall back to the
                            ! sub-step's starting composition rather than dividing by a vanishing sum.
                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = 1, num_species
                                Ys(eqn) = y0(eqn)
                            end do
                        end if
                        call get_temperature(energy, T, Ys, .true., T_new)
                        T = T_new
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = eqn_idx%species%beg, eqn_idx%species%end
                        q_cons_vf(eqn)%sf(x, y, z) = rho*Ys(eqn - eqn_idx%species%beg + 1)
                    end do
                    q_T_sf%sf(x, y, z) = T
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_chemistry_reaction_substep

    !> Compute species mass diffusion fluxes at cell interfaces using mixture-averaged diffusivities.
    subroutine s_compute_chemistry_diffusion_flux(idir, q_prim_qp, flux_src_vf, irx, iry, irz, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_qp
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in)                      :: irx, iry, irz
        integer, intent(in)                                    :: idir
        type(scalar_field), intent(in)                         :: q_T_sf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10) :: Xs_L, Xs_R, Xs_cell, Ys_L, Ys_R, Ys_cell
            real(wp), dimension(10) :: mass_diffusivities_mixavg1, mass_diffusivities_mixavg2
            real(wp), dimension(10) :: mass_diffusivities_mixavg_Cell, dXk_dxi, h_l, h_r, h_k
            real(wp), dimension(10) :: Mass_Diffu_Flux, dYk_dxi
        #:else
            real(wp), dimension(num_species) :: Xs_L, Xs_R, Xs_cell, Ys_L, Ys_R, Ys_cell
            real(wp), dimension(num_species) :: mass_diffusivities_mixavg1, mass_diffusivities_mixavg2
            real(wp), dimension(num_species) :: mass_diffusivities_mixavg_Cell, dXk_dxi, h_l, h_r, h_k
            real(wp), dimension(num_species) :: Mass_Diffu_Flux, dYk_dxi
        #:endif

        real(wp)              :: Mass_Diffu_Energy
        real(wp)              :: MW_L, MW_R, MW_cell, T_L, T_R, P_L, P_R, rho_L, rho_R, rho_cell, rho_Vic
        real(wp)              :: lambda_L, lambda_R, lambda_Cell, dT_dxi, grid_spacing
        real(wp)              :: Cp_L, Cp_R
        real(wp)              :: diffusivity_L, diffusivity_R, diffusivity_cell
        real(wp)              :: hmix_L, hmix_R, dh_dxi
        integer               :: x, y, z, i, n, eqn
        integer, dimension(3) :: offsets

        isc1 = irx; isc2 = iry; isc3 = irz

        $:GPU_UPDATE(device='[isc1, isc2, isc3]')

        if (chemistry) then
            ! Set offsets based on direction using array indexing
            offsets = 0
            offsets(idir) = 1
            ! Model 1: Mixture-Average Transport
            if (chem_params%transport_model == 1) then
                ! Note: Added 'i' and 'eqn' to private list.
                $:GPU_PARALLEL_LOOP(collapse=3,  private='[x, y, z, i, eqn, Ys_L, Ys_R, Ys_cell, Xs_L, Xs_R, &
                                    & mass_diffusivities_mixavg1, mass_diffusivities_mixavg2, mass_diffusivities_mixavg_Cell, &
                                    & h_l, h_r, Xs_cell, h_k, dXk_dxi, Mass_Diffu_Flux, Mass_Diffu_Energy, MW_L, MW_R, MW_cell, &
                                    & T_L, T_R, P_L, P_R, rho_L, rho_R, rho_cell, rho_Vic, lambda_L, lambda_R, lambda_Cell, &
                                    & dT_dxi, grid_spacing]', copyin='[offsets]')
                do z = isc3%beg, isc3%end
                    do y = isc2%beg, isc2%end
                        do x = isc1%beg, isc1%end
                            ! Calculate grid spacing using direction-based indexing
                            select case (idir)
                            case (1)
                                grid_spacing = x_cc(x + 1) - x_cc(x)
                            case (2)
                                grid_spacing = y_cc(y + 1) - y_cc(y)
                            case (3)
                                grid_spacing = z_cc(z + 1) - z_cc(z)
                            end select

                            ! Extract species mass fractions
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                Ys_L(i - eqn_idx%species%beg + 1) = q_prim_qp(i)%sf(x, y, z)
                                Ys_R(i - eqn_idx%species%beg + 1) = q_prim_qp(i)%sf(x + offsets(1), y + offsets(2), z + offsets(3))
                                Ys_cell(i - eqn_idx%species%beg + 1) = 0.5_wp*(Ys_L(i - eqn_idx%species%beg + 1) + Ys_R(i &
                                        & - eqn_idx%species%beg + 1))
                            end do

                            ! Calculate molecular weights and mole fractions
                            call get_mixture_molecular_weight(Ys_L, MW_L)
                            call get_mixture_molecular_weight(Ys_R, MW_R)
                            MW_cell = 0.5_wp*(MW_L + MW_R)

                            call get_mole_fractions(MW_L, Ys_L, Xs_L)
                            call get_mole_fractions(MW_R, Ys_R, Xs_R)

                            P_L = q_prim_qp(eqn_idx%E)%sf(x, y, z)
                            P_R = q_prim_qp(eqn_idx%E)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            rho_L = q_prim_qp(1)%sf(x, y, z)
                            rho_R = q_prim_qp(1)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            T_L = q_T_sf%sf(x, y, z)
                            T_R = q_T_sf%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            rho_cell = 0.5_wp*(rho_L + rho_R)
                            dT_dxi = (T_R - T_L)/grid_spacing

                            ! Get transport properties
                            call get_species_mass_diffusivities_mixavg(P_L, T_L, Ys_L, mass_diffusivities_mixavg1)
                            call get_species_mass_diffusivities_mixavg(P_R, T_R, Ys_R, mass_diffusivities_mixavg2)

                            call get_mixture_thermal_conductivity_mixavg(T_L, Ys_L, lambda_L)
                            call get_mixture_thermal_conductivity_mixavg(T_R, Ys_R, lambda_R)

                            call get_species_enthalpies_rt(T_L, h_l)
                            call get_species_enthalpies_rt(T_R, h_r)

                            ! Calculate species properties and gradients
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                h_l(i - eqn_idx%species%beg + 1) = h_l(i - eqn_idx%species%beg + 1) &
                                    & *gas_constant*T_L/molecular_weights(i - eqn_idx%species%beg + 1)
                                h_r(i - eqn_idx%species%beg + 1) = h_r(i - eqn_idx%species%beg + 1) &
                                    & *gas_constant*T_R/molecular_weights(i - eqn_idx%species%beg + 1)
                                Xs_cell(i - eqn_idx%species%beg + 1) = 0.5_wp*(Xs_L(i - eqn_idx%species%beg + 1) + Xs_R(i &
                                        & - eqn_idx%species%beg + 1))
                                h_k(i - eqn_idx%species%beg + 1) = 0.5_wp*(h_l(i - eqn_idx%species%beg + 1) + h_r(i &
                                    & - eqn_idx%species%beg + 1))
                                dXk_dxi(i - eqn_idx%species%beg + 1) = (Xs_R(i - eqn_idx%species%beg + 1) - Xs_L(i &
                                        & - eqn_idx%species%beg + 1))/grid_spacing
                            end do

                            ! Calculate mixture-averaged diffusivities
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                mass_diffusivities_mixavg_Cell(i - eqn_idx%species%beg + 1) = (mass_diffusivities_mixavg2(i &
                                                               & - eqn_idx%species%beg + 1) + mass_diffusivities_mixavg1(i &
                                                               & - eqn_idx%species%beg + 1))/2.0_wp
                            end do

                            lambda_Cell = 0.5_wp*(lambda_R + lambda_L)

                            ! Calculate mass diffusion fluxes
                            rho_Vic = 0.0_wp
                            Mass_Diffu_Energy = 0.0_wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = eqn_idx%species%beg, eqn_idx%species%end
                                Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1) = rho_cell*mass_diffusivities_mixavg_Cell(eqn &
                                                & - eqn_idx%species%beg + 1)*molecular_weights(eqn - eqn_idx%species%beg + 1) &
                                                & /MW_cell*dXk_dxi(eqn - eqn_idx%species%beg + 1)
                                rho_Vic = rho_Vic + Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1)
                                Mass_Diffu_Energy = Mass_Diffu_Energy + h_k(eqn - eqn_idx%species%beg + 1)*Mass_Diffu_Flux(eqn &
                                    & - eqn_idx%species%beg + 1)
                            end do

                            ! Apply corrections for mass conservation
                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = eqn_idx%species%beg, eqn_idx%species%end
                                Mass_Diffu_Energy = Mass_Diffu_Energy - h_k(eqn - eqn_idx%species%beg + 1)*Ys_cell(eqn &
                                    & - eqn_idx%species%beg + 1)*rho_Vic
                                Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1) = Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1) &
                                                & - rho_Vic*Ys_cell(eqn - eqn_idx%species%beg + 1)
                            end do

                            ! Add thermal conduction contribution
                            Mass_Diffu_Energy = lambda_Cell*dT_dxi + Mass_Diffu_Energy

                            ! Update flux arrays
                            flux_src_vf(eqn_idx%E)%sf(x, y, z) = flux_src_vf(eqn_idx%E)%sf(x, y, z) - Mass_Diffu_Energy

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = eqn_idx%species%beg, eqn_idx%species%end
                                flux_src_vf(eqn)%sf(x, y, z) = flux_src_vf(eqn)%sf(x, y, &
                                            & z) - Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                ! Model 2: Unity Lewis Number
            else if (chem_params%transport_model == 2) then
                ! Note: Added ALL scalars and 'i'/'eqn' to private list to prevent race conditions.
                $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, i, eqn, Ys_L, Ys_R, Ys_cell, dYk_dxi, Mass_Diffu_Flux, &
                                    & grid_spacing, MW_L, MW_R, MW_cell, P_L, P_R, rho_L, rho_R, rho_cell, T_L, T_R, Cp_L, Cp_R, &
                                    & hmix_L, hmix_R, dh_dxi, lambda_L, lambda_R, lambda_Cell, diffusivity_L, diffusivity_R, &
                                    & diffusivity_cell, Mass_Diffu_Energy]', copyin='[offsets]')
                do z = isc3%beg, isc3%end
                    do y = isc2%beg, isc2%end
                        do x = isc1%beg, isc1%end
                            ! Calculate grid spacing using direction-based indexing
                            select case (idir)
                            case (1)
                                grid_spacing = x_cc(x + 1) - x_cc(x)
                            case (2)
                                grid_spacing = y_cc(y + 1) - y_cc(y)
                            case (3)
                                grid_spacing = z_cc(z + 1) - z_cc(z)
                            end select

                            ! Extract species mass fractions
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                Ys_L(i - eqn_idx%species%beg + 1) = q_prim_qp(i)%sf(x, y, z)
                                Ys_R(i - eqn_idx%species%beg + 1) = q_prim_qp(i)%sf(x + offsets(1), y + offsets(2), z + offsets(3))
                                Ys_cell(i - eqn_idx%species%beg + 1) = 0.5_wp*(Ys_L(i - eqn_idx%species%beg + 1) + Ys_R(i &
                                        & - eqn_idx%species%beg + 1))
                            end do

                            ! Calculate molecular weights and mole fractions
                            call get_mixture_molecular_weight(Ys_L, MW_L)
                            call get_mixture_molecular_weight(Ys_R, MW_R)
                            MW_cell = 0.5_wp*(MW_L + MW_R)

                            P_L = q_prim_qp(eqn_idx%E)%sf(x, y, z)
                            P_R = q_prim_qp(eqn_idx%E)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            rho_L = q_prim_qp(1)%sf(x, y, z)
                            rho_R = q_prim_qp(1)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            T_L = q_T_sf%sf(x, y, z)
                            T_R = q_T_sf%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            rho_cell = 0.5_wp*(rho_L + rho_R)

                            call get_mixture_specific_heat_cp_mass(T_L, Ys_L, Cp_L)
                            call get_mixture_specific_heat_cp_mass(T_R, Ys_R, Cp_R)
                            call get_mixture_enthalpy_mass(T_L, Ys_L, hmix_L)
                            call get_mixture_enthalpy_mass(T_R, Ys_R, hmix_R)
                            dh_dxi = (hmix_R - hmix_L)/grid_spacing

                            ! Get transport properties
                            call get_mixture_thermal_conductivity_mixavg(T_L, Ys_L, lambda_L)
                            call get_mixture_thermal_conductivity_mixavg(T_R, Ys_R, lambda_R)

                            ! Calculate species properties and gradients
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%species%beg, eqn_idx%species%end
                                dYk_dxi(i - eqn_idx%species%beg + 1) = (Ys_R(i - eqn_idx%species%beg + 1) - Ys_L(i &
                                        & - eqn_idx%species%beg + 1))/grid_spacing
                            end do

                            ! Calculate mixture-averaged diffusivities
                            diffusivity_L = lambda_L/rho_L/Cp_L
                            diffusivity_R = lambda_R/rho_R/Cp_R

                            lambda_Cell = 0.5_wp*(lambda_R + lambda_L)
                            diffusivity_cell = 0.5_wp*(diffusivity_R + diffusivity_L)

                            ! Calculate mass diffusion fluxes
                            Mass_Diffu_Energy = 0.0_wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = eqn_idx%species%beg, eqn_idx%species%end
                                Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1) = rho_cell*diffusivity_cell*dYk_dxi(eqn &
                                                & - eqn_idx%species%beg + 1)
                            end do
                            Mass_Diffu_Energy = rho_cell*diffusivity_cell*dh_dxi

                            ! Update flux arrays
                            flux_src_vf(eqn_idx%E)%sf(x, y, z) = flux_src_vf(eqn_idx%E)%sf(x, y, z) - Mass_Diffu_Energy

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = eqn_idx%species%beg, eqn_idx%species%end
                                flux_src_vf(eqn)%sf(x, y, z) = flux_src_vf(eqn)%sf(x, y, &
                                            & z) - Mass_Diffu_Flux(eqn - eqn_idx%species%beg + 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_compute_chemistry_diffusion_flux

end module m_chemistry
