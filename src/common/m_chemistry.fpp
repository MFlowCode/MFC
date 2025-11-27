!!>
!! @file   m_chemistry.f90
!! @brief  Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

module m_chemistry

    use m_thermochem, only: &
        num_species, molecular_weights, get_temperature, get_net_production_rates, &
        get_mole_fractions, get_species_binary_mass_diffusivities, &
        get_species_mass_diffusivities_mixavg, gas_constant, get_mixture_molecular_weight, &
        get_mixture_energy_mass, get_mixture_thermal_conductivity_mixavg, get_species_enthalpies_rt, &
        get_mixture_viscosity_mixavg

    use m_global_parameters

    implicit none

    #:block DEF_AMD
        real(wp) :: molecular_weights_nonparameter(10) = &
                    (/2.016, 1.008, 15.999, 31.998, 17.007, 18.015, 33.006, &
                      34.014, 39.95, 28.014/)
        $:GPU_DECLARE(create='[molecular_weights_nonparameter]')
    #:endblock DEF_AMD

    type(int_bounds_info) :: isc1, isc2, isc3
    $:GPU_DECLARE(create='[isc1, isc2, isc3]')
    integer, dimension(3) :: offsets
    $:GPU_DECLARE(create='[offsets]')

contains

    subroutine compute_viscosity_and_inversion(T_L, Ys_L, T_R, Ys_R, Re_L, Re_R)

        $:GPU_ROUTINE(function_name='compute_viscosity_and_inversion',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(inout) :: T_L, T_R, Re_L, Re_R
        real(wp), dimension(num_species), intent(inout) :: Ys_R, Ys_L

        call get_mixture_viscosity_mixavg(T_L, Ys_L, Re_L)
        call get_mixture_viscosity_mixavg(T_R, Ys_R, Re_R)
        Re_L = 1.0_wp/Re_L
        Re_R = 1.0_wp/Re_R

    end subroutine compute_viscosity_and_inversion

    subroutine s_compute_q_T_sf(q_T_sf, q_cons_vf, bounds)

        ! Initialize the temperature field at the start of the simulation to
        ! reasonable values. Temperature is computed the regular way using the
        ! conservative variables.

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, eqn
        real(wp) :: energy, T_in
        real(wp), dimension(num_species) :: Ys

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = &
                            q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    end do

                    ! e = E - 1/2*|u|^2
                    ! cons. E_idx     = \rho E
                    ! cons. contxb    = \rho         (1-fluid model)
                    ! cons. momxb + i = \rho u_i
                    energy = q_cons_vf(E_idx)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    do eqn = momxb, momxe
                        energy = energy - &
                                 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z))**2._wp
                    end do

                    T_in = real(q_T_sf%sf(x, y, z), kind=wp)
                    call get_temperature(energy, dflt_T_guess, Ys, .true., T_in)
                    q_T_sf%sf(x, y, z) = T_in

                end do
            end do
        end do

    end subroutine s_compute_q_T_sf

    subroutine s_compute_T_from_primitives(q_T_sf, q_prim_vf, bounds)

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, i
        real(wp), dimension(num_species) :: Ys
        real(wp) :: mix_mol_weight

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    do i = chemxb, chemxe
                        Ys(i - chemxb + 1) = q_prim_vf(i)%sf(x, y, z)
                    end do

                    call get_mixture_molecular_weight(Ys, mix_mol_weight)
                    q_T_sf%sf(x, y, z) = q_prim_vf(E_idx)%sf(x, y, z)*mix_mol_weight/(gas_constant*q_prim_vf(1)%sf(x, y, z))
                end do
            end do
        end do

    end subroutine s_compute_T_from_primitives

    subroutine s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp, q_T_sf, q_prim_qp, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_qp, q_prim_qp
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z
        integer :: eqn
        real(wp) :: T
        real(wp) :: rho, omega_m
        real(wp), dimension(num_species) :: Ys
        real(wp), dimension(num_species) :: omega

        $:GPU_PARALLEL_LOOP(collapse=3, private='[Ys, omega, eqn, T, rho, omega, omega_m]', copyin='[bounds]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                    end do

                    rho = q_cons_qp(contxe)%sf(x, y, z)
                    T = q_T_sf%sf(x, y, z)

                    call get_net_production_rates(rho, T, Ys, omega)

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = chemxb, chemxe
                        #:block UNDEF_AMD
                            omega_m = molecular_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)
                        #:endblock UNDEF_AMD
                        #:block DEF_AMD
                            omega_m = molecular_weights_nonparameter(eqn - chemxb + 1)*omega(eqn - chemxb + 1)
                        #:endblock DEF_AMD
                        rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                    end do

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_chemistry_reaction_flux

    subroutine s_compute_chemistry_diffusion_flux(idir, q_prim_qp, flux_src_vf, irx, iry, irz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_qp
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz

        integer, intent(in) :: idir

        real(wp), dimension(num_species) :: Xs_L, Xs_R, Xs_cell, Ys_L, Ys_R, Ys_cell
        real(wp), dimension(num_species) :: mass_diffusivities_mixavg1, mass_diffusivities_mixavg2
        real(wp), dimension(num_species) :: mass_diffusivities_mixavg_Cell, dXk_dxi, h_l, h_r, h_k
        real(wp), dimension(num_species) :: Mass_Diffu_Flux
        real(wp) :: Mass_Diffu_Energy
        real(wp) :: MW_L, MW_R, MW_cell, Rgas_L, Rgas_R, T_L, T_R, P_L, P_R, rho_L, rho_R, rho_cell, rho_Vic
        real(wp) :: lambda_L, lambda_R, lambda_Cell, dT_dxi, grid_spacing

        integer :: x, y, z, i, n, eqn
        integer, dimension(3) :: offsets

        isc1 = irx; isc2 = iry; isc3 = irz

        $:GPU_UPDATE(device='[isc1,isc2,isc3]')

        if (chemistry) then
            ! Set offsets based on direction using array indexing
            offsets = 0
            offsets(idir) = 1
            #:block UNDEF_AMD
                $:GPU_PARALLEL_LOOP(collapse=3,  private='[x,y,z,Ys_L, Ys_R, Ys_cell, Xs_L, Xs_R, mass_diffusivities_mixavg1, mass_diffusivities_mixavg2, mass_diffusivities_mixavg_Cell, h_l, h_r, Xs_cell, h_k, dXk_dxi,Mass_Diffu_Flux, Mass_Diffu_Energy, MW_L, MW_R, MW_cell, Rgas_L, Rgas_R, T_L, T_R, P_L, P_R, rho_L, rho_R, rho_cell, rho_Vic, lambda_L, lambda_R, lambda_Cell, dT_dxi, grid_spacing]', copyin='[offsets]')
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
                            do i = chemxb, chemxe
                                Ys_L(i - chemxb + 1) = q_prim_qp(i)%sf(x, y, z)
                                Ys_R(i - chemxb + 1) = q_prim_qp(i)%sf(x + offsets(1), y + offsets(2), z + offsets(3))
                                Ys_cell(i - chemxb + 1) = 0.5_wp*(Ys_L(i - chemxb + 1) + Ys_R(i - chemxb + 1))
                            end do

                            ! Calculate molecular weights and mole fractions
                            call get_mixture_molecular_weight(Ys_L, MW_L)
                            call get_mixture_molecular_weight(Ys_R, MW_R)
                            MW_cell = 0.5_wp*(MW_L + MW_R)

                            call get_mole_fractions(MW_L, Ys_L, Xs_L)
                            call get_mole_fractions(MW_R, Ys_R, Xs_R)

                            ! Calculate gas constants and thermodynamic properties
                            Rgas_L = gas_constant/MW_L
                            Rgas_R = gas_constant/MW_R

                            P_L = q_prim_qp(E_idx)%sf(x, y, z)
                            P_R = q_prim_qp(E_idx)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            rho_L = q_prim_qp(1)%sf(x, y, z)
                            rho_R = q_prim_qp(1)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                            T_L = P_L/rho_L/Rgas_L
                            T_R = P_R/rho_R/Rgas_R

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
                            do i = chemxb, chemxe
                                h_l(i - chemxb + 1) = h_l(i - chemxb + 1)*gas_constant*T_L/molecular_weights(i - chemxb + 1)
                                h_r(i - chemxb + 1) = h_r(i - chemxb + 1)*gas_constant*T_R/molecular_weights(i - chemxb + 1)
                                Xs_cell(i - chemxb + 1) = 0.5_wp*(Xs_L(i - chemxb + 1) + Xs_R(i - chemxb + 1))
                                h_k(i - chemxb + 1) = 0.5_wp*(h_l(i - chemxb + 1) + h_r(i - chemxb + 1))
                                dXk_dxi(i - chemxb + 1) = (Xs_R(i - chemxb + 1) - Xs_L(i - chemxb + 1))/grid_spacing
                            end do

                            ! Calculate mixture-averaged diffusivities
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = chemxb, chemxe
                                mass_diffusivities_mixavg_Cell(i - chemxb + 1) = &
                                    (mass_diffusivities_mixavg2(i - chemxb + 1) + mass_diffusivities_mixavg1(i - chemxb + 1))/2.0_wp
                            end do

                            lambda_Cell = 0.5_wp*(lambda_R + lambda_L)

                            ! Calculate mass diffusion fluxes
                            rho_Vic = 0.0_wp
                            Mass_Diffu_Energy = 0.0_wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = chemxb, chemxe
                                Mass_Diffu_Flux(eqn - chemxb + 1) = rho_cell*mass_diffusivities_mixavg_Cell(eqn - chemxb + 1)* &
                                                                    molecular_weights(eqn - chemxb + 1)/MW_cell*dXk_dxi(eqn - chemxb + 1)
                                rho_Vic = rho_Vic + Mass_Diffu_Flux(eqn - chemxb + 1)
                                Mass_Diffu_Energy = Mass_Diffu_Energy + h_k(eqn - chemxb + 1)*Mass_Diffu_Flux(eqn - chemxb + 1)
                            end do

                            ! Apply corrections for mass conservation
                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = chemxb, chemxe
                                Mass_Diffu_Energy = Mass_Diffu_Energy - h_k(eqn - chemxb + 1)*Ys_cell(eqn - chemxb + 1)*rho_Vic
                                Mass_Diffu_Flux(eqn - chemxb + 1) = Mass_Diffu_Flux(eqn - chemxb + 1) - rho_Vic*Ys_cell(eqn - chemxb + 1)
                            end do

                            ! Add thermal conduction contribution
                            Mass_Diffu_Energy = lambda_Cell*dT_dxi + Mass_Diffu_Energy

                            ! Update flux arrays
                            flux_src_vf(E_idx)%sf(x, y, z) = flux_src_vf(E_idx)%sf(x, y, z) - Mass_Diffu_Energy

                            $:GPU_LOOP(parallelism='[seq]')
                            do eqn = chemxb, chemxe
                                flux_src_vf(eqn)%sf(x, y, z) = flux_src_vf(eqn)%sf(x, y, z) - Mass_diffu_Flux(eqn - chemxb + 1)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endblock UNDEF_AMD
        end if

    end subroutine s_compute_chemistry_diffusion_flux

end module m_chemistry
