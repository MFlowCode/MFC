! MFC-owned thermochemistry; derived from Pyrometheus 1.1.1 (MIT).


${gpu_routine}

module ${module_name}

    implicit none

    integer, parameter :: sp = selected_real_kind(6,37)   ! Single precision
    integer, parameter :: dp = selected_real_kind(15,307) ! Double precision

    integer, parameter :: num_elements = ${sol.n_elements}
    integer, parameter :: num_species = ${sol.n_species}
    integer, parameter :: num_reactions = ${sol.n_reactions}
    integer, parameter :: num_falloff = ${len(falloff_reactions)}
    ${real_type}, parameter :: one_atm = ${float_to_fortran(ct.one_atm)}
    ${real_type}, parameter :: gas_constant = ${float_to_fortran(ct.gas_constant)}
    ${real_type}, parameter :: molecular_weights(${sol.n_species}) = &
        (/ ${str_np(sol.molecular_weights)} /)
    ${real_type}, parameter :: inv_molecular_weights(${sol.n_species}) = &
        (/ ${str_np(1/sol.molecular_weights)} /)

    character(len=12), parameter :: species_names(${sol.n_species}) = &
        (/ ${", ".join('"'+'{0: <12}'.format(s)+'"' for s in sol.species_names)} /)

    character(len=4), parameter :: element_names(${sol.n_elements}) = &
        (/ ${", ".join('"'+'{0: <4}'.format(e)+'"' for e in sol.element_names)} /)

contains

    subroutine get_species_name(sp_index, sp_name)

        integer, intent(in) :: sp_index
        character(len=*), intent(out) :: sp_name

        sp_name = species_names(sp_index)

    end subroutine get_species_name

    subroutine get_species_index(sp_name, sp_index)

        character(len=*), intent(in) :: sp_name
        integer, intent(out) :: sp_index

        integer :: idx

        sp_index = 0
        loop:do idx = 1, num_species
            if(trim(adjustl(sp_name)) .eq. trim(species_names(idx))) then
                sp_index = idx
                exit loop
            end if
        end do loop

    end subroutine get_species_index

    subroutine get_element_index(el_name, el_index)

        character(len=*), intent(in) :: el_name
        integer, intent(out) :: el_index

        integer :: idx

        el_index = 0
        loop:do idx = 1, num_elements
            if(trim(adjustl(el_name)) .eq. trim(element_names(idx))) then
                el_index = idx
                exit loop
            end if
        end do loop

    end subroutine get_element_index

    subroutine get_specific_gas_constant(mass_fractions, specific_gas_constant)

        GPU_ROUTINE(get_specific_gas_constant)

        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: specific_gas_constant

        specific_gas_constant = gas_constant * ( &
                %for i in range(sol.n_species):
                    + inv_molecular_weights(${i+1})*mass_fractions(${i+1}) &
                %endfor
                )

    end subroutine get_specific_gas_constant

    subroutine get_density(pressure, temperature, mass_fractions, density)

        GPU_ROUTINE(get_density)

        ${real_type}, intent(in) :: pressure
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: density

        ${real_type} :: mix_mol_weight

        call get_mixture_molecular_weight(mass_fractions, mix_mol_weight)
        density = pressure * mix_mol_weight / (gas_constant * temperature)

    end subroutine get_density

    subroutine get_pressure(density, temperature, mass_fractions, pressure)

        GPU_ROUTINE(get_pressure)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: pressure

        ${real_type} :: mix_mol_weight

        call get_mixture_molecular_weight(mass_fractions, mix_mol_weight)
        pressure = density * gas_constant * temperature / mix_mol_weight

    end subroutine get_pressure

    subroutine get_mixture_molecular_weight(mass_fractions, mix_mol_weight)

        GPU_ROUTINE(get_mixture_molecular_weight)

        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: mix_mol_weight

        mix_mol_weight = 1.0d0 / ( &
                %for i in range(sol.n_species):
                    + inv_molecular_weights(${i+1})*mass_fractions(${i+1}) &
                %endfor
                )

    end subroutine get_mixture_molecular_weight

    subroutine get_concentrations(density, mass_fractions, concentrations)

        GPU_ROUTINE(get_concentrations)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: concentrations

        %for i in range(sol.n_species):
            concentrations(${i+1}) = density * &
                inv_molecular_weights(${i+1}) * mass_fractions(${i+1})
        %endfor

    end subroutine get_concentrations

    subroutine get_mole_fractions(mix_mol_weight, mass_fractions, mole_fractions)

        GPU_ROUTINE(get_mole_fractions)

        ${real_type}, intent(in) :: mix_mol_weight
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: mole_fractions

        %for i in range(sol.n_species):
            mole_fractions(${i+1}) = inv_molecular_weights(${i+1}) * &
                mass_fractions(${i+1}) * mix_mol_weight
        %endfor

    end subroutine get_mole_fractions

    subroutine get_mass_averaged_property(&
        & mass_fractions, spec_property, mix_property)

        GPU_ROUTINE(get_mass_averaged_property)

        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(in), dimension(${sol.n_species}) :: spec_property
        ${real_type}, intent(out) :: mix_property

        mix_property =  ( &
            %for i in range(sol.n_species):
                + inv_molecular_weights(${i+1})*mass_fractions(${i+1}) &
                *spec_property(${i+1}) &
            %endfor
        )

    end subroutine get_mass_averaged_property

    subroutine get_mixture_specific_heat_cp_mass(temperature, mass_fractions, cp_mix)

        GPU_ROUTINE(get_mixture_specific_heat_cp_mass)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: cp_mix

        ${real_type}, dimension(${sol.n_species}) :: cp0_r

        call get_species_specific_heats_r(temperature, cp0_r)
        call get_mass_averaged_property(mass_fractions, cp0_r, cp_mix)
        cp_mix = cp_mix * gas_constant

    end subroutine get_mixture_specific_heat_cp_mass

    subroutine get_mixture_specific_heat_cv_mass(temperature, mass_fractions, cv_mix)

        GPU_ROUTINE(get_mixture_specific_heat_cv_mass)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: cv_mix

        ${real_type}, dimension(${sol.n_species}) :: cp0_r

        call get_species_specific_heats_r(temperature, cp0_r)

        %for i in range(sol.n_species):
            cp0_r(${i+1}) = cp0_r(${i+1}) - 1.d0
        %endfor

        call get_mass_averaged_property(mass_fractions, cp0_r, cv_mix)
        cv_mix = cv_mix * gas_constant

    end subroutine get_mixture_specific_heat_cv_mass

    subroutine get_mixture_enthalpy_mass(temperature, mass_fractions, h_mix)

        GPU_ROUTINE(get_mixture_enthalpy_mass)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: h_mix

        ${real_type}, dimension(${sol.n_species}) :: h0_rt

        call get_species_enthalpies_rt(temperature, h0_rt)
        call get_mass_averaged_property(mass_fractions, h0_rt, h_mix)
        h_mix = h_mix * gas_constant * temperature

    end subroutine get_mixture_enthalpy_mass

    subroutine get_mixture_energy_mass(temperature, mass_fractions, e_mix)

        GPU_ROUTINE(get_mixture_energy_mass)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: e_mix

        ${real_type}, dimension(${sol.n_species}) :: h0_rt

        call get_species_enthalpies_rt(temperature, h0_rt)

        %for i in range(sol.n_species):
            h0_rt(${i+1}) = h0_rt(${i+1}) - 1.d0
        %endfor

        call get_mass_averaged_property(mass_fractions, h0_rt, e_mix)
        e_mix = e_mix * gas_constant * temperature

    end subroutine get_mixture_energy_mass

    subroutine get_species_specific_heats_r(temperature, cp0_r)

        GPU_ROUTINE(get_species_specific_heats_r)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: cp0_r

        %for i, sp in enumerate(sol.species()):
        cp0_r(${i+1}) = ${cgm(ce.poly_to_expr(sp.thermo, "temperature"))}
        %endfor

    end subroutine get_species_specific_heats_r

    subroutine get_species_enthalpies_rt(temperature, h0_rt)

        GPU_ROUTINE(get_species_enthalpies_rt)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: h0_rt

        %for i, sp in enumerate(sol.species()):
        h0_rt(${i+1}) = ${cgm(ce.poly_to_enthalpy_expr(sp.thermo, "temperature"))}
        %endfor

    end subroutine get_species_enthalpies_rt

    subroutine get_species_entropies_r(temperature, s0_r)

        GPU_ROUTINE(get_species_entropies_r)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: s0_r

        %for i, sp in enumerate(sol.species()):
        s0_r(${i+1}) = ${cgm(ce.poly_to_entropy_expr(sp.thermo, "temperature"))}
        %endfor

    end subroutine get_species_entropies_r

    subroutine get_species_gibbs_rt(temperature, g0_rt)

        GPU_ROUTINE(get_species_gibbs_rt)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: g0_rt

        ${real_type}, dimension(${sol.n_species}) :: h0_rt
        ${real_type}, dimension(${sol.n_species}) :: s0_r

        call get_species_enthalpies_rt(temperature, h0_rt)
        call get_species_entropies_r(temperature, s0_r)

        %for i in range(sol.n_species):
            g0_rt(${i+1}) = h0_rt(${i+1}) - s0_r(${i+1})
        %endfor

    end subroutine get_species_gibbs_rt

    subroutine get_equilibrium_constants(temperature, k_eq)

        GPU_ROUTINE(get_equilibrium_constants)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: k_eq

        ${real_type} :: rt
        ${real_type} :: c0

        ${real_type}, dimension(${sol.n_species}) :: g0_rt

        rt = gas_constant * temperature
        c0 = log(one_atm/rt)

        call get_species_gibbs_rt(temperature, g0_rt)

        %for i, react in enumerate(sol.reactions()):
        %if react.reversible:
        k_eq(${i+1}) = ${cgm(
            ce.equilibrium_constants_expr(sol, i, Variable("g0_rt")))}
        %else:
        k_eq(${i+1}) = -0.1d0*temperature
        %endif
        %endfor

    end subroutine get_equilibrium_constants

    subroutine get_temperature( &
        & enthalpy_or_energy, t_guess, mass_fractions, do_energy, temperature)

        GPU_ROUTINE(get_temperature)

        logical, intent(in) :: do_energy
        ${real_type}, intent(in)  :: enthalpy_or_energy
        ${real_type}, intent(in)  :: t_guess
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: temperature

        integer :: iter
        integer,      parameter :: num_iter = 500
        ${real_type}, parameter :: tol = 1.0d-06

        ${real_type} :: iter_temp
        ${real_type} :: iter_energy
        ${real_type} :: iter_energy_deriv
        ${real_type} :: iter_rhs
        ${real_type} :: iter_deriv

        iter_rhs = 0.d0
        iter_deriv = 1.d0
        iter_temp = t_guess

        do iter = 1, num_iter
            if(do_energy) then
                call get_mixture_specific_heat_cv_mass(&
                    & iter_temp, mass_fractions, iter_energy_deriv)
                call get_mixture_energy_mass(iter_temp, mass_fractions, iter_energy)
            else
                call get_mixture_specific_heat_cp_mass(&
                    & iter_temp, mass_fractions, iter_energy_deriv)
                call get_mixture_enthalpy_mass(&
                    & iter_temp, mass_fractions, iter_energy)
            endif
            iter_rhs = enthalpy_or_energy - iter_energy
            iter_deriv = (-1.d0)*iter_energy_deriv
            iter_temp = iter_temp - iter_rhs / iter_deriv
            if(abs(iter_rhs/iter_deriv) .lt. tol) exit
        end do

        temperature = iter_temp

    end subroutine get_temperature

    %if falloff_reactions:
    subroutine get_falloff_rates(temperature, concentrations, k_fwd)

        GPU_ROUTINE(get_falloff_rates)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: concentrations
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: k_fwd

        ${real_type}, dimension(${len(falloff_reactions)}) :: k_high
        ${real_type}, dimension(${len(falloff_reactions)}) :: k_low
        ${real_type}, dimension(${len(falloff_reactions)}) :: reduced_pressure
        ${real_type}, dimension(${len(falloff_reactions)}) :: falloff_center
        ${real_type}, dimension(${len(falloff_reactions)}) :: falloff_factor
        ${real_type}, dimension(${len(falloff_reactions)}) :: falloff_function

        %for i, (_, react) in enumerate(falloff_reactions):
        k_high(${i+1}) = ${cgm(ce.rate_coefficient_expr(
                                react.rate.high_rate,
                                Variable("temperature")))}
        %endfor

        %for i, (_, react) in enumerate(falloff_reactions):
        k_low(${i+1}) = ${cgm(ce.rate_coefficient_expr(
                                react.rate.low_rate,
                                Variable("temperature")))}
        %endfor

        %for i, (_, react) in enumerate(falloff_reactions):
        reduced_pressure(${i+1}) = (${cgm(
            ce.third_body_efficiencies_expr(sol,
                                            react,
                                            Variable("concentrations")))})*k_low(${i+1})/k_high(${i+1})
        %endfor

        %for i, (_, react) in enumerate(falloff_reactions):
        falloff_center(${i+1}) = ${cgm(ce.troe_falloff_center_expr(
            react, Variable("temperature")))}
        %endfor

        %for i, (_, react) in enumerate(falloff_reactions):
        falloff_factor(${i+1}) = ${cgm(ce.troe_falloff_factor_expr(react, i,
            Variable("reduced_pressure"), Variable("falloff_center")))}
        %endfor

        %for i, (_, react) in enumerate(falloff_reactions):
        falloff_function(${i+1}) = ${cgm(ce.falloff_function_expr(
            react, i,
            Variable("falloff_factor"),
            Variable("falloff_center")))}
        %endfor

        %for i, (j, react) in enumerate(falloff_reactions):
        k_fwd(${j+1}) = k_high(${i+1})*falloff_function(${i+1}) * &
            reduced_pressure(${i+1})/(1.d0 + reduced_pressure(${i+1}))
        %endfor

    end subroutine get_falloff_rates

    %endif
    subroutine get_fwd_rate_coefficients(temperature, concentrations, k_fwd)

        GPU_ROUTINE(get_fwd_rate_coefficients)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: concentrations
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: k_fwd

        %if falloff_reactions:
        ${real_type}, dimension(${len(falloff_reactions)}) :: k_falloff
        %endif

        %for i, react in enumerate(sol.reactions()):
        %if i in falloff_indices:
        k_fwd(${i+1}) = 0.d0
        %else:
        k_fwd(${i+1}) = ${cgm(ce.rate_coefficient_expr(react.rate,
                            Variable("temperature")))}
        %endif
        %endfor

        %for j, react in three_body_reactions:
        k_fwd(${j+1}) = k_fwd(${j+1}) * ( &
            ${cgm(ce.third_body_efficiencies_expr(
            sol, react, Variable("concentrations")))})
        %endfor

        %if falloff_reactions:
        call get_falloff_rates(temperature, concentrations, k_fwd)
        %endif

    end subroutine get_fwd_rate_coefficients

    subroutine get_net_rates_of_progress(temperature, concentrations, r_net)

        GPU_ROUTINE(get_net_rates_of_progress)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: concentrations
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: r_net

        ${real_type}, dimension(${sol.n_reactions}) :: k_fwd
        ${real_type}, dimension(${sol.n_reactions}) :: log_k_eq

        call get_fwd_rate_coefficients(temperature, concentrations, k_fwd)
        call get_equilibrium_constants(temperature, log_k_eq)
        %for i in range(sol.n_reactions):
        r_net(${i+1}) = ${cgm(ce.rate_of_progress_expr(sol, i,
                        Variable("concentrations"),
                        Variable("k_fwd"), Variable("log_k_eq")))}
        %endfor

    end subroutine get_net_rates_of_progress

    subroutine get_net_production_rates(density, temperature, mass_fractions, omega)

        GPU_ROUTINE(get_net_production_rates)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: omega

        ${real_type}, dimension(${sol.n_species})   :: concentrations
        ${real_type}, dimension(${sol.n_reactions}) :: r_net

        call get_concentrations(density, mass_fractions, concentrations)
        call get_net_rates_of_progress(temperature, concentrations, r_net)

        %for i, sp in enumerate(sol.species()):
        omega(${i+1}) = ${cgm(ce.production_rate_expr(sol,
            sp.name, Variable("r_net")))}
        %endfor

    end subroutine get_net_production_rates

    subroutine get_fwd_rates_of_progress(temperature, concentrations, r_fwd)

        GPU_ROUTINE(get_fwd_rates_of_progress)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: concentrations
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: r_fwd

        ${real_type}, dimension(${sol.n_reactions}) :: k_fwd

        call get_fwd_rate_coefficients(temperature, concentrations, k_fwd)
        %for i in range(sol.n_reactions):
        r_fwd(${i+1}) = ${cgm(ce.fwd_rate_of_progress_expr(sol, i,
                        Variable("concentrations"), Variable("k_fwd")))}
        %endfor

    end subroutine get_fwd_rates_of_progress

    subroutine get_rev_rates_of_progress(temperature, concentrations, r_rev)

        GPU_ROUTINE(get_rev_rates_of_progress)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: concentrations
        ${real_type}, intent(out), dimension(${sol.n_reactions}) :: r_rev

        ${real_type}, dimension(${sol.n_reactions}) :: k_fwd
        ${real_type}, dimension(${sol.n_reactions}) :: log_k_eq

        call get_fwd_rate_coefficients(temperature, concentrations, k_fwd)
        call get_equilibrium_constants(temperature, log_k_eq)
        %for i in range(sol.n_reactions):
        r_rev(${i+1}) = ${cgm(ce.rev_rate_of_progress_expr(sol, i,
                        Variable("concentrations"),
                        Variable("k_fwd"), Variable("log_k_eq")))}
        %endfor

    end subroutine get_rev_rates_of_progress

    subroutine get_creation_rates(density, temperature, mass_fractions, cdot)

        GPU_ROUTINE(get_creation_rates)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: cdot

        ${real_type}, dimension(${sol.n_species})   :: concentrations
        ${real_type}, dimension(${sol.n_reactions}) :: r_fwd, r_rev

        call get_concentrations(density, mass_fractions, concentrations)
        call get_fwd_rates_of_progress(temperature, concentrations, r_fwd)
        call get_rev_rates_of_progress(temperature, concentrations, r_rev)

        %for i, sp in enumerate(sol.species()):
        cdot(${i+1}) = ${cgm(ce.creation_rate_expr(sol, sp.name,
            Variable("r_fwd"), Variable("r_rev")))}
        %endfor

    end subroutine get_creation_rates

    subroutine get_destruction_rates(density, temperature, mass_fractions, ddot)

        GPU_ROUTINE(get_destruction_rates)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: ddot

        ${real_type}, dimension(${sol.n_species})   :: concentrations
        ${real_type}, dimension(${sol.n_reactions}) :: r_fwd, r_rev

        call get_concentrations(density, mass_fractions, concentrations)
        call get_fwd_rates_of_progress(temperature, concentrations, r_fwd)
        call get_rev_rates_of_progress(temperature, concentrations, r_rev)

        %for i, sp in enumerate(sol.species()):
        ddot(${i+1}) = ${cgm(ce.destruction_rate_expr(sol, sp.name,
            Variable("r_fwd"), Variable("r_rev")))}
        %endfor

    end subroutine get_destruction_rates

    subroutine get_creation_destruction_rates(density, temperature, &
        mass_fractions, cdot, ddot)

        GPU_ROUTINE(get_creation_destruction_rates)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in),  dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: cdot, ddot

        ${real_type}, dimension(${sol.n_species})   :: concentrations
        ${real_type}, dimension(${sol.n_reactions}) :: r_fwd, r_rev

        call get_concentrations(density, mass_fractions, concentrations)
        call get_fwd_rates_of_progress(temperature, concentrations, r_fwd)
        call get_rev_rates_of_progress(temperature, concentrations, r_rev)

        %for i, sp in enumerate(sol.species()):
        cdot(${i+1}) = ${cgm(ce.creation_rate_expr(sol, sp.name,
            Variable("r_fwd"), Variable("r_rev")))}
        ddot(${i+1}) = ${cgm(ce.destruction_rate_expr(sol, sp.name,
            Variable("r_fwd"), Variable("r_rev")))}
        %endfor

    end subroutine get_creation_destruction_rates

    subroutine get_species_viscosities(temperature, viscosities)

        GPU_ROUTINE(get_species_viscosities)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: viscosities

        %for sp in range(sol.n_species):
        viscosities(${sp+1}) = ${cgm(ce.viscosity_polynomial_expr(
            sol.get_viscosity_polynomial(sp),
            Variable("temperature")))}
        %endfor

    end subroutine get_species_viscosities

    subroutine get_species_thermal_conductivities(temperature, conductivities)

        GPU_ROUTINE(get_species_thermal_conductivities)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}) :: conductivities

        %for sp in range(sol.n_species):
        conductivities(${sp+1}) = ${cgm(ce.conductivity_polynomial_expr(
            sol.get_thermal_conductivity_polynomial(sp),
            Variable("temperature")))}
        %endfor

    end subroutine get_species_thermal_conductivities

    subroutine get_species_binary_mass_diffusivities(temperature, diffusivities)

        GPU_ROUTINE(get_species_binary_mass_diffusivities)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(out), dimension(${sol.n_species}, ${sol.n_species})&
            :: diffusivities

        %for i in range(sol.n_species):
        %for j in range(sol.n_species):
        diffusivities(${i + 1}, ${j + 1}) = ${cgm(ce.diffusivity_polynomial_expr(
            sol.get_binary_diff_coeffs_polynomial(i, j),
            Variable("temperature")))}
        %endfor
        %endfor

    end subroutine get_species_binary_mass_diffusivities

    subroutine get_mixture_viscosity_mixavg(&
        temperature, mass_fractions, mixture_viscosity_mixavg)

        GPU_ROUTINE(get_mixture_viscosity_mixavg)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: mixture_viscosity_mixavg

        ${real_type} :: mix_mol_weight
        ${real_type}, dimension(${sol.n_species}) :: &
            mole_fractions, viscosities, mix_rule_f

        call get_mixture_molecular_weight(mass_fractions, mix_mol_weight)
        call get_mole_fractions(mix_mol_weight, mass_fractions, mole_fractions)
        call get_species_viscosities(temperature, viscosities)

        %for sp in range(sol.n_species):
        mix_rule_f(${sp + 1}) = ${cgm(ce.viscosity_mixture_rule_wilke_expr(sol, sp,
            Variable("mole_fractions"), Variable("viscosities")))}
        %endfor

        mixture_viscosity_mixavg = sum(mole_fractions*viscosities/mix_rule_f)

    end subroutine get_mixture_viscosity_mixavg

    subroutine get_mixture_thermal_conductivity_mixavg(temperature, &
        mass_fractions, mixture_thermal_conductivity_mixavg)

        GPU_ROUTINE(get_mixture_thermal_conductivity_mixavg)

        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out) :: mixture_thermal_conductivity_mixavg

        ${real_type} :: mix_mol_weight
        ${real_type}, dimension(${sol.n_species}) :: mole_fractions, conductivities

        call get_mixture_molecular_weight(mass_fractions, mix_mol_weight)
        call get_mole_fractions(mix_mol_weight, mass_fractions, mole_fractions)
        call get_species_thermal_conductivities(temperature, conductivities)

        mixture_thermal_conductivity_mixavg = 0.5*(&
            sum(mole_fractions*conductivities) + &
            1/sum(mole_fractions/conductivities))

    end subroutine get_mixture_thermal_conductivity_mixavg

    subroutine get_species_mass_diffusivities_mixavg(&
        pressure, temperature, mass_fractions, mass_diffusivities_mixavg)

        GPU_ROUTINE(get_species_mass_diffusivities_mixavg)

        ${real_type}, intent(in) :: pressure, temperature
        ${real_type}, intent(in), dimension(${sol.n_species}) :: mass_fractions
        ${real_type}, intent(out), dimension(${sol.n_species}) :: &
            mass_diffusivities_mixavg

        ${real_type} :: mix_mol_weight
        ${real_type}, dimension(${sol.n_species}) :: mole_fractions, x_sum, denom
        ${real_type}, dimension(${sol.n_species}, ${sol.n_species}) :: bdiff_ij

        call get_mixture_molecular_weight(mass_fractions, mix_mol_weight)
        call get_mole_fractions(mix_mol_weight, mass_fractions, mole_fractions)
        call get_species_binary_mass_diffusivities(temperature, bdiff_ij)

        %for sp in range(sol.n_species):
        x_sum(${sp + 1}) = ${cgm(ce.diffusivity_mixture_rule_denom_expr(
                sol, sp, Variable("mole_fractions"), Variable("bdiff_ij")))}
        %endfor

        %for sp in range(sol.n_species):
        denom(${sp + 1}) = x_sum(${sp + 1}) - &
            mole_fractions(${sp + 1})/bdiff_ij(${sp + 1}, ${sp + 1})
        %endfor

        %for sp in range(sol.n_species):
        if (denom(${sp + 1}) .gt. 0d0) then
        mass_diffusivities_mixavg(${sp + 1}) = &
            (mix_mol_weight - &
                mole_fractions(${sp + 1})*molecular_weights(${sp + 1}))&
            /(pressure * mix_mol_weight * denom(${sp + 1}))
        else
        mass_diffusivities_mixavg(${sp + 1}) = &
            bdiff_ij(${sp + 1}, ${sp + 1}) / pressure
        end if
        %endfor

    end subroutine get_species_mass_diffusivities_mixavg

end module ${module_name}