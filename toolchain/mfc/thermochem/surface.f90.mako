! Generated from a Cantera interface mechanism by MFC. Do not edit.

${gpu_routine}

module ${module_name}

    use m_thermochem, only: ${kind}, num_species, gas_constant, get_concentrations, get_species_enthalpies_rt

    implicit none

    private
    public :: num_surface_reactions
    public :: get_surface_rates_of_progress
    public :: get_surface_net_production_rates
    public :: get_surface_reaction_heat_flux

    integer, parameter :: num_surface_reactions = ${len(reactions)}

contains

    !> Rates of progress [kmol/m^2/s] of the heterogeneous reactions.
    subroutine get_surface_rates_of_progress(density, temperature, mass_fractions, rates_of_progress)

        GPU_ROUTINE(get_surface_rates_of_progress)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(num_species) :: mass_fractions
        ${real_type}, intent(out), dimension(max(num_surface_reactions, 1)) :: rates_of_progress

        %if reactions:
        ${real_type}, dimension(num_species) :: concentrations

        %endif
        rates_of_progress = ${float_to_fortran(0)}
        %if reactions:
        call get_concentrations(density, mass_fractions, concentrations)
        %endif
        %for i, reaction in enumerate(reactions):
        ! ${reaction.equation}
        rates_of_progress(${i + 1}) = ${cgm(rate_of_progress(reaction))}
        %endfor

    end subroutine get_surface_rates_of_progress

    !> Net molar production rates [kmol/m^2/s] of gas species at the surface.
    subroutine get_surface_net_production_rates(density, temperature, mass_fractions, omega_s)

        GPU_ROUTINE(get_surface_net_production_rates)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(num_species) :: mass_fractions
        ${real_type}, intent(out), dimension(num_species) :: omega_s

        ${real_type}, dimension(max(num_surface_reactions, 1)) :: rates_of_progress

        call get_surface_rates_of_progress(density, temperature, mass_fractions, rates_of_progress)
        omega_s = ${float_to_fortran(0)}
        %for k, expr in production_rates:
        omega_s(${k + 1}) = ${cgm(expr)}
        %endfor

    end subroutine get_surface_net_production_rates

    !> Heat released at the surface [W/m^2]; positive for exothermic reactions.
    subroutine get_surface_reaction_heat_flux(density, temperature, mass_fractions, q_rxn)

        GPU_ROUTINE(get_surface_reaction_heat_flux)

        ${real_type}, intent(in) :: density
        ${real_type}, intent(in) :: temperature
        ${real_type}, intent(in), dimension(num_species) :: mass_fractions
        ${real_type}, intent(out) :: q_rxn

        ${real_type}, dimension(max(num_surface_reactions, 1)) :: rates_of_progress
        ${real_type}, dimension(num_species) :: h_rt
        ${real_type}, dimension(${max(len(bulk), 1)}) :: bulk_h_rt

        call get_surface_rates_of_progress(density, temperature, mass_fractions, rates_of_progress)
        call get_species_enthalpies_rt(temperature, h_rt)
        bulk_h_rt = ${float_to_fortran(0)}
        %for j, (name, h_rt_expr) in enumerate(bulk):
        ! ${name}
        bulk_h_rt(${j + 1}) = ${cgm(h_rt_expr)}
        %endfor

        q_rxn = ${float_to_fortran(0)}
        %for i, reaction in enumerate(reactions):
        q_rxn = q_rxn - (${cgm(reaction_enthalpy_rt(reaction))}) * rates_of_progress(${i + 1})
        %endfor
        q_rxn = q_rxn*gas_constant*temperature

    end subroutine get_surface_reaction_heat_flux

end module ${module_name}
