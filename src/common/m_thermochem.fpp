#:include 'case.fpp'

module m_thermochem

    #:if chemistry
        use m_pyrometheus, only: &
            num_species, species_names, gas_constant, mol_weights, &
            get_temperature, get_net_production_rates, get_pressure, &
            get_mixture_molecular_weight, get_mixture_energy_mass, &
            get_temperature
    #:endif

    implicit none

    #:if not chemistry
        integer, parameter :: num_species = 0
        character(len=:), allocatable, dimension(:) :: species_names
    #:endif

end module m_thermochem
