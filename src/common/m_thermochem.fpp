#:include 'case.fpp'

module m_thermochem

    #:if chemistry
        use m_pyrometheus
    #:else
        integer, parameter :: num_species = 0
        character(len=:), allocatable, dimension(:) :: species_names
    #:endif

end module m_thermochem
