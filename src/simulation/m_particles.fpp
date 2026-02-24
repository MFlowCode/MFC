!>
!! @file m_particles.f90
!! @brief Contains module m_particles

#:include 'macros.fpp'

!> @brief This module currently contains nothing. Placeholder for any solid particle routines.
module m_particles

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_particles_EL_kernels

    implicit none

contains

end module m_particles
