!>
!! @file m_hyperelastic.f90
!! @brief Contains module m_hyperelastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hyperelastic model
module m_hyperelastic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_hyperelastic_module, &
 s_calculate_cauchy_stress

contains

    subroutine s_initialize_hyperelastic_module()

    end subroutine s_initialize_hyperelastic_module

        !type(int_bounds_info), optional, intent(IN) :: ix, iy, iz

        !do l = izb, ize
        !    do k = iyb, iye
        !        do j = ixb, ixe

    subroutine s_calculate_cauchy_from_btensor(q_btensor_vf, ix, iy, iz, q_prim_vf)
        type(scalar_field), dimension(num_dims**2 + 1), intent(IN) :: btensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: sigma
        integer, intent(IN) :: j, k, l

        real(kind(0d0)), dimension(num_dims**2) :: tensorb, devbtensor
        real(kind(0d0)) :: jacobian
        integer :: i !< Generic loop iterators

        ! extracting the nxn tensor for the calculation
        do i = 1, num_dims**2
            tensorb(i) = btensor(i)%sf(j, k, l)
        end do
        jacobian = btensor(num_dims**2 + 1)%sf(j, k, l)
        call s_calculate_deviatoric(tensorb, devbtensor)
        sigma(:) = devbtensor(:)/jacobian

    end subroutine s_calculate_cauchy_stress

    function f_trace(symtensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: symtensor
        real(kind(0d0)) ::  f_trace

        f_trace  = symtensor(1)
        if (num_dims == 2) then
            f_trace = symtensor(1)+symtensor(3)
        else
            f_trace = symtensor(1)+symtensor(4)+symtensor(6)
        endif
    end function f_trace

     subroutine s_calculate_deviatoric(symtensor, devtensor)
        real(kind(0d0)), dimension(num_dims*2 + 1), intent(IN) :: symtensor
        real(kind(0d0)), dimension(num_dims*2), intent(OUT) :: devtensor
        real(kind(0d0)) :: trace
        devtensor = symtensor
        trace = f_trace(symtensor)
        devtensor(1) = symtensor(1) - (1d0/3d0)*trace      
        if (num_dims == 2) then
             devtensor(3) = symtensor(3) - (1d0/3d0)*trace
        else
             devtensor(4) = symtensor(4) - (1d0/3d0)*trace
             devtensor(6) = symtensor(6) - (1d0/3d0)*trace
        end if     
     end subroutine s_calculate_deviatoric

end module m_hyperelastic
