!!>
!! @file   m_chemistry.f90
!! @brief  Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

module m_chemistry

    use ieee_arithmetic

    use m_mpi_proxy
    use m_thermochem, only: &
        num_species, mol_weights, get_net_production_rates

    use m_global_parameters

    use m_finite_differences

    implicit none

    type(scalar_field), private :: grads(1:3)

    !$acc declare create(grads)

contains

    subroutine s_initialize_chemistry_module

        integer :: i

        do i = 1, 3
            @:ALLOCATE(grads(i)%sf(&
                & idwbuff(1)%beg:idwbuff(1)%end, &
                & idwbuff(2)%beg:idwbuff(2)%end, &
                & idwbuff(3)%beg:idwbuff(3)%end))
        end do

        !$acc kernels
        do i = 1, num_dims
            grads(i)%sf(:, :, :) = 0.0_wp
        end do
        !$acc end kernels

    end subroutine s_initialize_chemistry_module

    subroutine s_finalize_chemistry_module

        deallocate (grads(1)%sf, grads(2)%sf, grads(3)%sf)

    end subroutine s_finalize_chemistry_module

    subroutine s_compute_chemistry_advection_flux(flux_n, rhs_vf)

        type(vector_field), dimension(:), intent(IN) :: flux_n
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        type(int_bounds_info) :: ix, iy, iz

        integer :: x, y, z
        integer :: eqn

        real(kind(0d0)) :: flux_x, flux_y, flux_z

        #:for num_dims in range(1, 4)
            if (num_dims == ${num_dims}$) then
                !$acc parallel loop collapse(4) gang vector default(present) &
                !$acc private(flux_x, flux_y, flux_z)
                do z = idwint(3)%beg, idwint(3)%end
                    do y = idwint(2)%beg, idwint(2)%end
                        do x = idwint(1)%beg, idwint(1)%end
                            do eqn = chemxb, chemxe
                                ! \nabla \cdot (F)
                                flux_x = (flux_n(1)%vf(eqn)%sf(x - 1, y, z) - &
                                          flux_n(1)%vf(eqn)%sf(x, y, z))/dx(x)

                                #:if num_dims >= 2
                                    flux_y = (flux_n(2)%vf(eqn)%sf(x, y - 1, z) - &
                                              flux_n(2)%vf(eqn)%sf(x, y, z))/dy(y)
                                #:else
                                    flux_y = 0d0
                                #:endif

                                #:if num_dims == 3
                                    flux_z = (flux_n(3)%vf(eqn)%sf(x, y, z - 1) - &
                                              flux_n(3)%vf(eqn)%sf(x, y, z))/dz(z)
                                #:else
                                    flux_z = 0d0
                                #:endif

                                rhs_vf(eqn)%sf(x, y, z) = flux_x + flux_y + flux_z
                            end do

                            rhs_vf(T_idx)%sf(x, y, z) = 0d0
                        end do
                    end do
                end do
            end if
        #:endfor

    end subroutine s_compute_chemistry_advection_flux

    subroutine s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp, q_prim_qp)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf, q_cons_qp, q_prim_qp

        integer :: i

        integer :: x, y, z
        integer :: eqn

        real(wp) :: T
        integer :: o
        real(wp) :: dyn_pres
        real(wp) :: E


        real(wp) :: rho, omega_m
        real(wp), dimension(num_species) :: Ys
        real(wp), dimension(num_species) :: omega
        real(wp) :: cp_mix

        if (chemistry) then
            !$acc parallel loop collapse(3) gang vector default(present) &
            !$acc private(Ys, omega)
            do z = idwint(3)%beg, idwint(3)%end
                do y = idwint(2)%beg, idwint(2)%end
                    do x = idwint(1)%beg, idwint(1)%end

                        !$acc loop seq
                        do eqn = chemxb, chemxe
                            Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                        end do

                        rho = q_cons_qp(contxe)%sf(x, y, z)
                        T = q_prim_qp(T_idx)%sf(x, y, z)

                        call get_net_production_rates(rho, T, Ys, omega)

                        !$acc loop seq
                        do eqn = chemxb, chemxe

                            omega_m = mol_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)

                            rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                        end do

                    end do
                end do
            end do

        else

            @:ASSERT(.false., "Chemistry is not enabled")

        end if

    end subroutine s_compute_chemistry_reaction_flux

end module m_chemistry
