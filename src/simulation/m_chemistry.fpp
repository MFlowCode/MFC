!!>
!! @file   m_chemistry.f90
!! @brief  Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

module m_chemistry

    use ieee_arithmetic

    use m_mpi_proxy
    use m_thermochem
    use m_global_parameters
    use m_finite_differences

    implicit none

    type(int_bounds_info), private :: ix, iy, iz
    type(scalar_field), private :: grads(1:3)

    !$acc declare create(ix, iy, iz)
    !$acc declare create(grads)

contains

    subroutine s_initialize_chemistry_module

        integer :: i

        ix%beg = -buff_size
        if (n > 0) then; iy%beg = -buff_size; else; iy%beg = 0; end if
        if (p > 0) then; iz%beg = -buff_size; else; iz%beg = 0; end if

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

        !$acc update device(ix, iy, iz)

        do i = 1, 3
            @:ALLOCATE(grads(i)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        !$acc kernels
        do i = 1, num_dims
            grads(i)%sf(:, :, :) = 0.0d0
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
                do x = 0, m
                    do y = 0, n
                        do z = 0, p
                            do eqn = chemxb, chemxe
                                ! \nabla \cdot (F)
                                flux_x = (flux_n(1)%vf(eqn)%sf(x - 1, y, z) - &
                                          flux_n(1)%vf(eqn)%sf(x, y, z))/dx(x)

                                #:if num_dims >= 2
                                    flux_x = (flux_n(2)%vf(eqn)%sf(x, y - 1, z) - &
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

        real(kind(0d0)) :: T
        integer :: o
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: E

        real(kind(0d0)) :: rho, omega_m, omega_T
        real(kind(0d0)), dimension(num_species) :: Ys
        real(kind(0d0)), dimension(num_species) :: omega
        real(kind(0d0)), dimension(num_species) :: enthalpies
        real(kind(0d0)) :: cp_mix

        #:if chemistry

            !$acc parallel loop collapse(3) private(rho)
            do x = 0, m
                do y = 0, n
                    do z = 0, p

                        rho = 0d0
                        do eqn = chemxb, chemxe
                            rho = rho + q_cons_qp(eqn)%sf(x, y, z)
                        end do

                        do eqn = chemxb, chemxe
                            Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                        end do

                        dyn_pres = 0d0
                        do i = momxb, momxe
                            dyn_pres = dyn_pres + q_cons_qp(i)%sf(x, y, z)* &
                                       q_prim_qp(i)%sf(x, y, z)/2d0
                        end do

                        call get_temperature(q_cons_qp(E_idx)%sf(x, y, z)/rho - dyn_pres/rho, &
                            & 1200d0, Ys, .true., T)

                        call get_net_production_rates(rho, T, Ys, omega)

                        call get_species_enthalpies_rt(T, enthalpies)

                        do eqn = chemxb, chemxe

                            omega_m = mol_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)
                            omega_T = omega_m*enthalpies(eqn - chemxb + 1)*gas_constant*T

                            rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                        end do

                    end do
                end do
            end do

        #:else

            @:ASSERT(.false., "Chemistry is not enabled")

        #:endif

    end subroutine s_compute_chemistry_reaction_flux

end module m_chemistry
