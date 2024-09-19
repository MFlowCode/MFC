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

    #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
        subroutine s_compute_chemistry_rhs_${XYZ}$ (flux_n, rhs_vf, flux_src_vf, q_prim_vf)

            type(vector_field), dimension(:), intent(IN) :: flux_n
            type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf, flux_src_vf, q_prim_vf
            type(int_bounds_info) :: ix, iy, iz

            integer :: x, y, z
            integer :: eqn

            integer, parameter :: mx = ${1 if NORM_DIR == 1 else 0}$
            integer, parameter :: my = ${1 if NORM_DIR == 2 else 0}$
            integer, parameter :: mz = ${1 if NORM_DIR == 3 else 0}$

            !$acc parallel loop collapse(4) present(rhs_vf, flux_n)
            do x = 0, m
                do y = 0, n
                    do z = 0, p

                        do eqn = chemxb, chemxe

                            ! \nabla \cdot (F)
                            rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + &
                                                      (flux_n(${NORM_DIR}$)%vf(eqn)%sf(x - mx, y - my, z - mz) - &
                                                       flux_n(${NORM_DIR}$)%vf(eqn)%sf(x, y, z))/d${XYZ}$ (${XYZ}$)

                        end do

                    end do
                end do
            end do

        end subroutine s_compute_chemistry_rhs_${XYZ}$
    #:endfor

    subroutine s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp, q_prim_qp)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf, q_cons_qp, q_prim_qp

        integer :: i

        integer :: x, y, z
        integer :: eqn

        real(kind(0d0)) :: T
        integer :: o
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: E

        real(kind(0d0)) :: rho
        real(kind(1.d0)), dimension(num_species) :: Ys
        real(kind(1.d0)), dimension(num_species) :: omega
        real(kind(0d0)), dimension(num_species) :: enthalpies
        real(kind(0d0)) :: cp_mix

        #:if chemistry

            !$acc parallel loop collapse(4) private(rho)
            do x = 0, m
                do y = 0, n
                    do z = 0, p

                        ! Maybe use q_prim_vf instead?
                        rho = 0d0
                        do eqn = chemxb, chemxe
                            rho = rho + q_cons_qp(eqn)%sf(x, y, z)
                        end do

                        do eqn = chemxb, chemxe
                            Ys(eqn - chemxb + 1) = q_cons_qp(eqn)%sf(x, y, z)/rho
                        end do

                        dyn_pres = 0d0

                        do i = momxb, momxe
                            dyn_pres = dyn_pres + rho*q_cons_qp(i)%sf(x, y, z)* &
                                       q_cons_qp(i)%sf(x, y, z)/2d0
                        end do

                        call get_temperature(.true., q_cons_qp(E_idx)%sf(x, y, z) - dyn_pres, &
                            & q_prim_qp(tempxb)%sf(x, y, z), Ys, T)

                        call get_net_production_rates(rho, T, Ys, omega)

                        q_cons_qp(tempxb)%sf(x, y, z) = T
                        q_prim_qp(tempxb)%sf(x, y, z) = T

                        !print*, x, y, z, T, rho, Ys, omega, q_cons_qp(E_idx)%sf(x, y, z), dyn_pres

                        do eqn = chemxb, chemxe

                            rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + &
                                                      mol_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)

                        end do

                    end do
                end do
            end do

        #:else

            @:ASSERT(.false., "Chemistry is not enabled")

        #:endif

    end subroutine s_compute_chemistry_reaction_flux

    subroutine s_chemistry_normalize_cons(q_cons_qp)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_qp

        integer :: x, y, z
        integer :: eqn

        !$acc parallel loop collapse(4)
        do x = ix%beg, ix%end
            do y = iy%beg, iy%end
                do z = iz%beg, iz%end

                    do eqn = chemxb, chemxe
                        q_cons_qp(eqn)%sf(x, y, z) = max(0d0, q_cons_qp(eqn)%sf(x, y, z))
                    end do

                end do
            end do
        end do

    end subroutine s_chemistry_normalize_cons

end module m_chemistry
