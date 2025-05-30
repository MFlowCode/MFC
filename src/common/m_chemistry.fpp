!!>
!! @file   m_chemistry.f90
!! @brief  Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

module m_chemistry

    use m_thermochem, only: &
        num_species, molecular_weights, get_temperature, get_net_production_rates, &
        gas_constant, get_mixture_molecular_weight

    use m_global_parameters

    implicit none

contains

    subroutine s_compute_q_T_sf(q_T_sf, q_cons_vf, bounds)

        ! Initialize the temperature field at the start of the simulation to
        ! reasonable values. Temperature is computed the regular way using the
        ! conservative variables.

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, eqn
        real(wp) :: energy
        real(wp), dimension(num_species) :: Ys

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    !$acc loop seq
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = &
                            q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    end do

                    ! e = E - 1/2*|u|^2
                    ! cons. E_idx     = \rho E
                    ! cons. contxb    = \rho         (1-fluid model)
                    ! cons. momxb + i = \rho u_i
                    energy = q_cons_vf(E_idx)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    !$acc loop seq
                    do eqn = momxb, momxe
                        energy = energy - &
                                 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z))**2._wp
                    end do

                    call get_temperature(energy, dflt_T_guess, Ys, .true., q_T_sf%sf(x, y, z))
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
                    !$acc loop seq
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

        !$acc parallel loop collapse(3) gang vector default(present) &
        !$acc private(Ys, omega)
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end

                    !$acc loop seq
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                    end do

                    rho = q_cons_qp(contxe)%sf(x, y, z)
                    T = q_T_sf%sf(x, y, z)

                    call get_net_production_rates(rho, T, Ys, omega)

                    !$acc loop seq
                    do eqn = chemxb, chemxe

                        omega_m = molecular_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)

                        rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                    end do

                end do
            end do
        end do

    end subroutine s_compute_chemistry_reaction_flux

end module m_chemistry
