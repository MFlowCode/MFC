!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion
!#:include 'macros.fpp'

!> @brief This module consists of subroutines used in the calculation of matrix
!!              operations for the reference map tensor

module m_rmt_tensor_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_calculate_btensor, &
 f_elastic_energy

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), grad_xi)
    !$acc declare link(grad_xi)
#else
    real(kind(0d0)), allocatable, dimension(:) :: grad_xi
    ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
    ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
    ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
    !$acc declare create(grad_xi)
#endif

contains

    subroutine s_calculate_btensor(q_prim_vf, j, k, l, btensor)
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(num_dims*(num_dims+1)/2 + 1), intent(OUT) :: btensor
        integer, intent(IN) :: j, k, l
        real(kind(0d0)), dimension(num_dims**2) :: grad_xi, ftensor, tensorb
        integer :: i

        ! calculate the grad_xi, grad_xi is a nxn tensor
        call s_compute_grad_xi(q_prim_vf, j, k, l, grad_xi)
        ! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        call s_calculate_ainverse(grad_xi,ftensor)
        ! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        call s_calculate_atransposea(ftensor,tensorb)
        ! btensor is symmetric, save the data space
        ! 1: 1D, 3: 2D, 6: 3D
        btensor(1)%sf(j,k,l) = tensorb(1)
        if (num_dims > 1) then ! 2D
           btensor(2)%sf(j,k,l) = tensorb(2)
           btensor(3)%sf(j,k,l) = tensorb(4)
        end if
        if (num_dims > 2) then ! 3D
           btensor(3)%sf(j,k,l) = tensorb(3)
           btensor(4)%sf(j,k,l) = tensorb(5)
           btensor(5)%sf(j,k,l) = tensorb(6)
           btensor(6)%sf(j,k,l) = tensorb(9)
        end if
        ! store the determinant at the last entry of the btensor sf
        btensor(b_size)%sf(j,k,l) = f_determinant(ftensor)

    do i = 1, size(btensor)
      print*, 'btensor(', i, ')%sf(', j, ',', k, ',', l, ') = ', btensor(i)%sf(j,k,l)
    end do   
     
    end subroutine s_calculate_btensor

    function f_determinant(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)) :: f_determinant

        if (num_dims == 1) then
            f_determinant = tensor(1)
        elseif (num_dims == 2) then
            f_determinant = tensor(1)*tensor(4) - tensor(2)*tensor(3)
        else
            f_determinant = tensor(1)*(tensor(5)*tensor(9) - tensor(6)*tensor(8)) &
                            - tensor(2)*(tensor(4)*tensor(9) - tensor(6)*tensor(7)) &
                            + tensor(3)*(tensor(4)*tensor(8) - tensor(5)*tensor(7))
        end if
        ! error checking
        if (f_determinant == 0) then
            print *, 'ERROR: Determinant was zero'
            call s_mpi_abort()
        end if
    end function f_determinant

    subroutine s_calculate_atransposea(tensor, ata)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ata

        ata(1) = tensor(1)**2
        if (num_dims == 2) then
            ata(1) = ata(1) + tensor(3)**2
            ata(2) = tensor(1)*tensor(2) + tensor(3)*tensor(4)
            ata(3) = ata(2)
            ata(4) = tensor(2)**2 + tensor(4)**2
        elseif (num_dims == 3) then
            ata(1) = ata(1) + tensor(4)**2 + tensor(7)**2
            ata(5) = tensor(2) + tensor(5)**2 + tensor(8)**2
            ata(9) = tensor(3) + tensor(6)**2 + tensor(9)**2
            ata(2) = tensor(1)*tensor(2) + tensor(4)*tensor(5) + tensor(7)*tensor(8)
            ata(3) = tensor(1)*tensor(3) + tensor(4)*tensor(6) + tensor(7)*tensor(9)
            ata(6) = tensor(2)*tensor(3) + tensor(5)*tensor(6) + tensor(8)*tensor(9)
            ata(4) = ata(2)
            ata(7) = ata(3)
            ata(8) = ata(4)
        end if
    end subroutine s_calculate_atransposea

    subroutine s_calculate_adjointa(tensor, dja)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: dja

        if (num_dims == 1) then
            dja(1) = 1
        elseif (num_dims == 2) then
            dja(1) = tensor(4)
            dja(2) = -tensor(3)
            dja(3) = -tensor(2)
            dja(4) = tensor(1)
        elseif (num_dims == 3) then
            dja(1) = tensor(5)*tensor(9) - tensor(6)*tensor(8)
            dja(2) = -(tensor(2)*tensor(9) - tensor(3)*tensor(8))
            dja(3) = tensor(2)*tensor(6) - tensor(3)*tensor(5)
            dja(4) = -(tensor(4)*tensor(9) - tensor(6)*tensor(7))
            dja(5) = tensor(1)*tensor(9) - tensor(3)*tensor(7)
            dja(6) = -(tensor(1)*tensor(6) - tensor(4)*tensor(3))
            dja(7) = tensor(4)*tensor(8) - tensor(5)*tensor(7)
            dja(8) = -(tensor(1)*tensor(8) - tensor(2)*tensor(7))
            dja(9) = tensor(1)*tensor(5) - tensor(2)*tensor(4)
        end if
    end subroutine s_calculate_adjointa

    subroutine s_calculate_ainverse(tensor, ainv)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ainv
        real(kind(0d0)), dimension(num_dims**2) :: dja
        real(kind(0d0)) :: det

        call s_calculate_adjointa(tensor, dja)
        det = f_determinant(tensor)
        ainv(:) = dja(:)/det
    end subroutine s_calculate_ainverse

    ! neo-Hookean only at this time, will need to be changed later
    function f_elastic_energy(btensor, j, k, l)
        type(scalar_field), dimension(b_size), intent(IN) :: btensor
        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: invariant1, f_elastic_energy

        invariant1 = btensor(1)%sf(j,k,l)

        if (num_dims == 2) then
           invariant1 = invariant1 + btensor(3)%sf(j,k,l)
        elseif (num_dims == 3) then 
           invariant1 = invariant1 + btensor(4)%sf(j,k,l) + btensor(6)%sf(j,k,l)
        end if 
        ! compute the invariant without the elastic modulus
        f_elastic_energy = 0.5d0*(invariant1 - 3)/btensor(b_size)%sf(j, k, l)
    end function f_elastic_energy

    subroutine s_compute_grad_xi(q_prim_vf, j, k, l, grad_xi)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_dims**2), intent(INOUT) :: grad_xi
        integer, intent(IN) :: j, k, l

        ! dxix/dx
        grad_xi(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                       - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                       + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                       - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                       /(12d0*(x_cb(j) - x_cb(j-1)))
                       !/(12d0*dx(j))

        if (num_dims > 1) then
              ! dxiy / dx 
              grad_xi(2) = &
                    (q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                    - 8d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                    + 8d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                    - q_prim_vf(xibeg + 1)%sf(j + 2, k, l)) &
                    /(12d0*(x_cb(j) - x_cb(j-1)))
                    !/(12d0*dx(j))
              ! dxix / dy
              grad_xi(3) = &
                   (q_prim_vf(xibeg)%sf(j, k - 2, l) &
                    - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                    + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                    - q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                    /(12d0*(y_cb(k) - y_cb(k-1)))
                    !/(12d0*dy(k))
              ! dxiy / dy
              grad_xi(4) = &
                    (q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                    - 8d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                    + 8d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                    - q_prim_vf(xibeg + 1)%sf(j, k + 2, l)) &
                    /(12d0*(y_cb(k) - y_cb(k-1)))
                    !/(12d0*dy(k))
        end if 
        ! 3D
        if (num_dims > 2) then
               ! using results from upper if statement to map form 2x2 to 3x3 tensor
               grad_xi(5) = grad_xi(4)
               grad_xi(4) = grad_xi(3)
               ! dxix / dz
               grad_xi(3) = &
                  (q_prim_vf(xibeg)%sf(j, k, l - 2) &
                  - 8d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                  + 8d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                  - q_prim_vf(xibeg)%sf(j, k, l + 2)) &
                  /(12d0*(z_cb(l) - z_cb(l-1)))
                  !/(12d0*dz(l))
               ! dxiy / dz
               grad_xi(6) = &
                  (q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                  - 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                  + 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                  - q_prim_vf(xibeg + 1)%sf(j, k, l + 2)) &
                  /(12d0*(z_cb(l) - z_cb(l-1)))
                  !/(12d0*dz(l))
               ! dxiz / dx
               grad_xi(7) = &
                  (q_prim_vf(xiend)%sf(j - 2, k, l) &
                  - 8d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                  + 8d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                  - q_prim_vf(xiend)%sf(j + 2, k, l)) &
                  /(12d0*(x_cb(j) - x_cb(j-1)))
                  !/(12d0*dx(j))
               ! dxiz / dy
               grad_xi(8) = &
                  (q_prim_vf(xiend)%sf(j, k - 2, l) &
                  - 8d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                  + 8d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                  - q_prim_vf(xiend)%sf(j, k + 2, l)) &
                  /(12d0*(y_cb(k) - y_cb(k-1)))
                  !/(12d0*dy(k))
               ! dxiz / dz
               grad_xi(9) = &
                  (q_prim_vf(xiend)%sf(j, k, l - 2) &
                  - 8d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                  + 8d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                  - q_prim_vf(xiend)%sf(j, k, l + 2)) &
                  /(12d0*(z_cb(l) - z_cb(l-1)))
                  !/(12d0*dz(l))
       end if

    end subroutine s_compute_grad_xi

end module m_rmt_tensor_calc
