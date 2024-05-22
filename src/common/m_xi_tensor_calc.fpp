!>
!! @file m_xi_tensor_calc.f90
!! @brief Contains module m_xi_tensor_calc

!> @brief This module consists of subroutines used in the calculation of matrix
!!              operations for the reference map tensor

module m_xi_tensor_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================

    implicit none

    private; public :: s_compute_grad_xi, f_elastic_energy
 !s_calculate_ainverse, &
 !s_calculate_atransposea, &
 !f_determinant, &
 !s_compute_grad_xi

contains

    function f_determinant(tensor)
        !$acc routine seq
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
            print *, 'f_determinant :: ', f_determinant
            print *, 'ERROR: Determinant was zero'
            stop
        end if
    end function f_determinant

    subroutine s_calculate_atransposea(tensor, ata)
        !$acc routine seq
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
        !$acc routine seq
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
        !$acc routine seq
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ainv
        real(kind(0d0)), dimension(num_dims**2) :: dja
        real(kind(0d0)) :: det

        call s_calculate_adjointa(tensor, dja)
        ainv(:) = dja(:)/f_determinant(tensor)

    end subroutine s_calculate_ainverse

    ! neo-Hookean only at this time, will need to be changed later
    function f_elastic_energy(btensor, j, k, l)
        !$acc routine seq
        type(scalar_field), dimension(b_size), intent(IN) :: btensor
        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: invariant1, f_elastic_energy

        invariant1 = btensor(1)%sf(j, k, l)
        f_elastic_energy = 0d0
        if (num_dims == 2) then
            invariant1 = invariant1 + btensor(3)%sf(j, k, l)
        elseif (num_dims == 3) then
            invariant1 = invariant1 + btensor(4)%sf(j, k, l) + btensor(6)%sf(j, k, l)
        end if
        ! compute the invariant without the elastic modulus
        f_elastic_energy = 0.5d0*(invariant1 - 3)/btensor(b_size)%sf(j, k, l)
    end function f_elastic_energy

    subroutine s_compute_grad_xi(q_prim_vf, j, k, l, grad_xi, tensora, tensorc)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_dims**2+1), intent(OUT) :: grad_xi 
        real(kind(0d0)), dimension(num_dims**2+1), intent(OUT) :: tensora
        real(kind(0d0)), dimension(num_dims**2+1), intent(OUT) :: tensorc

        integer, intent(IN) :: j, k, l
        integer :: i

        if(j == 0) then
           ! dxix/dx
           grad_xi(1) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                      - 36d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                      + 16d0*q_prim_vf(xibeg)%sf(j + 3, k, l) &
                      -  3d0*q_prim_vf(xibeg)%sf(j + 4, k, l) ) &
                     /(12d0*(x_cb(j+1) - x_cb(j)))
        else if (j == 1) then
           ! dxix/dx
           grad_xi(1) = (-3d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                      - 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                      -  6d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                      +      q_prim_vf(xibeg)%sf(j + 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == m - 1) then
           ! dxix/dx
           grad_xi(1) = (3d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                      + 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                      +  6d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                      -      q_prim_vf(xibeg)%sf(j - 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == m) then
           ! dxix/dx
           grad_xi(1) =(25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                      + 36d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                      - 16d0*q_prim_vf(xibeg)%sf(j - 3, k, l) &
                      +  3d0*q_prim_vf(xibeg)%sf(j - 4, k, l) ) &
                     /(12d0*(x_cb(j) - x_cb(j-1)))
        else
           ! dxix/dx
           grad_xi(1) = (   q_prim_vf(xibeg)%sf(j - 2, k, l) &
                      - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                      + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                      -     q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
        end if  

        if (num_dims > 1) then
          if(j == 0) then
             ! dxiy / dx
             grad_xi(2) = (-25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg+1)%sf(j + 1, k, l) &
                      - 36d0*q_prim_vf(xibeg+1)%sf(j + 2, k, l) &
                      + 16d0*q_prim_vf(xibeg+1)%sf(j + 3, k, l) &
                      -  3d0*q_prim_vf(xibeg+1)%sf(j + 4, k, l) ) &
                     /(12d0*(x_cb(j + 1) - x_cb(j)))
          else if (j == 1) then
             ! dxiy / dx
             grad_xi(2) = (-3d0*q_prim_vf(xibeg+1)%sf(j - 1, k, l) &
                      - 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg+1)%sf(j + 1, k, l) &
                      -  6d0*q_prim_vf(xibeg+1)%sf(j + 2, k, l) &
                      +      q_prim_vf(xibeg+1)%sf(j + 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else if (j == m - 1) then
             ! dxiy / dx
             grad_xi(2) = (3d0*q_prim_vf(xibeg+1)%sf(j + 1, k, l) &
                      + 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg+1)%sf(j - 1, k, l) &
                      +  6d0*q_prim_vf(xibeg+1)%sf(j - 2, k, l) &
                      -      q_prim_vf(xibeg+1)%sf(j - 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else if (j == m) then
             ! dxiy / dx
             grad_xi(2) = (25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg+1)%sf(j - 1, k, l) &
                      + 36d0*q_prim_vf(xibeg+1)%sf(j - 2, k, l) &
                      - 16d0*q_prim_vf(xibeg+1)%sf(j - 3, k, l) &
                      +  3d0*q_prim_vf(xibeg+1)%sf(j - 4, k, l) ) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else
             ! dxiy / dx
             grad_xi(2) = ( q_prim_vf(xibeg+1)%sf(j - 2, k, l) &
                      - 8d0*q_prim_vf(xibeg+1)%sf(j - 1, k, l) &
                      + 8d0*q_prim_vf(xibeg+1)%sf(j + 1, k, l) &
                      -     q_prim_vf(xibeg+1)%sf(j + 2, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          end if  

          if(k == 0) then
             ! dxix / dy
             grad_xi(3) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                      - 36d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                      + 16d0*q_prim_vf(xibeg)%sf(j, k + 3, l) &
                      -  3d0*q_prim_vf(xibeg)%sf(j, k + 4, l) ) &
                     /(12d0*(y_cb(k+1) - y_cb(k)))
          else if (k == 1) then
             ! dxix / dy
             grad_xi(3) = (-3d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                      - 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                      -  6d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                      +      q_prim_vf(xibeg)%sf(j, k + 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n - 1) then
             ! dxix / dy
             grad_xi(3) = (3d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                      + 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                      +  6d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                      -      q_prim_vf(xibeg)%sf(j, k - 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n) then
             ! dxix / dy
             grad_xi(3) =(25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                      + 36d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                      - 16d0*q_prim_vf(xibeg)%sf(j, k - 3, l) &
                      +  3d0*q_prim_vf(xibeg)%sf(j, k - 4, l) ) &
                     /(12d0*(y_cb(j) - y_cb(j-1)))
          else
             ! dxix / dy
             grad_xi(3) = ( q_prim_vf(xibeg)%sf(j, k - 2, l) &
                      - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                      + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                      -     q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          end if  

          if(k == 0) then
             ! dxiy / dy
             grad_xi(4) = (-25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg+1)%sf(j, k + 1, l) &
                      - 36d0*q_prim_vf(xibeg+1)%sf(j, k + 2, l) &
                      + 16d0*q_prim_vf(xibeg+1)%sf(j, k + 3, l) &
                      -  3d0*q_prim_vf(xibeg+1)%sf(j, k + 4, l) ) &
                     /(12d0*(y_cb(k+1) - y_cb(k)))
          else if (k == 1) then
             ! dxiy / dy
             grad_xi(4) = (-3d0*q_prim_vf(xibeg+1)%sf(j, k - 1, l) &
                      - 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg+1)%sf(j, k + 1, l) &
                      -  6d0*q_prim_vf(xibeg+1)%sf(j, k + 2, l) &
                      +      q_prim_vf(xibeg+1)%sf(j, k + 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n - 1) then
             ! dxiy / dy
             grad_xi(4) = (3d0*q_prim_vf(xibeg+1)%sf(j, k + 1, l) &
                      + 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg+1)%sf(j, k - 1, l) &
                      +  6d0*q_prim_vf(xibeg+1)%sf(j, k - 2, l) &
                      -      q_prim_vf(xibeg+1)%sf(j, k - 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n) then
             ! dxiy / dy
             grad_xi(4) =(25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg+1)%sf(j, k - 1, l) &
                      + 36d0*q_prim_vf(xibeg+1)%sf(j, k - 2, l) &
                      - 16d0*q_prim_vf(xibeg+1)%sf(j, k - 3, l) &
                      +  3d0*q_prim_vf(xibeg+1)%sf(j, k - 4, l) ) &
                     /(12d0*(y_cb(j) - y_cb(j-1)))
          else
             ! dxiy / dy
             grad_xi(4) = ( q_prim_vf(xibeg+1)%sf(j, k - 2, l) &
                      - 8d0*q_prim_vf(xibeg+1)%sf(j, k - 1, l) &
                      + 8d0*q_prim_vf(xibeg+1)%sf(j, k + 1, l) &
                      -     q_prim_vf(xibeg+1)%sf(j, k + 2, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          end if  

        end if
        ! 3D
        if (num_dims > 2) then
            ! using results from upper if statement to map form 2x2 to 3x3 tensor
            grad_xi(5) = grad_xi(4)
            grad_xi(4) = grad_xi(3)

          if(l == 0) then
             ! dxix / dz
             grad_xi(3) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg)%sf(j, k , l + 1) &
                      - 36d0*q_prim_vf(xibeg)%sf(j, k , l + 2) &
                      + 16d0*q_prim_vf(xibeg)%sf(j, k , l + 3) &
                      -  3d0*q_prim_vf(xibeg)%sf(j, k , l + 4) ) &
                     /(12d0*(z_cb(k+1) - z_cb(k)))
          else if (l == 1) then
             ! dxix / dz
             grad_xi(3) = (-3d0*q_prim_vf(xibeg)%sf(j, k , l - 1) &
                      - 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg)%sf(j, k , l + 1) &
                      -  6d0*q_prim_vf(xibeg)%sf(j, k , l + 2) &
                      +      q_prim_vf(xibeg)%sf(j, k , l + 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p - 1) then
             ! dxix / dz
             grad_xi(3) = (3d0*q_prim_vf(xibeg)%sf(j, k , l + 1) &
                      + 10d0*q_prim_vf(xibeg)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg)%sf(j, k , l - 1) &
                      +  6d0*q_prim_vf(xibeg)%sf(j, k , l - 2) &
                      -      q_prim_vf(xibeg)%sf(j, k , l - 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p) then
             ! dxix / dz
             grad_xi(3) =(25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg)%sf(j, k , l - 1) &
                      + 36d0*q_prim_vf(xibeg)%sf(j, k , l - 2) &
                      - 16d0*q_prim_vf(xibeg)%sf(j, k , l - 3) &
                      +  3d0*q_prim_vf(xibeg)%sf(j, k , l - 4) ) &
                     /(12d0*(z_cb(j) - z_cb(j-1)))
          else
             ! dxix / dz
             grad_xi(3) = ( q_prim_vf(xibeg)%sf(j, k , l - 2) &
                      - 8d0*q_prim_vf(xibeg)%sf(j, k , l - 1) &
                      + 8d0*q_prim_vf(xibeg)%sf(j, k , l + 1) &
                      -     q_prim_vf(xibeg)%sf(j, k , l + 2)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          end if  

          if(l == 0) then
             ! dxiy / dz
             grad_xi(6) = (-25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xibeg+1)%sf(j, k , l + 1) &
                      - 36d0*q_prim_vf(xibeg+1)%sf(j, k , l + 2) &
                      + 16d0*q_prim_vf(xibeg+1)%sf(j, k , l + 3) &
                      -  3d0*q_prim_vf(xibeg+1)%sf(j, k , l + 4) ) &
                     /(12d0*(z_cb(k+1) - z_cb(k)))
          else if (l == 1) then
             ! dxiy / dz
             grad_xi(6) = (-3d0*q_prim_vf(xibeg+1)%sf(j, k , l - 1) &
                      - 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xibeg+1)%sf(j, k , l + 1) &
                      -  6d0*q_prim_vf(xibeg+1)%sf(j, k , l + 2) &
                      +      q_prim_vf(xibeg+1)%sf(j, k , l + 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p - 1) then
             ! dxiy / dz
             grad_xi(6) = (3d0*q_prim_vf(xibeg+1)%sf(j, k , l + 1) &
                      + 10d0*q_prim_vf(xibeg+1)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xibeg+1)%sf(j, k , l - 1) &
                      +  6d0*q_prim_vf(xibeg+1)%sf(j, k , l - 2) &
                      -      q_prim_vf(xibeg+1)%sf(j, k , l - 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p) then
             ! dxiy / dz
             grad_xi(6) =(25d0*q_prim_vf(xibeg+1)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xibeg+1)%sf(j, k , l - 1) &
                      + 36d0*q_prim_vf(xibeg+1)%sf(j, k , l - 2) &
                      - 16d0*q_prim_vf(xibeg+1)%sf(j, k , l - 3) &
                      +  3d0*q_prim_vf(xibeg+1)%sf(j, k , l - 4) ) &
                     /(12d0*(z_cb(j) - z_cb(j-1)))
          else
             ! dxiy / dz
             grad_xi(6) = ( q_prim_vf(xibeg+1)%sf(j, k , l - 2) &
                      - 8d0*q_prim_vf(xibeg+1)%sf(j, k , l - 1) &
                      + 8d0*q_prim_vf(xibeg+1)%sf(j, k , l + 1) &
                      -     q_prim_vf(xibeg+1)%sf(j, k , l + 2)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          end if  

          if(j == 0) then
             ! dxiz / dx
             grad_xi(7) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                      - 36d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                      + 16d0*q_prim_vf(xiend)%sf(j + 3, k, l) &
                      -  3d0*q_prim_vf(xiend)%sf(j + 4, k, l) ) &
                     /(12d0*(x_cb(j + 1) - x_cb(j)))
          else if (j == 1) then
             ! dxiz / dx
             grad_xi(7) = (-3d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                      - 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                      -  6d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                      +      q_prim_vf(xiend)%sf(j + 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else if (j == m - 1) then
             ! dxiz / dx
             grad_xi(7) = (3d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                      + 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                      +  6d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                      -      q_prim_vf(xiend)%sf(j - 3, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else if (j == m) then
             ! dxiz / dx
             grad_xi(7) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                      + 36d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                      - 16d0*q_prim_vf(xiend)%sf(j - 3, k, l) &
                      +  3d0*q_prim_vf(xiend)%sf(j - 4, k, l) ) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          else
             ! dxiz / dx
             grad_xi(7) = ( q_prim_vf(xiend)%sf(j - 2, k, l) &
                      - 8d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                      + 8d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                      -     q_prim_vf(xiend)%sf(j + 2, k, l)) &
                     /(12d0*(x_cb(j) - x_cb(j - 1)))
          end if  

          if(k == 0) then
             ! dxiz / dy
             grad_xi(8) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                      - 36d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                      + 16d0*q_prim_vf(xiend)%sf(j, k + 3, l) &
                      -  3d0*q_prim_vf(xiend)%sf(j, k + 4, l) ) &
                     /(12d0*(y_cb(k+1) - y_cb(k)))
          else if (k == 1) then
             ! dxiz / dy
             grad_xi(8) = (-3d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                      - 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                      -  6d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                      +      q_prim_vf(xiend)%sf(j, k + 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n - 1) then
             ! dxiz / dy
             grad_xi(8) = (3d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                      + 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                      +  6d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                      -      q_prim_vf(xiend)%sf(j, k - 3, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          else if (k == n) then
             ! dxiz / dy
             grad_xi(8) =(25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                      + 36d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                      - 16d0*q_prim_vf(xiend)%sf(j, k - 3, l) &
                      +  3d0*q_prim_vf(xiend)%sf(j, k - 4, l) ) &
                     /(12d0*(y_cb(j) - y_cb(j-1)))
          else
             ! dxiz / dy
             grad_xi(8) = ( q_prim_vf(xiend)%sf(j, k - 2, l) &
                      - 8d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                      + 8d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                      -     q_prim_vf(xiend)%sf(j, k + 2, l)) &
                     /(12d0*(y_cb(j) - y_cb(j - 1)))
          end if  

          if(l == 0) then
             ! dxiz / dz
             grad_xi(9) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      + 48d0*q_prim_vf(xiend)%sf(j, k , l + 1) &
                      - 36d0*q_prim_vf(xiend)%sf(j, k , l + 2) &
                      + 16d0*q_prim_vf(xiend)%sf(j, k , l + 3) &
                      -  3d0*q_prim_vf(xiend)%sf(j, k , l + 4) ) &
                     /(12d0*(z_cb(k+1) - z_cb(k)))
          else if (l == 1) then
             ! dxiz / dz
             grad_xi(9) = (-3d0*q_prim_vf(xiend)%sf(j, k , l - 1) &
                      - 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      + 18d0*q_prim_vf(xiend)%sf(j, k , l + 1) &
                      -  6d0*q_prim_vf(xiend)%sf(j, k , l + 2) &
                      +      q_prim_vf(xiend)%sf(j, k , l + 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p - 1) then
             ! dxiz / dz
             grad_xi(9) = (3d0*q_prim_vf(xiend)%sf(j, k , l + 1) &
                      + 10d0*q_prim_vf(xiend)%sf(j,k,l) &
                      - 18d0*q_prim_vf(xiend)%sf(j, k , l - 1) &
                      +  6d0*q_prim_vf(xiend)%sf(j, k , l - 2) &
                      -      q_prim_vf(xiend)%sf(j, k , l - 3)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          else if (l == p) then
             ! dxiz / dz
             grad_xi(9) =(25d0*q_prim_vf(xiend)%sf(j, k, l) &
                      - 48d0*q_prim_vf(xiend)%sf(j, k , l - 1) &
                      + 36d0*q_prim_vf(xiend)%sf(j, k , l - 2) &
                      - 16d0*q_prim_vf(xiend)%sf(j, k , l - 3) &
                      +  3d0*q_prim_vf(xiend)%sf(j, k , l - 4) ) &
                     /(12d0*(z_cb(j) - z_cb(j-1)))
          else
             ! dxiz / dz
             grad_xi(9) = ( q_prim_vf(xiend)%sf(j, k , l - 2) &
                      - 8d0*q_prim_vf(xiend)%sf(j, k , l - 1) &
                      + 8d0*q_prim_vf(xiend)%sf(j, k , l + 1) &
                      -     q_prim_vf(xiend)%sf(j, k , l + 2)) &
                     /(12d0*(z_cb(j) - z_cb(j - 1)))
          end if  
        end if

    !if(proc_rank == 0) then
    !    do i = 1, num_dims**2
    !        print *, "i :: ",i,", grad_xi :: ",grad_xi(i)
    !    end do
    !end if

    end subroutine s_compute_grad_xi

end module m_xi_tensor_calc
