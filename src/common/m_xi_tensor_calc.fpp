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

    private; public ::  f_elastic_energy, &
 s_compute_grad_xi

contains

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
        f_elastic_energy = 0.5d0*(invariant1 - 3.0d0)/btensor(b_size)%sf(j, k, l)
    end function f_elastic_energy

    subroutine s_compute_grad_xi(q_prim_vf, j, k, l, ftensor, grad_xi, tensorb)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_dims**2+1), intent(INOUT) :: ftensor
        real(kind(0d0)), dimension(num_dims**2+1), intent(INOUT) :: grad_xi, tensorb
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: determinant
        integer :: i

        !print *, 'xibeg :: ',xibeg,', qprim :: ',q_prim_vf(xibeg)%sf(j,k,l)

        ! grad_xi definition / organization
        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
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

        !do i = 1, num_dims**2
        !   print *, 'i :: ',i,', grad_xi :: ',grad_xi(i)
        !end do 

        ! calculating the adjoint of grad_xi tensor in preparation for 
        ! calculating the inverse of the tensor
        if (num_dims == 1) then
            ftensor(1) = 1
        elseif (num_dims == 2) then
            ftensor(1) = grad_xi(4)
            ftensor(2) = -grad_xi(3)
            ftensor(3) = -grad_xi(2)
            ftensor(4) = grad_xi(1)
        elseif (num_dims == 3) then
            ftensor(1) = grad_xi(5)*grad_xi(9) - grad_xi(6)*grad_xi(8)
            ftensor(2) = -(grad_xi(2)*grad_xi(9) - grad_xi(3)*grad_xi(8))
            ftensor(3) = grad_xi(2)*grad_xi(6) - grad_xi(3)*grad_xi(5)
            ftensor(4) = -(grad_xi(4)*grad_xi(9) - grad_xi(6)*grad_xi(7))
            ftensor(5) = grad_xi(1)*grad_xi(9) - grad_xi(3)*grad_xi(7)
            ftensor(6) = -(grad_xi(1)*grad_xi(6) - grad_xi(4)*grad_xi(3))
            ftensor(7) = grad_xi(4)*grad_xi(8) - grad_xi(5)*grad_xi(7)
            ftensor(8) = -(grad_xi(1)*grad_xi(8) - grad_xi(2)*grad_xi(7))
            ftensor(9) = grad_xi(1)*grad_xi(5) - grad_xi(2)*grad_xi(4)
        end if

        ! calculating the determinant of the grad_xi tensor 
        if (num_dims == 1) then
            determinant = grad_xi(1)
        elseif (num_dims == 2) then
            determinant = grad_xi(1)*grad_xi(4) - grad_xi(2)*grad_xi(3)
        else
            determinant = grad_xi(1)*(grad_xi(5)*grad_xi(9) - grad_xi(6)*grad_xi(8)) &
                            - grad_xi(2)*(grad_xi(4)*grad_xi(9) - grad_xi(6)*grad_xi(7)) &
                            + grad_xi(3)*(grad_xi(4)*grad_xi(8) - grad_xi(5)*grad_xi(7))
        end if

        ! error checking
        if (determinant == 0) then
            if(proc_rank == 0) then
            print *, 'determinant :: ', determinant
            !print *, 'ERROR: Determinant was zero'
            !stop
            end if
        end if
        ! calculating the inverse and saving it in tensorb, which is F tensor
        tensorb(:) = ftensor(:)/determinant

        ! calculating F transpose F 
        ftensor(1) = tensorb(1)**2
        if (num_dims == 2) then
            ftensor(1) = ftensor(1) + tensorb(3)**2
            ftensor(2) = tensorb(1)*tensorb(2) + tensorb(3)*tensorb(4)
            ftensor(3) = ftensor(2)
            ftensor(4) = tensorb(2)**2 + tensorb(4)**2
        elseif (num_dims == 3) then
            ftensor(1) = ftensor(1) + tensorb(4)**2 + tensorb(7)**2
            ftensor(5) = tensorb(2) + tensorb(5)**2 + tensorb(8)**2
            ftensor(9) = tensorb(3) + tensorb(6)**2 + tensorb(9)**2
            ftensor(2) = tensorb(1)*tensorb(2) + tensorb(4)*tensorb(5) + tensorb(7)*tensorb(8)
            ftensor(3) = tensorb(1)*tensorb(3) + tensorb(4)*tensorb(6) + tensorb(7)*tensorb(9)
            ftensor(6) = tensorb(2)*tensorb(3) + tensorb(5)*tensorb(6) + tensorb(8)*tensorb(9)
            ftensor(4) = ftensor(2)
            ftensor(7) = ftensor(3)
            ftensor(8) = ftensor(4)
        end if

        ! calculating the determinant of the F tensor and storing in last entry of ftensor
        if (num_dims == 1) then ! 1D
            ftensor(num_dims**2+1) = tensorb(1)
        elseif (num_dims == 2) then ! 2D
            ftensor(num_dims**2+1) = tensorb(1)*tensorb(4) - tensorb(2)*tensorb(3)
        else ! 3D
            ftensor(num_dims**2+1) = tensorb(1)*(tensorb(5)*tensorb(9) - tensorb(6)*tensorb(8)) &
                            - tensorb(2)*(tensorb(4)*tensorb(9) - tensorb(6)*tensorb(7)) &
                            + tensorb(3)*(tensorb(4)*tensorb(8) - tensorb(5)*tensorb(7))
        end if

    end subroutine s_compute_grad_xi

end module m_xi_tensor_calc
