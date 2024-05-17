    subroutine s_compute_grad_xi(q_prim_vf, j, k, l, grad_xi)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(num_dims**2), intent(INOUT) :: grad_xi
        integer, intent(IN) :: j, k, l

        ! dxix/dx
        grad_xi(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                       - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                       + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                       - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                       /(12d0*dx(k))

        if (num_dims > 1) then
              ! dxiy / dx 
              grad_xi(2) = &
                    (q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                    - 8d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                    + 8d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                    - q_prim_vf(xibeg + 1)%sf(j + 2, k, l)) &
                    /(12d0*dx(k))
              ! dxix / dy
              grad_xi(3) = &
                   (q_prim_vf(xibeg)%sf(j, k - 2, l) &
                    - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                    + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                    - q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                    /(12d0*dy(l))
              ! dxiy / dy
              grad_xi(4) = &
                    (q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                    - 8d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                    + 8d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                    - q_prim_vf(xibeg + 1)%sf(j, k + 2, l)) &
                    /(12d0*dy(l))
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
                  /(12d0*dz(q))
               ! dxiy / dz
               grad_xi(6) = &
                  (q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                  - 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                  + 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                  - q_prim_vf(xibeg + 1)%sf(j, k, l + 2)) &
                  /(12d0*dz(q))
               ! dxiz / dx
               grad_xi(7) = &
                  (q_prim_vf(xiend)%sf(j - 2, k, l) &
                  - 8d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                  + 8d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                  - q_prim_vf(xiend)%sf(j + 2, k, l)) &
                  /(12d0*dx(k))
               ! dxiz / dy
               grad_xi(8) = &
                  (q_prim_vf(xiend)%sf(j, k - 2, l) &
                  - 8d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                  + 8d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                  - q_prim_vf(xiend)%sf(j, k + 2, l)) &
                  /(12d0*dy(l))
               ! dxiz / dz
               grad_xi(9) = &
                  (q_prim_vf(xiend)%sf(j, k, l - 2) &
                  - 8d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                  + 8d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                  - q_prim_vf(xiend)%sf(j, k, l + 2)) &
                  /(12d0*dz(q))
       end if
    end subroutine s_compute_grad_xi

end module m_rmt_tensor_calc
