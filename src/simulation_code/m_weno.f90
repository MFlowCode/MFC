!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names 
!!  are listed in the copyright file included with this source 
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published 
!!  by the Free Software Foundation, either version 3 of the license 
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!  
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).  
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_weno.f90
!! @brief Contains module m_weno
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  Weighted essentially non-oscillatory (WENO) reconstruction scheme
!!              that is supplemented with monotonicity preserving bounds (MPWENO)
!!              and a mapping function that boosts the accuracy of the non-linear
!!              weights (WENOM). MPWENO, see Balsara and Shu (2000), prevents the
!!              reconstructed values to lay outside the range set by the stencil,
!!              while WENOM, see Henrick et al. (2005), recovers the formal order
!!              of accuracy of the reconstruction at critical points. Please note
!!              that the basic WENO approach is implemented according to the work
!!              of Jiang and Shu (1996).
MODULE m_weno
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_weno_module, s_weno, s_finalize_weno_module
    
    !> @name The cell-average variables that will be WENO-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the WENO procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed WENO-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{
    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: v_rs_wsL, v_rs_wsR
    !> @}
  
    !> @name Left and right WENO-reconstructed values of the cell-average variables.
    !! Note that the reshaped property of the variables from which these were
    !! obtained, v_rs_wsL and v_rs_wsR, is initially kept. Once the reshaping
    !! is undone, the reconstructed values are moved into vL_vf and vR_vf.
    !> @{
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: vL_rs_vf, vR_rs_vf 
    !> @}
     
    
    ! WENO Coefficients ========================================================
   

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at
    !! the left and right quadrature points (QP), in the x-, y- and z-directions.
    !! Note that the first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the cell-location in the relevant coordinate direction.
    !> @{
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbL_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbL_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbL_z
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpL_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpL_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpL_z
    
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbR_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbR_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_cbR_z
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpR_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpR_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: poly_coef_qpR_z
    
    REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: poly_coef_L => NULL()
    REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: poly_coef_R => NULL()
    !> @}

    !> @name The ideal weights at the left and the right cell-boundaries and at the
    !! left and the right quadrature points, in x-, y- and z-directions. Note
    !! that the first dimension of the array identifies the weight, while the
    !! last denotes the cell-location in the relevant coordinate direction.
    !> @{
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbL_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbL_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbL_z
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpL_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpL_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpL_z
    
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbR_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbR_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_cbR_z
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpR_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpR_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:) :: d_qpR_z
    
    REAL(KIND(0d0)), POINTER, DIMENSION(:,:) :: d_L => NULL()
    REAL(KIND(0d0)), POINTER, DIMENSION(:,:) :: d_R => NULL()
    !> @}

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note
    !! that the first array dimension identifies the smoothness indicator, the
    !! second identifies the position of its coefficients and the last denotes
    !! the cell-location in the relevant coordinate direction.
    !> @{
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: beta_coef_x
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: beta_coef_y
    REAL(KIND(0d0)), TARGET , ALLOCATABLE, DIMENSION(:,:,:) :: beta_coef_z
    
    REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: beta_coef => NULL()
    !> @}

    ! END: WENO Coefficients ===================================================
    
    

    INTEGER :: v_size !< Number of WENO-reconstructed cell-average variables
    
    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    TYPE(bounds_info) :: is1,is2,is3
    !> @}
    
    CONTAINS
        
        
        
        
        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.        
        SUBROUTINE s_initialize_weno_module() ! --------------------------------

            TYPE(bounds_info) :: ix,iy,iz !< Indical bounds in the x-, y- and z-directions
            
            
            IF(weno_order == 1) RETURN
            
            
            ! Allocating WENO-stencil for the variables to be WENO-reconstructed
            ALLOCATE(v_rs_wsL(-weno_polyn:weno_polyn))
            ALLOCATE(v_rs_wsR(-weno_polyn:weno_polyn))
            
            
            ! Allocating/Computing WENO Coefficients in x-direction ============
            ix%beg = -buff_size + weno_polyn; ix%end = m - ix%beg
            
            ALLOCATE(poly_coef_cbL_x(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      ix%beg :    ix%end    ))
            ALLOCATE(poly_coef_cbR_x(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      ix%beg :    ix%end    ))
            
            ALLOCATE(d_cbL_x(0:weno_polyn,ix%beg:ix%end))
            ALLOCATE(d_cbR_x(0:weno_polyn,ix%beg:ix%end))
            
            ALLOCATE(beta_coef_x(    0   :    weno_polyn   , &
                                     0   : 2*(weno_polyn-1), &
                                  ix%beg :      ix%end      ))
            
            CALL s_compute_weno_coefficients(1,1,ix)
            
            IF(commute_err .OR. split_err) THEN
                
                ALLOCATE(poly_coef_qpL_x(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          ix%beg :    ix%end    ))
                ALLOCATE(poly_coef_qpR_x(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          ix%beg :    ix%end    ))
                
                ALLOCATE(d_qpL_x(0:weno_polyn,ix%beg:ix%end))
                ALLOCATE(d_qpR_x(0:weno_polyn,ix%beg:ix%end))
                
                CALL s_compute_weno_coefficients(1,2,ix)
                
            END IF
            ! ==================================================================
            
            
            ! Allocating/Computing WENO Coefficients in y-direction ============
            IF(n == 0) RETURN
            
            iy%beg = -buff_size + weno_polyn; iy%end = n - iy%beg
            
            ALLOCATE(poly_coef_cbL_y(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      iy%beg :    iy%end    ))
            ALLOCATE(poly_coef_cbR_y(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      iy%beg :    iy%end    ))
            
            ALLOCATE(d_cbL_y(0:weno_polyn,iy%beg:iy%end))
            ALLOCATE(d_cbR_y(0:weno_polyn,iy%beg:iy%end))
            
            ALLOCATE(beta_coef_y(    0   :    weno_polyn   , &
                                     0   : 2*(weno_polyn-1), &
                                  iy%beg :      iy%end      ))
            
            CALL s_compute_weno_coefficients(2,1,iy)
            
            IF(commute_err .OR. split_err) THEN
                
                ALLOCATE(poly_coef_qpL_y(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          iy%beg :    iy%end    ))
                ALLOCATE(poly_coef_qpR_y(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          iy%beg :    iy%end    ))
                
                ALLOCATE(d_qpL_y(0:weno_polyn,iy%beg:iy%end))
                ALLOCATE(d_qpR_y(0:weno_polyn,iy%beg:iy%end))
                
                CALL s_compute_weno_coefficients(2,2,iy)
                
            END IF
            ! ==================================================================
            
            
            ! Allocating/Computing WENO Coefficients in z-direction ============
            IF(p == 0) RETURN
            
            iz%beg = -buff_size + weno_polyn; iz%end = p - iz%beg
            
            ALLOCATE(poly_coef_cbL_z(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      iz%beg :    iz%end    ))
            ALLOCATE(poly_coef_cbR_z(    0   :  weno_polyn , &
                                         0   : weno_polyn-1, &
                                      iz%beg :    iz%end    ))
            
            ALLOCATE(d_cbL_z(0:weno_polyn,iz%beg:iz%end))
            ALLOCATE(d_cbR_z(0:weno_polyn,iz%beg:iz%end))
            
            ALLOCATE(beta_coef_z(    0   :    weno_polyn   , &
                                     0   : 2*(weno_polyn-1), &
                                  iz%beg :      iz%end      ))
            
            CALL s_compute_weno_coefficients(3,1,iz)
            
            IF(commute_err .OR. split_err) THEN
                
                ALLOCATE(poly_coef_qpL_z(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          iz%beg :    iz%end    ))
                ALLOCATE(poly_coef_qpR_z(    0   :  weno_polyn , &
                                             0   : weno_polyn-1, &
                                          iz%beg :    iz%end    ))
                
                ALLOCATE(d_qpL_z(0:weno_polyn,iz%beg:iz%end))
                ALLOCATE(d_qpR_z(0:weno_polyn,iz%beg:iz%end))
                
                CALL s_compute_weno_coefficients(3,2,iz)
                
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_initialize_weno_module ! ------------------------------
        
        
        
        !>  The purpose of this subroutine is to compute the grid
        !!      dependent coefficients of the WENO polynomials, ideal
        !!      weights and smoothness indicators, provided the order,
        !!      the coordinate direction and the location of the WENO
        !!      reconstruction.
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param weno_loc Location of the WENO reconstruction
        !! @param is Index bounds in the s-direction
        SUBROUTINE s_compute_weno_coefficients(weno_dir, weno_loc, is) ! -------
            
            
            INTEGER, INTENT(IN) :: weno_dir, weno_loc
            TYPE(bounds_info), INTENT(IN) :: is
            INTEGER :: s
            

            REAL(KIND(0d0)), POINTER, DIMENSION(:) :: s_cb => NULL() !<
            !! Cell-boundary locations in the s-direction
            

            TYPE(bounds_info) :: bc_s !< Boundary conditions (BC) in the s-direction
            

            INTEGER :: i !< Generic loop iterator
            
            
            ! Associating WENO coefficients pointers
            CALL s_associate_weno_coefficients_pointers(weno_dir, weno_loc)
            
            
            ! Determining the number of cells, the cell-boundary locations and
            ! the boundary conditions in the coordinate direction selected for
            ! the WENO reconstruction
            IF(weno_dir == 1) THEN
                s = m; s_cb => x_cb; bc_s = bc_x
            ELSEIF(weno_dir == 2) THEN
                s = n; s_cb => y_cb; bc_s = bc_y
            ELSE
                s = p; s_cb => z_cb; bc_s = bc_z
            END IF
            
            
            ! Computing WENO3 Coefficients =====================================
            IF(weno_order == 3) THEN
                
                IF(weno_loc == 1) THEN              ! Cell-boundaries
                    
                    DO i = is%beg-1, is%end-1
                        
                        poly_coef_R(0,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                               (s_cb( i )-s_cb(i+2))
                        poly_coef_R(1,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                               (s_cb(i-1)-s_cb(i+1))
                        
                        poly_coef_L(0,0,i+1) = -poly_coef_R(0,0,i+1)
                        poly_coef_L(1,0,i+1) = -poly_coef_R(1,0,i+1)
                        
                        d_R(0,i+1) = (s_cb(i-1)-s_cb(i+1)) / &
                                     (s_cb(i-1)-s_cb(i+2))
                        d_L(0,i+1) = (s_cb(i-1)-s_cb( i )) / &
                                     (s_cb(i-1)-s_cb(i+2))
                        
                        d_R(1,i+1) = 1d0 - d_R(0,i+1)
                        d_L(1,i+1) = 1d0 - d_L(0,i+1)
                        
                        beta_coef(0,0,i+1) = 4d0*(s_cb( i )-s_cb(i+1))**2d0 / &
                                                 (s_cb( i )-s_cb(i+2))**2d0
                        beta_coef(1,0,i+1) = 4d0*(s_cb( i )-s_cb(i+1))**2d0 / &
                                                 (s_cb(i-1)-s_cb(i+1))**2d0
                        
                    END DO
                    
                ELSE                                ! Quadrature points
                    
                    ! Only used for higher-order reconstruction
                    DO i = is%beg-1, is%end-1
                        
                        poly_coef_R(0,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                               (s_cb( i )-s_cb(i+2)) / &
                                                SQRT(3d0)
                        poly_coef_R(1,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                               (s_cb(i-1)-s_cb(i+1)) / &
                                                SQRT(3d0)
                        
                        poly_coef_L(0,0,i+1) = -poly_coef_R(0,0,i+1)
                        poly_coef_L(1,0,i+1) = -poly_coef_R(1,0,i+1)
                        
                        d_R(0,i+1) = &
                            ((s_cb( i )-s_cb(i-1))+(s_cb(i+1)-s_cb(i-1))) / &
                            ((s_cb(i+2)-s_cb(i-1))+(s_cb(i+2)-s_cb(i-1)))
                        
                        d_R(1,i+1) = 1d0 - d_R(0,i+1)
                        
                        d_L(0,i+1) = d_R(0,i+1); d_L(1,i+1) = d_R(1,i+1)
                        
                    END DO
                    
                END IF
                
                ! Modifying the ideal weights coefficients in the neighborhood
                ! of beginning and end Riemann state extrapolation BC to avoid
                ! any contributions from outside of the physical domain during
                ! the WENO reconstruction
                IF (null_weights) THEN
                    IF(bc_s%beg == -4) THEN
                        d_R(1,0) = 0d0; d_R(0,0) = 1d0
                        d_L(1,0) = 0d0; d_L(0,0) = 1d0
                    END IF
                    
                    IF(bc_s%end == -4) THEN
                        d_R(0,s) = 0d0; d_R(1,s) = 1d0
                        d_L(0,s) = 0d0; d_L(1,s) = 1d0
                    END IF
                END IF
            ! END: Computing WENO3 Coefficients ================================
            
            
            ! Computing WENO5 Coefficients =====================================
            ELSE
                
                IF(weno_loc == 1) THEN              ! Cell-boundaries
                    
                    DO i = is%beg-1, is%end-1
                        
                        poly_coef_R(0,0,i+1) = &
                            ((s_cb( i )-s_cb(i+1))*(s_cb(i+1)-s_cb(i+2))) / &
                            ((s_cb( i )-s_cb(i+3))*(s_cb(i+3)-s_cb(i+1)))
                        poly_coef_R(1,0,i+1) = &
                            ((s_cb(i-1)-s_cb(i+1))*(s_cb(i+1)-s_cb( i ))) / &
                            ((s_cb(i-1)-s_cb(i+2))*(s_cb(i+2)-s_cb( i )))
                        poly_coef_R(1,1,i+1) = &
                            ((s_cb( i )-s_cb(i+1))*(s_cb(i+1)-s_cb(i+2))) / &
                            ((s_cb(i-1)-s_cb(i+1))*(s_cb(i-1)-s_cb(i+2)))
                        poly_coef_R(2,1,i+1) = &
                            ((s_cb( i )-s_cb(i+1))*(s_cb(i+1)-s_cb(i-1))) / &
                            ((s_cb(i-2)-s_cb( i ))*(s_cb(i-2)-s_cb(i+1)))
                        poly_coef_L(0,0,i+1) = &
                            ((s_cb(i+1)-s_cb( i ))*(s_cb( i )-s_cb(i+2))) / &
                            ((s_cb( i )-s_cb(i+3))*(s_cb(i+3)-s_cb(i+1)))
                        poly_coef_L(1,0,i+1) = &
                            ((s_cb( i )-s_cb(i-1))*(s_cb( i )-s_cb(i+1))) / &
                            ((s_cb(i-1)-s_cb(i+2))*(s_cb( i )-s_cb(i+2)))
                        poly_coef_L(1,1,i+1) = &
                            ((s_cb(i+1)-s_cb( i ))*(s_cb( i )-s_cb(i+2))) / &
                            ((s_cb(i-1)-s_cb(i+1))*(s_cb(i-1)-s_cb(i+2)))
                        poly_coef_L(2,1,i+1) = &
                            ((s_cb(i-1)-s_cb( i ))*(s_cb( i )-s_cb(i+1))) / &
                            ((s_cb(i-2)-s_cb( i ))*(s_cb(i-2)-s_cb(i+1)))
                        
                        poly_coef_R(0,1,i+1) = &
                            ((s_cb( i )-s_cb(i+2))+(s_cb(i+1)-s_cb(i+3))) / &
                            ((s_cb( i )-s_cb(i+2))*(s_cb( i )-s_cb(i+3))) * &
                            ((s_cb( i )-s_cb(i+1)))
                        poly_coef_R(2,0,i+1) = &
                            ((s_cb(i-2)-s_cb(i+1))+(s_cb(i-1)-s_cb(i+1))) / &
                            ((s_cb(i-1)-s_cb(i+1))*(s_cb(i+1)-s_cb(i-2))) * &
                            ((s_cb(i+1)-s_cb( i )))
                        poly_coef_L(0,1,i+1) = &
                            ((s_cb( i )-s_cb(i+2))+(s_cb( i )-s_cb(i+3))) / &
                            ((s_cb( i )-s_cb(i+2))*(s_cb( i )-s_cb(i+3))) * &
                            ((s_cb(i+1)-s_cb( i )))
                        poly_coef_L(2,0,i+1) = &
                            ((s_cb(i-2)-s_cb( i ))+(s_cb(i-1)-s_cb(i+1))) / &
                            ((s_cb(i-2)-s_cb(i+1))*(s_cb(i+1)-s_cb(i-1))) * &
                            ((s_cb( i )-s_cb(i+1)))
                        
                        d_R(0,i+1) = &
                            ((s_cb(i-2)-s_cb(i+1))*(s_cb(i+1)-s_cb(i-1))) / &
                            ((s_cb(i-2)-s_cb(i+3))*(s_cb(i+3)-s_cb(i-1)))
                        d_R(2,i+1) = &
                            ((s_cb(i+1)-s_cb(i+2))*(s_cb(i+1)-s_cb(i+3))) / &
                            ((s_cb(i-2)-s_cb(i+2))*(s_cb(i-2)-s_cb(i+3)))
                        d_L(0,i+1) = &
                            ((s_cb(i-2)-s_cb( i ))*(s_cb( i )-s_cb(i-1))) / &
                            ((s_cb(i-2)-s_cb(i+3))*(s_cb(i+3)-s_cb(i-1)))
                        d_L(2,i+1) = &
                            ((s_cb( i )-s_cb(i+2))*(s_cb( i )-s_cb(i+3))) / &
                            ((s_cb(i-2)-s_cb(i+2))*(s_cb(i-2)-s_cb(i+3)))
                        
                        d_R(1,i+1) = 1d0 - d_R(0,i+1) - d_R(2,i+1)
                        d_L(1,i+1) = 1d0 - d_L(0,i+1) - d_L(2,i+1)
                        
                        beta_coef(0,0,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(10d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+(s_cb(i+1)-s_cb(i))*(s_cb(i+2)     - &
                            s_cb(i+1))+(s_cb(i+2)-s_cb(i+1))**2d0)/((s_cb(i) - &
                            s_cb(i+3))**2d0*(s_cb(i+1)-s_cb(i+3))**2d0)
                        beta_coef(0,1,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(19d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0-(s_cb(i+1)-s_cb(i))*(s_cb(i+3)     - &
                            s_cb(i+1))+2d0*(s_cb(i+2)-s_cb(i))*((s_cb(i+2)   - &
                            s_cb(i))+(s_cb(i+3)-s_cb(i+1))))/((s_cb(i)       - &
                            s_cb(i+2))*(s_cb(i)-s_cb(i+3))**2d0*(s_cb(i+3)   - &
                            s_cb(i+1)))
                        beta_coef(0,2,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(10d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+(s_cb(i+1)-s_cb(i))*((s_cb(i+2)    - &
                            s_cb(i))+(s_cb(i+3)-s_cb(i+1)))+((s_cb(i+2)      - &
                            s_cb(i))+(s_cb(i+3)-s_cb(i+1)))**2d0)/((s_cb(i)  - &
                            s_cb(i+2))**2d0*(s_cb(i)-s_cb(i+3))**2d0)
                        beta_coef(1,0,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(10d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+(s_cb(i)-s_cb(i-1))**2d0+(s_cb(i)  - &
                            s_cb(i-1))*(s_cb(i+1)-s_cb(i)))/((s_cb(i-1)      - &
                            s_cb(i+2))**2d0*(s_cb(i)-s_cb(i+2))**2d0)
                        beta_coef(1,1,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*((s_cb(i)           - &
                            s_cb(i+1))*((s_cb(i)-s_cb(i-1))+20d0*(s_cb(i+1)  - &
                            s_cb(i)))+(2d0*(s_cb(i)-s_cb(i-1))+(s_cb(i+1)    - &
                            s_cb(i)))*(s_cb(i+2)-s_cb(i)))/((s_cb(i+1)       - &
                            s_cb(i-1))*(s_cb(i-1)-s_cb(i+2))**2d0*(s_cb(i+2) - &
                            s_cb(i)))
                        beta_coef(1,2,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(10d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+(s_cb(i+1)-s_cb(i))*(s_cb(i+2)     - &
                            s_cb(i+1))+(s_cb(i+2)-s_cb(i+1))**2d0)           / &
                            ((s_cb(i-1)-s_cb(i+1))**2d0*(s_cb(i-1)           - &
                            s_cb(i+2))**2d0)
                        beta_coef(2,0,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(12d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+((s_cb(i)-s_cb(i-2))+(s_cb(i)      - &
                            s_cb(i-1)))**2d0+3d0*((s_cb(i)-s_cb(i-2))        + &
                            (s_cb(i)-s_cb(i-1)))*(s_cb(i+1)-s_cb(i)))        / &
                            ((s_cb(i-2)-s_cb(i+1))**2d0*(s_cb(i-1)           - &
                            s_cb(i+1))**2d0)
                        beta_coef(2,1,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(19d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+((s_cb(i)-s_cb(i-2))*(s_cb(i)      - &
                            s_cb(i+1)))+2d0*(s_cb(i+1)-s_cb(i-1))*((s_cb(i)  - &
                            s_cb(i-2))+(s_cb(i+1)-s_cb(i-1))))/((s_cb(i-2)   - &
                            s_cb(i))*(s_cb(i-2)-s_cb(i+1))**2d0*(s_cb(i+1)   - &
                            s_cb(i-1)))
                        beta_coef(2,2,i+1) = &
                            4d0*(s_cb(i)-s_cb(i+1))**2d0*(10d0*(s_cb(i+1)    - &
                            s_cb(i))**2d0+(s_cb(i)-s_cb(i-1))**2d0+(s_cb(i)  - &
                            s_cb(i-1))*(s_cb(i+1)-s_cb(i)))/((s_cb(i-2)      - &
                            s_cb(i))**2d0*(s_cb(i-2)-s_cb(i+1))**2d0)
                        
                    END DO
                    
                ELSE                                ! Quadrature points
                    ! Only used for higher-order reconstruction
                    
                    DO i = is%beg-1, is%end-1
                        
                        poly_coef_R(0,0,i+1) = &
                            ((s_cb( i )-s_cb(i+2))+(s_cb(i+1)-s_cb(i+2))) / &
                            ((s_cb( i )-s_cb(i+3))*(s_cb(i+3)-s_cb(i+1))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        poly_coef_R(1,0,i+1) = &
                            ((s_cb(i-1)-s_cb( i ))+(s_cb(i-1)-s_cb(i+1))) / &
                            ((s_cb(i-1)-s_cb(i+2))*(s_cb( i )-s_cb(i+2))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        poly_coef_R(1,1,i+1) = &
                            ((s_cb( i )-s_cb(i+2))+(s_cb(i+1)-s_cb(i+2))) / &
                            ((s_cb(i-1)-s_cb(i+1))*(s_cb(i-1)-s_cb(i+2))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        poly_coef_R(2,1,i+1) = &
                            ((s_cb( i )-s_cb(i-1))+(s_cb(i+1)-s_cb(i-1))) / &
                            ((s_cb(i-2)-s_cb( i ))*(s_cb(i-2)-s_cb(i+1))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        
                        poly_coef_R(0,1,i+1) = &
                            ((s_cb( i )-s_cb(i+2))+(s_cb( i )-s_cb(i+3))  + &
                             (s_cb( i )-s_cb(i+2))+(s_cb(i+1)-s_cb(i+3))) / &
                            ((s_cb( i )-s_cb(i+2))*(s_cb( i )-s_cb(i+3))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        poly_coef_R(2,0,i+1) = &
                            ((s_cb(i-2)-s_cb(i+1))+(s_cb(i-2)-s_cb(i+1))  + &
                             (s_cb(i-1)-s_cb( i ))+(s_cb(i-1)-s_cb(i+1))) / &
                            ((s_cb(i-2)-s_cb(i+1))*(s_cb(i-1)-s_cb(i+1))) * &
                            ((s_cb( i )-s_cb(i+1)))/(2d0*SQRT(3d0))
                        
                        poly_coef_L(0,0,i+1) = -poly_coef_R(0,0,i+1)
                        poly_coef_L(0,1,i+1) = -poly_coef_R(0,1,i+1)
                        poly_coef_L(1,0,i+1) = -poly_coef_R(1,0,i+1)
                        poly_coef_L(1,1,i+1) = -poly_coef_R(1,1,i+1)
                        poly_coef_L(2,0,i+1) = -poly_coef_R(2,0,i+1)
                        poly_coef_L(2,1,i+1) = -poly_coef_R(2,1,i+1)
                        
                        d_R(0,i+1) = &
                            ((s_cb(i+1)-s_cb(i))*(18d0*(s_cb(i-2)-s_cb(i))   * &
                            (s_cb(i)-s_cb(i-1))-6d0*((s_cb(i)-s_cb(i-2))     + &
                            (s_cb(i)-s_cb(i-1)))*(s_cb(i+1)-s_cb(i))         + &
                            SQRT(3d0)*(s_cb(i+1)-s_cb(i))**2d0)              + &
                            6d0*(6d0*(s_cb(i)-s_cb(i-2))*(s_cb(i)-s_cb(i-1)) + &
                            3d0*((s_cb(i)-s_cb(i-2))+(s_cb(i)-s_cb(i-1)))    * &
                            (s_cb(i+1)-s_cb(i))+(s_cb(i+1)-s_cb(i))**2d0)    * &
                            (s_cb(i+2)-s_cb(i)))/(18d0*((s_cb(i)-s_cb(i+2))  + &
                            (s_cb(i+1)-s_cb(i+2)))*(s_cb(i-1)-s_cb(i+3))     * &
                            (s_cb(i+3)-s_cb(i-2)))
                        d_R(2,i+1) = &
                            ((6d0+SQRT(3d0))*(s_cb(i)-s_cb(i+1))**3d0        + &
                            36d0*(s_cb(i)-s_cb(i-1))*(s_cb(i+2)-s_cb(i))     * &
                            (s_cb(i+3)-s_cb(i+1))-6d0*(s_cb(i+1)             - &
                            s_cb(i))**2d0*(2d0*(s_cb(i)-s_cb(i-1))           - &
                            2d0*(s_cb(i+2)-s_cb(i))+(s_cb(i+3)-s_cb(i+1)))   + &
                            18d0*(s_cb(i+1)-s_cb(i))*((s_cb(i)-s_cb(i-1))    * &
                            (s_cb(i+2)-s_cb(i))+((s_cb(i-1)-s_cb(i))         + &
                            (s_cb(i+2)-s_cb(i)))*(s_cb(i+3)-s_cb(i+1))))     / &
                            (18d0*((s_cb(i)-s_cb(i-1))+(s_cb(i+1)            - &
                            s_cb(i-1)))*(s_cb(i-2)-s_cb(i+2))*(s_cb(i-2)     - &
                            s_cb(i+3)))
                        d_L(0,i+1) = &
                            ((s_cb(i)-s_cb(i+1))*(18d0*(s_cb(i)-s_cb(i-2))   * &
                            (s_cb(i)-s_cb(i-1))+6d0*((s_cb(i)-s_cb(i-2))     + &
                            (s_cb(i)-s_cb(i-1)))*(s_cb(i+1)-s_cb(i))         + &
                            SQRT(3d0)*(s_cb(i+1)-s_cb(i))**2d0)              + &
                            6d0*(6d0*(s_cb(i)-s_cb(i-2))*(s_cb(i)-s_cb(i-1)) + &
                            3d0*((s_cb(i)-s_cb(i-2))+(s_cb(i)-s_cb(i-1)))    * &
                            (s_cb(i+1)-s_cb(i))+(s_cb(i+1)-s_cb(i))**2d0)    * &
                            (s_cb(i+2)-s_cb(i)))/(18d0*((s_cb(i)-s_cb(i+2))  + &
                            (s_cb(i+1)-s_cb(i+2)))*(s_cb(i-1)-s_cb(i+3))     * &
                            (s_cb(i+3)-s_cb(i-2)))
                        d_L(2,i+1) = &
                            (SQRT(3d0)*(s_cb(i+1)-s_cb(i))**3d0              + &
                            36d0*(s_cb(i+1)-s_cb(i-1))*(s_cb(i+2)-s_cb(i+1)) * &
                            (s_cb(i+3)-s_cb(i+1))-6d0*(s_cb(i+1)             - &
                            s_cb(i))**2d0*((s_cb(i-1)-s_cb(i+1))+(s_cb(i+2)  - &
                            s_cb(i+1))+(s_cb(i+3)-s_cb(i+1)))                + &
                            18d0*(s_cb(i+1)-s_cb(i))*((s_cb(i+2)-s_cb(i+1))  * &
                            (s_cb(i+1)-s_cb(i+3))+(s_cb(i+1)-s_cb(i-1))      * &
                            ((s_cb(i+2)-s_cb(i+1))+(s_cb(i+3)-s_cb(i+1)))))  / &
                            (18d0*((s_cb(i)-s_cb(i-1))+(s_cb(i+1)            - &
                            s_cb(i-1)))*(s_cb(i-2)-s_cb(i+2))*(s_cb(i-2)     - &
                            s_cb(i+3)))
                        
                        d_R(1,i+1) = 1d0 - d_R(0,i+1) - d_R(2,i+1)
                        d_L(1,i+1) = 1d0 - d_L(0,i+1) - d_L(2,i+1)
                        
                    END DO
                    
                END IF
                
                ! Modifying the ideal weights coefficients in the neighborhood
                ! of beginning and end Riemann state extrapolation BC to avoid
                ! any contributions from outside of the physical domain during
                ! the WENO reconstruction
                IF (null_weights) THEN
                    IF(bc_s%beg == -4) THEN
                        d_R(1:2, 0 ) = 0d0; d_R(0, 0 ) = 1d0
                        d_L(1:2, 0 ) = 0d0; d_L(0, 0 ) = 1d0
                        d_R( 2 , 1 ) = 0d0; d_R(:, 1 ) = d_R(:, 1 )/SUM(d_R(:, 1 ))
                        d_L( 2 , 1 ) = 0d0; d_L(:, 1 ) = d_L(:, 1 )/SUM(d_L(:, 1 ))
                    END IF
                    
                    IF(bc_s%end == -4) THEN
                        d_R( 0 ,s-1) = 0d0; d_R(:,s-1) = d_R(:,s-1)/SUM(d_R(:,s-1))
                        d_L( 0 ,s-1) = 0d0; d_L(:,s-1) = d_L(:,s-1)/SUM(d_L(:,s-1))
                        d_R(0:1, s ) = 0d0; d_R(2, s ) = 1d0
                        d_L(0:1, s ) = 0d0; d_L(2, s ) = 1d0
                    END IF
                END IF
                
            END IF
            ! END: Computing WENO5 Coefficients ================================
            
            
            ! Nullifying WENO coefficients and cell-boundary locations pointers
            NULLIFY(poly_coef_L, poly_coef_R, d_L, d_R, beta_coef, s_cb)
            
            
        END SUBROUTINE s_compute_weno_coefficients ! ---------------------------
        
        
        
        
        !>  The purpose of the procedure is to associate the WENO
        !!      coefficients' pointers with their appropriate targets,
        !!      based on the coordinate direction and the location of
        !!      the WENO reconstruction.      
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param weno_loc Location of the WENO reconstruction
        SUBROUTINE s_associate_weno_coefficients_pointers(weno_dir, weno_loc)

            
            INTEGER, INTENT(IN) :: weno_dir, weno_loc
            
            
            ! Associating WENO Coefficients in x-direction =====================
            IF(weno_dir == 1) THEN
                
                IF(weno_loc == 1) THEN                  ! Cell-boundaries
                    poly_coef_L => poly_coef_cbL_x
                    poly_coef_R => poly_coef_cbR_x
                    d_L => d_cbL_x
                    d_R => d_cbR_x
                ELSE                                    ! Quadrature points
                    poly_coef_L => poly_coef_qpL_x
                    poly_coef_R => poly_coef_qpR_x
                    d_L => d_qpL_x
                    d_R => d_qpR_x
                END IF
                
                beta_coef => beta_coef_x
                
            ! ==================================================================
            
            
            ! Associating WENO Coefficients in y-direction =====================
            ELSEIF(weno_dir == 2) THEN
                
                IF(weno_loc == 1) THEN                  ! Cell-boundaries
                    poly_coef_L => poly_coef_cbL_y
                    poly_coef_R => poly_coef_cbR_y
                    d_L => d_cbL_y
                    d_R => d_cbR_y
                ELSE                                    ! Quadrature points
                    poly_coef_L => poly_coef_qpL_y
                    poly_coef_R => poly_coef_qpR_y
                    d_L => d_qpL_y
                    d_R => d_qpR_y
                END IF
                
                beta_coef => beta_coef_y
                
            ! ==================================================================
            
            
            ! Associating WENO Coefficients in z-direction =====================
            ELSE
                
                IF(weno_loc == 1) THEN                  ! Cell-boundaries
                    poly_coef_L => poly_coef_cbL_z
                    poly_coef_R => poly_coef_cbR_z
                    d_L => d_cbL_z
                    d_R => d_cbR_z
                ELSE                                    ! Quadrature points
                    poly_coef_L => poly_coef_qpL_z
                    poly_coef_R => poly_coef_qpR_z
                    d_L => d_qpL_z
                    d_R => d_qpR_z
                END IF
                
                beta_coef => beta_coef_z
                
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_associate_weno_coefficients_pointers ! ----------------
        
        
        
        !>  WENO reconstruction that is improved with monotonicity
        !!      preserving bounds (MPWENO) and a mapping function that
        !!      boosts the accuracy of the non-linear weights (WENOM).
        !!      MPWENO, Balsara and Shu (2000), prevents reconstructed
        !!      values to reside outside the range set by the stencil,
        !!      while WENOM, Henrick et al. (2005), recovers the order
        !!      of accuracy of the reconstruction for critical points.
        !!      Notice that the basic WENO scheme is implemented based
        !!      on the work of Jiang and Shu (1996).
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param cd_vars Characteristic decomposition state variables type
        !! @param norm_dir Characteristic decommposition coordinate direction
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param weno_loc Coordinate location of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
        SUBROUTINE s_weno( v_vf, vL_vf, vR_vf, cd_vars,  & ! -------------------
                           norm_dir, weno_dir, weno_loc, &
                           ix,iy,iz                      )

            
            TYPE(scalar_field), DIMENSION(:), INTENT(IN) :: v_vf
            TYPE(scalar_field), DIMENSION(:), INTENT(INOUT) :: vL_vf, vR_vf
            INTEGER, INTENT(IN) :: cd_vars
            INTEGER, INTENT(IN) :: norm_dir
            INTEGER, INTENT(IN) :: weno_dir, weno_loc
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            REAL(KIND(0d0)), DIMENSION(-weno_polyn:weno_polyn-1) :: dvd !<
            !! Newton divided differences
            
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn) ::  poly_L,  poly_R  !< Left/right polynominals
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn) :: alpha_L, alpha_R  !< Left/right nonlinear weights
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn) :: omega_L, omega_R  !< Left/right nonlinear weights
            

            REAL(KIND(0d0)), DIMENSION(0:weno_polyn) :: beta !< Smoothness indicators


            REAL(KIND(0d0)), DIMENSION(-weno_polyn:weno_polyn) :: scaling_stencil, scaled_vars
            REAL(KIND(0d0)) :: min_u, max_u 
            

            INTEGER :: i,j,k,l,q !< Generic loop iterators
            
            
            ! Reshaping and/or projecting onto characteristic fields inputted
            ! data and in addition associating the WENO coefficients pointers
            IF(weno_order /= 1) THEN
                CALL s_initialize_weno( v_vf, vL_vf, vR_vf, cd_vars, &
                                        norm_dir, weno_dir, ix,iy,iz )
                CALL s_associate_weno_coefficients_pointers(weno_dir, weno_loc)
            END IF
            
            
            ! WENO1 ============================================================
            IF(weno_order == 1) THEN
                
                DO i = 1, UBOUND(v_vf,1)
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                                vL_vf(i)%sf(j,k,l) = v_vf(i)%sf(j,k,l)
                                vR_vf(i)%sf(j,k,l) = v_vf(i)%sf(j,k,l)
                            END DO
                        END DO
                    END DO
                END DO
                
            ! ==================================================================
            
            
            ! WENO3 ============================================================
            ELSEIF(weno_order == 3) THEN
                
                DO i = 1, v_size
                    DO l = is3%beg, is3%end
                        DO k = is2%beg, is2%end
                            DO j = is1%beg, is1%end
                                ! reconstruct from left side

                                dvd( 0) = v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL( 0)%vf(i)%sf(j,k,l)
                                dvd(-1) = v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL(-1)%vf(i)%sf(j,k,l)
                                
                                ! poly_coef_R(0,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                !                        (s_cb( i )-s_cb(i+2)) = 1/2
                                ! poly_coef_R(1,0,i+1) = (s_cb( i )-s_cb(i+1)) / &
                                !                        (s_cb(i-1)-s_cb(i+1)) = 1/2
                                
                                ! poly_coef_L(0,0,i+1) = -poly_coef_R(0,0,i+1) = -1/2
                                ! poly_coef_L(1,0,i+1) = -poly_coef_R(1,0,i+1) = -1/2
 
                                ! so: poly_L[0] = v[j]-(1/2)(v[j+1]-v[j])
                                ! so: poly_L[1] = v[j]-(1/2)(v[j]-v[j-1])
                                poly_L(0) = v_rs_wsL(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_L(0,0,j)*dvd( 0)
                                poly_L(1) = v_rs_wsL(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_L(1,0,j)*dvd(-1)
                                
                                ! beta_coef(0,0,i+1) = 4d0*(s_cb( i )-s_cb(i+1))**2d0 / &
                                !                          (s_cb( i )-s_cb(i+2))**2d0
                                ! beta_coef(1,0,i+1) = 4d0*(s_cb( i )-s_cb(i+1))**2d0 / &
                                !                          (s_cb(i-1)-s_cb(i+1))**2d0

                                ! so for uniform mesh: beta_coef(0,0,j) = 4*((-dx)**2) / (-2dx)**2 = 1
                                ! so for uniform mesh: beta_coef(1,0,j) = 4*((-dx)**2) / (-2dx)**2 = 1

                                ! so: beta[0] = (v[j+1]-v[j])**2 + weno_eps
                                ! so: beta[1] = (v[j]-v[j-1])**2 + weno_eps
                                beta(0) = beta_coef(0,0,j)*dvd( 0)*dvd( 0) &
                                        + weno_eps
                                beta(1) = beta_coef(1,0,j)*dvd(-1)*dvd(-1) &
                                        + weno_eps
                                
                                ! d_L(0,i+1) = (s_cb(i-1)-s_cb( i )) / &
                                !              (s_cb(i-1)-s_cb(i+2))
                                !            = (-dx)/(-3dx) = 1/3
                                ! d_L(1,i+1) = 1d0 - d_L(0,i+1)
                                !            = 1-1/3 = 2/3

                                ! so: alpha_L[0] = (1/3)/( (v[j+1]-v[j])^2 + weno_eps )^2
                                ! so: alpha_L[1] = (2/3)/( (v[j]-v[j-1])^2 + weno_eps )^2
                                alpha_L = d_L(:,j)/(beta*beta)
                                
                                ! so: omega_L[0] = (1/3)/( (v[j+1]-v[j])^2 + weno_eps )^2
                                !                  --------------------------------------
                                ! (2/3)/( (v[j]-v[j-1])^2 + weno_eps )^2 + (1/3)/( (v[j+1]-v[j])^2 + weno_eps )^2

                                ! so: omega_L[1] = (2/3)/( (v[j]-v[j-1])^2 + weno_eps )^2
                                !                  --------------------------------------
                                ! (2/3)/( (v[j]-v[j-1])^2 + weno_eps )^2 + (1/3)/( (v[j+1]-v[j])^2 + weno_eps )^2
                                omega_L = alpha_L/SUM(alpha_L)


                                ! reconstruct from right side
                                dvd( 0) = v_rs_wsR( 1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR( 0)%vf(i)%sf(j,k,l)
                                dvd(-1) = v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR(-1)%vf(i)%sf(j,k,l)

                                ! poly_coef_R(0,0,j) = -dx/(-2 dx) = 1/2
                                ! poly_coef_R(1,0,j) = -dx/(-2 dx) = 1/2

                                poly_R(0) = v_rs_wsR(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_R(0,0,j)*dvd( 0)
                                poly_R(1) = v_rs_wsR(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_R(1,0,j)*dvd(-1)
                                
                                beta(0) = beta_coef(0,0,j)*dvd( 0)*dvd( 0) &
                                        + weno_eps
                                beta(1) = beta_coef(1,0,j)*dvd(-1)*dvd(-1) &
                                        + weno_eps
                                

                                ! d_R(0,i+1) = (s_cb(i-1)-s_cb(i+1)) / &
                                !              (s_cb(i-1)-s_cb(i+2))
                                !            = 2/3
                                ! d_R(1,i+1) = 1d0 - d_R(0,i+1)
                                !            = 1-2/3 = 1/3

                                alpha_R = d_R(:,j)/(beta*beta)
                                
                                omega_R = alpha_R/SUM(alpha_R)
                                
                                IF(mapped_weno) THEN
                                    CALL s_map_nonlinear_weights( d_L(:,j), &
                                                                   alpha_L, &
                                                                   omega_L  )
                                    CALL s_map_nonlinear_weights( d_R(:,j), &
                                                                   alpha_R, &
                                                                   omega_R  )
                                END IF
                                
                                vL_rs_vf(i)%sf(j,k,l) = SUM(omega_L*poly_L)
                                vR_rs_vf(i)%sf(j,k,l) = SUM(omega_R*poly_R)
                                
                            END DO
                        END DO
                    END DO
                END DO
                
            ! END: WENO3 =======================================================
            
            
            ! WENO5 ============================================================
            ELSE
                
                DO i = 1, v_size
                    DO l = is3%beg, is3%end
                        DO k = is2%beg, is2%end
                            DO j = is1%beg, is1%end
                                
                                dvd( 1) = v_rs_wsL( 2)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL( 1)%vf(i)%sf(j,k,l)
                                dvd( 0) = v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL( 0)%vf(i)%sf(j,k,l)
                                dvd(-1) = v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL(-1)%vf(i)%sf(j,k,l)
                                dvd(-2) = v_rs_wsL(-1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsL(-2)%vf(i)%sf(j,k,l)
                                
                                poly_L(0) = v_rs_wsL(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_L(0,0,j)*dvd( 1)  &
                                          + poly_coef_L(0,1,j)*dvd( 0)
                                poly_L(1) = v_rs_wsL(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_L(1,0,j)*dvd( 0)  &
                                          + poly_coef_L(1,1,j)*dvd(-1)
                                poly_L(2) = v_rs_wsL(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_L(2,0,j)*dvd(-1)  &
                                          + poly_coef_L(2,1,j)*dvd(-2)
                                
                                beta(0) = beta_coef(0,0,j)*dvd( 1)*dvd( 1) &
                                        + beta_coef(0,1,j)*dvd( 1)*dvd( 0) &
                                        + beta_coef(0,2,j)*dvd( 0)*dvd( 0) &
                                        + weno_eps
                                beta(1) = beta_coef(1,0,j)*dvd( 0)*dvd( 0) &
                                        + beta_coef(1,1,j)*dvd( 0)*dvd(-1) &
                                        + beta_coef(1,2,j)*dvd(-1)*dvd(-1) &
                                        + weno_eps
                                beta(2) = beta_coef(2,0,j)*dvd(-1)*dvd(-1) &
                                        + beta_coef(2,1,j)*dvd(-1)*dvd(-2) &
                                        + beta_coef(2,2,j)*dvd(-2)*dvd(-2) &
                                        + weno_eps
                                
                                alpha_L = d_L(:,j)/(beta*beta)
                                
                                omega_L = alpha_L/SUM(alpha_L)
                                
                                dvd( 1) = v_rs_wsR( 2)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR( 1)%vf(i)%sf(j,k,l)
                                dvd( 0) = v_rs_wsR( 1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR( 0)%vf(i)%sf(j,k,l)
                                dvd(-1) = v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR(-1)%vf(i)%sf(j,k,l)
                                dvd(-2) = v_rs_wsR(-1)%vf(i)%sf(j,k,l) &
                                        - v_rs_wsR(-2)%vf(i)%sf(j,k,l)
                                
                                poly_R(0) = v_rs_wsR(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_R(0,0,j)*dvd( 1)  &
                                          + poly_coef_R(0,1,j)*dvd( 0)
                                poly_R(1) = v_rs_wsR(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_R(1,0,j)*dvd( 0)  &
                                          + poly_coef_R(1,1,j)*dvd(-1)
                                poly_R(2) = v_rs_wsR(0)%vf(i)%sf(j,k,l) &
                                          + poly_coef_R(2,0,j)*dvd(-1)  &
                                          + poly_coef_R(2,1,j)*dvd(-2)
                                
                                beta(0) = beta_coef(0,0,j)*dvd( 1)*dvd( 1) &
                                        + beta_coef(0,1,j)*dvd( 1)*dvd( 0) &
                                        + beta_coef(0,2,j)*dvd( 0)*dvd( 0) &
                                        + weno_eps
                                beta(1) = beta_coef(1,0,j)*dvd( 0)*dvd( 0) &
                                        + beta_coef(1,1,j)*dvd( 0)*dvd(-1) &
                                        + beta_coef(1,2,j)*dvd(-1)*dvd(-1) &
                                        + weno_eps
                                beta(2) = beta_coef(2,0,j)*dvd(-1)*dvd(-1) &
                                        + beta_coef(2,1,j)*dvd(-1)*dvd(-2) &
                                        + beta_coef(2,2,j)*dvd(-2)*dvd(-2) &
                                        + weno_eps
                                
                                alpha_R = d_R(:,j)/(beta*beta)
                                
                                omega_R = alpha_R/SUM(alpha_R)
                                
                                IF(mapped_weno) THEN
                                    CALL s_map_nonlinear_weights( d_L(:,j), &
                                                                   alpha_L, &
                                                                   omega_L  )
                                    CALL s_map_nonlinear_weights( d_R(:,j), &
                                                                   alpha_R, &
                                                                   omega_R  )
                                END IF
                                
                                vL_rs_vf(i)%sf(j,k,l) = SUM(omega_L*poly_L)
                                vR_rs_vf(i)%sf(j,k,l) = SUM(omega_R*poly_R)
                                
                                IF(mp_weno .AND. weno_loc == 1) THEN
                                    CALL s_preserve_monotonicity(i,j,k,l)
                                END IF
                                
                            END DO
                        END DO
                    END DO
                END DO
                
            END IF
            ! END: WENO5 =======================================================
            
            
            ! Reshaping and/or projecting onto physical fields the outputted
            ! data, as well as disassociating the WENO coefficients pointers
            IF(weno_order /= 1) THEN
                CALL s_finalize_weno(vL_vf, vR_vf, cd_vars, weno_dir, ix,iy,iz)
                NULLIFY(poly_coef_L, poly_coef_R, d_L, d_R, beta_coef)
            END IF
            
            
        END SUBROUTINE s_weno ! ------------------------------------------------
        
        
        
        
        !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      WENO reconstruction.      
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param cd_vars Characteristic decomposition state variables type
        !! @param norm_dir Characteristic decommposition coordinate direction
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
        SUBROUTINE s_initialize_weno( v_vf, vL_vf, vR_vf, cd_vars, & ! ---------
                                      norm_dir, weno_dir, ix,iy,iz )

            
            TYPE(scalar_field), DIMENSION(:), INTENT(IN) :: v_vf
            TYPE(scalar_field), DIMENSION(:), INTENT(INOUT) :: vL_vf, vR_vf
            INTEGER, INTENT(IN) :: cd_vars
            INTEGER, INTENT(IN) :: norm_dir
            INTEGER, INTENT(IN) :: weno_dir
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            
            ! Determining the number of cell-average variables which will be
            ! WENO-reconstructed and mapping their indical bounds in the x-,
            ! y- and z-directions to those in the s1-, s2- and s3-directions
            ! as to reshape the inputted data in the coordinate direction of
            ! the WENO reconstruction
            v_size = UBOUND(v_vf,1)
            
            IF(weno_dir == 1) THEN
                is1 = ix; is2 = iy; is3 = iz
            ELSEIF(weno_dir == 2) THEN
                is1 = iy; is2 = ix; is3 = iz
            ELSE
                is1 = iz; is2 = iy; is3 = ix
            END IF
            
            
            ! Allocating the cell-average variables, which are reshaped, and/or
            ! characteristically decomposed, in the coordinate direction of the
            ! WENO reconstruction
            DO i = -weno_polyn, weno_polyn
                
                ALLOCATE(v_rs_wsL(i)%vf(1:v_size), v_rs_wsR(i)%vf(1:v_size))
                
                DO j = 1, v_size
                    
                    ALLOCATE(v_rs_wsL(i)%vf(j)%sf( is1%beg:is1%end, &
                                                   is2%beg:is2%end, &
                                                   is3%beg:is3%end ))
                    
                    IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                        ALLOCATE(v_rs_wsR(i)%vf(j)%sf( is1%beg:is1%end, &
                                                       is2%beg:is2%end, &
                                                       is3%beg:is3%end ))
                    ELSE
                        v_rs_wsR(i)%vf(j)%sf => v_rs_wsL(i)%vf(j)%sf
                    END IF
                    
                END DO
                
            END DO
            
            
            ! Allocating the left and right WENO reconstructions of the cell-
            ! average variables that are reshaped, and/or characteristically
            ! decomposed, in the coordinate direction of WENO reconstruction
            ALLOCATE(vL_rs_vf(1:v_size), vR_rs_vf(1:v_size))
            
            IF(weno_dir == 1) THEN
                DO i = 1, v_size
                    vL_rs_vf(i)%sf => vL_vf(i)%sf
                    vR_rs_vf(i)%sf => vR_vf(i)%sf
                END DO
            ELSE
                DO i = 1, v_size
                    ALLOCATE(vL_rs_vf(i)%sf( is1%beg:is1%end, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end ))
                    ALLOCATE(vR_rs_vf(i)%sf( is1%beg:is1%end, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end ))
                END DO
            END IF
            
            
            ! Reshaping/Projecting onto Characteristic Fields in x-direction ===
            IF(weno_dir == 1) THEN
                
                DO i = -weno_polyn, weno_polyn
                    DO j = 1, v_size
                        DO k = ix%beg, ix%end
                            v_rs_wsL(i)%vf(j)%sf(k,:,:) = &
                                 v_vf(j)%sf(i+k,iy%beg:iy%end,iz%beg:iz%end)
                        END DO
                    END DO
                END DO
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    DO i = -weno_polyn, weno_polyn
                        DO j = 1, v_size
                            v_rs_wsR(i)%vf(j)%sf = v_rs_wsL(i)%vf(j)%sf
                        END DO
                    END DO
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsL, norm_dir, ix,iy,iz,-1 )
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsR, norm_dir, ix,iy,iz, 0 )
                    ELSE
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsL, norm_dir, ix,iy,iz,-1 )
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsR, norm_dir, ix,iy,iz, 0 )
                    END IF
                    
                END IF
                
            ! ==================================================================
            
            
            ! Reshaping/Projecting onto Characteristic Fields in y-direction ===
            ELSEIF(weno_dir == 2) THEN
                
                DO i = -weno_polyn, weno_polyn
                    DO j = 1, v_size
                        DO k = ix%beg, ix%end
                            DO l = iy%beg, iy%end
                                v_rs_wsL(i)%vf(j)%sf(l,k,:) = &
                                         v_vf(j)%sf(k,i+l,iz%beg:iz%end)
                            END DO
                        END DO
                    END DO
                END DO
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    DO i = -weno_polyn, weno_polyn
                        DO j = 1, v_size
                            v_rs_wsR(i)%vf(j)%sf = v_rs_wsL(i)%vf(j)%sf
                        END DO
                    END DO
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsL, norm_dir, iy,ix,iz,-1 )
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsR, norm_dir, iy,ix,iz, 0 )
                    ELSE
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsL, norm_dir, iy,ix,iz,-1 )
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsR, norm_dir, iy,ix,iz, 0 )
                    END IF
                    
                END IF
                
            ! ==================================================================
            
            
            ! Reshaping/Projecting onto Characteristic Fields in z-direction ===
            ELSE
                
                DO i = -weno_polyn, weno_polyn
                    DO j = 1, v_size
                        DO k = ix%beg,ix%end
                            DO l = iz%beg,iz%end
                                v_rs_wsL(i)%vf(j)%sf(l,:,k) = &
                                         v_vf(j)%sf(k,iy%beg:iy%end,i+l)
                            END DO
                        END DO
                    END DO
                END DO
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    DO i = -weno_polyn, weno_polyn
                        DO j = 1, v_size
                            v_rs_wsR(i)%vf(j)%sf = v_rs_wsL(i)%vf(j)%sf
                        END DO
                    END DO
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsL, norm_dir, iz,iy,ix,-1 )
                      CALL s_convert_conservative_to_characteristic_variables( &
                                               v_rs_wsR, norm_dir, iz,iy,ix, 0 )
                    ELSE
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsL, norm_dir, iz,iy,ix,-1 )
                      CALL s_convert_primitive_to_characteristic_variables(    &
                                               v_rs_wsR, norm_dir, iz,iy,ix, 0 )
                    END IF
                    
                END IF
                
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_initialize_weno ! -------------------------------------
        
        
        
        
        !>  The goal of this procedure is to map the nonlinear WENO
        !!      weights to the more accurate nonlinear WENOM weights in
        !!      order to reinstate the optimal order of accuracy of the
        !!      reconstruction in the proximity of critical points, see
        !!      Henrick et al. (2005).
        !!  @param d_K Cell boundary pointer
        !!  @param alpha_K ideal weights
        !!  @param omega_K nonlinear weights
        SUBROUTINE s_map_nonlinear_weights(d_K, alpha_K, omega_K) ! ------------
            
            
            ! Ideal and nonlinear weights
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn), INTENT(IN)    ::     d_K
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn), INTENT(INOUT) :: alpha_K
            REAL(KIND(0d0)), DIMENSION(0:weno_polyn), INTENT(INOUT) :: omega_K
            
            
            ! Mapping the WENO nonlinear weights to the WENOM nonlinear weights
            IF(MINVAL(d_K) == 0d0 .OR. MAXVAL(d_K) == 1d0) RETURN
            
            alpha_K = (d_K*(1d0 + d_K - 3d0*omega_K) + omega_K**2d0) &
                    * (omega_K/(d_K**2d0 + omega_K*(1d0 - 2d0*d_K)))
            
            omega_K = alpha_K/SUM(alpha_K)
            
            
        END SUBROUTINE s_map_nonlinear_weights ! -------------------------------
        
        
        
        !>  The goal of this subroutine is to ensure that the WENO
        !!      reconstruction is monotonic. The latter is achieved by
        !!      enforcing monotonicity preserving bounds of Suresh and
        !!      Huynh (1997). The resulting MPWENO reconstruction, see
        !!      Balsara and Shu (2000), ensures that the reconstructed
        !!      values do not reside outside the range spanned by WENO
        !!      stencil.
        !!  @param i Equation number
        !!  @param j First-coordinate cell index
        !!  @param k Second-coordinate cell index
        !!  @param l Third-coordinate cell index
        SUBROUTINE s_preserve_monotonicity(i,j,k,l) ! --------------------------

            
            INTEGER, INTENT(IN) :: i,j,k,l
            

            REAL(KIND(0d0)), DIMENSION(-1:1) :: d !< Curvature measures at the zone centers
            

            REAL(KIND(0d0)) :: d_MD, d_LC !<
            !! Median (md) curvature and large curvature (LC) measures
            
            ! The left and right upper bounds (UL), medians, large curvatures,
            ! minima, and maxima of the WENO-reconstructed values of the cell-
            ! average variables.
            REAL(KIND(0d0)) :: vL_UL , vR_UL
            REAL(KIND(0d0)) :: vL_MD , vR_MD
            REAL(KIND(0d0)) :: vL_LC , vR_LC
            REAL(KIND(0d0)) :: vL_min, vR_min
            REAL(KIND(0d0)) :: vL_max, vR_max
            
            REAL(KIND(0d0)), PARAMETER :: alpha = 2d0 !>
            !! Determines the maximum Courant–Friedrichs–Lewy (CFL) number that
            !! may be utilized with the scheme. In theory, for stability, a CFL
            !! number less than 1/(1+alpha) is necessary. The default value for
            !! alpha is 2.
            


            REAL(KIND(0d0)), PARAMETER :: beta = 4d0/3d0 !<
            !! Determines the amount of freedom available from utilizing a large
            !! value for the local curvature. The default value for beta is 4/3.
            
            
            ! Left Monotonicity Preserving Bound ===============================
            d(-1) = v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                  + v_rs_wsL(-2)%vf(i)%sf(j,k,l) &
                  - v_rs_wsL(-1)%vf(i)%sf(j,k,l) &
                  * 2d0
            d( 0) = v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                  + v_rs_wsL(-1)%vf(i)%sf(j,k,l) &
                  - v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                  * 2d0
            d( 1) = v_rs_wsL( 2)%vf(i)%sf(j,k,l) &
                  + v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                  - v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                  * 2d0
            
            d_MD = (SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,4d0*d(0)-d(-1))) &
                 *  ABS((SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,d(-1)))     &
                       *(SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,d( 0))))    &
                 *  MIN(ABS(4d0*d(-1)-d( 0)),ABS(d(-1)),                 &
                        ABS(4d0*d( 0)-d(-1)),ABS(d( 0)))/8d0
            
            d_LC = (SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,4d0*d(1)-d( 0))) &
                 *  ABS((SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,d( 0)))     &
                       *(SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,d( 1))))    &
                 *  MIN(ABS(4d0*d( 0)-d( 1)),ABS(d( 0)),                 &
                        ABS(4d0*d( 1)-d( 0)),ABS(d( 1)))/8d0
            
            vL_UL =   v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                  - ( v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                    - v_rs_wsL( 0)%vf(i)%sf(j,k,l) )*alpha
            
            vL_MD = ( v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                    + v_rs_wsL(-1)%vf(i)%sf(j,k,l) &
                    - d_MD                         )*5d-1
            
            vL_LC =   v_rs_wsL( 0)%vf(i)%sf(j,k,l) &
                  - ( v_rs_wsL( 1)%vf(i)%sf(j,k,l) &
                    - v_rs_wsL( 0)%vf(i)%sf(j,k,l) )*5d-1 + beta*d_LC
            
            vL_min = MAX(MIN(v_rs_wsL( 0)%vf(i)%sf(j,k,l), &
                             v_rs_wsL(-1)%vf(i)%sf(j,k,l), &
                             vL_MD                      ), &
                         MIN(v_rs_wsL( 0)%vf(i)%sf(j,k,l), &
                             vL_UL                       , &
                             vL_LC                       ) )
            
            vL_max = MIN(MAX(v_rs_wsL( 0)%vf(i)%sf(j,k,l), &
                             v_rs_wsL(-1)%vf(i)%sf(j,k,l), &
                             vL_MD                      ), &
                         MAX(v_rs_wsL( 0)%vf(i)%sf(j,k,l), &
                             vL_UL                       , &
                             vL_LC                       ) )
            
            vL_rs_vf(i)%sf(j,k,l) = vL_rs_vf(i)%sf(j,k,l) &
                                  + (SIGN(5d-1,vL_min-vL_rs_vf(i)%sf(j,k,l))  &
                                    +SIGN(5d-1,vL_max-vL_rs_vf(i)%sf(j,k,l))) &
                                    * MIN( ABS(vL_min-vL_rs_vf(i)%sf(j,k,l)), &
                                           ABS(vL_max-vL_rs_vf(i)%sf(j,k,l)))
            ! END: Left Monotonicity Preserving Bound ==========================
            
            
            ! Right Monotonicity Preserving Bound ==============================
            d(-1) = v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                  + v_rs_wsR(-2)%vf(i)%sf(j,k,l) &
                  - v_rs_wsR(-1)%vf(i)%sf(j,k,l)*2d0
            d( 0) = v_rs_wsR( 1)%vf(i)%sf(j,k,l) &
                  + v_rs_wsR(-1)%vf(i)%sf(j,k,l) &
                  - v_rs_wsR( 0)%vf(i)%sf(j,k,l)*2d0
            d( 1) = v_rs_wsR( 2)%vf(i)%sf(j,k,l) &
                  + v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                  - v_rs_wsR( 1)%vf(i)%sf(j,k,l)*2d0
            
            d_MD = (SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,4d0*d(1)-d( 0))) &
                 *  ABS((SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,d( 0)))     &
                       *(SIGN(1d0,4d0*d( 0)-d(1)) + SIGN(1d0,d( 1))))    &
                 *  MIN(ABS(4d0*d( 0)-d( 1)),ABS(d( 0)),                 &
                        ABS(4d0*d( 1)-d( 0)),ABS(d( 1)))/8d0
            
            d_LC = (SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,4d0*d(0)-d(-1))) &
                 *  ABS((SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,d(-1)))     &
                       *(SIGN(1d0,4d0*d(-1)-d(0)) + SIGN(1d0,d( 0))))    &
                 *  MIN(ABS(4d0*d(-1)-d( 0)),ABS(d(-1)),                 &
                        ABS(4d0*d( 0)-d(-1)),ABS(d( 0)))/8d0
            
            vR_UL =   v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                  + ( v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                    - v_rs_wsR(-1)%vf(i)%sf(j,k,l) )*alpha
            
            vR_MD = ( v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                    + v_rs_wsR( 1)%vf(i)%sf(j,k,l) &
                    - d_MD                         )*5d-1
            
            vR_LC =   v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                  + ( v_rs_wsR( 0)%vf(i)%sf(j,k,l) &
                    - v_rs_wsR(-1)%vf(i)%sf(j,k,l) )*5d-1 + beta*d_LC
            
            vR_min = MAX(MIN(v_rs_wsR(0)%vf(i)%sf(j,k,l), &
                             v_rs_wsR(1)%vf(i)%sf(j,k,l), &
                             vR_MD                     ), &
                         MIN(v_rs_wsR(0)%vf(i)%sf(j,k,l), &
                             vR_UL                      , &
                             vR_LC                      ) )
            
            vR_max = MIN(MAX(v_rs_wsR(0)%vf(i)%sf(j,k,l), &
                             v_rs_wsR(1)%vf(i)%sf(j,k,l), &
                             vR_MD                     ), &
                         MAX(v_rs_wsR(0)%vf(i)%sf(j,k,l), &
                             vR_UL                      , &
                             vR_LC                      ) )
            
            vR_rs_vf(i)%sf(j,k,l) = vR_rs_vf(i)%sf(j,k,l) &
                                  + (SIGN(5d-1,vR_min-vR_rs_vf(i)%sf(j,k,l))  &
                                    +SIGN(5d-1,vR_max-vR_rs_vf(i)%sf(j,k,l))) &
                                    * MIN( ABS(vR_min-vR_rs_vf(i)%sf(j,k,l)), &
                                           ABS(vR_max-vR_rs_vf(i)%sf(j,k,l)))
            ! END: Right Monotonicity Preserving Bound =========================
            
            
        END SUBROUTINE s_preserve_monotonicity ! -------------------------------
        
        
        
        
        !>  Deallocation and/or disassociation procedures that are
        !!      necessary in order to finalize the WENO reconstruction
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param cd_vars Characteristic decomposition state variables type
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
        SUBROUTINE s_finalize_weno( vL_vf, vR_vf, cd_vars, & ! -----------------
                                    weno_dir, ix,iy,iz     )
            
            TYPE(scalar_field), DIMENSION(:), INTENT(INOUT) :: vL_vf, vR_vf
            INTEGER, INTENT(IN) :: cd_vars
            INTEGER, INTENT(IN) :: weno_dir
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: i,j,k !< Generic loop iterators
            
            ! Reshaping/Projecting onto Physical Fields in x-direction =========
            IF(weno_dir == 1) THEN
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vL_rs_vf, ix,iy,iz,-1 )
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vR_rs_vf, ix,iy,iz, 0 )
                    ELSE
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vL_rs_vf, ix,iy,iz,-1 )
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vR_rs_vf, ix,iy,iz, 0 )
                    END IF
                    
                END IF
                
            ! ==================================================================
            
            
            ! Reshaping/Projecting onto Physical Fields in y-direction =========
            ELSEIF(weno_dir == 2) THEN
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vL_rs_vf, iy,ix,iz,-1 )
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vR_rs_vf, iy,ix,iz, 0 )
                    ELSE
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vL_rs_vf, iy,ix,iz,-1 )
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vR_rs_vf, iy,ix,iz, 0 )
                    END IF
                    
                END IF
                
                DO i = 1, v_size
                    DO j = ix%beg, ix%end
                        DO k = iy%beg, iy%end
                          vL_vf(i)%sf(j,k,iz%beg:iz%end) = vL_rs_vf(i)%sf(k,j,:)
                          vR_vf(i)%sf(j,k,iz%beg:iz%end) = vR_rs_vf(i)%sf(k,j,:)
                        END DO
                    END DO
                END DO
                
            ! ==================================================================
            
            
            ! Reshaping/Projecting onto Physical Fields in z-direction =========
            ELSE
                
                IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                    
                    IF(cd_vars == 1) THEN
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vL_rs_vf, iz,iy,ix,-1 )
                      CALL s_convert_characteristic_to_conservative_variables( &
                                                         vR_rs_vf, iz,iy,ix, 0 )
                    ELSE
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vL_rs_vf, iz,iy,ix,-1 )
                      CALL s_convert_characteristic_to_primitive_variables(    &
                                                         vR_rs_vf, iz,iy,ix, 0 )
                    END IF
                    
                END IF
                
                DO i = 1, v_size
                    DO j = ix%beg, ix%end
                        DO k = iz%beg, iz%end
                          vL_vf(i)%sf(j,iy%beg:iy%end,k) = vL_rs_vf(i)%sf(k,:,j)
                          vR_vf(i)%sf(j,iy%beg:iy%end,k) = vR_rs_vf(i)%sf(k,:,j)
                        END DO
                    END DO
                END DO
                
            END IF
            ! ==================================================================
            
            
            ! Deallocating the cell-average variables that were reshaped and/or
            ! characteristically decomposed in the coordinate direction of WENO
            ! reconstruction
            DO i = -weno_polyn, weno_polyn
                
                DO j = 1, v_size
                    
                    DEALLOCATE(v_rs_wsL(i)%vf(j)%sf)
                    
                    IF(char_decomp .AND. cd_vars /= dflt_int) THEN
                        DEALLOCATE(v_rs_wsR(i)%vf(j)%sf)
                    ELSE
                        v_rs_wsR(i)%vf(j)%sf => NULL()
                    END IF
                    
                END DO
                
                DEALLOCATE(v_rs_wsL(i)%vf, v_rs_wsR(i)%vf)
                
            END DO
            
            
            ! Deallocating the left and right WENO reconstructions of the cell-
            ! average variables which were reshaped, and/or characteristically
            ! decomposed, in the coordinate direction of WENO reconstruction
            IF(weno_dir == 1) THEN
                DO i = 1, v_size
                    vL_rs_vf(i)%sf => NULL()
                    vR_rs_vf(i)%sf => NULL()
                END DO
            ELSE
                DO i = 1, v_size
                    DEALLOCATE(vL_rs_vf(i)%sf)
                    DEALLOCATE(vR_rs_vf(i)%sf)
                END DO
            END IF
            
            DEALLOCATE(vL_rs_vf, vR_rs_vf)
            
            
        END SUBROUTINE s_finalize_weno ! ---------------------------------------
        
        
        !>  Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_weno_module() ! ----------------------------------

            
            
            IF(weno_order == 1) RETURN
            
            
            ! Deallocating the WENO-stencil of the WENO-reconstructed variables
            DEALLOCATE(v_rs_wsL, v_rs_wsR)
            
            
            ! Deallocating WENO coefficients in x-direction ====================
            DEALLOCATE(poly_coef_cbL_x, poly_coef_cbR_x)
            DEALLOCATE(d_cbL_x, d_cbR_x)
            
            IF(commute_err .OR. split_err) THEN
                DEALLOCATE(poly_coef_qpL_x, poly_coef_qpR_x)
                DEALLOCATE(d_qpL_x, d_qpR_x)
            END IF
            
            DEALLOCATE(beta_coef_x)
            ! ==================================================================
            
            
            ! Deallocating WENO coefficients in y-direction ====================
            IF(n == 0) RETURN
            
            DEALLOCATE(poly_coef_cbL_y, poly_coef_cbR_y)
            DEALLOCATE(d_cbL_y, d_cbR_y)
            
            IF(commute_err .OR. split_err) THEN
                DEALLOCATE(poly_coef_qpL_y, poly_coef_qpR_y)
                DEALLOCATE(d_qpL_y, d_qpR_y)
            END IF
            
            DEALLOCATE(beta_coef_y)
            ! ==================================================================
            
            
            ! Deallocating WENO coefficients in z-direction ====================
            IF(p == 0) RETURN
            
            DEALLOCATE(poly_coef_cbL_z, poly_coef_cbR_z)
            DEALLOCATE(d_cbL_z, d_cbR_z)
            
            IF(commute_err .OR. split_err) THEN
                DEALLOCATE(poly_coef_qpL_z, poly_coef_qpR_z)
                DEALLOCATE(d_qpL_z, d_qpR_z)
            END IF
            
            DEALLOCATE(beta_coef_z)
            ! ==================================================================
            
            
        END SUBROUTINE s_finalize_weno_module ! --------------------------------
        
        
        
        
        
END MODULE m_weno
