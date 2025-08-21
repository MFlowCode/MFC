!>
!! @file m_weno.f90
!! @brief Contains module m_weno
#:include 'macros.fpp'

!> @brief  Weighted essentially non-oscillatory (WENO) reconstruction scheme
!!              that is supplemented with monotonicity preserving bounds (MPWENO)
!!              and a mapping function that boosts the accuracy of the non-linear
!!              weights (WENOM). MPWENO, see Balsara and Shu (2000), prevents the
!!              reconstructed values to lay outside the range set by the stencil,
!!              while WENOM, see Henrick et al. (2005), recovers the formal order
!!              of accuracy of the reconstruction at critical points. Please note
!!              that the basic WENO approach is implemented according to the work
!!              of Jiang and Shu (1996). WENO-Z, which is less dissipative than
!!              WENO-JS and WENO-M, is implemented according to the work of
!!              Borges, et al. (2008). TENO, which is even less dissipative than
!!              WENO-Z but is less robust, is implemented according to the work
!!              of Fu et al. (2016).
module m_weno

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

#ifdef MFC_OPENACC
    use openacc
#endif

    use m_mpi_proxy

    use m_muscl !< For Interface Compression

    private; public :: s_initialize_weno_module, s_initialize_weno, s_finalize_weno_module, s_weno

    !> @name The cell-average variables that will be WENO-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the WENO procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed WENO-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{
    real(wp), allocatable, dimension(:, :, :, :) :: v_rs_ws_x, v_rs_ws_y, v_rs_ws_z
    !> @}
    $:GPU_DECLARE(create='[v_rs_ws_x,v_rs_ws_y,v_rs_ws_z]')

    ! WENO Coefficients

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at
    !! the left and right quadrature points (QP), in the x-, y- and z-directions.
    !! Note that the first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_x
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_y
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_z
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_x
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_y
    real(wp), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_z
    !> @}
    $:GPU_DECLARE(create='[poly_coef_cbL_x,poly_coef_cbL_y,poly_coef_cbL_z]')
    $:GPU_DECLARE(create='[poly_coef_cbR_x,poly_coef_cbR_y,poly_coef_cbR_z]')

    !> @name The ideal weights at the left and the right cell-boundaries and at the
    !! left and the right quadrature points, in x-, y- and z-directions. Note
    !! that the first dimension of the array identifies the weight, while the
    !! last denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(wp), target, allocatable, dimension(:, :) :: d_cbL_x
    real(wp), target, allocatable, dimension(:, :) :: d_cbL_y
    real(wp), target, allocatable, dimension(:, :) :: d_cbL_z

    real(wp), target, allocatable, dimension(:, :) :: d_cbR_x
    real(wp), target, allocatable, dimension(:, :) :: d_cbR_y
    real(wp), target, allocatable, dimension(:, :) :: d_cbR_z
    !> @}
    $:GPU_DECLARE(create='[d_cbL_x,d_cbL_y,d_cbL_z,d_cbR_x,d_cbR_y,d_cbR_z]')

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note
    !! that the first array dimension identifies the smoothness indicator, the
    !! second identifies the position of its coefficients and the last denotes
    !! the cell-location in the relevant coordinate direction.
    !> @{
    real(wp), target, allocatable, dimension(:, :, :) :: beta_coef_x
    real(wp), target, allocatable, dimension(:, :, :) :: beta_coef_y
    real(wp), target, allocatable, dimension(:, :, :) :: beta_coef_z
    !> @}
    $:GPU_DECLARE(create='[beta_coef_x,beta_coef_y,beta_coef_z]')

    ! END: WENO Coefficients

    integer :: v_size !< Number of WENO-reconstructed cell-average variables
    $:GPU_DECLARE(create='[v_size]')

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1_weno, is2_weno, is3_weno
#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
    $:GPU_DECLARE(create='[is1_weno,is2_weno,is3_weno]')
#endif
    !
    !> @}

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_weno_module

        if (weno_order == 1) return

        ! Allocating/Computing WENO Coefficients in x-direction
        is1_weno%beg = -buff_size; is1_weno%end = m - is1_weno%beg
        if (n == 0) then
            is2_weno%beg = 0
        else
            is2_weno%beg = -buff_size; 
        end if

        is2_weno%end = n - is2_weno%beg

        if (p == 0) then
            is3_weno%beg = 0
        else
            is3_weno%beg = -buff_size
        end if

        is3_weno%end = p - is3_weno%beg

        @:ALLOCATE(poly_coef_cbL_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_x(0:weno_num_stencils, is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_x(0:weno_num_stencils, is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn*(weno_polyn + 1)/2 - 1))
        ! Number of cross terms for dvd = (k-1)(k-1+1)/2, where weno_polyn = k-1
        ! Note: k-1 not k because we are using value differences (dvd) not the values themselves

        call s_compute_weno_coefficients(1, is1_weno)

        @:ALLOCATE(v_rs_ws_x(is1_weno%beg:is1_weno%end, &
            is2_weno%beg:is2_weno%end, is3_weno%beg:is3_weno%end, 1:sys_size))

        ! Allocating/Computing WENO Coefficients in y-direction
        if (n == 0) return

        is2_weno%beg = -buff_size; is2_weno%end = n - is2_weno%beg
        is1_weno%beg = -buff_size; is1_weno%end = m - is1_weno%beg

        if (p == 0) then
            is3_weno%beg = 0
        else
            is3_weno%beg = -buff_size
        end if

        is3_weno%end = p - is3_weno%beg

        @:ALLOCATE(poly_coef_cbL_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_y(0:weno_num_stencils, is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_y(0:weno_num_stencils, is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn*(weno_polyn + 1)/2 - 1))

        call s_compute_weno_coefficients(2, is2_weno)

        @:ALLOCATE(v_rs_ws_y(is2_weno%beg:is2_weno%end, &
            is1_weno%beg:is1_weno%end, is3_weno%beg:is3_weno%end, 1:sys_size))

        ! Allocating/Computing WENO Coefficients in z-direction
        if (p == 0) return

        is2_weno%beg = -buff_size; is2_weno%end = n - is2_weno%beg
        is1_weno%beg = -buff_size; is1_weno%end = m - is1_weno%beg
        is3_weno%beg = -buff_size; is3_weno%end = p - is3_weno%beg

        @:ALLOCATE(poly_coef_cbL_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_z(0:weno_num_stencils, is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_z(0:weno_num_stencils, is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, &
            0:weno_polyn*(weno_polyn + 1)/2 - 1))

        call s_compute_weno_coefficients(3, is3_weno)

        @:ALLOCATE(v_rs_ws_z(is3_weno%beg:is3_weno%end, &
            is2_weno%beg:is2_weno%end, is1_weno%beg:is1_weno%end, 1:sys_size))

    end subroutine s_initialize_weno_module

    !>  The purpose of this subroutine is to compute the grid
        !!      dependent coefficients of the WENO polynomials, ideal
        !!      weights and smoothness indicators, provided the order,
        !!      the coordinate direction and the location of the WENO
        !!      reconstruction.
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param is Index bounds in the s-direction
    subroutine s_compute_weno_coefficients(weno_dir, is)

        integer, intent(in) :: weno_dir
        type(int_bounds_info), intent(in) :: is
        integer :: s

        real(wp), pointer, dimension(:) :: s_cb => null() !<
            !! Cell-boundary locations in the s-direction

        type(int_bounds_info) :: bc_s !< Boundary conditions (BC) in the s-direction

        integer :: i !< Generic loop iterator

        real(wp) :: w(1:8) ! Intermediate var for ideal weights: s_cb across overall stencil
        real(wp) :: y(1:4) ! Intermediate var for poly & beta: diff(s_cb) across sub-stencil

        ! Determining the number of cells, the cell-boundary locations and
        ! the boundary conditions in the coordinate direction selected for
        ! the WENO reconstruction
        if (weno_dir == 1) then
            s = m; s_cb => x_cb; bc_s = bc_x
        elseif (weno_dir == 2) then
            s = n; s_cb => y_cb; bc_s = bc_y
        else
            s = p; s_cb => z_cb; bc_s = bc_z
        end if

        #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            ! Computing WENO3 Coefficients
            if (weno_dir == ${WENO_DIR}$) then
                if (weno_order == 3) then
                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = (s_cb(i) - s_cb(i + 1))/ &
                                                              (s_cb(i) - s_cb(i + 2))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = (s_cb(i) - s_cb(i + 1))/ &
                                                              (s_cb(i - 1) - s_cb(i + 1))

                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 0, 0)
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 1, 0)

                        d_cbR_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i + 1))/ &
                                                   (s_cb(i - 1) - s_cb(i + 2))
                        d_cbL_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i))/ &
                                                   (s_cb(i - 1) - s_cb(i + 2))

                        d_cbR_${XYZ}$ (1, i + 1) = 1._wp - d_cbR_${XYZ}$ (0, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1._wp - d_cbL_${XYZ}$ (0, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp/ &
                                                          (s_cb(i) - s_cb(i + 2))**2._wp
                        beta_coef_${XYZ}$ (i + 1, 1, 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp/ &
                                                          (s_cb(i - 1) - s_cb(i + 1))**2._wp

                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood
                    ! of beginning and end Riemann state extrapolation BC to avoid
                    ! any contributions from outside of the physical domain during
                    ! the WENO reconstruction
                    if (null_weights) then
                        if (bc_s%beg == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (1, 0) = 0._wp; d_cbR_${XYZ}$ (0, 0) = 1._wp
                            d_cbL_${XYZ}$ (1, 0) = 0._wp; d_cbL_${XYZ}$ (0, 0) = 1._wp
                        end if

                        if (bc_s%end == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (0, s) = 0._wp; d_cbR_${XYZ}$ (1, s) = 1._wp
                            d_cbL_${XYZ}$ (0, s) = 0._wp; d_cbL_${XYZ}$ (1, s) = 1._wp
                        end if
                    end if
                    ! END: Computing WENO3 Coefficients

                    ! Computing WENO5 Coefficients
                elseif (weno_order == 5) then

                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                            ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i)))/ &
                            ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i + 2) - s_cb(i)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 1) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, 1) = &
                            ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = &
                            ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                            ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = &
                            ((s_cb(i) - s_cb(i - 1))*(s_cb(i) - s_cb(i + 1)))/ &
                            ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 1) = &
                            ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, 1) = &
                            ((s_cb(i - 1) - s_cb(i))*(s_cb(i) - s_cb(i + 1)))/ &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 1) = &
                            ((s_cb(i) - s_cb(i + 2)) + (s_cb(i + 1) - s_cb(i + 3)))/ &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                            ((s_cb(i) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, 0) = &
                            ((s_cb(i - 2) - s_cb(i + 1)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                            ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 2)))* &
                            ((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 1) = &
                            ((s_cb(i) - s_cb(i + 2)) + (s_cb(i) - s_cb(i + 3)))/ &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                            ((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, 0) = &
                            ((s_cb(i - 2) - s_cb(i)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))* &
                            ((s_cb(i) - s_cb(i + 1)))

                        d_cbR_${XYZ}$ (0, i + 1) = &
                            ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbR_${XYZ}$ (2, i + 1) = &
                            ((s_cb(i + 1) - s_cb(i + 2))*(s_cb(i + 1) - s_cb(i + 3)))/ &
                            ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))
                        d_cbL_${XYZ}$ (0, i + 1) = &
                            ((s_cb(i - 2) - s_cb(i))*(s_cb(i) - s_cb(i - 1)))/ &
                            ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbL_${XYZ}$ (2, i + 1) = &
                            ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))/ &
                            ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))

                        d_cbR_${XYZ}$ (1, i + 1) = 1._wp - d_cbR_${XYZ}$ (0, i + 1) - d_cbR_${XYZ}$ (2, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1._wp - d_cbL_${XYZ}$ (0, i + 1) - d_cbL_${XYZ}$ (2, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                                                                     s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2._wp)/((s_cb(i) - &
                                                                                                                                                                          s_cb(i + 3))**2._wp*(s_cb(i + 1) - s_cb(i + 3))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 0, 1) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(19._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - &
                                                                                                                     s_cb(i + 1)) + 2._wp*(s_cb(i + 2) - s_cb(i))*((s_cb(i + 2) - &
                                                                                                                                                                    s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - &
                                                                                                                                                                                                               s_cb(i + 2))*(s_cb(i) - s_cb(i + 3))**2._wp*(s_cb(i + 3) - &
                                                                                                                                                                                                                                                            s_cb(i + 1)))

                        beta_coef_${XYZ}$ (i + 1, 0, 2) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - &
                                                                                                                      s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) + ((s_cb(i + 2) - &
                                                                                                                                                                  s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2._wp)/((s_cb(i) - &
                                                                                                                                                                                                                    s_cb(i + 2))**2._wp*(s_cb(i) - s_cb(i + 3))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 1, 0) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + (s_cb(i) - s_cb(i - 1))**2._wp + (s_cb(i) - &
                                                                                                                              s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) - &
                                                                                                                                                                      s_cb(i + 2))**2._wp*(s_cb(i) - s_cb(i + 2))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 1, 1) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*((s_cb(i) - &
                                                                   s_cb(i + 1))*((s_cb(i) - s_cb(i - 1)) + 20._wp*(s_cb(i + 1) - &
                                                                                                                   s_cb(i))) + (2._wp*(s_cb(i) - s_cb(i - 1)) + (s_cb(i + 1) - &
                                                                                                                                                                 s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - &
                                                                                                                                                                                                      s_cb(i - 1))*(s_cb(i - 1) - s_cb(i + 2))**2._wp*(s_cb(i + 2) - &
                                                                                                                                                                                                                                                       s_cb(i)))

                        beta_coef_${XYZ}$ (i + 1, 1, 2) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                                                                     s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2._wp)/ &
                            ((s_cb(i - 1) - s_cb(i + 1))**2._wp*(s_cb(i - 1) - &
                                                                 s_cb(i + 2))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 2, 0) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(12._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - &
                                                                                                                        s_cb(i - 1)))**2._wp + 3._wp*((s_cb(i) - s_cb(i - 2)) + &
                                                                                                                                                      (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/ &
                            ((s_cb(i - 2) - s_cb(i + 1))**2._wp*(s_cb(i - 1) - &
                                                                 s_cb(i + 1))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 2, 1) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(19._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - &
                                                                                                                      s_cb(i + 1))) + 2._wp*(s_cb(i + 1) - s_cb(i - 1))*((s_cb(i) - &
                                                                                                                                                                          s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) - &
                                                                                                                                                                                                                         s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2._wp*(s_cb(i + 1) - &
                                                                                                                                                                                                                                                                      s_cb(i - 1)))

                        beta_coef_${XYZ}$ (i + 1, 2, 2) = &
                            4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - &
                                                                          s_cb(i))**2._wp + (s_cb(i) - s_cb(i - 1))**2._wp + (s_cb(i) - &
                                                                                                                              s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) - &
                                                                                                                                                                      s_cb(i))**2._wp*(s_cb(i - 2) - s_cb(i + 1))**2._wp)

                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood
                    ! of beginning and end Riemann state extrapolation BC to avoid
                    ! any contributions from outside of the physical domain during
                    ! the WENO reconstruction
                    if (null_weights) then
                        if (bc_s%beg == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (1:2, 0) = 0._wp; d_cbR_${XYZ}$ (0, 0) = 1._wp
                            d_cbL_${XYZ}$ (1:2, 0) = 0._wp; d_cbL_${XYZ}$ (0, 0) = 1._wp
                            d_cbR_${XYZ}$ (2, 1) = 0._wp; d_cbR_${XYZ}$ (:, 1) = d_cbR_${XYZ}$ (:, 1)/sum(d_cbR_${XYZ}$ (:, 1))
                            d_cbL_${XYZ}$ (2, 1) = 0._wp; d_cbL_${XYZ}$ (:, 1) = d_cbL_${XYZ}$ (:, 1)/sum(d_cbL_${XYZ}$ (:, 1))
                        end if

                        if (bc_s%end == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (0, s - 1) = 0._wp; d_cbR_${XYZ}$ (:, s - 1) = d_cbR_${XYZ}$ (:, s - 1)/sum(d_cbR_${XYZ}$ (:, s - 1))
                            d_cbL_${XYZ}$ (0, s - 1) = 0._wp; d_cbL_${XYZ}$ (:, s - 1) = d_cbL_${XYZ}$ (:, s - 1)/sum(d_cbL_${XYZ}$ (:, s - 1))
                            d_cbR_${XYZ}$ (0:1, s) = 0._wp; d_cbR_${XYZ}$ (2, s) = 1._wp
                            d_cbL_${XYZ}$ (0:1, s) = 0._wp; d_cbL_${XYZ}$ (2, s) = 1._wp
                        end if
                    end if

                else ! WENO7

                    if (.not. teno) then

                        do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                            ! Reference: Shu (1997) "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws"
                            ! Equation 2.20: Polynomial Coefficients (poly_coef_cb)
                            ! Equation 2.61: Smoothness Indicators (beta_coef)
                            ! To reduce computational cost, we leverage the fact that all polynomial coefficients in a stencil sum to 1
                            ! and compute the polynomial coefficients (poly_coef_cb) for the cell value differences (dvd) instead of the values themselves.
                            ! The computation of coefficients is further simplified by using grid spacing (y or w) rather than the grid locations (s_cb) directly.
                            ! Ideal weights (d_cb) are obtained by comparing the grid location coefficients of the polynomial coefficients.
                            ! The smoothness indicators (beta_coef) are calculated through numerical differentiation and integration of each cross term of the polynomial coefficients,
                            ! using the cell value differences (dvd) instead of the values themselves.
                            ! While the polynomial coefficients sum to 1, the derivative of 1 is 0, which means it does not create additional cross terms in the smoothness indicators.

                            w = s_cb(i - 3:i + 4) - s_cb(i) ! Offset using s_cb(i) to reduce floating point error
                            d_cbR_${XYZ}$ (0, i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))) !&
                            d_cbR_${XYZ}$ (1, i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1)*w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7)*w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))) !&
                            d_cbR_${XYZ}$ (2, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2)*w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7)*w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))*(w(3) - w(8))) !&
                            d_cbR_${XYZ}$ (3, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8))*(w(3) - w(8))) !&

                            w = s_cb(i + 4:i - 3:-1) - s_cb(i)
                            d_cbL_${XYZ}$ (0, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8))*(w(3) - w(8))) !&
                            d_cbL_${XYZ}$ (1, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2)*w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7)*w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))*(w(3) - w(8))) !&
                            d_cbL_${XYZ}$ (2, i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1)*w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7)*w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))) !&
                            d_cbL_${XYZ}$ (3, i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))) !&
                            ! Note: Left has the reversed order of both points and coefficients compared to the right

                            y = s_cb(i + 1:i + 4) - s_cb(i:i + 3)
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = (y(1)*y(2)*(y(2) + y(3)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, 1) = -(y(1)*y(2)*(3*y(2)**2 + 6*y(2)*y(3) + 3*y(2)*y(4) + 2*y(1)*y(2) + 3*y(3)**2 + 3*y(3)*y(4) + 2*y(1)*y(3) + y(4)**2 + y(1)*y(4)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, 2) = (y(1)*(y(1)**2 + 3*y(1)*y(2) + 2*y(1)*y(3) + y(4)*y(1) + 3*y(2)**2 + 4*y(2)*y(3) + 2*y(4)*y(2) + y(3)**2 + y(4)*y(3)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i:i + 3) - s_cb(i - 1:i + 2)
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = -(y(2)*y(3)*(y(1) + y(2)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, 1) = (y(2)*(y(1) + y(2))*(y(2)**2 + 4*y(2)*y(3) + 2*y(2)*y(4) + y(1)*y(2) + 3*y(3)**2 + 3*y(3)*y(4) + 2*y(1)*y(3) + y(4)**2 + y(1)*y(4)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, 2) = (y(2)*y(3)*(y(3) + y(4)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i - 1:i + 2) - s_cb(i - 2:i + 1)
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, 0) = (y(3)*(y(2) + y(3))*(y(1) + y(2) + y(3)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, 1) = (y(3)*y(4)*(y(1)**2 + 3*y(1)*y(2) + 3*y(1)*y(3) + y(4)*y(1) + 3*y(2)**2 + 6*y(2)*y(3) + 2*y(4)*y(2) + 3*y(3)**2 + 2*y(4)*y(3)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, 2) = -(y(3)*y(4)*(y(2) + y(3)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i - 2:i + 1) - s_cb(i - 3:i)
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, 0) = (y(4)*(y(2)**2 + 4*y(2)*y(3) + 4*y(2)*y(4) + y(1)*y(2) + 3*y(3)**2 + 6*y(3)*y(4) + 2*y(1)*y(3) + 3*y(4)**2 + 2*y(1)*y(4)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, 1) = -(y(4)*(y(3) + y(4))*(y(1)**2 + 3*y(1)*y(2) + 3*y(1)*y(3) + 2*y(1)*y(4) + 3*y(2)**2 + 6*y(2)*y(3) + 4*y(2)*y(4) + 3*y(3)**2 + 4*y(3)*y(4) + y(4)**2))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, 2) = (y(4)*(y(3) + y(4))*(y(2) + y(3) + y(4)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i + 1:i - 2:-1) - s_cb(i:i - 3:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, 2) = (y(1)*y(2)*(y(2) + y(3)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, 1) = -(y(1)*y(2)*(3*y(2)**2 + 6*y(2)*y(3) + 3*y(2)*y(4) + 2*y(1)*y(2) + 3*y(3)**2 + 3*y(3)*y(4) + 2*y(1)*y(3) + y(4)**2 + y(1)*y(4)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, 0) = (y(1)*(y(1)**2 + 3*y(1)*y(2) + 2*y(1)*y(3) + y(4)*y(1) + 3*y(2)**2 + 4*y(2)*y(3) + 2*y(4)*y(2) + y(3)**2 + y(4)*y(3)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i + 2:i - 1:-1) - s_cb(i + 1:i - 2:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, 2) = -(y(2)*y(3)*(y(1) + y(2)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, 1) = (y(2)*(y(1) + y(2))*(y(2)**2 + 4*y(2)*y(3) + 2*y(2)*y(4) + y(1)*y(2) + 3*y(3)**2 + 3*y(3)*y(4) + 2*y(1)*y(3) + y(4)**2 + y(1)*y(4)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, 0) = (y(2)*y(3)*(y(3) + y(4)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i + 3:i:-1) - s_cb(i + 2:i - 1:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, 2) = (y(3)*(y(2) + y(3))*(y(1) + y(2) + y(3)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, 1) = (y(3)*y(4)*(y(1)**2 + 3*y(1)*y(2) + 3*y(1)*y(3) + y(4)*y(1) + 3*y(2)**2 + 6*y(2)*y(3) + 2*y(4)*y(2) + 3*y(3)**2 + 2*y(4)*y(3)))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = -(y(3)*y(4)*(y(2) + y(3)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            y = s_cb(i + 4:i + 1:-1) - s_cb(i + 3:i:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, 2) = (y(4)*(y(2)**2 + 4*y(2)*y(3) + 4*y(2)*y(4) + y(1)*y(2) + 3*y(3)**2 + 6*y(3)*y(4) + 2*y(1)*y(3) + 3*y(4)**2 + 2*y(1)*y(4)))/((y(3) + y(4))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, 1) = -(y(4)*(y(3) + y(4))*(y(1)**2 + 3*y(1)*y(2) + 3*y(1)*y(3) + 2*y(1)*y(4) + 3*y(2)**2 + 6*y(2)*y(3) + 4*y(2)*y(4) + 3*y(3)**2 + 4*y(3)*y(4) + y(4)**2))/((y(2) + y(3))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))) !&
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = (y(4)*(y(3) + y(4))*(y(2) + y(3) + y(4)))/((y(1) + y(2))*(y(1) + y(2) + y(3))*(y(1) + y(2) + y(3) + y(4))) !&

                            poly_coef_cbL_${XYZ}$ (i + 1, :, :) = -poly_coef_cbL_${XYZ}$ (i + 1, :, :)
                            ! Note: negative sign as the direction of taking the difference (dvd) is reversed

                            y = s_cb(i - 2:i + 1) - s_cb(i - 3:i)
                            beta_coef_${XYZ}$ (i + 1, 3, 0) = (4*y(4)**2*(5*y(1)**2*y(2)**2 + 20*y(1)**2*y(2)*y(3) + 15*y(1)**2*y(2)*y(4) + 20*y(1)**2*y(3)**2 + 30*y(1)**2*y(3)*y(4) + 60*y(1)**2*y(4)**2 + 10*y(1)*y(2)**3 + 60*y(1)*y(2)**2*y(3) + 45*y(1)*y(2)**2*y(4) + 110*y(1)*y(2)*y(3)**2 + 165*y(1)*y(2)*y(3)*y(4) & !&
                                                              + 260*y(1)*y(2)*y(4)**2 + 60*y(1)*y(3)**3 + 135*y(1)*y(3)**2*y(4) + 400*y(1)*y(3)*y(4)**2 + 225*y(1)*y(4)**3 + 5*y(2)**4 + 40*y(2)**3*y(3) + 30*y(2)**3*y(4) + 110*y(2)**2*y(3)**2 + 165*y(2)**2*y(3)*y(4) + 260*y(2)**2*y(4)**2 + 120*y(2)*y(3)**3 & !&
                                                              + 270*y(2)*y(3)**2*y(4) + 800*y(2)*y(3)*y(4)**2 + 450*y(2)*y(4)**3 + 45*y(3)**4 + 135*y(3)**3*y(4) + 600*y(3)**2*y(4)**2 + 675*y(3)*y(4)**3 + 996*y(4)**4))/(5*(y(3) + y(4))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 3, 1) = -(4*y(4)**2*(10*y(1)**3*y(2)*y(3) + 5*y(1)**3*y(2)*y(4) + 20*y(1)**3*y(3)**2 + 25*y(1)**3*y(3)*y(4) + 105*y(1)**3*y(4)**2 + 40*y(1)**2*y(2)**2*y(3) + 20*y(1)**2*y(2)**2*y(4) + 130*y(1)**2*y(2)*y(3)**2 + 155*y(1)**2*y(2)*y(3)*y(4) + 535*y(1)**2*y(2)*y(4)**2 & !&
                                                              + 90*y(1)**2*y(3)**3 + 165*y(1)**2*y(3)**2*y(4) + 790*y(1)**2*y(3)*y(4)**2 + 415*y(1)**2*y(4)**3 + 60*y(1)*y(2)**3*y(3) + 30*y(1)*y(2)**3*y(4) + 270*y(1)*y(2)**2*y(3)**2 + 315*y(1)*y(2)**2*y(3)*y(4) + 975*y(1)*y(2)**2*y(4)**2 + 360*y(1)*y(2)*y(3)**3 & !&
                                                              + 645*y(1)*y(2)*y(3)**2*y(4) + 2850*y(1)*y(2)*y(3)*y(4)**2 + 1460*y(1)*y(2)*y(4)**3 + 150*y(1)*y(3)**4 + 360*y(1)*y(3)**3*y(4) + 2000*y(1)*y(3)**2*y(4)**2 + 2005*y(1)*y(3)*y(4)**3 + 2077*y(1)*y(4)**4 + 30*y(2)**4*y(3) + 15*y(2)**4*y(4) + 180*y(2)**3*y(3)**2 & !&
                                                              + 210*y(2)**3*y(3)*y(4) + 650*y(2)**3*y(4)**2 + 360*y(2)**2*y(3)**3 + 645*y(2)**2*y(3)**2*y(4) + 2850*y(2)**2*y(3)*y(4)**2 + 1460*y(2)**2*y(4)**3 + 300*y(2)*y(3)**4 + 720*y(2)*y(3)**3*y(4) + 4000*y(2)*y(3)**2*y(4)**2 + 4010*y(2)*y(3)*y(4)**3 + 4154*y(2)*y(4)**4 & !&
                                                              + 90*y(3)**5 + 270*y(3)**4*y(4) + 1800*y(3)**3*y(4)**2 + 2655*y(3)**2*y(4)**3 + 4464*y(3)*y(4)**4 + 1767*y(4)**5))/(5*(y(2) + y(3))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 3, 2) = (4*y(4)**2*(10*y(2)**3*y(3) + 5*y(2)**3*y(4) + 50*y(2)**2*y(3)**2 + 60*y(2)**2*y(3)*y(4) + 10*y(1)*y(2)**2*y(3) + 215*y(2)**2*y(4)**2 + 5*y(1)*y(2)**2*y(4) + 70*y(2)*y(3)**3 + 130*y(2)*y(3)**2*y(4) + 30*y(1)*y(2)*y(3)**2 + 775*y(2)*y(3)*y(4)**2 & !&
                                                              + 35*y(1)*y(2)*y(3)*y(4) + 415*y(2)*y(4)**3 + 110*y(1)*y(2)*y(4)**2 + 30*y(3)**4 + 75*y(3)**3*y(4) + 20*y(1)*y(3)**3 + 665*y(3)**2*y(4)**2 + 35*y(1)*y(3)**2*y(4) + 725*y(3)*y(4)**3 + 220*y(1)*y(3)*y(4)**2 + 1767*y(4)**4 + 105*y(1)*y(4)**3)) & !&
                                                              /(5*(y(1) + y(2))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 3, 3) = (4*y(4)**2*(5*y(1)**4*y(3)**2 + 5*y(1)**4*y(3)*y(4) + 50*y(1)**4*y(4)**2 + 30*y(1)**3*y(2)*y(3)**2 + 30*y(1)**3*y(2)*y(3)*y(4) + 300*y(1)**3*y(2)*y(4)**2 + 30*y(1)**3*y(3)**3 + 45*y(1)**3*y(3)**2*y(4) + 415*y(1)**3*y(3)*y(4)**2 + 200*y(1)**3*y(4)**3 & !&
                                                              + 75*y(1)**2*y(2)**2*y(3)**2 + 75*y(1)**2*y(2)**2*y(3)*y(4) + 750*y(1)**2*y(2)**2*y(4)**2 + 150*y(1)**2*y(2)*y(3)**3 + 225*y(1)**2*y(2)*y(3)**2*y(4) + 2075*y(1)**2*y(2)*y(3)*y(4)**2 + 1000*y(1)**2*y(2)*y(4)**3 + 75*y(1)**2*y(3)**4 + 150*y(1)**2*y(3)**3*y(4) & !&
                                                              + 1390*y(1)**2*y(3)**2*y(4)**2 + 1315*y(1)**2*y(3)*y(4)**3 + 1081*y(1)**2*y(4)**4 + 90*y(1)*y(2)**3*y(3)**2 + 90*y(1)*y(2)**3*y(3)*y(4) + 900*y(1)*y(2)**3*y(4)**2 + 270*y(1)*y(2)**2*y(3)**3 + 405*y(1)*y(2)**2*y(3)**2*y(4) + 3735*y(1)*y(2)**2*y(3)*y(4)**2 & !&
                                                              + 1800*y(1)*y(2)**2*y(4)**3 + 270*y(1)*y(2)*y(3)**4 + 540*y(1)*y(2)*y(3)**3*y(4) + 5025*y(1)*y(2)*y(3)**2*y(4)**2 + 4755*y(1)*y(2)*y(3)*y(4)**3 + 4224*y(1)*y(2)*y(4)**4 + 90*y(1)*y(3)**5 + 225*y(1)*y(3)**4*y(4) + 2190*y(1)*y(3)**3*y(4)**2 + 3060*y(1)*y(3)**2*y(4)**3 & !&
                                                              + 4529*y(1)*y(3)*y(4)**4 + 1762*y(1)*y(4)**5 + 45*y(2)**4*y(3)**2 + 45*y(2)**4*y(3)*y(4) + 450*y(2)**4*y(4)**2 + 180*y(2)**3*y(3)**3 + 270*y(2)**3*y(3)**2*y(4) + 2490*y(2)**3*y(3)*y(4)**2 + 1200*y(2)**3*y(4)**3 + 270*y(2)**2*y(3)**4 + 540*y(2)**2*y(3)**3*y(4) & !&
                                                              + 5025*y(2)**2*y(3)**2*y(4)**2 + 4755*y(2)**2*y(3)*y(4)**3 + 4224*y(2)**2*y(4)**4 + 180*y(2)*y(3)**5 + 450*y(2)*y(3)**4*y(4) + 4380*y(2)*y(3)**3*y(4)**2 + 6120*y(2)*y(3)**2*y(4)**3 + 9058*y(2)*y(3)*y(4)**4 + 3524*y(2)*y(4)**5 + 45*y(3)**6 + 135*y(3)**5*y(4) & !&
                                                              + 1395*y(3)**4*y(4)**2 + 2565*y(3)**3*y(4)**3 + 4884*y(3)**2*y(4)**4 + 3624*y(3)*y(4)**5 + 831*y(4)**6))/(5*(y(2) + y(3))**2*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 3, 4) = -(4*y(4)**2*(10*y(1)**2*y(2)*y(3)**2 + 10*y(1)**2*y(2)*y(3)*y(4) + 100*y(1)**2*y(2)*y(4)**2 + 10*y(1)**2*y(3)**3 + 15*y(1)**2*y(3)**2*y(4) + 205*y(1)**2*y(3)*y(4)**2 + 100*y(1)**2*y(4)**3 + 30*y(1)*y(2)**2*y(3)**2 + 30*y(1)*y(2)**2*y(3)*y(4) & !&
                                                              + 300*y(1)*y(2)**2*y(4)**2 + 60*y(1)*y(2)*y(3)**3 + 90*y(1)*y(2)*y(3)**2*y(4) + 1030*y(1)*y(2)*y(3)*y(4)**2 + 500*y(1)*y(2)*y(4)**3 + 30*y(1)*y(3)**4 + 60*y(1)*y(3)**3*y(4) + 835*y(1)*y(3)**2*y(4)**2 + 805*y(1)*y(3)*y(4)**3 + 1762*y(1)*y(4)**4 + 30*y(2)**3*y(3)**2 & !&
                                                              + 30*y(2)**3*y(3)*y(4) + 300*y(2)**3*y(4)**2 + 90*y(2)**2*y(3)**3 + 135*y(2)**2*y(3)**2*y(4) + 1445*y(2)**2*y(3)*y(4)**2 + 700*y(2)**2*y(4)**3 + 90*y(2)*y(3)**4 + 180*y(2)*y(3)**3*y(4) + 2205*y(2)*y(3)**2*y(4)**2 + 2115*y(2)*y(3)*y(4)**3 + 3624*y(2)*y(4)**4 & !&
                                                              + 30*y(3)**5 + 75*y(3)**4*y(4) + 1060*y(3)**3*y(4)**2 + 1515*y(3)**2*y(4)**3 + 3824*y(3)*y(4)**4 + 1662*y(4)**5))/(5*(y(1) + y(2))*(y(2) + y(3))*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 3, 5) = (4*y(4)**2*(5*y(2)**2*y(3)**2 + 5*y(2)**2*y(3)*y(4) + 50*y(2)**2*y(4)**2 + 10*y(2)*y(3)**3 + 15*y(2)*y(3)**2*y(4) + 205*y(2)*y(3)*y(4)**2 + 100*y(2)*y(4)**3 + 5*y(3)**4 + 10*y(3)**3*y(4) + 205*y(3)**2*y(4)**2 + 200*y(3)*y(4)**3 + 831*y(4)**4))/(5*(y(1) & !&
                                                              + y(2))**2*(y(1) + y(2) + y(3))**2*(y(1) + y(2) + y(3) + y(4))**2) !&

                            y = s_cb(i - 1:i + 2) - s_cb(i - 2:i + 1)
                            beta_coef_${XYZ}$ (i + 1, 2, 0) = (4*y(3)**2*(5*y(1)**2*y(2)**2 + 5*y(1)**2*y(2)*y(3) + 50*y(1)**2*y(3)**2 + 10*y(1)*y(2)**3 + 15*y(1)*y(2)**2*y(3) + 205*y(1)*y(2)*y(3)**2 + 100*y(1)*y(3)**3 + 5*y(2)**4 + 10*y(2)**3*y(3) + 205*y(2)**2*y(3)**2 + 200*y(2)*y(3)**3 + 831*y(3)**4))/(5*(y(3) & !&
                                                              + y(4))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 2, 1) = (4*y(3)**2*(5*y(1)**3*y(2)*y(3) + 10*y(1)**3*y(2)*y(4) - 95*y(1)**3*y(3)**2 + 5*y(1)**3*y(3)*y(4) + 20*y(1)**2*y(2)**2*y(3) + 40*y(1)**2*y(2)**2*y(4) - 465*y(1)**2*y(2)*y(3)**2 + 55*y(1)**2*y(2)*y(3)*y(4) + 10*y(1)**2*y(2)*y(4)**2 - 285*y(1)**2*y(3)**3 & !&
                                                              + 20*y(1)**2*y(3)**2*y(4) + 5*y(1)**2*y(3)*y(4)**2 + 30*y(1)*y(2)**3*y(3) + 60*y(1)*y(2)**3*y(4) - 825*y(1)*y(2)**2*y(3)**2 + 135*y(1)*y(2)**2*y(3)*y(4) + 30*y(1)*y(2)**2*y(4)**2 - 1040*y(1)*y(2)*y(3)**3 + 100*y(1)*y(2)*y(3)**2*y(4) + 35*y(1)*y(2)*y(3)*y(4)**2 & !&
                                                              - 1847*y(1)*y(3)**4 + 125*y(1)*y(3)**3*y(4) + 110*y(1)*y(3)**2*y(4)**2 + 15*y(2)**4*y(3) + 30*y(2)**4*y(4) - 550*y(2)**3*y(3)**2 + 90*y(2)**3*y(3)*y(4) + 20*y(2)**3*y(4)**2 - 1040*y(2)**2*y(3)**3 + 100*y(2)**2*y(3)**2*y(4) + 35*y(2)**2*y(3)*y(4)**2 & !&
                                                              - 3694*y(2)*y(3)**4 + 250*y(2)*y(3)**3*y(4) + 220*y(2)*y(3)**2*y(4)**2 - 3219*y(3)**5 - 1452*y(3)**4*y(4) + 105*y(3)**3*y(4)**2))/(5*(y(2) + y(3))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 2, 2) = -(4*y(3)**2*(5*y(2)**3*y(3) - 95*y(2)*y(3)**3 - 190*y(2)**2*y(3)**2 + 10*y(2)**3*y(4) + 100*y(3)**3*y(4) - 1562*y(3)**4 - 95*y(1)*y(2)*y(3)**2 + 5*y(1)*y(2)**2*y(3) + 10*y(1)*y(2)**2*y(4) + 100*y(1)*y(3)**2*y(4) + 205*y(2)*y(3)**2*y(4) & !&
                                                              + 15*y(2)**2*y(3)*y(4) + 10*y(1)*y(2)*y(3)*y(4)))/(5*(y(1) + y(2))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 2, 3) = (4*y(3)**2*(50*y(1)**4*y(3)**2 + 5*y(1)**4*y(3)*y(4) + 5*y(1)**4*y(4)**2 + 300*y(1)**3*y(2)*y(3)**2 + 30*y(1)**3*y(2)*y(3)*y(4) + 30*y(1)**3*y(2)*y(4)**2 + 200*y(1)**3*y(3)**3 + 25*y(1)**3*y(3)**2*y(4) + 35*y(1)**3*y(3)*y(4)**2 + 10*y(1)**3*y(4)**3 & !&
                                                              + 750*y(1)**2*y(2)**2*y(3)**2 + 75*y(1)**2*y(2)**2*y(3)*y(4) + 75*y(1)**2*y(2)**2*y(4)**2 + 1000*y(1)**2*y(2)*y(3)**3 + 125*y(1)**2*y(2)*y(3)**2*y(4) + 175*y(1)**2*y(2)*y(3)*y(4)**2 + 50*y(1)**2*y(2)*y(4)**3 + 1081*y(1)**2*y(3)**4 - 50*y(1)**2*y(3)**3*y(4) & !&
                                                              - 10*y(1)**2*y(3)**2*y(4)**2 + 45*y(1)**2*y(3)*y(4)**3 + 5*y(1)**2*y(4)**4 + 900*y(1)*y(2)**3*y(3)**2 + 90*y(1)*y(2)**3*y(3)*y(4) + 90*y(1)*y(2)**3*y(4)**2 + 1800*y(1)*y(2)**2*y(3)**3 + 225*y(1)*y(2)**2*y(3)**2*y(4) + 315*y(1)*y(2)**2*y(3)*y(4)**2 & !&
                                                              + 90*y(1)*y(2)**2*y(4)**3 + 4224*y(1)*y(2)*y(3)**4 - 120*y(1)*y(2)*y(3)**3*y(4) + 25*y(1)*y(2)*y(3)**2*y(4)**2 + 165*y(1)*y(2)*y(3)*y(4)**3 + 20*y(1)*y(2)*y(4)**4 + 3324*y(1)*y(3)**5 + 1407*y(1)*y(3)**4*y(4) - 100*y(1)*y(3)**3*y(4)**2 + 70*y(1)*y(3)**2*y(4)**3 & !&
                                                              + 15*y(1)*y(3)*y(4)**4 + 450*y(2)**4*y(3)**2 + 45*y(2)**4*y(3)*y(4) + 45*y(2)**4*y(4)**2 + 1200*y(2)**3*y(3)**3 + 150*y(2)**3*y(3)**2*y(4) + 210*y(2)**3*y(3)*y(4)**2 + 60*y(2)**3*y(4)**3 + 4224*y(2)**2*y(3)**4 - 120*y(2)**2*y(3)**3*y(4) + 25*y(2)**2*y(3)**2*y(4)**2 & !&
                                                              + 165*y(2)**2*y(3)*y(4)**3 + 20*y(2)**2*y(4)**4 + 6648*y(2)*y(3)**5 + 2814*y(2)*y(3)**4*y(4) - 200*y(2)*y(3)**3*y(4)**2 + 140*y(2)*y(3)**2*y(4)**3 + 30*y(2)*y(3)*y(4)**4 + 3174*y(3)**6 + 3039*y(3)**5*y(4) + 771*y(3)**4*y(4)**2 + 135*y(3)**3*y(4)**3 + 60*y(3)**2*y(4)**4)) & !&
                                                              /(5*(y(2) + y(3))**2*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 2, 4) = -(4*y(3)**2*(100*y(1)**2*y(2)*y(3)**2 + 10*y(1)**2*y(2)*y(3)*y(4) + 10*y(1)**2*y(2)*y(4)**2 - 95*y(1)**2*y(3)**2*y(4) + 5*y(1)**2*y(3)*y(4)**2 + 300*y(1)*y(2)**2*y(3)**2 + 30*y(1)*y(2)**2*y(3)*y(4) + 30*y(1)*y(2)**2*y(4)**2 + 200*y(1)*y(2)*y(3)**3 & !&
                                                              - 260*y(1)*y(2)*y(3)**2*y(4) + 50*y(1)*y(2)*y(3)*y(4)**2 + 10*y(1)*y(2)*y(4)**3 + 1562*y(1)*y(3)**4 - 190*y(1)*y(3)**3*y(4) + 15*y(1)*y(3)**2*y(4)**2 + 5*y(1)*y(3)*y(4)**3 + 300*y(2)**3*y(3)**2 + 30*y(2)**3*y(3)*y(4) + 30*y(2)**3*y(4)**2 + 400*y(2)**2*y(3)**3 & !&
                                                              - 235*y(2)**2*y(3)**2*y(4) + 85*y(2)**2*y(3)*y(4)**2 + 20*y(2)**2*y(4)**3 + 3224*y(2)*y(3)**4 - 460*y(2)*y(3)**3*y(4) - 35*y(2)*y(3)**2*y(4)**2 + 25*y(2)*y(3)*y(4)**3 + 3124*y(3)**5 + 1467*y(3)**4*y(4) + 110*y(3)**3*y(4)**2 + 105*y(3)**2*y(4)**3)) & !&
                                                              /(5*(y(1) + y(2))*(y(2) + y(3))*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 2, 5) = (4*y(3)**2*(50*y(2)**2*y(3)**2 + 5*y(2)**2*y(3)*y(4) + 5*y(2)**2*y(4)**2 - 95*y(2)*y(3)**2*y(4) + 5*y(2)*y(3)*y(4)**2 + 781*y(3)**4 + 50*y(3)**2*y(4)**2))/(5*(y(1) + y(2))**2*(y(1) + y(2) + y(3))**2*(y(1) + y(2) + y(3) + y(4))**2) !&

                            y = s_cb(i:i + 3) - s_cb(i - 1:i + 2)
                            beta_coef_${XYZ}$ (i + 1, 1, 0) = (4*y(2)**2*(50*y(1)**2*y(2)**2 + 5*y(1)**2*y(2)*y(3) + 5*y(1)**2*y(3)**2 - 95*y(1)*y(2)**2*y(3) + 5*y(1)*y(2)*y(3)**2 + 781*y(2)**4 + 50*y(2)**2*y(3)**2))/(5*(y(3) + y(4))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 1, 1) = -(4*y(2)**2*(105*y(1)**3*y(2)**2 + 25*y(1)**3*y(2)*y(3) + 5*y(1)**3*y(2)*y(4) + 20*y(1)**3*y(3)**2 + 10*y(1)**3*y(3)*y(4) + 110*y(1)**2*y(2)**3 - 35*y(1)**2*y(2)**2*y(3) + 15*y(1)**2*y(2)**2*y(4) + 85*y(1)**2*y(2)*y(3)**2 + 50*y(1)**2*y(2)*y(3)*y(4) & !&
                                                              + 5*y(1)**2*y(2)*y(4)**2 + 30*y(1)**2*y(3)**3 + 30*y(1)**2*y(3)**2*y(4) + 10*y(1)**2*y(3)*y(4)**2 + 1467*y(1)*y(2)**4 - 460*y(1)*y(2)**3*y(3) - 190*y(1)*y(2)**3*y(4) - 235*y(1)*y(2)**2*y(3)**2 - 260*y(1)*y(2)**2*y(3)*y(4) - 95*y(1)*y(2)**2*y(4)**2 & !&
                                                              + 30*y(1)*y(2)*y(3)**3 + 30*y(1)*y(2)*y(3)**2*y(4) + 10*y(1)*y(2)*y(3)*y(4)**2 + 3124*y(2)**5 + 3224*y(2)**4*y(3) + 1562*y(2)**4*y(4) + 400*y(2)**3*y(3)**2 + 200*y(2)**3*y(3)*y(4) + 300*y(2)**2*y(3)**3 + 300*y(2)**2*y(3)**2*y(4) + 100*y(2)**2*y(3)*y(4)**2)) & !&
                                                              /(5*(y(2) + y(3))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 1, 2) = -(4*y(2)**2*(100*y(1)*y(2)**3 - 190*y(2)**2*y(3)**2 + 10*y(1)*y(3)**3 + 5*y(2)*y(3)**3 - 95*y(2)**3*y(3) - 1562*y(2)**4 + 15*y(1)*y(2)*y(3)**2 + 205*y(1)*y(2)**2*y(3) + 100*y(1)*y(2)**2*y(4) + 10*y(1)*y(3)**2*y(4) + 5*y(2)*y(3)**2*y(4) - 95*y(2)**2*y(3)*y(4) & !&
                                                              + 10*y(1)*y(2)*y(3)*y(4)))/(5*(y(1) + y(2))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 1, 3) = (4*y(2)**2*(60*y(1)**4*y(2)**2 + 30*y(1)**4*y(2)*y(3) + 15*y(1)**4*y(2)*y(4) + 20*y(1)**4*y(3)**2 + 20*y(1)**4*y(3)*y(4) + 5*y(1)**4*y(4)**2 + 135*y(1)**3*y(2)**3 + 140*y(1)**3*y(2)**2*y(3) + 70*y(1)**3*y(2)**2*y(4) + 165*y(1)**3*y(2)*y(3)**2 & !&
                                                              + 165*y(1)**3*y(2)*y(3)*y(4) + 45*y(1)**3*y(2)*y(4)**2 + 60*y(1)**3*y(3)**3 + 90*y(1)**3*y(3)**2*y(4) + 50*y(1)**3*y(3)*y(4)**2 + 10*y(1)**3*y(4)**3 + 771*y(1)**2*y(2)**4 - 200*y(1)**2*y(2)**3*y(3) - 100*y(1)**2*y(2)**3*y(4) + 25*y(1)**2*y(2)**2*y(3)**2 & !&
                                                              + 25*y(1)**2*y(2)**2*y(3)*y(4) - 10*y(1)**2*y(2)**2*y(4)**2 + 210*y(1)**2*y(2)*y(3)**3 + 315*y(1)**2*y(2)*y(3)**2*y(4) + 175*y(1)**2*y(2)*y(3)*y(4)**2 + 35*y(1)**2*y(2)*y(4)**3 + 45*y(1)**2*y(3)**4 + 90*y(1)**2*y(3)**3*y(4) + 75*y(1)**2*y(3)**2*y(4)**2 & !&
                                                              + 30*y(1)**2*y(3)*y(4)**3 + 5*y(1)**2*y(4)**4 + 3039*y(1)*y(2)**5 + 2814*y(1)*y(2)**4*y(3) + 1407*y(1)*y(2)**4*y(4) - 120*y(1)*y(2)**3*y(3)**2 - 120*y(1)*y(2)**3*y(3)*y(4) - 50*y(1)*y(2)**3*y(4)**2 + 150*y(1)*y(2)**2*y(3)**3 + 225*y(1)*y(2)**2*y(3)**2*y(4) & !&
                                                              + 125*y(1)*y(2)**2*y(3)*y(4)**2 + 25*y(1)*y(2)**2*y(4)**3 + 45*y(1)*y(2)*y(3)**4 + 90*y(1)*y(2)*y(3)**3*y(4) + 75*y(1)*y(2)*y(3)**2*y(4)**2 + 30*y(1)*y(2)*y(3)*y(4)**3 + 5*y(1)*y(2)*y(4)**4 + 3174*y(2)**6 + 6648*y(2)**5*y(3) + 3324*y(2)**5*y(4) & !&
                                                              + 4224*y(2)**4*y(3)**2 + 4224*y(2)**4*y(3)*y(4) + 1081*y(2)**4*y(4)**2 + 1200*y(2)**3*y(3)**3 + 1800*y(2)**3*y(3)**2*y(4) + 1000*y(2)**3*y(3)*y(4)**2 + 200*y(2)**3*y(4)**3 + 450*y(2)**2*y(3)**4 + 900*y(2)**2*y(3)**3*y(4) + 750*y(2)**2*y(3)**2*y(4)**2 & !&
                                                              + 300*y(2)**2*y(3)*y(4)**3 + 50*y(2)**2*y(4)**4))/(5*(y(2) + y(3))**2*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 1, 4) = (4*y(2)**2*(105*y(1)**2*y(2)**3 + 220*y(1)**2*y(2)**2*y(3) + 110*y(1)**2*y(2)**2*y(4) + 35*y(1)**2*y(2)*y(3)**2 + 35*y(1)**2*y(2)*y(3)*y(4) + 5*y(1)**2*y(2)*y(4)**2 + 20*y(1)**2*y(3)**3 + 30*y(1)**2*y(3)**2*y(4) + 10*y(1)**2*y(3)*y(4)**2 - 1452*y(1)*y(2)**4 & !&
                                                              + 250*y(1)*y(2)**3*y(3) + 125*y(1)*y(2)**3*y(4) + 100*y(1)*y(2)**2*y(3)**2 + 100*y(1)*y(2)**2*y(3)*y(4) + 20*y(1)*y(2)**2*y(4)**2 + 90*y(1)*y(2)*y(3)**3 + 135*y(1)*y(2)*y(3)**2*y(4) + 55*y(1)*y(2)*y(3)*y(4)**2 + 5*y(1)*y(2)*y(4)**3 + 30*y(1)*y(3)**4 & !&
                                                              + 60*y(1)*y(3)**3*y(4) + 40*y(1)*y(3)**2*y(4)**2 + 10*y(1)*y(3)*y(4)**3 - 3219*y(2)**5 - 3694*y(2)**4*y(3) - 1847*y(2)**4*y(4) - 1040*y(2)**3*y(3)**2 - 1040*y(2)**3*y(3)*y(4) - 285*y(2)**3*y(4)**2 - 550*y(2)**2*y(3)**3 - 825*y(2)**2*y(3)**2*y(4) & !&
                                                              - 465*y(2)**2*y(3)*y(4)**2 - 95*y(2)**2*y(4)**3 + 15*y(2)*y(3)**4 + 30*y(2)*y(3)**3*y(4) + 20*y(2)*y(3)**2*y(4)**2 + 5*y(2)*y(3)*y(4)**3))/(5*(y(1) + y(2))*(y(2) + y(3))*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 1, 5) = (4*y(2)**2*(831*y(2)**4 + 200*y(2)**3*y(3) + 100*y(2)**3*y(4) + 205*y(2)**2*y(3)**2 + 205*y(2)**2*y(3)*y(4) + 50*y(2)**2*y(4)**2 + 10*y(2)*y(3)**3 + 15*y(2)*y(3)**2*y(4) + 5*y(2)*y(3)*y(4)**2 + 5*y(3)**4 + 10*y(3)**3*y(4) + 5*y(3)**2*y(4)**2))/(5*(y(1) & !&
                                                              + y(2))**2*(y(1) + y(2) + y(3))**2*(y(1) + y(2) + y(3) + y(4))**2) !&

                            y = s_cb(i + 1:i + 4) - s_cb(i:i + 3)
                            beta_coef_${XYZ}$ (i + 1, 0, 0) = (4*y(1)**2*(831*y(1)**4 + 200*y(1)**3*y(2) + 100*y(1)**3*y(3) + 205*y(1)**2*y(2)**2 + 205*y(1)**2*y(2)*y(3) + 50*y(1)**2*y(3)**2 + 10*y(1)*y(2)**3 + 15*y(1)*y(2)**2*y(3) + 5*y(1)*y(2)*y(3)**2 + 5*y(2)**4 + 10*y(2)**3*y(3) + 5*y(2)**2*y(3)**2))/(5*(y(3) & !&
                                                              + y(4))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 0, 1) = -(4*y(1)**2*(1662*y(1)**5 + 3824*y(1)**4*y(2) + 3624*y(1)**4*y(3) + 1762*y(1)**4*y(4) + 1515*y(1)**3*y(2)**2 + 2115*y(1)**3*y(2)*y(3) + 805*y(1)**3*y(2)*y(4) + 700*y(1)**3*y(3)**2 + 500*y(1)**3*y(3)*y(4) + 100*y(1)**3*y(4)**2 + 1060*y(1)**2*y(2)**3 & !&
                                                              + 2205*y(1)**2*y(2)**2*y(3) + 835*y(1)**2*y(2)**2*y(4) + 1445*y(1)**2*y(2)*y(3)**2 + 1030*y(1)**2*y(2)*y(3)*y(4) + 205*y(1)**2*y(2)*y(4)**2 + 300*y(1)**2*y(3)**3 + 300*y(1)**2*y(3)**2*y(4) + 100*y(1)**2*y(3)*y(4)**2 + 75*y(1)*y(2)**4 + 180*y(1)*y(2)**3*y(3) & !&
                                                              + 60*y(1)*y(2)**3*y(4) + 135*y(1)*y(2)**2*y(3)**2 + 90*y(1)*y(2)**2*y(3)*y(4) + 15*y(1)*y(2)**2*y(4)**2 + 30*y(1)*y(2)*y(3)**3 + 30*y(1)*y(2)*y(3)**2*y(4) + 10*y(1)*y(2)*y(3)*y(4)**2 + 30*y(2)**5 + 90*y(2)**4*y(3) + 30*y(2)**4*y(4) + 90*y(2)**3*y(3)**2 & !&
                                                              + 60*y(2)**3*y(3)*y(4) + 10*y(2)**3*y(4)**2 + 30*y(2)**2*y(3)**3 + 30*y(2)**2*y(3)**2*y(4) + 10*y(2)**2*y(3)*y(4)**2))/(5*(y(2) + y(3))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 0, 2) = (4*y(1)**2*(1767*y(1)**4 + 725*y(1)**3*y(2) + 415*y(1)**3*y(3) + 105*y(4)*y(1)**3 + 665*y(1)**2*y(2)**2 + 775*y(1)**2*y(2)*y(3) + 220*y(4)*y(1)**2*y(2) + 215*y(1)**2*y(3)**2 + 110*y(4)*y(1)**2*y(3) + 75*y(1)*y(2)**3 + 130*y(1)*y(2)**2*y(3) + 35*y(4)*y(1)*y(2)**2 & !&
                                                              + 60*y(1)*y(2)*y(3)**2 + 35*y(4)*y(1)*y(2)*y(3) + 5*y(1)*y(3)**3 + 5*y(4)*y(1)*y(3)**2 + 30*y(2)**4 + 70*y(2)**3*y(3) + 20*y(4)*y(2)**3 + 50*y(2)**2*y(3)**2 + 30*y(4)*y(2)**2*y(3) + 10*y(2)*y(3)**3 + 10*y(4)*y(2)*y(3)**2)) & !&
                                                              /(5*(y(1) + y(2))*(y(3) + y(4))*(y(1) + y(2) + y(3))*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 0, 3) = (4*y(1)**2*(831*y(1)**6 + 3624*y(1)**5*y(2) + 3524*y(1)**5*y(3) + 1762*y(1)**5*y(4) + 4884*y(1)**4*y(2)**2 + 9058*y(1)**4*y(2)*y(3) + 4529*y(1)**4*y(2)*y(4) + 4224*y(1)**4*y(3)**2 + 4224*y(1)**4*y(3)*y(4) + 1081*y(1)**4*y(4)**2 + 2565*y(1)**3*y(2)**3 & !&
                                                              + 6120*y(1)**3*y(2)**2*y(3) + 3060*y(1)**3*y(2)**2*y(4) + 4755*y(1)**3*y(2)*y(3)**2 + 4755*y(1)**3*y(2)*y(3)*y(4) + 1315*y(1)**3*y(2)*y(4)**2 + 1200*y(1)**3*y(3)**3 + 1800*y(1)**3*y(3)**2*y(4) + 1000*y(1)**3*y(3)*y(4)**2 + 200*y(1)**3*y(4)**3 + 1395*y(1)**2*y(2)**4 & !&
                                                              + 4380*y(1)**2*y(2)**3*y(3) + 2190*y(1)**2*y(2)**3*y(4) + 5025*y(1)**2*y(2)**2*y(3)**2 + 5025*y(1)**2*y(2)**2*y(3)*y(4) + 1390*y(1)**2*y(2)**2*y(4)**2 + 2490*y(1)**2*y(2)*y(3)**3 + 3735*y(1)**2*y(2)*y(3)**2*y(4) + 2075*y(1)**2*y(2)*y(3)*y(4)**2 & !&
                                                              + 415*y(1)**2*y(2)*y(4)**3 + 450*y(1)**2*y(3)**4 + 900*y(1)**2*y(3)**3*y(4) + 750*y(1)**2*y(3)**2*y(4)**2 + 300*y(1)**2*y(3)*y(4)**3 + 50*y(1)**2*y(4)**4 + 135*y(1)*y(2)**5 + 450*y(1)*y(2)**4*y(3) + 225*y(1)*y(2)**4*y(4) + 540*y(1)*y(2)**3*y(3)**2 & !&
                                                              + 540*y(1)*y(2)**3*y(3)*y(4) + 150*y(1)*y(2)**3*y(4)**2 + 270*y(1)*y(2)**2*y(3)**3 + 405*y(1)*y(2)**2*y(3)**2*y(4) + 225*y(1)*y(2)**2*y(3)*y(4)**2 + 45*y(1)*y(2)**2*y(4)**3 + 45*y(1)*y(2)*y(3)**4 + 90*y(1)*y(2)*y(3)**3*y(4) + 75*y(1)*y(2)*y(3)**2*y(4)**2 & !&
                                                              + 30*y(1)*y(2)*y(3)*y(4)**3 + 5*y(1)*y(2)*y(4)**4 + 45*y(2)**6 + 180*y(2)**5*y(3) + 90*y(2)**5*y(4) + 270*y(2)**4*y(3)**2 + 270*y(2)**4*y(3)*y(4) + 75*y(2)**4*y(4)**2 + 180*y(2)**3*y(3)**3 + 270*y(2)**3*y(3)**2*y(4) + 150*y(2)**3*y(3)*y(4)**2 + 30*y(2)**3*y(4)**3 & !&
                                                              + 45*y(2)**2*y(3)**4 + 90*y(2)**2*y(3)**3*y(4) + 75*y(2)**2*y(3)**2*y(4)**2 + 30*y(2)**2*y(3)*y(4)**3 + 5*y(2)**2*y(4)**4))/(5*(y(2) + y(3))**2*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))**2*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 0, 4) = -(4*y(1)**2*(1767*y(1)**5 + 4464*y(1)**4*y(2) + 4154*y(1)**4*y(3) + 2077*y(1)**4*y(4) + 2655*y(1)**3*y(2)**2 + 4010*y(1)**3*y(2)*y(3) + 2005*y(1)**3*y(2)*y(4) + 1460*y(1)**3*y(3)**2 + 1460*y(1)**3*y(3)*y(4) + 415*y(1)**3*y(4)**2 + 1800*y(1)**2*y(2)**3 & !&
                                                              + 4000*y(1)**2*y(2)**2*y(3) + 2000*y(1)**2*y(2)**2*y(4) + 2850*y(1)**2*y(2)*y(3)**2 + 2850*y(1)**2*y(2)*y(3)*y(4) + 790*y(1)**2*y(2)*y(4)**2 + 650*y(1)**2*y(3)**3 + 975*y(1)**2*y(3)**2*y(4) + 535*y(1)**2*y(3)*y(4)**2 + 105*y(1)**2*y(4)**3 + 270*y(1)*y(2)**4 & !&
                                                              + 720*y(1)*y(2)**3*y(3) + 360*y(1)*y(2)**3*y(4) + 645*y(1)*y(2)**2*y(3)**2 + 645*y(1)*y(2)**2*y(3)*y(4) + 165*y(1)*y(2)**2*y(4)**2 + 210*y(1)*y(2)*y(3)**3 + 315*y(1)*y(2)*y(3)**2*y(4) + 155*y(1)*y(2)*y(3)*y(4)**2 + 25*y(1)*y(2)*y(4)**3 + 15*y(1)*y(3)**4 & !&
                                                              + 30*y(1)*y(3)**3*y(4) + 20*y(1)*y(3)**2*y(4)**2 + 5*y(1)*y(3)*y(4)**3 + 90*y(2)**5 + 300*y(2)**4*y(3) + 150*y(2)**4*y(4) + 360*y(2)**3*y(3)**2 + 360*y(2)**3*y(3)*y(4) + 90*y(2)**3*y(4)**2 + 180*y(2)**2*y(3)**3 + 270*y(2)**2*y(3)**2*y(4) + 130*y(2)**2*y(3)*y(4)**2 & !&
                                                              + 20*y(2)**2*y(4)**3 + 30*y(2)*y(3)**4 + 60*y(2)*y(3)**3*y(4) + 40*y(2)*y(3)**2*y(4)**2 + 10*y(2)*y(3)*y(4)**3))/(5*(y(1) + y(2))*(y(2) + y(3))*(y(1) + y(2) + y(3))**2*(y(2) + y(3) + y(4))*(y(1) + y(2) + y(3) + y(4))**2) !&
                            beta_coef_${XYZ}$ (i + 1, 0, 5) = (4*y(1)**2*(996*y(1)**4 + 675*y(1)**3*y(2) + 450*y(1)**3*y(3) + 225*y(1)**3*y(4) + 600*y(1)**2*y(2)**2 + 800*y(1)**2*y(2)*y(3) + 400*y(1)**2*y(2)*y(4) + 260*y(1)**2*y(3)**2 + 260*y(1)**2*y(3)*y(4) + 60*y(1)**2*y(4)**2 + 135*y(1)*y(2)**3 + 270*y(1)*y(2)**2*y(3) & !&
                                                              + 135*y(1)*y(2)**2*y(4) + 165*y(1)*y(2)*y(3)**2 + 165*y(1)*y(2)*y(3)*y(4) + 30*y(1)*y(2)*y(4)**2 + 30*y(1)*y(3)**3 + 45*y(1)*y(3)**2*y(4) + 15*y(1)*y(3)*y(4)**2 + 45*y(2)**4 + 120*y(2)**3*y(3) + 60*y(2)**3*y(4) + 110*y(2)**2*y(3)**2 + 110*y(2)**2*y(3)*y(4) & !&
                                                              + 20*y(2)**2*y(4)**2 + 40*y(2)*y(3)**3 + 60*y(2)*y(3)**2*y(4) + 20*y(2)*y(3)*y(4)**2 + 5*y(3)**4 + 10*y(3)**3*y(4) + 5*y(3)**2*y(4)**2))/(5*(y(1) + y(2))**2*(y(1) + y(2) + y(3))**2*(y(1) + y(2) + y(3) + y(4))**2) !&

                        end do

                    else ! TENO (only supports uniform grid)
                        ! (Fu, et al., 2016) Table 2 (for right flux)
                        d_cbL_${XYZ}$ (0, :) = 18._wp/35._wp
                        d_cbL_${XYZ}$ (1, :) = 3._wp/35._wp
                        d_cbL_${XYZ}$ (2, :) = 9._wp/35._wp
                        d_cbL_${XYZ}$ (3, :) = 1._wp/35._wp
                        d_cbL_${XYZ}$ (4, :) = 4._wp/35._wp

                        d_cbR_${XYZ}$ (0, :) = 18._wp/35._wp
                        d_cbR_${XYZ}$ (1, :) = 9._wp/35._wp
                        d_cbR_${XYZ}$ (2, :) = 3._wp/35._wp
                        d_cbR_${XYZ}$ (3, :) = 4._wp/35._wp
                        d_cbR_${XYZ}$ (4, :) = 1._wp/35._wp

                    end if
                end if

            end if
        #:endfor

        if (weno_dir == 1) then
            $:GPU_UPDATE(device='[poly_coef_cbL_x,poly_coef_cbR_x,d_cbL_x,d_cbR_x,beta_coef_x]')
        elseif (weno_dir == 2) then
            $:GPU_UPDATE(device='[poly_coef_cbL_y,poly_coef_cbR_y,d_cbL_y,d_cbR_y,beta_coef_y]')
        else
            $:GPU_UPDATE(device='[poly_coef_cbL_z,poly_coef_cbR_z,d_cbL_z,d_cbR_z,beta_coef_z]')
        end if

        ! Nullifying WENO coefficients and cell-boundary locations pointers

        nullify (s_cb)

    end subroutine s_compute_weno_coefficients

    subroutine s_weno(v_vf, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, &
                      weno_dir, &
                      is1_weno_d, is2_weno_d, is3_weno_d)

        type(scalar_field), dimension(1:), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in) :: weno_dir
        type(int_bounds_info), intent(in) :: is1_weno_d, is2_weno_d, is3_weno_d

        real(wp), dimension(-weno_polyn:weno_polyn - 1) :: dvd
        real(wp), dimension(0:weno_num_stencils) :: poly
        real(wp), dimension(0:weno_num_stencils) :: alpha
        real(wp), dimension(0:weno_num_stencils) :: omega
        real(wp), dimension(0:weno_num_stencils) :: beta
        real(wp), dimension(0:weno_num_stencils) :: delta
        real(wp), dimension(-3:3) :: v ! temporary field value array for clarity (WENO7 only)
        real(wp) :: tau

        integer :: i, j, k, l

        is1_weno = is1_weno_d
        is2_weno = is2_weno_d
        is3_weno = is3_weno_d

        $:GPU_UPDATE(device='[is1_weno,is2_weno,is3_weno]')

        if (weno_order /= 1) then
            call s_initialize_weno(v_vf, &
                                   weno_dir)
        end if

        if (weno_order == 1) then
            if (weno_dir == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                vL_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                vR_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
            else if (weno_dir == 2) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                vL_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                vR_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
            else if (weno_dir == 3) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                vL_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                vR_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
            end if
        elseif (weno_order == 3) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[beta,dvd,poly,omega,alpha,tau]')
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                do i = 1, v_size
                                    ! reconstruct from left side

                                    dvd(0) = v_rs_ws_${XYZ}$ (j + 1, k, l, i) &
                                             - v_rs_ws_${XYZ}$ (j, k, l, i)
                                    dvd(-1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              - v_rs_ws_${XYZ}$ (j - 1, k, l, i)

                                    poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbL_${XYZ}$ (j, 0, 0)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbL_${XYZ}$ (j, 1, 0)*dvd(-1)

                                    beta(0) = beta_coef_${XYZ}$ (j, 0, 0)*dvd(0)*dvd(0) &
                                              + weno_eps
                                    beta(1) = beta_coef_${XYZ}$ (j, 1, 0)*dvd(-1)*dvd(-1) &
                                              + weno_eps

                                    if (wenojs) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1._wp + d_cbL_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Borges, et al. (2008)

                                        tau = abs(beta(1) - beta(0))
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1._wp + tau/beta)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                    ! reconstruct from right side

                                    poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbR_${XYZ}$ (j, 0, 0)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbR_${XYZ}$ (j, 1, 0)*dvd(-1)

                                    if (wenojs) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1._wp + d_cbR_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then

                                        alpha = d_cbR_${XYZ}$ (:, j)*(1._wp + tau/beta)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                end do
                            end do
                        end do
                    end do
                end if
            #:endfor
        elseif (weno_order == 5) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[dvd,poly,beta,alpha,omega,tau,delta]')
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, v_size
                                    ! reconstruct from left side

                                    dvd(1) = v_rs_ws_${XYZ}$ (j + 2, k, l, i) &
                                             - v_rs_ws_${XYZ}$ (j + 1, k, l, i)
                                    dvd(0) = v_rs_ws_${XYZ}$ (j + 1, k, l, i) &
                                             - v_rs_ws_${XYZ}$ (j, k, l, i)
                                    dvd(-1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              - v_rs_ws_${XYZ}$ (j - 1, k, l, i)
                                    dvd(-2) = v_rs_ws_${XYZ}$ (j - 1, k, l, i) &
                                              - v_rs_ws_${XYZ}$ (j - 2, k, l, i)

                                    poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbL_${XYZ}$ (j, 0, 0)*dvd(1) &
                                              + poly_coef_cbL_${XYZ}$ (j, 0, 1)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbL_${XYZ}$ (j, 1, 0)*dvd(0) &
                                              + poly_coef_cbL_${XYZ}$ (j, 1, 1)*dvd(-1)
                                    poly(2) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbL_${XYZ}$ (j, 2, 0)*dvd(-1) &
                                              + poly_coef_cbL_${XYZ}$ (j, 2, 1)*dvd(-2)

                                    beta(0) = beta_coef_${XYZ}$ (j, 0, 0)*dvd(1)*dvd(1) &
                                              + beta_coef_${XYZ}$ (j, 0, 1)*dvd(1)*dvd(0) &
                                              + beta_coef_${XYZ}$ (j, 0, 2)*dvd(0)*dvd(0) &
                                              + weno_eps
                                    beta(1) = beta_coef_${XYZ}$ (j, 1, 0)*dvd(0)*dvd(0) &
                                              + beta_coef_${XYZ}$ (j, 1, 1)*dvd(0)*dvd(-1) &
                                              + beta_coef_${XYZ}$ (j, 1, 2)*dvd(-1)*dvd(-1) &
                                              + weno_eps
                                    beta(2) = beta_coef_${XYZ}$ (j, 2, 0)*dvd(-1)*dvd(-1) &
                                              + beta_coef_${XYZ}$ (j, 2, 1)*dvd(-1)*dvd(-2) &
                                              + beta_coef_${XYZ}$ (j, 2, 2)*dvd(-2)*dvd(-2) &
                                              + weno_eps

                                    if (wenojs) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1._wp + d_cbL_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Borges, et al. (2008)

                                        tau = abs(beta(2) - beta(0))                   ! Equation 25
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1._wp + tau/beta)  ! Equation 28 (note: weno_eps was already added to beta)

                                    elseif (teno) then
                                        ! Fu, et al. (2016)
                                        ! Fu''s code: https://dx.doi.org/10.13140/RG.2.2.36250.34247
                                        tau = abs(beta(2) - beta(0))
                                        alpha = 1._wp + tau/beta                    ! Equation 22 (reuse alpha as gamma; pick C=1 & q=6)
                                        alpha = (alpha*alpha*alpha)**2._wp          ! Equation 22 cont. (some CPU compilers cannot optimize x**6.0)
                                        omega = alpha/sum(alpha)                    ! Equation 25 (reuse omega as xi)
                                        delta = merge(0._wp, 1._wp, omega < teno_CT)! Equation 26
                                        alpha = delta*d_cbL_${XYZ}$ (:, j)          ! Equation 27

                                    end if

                                    omega = alpha/sum(alpha)

                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                    ! reconstruct from right side

                                    poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbR_${XYZ}$ (j, 0, 0)*dvd(1) &
                                              + poly_coef_cbR_${XYZ}$ (j, 0, 1)*dvd(0)
                                    poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbR_${XYZ}$ (j, 1, 0)*dvd(0) &
                                              + poly_coef_cbR_${XYZ}$ (j, 1, 1)*dvd(-1)
                                    poly(2) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                              + poly_coef_cbR_${XYZ}$ (j, 2, 0)*dvd(-1) &
                                              + poly_coef_cbR_${XYZ}$ (j, 2, 1)*dvd(-2)

                                    if (wenojs) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1._wp + d_cbR_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then

                                        alpha = d_cbR_${XYZ}$ (:, j)*(1._wp + tau/beta)

                                    elseif (teno) then
                                        alpha = delta*d_cbR_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                end do
                            end do
                        end do
                    end do

                    if (mp_weno) then
                        call s_preserve_monotonicity(v_rs_ws_${XYZ}$, vL_rs_vf_${XYZ}$, &
                                                     vR_rs_vf_${XYZ}$)
                    end if
                end if
            #:endfor
        elseif (weno_order == 7) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[poly,beta,alpha,omega,tau,delta,dvd,v]')
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, v_size

                                    if (teno) v = v_rs_ws_${XYZ}$ (j - 3:j + 3, k, l, i) ! temporary field value array for clarity

                                    if (.not. teno) then
                                        dvd(2) = v_rs_ws_${XYZ}$ (j + 3, k, l, i) &
                                                 - v_rs_ws_${XYZ}$ (j + 2, k, l, i)
                                        dvd(1) = v_rs_ws_${XYZ}$ (j + 2, k, l, i) &
                                                 - v_rs_ws_${XYZ}$ (j + 1, k, l, i)
                                        dvd(0) = v_rs_ws_${XYZ}$ (j + 1, k, l, i) &
                                                 - v_rs_ws_${XYZ}$ (j, k, l, i)
                                        dvd(-1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  - v_rs_ws_${XYZ}$ (j - 1, k, l, i)
                                        dvd(-2) = v_rs_ws_${XYZ}$ (j - 1, k, l, i) &
                                                  - v_rs_ws_${XYZ}$ (j - 2, k, l, i)
                                        dvd(-3) = v_rs_ws_${XYZ}$ (j - 2, k, l, i) &
                                                  - v_rs_ws_${XYZ}$ (j - 3, k, l, i)

                                        poly(3) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 0, 0)*dvd(2) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 0, 1)*dvd(1) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 0, 2)*dvd(0)
                                        poly(2) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 1, 0)*dvd(1) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 1, 1)*dvd(0) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 1, 2)*dvd(-1)
                                        poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 2, 0)*dvd(0) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 2, 1)*dvd(-1) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 2, 2)*dvd(-2)
                                        poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 3, 0)*dvd(-1) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 3, 1)*dvd(-2) &
                                                  + poly_coef_cbL_${XYZ}$ (j, 3, 2)*dvd(-3)

                                    else
                                        ! (Fu, et al., 2016) Table 1
                                        ! Note: Unlike TENO5, TENO7 stencils differ from WENO7 stencils
                                        ! See Figure 2 (right) for right-sided flux (at i+1/2)
                                        ! Here we need the left-sided flux, so we flip the weights with respect to the x=i point
                                        ! But we need to keep the stencil order to reuse the beta coefficients
                                        poly(0) = ( 2._wp*v(-1) +  5._wp*v( 0) -  1._wp*v( 1)) / 6._wp !&
                                        poly(1) = (11._wp*v( 0) -  7._wp*v( 1) +  2._wp*v( 2)) / 6._wp !&
                                        poly(2) = (-1._wp*v(-2) +  5._wp*v(-1) +  2._wp*v( 0)) / 6._wp !&
                                        poly(3) = (25._wp*v( 0) - 23._wp*v( 1) + 13._wp*v( 2) - 3._wp*v( 3)) / 12._wp !&
                                        poly(4) = ( 1._wp*v(-3) -  5._wp*v(-2) + 13._wp*v(-1) + 3._wp*v( 0)) / 12._wp !&
                                    end if

                                    if (.not. teno) then

                                        beta(3) = beta_coef_${XYZ}$ (j, 0, 0)*dvd(2)*dvd(2) &
                                                  + beta_coef_${XYZ}$ (j, 0, 1)*dvd(2)*dvd(1) &
                                                  + beta_coef_${XYZ}$ (j, 0, 2)*dvd(2)*dvd(0) &
                                                  + beta_coef_${XYZ}$ (j, 0, 3)*dvd(1)*dvd(1) &
                                                  + beta_coef_${XYZ}$ (j, 0, 4)*dvd(1)*dvd(0) &
                                                  + beta_coef_${XYZ}$ (j, 0, 5)*dvd(0)*dvd(0) &
                                                  + weno_eps

                                        beta(2) = beta_coef_${XYZ}$ (j, 1, 0)*dvd(1)*dvd(1) &
                                                  + beta_coef_${XYZ}$ (j, 1, 1)*dvd(1)*dvd(0) &
                                                  + beta_coef_${XYZ}$ (j, 1, 2)*dvd(1)*dvd(-1) &
                                                  + beta_coef_${XYZ}$ (j, 1, 3)*dvd(0)*dvd(0) &
                                                  + beta_coef_${XYZ}$ (j, 1, 4)*dvd(0)*dvd(-1) &
                                                  + beta_coef_${XYZ}$ (j, 1, 5)*dvd(-1)*dvd(-1) &
                                                  + weno_eps

                                        beta(1) = beta_coef_${XYZ}$ (j, 2, 0)*dvd(0)*dvd(0) &
                                                  + beta_coef_${XYZ}$ (j, 2, 1)*dvd(0)*dvd(-1) &
                                                  + beta_coef_${XYZ}$ (j, 2, 2)*dvd(0)*dvd(-2) &
                                                  + beta_coef_${XYZ}$ (j, 2, 3)*dvd(-1)*dvd(-1) &
                                                  + beta_coef_${XYZ}$ (j, 2, 4)*dvd(-1)*dvd(-2) &
                                                  + beta_coef_${XYZ}$ (j, 2, 5)*dvd(-2)*dvd(-2) &
                                                  + weno_eps

                                        beta(0) = beta_coef_${XYZ}$ (j, 3, 0)*dvd(-1)*dvd(-1) &
                                                  + beta_coef_${XYZ}$ (j, 3, 1)*dvd(-1)*dvd(-2) &
                                                  + beta_coef_${XYZ}$ (j, 3, 2)*dvd(-1)*dvd(-3) &
                                                  + beta_coef_${XYZ}$ (j, 3, 3)*dvd(-2)*dvd(-2) &
                                                  + beta_coef_${XYZ}$ (j, 3, 4)*dvd(-2)*dvd(-3) &
                                                  + beta_coef_${XYZ}$ (j, 3, 5)*dvd(-3)*dvd(-3) &
                                                  + weno_eps

                                    else ! TENO
                                        ! High-Order Low-Dissipation Targeted ENO Schemes for Ideal Magnetohydrodynamics (Fu & Tang, 2019) Section 3.2
                                        beta(0) = 13._wp/12._wp*(v(-1) - 2._wp*v( 0) + v( 1))**2._wp + ((    v(-1)             -     v( 1))**2._wp)/4._wp + weno_eps !&
                                        beta(1) = 13._wp/12._wp*(v( 0) - 2._wp*v( 1) + v( 2))**2._wp + ((3._wp*v( 0) - 4._wp*v( 1) +     v( 2))**2._wp)/4._wp + weno_eps !&
                                        beta(2) = 13._wp/12._wp*(v(-2) - 2._wp*v(-1) + v( 0))**2._wp + ((    v(-2) - 4._wp*v(-1) + 3._wp*v( 0))**2._wp)/4._wp + weno_eps !&

                                        beta(3) = ( v( 0)*(2107._wp*v( 0) - 9402._wp*v( 1) +  7042._wp*v( 2) - 1854._wp*v( 3)) & !&
                                                  + v( 1)*(              11003._wp*v( 1) - 17246._wp*v( 2) + 4642._wp*v( 3)) & !&
                                                  + v( 2)*(                               7043._wp*v( 2) - 3882._wp*v( 3)) & !&
                                                  + v( 3)*(                                               547._wp*v( 3)) ) / 240._wp & !&
                                                  + weno_eps !&

                                        beta(4) = ( v(-3)*(547._wp*v(-3) - 3882._wp*v(-2) +  4642._wp*v(-1) - 1854._wp*v( 0)) & !&
                                                  + v(-2)*(              7043._wp*v(-2) - 17246._wp*v(-1) + 7042._wp*v( 0)) & !&
                                                  + v(-1)*(                             11003._wp*v(-1) - 9402._wp*v( 0)) & !&
                                                  + v( 0)*(                                             2107._wp*v( 0)) ) / 240._wp & !&
                                                  + weno_eps !&
                                    end if

                                    if (wenojs) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1._wp + d_cbL_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Castro, et al. (2010)
                                        ! Don & Borges (2013) also helps
                                        tau = abs(beta(3) - beta(0)) ! Equation 50
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1._wp + (tau/beta)**wenoz_q) ! q = 2,3,4 for stability

                                    elseif (teno) then
                                        tau = abs(beta(4) - beta(3)) ! Note the reordering of stencils
                                        alpha = 1._wp + tau/beta
                                        alpha = (alpha*alpha*alpha)**2._wp ! some CPU compilers cannot optimize x**6.0
                                        omega = alpha/sum(alpha)
                                        delta = merge(0._wp, 1._wp, omega < teno_CT)
                                        alpha = delta*d_cbL_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                    if (.not. teno) then
                                        poly(3) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 0, 0)*dvd(2) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 0, 1)*dvd(1) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 0, 2)*dvd(0)
                                        poly(2) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 1, 0)*dvd(1) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 1, 1)*dvd(0) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 1, 2)*dvd(-1)
                                        poly(1) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 2, 0)*dvd(0) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 2, 1)*dvd(-1) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 2, 2)*dvd(-2)
                                        poly(0) = v_rs_ws_${XYZ}$ (j, k, l, i) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 3, 0)*dvd(-1) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 3, 1)*dvd(-2) &
                                                  + poly_coef_cbR_${XYZ}$ (j, 3, 2)*dvd(-3)
                                    else
                                        poly(0) = (-1._wp*v(-1) +  5._wp*v( 0) +  2._wp*v( 1)) / 6._wp !&
                                        poly(1) = ( 2._wp*v( 0) +  5._wp*v( 1) -  1._wp*v( 2)) / 6._wp !&
                                        poly(2) = ( 2._wp*v(-2) -  7._wp*v(-1) + 11._wp*v( 0)) / 6._wp !&
                                        poly(3) = ( 3._wp*v( 0) + 13._wp*v( 1) -  5._wp*v( 2) +  1._wp*v( 3)) / 12._wp !&
                                        poly(4) = (-3._wp*v(-3) + 13._wp*v(-2) - 23._wp*v(-1) + 25._wp*v( 0)) / 12._wp !&
                                    end if

                                    if (wenojs) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1._wp + d_cbR_${XYZ}$ (:, j) - 3._wp*omega) + omega**2._wp) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2._wp + omega*(1._wp - 2._wp*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        alpha = d_cbR_${XYZ}$ (:, j)*(1._wp + (tau/beta)**wenoz_q)

                                    elseif (teno) then
                                        alpha = delta*d_cbR_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                end do
                            end do
                        end do
                    end do

                end if
            #:endfor
        end if

        if (int_comp) then
            call s_interface_compression(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, &
                                         vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, &
                                         weno_dir, is1_weno_d, is2_weno_d, is3_weno_d)
        end if
    end subroutine s_weno

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      WENO reconstruction.
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param norm_dir Characteristic decommposition coordinate direction
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param is1_weno Index bounds in first coordinate direction
        !! @param is2_weno Index bounds in second coordinate direction
        !! @param is3_weno Index bounds in third coordinate direction
    subroutine s_initialize_weno(v_vf, &
                                 weno_dir)

        type(scalar_field), dimension(:), intent(IN) :: v_vf

        integer, intent(IN) :: weno_dir

        integer :: j, k, l, q

        ! Determining the number of cell-average variables which will be
        ! WENO-reconstructed and mapping their indical bounds in the x-,
        ! y- and z-directions to those in the s1-, s2- and s3-directions
        ! as to reshape the inputted data in the coordinate direction of
        ! the WENO reconstruction
        v_size = ubound(v_vf, 1)
        $:GPU_UPDATE(device='[v_size]')

        if (weno_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do j = 1, v_size
                do q = is3_weno%beg, is3_weno%end
                    do l = is2_weno%beg, is2_weno%end
                        do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                            v_rs_ws_x(k, l, q, j) = v_vf(j)%sf(k, l, q)
                        end do
                    end do
                end do
            end do
        end if

        ! Reshaping/Projecting onto Characteristic Fields in y-direction
        if (n == 0) return

        if (weno_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do j = 1, v_size
                do q = is3_weno%beg, is3_weno%end
                    do l = is2_weno%beg, is2_weno%end
                        do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                            v_rs_ws_y(k, l, q, j) = v_vf(j)%sf(l, k, q)
                        end do
                    end do
                end do
            end do
        end if

        ! Reshaping/Projecting onto Characteristic Fields in z-direction
        if (p == 0) return

        if (weno_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do j = 1, v_size
                do q = is3_weno%beg, is3_weno%end
                    do l = is2_weno%beg, is2_weno%end
                        do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                            v_rs_ws_z(k, l, q, j) = v_vf(j)%sf(q, l, k)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_initialize_weno

    !>  The goal of this subroutine is to ensure that the WENO
        !!      reconstruction is monotonic. The latter is achieved by
        !!      enforcing monotonicity preserving bounds of Suresh and
        !!      Huynh (1997). The resulting MPWENO reconstruction, see
        !!      Balsara and Shu (2000), ensures that the reconstructed
        !!      values do not reside outside the range spanned by WENO
        !!      stencil.
        !!  @param i Equation number
        !!  @param j First-coordinate cell index
        !!  @param k Secone-coordinate cell index
        !!  @param l Thire-coordinate cell index
    pure subroutine s_preserve_monotonicity(v_rs_ws, vL_rs_vf, vR_rs_vf)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(IN) :: v_rs_ws
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(INOUT) :: vL_rs_vf, vR_rs_vf

        integer :: i, j, k, l

        real(wp), dimension(-1:1) :: d !< Curvature measures at the zone centers

        real(wp) :: d_MD, d_LC !<
            !! Median (md) curvature and large curvature (LC) measures

        ! The left and right upper bounds (UL), medians, large curvatures,
        ! minima, and maxima of the WENO-reconstructed values of the cell-
        ! average variables.
        real(wp) :: vL_UL, vR_UL
        real(wp) :: vL_MD, vR_MD
        real(wp) :: vL_LC, vR_LC
        real(wp) :: vL_min, vR_min
        real(wp) :: vL_max, vR_max

        real(wp), parameter :: alpha = 2._wp !>
            !! Determines the maximum Courant–Friedrichs–Lewy (CFL) number that
            !! may be utilized with the scheme. In theory, for stability, a CFL
            !! number less than 1/(1+alpha) is necessary. The default value for
            !! alpha is 2.

        real(wp), parameter :: beta = 4._wp/3._wp !<
            !! Determines the amount of freedom available from utilizing a large
            !! value for the local curvature. The default value for beta is 4/3.

        real(wp), parameter :: alpha_mp = 2._wp
        real(wp), parameter :: beta_mp = 4._wp/3._wp

        $:GPU_PARALLEL_LOOP(collapse=4,private='[d]')
        do l = is3_weno%beg, is3_weno%end
            do k = is2_weno%beg, is2_weno%end
                do j = is1_weno%beg, is1_weno%end
                    do i = 1, v_size
                        d(-1) = v_rs_ws(j, k, l, i) &
                                + v_rs_ws(j - 2, k, l, i) &
                                - v_rs_ws(j - 1, k, l, i) &
                                *2._wp
                        d(0) = v_rs_ws(j + 1, k, l, i) &
                               + v_rs_ws(j - 1, k, l, i) &
                               - v_rs_ws(j, k, l, i) &
                               *2._wp
                        d(1) = v_rs_ws(j + 2, k, l, i) &
                               + v_rs_ws(j, k, l, i) &
                               - v_rs_ws(j + 1, k, l, i) &
                               *2._wp

                        d_MD = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1))) &
                               *abs((sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1))) &
                                    *(sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(0)))) &
                               *min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), &
                                    abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                        d_LC = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0))) &
                               *abs((sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(0))) &
                                    *(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(1)))) &
                               *min(abs(4._wp*d(0) - d(1)), abs(d(0)), &
                                    abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                        vL_UL = v_rs_ws(j, k, l, i) &
                                - (v_rs_ws(j + 1, k, l, i) &
                                   - v_rs_ws(j, k, l, i))*alpha_mp

                        vL_MD = (v_rs_ws(j, k, l, i) &
                                 + v_rs_ws(j - 1, k, l, i) &
                                 - d_MD)*5.e-1_wp

                        vL_LC = v_rs_ws(j, k, l, i) &
                                - (v_rs_ws(j + 1, k, l, i) &
                                   - v_rs_ws(j, k, l, i))*5.e-1_wp + beta_mp*d_LC

                        vL_min = max(min(v_rs_ws(j, k, l, i), &
                                         v_rs_ws(j - 1, k, l, i), &
                                         vL_MD), &
                                     min(v_rs_ws(j, k, l, i), &
                                         vL_UL, &
                                         vL_LC))

                        vL_max = min(max(v_rs_ws(j, k, l, i), &
                                         v_rs_ws(j - 1, k, l, i), &
                                         vL_MD), &
                                     max(v_rs_ws(j, k, l, i), &
                                         vL_UL, &
                                         vL_LC))

                        vL_rs_vf(j, k, l, i) = vL_rs_vf(j, k, l, i) &
                                               + (sign(5.e-1_wp, vL_min - vL_rs_vf(j, k, l, i)) &
                                                  + sign(5.e-1_wp, vL_max - vL_rs_vf(j, k, l, i))) &
                                               *min(abs(vL_min - vL_rs_vf(j, k, l, i)), &
                                                    abs(vL_max - vL_rs_vf(j, k, l, i)))
                        ! END: Left Monotonicity Preserving Bound

                        ! Right Monotonicity Preserving Bound
                        d(-1) = v_rs_ws(j, k, l, i) &
                                + v_rs_ws(j - 2, k, l, i) &
                                - v_rs_ws(j - 1, k, l, i) &
                                *2._wp
                        d(0) = v_rs_ws(j + 1, k, l, i) &
                               + v_rs_ws(j - 1, k, l, i) &
                               - v_rs_ws(j, k, l, i) &
                               *2._wp
                        d(1) = v_rs_ws(j + 2, k, l, i) &
                               + v_rs_ws(j, k, l, i) &
                               - v_rs_ws(j + 1, k, l, i) &
                               *2._wp

                        d_MD = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0))) &
                               *abs((sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(0))) &
                                    *(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(1)))) &
                               *min(abs(4._wp*d(0) - d(1)), abs(d(0)), &
                                    abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                        d_LC = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1))) &
                               *abs((sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1))) &
                                    *(sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(0)))) &
                               *min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), &
                                    abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                        vR_UL = v_rs_ws(j, k, l, i) &
                                + (v_rs_ws(j, k, l, i) &
                                   - v_rs_ws(j - 1, k, l, i))*alpha_mp

                        vR_MD = (v_rs_ws(j, k, l, i) &
                                 + v_rs_ws(j + 1, k, l, i) &
                                 - d_MD)*5.e-1_wp

                        vR_LC = v_rs_ws(j, k, l, i) &
                                + (v_rs_ws(j, k, l, i) &
                                   - v_rs_ws(j - 1, k, l, i))*5.e-1_wp + beta_mp*d_LC

                        vR_min = max(min(v_rs_ws(j, k, l, i), &
                                         v_rs_ws(j + 1, k, l, i), &
                                         vR_MD), &
                                     min(v_rs_ws(j, k, l, i), &
                                         vR_UL, &
                                         vR_LC))

                        vR_max = min(max(v_rs_ws(j, k, l, i), &
                                         v_rs_ws(j + 1, k, l, i), &
                                         vR_MD), &
                                     max(v_rs_ws(j, k, l, i), &
                                         vR_UL, &
                                         vR_LC))

                        vR_rs_vf(j, k, l, i) = vR_rs_vf(j, k, l, i) &
                                               + (sign(5.e-1_wp, vR_min - vR_rs_vf(j, k, l, i)) &
                                                  + sign(5.e-1_wp, vR_max - vR_rs_vf(j, k, l, i))) &
                                               *min(abs(vR_min - vR_rs_vf(j, k, l, i)), &
                                                    abs(vR_max - vR_rs_vf(j, k, l, i)))
                        ! END: Right Monotonicity Preserving Bound
                    end do
                end do
            end do
        end do

    end subroutine s_preserve_monotonicity

    !>  Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_weno_module()

        if (weno_order == 1) return

        ! Deallocating the WENO-stencil of the WENO-reconstructed variables

        !deallocate(vL_rs_vf_x, vR_rs_vf_x)
        @:DEALLOCATE(v_rs_ws_x)

        ! Deallocating WENO coefficients in x-direction
        @:DEALLOCATE(poly_coef_cbL_x, poly_coef_cbR_x)
        @:DEALLOCATE(d_cbL_x, d_cbR_x)
        @:DEALLOCATE(beta_coef_x)

        ! Deallocating WENO coefficients in y-direction
        if (n == 0) return

        !deallocate(vL_rs_vf_y, vR_rs_vf_y)
        @:DEALLOCATE(v_rs_ws_y)

        @:DEALLOCATE(poly_coef_cbL_y, poly_coef_cbR_y)
        @:DEALLOCATE(d_cbL_y, d_cbR_y)
        @:DEALLOCATE(beta_coef_y)

        ! Deallocating WENO coefficients in z-direction
        if (p == 0) return

        !deallocate(vL_rs_vf_z, vR_rs_vf_z)
        @:DEALLOCATE(v_rs_ws_z)

        @:DEALLOCATE(poly_coef_cbL_z, poly_coef_cbR_z)
        @:DEALLOCATE(d_cbL_z, d_cbR_z)
        @:DEALLOCATE(beta_coef_z)

    end subroutine s_finalize_weno_module

end module m_weno
