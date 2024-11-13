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
    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

#ifdef MFC_OPENACC
    use openacc
#endif

    use m_mpi_proxy
    ! ==========================================================================

    !implicit none

    private; public :: s_initialize_weno_module, s_initialize_weno, s_finalize_weno_module, s_weno

    !> @name The cell-average variables that will be WENO-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the WENO procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed WENO-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{

    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: v_rs_ws_x, v_rs_ws_y, v_rs_ws_z
    !> @}

    ! WENO Coefficients ========================================================

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at
    !! the left and right quadrature points (QP), in the x-, y- and z-directions.
    !! Note that the first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_z

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_z

    !    real(kind(0d0)), pointer, dimension(:, :, :) :: poly_coef_L => null()
    !    real(kind(0d0)), pointer, dimension(:, :, :) :: poly_coef_R => null()
    !> @}

    !> @name The ideal weights at the left and the right cell-boundaries and at the
    !! left and the right quadrature points, in x-, y- and z-directions. Note
    !! that the first dimension of the array identifies the weight, while the
    !! last denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_y
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_z

    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_y
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_z
!    real(kind(0d0)), pointer, dimension(:, :) :: d_L => null()
!    real(kind(0d0)), pointer, dimension(:, :) :: d_R => null()
    !> @}

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note
    !! that the first array dimension identifies the smoothness indicator, the
    !! second identifies the position of its coefficients and the last denotes
    !! the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_z
!    real(kind(0d0)), pointer, dimension(:, :, :) :: beta_coef => null()
    !> @}

    ! TODO GPU and improve implementation
    real(kind(0d0)), dimension(0:3,4) :: CpL
    real(kind(0d0)), dimension(0:3,4) :: CpR
    real(kind(0d0)), dimension(0:3,10) :: Cb

    real(kind(0d0)) :: w(1:8) ! Intermediate variables for overall stencil for ideal weights
    real(kind(0d0)) :: x(1:5) ! Intermediate variables for sub-stencils for poly and beta

    ! END: WENO Coefficients ===================================================

    integer :: v_size !< Number of WENO-reconstructed cell-average variables

    !$acc declare create(v_size)

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1_weno, is2_weno, is3_weno
    !$acc declare create(is1_weno, is2_weno, is3_weno)
    !
    !> @}

    real(kind(0d0)) :: test
    !$acc declare create(test)

    !$acc declare create( &
    !$acc                v_rs_ws_x, v_rs_ws_y, v_rs_ws_z, &
    !$acc                poly_coef_cbL_x,poly_coef_cbL_y,poly_coef_cbL_z, &
    !$acc                poly_coef_cbR_x,poly_coef_cbR_y,poly_coef_cbR_z,d_cbL_x,       &
    !$acc                d_cbL_y,d_cbL_z,d_cbR_x,d_cbR_y,d_cbR_z,beta_coef_x,beta_coef_y,beta_coef_z)

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_weno_module

        if (weno_order == 1) return

        ! Allocating/Computing WENO Coefficients in x-direction ============
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
            0:2*(weno_polyn - 1)))

        call s_compute_weno_coefficients(1, is1_weno)

        @:ALLOCATE(v_rs_ws_x(is1_weno%beg:is1_weno%end, &
            is2_weno%beg:is2_weno%end, is3_weno%beg:is3_weno%end, 1:sys_size))

        ! ==================================================================

        ! Allocating/Computing WENO Coefficients in y-direction ============
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
            0:2*(weno_polyn - 1)))

        call s_compute_weno_coefficients(2, is2_weno)

        @:ALLOCATE(v_rs_ws_y(is2_weno%beg:is2_weno%end, &
            is1_weno%beg:is1_weno%end, is3_weno%beg:is3_weno%end, 1:sys_size))

        ! ==================================================================

        ! Allocating/Computing WENO Coefficients in z-direction ============
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
            0:2*(weno_polyn - 1)))

        call s_compute_weno_coefficients(3, is3_weno)

        @:ALLOCATE(v_rs_ws_z(is3_weno%beg:is3_weno%end, &
            is2_weno%beg:is2_weno%end, is1_weno%beg:is1_weno%end, 1:sys_size))

        ! ==================================================================

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

        real(kind(0d0)), pointer, dimension(:) :: s_cb => null() !<
            !! Cell-boundary locations in the s-direction

        type(int_bounds_info) :: bc_s !< Boundary conditions (BC) in the s-direction

        integer :: i !< Generic loop iterator

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
            ! Computing WENO3 Coefficients =====================================
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

                        d_cbR_${XYZ}$ (1, i + 1) = 1d0 - d_cbR_${XYZ}$ (0, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1d0 - d_cbL_${XYZ}$ (0, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                                          (s_cb(i) - s_cb(i + 2))**2d0
                        beta_coef_${XYZ}$ (i + 1, 1, 0) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                                          (s_cb(i - 1) - s_cb(i + 1))**2d0

                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood
                    ! of beginning and end Riemann state extrapolation BC to avoid
                    ! any contributions from outside of the physical domain during
                    ! the WENO reconstruction
                    if (null_weights) then
                        if (bc_s%beg == -4) then
                            d_cbR_${XYZ}$ (1, 0) = 0d0; d_cbR_${XYZ}$ (0, 0) = 1d0
                            d_cbL_${XYZ}$ (1, 0) = 0d0; d_cbL_${XYZ}$ (0, 0) = 1d0
                        end if

                        if (bc_s%end == -4) then
                            d_cbR_${XYZ}$ (0, s) = 0d0; d_cbR_${XYZ}$ (1, s) = 1d0
                            d_cbL_${XYZ}$ (0, s) = 0d0; d_cbL_${XYZ}$ (1, s) = 1d0
                        end if
                    end if
                    ! END: Computing WENO3 Coefficients ================================

                    ! Computing WENO5 Coefficients =====================================
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

                        d_cbR_${XYZ}$ (1, i + 1) = 1d0 - d_cbR_${XYZ}$ (0, i + 1) - d_cbR_${XYZ}$ (2, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1d0 - d_cbL_${XYZ}$ (0, i + 1) - d_cbL_${XYZ}$ (2, i + 1)

                        beta_coef_${XYZ}$ (i + 1, 0, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                                                             s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/((s_cb(i) - &
                                                                                                                                                                s_cb(i + 3))**2d0*(s_cb(i + 1) - s_cb(i + 3))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 0, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - &
                                                                                                             s_cb(i + 1)) + 2d0*(s_cb(i + 2) - s_cb(i))*((s_cb(i + 2) - &
                                                                                                                                                          s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - &
                                                                                                                                                                                                     s_cb(i + 2))*(s_cb(i) - s_cb(i + 3))**2d0*(s_cb(i + 3) - &
                                                                                                                                                                                                                                                s_cb(i + 1)))

                        beta_coef_${XYZ}$ (i + 1, 0, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - &
                                                                                                              s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) + ((s_cb(i + 2) - &
                                                                                                                                                          s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2d0)/((s_cb(i) - &
                                                                                                                                                                                                          s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 3))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 1, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                                                                                                                    s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) - &
                                                                                                                                                            s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 2))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 1, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*((s_cb(i) - &
                                                               s_cb(i + 1))*((s_cb(i) - s_cb(i - 1)) + 20d0*(s_cb(i + 1) - &
                                                                                                             s_cb(i))) + (2d0*(s_cb(i) - s_cb(i - 1)) + (s_cb(i + 1) - &
                                                                                                                                                         s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - &
                                                                                                                                                                                              s_cb(i - 1))*(s_cb(i - 1) - s_cb(i + 2))**2d0*(s_cb(i + 2) - &
                                                                                                                                                                                                                                             s_cb(i)))

                        beta_coef_${XYZ}$ (i + 1, 1, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                                                                                                             s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/ &
                            ((s_cb(i - 1) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                               s_cb(i + 2))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 2, 0) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(12d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - &
                                                                                                                s_cb(i - 1)))**2d0 + 3d0*((s_cb(i) - s_cb(i - 2)) + &
                                                                                                                                          (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/ &
                            ((s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                               s_cb(i + 1))**2d0)

                        beta_coef_${XYZ}$ (i + 1, 2, 1) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - &
                                                                                                              s_cb(i + 1))) + 2d0*(s_cb(i + 1) - s_cb(i - 1))*((s_cb(i) - &
                                                                                                                                                                s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) - &
                                                                                                                                                                                                               s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i + 1) - &
                                                                                                                                                                                                                                                          s_cb(i - 1)))

                        beta_coef_${XYZ}$ (i + 1, 2, 2) = &
                            4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                                                                    s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                                                                                                                    s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) - &
                                                                                                                                                            s_cb(i))**2d0*(s_cb(i - 2) - s_cb(i + 1))**2d0)

                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood
                    ! of beginning and end Riemann state extrapolation BC to avoid
                    ! any contributions from outside of the physical domain during
                    ! the WENO reconstruction
                    if (null_weights) then
                        if (bc_s%beg == -4) then
                            d_cbR_${XYZ}$ (1:2, 0) = 0d0; d_cbR_${XYZ}$ (0, 0) = 1d0
                            d_cbL_${XYZ}$ (1:2, 0) = 0d0; d_cbL_${XYZ}$ (0, 0) = 1d0
                            d_cbR_${XYZ}$ (2, 1) = 0d0; d_cbR_${XYZ}$ (:, 1) = d_cbR_${XYZ}$ (:, 1)/sum(d_cbR_${XYZ}$ (:, 1))
                            d_cbL_${XYZ}$ (2, 1) = 0d0; d_cbL_${XYZ}$ (:, 1) = d_cbL_${XYZ}$ (:, 1)/sum(d_cbL_${XYZ}$ (:, 1))
                        end if

                        if (bc_s%end == -4) then
                            d_cbR_${XYZ}$ (0, s - 1) = 0d0; d_cbR_${XYZ}$ (:, s - 1) = d_cbR_${XYZ}$ (:, s - 1)/sum(d_cbR_${XYZ}$ (:, s - 1))
                            d_cbL_${XYZ}$ (0, s - 1) = 0d0; d_cbL_${XYZ}$ (:, s - 1) = d_cbL_${XYZ}$ (:, s - 1)/sum(d_cbL_${XYZ}$ (:, s - 1))
                            d_cbR_${XYZ}$ (0:1, s) = 0d0; d_cbR_${XYZ}$ (2, s) = 1d0
                            d_cbL_${XYZ}$ (0:1, s) = 0d0; d_cbL_${XYZ}$ (2, s) = 1d0
                        end if
                    end if

                else ! WENO7
                    ! Note: WENO7 only supports uniform grid
                    if (.not. teno) then
                        ! ! (Balsara & Shu, 2000) Page 11 Section III.a
                        ! d_cbL_${XYZ}$ (0, :) = 4d0/35d0
                        ! d_cbL_${XYZ}$ (1, :) = 18d0/35d0
                        ! d_cbL_${XYZ}$ (2, :) = 12d0/35d0
                        ! d_cbL_${XYZ}$ (3, :) = 1d0/35d0

                        ! d_cbR_${XYZ}$ (0, :) = 1d0/35d0
                        ! d_cbR_${XYZ}$ (1, :) = 12d0/35d0
                        ! d_cbR_${XYZ}$ (2, :) = 18d0/35d0
                        ! d_cbR_${XYZ}$ (3, :) = 4d0/35d0

                        do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn
                            w = s_cb(i + 4:i - 3:-1) ! Left has the reversed order of both points and coefficients compared to the right
                            d_cbL_${XYZ}$ (0, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8))*(w(3) - w(8)))
                            d_cbL_${XYZ}$ (1, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2)*w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7)*w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))*(w(3) - w(8)))
                            d_cbL_${XYZ}$ (2, i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1)*w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7)*w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8)))
                            d_cbL_${XYZ}$ (3, i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8)))

                            w = s_cb(i - 3:i + 4)
                            d_cbR_${XYZ}$ (0, i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8)))
                            d_cbR_${XYZ}$ (1, i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1)*w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7)*w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8)))
                            d_cbR_${XYZ}$ (2, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2)*w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7)*w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8))*(w(3) - w(8)))
                            d_cbR_${XYZ}$ (3, i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8))*(w(3) - w(8)))
                        end do

                        ! CpL(0,1) = 1d0/12d0
                        ! CpL(0,2) = -5d0/12d0
                        ! CpL(0,3) = 13d0/12d0
                        ! CpL(0,4) = 3d0/12d0

                        ! CpL(1,1) = -1d0/12d0
                        ! CpL(1,2) = 7d0/12d0
                        ! CpL(1,3) = 7d0/12d0
                        ! CpL(1,4) = -1d0/12d0

                        ! CpL(2,1) = 3d0/12d0
                        ! CpL(2,2) = 13d0/12d0
                        ! CpL(2,3) = -5d0/12d0
                        ! CpL(2,4) = 1d0/12d0

                        ! CpL(3,1) = 25d0/12d0
                        ! CpL(3,2) = -23d0/12d0
                        ! CpL(3,3) = 13d0/12d0
                        ! CpL(3,4) = -3d0/12d0


                        ! CpR(0,1) = -3d0/12d0
                        ! CpR(0,2) = 13d0/12d0
                        ! CpR(0,3) = -23d0/12d0
                        ! CpR(0,4) = 25d0/12d0

                        ! CpR(1,1) = 1d0/12d0
                        ! CpR(1,2) = -5d0/12d0
                        ! CpR(1,3) = 13d0/12d0
                        ! CpR(1,4) = 3d0/12d0

                        ! CpR(2,1) = -1d0/12d0
                        ! CpR(2,2) = 7d0/12d0
                        ! CpR(2,3) = 7d0/12d0
                        ! CpR(2,4) = -1d0/12d0

                        ! CpR(3,1) = 3d0/12d0
                        ! CpR(3,2) = 13d0/12d0
                        ! CpR(3,3) = -5d0/12d0
                        ! CpR(3,4) = 1d0/12d0


                        ! Cb(0,1) = 547d0/240d0
                        ! Cb(0,2) = -3882d0/240d0
                        ! Cb(0,3) = 4642d0/240d0
                        ! Cb(0,4) = -1854d0/240d0
                        ! Cb(0,5) = 7043d0/240d0
                        ! Cb(0,6) = -17246d0/240d0
                        ! Cb(0,7) = 7042d0/240d0
                        ! Cb(0,8) = 11003d0/240d0
                        ! Cb(0,9) = -9402d0/240d0
                        ! Cb(0,10) = 2107d0/240d0

                        ! Cb(1,1) = 267d0/240d0
                        ! Cb(1,2) = -1642d0/240d0
                        ! Cb(1,3) = 1602d0/240d0
                        ! Cb(1,4) = -494d0/240d0
                        ! Cb(1,5) = 2843d0/240d0
                        ! Cb(1,6) = -5966d0/240d0
                        ! Cb(1,7) = 1922d0/240d0
                        ! Cb(1,8) = 3443d0/240d0
                        ! Cb(1,9) = -2522d0/240d0
                        ! Cb(1,10) = 547d0/240d0

                        ! Cb(2,1) = 547d0/240d0
                        ! Cb(2,2) = -2522d0/240d0
                        ! Cb(2,3) = 1922d0/240d0
                        ! Cb(2,4) = -494d0/240d0
                        ! Cb(2,5) = 3443d0/240d0
                        ! Cb(2,6) = -5966d0/240d0
                        ! Cb(2,7) = 1602d0/240d0
                        ! Cb(2,8) = 2843d0/240d0
                        ! Cb(2,9) = -1642d0/240d0
                        ! Cb(2,10) = 267d0/240d0

                        ! Cb(3,1) = 2107d0/240d0
                        ! Cb(3,2) = -9402d0/240d0
                        ! Cb(3,3) = 7042d0/240d0
                        ! Cb(3,4) = -1854d0/240d0
                        ! Cb(3,5) = 11003d0/240d0
                        ! Cb(3,6) = -17246d0/240d0
                        ! Cb(3,7) = 4642d0/240d0
                        ! Cb(3,8) = 7043d0/240d0
                        ! Cb(3,9) = -3882d0/240d0
                        ! Cb(3,10) = 547d0/240d0

                        do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn

                            ! polynomial coefficients (right)

                            x = x_cb(i-3:i+1)
                            x = x - x(1)

                            CpR(0,1) = -((x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5)))

                            CpR(0,2) = ((x(2) - x(3))*((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)) + (x(1) - x(5))*(x(2) - x(5))*(x(4) - x(5)) + (x(1) - x(5))*(x(3) - x(5))*(x(4) - x(5)) + (x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))) + ((x(1) - x(5))*(x(2) - x(5))*(x(4) - x(5)))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5)))

                            CpR(0,3) = ((x(3) - x(4))*((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)) + (x(1) - x(5))*(x(2) - x(5))*(x(4) - x(5)) + (x(1) - x(5))*(x(3) - x(5))*(x(4) - x(5)) + (x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))) - ((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5)))

                            CpR(0,4) = ((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)) + (x(1) - x(5))*(x(2) - x(5))*(x(4) - x(5)) + (x(1) - x(5))*(x(3) - x(5))*(x(4) - x(5)) + (x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))

                            x = x_cb(i-2:i+2)
                            x = x - x(1)

                            CpR(1,1) = ((x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5)))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5)))

                            CpR(1,2) = ((x(2) - x(3))*((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5)) - (x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4)) + (x(1) - x(4))*(x(3) - x(4))*(x(4) - x(5)) + (x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) - ((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5)))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) + ((x(1) - x(4))*(x(2) - x(3))*(x(2) - x(4))*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))

                            CpR(1,3) = ((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5)) - (x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4)) + (x(1) - x(4))*(x(3) - x(4))*(x(4) - x(5)) + (x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) + ((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))**2)/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))

                            CpR(1,4) = ((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))

                            x = x_cb(i-1:i+3)
                            x = x - x(1)

                            CpR(2,1) = -((x(2) - x(3))*(x(3) - x(4))*(x(3) - x(5)))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5)))

                            CpR(2,2) = ((x(1) - x(3))*(x(2) - x(3))**2*(x(3) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) - ((x(1) - x(3))*(x(2) - x(3))*(x(3) - x(4)) + (x(1) - x(3))*(x(2) - x(3))*(x(3) - x(5)) - (x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5)) - (x(2) - x(3))*(x(3) - x(4))*(x(3) - x(5)))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(1) - x(3))*(x(2) - x(3))**2*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))

                            CpR(2,3) = ((x(1) - x(3))*(x(2) - x(3))*(x(3) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(1) - x(3))*(x(2) - x(3))*(x(3) - x(4))**2)/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))

                            CpR(2,4) = -((x(1) - x(3))*(x(2) - x(3))*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))

                            x = x_cb(i:i+4)
                            x = x - x(1)

                            CpR(3,1) = ((x(2) - x(3))*(x(2) - x(4))*(x(2) - x(5)))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5)))
                            
                            CpR(3,2) = ((x(1) - x(2))*(x(2) - x(4))*(x(2) - x(5)))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(1) - x(2))*(x(2) - x(3))**2*(x(2) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(1) - x(2))*(x(2) - x(3))**2*(x(2) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))
                            
                            CpR(3,3) = ((x(1) - x(2))*(x(2) - x(3))*(x(2) - x(4))*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))) - ((x(1) - x(2))*(x(2) - x(3))*(x(2) - x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5)))

                            CpR(3,4) = ((x(1) - x(2))*(x(2) - x(3))*(x(2) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5)))

                            ! polynomial coefficients (left)

                            CpL(0,1) = CpR(3,4)
                            CpL(0,2) = CpR(3,3)
                            CpL(0,3) = CpR(3,2)
                            CpL(0,4) = CpR(3,1)

                            CpL(1,1) = CpR(2,4)
                            CpL(1,2) = CpR(2,3)
                            CpL(1,3) = CpR(2,2)
                            CpL(1,4) = CpR(2,1)

                            CpL(2,1) = CpR(1,4)
                            CpL(2,2) = CpR(1,3)
                            CpL(2,3) = CpR(1,2)
                            CpL(2,4) = CpR(1,1)

                            CpL(3,1) = CpR(0,4)
                            CpL(3,2) = CpR(0,3)
                            CpL(3,3) = CpR(0,2)
                            CpL(3,4) = CpR(0,1)

                            ! beta coefficients

                            x = x_cb(i-3:i+1)
                            x = x - x(1)

                            Cb(0,1) = (4*(x(4) - x(5))**2*(5*x(2)**2*x(3)**2 - 5*x(2)**2*x(3)*x(4) - 5*x(2)**2*x(3)*x(5) + 50*x(2)**2*x(4)**2 - 95*x(2)**2*x(4)*x(5) + 50*x(2)**2*x(5)**2 - 5*x(2)*x(3)**2*x(4) - 5*x(2)*x(3)**2*x(5) + 105*x(2)*x(3)*x(4)**2 - 190*x(2)*x(3)*x(4)*x(5) + 105*x(2)*x(3)*x(5)**2 - 100*x(2)*x(4)**3 + 95*x(2)*x(4)**2*x(5) + 95*x(2)*x(4)*x(5)**2 - 100*x(2)*x(5)**3 + 50*x(3)**2*x(4)**2 - 95*x(3)**2*x(4)*x(5) + 50*x(3)**2*x(5)**2 - 100*x(3)*x(4)**3 + 95*x(3)*x(4)**2*x(5) + 95*x(3)*x(4)*x(5)**2 - 100*x(3)*x(5)**3 + 831*x(4)**4 - 3124*x(4)**3*x(5) + 4591*x(4)**2*x(5)**2 - 3124*x(4)*x(5)**3 + 831*x(5)**4))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2)

                            Cb(0,2) = (4*(x(4) - x(5))**2*(- 10*x(1)**3*x(2)*x(3)**2 + 10*x(1)**3*x(2)*x(3)*x(4) + 10*x(1)**3*x(2)*x(3)*x(5) - 100*x(1)**3*x(2)*x(4)**2 + 190*x(1)**3*x(2)*x(4)*x(5) - 100*x(1)**3*x(2)*x(5)**2 + 5*x(1)**3*x(3)**2*x(4) + 5*x(1)**3*x(3)**2*x(5) - 105*x(1)**3*x(3)*x(4)**2 + 190*x(1)**3*x(3)*x(4)*x(5) - 105*x(1)**3*x(3)*x(5)**2 + 100*x(1)**3*x(4)**3 - 95*x(1)**3*x(4)**2*x(5) - 95*x(1)**3*x(4)*x(5)**2 + 100*x(1)**3*x(5)**3 - 10*x(1)**2*x(2)**2*x(3)**2 + 10*x(1)**2*x(2)**2*x(3)*x(4) + 10*x(1)**2*x(2)**2*x(3)*x(5) - 100*x(1)**2*x(2)**2*x(4)**2 + 190*x(1)**2*x(2)**2*x(4)*x(5) - 100*x(1)**2*x(2)**2*x(5)**2 + 10*x(1)**2*x(2)*x(3)**3 + 10*x(1)**2*x(2)*x(3)**2*x(4) + 10*x(1)**2*x(2)*x(3)**2*x(5) - 120*x(1)**2*x(2)*x(3)*x(4)**2 + 170*x(1)**2*x(2)*x(3)*x(4)*x(5) - 120*x(1)**2*x(2)*x(3)*x(5)**2 + 300*x(1)**2*x(2)*x(4)**3 - 280*x(1)**2*x(2)*x(4)**2*x(5) - 280*x(1)**2*x(2)*x(4)*x(5)**2 + 300*x(1)**2*x(2)*x(5)**3 - 5*x(1)**2*x(3)**3*x(4) - 5*x(1)**2*x(3)**3*x(5) - 10*x(1)**2*x(3)**2*x(4)*x(5) + 205*x(1)**2*x(3)*x(4)**3 - 180*x(1)**2*x(3)*x(4)**2*x(5) - 180*x(1)**2*x(3)*x(4)*x(5)**2 + 205*x(1)**2*x(3)*x(5)**3 - 1762*x(1)**2*x(4)**4 + 6243*x(1)**2*x(4)**3*x(5) - 8992*x(1)**2*x(4)**2*x(5)**2 + 6243*x(1)**2*x(4)*x(5)**3 - 1762*x(1)**2*x(5)**4 - 10*x(1)*x(2)**3*x(3)**2 + 10*x(1)*x(2)**3*x(3)*x(4) + 10*x(1)*x(2)**3*x(3)*x(5) - 100*x(1)*x(2)**3*x(4)**2 + 190*x(1)*x(2)**3*x(4)*x(5) - 100*x(1)*x(2)**3*x(5)**2 + 10*x(1)*x(2)**2*x(3)**3 + 10*x(1)*x(2)**2*x(3)**2*x(4) + 10*x(1)*x(2)**2*x(3)**2*x(5) - 120*x(1)*x(2)**2*x(3)*x(4)**2 + 170*x(1)*x(2)**2*x(3)*x(4)*x(5) - 120*x(1)*x(2)**2*x(3)*x(5)**2 + 300*x(1)*x(2)**2*x(4)**3 - 280*x(1)*x(2)**2*x(4)**2*x(5) - 280*x(1)*x(2)**2*x(4)*x(5)**2 + 300*x(1)*x(2)**2*x(5)**3 - 20*x(1)*x(2)*x(3)**3*x(4) - 20*x(1)*x(2)*x(3)**3*x(5) + 110*x(1)*x(2)*x(3)**2*x(4)**2 - 200*x(1)*x(2)*x(3)**2*x(4)*x(5) + 110*x(1)*x(2)*x(3)**2*x(5)**2 + 110*x(1)*x(2)*x(3)*x(4)**3 - 70*x(1)*x(2)*x(3)*x(4)**2*x(5) - 70*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 110*x(1)*x(2)*x(3)*x(5)**3 - 1862*x(1)*x(2)*x(4)**4 + 6138*x(1)*x(2)*x(4)**3*x(5) - 8612*x(1)*x(2)*x(4)**2*x(5)**2 + 6138*x(1)*x(2)*x(4)*x(5)**3 - 1862*x(1)*x(2)*x(5)**4 + 105*x(1)*x(3)**3*x(4)**2 - 180*x(1)*x(3)**3*x(4)*x(5) + 105*x(1)*x(3)**3*x(5)**2 - 205*x(1)*x(3)**2*x(4)**3 + 190*x(1)*x(3)**2*x(4)**2*x(5) + 190*x(1)*x(3)**2*x(4)*x(5)**2 - 205*x(1)*x(3)**2*x(5)**3 + 1562*x(1)*x(3)*x(4)**4 - 6358*x(1)*x(3)*x(4)**3*x(5) + 9562*x(1)*x(3)*x(4)**2*x(5)**2 - 6358*x(1)*x(3)*x(4)*x(5)**3 + 1562*x(1)*x(3)*x(5)**4 + 1662*x(1)*x(4)**5 - 4486*x(1)*x(4)**4*x(5) + 2839*x(1)*x(4)**3*x(5)**2 + 2839*x(1)*x(4)**2*x(5)**3 - 4486*x(1)*x(4)*x(5)**4 + 1662*x(1)*x(5)**5 - 10*x(2)**4*x(3)**2 + 10*x(2)**4*x(3)*x(4) + 10*x(2)**4*x(3)*x(5) - 100*x(2)**4*x(4)**2 + 190*x(2)**4*x(4)*x(5) - 100*x(2)**4*x(5)**2 + 10*x(2)**3*x(3)**3 + 10*x(2)**3*x(3)**2*x(4) + 10*x(2)**3*x(3)**2*x(5) - 120*x(2)**3*x(3)*x(4)**2 + 170*x(2)**3*x(3)*x(4)*x(5) - 120*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 - 280*x(2)**3*x(4)**2*x(5) - 280*x(2)**3*x(4)*x(5)**2 + 300*x(2)**3*x(5)**3 - 20*x(2)**2*x(3)**3*x(4) - 20*x(2)**2*x(3)**3*x(5) + 110*x(2)**2*x(3)**2*x(4)**2 - 200*x(2)**2*x(3)**2*x(4)*x(5) + 110*x(2)**2*x(3)**2*x(5)**2 + 110*x(2)**2*x(3)*x(4)**3 - 70*x(2)**2*x(3)*x(4)**2*x(5) - 70*x(2)**2*x(3)*x(4)*x(5)**2 + 110*x(2)**2*x(3)*x(5)**3 - 1862*x(2)**2*x(4)**4 + 6138*x(2)**2*x(4)**3*x(5) - 8612*x(2)**2*x(4)**2*x(5)**2 + 6138*x(2)**2*x(4)*x(5)**3 - 1862*x(2)**2*x(5)**4 + 110*x(2)*x(3)**3*x(4)**2 - 160*x(2)*x(3)**3*x(4)*x(5) + 110*x(2)*x(3)**3*x(5)**2 - 310*x(2)*x(3)**2*x(4)**3 + 270*x(2)*x(3)**2*x(4)**2*x(5) + 270*x(2)*x(3)**2*x(4)*x(5)**2 - 310*x(2)*x(3)**2*x(5)**3 + 1662*x(2)*x(3)*x(4)**4 - 6358*x(2)*x(3)*x(4)**3*x(5) + 9372*x(2)*x(3)*x(4)**2*x(5)**2 - 6358*x(2)*x(3)*x(4)*x(5)**3 + 1662*x(2)*x(3)*x(5)**4 + 1662*x(2)*x(4)**5 - 4386*x(2)*x(4)**4*x(5) + 2744*x(2)*x(4)**3*x(5)**2 + 2744*x(2)*x(4)**2*x(5)**3 - 4386*x(2)*x(4)*x(5)**4 + 1662*x(2)*x(5)**5 - 100*x(3)**3*x(4)**3 + 85*x(3)**3*x(4)**2*x(5) + 85*x(3)**3*x(4)*x(5)**2 - 100*x(3)**3*x(5)**3 + 200*x(3)**2*x(4)**4 + 15*x(3)**2*x(4)**3*x(5) - 380*x(3)**2*x(4)**2*x(5)**2 + 15*x(3)**2*x(4)*x(5)**3 + 200*x(3)**2*x(5)**4 - 1662*x(3)*x(4)**5 + 4686*x(3)*x(4)**4*x(5) - 3029*x(3)*x(4)**3*x(5)**2 - 3029*x(3)*x(4)**2*x(5)**3 + 4686*x(3)*x(4)*x(5)**4 - 1662*x(3)*x(5)**5 - 1662*x(4)**5*x(5) + 6248*x(4)**4*x(5)**2 - 9182*x(4)**3*x(5)**3 + 6248*x(4)**2*x(5)**4 - 1662*x(4)*x(5)**5))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(0,3) = (4*(x(4) - x(5))**2*(10*x(1)**2*x(2)**3*x(3) - 5*x(1)**2*x(2)**3*x(4) - 5*x(1)**2*x(2)**3*x(5) + 10*x(1)**2*x(2)**2*x(3)**2 - 25*x(1)**2*x(2)**2*x(3)*x(4) - 25*x(1)**2*x(2)**2*x(3)*x(5) + 110*x(1)**2*x(2)**2*x(4)**2 - 180*x(1)**2*x(2)**2*x(4)*x(5) + 110*x(1)**2*x(2)**2*x(5)**2 + 10*x(1)**2*x(2)*x(3)**3 - 25*x(1)**2*x(2)*x(3)**2*x(4) - 25*x(1)**2*x(2)*x(3)**2*x(5) + 220*x(1)**2*x(2)*x(3)*x(4)**2 - 340*x(1)**2*x(2)*x(3)*x(4)*x(5) + 220*x(1)**2*x(2)*x(3)*x(5)**2 - 205*x(1)**2*x(2)*x(4)**3 + 175*x(1)**2*x(2)*x(4)**2*x(5) + 175*x(1)**2*x(2)*x(4)*x(5)**2 - 205*x(1)**2*x(2)*x(5)**3 - 5*x(1)**2*x(3)**3*x(4) - 5*x(1)**2*x(3)**3*x(5) + 110*x(1)**2*x(3)**2*x(4)**2 - 180*x(1)**2*x(3)**2*x(4)*x(5) + 110*x(1)**2*x(3)**2*x(5)**2 - 205*x(1)**2*x(3)*x(4)**3 + 175*x(1)**2*x(3)*x(4)**2*x(5) + 175*x(1)**2*x(3)*x(4)*x(5)**2 - 205*x(1)**2*x(3)*x(5)**3 + 100*x(1)**2*x(4)**4 + 10*x(1)**2*x(4)**3*x(5) - 190*x(1)**2*x(4)**2*x(5)**2 + 10*x(1)**2*x(4)*x(5)**3 + 100*x(1)**2*x(5)**4 + 10*x(1)*x(2)**3*x(3)**2 - 20*x(1)*x(2)**3*x(3)*x(4) - 20*x(1)*x(2)**3*x(3)*x(5) + 105*x(1)*x(2)**3*x(4)**2 - 180*x(1)*x(2)**3*x(4)*x(5) + 105*x(1)*x(2)**3*x(5)**2 + 10*x(1)*x(2)**2*x(3)**3 - 40*x(1)*x(2)**2*x(3)**2*x(4) - 40*x(1)*x(2)**2*x(3)**2*x(5) + 345*x(1)*x(2)**2*x(3)*x(4)**2 - 500*x(1)*x(2)**2*x(3)*x(4)*x(5) + 345*x(1)*x(2)**2*x(3)*x(5)**2 - 410*x(1)*x(2)**2*x(4)**3 + 350*x(1)*x(2)**2*x(4)**2*x(5) + 350*x(1)*x(2)**2*x(4)*x(5)**2 - 410*x(1)*x(2)**2*x(5)**3 - 20*x(1)*x(2)*x(3)**3*x(4) - 20*x(1)*x(2)*x(3)**3*x(5) + 345*x(1)*x(2)*x(3)**2*x(4)**2 - 500*x(1)*x(2)*x(3)**2*x(4)*x(5) + 345*x(1)*x(2)*x(3)**2*x(5)**2 - 830*x(1)*x(2)*x(3)*x(4)**3 + 670*x(1)*x(2)*x(3)*x(4)**2*x(5) + 670*x(1)*x(2)*x(3)*x(4)*x(5)**2 - 830*x(1)*x(2)*x(3)*x(5)**3 + 2067*x(1)*x(2)*x(4)**4 - 6208*x(1)*x(2)*x(4)**3*x(5) + 8452*x(1)*x(2)*x(4)**2*x(5)**2 - 6208*x(1)*x(2)*x(4)*x(5)**3 + 2067*x(1)*x(2)*x(5)**4 + 105*x(1)*x(3)**3*x(4)**2 - 180*x(1)*x(3)**3*x(4)*x(5) + 105*x(1)*x(3)**3*x(5)**2 - 410*x(1)*x(3)**2*x(4)**3 + 350*x(1)*x(3)**2*x(4)**2*x(5) + 350*x(1)*x(3)**2*x(4)*x(5)**2 - 410*x(1)*x(3)**2*x(5)**3 + 2067*x(1)*x(3)*x(4)**4 - 6208*x(1)*x(3)*x(4)**3*x(5) + 8452*x(1)*x(3)*x(4)**2*x(5)**2 - 6208*x(1)*x(3)*x(4)*x(5)**3 + 2067*x(1)*x(3)*x(5)**4 - 1762*x(1)*x(4)**5 + 4476*x(1)*x(4)**4*x(5) - 2754*x(1)*x(4)**3*x(5)**2 - 2754*x(1)*x(4)**2*x(5)**3 + 4476*x(1)*x(4)*x(5)**4 - 1762*x(1)*x(5)**5 + 10*x(2)**3*x(3)**3 - 20*x(2)**3*x(3)**2*x(4) - 20*x(2)**3*x(3)**2*x(5) + 110*x(2)**3*x(3)*x(4)**2 - 160*x(2)**3*x(3)*x(4)*x(5) + 110*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 + 85*x(2)**3*x(4)**2*x(5) + 85*x(2)**3*x(4)*x(5)**2 - 100*x(2)**3*x(5)**3 - 20*x(2)**2*x(3)**3*x(4) - 20*x(2)**2*x(3)**3*x(5) + 240*x(2)**2*x(3)**2*x(4)**2 - 320*x(2)**2*x(3)**2*x(4)*x(5) + 240*x(2)**2*x(3)**2*x(5)**2 - 520*x(2)**2*x(3)*x(4)**3 + 405*x(2)**2*x(3)*x(4)**2*x(5) + 405*x(2)**2*x(3)*x(4)*x(5)**2 - 520*x(2)**2*x(3)*x(5)**3 + 300*x(2)**2*x(4)**4 + 30*x(2)**2*x(4)**3*x(5) - 550*x(2)**2*x(4)**2*x(5)**2 + 30*x(2)**2*x(4)*x(5)**3 + 300*x(2)**2*x(5)**4 + 110*x(2)*x(3)**3*x(4)**2 - 160*x(2)*x(3)**3*x(4)*x(5) + 110*x(2)*x(3)**3*x(5)**2 - 520*x(2)*x(3)**2*x(4)**3 + 405*x(2)*x(3)**2*x(4)**2*x(5) + 405*x(2)*x(3)**2*x(4)*x(5)**2 - 520*x(2)*x(3)**2*x(5)**3 + 2272*x(2)*x(3)*x(4)**4 - 6178*x(2)*x(3)*x(4)**3*x(5) + 8122*x(2)*x(3)*x(4)**2*x(5)**2 - 6178*x(2)*x(3)*x(4)*x(5)**3 + 2272*x(2)*x(3)*x(5)**4 - 1862*x(2)*x(4)**5 + 4371*x(2)*x(4)**4*x(5) - 2579*x(2)*x(4)**3*x(5)**2 - 2579*x(2)*x(4)**2*x(5)**3 + 4371*x(2)*x(4)*x(5)**4 - 1862*x(2)*x(5)**5 - 100*x(3)**3*x(4)**3 + 85*x(3)**3*x(4)**2*x(5) + 85*x(3)**3*x(4)*x(5)**2 - 100*x(3)**3*x(5)**3 + 300*x(3)**2*x(4)**4 + 30*x(3)**2*x(4)**3*x(5) - 550*x(3)**2*x(4)**2*x(5)**2 + 30*x(3)**2*x(4)*x(5)**3 + 300*x(3)**2*x(5)**4 - 1862*x(3)*x(4)**5 + 4371*x(3)*x(4)**4*x(5) - 2579*x(3)*x(4)**3*x(5)**2 - 2579*x(3)*x(4)**2*x(5)**3 + 4371*x(3)*x(4)*x(5)**4 - 1862*x(3)*x(5)**5 + 1662*x(4)**6 - 4486*x(4)**5*x(5) + 4606*x(4)**4*x(5)**2 - 3504*x(4)**3*x(5)**3 + 4606*x(4)**2*x(5)**4 - 4486*x(4)*x(5)**5 + 1662*x(5)**6))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(0,4) = (4*(x(4) - x(5))**2*(- 10*x(2)**2*x(3)**2 + 10*x(2)**2*x(3)*x(4) + 20*x(2)**2*x(3)*x(5) - 10*x(1)*x(2)**2*x(3) - 100*x(2)**2*x(4)**2 + 185*x(2)**2*x(4)*x(5) + 5*x(1)*x(2)**2*x(4) - 105*x(2)**2*x(5)**2 + 5*x(1)*x(2)**2*x(5) + 10*x(2)*x(3)**2*x(4) + 20*x(2)*x(3)**2*x(5) - 10*x(1)*x(2)*x(3)**2 - 210*x(2)*x(3)*x(4)**2 + 365*x(2)*x(3)*x(4)*x(5) + 15*x(1)*x(2)*x(3)*x(4) - 235*x(2)*x(3)*x(5)**2 + 25*x(1)*x(2)*x(3)*x(5) + 200*x(2)*x(4)**3 - 85*x(2)*x(4)**2*x(5) - 105*x(1)*x(2)*x(4)**2 - 375*x(2)*x(4)*x(5)**2 + 185*x(1)*x(2)*x(4)*x(5) + 310*x(2)*x(5)**3 - 110*x(1)*x(2)*x(5)**2 - 100*x(3)**2*x(4)**2 + 185*x(3)**2*x(4)*x(5) + 5*x(1)*x(3)**2*x(4) - 105*x(3)**2*x(5)**2 + 5*x(1)*x(3)**2*x(5) + 200*x(3)*x(4)**3 - 85*x(3)*x(4)**2*x(5) - 105*x(1)*x(3)*x(4)**2 - 375*x(3)*x(4)*x(5)**2 + 185*x(1)*x(3)*x(4)*x(5) + 310*x(3)*x(5)**3 - 110*x(1)*x(3)*x(5)**2 - 1662*x(4)**4 + 6148*x(4)**3*x(5) + 100*x(1)*x(4)**3 - 9092*x(4)**2*x(5)**2 - 90*x(1)*x(4)**2*x(5) + 6343*x(4)*x(5)**3 - 95*x(1)*x(4)*x(5)**2 - 1767*x(5)**4 + 105*x(1)*x(5)**3))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(0,5) = x(4)*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(5)*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + x(4)*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(5)*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(4)**3*(((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (24*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))))*(x(4) - x(5)))/3 - (x(5)**3*(((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (24*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))))*(x(4) - x(5)))/3 - x(4)**2*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(5)**2*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + (576*(x(4) - x(5))**6*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) + (192*x(4)**3*(x(4) - x(5))**3*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (192*x(5)**3*(x(4) - x(5))**3*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) + (144*x(4)**5*(x(4) - x(5))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (144*x(5)**5*(x(4) - x(5))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (6*x(4)**4*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) + (6*x(5)**4*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) - (24*x(4)**2*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) + (24*x(5)**2*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(0,6) = (4*(x(4) - x(5))**2*(- 10*x(1)**5*x(2)**2*x(3) + 5*x(1)**5*x(2)**2*x(4) + 5*x(1)**5*x(2)**2*x(5) - 10*x(1)**5*x(2)*x(3)**2 + 20*x(1)**5*x(2)*x(3)*x(4) + 20*x(1)**5*x(2)*x(3)*x(5) - 105*x(1)**5*x(2)*x(4)**2 + 180*x(1)**5*x(2)*x(4)*x(5) - 105*x(1)**5*x(2)*x(5)**2 - 10*x(1)**5*x(3)**3 + 20*x(1)**5*x(3)**2*x(4) + 20*x(1)**5*x(3)**2*x(5) - 110*x(1)**5*x(3)*x(4)**2 + 160*x(1)**5*x(3)*x(4)*x(5) - 110*x(1)**5*x(3)*x(5)**2 + 100*x(1)**5*x(4)**3 - 85*x(1)**5*x(4)**2*x(5) - 85*x(1)**5*x(4)*x(5)**2 + 100*x(1)**5*x(5)**3 - 10*x(1)**4*x(2)**3*x(3) + 5*x(1)**4*x(2)**3*x(4) + 5*x(1)**4*x(2)**3*x(5) - 10*x(1)**4*x(2)**2*x(3)**2 + 50*x(1)**4*x(2)**2*x(3)*x(4) + 50*x(1)**4*x(2)**2*x(3)*x(5) - 220*x(1)**4*x(2)**2*x(4)**2 + 350*x(1)**4*x(2)**2*x(4)*x(5) - 220*x(1)**4*x(2)**2*x(5)**2 - 10*x(1)**4*x(2)*x(3)**3 + 50*x(1)**4*x(2)*x(3)**2*x(4) + 50*x(1)**4*x(2)*x(3)**2*x(5) - 370*x(1)**4*x(2)*x(3)*x(4)**2 + 440*x(1)**4*x(2)*x(3)*x(4)*x(5) - 370*x(1)**4*x(2)*x(3)*x(5)**2 + 615*x(1)**4*x(2)*x(4)**3 - 510*x(1)**4*x(2)*x(4)**2*x(5) - 510*x(1)**4*x(2)*x(4)*x(5)**2 + 615*x(1)**4*x(2)*x(5)**3 + 10*x(1)**4*x(3)**4 + 10*x(1)**4*x(3)**3*x(4) + 10*x(1)**4*x(3)**3*x(5) - 150*x(1)**4*x(3)**2*x(4)**2 + 120*x(1)**4*x(3)**2*x(4)*x(5) - 150*x(1)**4*x(3)**2*x(5)**2 + 530*x(1)**4*x(3)*x(4)**3 - 370*x(1)**4*x(3)*x(4)**2*x(5) - 370*x(1)**4*x(3)*x(4)*x(5)**2 + 530*x(1)**4*x(3)*x(5)**3 - 400*x(1)**4*x(4)**4 - 45*x(1)**4*x(4)**3*x(5) + 720*x(1)**4*x(4)**2*x(5)**2 - 45*x(1)**4*x(4)*x(5)**3 - 400*x(1)**4*x(5)**4 - 10*x(1)**3*x(2)**4*x(3) + 5*x(1)**3*x(2)**4*x(4) + 5*x(1)**3*x(2)**4*x(5) - 10*x(1)**3*x(2)**3*x(3)**2 + 50*x(1)**3*x(2)**3*x(3)*x(4) + 50*x(1)**3*x(2)**3*x(3)*x(5) - 220*x(1)**3*x(2)**3*x(4)**2 + 350*x(1)**3*x(2)**3*x(4)*x(5) - 220*x(1)**3*x(2)**3*x(5)**2 - 10*x(1)**3*x(2)**2*x(3)**3 + 50*x(1)**3*x(2)**2*x(3)**2*x(4) + 50*x(1)**3*x(2)**2*x(3)**2*x(5) - 500*x(1)**3*x(2)**2*x(3)*x(4)**2 + 550*x(1)**3*x(2)**2*x(3)*x(4)*x(5) - 500*x(1)**3*x(2)**2*x(3)*x(5)**2 + 930*x(1)**3*x(2)**2*x(4)**3 - 750*x(1)**3*x(2)**2*x(4)**2*x(5) - 750*x(1)**3*x(2)**2*x(4)*x(5)**2 + 930*x(1)**3*x(2)**2*x(5)**3 + 20*x(1)**3*x(2)*x(3)**4 - 10*x(1)**3*x(2)*x(3)**3*x(4) - 10*x(1)**3*x(2)*x(3)**3*x(5) - 170*x(1)**3*x(2)*x(3)**2*x(4)**2 + 70*x(1)**3*x(2)*x(3)**2*x(4)*x(5) - 170*x(1)**3*x(2)*x(3)**2*x(5)**2 + 1190*x(1)**3*x(2)*x(3)*x(4)**3 - 750*x(1)**3*x(2)*x(3)*x(4)**2*x(5) - 750*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 1190*x(1)**3*x(2)*x(3)*x(5)**3 - 2877*x(1)**3*x(2)*x(4)**4 + 5998*x(1)**3*x(2)*x(4)**3*x(5) - 6852*x(1)**3*x(2)*x(4)**2*x(5)**2 + 5998*x(1)**3*x(2)*x(4)*x(5)**3 - 2877*x(1)**3*x(2)*x(5)**4 - 30*x(1)**3*x(3)**4*x(4) - 30*x(1)**3*x(3)**4*x(5) + 130*x(1)**3*x(3)**3*x(4)**2 - 170*x(1)**3*x(3)**3*x(4)*x(5) + 130*x(1)**3*x(3)**3*x(5)**2 + 130*x(1)**3*x(3)**2*x(4)**3 - 10*x(1)**3*x(3)**2*x(4)**2*x(5) - 10*x(1)**3*x(3)**2*x(4)*x(5)**2 + 130*x(1)**3*x(3)**2*x(5)**3 - 2392*x(1)**3*x(3)*x(4)**4 + 5998*x(1)**3*x(3)*x(4)**3*x(5) - 7872*x(1)**3*x(3)*x(4)**2*x(5)**2 + 5998*x(1)**3*x(3)*x(4)*x(5)**3 - 2392*x(1)**3*x(3)*x(5)**4 + 2162*x(1)**3*x(4)**5 - 3941*x(1)**3*x(4)**4*x(5) + 1974*x(1)**3*x(4)**3*x(5)**2 + 1974*x(1)**3*x(4)**2*x(5)**3 - 3941*x(1)**3*x(4)*x(5)**4 + 2162*x(1)**3*x(5)**5 - 10*x(1)**2*x(2)**5*x(3) + 5*x(1)**2*x(2)**5*x(4) + 5*x(1)**2*x(2)**5*x(5) - 10*x(1)**2*x(2)**4*x(3)**2 + 50*x(1)**2*x(2)**4*x(3)*x(4) + 50*x(1)**2*x(2)**4*x(3)*x(5) - 220*x(1)**2*x(2)**4*x(4)**2 + 350*x(1)**2*x(2)**4*x(4)*x(5) - 220*x(1)**2*x(2)**4*x(5)**2 - 10*x(1)**2*x(2)**3*x(3)**3 + 50*x(1)**2*x(2)**3*x(3)**2*x(4) + 50*x(1)**2*x(2)**3*x(3)**2*x(5) - 500*x(1)**2*x(2)**3*x(3)*x(4)**2 + 550*x(1)**2*x(2)**3*x(3)*x(4)*x(5) - 500*x(1)**2*x(2)**3*x(3)*x(5)**2 + 930*x(1)**2*x(2)**3*x(4)**3 - 750*x(1)**2*x(2)**3*x(4)**2*x(5) - 750*x(1)**2*x(2)**3*x(4)*x(5)**2 + 930*x(1)**2*x(2)**3*x(5)**3 + 30*x(1)**2*x(2)**2*x(3)**4 - 30*x(1)**2*x(2)**2*x(3)**3*x(4) - 30*x(1)**2*x(2)**2*x(3)**3*x(5) - 60*x(1)**2*x(2)**2*x(3)**2*x(4)**2 - 90*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) - 60*x(1)**2*x(2)**2*x(3)**2*x(5)**2 + 1300*x(1)**2*x(2)**2*x(3)*x(4)**3 - 780*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) - 780*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 + 1300*x(1)**2*x(2)**2*x(3)*x(5)**3 - 3182*x(1)**2*x(2)**2*x(4)**4 + 5848*x(1)**2*x(2)**2*x(4)**3*x(5) - 6132*x(1)**2*x(2)**2*x(4)**2*x(5)**2 + 5848*x(1)**2*x(2)**2*x(4)*x(5)**3 - 3182*x(1)**2*x(2)**2*x(5)**4 - 60*x(1)**2*x(2)*x(3)**4*x(4) - 60*x(1)**2*x(2)*x(3)**4*x(5) + 390*x(1)**2*x(2)*x(3)**3*x(4)**2 - 450*x(1)**2*x(2)*x(3)**3*x(4)*x(5) + 390*x(1)**2*x(2)*x(3)**3*x(5)**2 - 290*x(1)**2*x(2)*x(3)**2*x(4)**3 + 330*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) + 330*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 - 290*x(1)**2*x(2)*x(3)**2*x(5)**3 - 2812*x(1)**2*x(2)*x(3)*x(4)**4 + 5658*x(1)**2*x(2)*x(3)*x(4)**3*x(5) - 6912*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 + 5658*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 - 2812*x(1)**2*x(2)*x(3)*x(5)**4 + 5991*x(1)**2*x(2)*x(4)**5 - 12148*x(1)**2*x(2)*x(4)**4*x(5) + 6622*x(1)**2*x(2)*x(4)**3*x(5)**2 + 6622*x(1)**2*x(2)*x(4)**2*x(5)**3 - 12148*x(1)**2*x(2)*x(4)*x(5)**4 + 5991*x(1)**2*x(2)*x(5)**5 + 130*x(1)**2*x(3)**4*x(4)**2 - 110*x(1)**2*x(3)**4*x(4)*x(5) + 130*x(1)**2*x(3)**4*x(5)**2 - 550*x(1)**2*x(3)**3*x(4)**3 + 350*x(1)**2*x(3)**3*x(4)**2*x(5) + 350*x(1)**2*x(3)**3*x(4)*x(5)**2 - 550*x(1)**2*x(3)**3*x(5)**3 + 1972*x(1)**2*x(3)**2*x(4)**4 - 6338*x(1)**2*x(3)**2*x(4)**3*x(5) + 8832*x(1)**2*x(3)**2*x(4)**2*x(5)**2 - 6338*x(1)**2*x(3)**2*x(4)*x(5)**3 + 1972*x(1)**2*x(3)**2*x(5)**4 + 1972*x(1)**2*x(3)*x(4)**5 - 3816*x(1)**2*x(3)*x(4)**4*x(5) + 2144*x(1)**2*x(3)*x(4)**3*x(5)**2 + 2144*x(1)**2*x(3)*x(4)**2*x(5)**3 - 3816*x(1)**2*x(3)*x(4)*x(5)**4 + 1972*x(1)**2*x(3)*x(5)**5 - 3524*x(1)**2*x(4)**6 + 6695*x(1)**2*x(4)**5*x(5) - 2844*x(1)**2*x(4)**4*x(5)**2 - 1104*x(1)**2*x(4)**3*x(5)**3 - 2844*x(1)**2*x(4)**2*x(5)**4 + 6695*x(1)**2*x(4)*x(5)**5 - 3524*x(1)**2*x(5)**6 - 10*x(1)*x(2)**5*x(3)**2 + 20*x(1)*x(2)**5*x(3)*x(4) + 20*x(1)*x(2)**5*x(3)*x(5) - 105*x(1)*x(2)**5*x(4)**2 + 180*x(1)*x(2)**5*x(4)*x(5) - 105*x(1)*x(2)**5*x(5)**2 - 10*x(1)*x(2)**4*x(3)**3 + 50*x(1)*x(2)**4*x(3)**2*x(4) + 50*x(1)*x(2)**4*x(3)**2*x(5) - 370*x(1)*x(2)**4*x(3)*x(4)**2 + 440*x(1)*x(2)**4*x(3)*x(4)*x(5) - 370*x(1)*x(2)**4*x(3)*x(5)**2 + 615*x(1)*x(2)**4*x(4)**3 - 510*x(1)*x(2)**4*x(4)**2*x(5) - 510*x(1)*x(2)**4*x(4)*x(5)**2 + 615*x(1)*x(2)**4*x(5)**3 + 20*x(1)*x(2)**3*x(3)**4 - 10*x(1)*x(2)**3*x(3)**3*x(4) - 10*x(1)*x(2)**3*x(3)**3*x(5) - 170*x(1)*x(2)**3*x(3)**2*x(4)**2 + 70*x(1)*x(2)**3*x(3)**2*x(4)*x(5) - 170*x(1)*x(2)**3*x(3)**2*x(5)**2 + 1190*x(1)*x(2)**3*x(3)*x(4)**3 - 750*x(1)*x(2)**3*x(3)*x(4)**2*x(5) - 750*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 + 1190*x(1)*x(2)**3*x(3)*x(5)**3 - 2877*x(1)*x(2)**3*x(4)**4 + 5998*x(1)*x(2)**3*x(4)**3*x(5) - 6852*x(1)*x(2)**3*x(4)**2*x(5)**2 + 5998*x(1)*x(2)**3*x(4)*x(5)**3 - 2877*x(1)*x(2)**3*x(5)**4 - 60*x(1)*x(2)**2*x(3)**4*x(4) - 60*x(1)*x(2)**2*x(3)**4*x(5) + 390*x(1)*x(2)**2*x(3)**3*x(4)**2 - 450*x(1)*x(2)**2*x(3)**3*x(4)*x(5) + 390*x(1)*x(2)**2*x(3)**3*x(5)**2 - 290*x(1)*x(2)**2*x(3)**2*x(4)**3 + 330*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) + 330*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 - 290*x(1)*x(2)**2*x(3)**2*x(5)**3 - 2812*x(1)*x(2)**2*x(3)*x(4)**4 + 5658*x(1)*x(2)**2*x(3)*x(4)**3*x(5) - 6912*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 + 5658*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 - 2812*x(1)*x(2)**2*x(3)*x(5)**4 + 5991*x(1)*x(2)**2*x(4)**5 - 12148*x(1)*x(2)**2*x(4)**4*x(5) + 6622*x(1)*x(2)**2*x(4)**3*x(5)**2 + 6622*x(1)*x(2)**2*x(4)**2*x(5)**3 - 12148*x(1)*x(2)**2*x(4)*x(5)**4 + 5991*x(1)*x(2)**2*x(5)**5 + 250*x(1)*x(2)*x(3)**4*x(4)**2 - 260*x(1)*x(2)*x(3)**4*x(4)*x(5) + 250*x(1)*x(2)*x(3)**4*x(5)**2 - 1090*x(1)*x(2)*x(3)**3*x(4)**3 + 710*x(1)*x(2)*x(3)**3*x(4)**2*x(5) + 710*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 - 1090*x(1)*x(2)*x(3)**3*x(5)**3 + 4254*x(1)*x(2)*x(3)**2*x(4)**4 - 12586*x(1)*x(2)*x(3)**2*x(4)**3*x(5) + 17154*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 - 12586*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 + 4254*x(1)*x(2)*x(3)**2*x(5)**4 + 310*x(1)*x(2)*x(3)*x(4)**5 + 1190*x(1)*x(2)*x(3)*x(4)**4*x(5) - 1110*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 - 1110*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 + 1190*x(1)*x(2)*x(3)*x(4)*x(5)**4 + 310*x(1)*x(2)*x(3)*x(5)**5 - 5286*x(1)*x(2)*x(4)**6 + 7442*x(1)*x(2)*x(4)**5*x(5) + 5096*x(1)*x(2)*x(4)**4*x(5)**2 - 15254*x(1)*x(2)*x(4)**3*x(5)**3 + 5096*x(1)*x(2)*x(4)**2*x(5)**4 + 7442*x(1)*x(2)*x(4)*x(5)**5 - 5286*x(1)*x(2)*x(5)**6 - 210*x(1)*x(3)**4*x(4)**3 + 120*x(1)*x(3)**4*x(4)**2*x(5) + 120*x(1)*x(3)**4*x(4)*x(5)**2 - 210*x(1)*x(3)**4*x(5)**3 + 720*x(1)*x(3)**3*x(4)**4 + 150*x(1)*x(3)**3*x(4)**3*x(5) - 1170*x(1)*x(3)**3*x(4)**2*x(5)**2 + 150*x(1)*x(3)**3*x(4)*x(5)**3 + 720*x(1)*x(3)**3*x(5)**4 - 3834*x(1)*x(3)**2*x(4)**5 + 8812*x(1)*x(3)**2*x(4)**4*x(5) - 5218*x(1)*x(3)**2*x(4)**3*x(5)**2 - 5218*x(1)*x(3)**2*x(4)**2*x(5)**3 + 8812*x(1)*x(3)**2*x(4)*x(5)**4 - 3834*x(1)*x(3)**2*x(5)**5 + 1662*x(1)*x(3)*x(4)**6 - 6558*x(1)*x(3)*x(4)**5*x(5) + 10804*x(1)*x(3)*x(4)**4*x(5)**2 - 11986*x(1)*x(3)*x(4)**3*x(5)**3 + 10804*x(1)*x(3)*x(4)**2*x(5)**4 - 6558*x(1)*x(3)*x(4)*x(5)**5 + 1662*x(1)*x(3)*x(5)**6 + 1662*x(1)*x(4)**7 - 962*x(1)*x(4)**6*x(5) - 4251*x(1)*x(4)**5*x(5)**2 + 3681*x(1)*x(4)**4*x(5)**3 + 3681*x(1)*x(4)**3*x(5)**4 - 4251*x(1)*x(4)**2*x(5)**5 - 962*x(1)*x(4)*x(5)**6 + 1662*x(1)*x(5)**7 - 10*x(2)**5*x(3)**3 + 20*x(2)**5*x(3)**2*x(4) + 20*x(2)**5*x(3)**2*x(5) - 110*x(2)**5*x(3)*x(4)**2 + 160*x(2)**5*x(3)*x(4)*x(5) - 110*x(2)**5*x(3)*x(5)**2 + 100*x(2)**5*x(4)**3 - 85*x(2)**5*x(4)**2*x(5) - 85*x(2)**5*x(4)*x(5)**2 + 100*x(2)**5*x(5)**3 + 10*x(2)**4*x(3)**4 + 10*x(2)**4*x(3)**3*x(4) + 10*x(2)**4*x(3)**3*x(5) - 150*x(2)**4*x(3)**2*x(4)**2 + 120*x(2)**4*x(3)**2*x(4)*x(5) - 150*x(2)**4*x(3)**2*x(5)**2 + 530*x(2)**4*x(3)*x(4)**3 - 370*x(2)**4*x(3)*x(4)**2*x(5) - 370*x(2)**4*x(3)*x(4)*x(5)**2 + 530*x(2)**4*x(3)*x(5)**3 - 400*x(2)**4*x(4)**4 - 45*x(2)**4*x(4)**3*x(5) + 720*x(2)**4*x(4)**2*x(5)**2 - 45*x(2)**4*x(4)*x(5)**3 - 400*x(2)**4*x(5)**4 - 30*x(2)**3*x(3)**4*x(4) - 30*x(2)**3*x(3)**4*x(5) + 130*x(2)**3*x(3)**3*x(4)**2 - 170*x(2)**3*x(3)**3*x(4)*x(5) + 130*x(2)**3*x(3)**3*x(5)**2 + 130*x(2)**3*x(3)**2*x(4)**3 - 10*x(2)**3*x(3)**2*x(4)**2*x(5) - 10*x(2)**3*x(3)**2*x(4)*x(5)**2 + 130*x(2)**3*x(3)**2*x(5)**3 - 2392*x(2)**3*x(3)*x(4)**4 + 5998*x(2)**3*x(3)*x(4)**3*x(5) - 7872*x(2)**3*x(3)*x(4)**2*x(5)**2 + 5998*x(2)**3*x(3)*x(4)*x(5)**3 - 2392*x(2)**3*x(3)*x(5)**4 + 2162*x(2)**3*x(4)**5 - 3941*x(2)**3*x(4)**4*x(5) + 1974*x(2)**3*x(4)**3*x(5)**2 + 1974*x(2)**3*x(4)**2*x(5)**3 - 3941*x(2)**3*x(4)*x(5)**4 + 2162*x(2)**3*x(5)**5 + 130*x(2)**2*x(3)**4*x(4)**2 - 110*x(2)**2*x(3)**4*x(4)*x(5) + 130*x(2)**2*x(3)**4*x(5)**2 - 550*x(2)**2*x(3)**3*x(4)**3 + 350*x(2)**2*x(3)**3*x(4)**2*x(5) + 350*x(2)**2*x(3)**3*x(4)*x(5)**2 - 550*x(2)**2*x(3)**3*x(5)**3 + 1972*x(2)**2*x(3)**2*x(4)**4 - 6338*x(2)**2*x(3)**2*x(4)**3*x(5) + 8832*x(2)**2*x(3)**2*x(4)**2*x(5)**2 - 6338*x(2)**2*x(3)**2*x(4)*x(5)**3 + 1972*x(2)**2*x(3)**2*x(5)**4 + 1972*x(2)**2*x(3)*x(4)**5 - 3816*x(2)**2*x(3)*x(4)**4*x(5) + 2144*x(2)**2*x(3)*x(4)**3*x(5)**2 + 2144*x(2)**2*x(3)*x(4)**2*x(5)**3 - 3816*x(2)**2*x(3)*x(4)*x(5)**4 + 1972*x(2)**2*x(3)*x(5)**5 - 3524*x(2)**2*x(4)**6 + 6695*x(2)**2*x(4)**5*x(5) - 2844*x(2)**2*x(4)**4*x(5)**2 - 1104*x(2)**2*x(4)**3*x(5)**3 - 2844*x(2)**2*x(4)**2*x(5)**4 + 6695*x(2)**2*x(4)*x(5)**5 - 3524*x(2)**2*x(5)**6 - 210*x(2)*x(3)**4*x(4)**3 + 120*x(2)*x(3)**4*x(4)**2*x(5) + 120*x(2)*x(3)**4*x(4)*x(5)**2 - 210*x(2)*x(3)**4*x(5)**3 + 720*x(2)*x(3)**3*x(4)**4 + 150*x(2)*x(3)**3*x(4)**3*x(5) - 1170*x(2)*x(3)**3*x(4)**2*x(5)**2 + 150*x(2)*x(3)**3*x(4)*x(5)**3 + 720*x(2)*x(3)**3*x(5)**4 - 3834*x(2)*x(3)**2*x(4)**5 + 8812*x(2)*x(3)**2*x(4)**4*x(5) - 5218*x(2)*x(3)**2*x(4)**3*x(5)**2 - 5218*x(2)*x(3)**2*x(4)**2*x(5)**3 + 8812*x(2)*x(3)**2*x(4)*x(5)**4 - 3834*x(2)*x(3)**2*x(5)**5 + 1662*x(2)*x(3)*x(4)**6 - 6558*x(2)*x(3)*x(4)**5*x(5) + 10804*x(2)*x(3)*x(4)**4*x(5)**2 - 11986*x(2)*x(3)*x(4)**3*x(5)**3 + 10804*x(2)*x(3)*x(4)**2*x(5)**4 - 6558*x(2)*x(3)*x(4)*x(5)**5 + 1662*x(2)*x(3)*x(5)**6 + 1662*x(2)*x(4)**7 - 962*x(2)*x(4)**6*x(5) - 4251*x(2)*x(4)**5*x(5)**2 + 3681*x(2)*x(4)**4*x(5)**3 + 3681*x(2)*x(4)**3*x(5)**4 - 4251*x(2)*x(4)**2*x(5)**5 - 962*x(2)*x(4)*x(5)**6 + 1662*x(2)*x(5)**7 + 100*x(3)**4*x(4)**4 + 20*x(3)**4*x(4)**3*x(5) - 150*x(3)**4*x(4)**2*x(5)**2 + 20*x(3)**4*x(4)*x(5)**3 + 100*x(3)**4*x(5)**4 - 300*x(3)**3*x(4)**5 - 340*x(3)**3*x(4)**4*x(5) + 490*x(3)**3*x(4)**3*x(5)**2 + 490*x(3)**3*x(4)**2*x(5)**3 - 340*x(3)**3*x(4)*x(5)**4 - 300*x(3)**3*x(5)**5 + 1862*x(3)**2*x(4)**6 - 2604*x(3)**2*x(4)**5*x(5) - 1792*x(3)**2*x(4)**4*x(5)**2 + 5378*x(3)**2*x(4)**3*x(5)**3 - 1792*x(3)**2*x(4)**2*x(5)**4 - 2604*x(3)**2*x(4)*x(5)**5 + 1862*x(3)**2*x(5)**6 - 1662*x(3)*x(4)**7 + 4586*x(3)*x(4)**6*x(5) - 4596*x(3)*x(4)**5*x(5)**2 + 1652*x(3)*x(4)**4*x(5)**3 + 1652*x(3)*x(4)**3*x(5)**4 - 4596*x(3)*x(4)**2*x(5)**5 + 4586*x(3)*x(4)*x(5)**6 - 1662*x(3)*x(5)**7 - 1662*x(4)**7*x(5) + 4486*x(4)**6*x(5)**2 - 4606*x(4)**5*x(5)**3 + 3504*x(4)**4*x(5)**4 - 4606*x(4)**3*x(5)**5 + 4486*x(4)**2*x(5)**6 - 1662*x(4)*x(5)**7))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(0,7) = -(4*(x(4) - x(5))**2*(- 10*x(1)**4*x(2)*x(3) + 5*x(1)**4*x(2)*x(4) + 5*x(1)**4*x(2)*x(5) - 10*x(1)**4*x(3)**2 + 10*x(1)**4*x(3)*x(4) + 20*x(1)**4*x(3)*x(5) - 100*x(1)**4*x(4)**2 + 185*x(1)**4*x(4)*x(5) - 105*x(1)**4*x(5)**2 - 10*x(1)**3*x(2)**2*x(3) + 5*x(1)**3*x(2)**2*x(4) + 5*x(1)**3*x(2)**2*x(5) - 10*x(1)**3*x(2)*x(3)**2 + 30*x(1)**3*x(2)*x(3)*x(4) + 50*x(1)**3*x(2)*x(3)*x(5) - 210*x(1)**3*x(2)*x(4)**2 + 360*x(1)**3*x(2)*x(4)*x(5) - 220*x(1)**3*x(2)*x(5)**2 + 10*x(1)**3*x(3)**3 + 10*x(1)**3*x(3)**2*x(4) + 10*x(1)**3*x(3)**2*x(5) - 120*x(1)**3*x(3)*x(4)**2 + 150*x(1)**3*x(3)*x(4)*x(5) - 150*x(1)**3*x(3)*x(5)**2 + 300*x(1)**3*x(4)**3 - 170*x(1)**3*x(4)**2*x(5) - 455*x(1)**3*x(4)*x(5)**2 + 415*x(1)**3*x(5)**3 - 10*x(1)**2*x(2)**3*x(3) + 5*x(1)**2*x(2)**3*x(4) + 5*x(1)**2*x(2)**3*x(5) - 10*x(1)**2*x(2)**2*x(3)**2 + 30*x(1)**2*x(2)**2*x(3)*x(4) + 50*x(1)**2*x(2)**2*x(3)*x(5) - 210*x(1)**2*x(2)**2*x(4)**2 + 360*x(1)**2*x(2)**2*x(4)*x(5) - 220*x(1)**2*x(2)**2*x(5)**2 + 20*x(1)**2*x(2)*x(3)**3 - 10*x(1)**2*x(2)*x(3)**2*x(5) - 130*x(1)**2*x(2)*x(3)*x(4)**2 + 110*x(1)**2*x(2)*x(3)*x(4)*x(5) - 170*x(1)**2*x(2)*x(3)*x(5)**2 + 505*x(1)**2*x(2)*x(4)**3 - 335*x(1)**2*x(2)*x(4)**2*x(5) - 620*x(1)**2*x(2)*x(4)*x(5)**2 + 630*x(1)**2*x(2)*x(5)**3 - 20*x(1)**2*x(3)**3*x(4) - 30*x(1)**2*x(3)**3*x(5) + 110*x(1)**2*x(3)**2*x(4)**2 - 190*x(1)**2*x(3)**2*x(4)*x(5) + 130*x(1)**2*x(3)**2*x(5)**2 + 110*x(1)**2*x(3)*x(4)**3 - 60*x(1)**2*x(3)*x(4)**2*x(5) - 30*x(1)**2*x(3)*x(4)*x(5)**2 + 130*x(1)**2*x(3)*x(5)**3 - 1862*x(1)**2*x(4)**4 + 5933*x(1)**2*x(4)**3*x(5) - 8447*x(1)**2*x(4)**2*x(5)**2 + 6303*x(1)**2*x(4)*x(5)**3 - 2077*x(1)**2*x(5)**4 - 10*x(1)*x(2)**4*x(3) + 5*x(1)*x(2)**4*x(4) + 5*x(1)*x(2)**4*x(5) - 10*x(1)*x(2)**3*x(3)**2 + 30*x(1)*x(2)**3*x(3)*x(4) + 50*x(1)*x(2)**3*x(3)*x(5) - 210*x(1)*x(2)**3*x(4)**2 + 360*x(1)*x(2)**3*x(4)*x(5) - 220*x(1)*x(2)**3*x(5)**2 + 20*x(1)*x(2)**2*x(3)**3 - 10*x(1)*x(2)**2*x(3)**2*x(5) - 130*x(1)*x(2)**2*x(3)*x(4)**2 + 110*x(1)*x(2)**2*x(3)*x(4)*x(5) - 170*x(1)*x(2)**2*x(3)*x(5)**2 + 505*x(1)*x(2)**2*x(4)**3 - 335*x(1)*x(2)**2*x(4)**2*x(5) - 620*x(1)*x(2)**2*x(4)*x(5)**2 + 630*x(1)*x(2)**2*x(5)**3 - 35*x(1)*x(2)*x(3)**3*x(4) - 45*x(1)*x(2)*x(3)**3*x(5) + 225*x(1)*x(2)*x(3)**2*x(4)**2 - 345*x(1)*x(2)*x(3)**2*x(4)*x(5) + 260*x(1)*x(2)*x(3)**2*x(5)**2 + 10*x(1)*x(2)*x(3)*x(4)**3 + 40*x(1)*x(2)*x(3)*x(4)**2*x(5) + 85*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 25*x(1)*x(2)*x(3)*x(5)**3 - 1962*x(1)*x(2)*x(4)**4 + 5818*x(1)*x(2)*x(4)**3*x(5) - 8077*x(1)*x(2)*x(4)**2*x(5)**2 + 6183*x(1)*x(2)*x(4)*x(5)**3 - 2182*x(1)*x(2)*x(5)**4 + 110*x(1)*x(3)**3*x(4)**2 - 145*x(1)*x(3)**3*x(4)*x(5) + 125*x(1)*x(3)**3*x(5)**2 - 310*x(1)*x(3)**2*x(4)**3 + 155*x(1)*x(3)**2*x(4)**2*x(5) + 425*x(1)*x(3)**2*x(4)*x(5)**2 - 440*x(1)*x(3)**2*x(5)**3 + 1662*x(1)*x(3)*x(4)**4 - 6258*x(1)*x(3)*x(4)**3*x(5) + 9272*x(1)*x(3)*x(4)**2*x(5)**2 - 6473*x(1)*x(3)*x(4)*x(5)**3 + 1767*x(1)*x(3)*x(5)**4 + 1662*x(1)*x(4)**5 - 4286*x(1)*x(4)**4*x(5) + 2859*x(1)*x(4)**3*x(5)**2 + 2374*x(1)*x(4)**2*x(5)**3 - 4266*x(1)*x(4)*x(5)**4 + 1767*x(1)*x(5)**5 - 10*x(2)**4*x(3)**2 + 10*x(2)**4*x(3)*x(4) + 20*x(2)**4*x(3)*x(5) - 100*x(2)**4*x(4)**2 + 185*x(2)**4*x(4)*x(5) - 105*x(2)**4*x(5)**2 + 10*x(2)**3*x(3)**3 + 10*x(2)**3*x(3)**2*x(4) + 10*x(2)**3*x(3)**2*x(5) - 120*x(2)**3*x(3)*x(4)**2 + 150*x(2)**3*x(3)*x(4)*x(5) - 150*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 - 170*x(2)**3*x(4)**2*x(5) - 455*x(2)**3*x(4)*x(5)**2 + 415*x(2)**3*x(5)**3 - 20*x(2)**2*x(3)**3*x(4) - 30*x(2)**2*x(3)**3*x(5) + 110*x(2)**2*x(3)**2*x(4)**2 - 190*x(2)**2*x(3)**2*x(4)*x(5) + 130*x(2)**2*x(3)**2*x(5)**2 + 110*x(2)**2*x(3)*x(4)**3 - 60*x(2)**2*x(3)*x(4)**2*x(5) - 30*x(2)**2*x(3)*x(4)*x(5)**2 + 130*x(2)**2*x(3)*x(5)**3 - 1862*x(2)**2*x(4)**4 + 5933*x(2)**2*x(4)**3*x(5) - 8447*x(2)**2*x(4)**2*x(5)**2 + 6303*x(2)**2*x(4)*x(5)**3 - 2077*x(2)**2*x(5)**4 + 110*x(2)*x(3)**3*x(4)**2 - 145*x(2)*x(3)**3*x(4)*x(5) + 125*x(2)*x(3)**3*x(5)**2 - 310*x(2)*x(3)**2*x(4)**3 + 155*x(2)*x(3)**2*x(4)**2*x(5) + 425*x(2)*x(3)**2*x(4)*x(5)**2 - 440*x(2)*x(3)**2*x(5)**3 + 1662*x(2)*x(3)*x(4)**4 - 6258*x(2)*x(3)*x(4)**3*x(5) + 9272*x(2)*x(3)*x(4)**2*x(5)**2 - 6473*x(2)*x(3)*x(4)*x(5)**3 + 1767*x(2)*x(3)*x(5)**4 + 1662*x(2)*x(4)**5 - 4286*x(2)*x(4)**4*x(5) + 2859*x(2)*x(4)**3*x(5)**2 + 2374*x(2)*x(4)**2*x(5)**3 - 4266*x(2)*x(4)*x(5)**4 + 1767*x(2)*x(5)**5 - 100*x(3)**3*x(4)**3 + 80*x(3)**3*x(4)**2*x(5) + 65*x(3)**3*x(4)*x(5)**2 - 105*x(3)**3*x(5)**3 + 200*x(3)**2*x(4)**4 + 120*x(3)**2*x(4)**3*x(5) - 455*x(3)**2*x(4)**2*x(5)**2 - 45*x(3)**2*x(4)*x(5)**3 + 310*x(3)**2*x(5)**4 - 1662*x(3)*x(4)**5 + 4586*x(3)*x(4)**4*x(5) - 3034*x(3)*x(4)**3*x(5)**2 - 2844*x(3)*x(4)**2*x(5)**3 + 4681*x(3)*x(4)*x(5)**4 - 1767*x(3)*x(5)**5 - 1662*x(4)**5*x(5) + 6148*x(4)**4*x(5)**2 - 9092*x(4)**3*x(5)**3 + 6343*x(4)**2*x(5)**4 - 1767*x(4)*x(5)**5))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(0,8) = (24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(4) - x(5))**6 + (x(4)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(4) - x(5)))/5 - (x(5)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(4) - x(5)))/5 + x(4)*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(5)*(x(4) - x(5))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(4)**3*(x(4) - x(5))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 - (x(5)**3*(x(4) - x(5))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 + (x(4)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(4) - x(5))**3)/3 - (x(5)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(4) - x(5))**3)/3 + x(4)*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(5)*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(4)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(5)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(4) - x(5))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - x(4)**2*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(5)**2*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - (x(4)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2 + (x(5)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(4) - x(5))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2

                            Cb(0,9) = (4*(x(4) - x(5))**2*(- 10*x(1)**3*x(2)**3 - 20*x(1)**3*x(2)**2*x(3) + 20*x(1)**3*x(2)**2*x(4) + 30*x(1)**3*x(2)**2*x(5) - 20*x(1)**3*x(2)*x(3)**2 + 35*x(1)**3*x(2)*x(3)*x(4) + 45*x(1)**3*x(2)*x(3)*x(5) - 110*x(1)**3*x(2)*x(4)**2 + 145*x(1)**3*x(2)*x(4)*x(5) - 125*x(1)**3*x(2)*x(5)**2 - 10*x(1)**3*x(3)**3 + 20*x(1)**3*x(3)**2*x(4) + 30*x(1)**3*x(3)**2*x(5) - 110*x(1)**3*x(3)*x(4)**2 + 145*x(1)**3*x(3)*x(4)*x(5) - 125*x(1)**3*x(3)*x(5)**2 + 100*x(1)**3*x(4)**3 - 80*x(1)**3*x(4)**2*x(5) - 65*x(1)**3*x(4)*x(5)**2 + 105*x(1)**3*x(5)**3 - 20*x(1)**2*x(2)**3*x(3) + 20*x(1)**2*x(2)**3*x(4) + 30*x(1)**2*x(2)**3*x(5) - 30*x(1)**2*x(2)**2*x(3)**2 + 75*x(1)**2*x(2)**2*x(3)*x(4) + 105*x(1)**2*x(2)**2*x(3)*x(5) - 240*x(1)**2*x(2)**2*x(4)**2 + 285*x(1)**2*x(2)**2*x(4)*x(5) - 285*x(1)**2*x(2)**2*x(5)**2 - 20*x(1)**2*x(2)*x(3)**3 + 75*x(1)**2*x(2)*x(3)**2*x(4) + 105*x(1)**2*x(2)*x(3)**2*x(5) - 480*x(1)**2*x(2)*x(3)*x(4)**2 + 555*x(1)**2*x(2)*x(3)*x(4)*x(5) - 555*x(1)**2*x(2)*x(3)*x(5)**2 + 520*x(1)**2*x(2)*x(4)**3 - 270*x(1)**2*x(2)*x(4)**2*x(5) - 510*x(1)**2*x(2)*x(4)*x(5)**2 + 670*x(1)**2*x(2)*x(5)**3 + 20*x(1)**2*x(3)**3*x(4) + 30*x(1)**2*x(3)**3*x(5) - 240*x(1)**2*x(3)**2*x(4)**2 + 285*x(1)**2*x(3)**2*x(4)*x(5) - 285*x(1)**2*x(3)**2*x(5)**2 + 520*x(1)**2*x(3)*x(4)**3 - 270*x(1)**2*x(3)*x(4)**2*x(5) - 510*x(1)**2*x(3)*x(4)*x(5)**2 + 670*x(1)**2*x(3)*x(5)**3 - 300*x(1)**2*x(4)**4 - 140*x(1)**2*x(4)**3*x(5) + 600*x(1)**2*x(4)**2*x(5)**2 + 5*x(1)**2*x(4)*x(5)**3 - 415*x(1)**2*x(5)**4 - 20*x(1)*x(2)**3*x(3)**2 + 35*x(1)*x(2)**3*x(3)*x(4) + 45*x(1)*x(2)**3*x(3)*x(5) - 110*x(1)*x(2)**3*x(4)**2 + 145*x(1)*x(2)**3*x(4)*x(5) - 125*x(1)*x(2)**3*x(5)**2 - 20*x(1)*x(2)**2*x(3)**3 + 75*x(1)*x(2)**2*x(3)**2*x(4) + 105*x(1)*x(2)**2*x(3)**2*x(5) - 480*x(1)*x(2)**2*x(3)*x(4)**2 + 555*x(1)*x(2)**2*x(3)*x(4)*x(5) - 555*x(1)*x(2)**2*x(3)*x(5)**2 + 520*x(1)*x(2)**2*x(4)**3 - 270*x(1)*x(2)**2*x(4)**2*x(5) - 510*x(1)*x(2)**2*x(4)*x(5)**2 + 670*x(1)*x(2)**2*x(5)**3 + 35*x(1)*x(2)*x(3)**3*x(4) + 45*x(1)*x(2)*x(3)**3*x(5) - 480*x(1)*x(2)*x(3)**2*x(4)**2 + 555*x(1)*x(2)*x(3)**2*x(4)*x(5) - 555*x(1)*x(2)*x(3)**2*x(5)**2 + 1155*x(1)*x(2)*x(3)*x(4)**3 - 585*x(1)*x(2)*x(3)*x(4)**2*x(5) - 1080*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 1470*x(1)*x(2)*x(3)*x(5)**3 - 2272*x(1)*x(2)*x(4)**4 + 5853*x(1)*x(2)*x(4)**3*x(5) - 7947*x(1)*x(2)*x(4)**2*x(5)**2 + 6338*x(1)*x(2)*x(4)*x(5)**3 - 2622*x(1)*x(2)*x(5)**4 - 110*x(1)*x(3)**3*x(4)**2 + 145*x(1)*x(3)**3*x(4)*x(5) - 125*x(1)*x(3)**3*x(5)**2 + 520*x(1)*x(3)**2*x(4)**3 - 270*x(1)*x(3)**2*x(4)**2*x(5) - 510*x(1)*x(3)**2*x(4)*x(5)**2 + 670*x(1)*x(3)**2*x(5)**3 - 2272*x(1)*x(3)*x(4)**4 + 5853*x(1)*x(3)*x(4)**3*x(5) - 7947*x(1)*x(3)*x(4)**2*x(5)**2 + 6338*x(1)*x(3)*x(4)*x(5)**3 - 2622*x(1)*x(3)*x(5)**4 + 1862*x(1)*x(4)**5 - 4166*x(1)*x(4)**4*x(5) + 2619*x(1)*x(4)**3*x(5)**2 + 2279*x(1)*x(4)**2*x(5)**3 - 4311*x(1)*x(4)*x(5)**4 + 2077*x(1)*x(5)**5 - 10*x(2)**3*x(3)**3 + 20*x(2)**3*x(3)**2*x(4) + 30*x(2)**3*x(3)**2*x(5) - 110*x(2)**3*x(3)*x(4)**2 + 145*x(2)**3*x(3)*x(4)*x(5) - 125*x(2)**3*x(3)*x(5)**2 + 100*x(2)**3*x(4)**3 - 80*x(2)**3*x(4)**2*x(5) - 65*x(2)**3*x(4)*x(5)**2 + 105*x(2)**3*x(5)**3 + 20*x(2)**2*x(3)**3*x(4) + 30*x(2)**2*x(3)**3*x(5) - 240*x(2)**2*x(3)**2*x(4)**2 + 285*x(2)**2*x(3)**2*x(4)*x(5) - 285*x(2)**2*x(3)**2*x(5)**2 + 520*x(2)**2*x(3)*x(4)**3 - 270*x(2)**2*x(3)*x(4)**2*x(5) - 510*x(2)**2*x(3)*x(4)*x(5)**2 + 670*x(2)**2*x(3)*x(5)**3 - 300*x(2)**2*x(4)**4 - 140*x(2)**2*x(4)**3*x(5) + 600*x(2)**2*x(4)**2*x(5)**2 + 5*x(2)**2*x(4)*x(5)**3 - 415*x(2)**2*x(5)**4 - 110*x(2)*x(3)**3*x(4)**2 + 145*x(2)*x(3)**3*x(4)*x(5) - 125*x(2)*x(3)**3*x(5)**2 + 520*x(2)*x(3)**2*x(4)**3 - 270*x(2)*x(3)**2*x(4)**2*x(5) - 510*x(2)*x(3)**2*x(4)*x(5)**2 + 670*x(2)*x(3)**2*x(5)**3 - 2272*x(2)*x(3)*x(4)**4 + 5853*x(2)*x(3)*x(4)**3*x(5) - 7947*x(2)*x(3)*x(4)**2*x(5)**2 + 6338*x(2)*x(3)*x(4)*x(5)**3 - 2622*x(2)*x(3)*x(5)**4 + 1862*x(2)*x(4)**5 - 4166*x(2)*x(4)**4*x(5) + 2619*x(2)*x(4)**3*x(5)**2 + 2279*x(2)*x(4)**2*x(5)**3 - 4311*x(2)*x(4)*x(5)**4 + 2077*x(2)*x(5)**5 + 100*x(3)**3*x(4)**3 - 80*x(3)**3*x(4)**2*x(5) - 65*x(3)**3*x(4)*x(5)**2 + 105*x(3)**3*x(5)**3 - 300*x(3)**2*x(4)**4 - 140*x(3)**2*x(4)**3*x(5) + 600*x(3)**2*x(4)**2*x(5)**2 + 5*x(3)**2*x(4)*x(5)**3 - 415*x(3)**2*x(5)**4 + 1862*x(3)*x(4)**5 - 4166*x(3)*x(4)**4*x(5) + 2619*x(3)*x(4)**3*x(5)**2 + 2279*x(3)*x(4)**2*x(5)**3 - 4311*x(3)*x(4)*x(5)**4 + 2077*x(3)*x(5)**5 - 1662*x(4)**6 + 4386*x(4)**5*x(5) - 4716*x(4)**4*x(5)**2 + 3669*x(4)**3*x(5)**3 - 4461*x(4)**2*x(5)**4 + 4371*x(4)*x(5)**5 - 1767*x(5)**6))/(5*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            Cb(0,10) = -(4*(x(4) - x(5))**2*(- 5*x(1)**2*x(2)**2 - 10*x(1)**2*x(2)*x(3) + 5*x(1)**2*x(2)*x(4) + 15*x(1)**2*x(2)*x(5) - 5*x(1)**2*x(3)**2 + 5*x(1)**2*x(3)*x(4) + 15*x(1)**2*x(3)*x(5) - 50*x(1)**2*x(4)**2 + 90*x(1)**2*x(4)*x(5) - 60*x(1)**2*x(5)**2 - 10*x(1)*x(2)**2*x(3) + 5*x(1)*x(2)**2*x(4) + 15*x(1)*x(2)**2*x(5) - 10*x(1)*x(2)*x(3)**2 + 15*x(1)*x(2)*x(3)*x(4) + 45*x(1)*x(2)*x(3)*x(5) - 105*x(1)*x(2)*x(4)**2 + 175*x(1)*x(2)*x(4)*x(5) - 140*x(1)*x(2)*x(5)**2 + 5*x(1)*x(3)**2*x(4) + 15*x(1)*x(3)**2*x(5) - 105*x(1)*x(3)*x(4)**2 + 175*x(1)*x(3)*x(4)*x(5) - 140*x(1)*x(3)*x(5)**2 + 100*x(1)*x(4)**3 + 10*x(1)*x(4)**2*x(5) - 275*x(1)*x(4)*x(5)**2 + 225*x(1)*x(5)**3 - 5*x(2)**2*x(3)**2 + 5*x(2)**2*x(3)*x(4) + 15*x(2)**2*x(3)*x(5) - 50*x(2)**2*x(4)**2 + 90*x(2)**2*x(4)*x(5) - 60*x(2)**2*x(5)**2 + 5*x(2)*x(3)**2*x(4) + 15*x(2)*x(3)**2*x(5) - 105*x(2)*x(3)*x(4)**2 + 175*x(2)*x(3)*x(4)*x(5) - 140*x(2)*x(3)*x(5)**2 + 100*x(2)*x(4)**3 + 10*x(2)*x(4)**2*x(5) - 275*x(2)*x(4)*x(5)**2 + 225*x(2)*x(5)**3 - 50*x(3)**2*x(4)**2 + 90*x(3)**2*x(4)*x(5) - 60*x(3)**2*x(5)**2 + 100*x(3)*x(4)**3 + 10*x(3)*x(4)**2*x(5) - 275*x(3)*x(4)*x(5)**2 + 225*x(3)*x(5)**3 - 831*x(4)**4 + 3024*x(4)**3*x(5) - 4551*x(4)**2*x(5)**2 + 3309*x(4)*x(5)**3 - 996*x(5)**4))/(5*(x(1) - x(5))**2*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            x = x_cb(i-2:i+2)
                            x = x - x(1)

                            Cb(1,1) = (4*(x(3) - x(4))**2*(50*x(2)**2*x(3)**2 - 95*x(2)**2*x(3)*x(4) - 5*x(2)**2*x(3)*x(5) + 50*x(2)**2*x(4)**2 - 5*x(2)**2*x(4)*x(5) + 5*x(2)**2*x(5)**2 - 100*x(2)*x(3)**3 + 95*x(2)*x(3)**2*x(4) + 105*x(2)*x(3)**2*x(5) + 95*x(2)*x(3)*x(4)**2 - 190*x(2)*x(3)*x(4)*x(5) - 5*x(2)*x(3)*x(5)**2 - 100*x(2)*x(4)**3 + 105*x(2)*x(4)**2*x(5) - 5*x(2)*x(4)*x(5)**2 + 831*x(3)**4 - 3124*x(3)**3*x(4) - 100*x(3)**3*x(5) + 4591*x(3)**2*x(4)**2 + 95*x(3)**2*x(4)*x(5) + 50*x(3)**2*x(5)**2 - 3124*x(3)*x(4)**3 + 95*x(3)*x(4)**2*x(5) - 95*x(3)*x(4)*x(5)**2 + 831*x(4)**4 - 100*x(4)**3*x(5) + 50*x(4)**2*x(5)**2))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2)

                            Cb(1,2) = (4*(x(3) - x(4))**2*(- 100*x(1)**3*x(2)*x(3)**2 + 190*x(1)**3*x(2)*x(3)*x(4) + 10*x(1)**3*x(2)*x(3)*x(5) - 100*x(1)**3*x(2)*x(4)**2 + 10*x(1)**3*x(2)*x(4)*x(5) - 10*x(1)**3*x(2)*x(5)**2 + 100*x(1)**3*x(3)**3 - 95*x(1)**3*x(3)**2*x(4) - 105*x(1)**3*x(3)**2*x(5) - 95*x(1)**3*x(3)*x(4)**2 + 190*x(1)**3*x(3)*x(4)*x(5) + 5*x(1)**3*x(3)*x(5)**2 + 100*x(1)**3*x(4)**3 - 105*x(1)**3*x(4)**2*x(5) + 5*x(1)**3*x(4)*x(5)**2 - 100*x(1)**2*x(2)**2*x(3)**2 + 190*x(1)**2*x(2)**2*x(3)*x(4) + 10*x(1)**2*x(2)**2*x(3)*x(5) - 100*x(1)**2*x(2)**2*x(4)**2 + 10*x(1)**2*x(2)**2*x(4)*x(5) - 10*x(1)**2*x(2)**2*x(5)**2 + 300*x(1)**2*x(2)*x(3)**3 - 280*x(1)**2*x(2)*x(3)**2*x(4) - 120*x(1)**2*x(2)*x(3)**2*x(5) - 280*x(1)**2*x(2)*x(3)*x(4)**2 + 170*x(1)**2*x(2)*x(3)*x(4)*x(5) + 10*x(1)**2*x(2)*x(3)*x(5)**2 + 300*x(1)**2*x(2)*x(4)**3 - 120*x(1)**2*x(2)*x(4)**2*x(5) + 10*x(1)**2*x(2)*x(4)*x(5)**2 + 10*x(1)**2*x(2)*x(5)**3 - 1762*x(1)**2*x(3)**4 + 6243*x(1)**2*x(3)**3*x(4) + 205*x(1)**2*x(3)**3*x(5) - 8992*x(1)**2*x(3)**2*x(4)**2 - 180*x(1)**2*x(3)**2*x(4)*x(5) + 6243*x(1)**2*x(3)*x(4)**3 - 180*x(1)**2*x(3)*x(4)**2*x(5) - 10*x(1)**2*x(3)*x(4)*x(5)**2 - 5*x(1)**2*x(3)*x(5)**3 - 1762*x(1)**2*x(4)**4 + 205*x(1)**2*x(4)**3*x(5) - 5*x(1)**2*x(4)*x(5)**3 - 100*x(1)*x(2)**3*x(3)**2 + 190*x(1)*x(2)**3*x(3)*x(4) + 10*x(1)*x(2)**3*x(3)*x(5) - 100*x(1)*x(2)**3*x(4)**2 + 10*x(1)*x(2)**3*x(4)*x(5) - 10*x(1)*x(2)**3*x(5)**2 + 300*x(1)*x(2)**2*x(3)**3 - 280*x(1)*x(2)**2*x(3)**2*x(4) - 120*x(1)*x(2)**2*x(3)**2*x(5) - 280*x(1)*x(2)**2*x(3)*x(4)**2 + 170*x(1)*x(2)**2*x(3)*x(4)*x(5) + 10*x(1)*x(2)**2*x(3)*x(5)**2 + 300*x(1)*x(2)**2*x(4)**3 - 120*x(1)*x(2)**2*x(4)**2*x(5) + 10*x(1)*x(2)**2*x(4)*x(5)**2 + 10*x(1)*x(2)**2*x(5)**3 - 1862*x(1)*x(2)*x(3)**4 + 6138*x(1)*x(2)*x(3)**3*x(4) + 110*x(1)*x(2)*x(3)**3*x(5) - 8612*x(1)*x(2)*x(3)**2*x(4)**2 - 70*x(1)*x(2)*x(3)**2*x(4)*x(5) + 110*x(1)*x(2)*x(3)**2*x(5)**2 + 6138*x(1)*x(2)*x(3)*x(4)**3 - 70*x(1)*x(2)*x(3)*x(4)**2*x(5) - 200*x(1)*x(2)*x(3)*x(4)*x(5)**2 - 20*x(1)*x(2)*x(3)*x(5)**3 - 1862*x(1)*x(2)*x(4)**4 + 110*x(1)*x(2)*x(4)**3*x(5) + 110*x(1)*x(2)*x(4)**2*x(5)**2 - 20*x(1)*x(2)*x(4)*x(5)**3 + 1662*x(1)*x(3)**5 - 4486*x(1)*x(3)**4*x(4) + 1562*x(1)*x(3)**4*x(5) + 2839*x(1)*x(3)**3*x(4)**2 - 6358*x(1)*x(3)**3*x(4)*x(5) - 205*x(1)*x(3)**3*x(5)**2 + 2839*x(1)*x(3)**2*x(4)**3 + 9562*x(1)*x(3)**2*x(4)**2*x(5) + 190*x(1)*x(3)**2*x(4)*x(5)**2 + 105*x(1)*x(3)**2*x(5)**3 - 4486*x(1)*x(3)*x(4)**4 - 6358*x(1)*x(3)*x(4)**3*x(5) + 190*x(1)*x(3)*x(4)**2*x(5)**2 - 180*x(1)*x(3)*x(4)*x(5)**3 + 1662*x(1)*x(4)**5 + 1562*x(1)*x(4)**4*x(5) - 205*x(1)*x(4)**3*x(5)**2 + 105*x(1)*x(4)**2*x(5)**3 - 100*x(2)**4*x(3)**2 + 190*x(2)**4*x(3)*x(4) + 10*x(2)**4*x(3)*x(5) - 100*x(2)**4*x(4)**2 + 10*x(2)**4*x(4)*x(5) - 10*x(2)**4*x(5)**2 + 300*x(2)**3*x(3)**3 - 280*x(2)**3*x(3)**2*x(4) - 120*x(2)**3*x(3)**2*x(5) - 280*x(2)**3*x(3)*x(4)**2 + 170*x(2)**3*x(3)*x(4)*x(5) + 10*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 - 120*x(2)**3*x(4)**2*x(5) + 10*x(2)**3*x(4)*x(5)**2 + 10*x(2)**3*x(5)**3 - 1862*x(2)**2*x(3)**4 + 6138*x(2)**2*x(3)**3*x(4) + 110*x(2)**2*x(3)**3*x(5) - 8612*x(2)**2*x(3)**2*x(4)**2 - 70*x(2)**2*x(3)**2*x(4)*x(5) + 110*x(2)**2*x(3)**2*x(5)**2 + 6138*x(2)**2*x(3)*x(4)**3 - 70*x(2)**2*x(3)*x(4)**2*x(5) - 200*x(2)**2*x(3)*x(4)*x(5)**2 - 20*x(2)**2*x(3)*x(5)**3 - 1862*x(2)**2*x(4)**4 + 110*x(2)**2*x(4)**3*x(5) + 110*x(2)**2*x(4)**2*x(5)**2 - 20*x(2)**2*x(4)*x(5)**3 + 1662*x(2)*x(3)**5 - 4386*x(2)*x(3)**4*x(4) + 1662*x(2)*x(3)**4*x(5) + 2744*x(2)*x(3)**3*x(4)**2 - 6358*x(2)*x(3)**3*x(4)*x(5) - 310*x(2)*x(3)**3*x(5)**2 + 2744*x(2)*x(3)**2*x(4)**3 + 9372*x(2)*x(3)**2*x(4)**2*x(5) + 270*x(2)*x(3)**2*x(4)*x(5)**2 + 110*x(2)*x(3)**2*x(5)**3 - 4386*x(2)*x(3)*x(4)**4 - 6358*x(2)*x(3)*x(4)**3*x(5) + 270*x(2)*x(3)*x(4)**2*x(5)**2 - 160*x(2)*x(3)*x(4)*x(5)**3 + 1662*x(2)*x(4)**5 + 1662*x(2)*x(4)**4*x(5) - 310*x(2)*x(4)**3*x(5)**2 + 110*x(2)*x(4)**2*x(5)**3 - 1662*x(3)**5*x(4) - 1662*x(3)**5*x(5) + 6248*x(3)**4*x(4)**2 + 4686*x(3)**4*x(4)*x(5) + 200*x(3)**4*x(5)**2 - 9182*x(3)**3*x(4)**3 - 3029*x(3)**3*x(4)**2*x(5) + 15*x(3)**3*x(4)*x(5)**2 - 100*x(3)**3*x(5)**3 + 6248*x(3)**2*x(4)**4 - 3029*x(3)**2*x(4)**3*x(5) - 380*x(3)**2*x(4)**2*x(5)**2 + 85*x(3)**2*x(4)*x(5)**3 - 1662*x(3)*x(4)**5 + 4686*x(3)*x(4)**4*x(5) + 15*x(3)*x(4)**3*x(5)**2 + 85*x(3)*x(4)**2*x(5)**3 - 1662*x(4)**5*x(5) + 200*x(4)**4*x(5)**2 - 100*x(4)**3*x(5)**3))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(1,3) = (4*(x(3) - x(4))**2*(- 5*x(1)**2*x(2)**3*x(3) - 5*x(1)**2*x(2)**3*x(4) + 10*x(1)**2*x(2)**3*x(5) + 105*x(1)**2*x(2)**2*x(3)**2 - 185*x(1)**2*x(2)**2*x(3)*x(4) - 10*x(1)**2*x(2)**2*x(3)*x(5) + 110*x(1)**2*x(2)**2*x(4)**2 - 20*x(1)**2*x(2)**2*x(4)*x(5) - 200*x(1)**2*x(2)*x(3)**2*x(4) - 10*x(1)**2*x(2)*x(3)**2*x(5) + 385*x(1)**2*x(2)*x(3)*x(4)**2 + 20*x(1)**2*x(2)*x(3)*x(5)**2 - 205*x(1)**2*x(2)*x(4)**3 + 10*x(1)**2*x(2)*x(4)**2*x(5) + 10*x(1)**2*x(2)*x(4)*x(5)**2 - 10*x(1)**2*x(2)*x(5)**3 - 100*x(1)**2*x(3)**4 + 195*x(1)**2*x(3)**3*x(4) + 205*x(1)**2*x(3)**3*x(5) - 385*x(1)**2*x(3)**2*x(4)*x(5) - 110*x(1)**2*x(3)**2*x(5)**2 - 195*x(1)**2*x(3)*x(4)**3 + 200*x(1)**2*x(3)*x(4)**2*x(5) + 185*x(1)**2*x(3)*x(4)*x(5)**2 + 5*x(1)**2*x(3)*x(5)**3 + 100*x(1)**2*x(4)**4 - 105*x(1)**2*x(4)**2*x(5)**2 + 5*x(1)**2*x(4)*x(5)**3 + 100*x(1)*x(2)**3*x(3)**2 - 185*x(1)*x(2)**3*x(3)*x(4) - 5*x(1)*x(2)**3*x(3)*x(5) + 105*x(1)*x(2)**3*x(4)**2 - 15*x(1)*x(2)**3*x(4)*x(5) - 100*x(1)*x(2)**2*x(3)**3 - 205*x(1)*x(2)**2*x(3)**2*x(4) - 5*x(1)*x(2)**2*x(3)**2*x(5) + 665*x(1)*x(2)**2*x(3)*x(4)**2 + 5*x(1)*x(2)**2*x(3)*x(4)*x(5) + 20*x(1)*x(2)**2*x(3)*x(5)**2 - 410*x(1)*x(2)**2*x(4)**3 + 30*x(1)*x(2)**2*x(4)**2*x(5) + 10*x(1)*x(2)**2*x(4)*x(5)**2 - 10*x(1)*x(2)**2*x(5)**3 + 1462*x(1)*x(2)*x(3)**4 - 5858*x(1)*x(2)*x(3)**3*x(4) + 210*x(1)*x(2)*x(3)**3*x(5) + 9382*x(1)*x(2)*x(3)**2*x(4)**2 - 380*x(1)*x(2)*x(3)**2*x(4)*x(5) - 110*x(1)*x(2)*x(3)**2*x(5)**2 - 7023*x(1)*x(2)*x(3)*x(4)**3 + 205*x(1)*x(2)*x(3)*x(4)**2*x(5) + 170*x(1)*x(2)*x(3)*x(4)*x(5)**2 - 10*x(1)*x(2)*x(3)*x(5)**3 + 2067*x(1)*x(2)*x(4)**4 - 15*x(1)*x(2)*x(4)**3*x(5) - 120*x(1)*x(2)*x(4)**2*x(5)**2 + 10*x(1)*x(2)*x(4)*x(5)**3 + 10*x(1)*x(2)*x(5)**4 + 1662*x(1)*x(3)**5 - 7810*x(1)*x(3)**4*x(4) - 1762*x(1)*x(3)**4*x(5) + 15235*x(1)*x(3)**3*x(4)**2 + 6238*x(1)*x(3)**3*x(4)*x(5) + 95*x(1)*x(3)**3*x(5)**2 - 15430*x(1)*x(3)**2*x(4)**3 - 8797*x(1)*x(3)**2*x(4)**2*x(5) + 15*x(1)*x(3)**2*x(4)*x(5)**2 + 10*x(1)*x(3)**2*x(5)**3 + 8105*x(1)*x(3)*x(4)**4 + 5853*x(1)*x(3)*x(4)**3*x(5) - 285*x(1)*x(3)*x(4)**2*x(5)**2 - 5*x(1)*x(3)*x(5)**4 - 1762*x(1)*x(4)**5 - 1562*x(1)*x(4)**4*x(5) + 205*x(1)*x(4)**3*x(5)**2 - 5*x(1)*x(4)*x(5)**4 + 100*x(2)**3*x(3)**3 - 290*x(2)**3*x(3)**2*x(4) - 110*x(2)**3*x(3)**2*x(5) + 290*x(2)**3*x(3)*x(4)**2 + 185*x(2)**3*x(3)*x(4)*x(5) + 20*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 - 95*x(2)**3*x(4)**2*x(5) + 10*x(2)**3*x(4)*x(5)**2 - 10*x(2)**3*x(5)**3 - 200*x(2)**2*x(3)**4 + 290*x(2)**2*x(3)**3*x(4) + 310*x(2)**2*x(3)**3*x(5) + 290*x(2)**2*x(3)**2*x(4)**2 - 375*x(2)**2*x(3)**2*x(4)*x(5) - 110*x(2)**2*x(3)**2*x(5)**2 - 680*x(2)**2*x(3)*x(4)**3 - 75*x(2)**2*x(3)*x(4)**2*x(5) + 170*x(2)**2*x(3)*x(4)*x(5)**2 - 10*x(2)**2*x(3)*x(5)**3 + 300*x(2)**2*x(4)**4 + 190*x(2)**2*x(4)**3*x(5) - 120*x(2)**2*x(4)**2*x(5)**2 + 10*x(2)**2*x(4)*x(5)**3 + 10*x(2)**2*x(5)**4 + 1662*x(2)*x(3)**5 - 7710*x(2)*x(3)**4*x(4) - 1662*x(2)*x(3)**4*x(5) + 15040*x(2)*x(3)**3*x(4)**2 + 6038*x(2)*x(3)**3*x(4)*x(5) - 110*x(2)*x(3)**3*x(5)**2 - 15430*x(2)*x(3)**2*x(4)**3 - 8792*x(2)*x(3)**2*x(4)**2*x(5) + 300*x(2)*x(3)**2*x(4)*x(5)**2 + 120*x(2)*x(3)**2*x(5)**3 + 8300*x(2)*x(3)*x(4)**4 + 6043*x(2)*x(3)*x(4)**3*x(5) - 300*x(2)*x(3)*x(4)**2*x(5)**2 - 170*x(2)*x(3)*x(4)*x(5)**3 - 10*x(2)*x(3)*x(5)**4 - 1862*x(2)*x(4)**5 - 1657*x(2)*x(4)**4*x(5) + 110*x(2)*x(4)**3*x(5)**2 + 110*x(2)*x(4)**2*x(5)**3 - 20*x(2)*x(4)*x(5)**4 - 1662*x(3)**5*x(4) - 1662*x(3)**5*x(5) + 7910*x(3)**4*x(4)**2 + 8010*x(3)**4*x(4)*x(5) + 1862*x(3)**4*x(5)**2 - 15430*x(3)**3*x(4)**3 - 15625*x(3)**3*x(4)**2*x(5) - 6533*x(3)**3*x(4)*x(5)**2 - 300*x(3)**3*x(5)**3 + 15430*x(3)**2*x(4)**4 + 15430*x(3)**2*x(4)**3*x(5) + 9087*x(3)**2*x(4)**2*x(5)**2 + 370*x(3)**2*x(4)*x(5)**3 + 100*x(3)**2*x(5)**4 - 7910*x(3)*x(4)**5 - 7715*x(3)*x(4)**4*x(5) - 5948*x(3)*x(4)**3*x(5)**2 + 85*x(3)*x(4)**2*x(5)**3 - 185*x(3)*x(4)*x(5)**4 + 1662*x(4)**6 + 1562*x(4)**5*x(5) + 1562*x(4)**4*x(5)**2 - 205*x(4)**3*x(5)**3 + 105*x(4)**2*x(5)**4))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(1,4) = -(4*(x(3) - x(4))**2*(100*x(2)**2*x(3)**2 + 100*x(2)**2*x(4)**2 + 9182*x(3)**2*x(4)**2 - 100*x(1)*x(3)**3 - 100*x(1)*x(4)**3 - 200*x(2)*x(3)**3 - 200*x(2)*x(4)**3 - 6248*x(3)*x(4)**3 - 6248*x(3)**3*x(4) - 100*x(3)**3*x(5) - 100*x(4)**3*x(5) + 1662*x(3)**4 + 1662*x(4)**4 + 105*x(1)*x(2)*x(3)**2 - 5*x(1)*x(2)**2*x(3) + 105*x(1)*x(2)*x(4)**2 - 5*x(1)*x(2)**2*x(4) + 95*x(1)*x(3)*x(4)**2 + 10*x(1)*x(2)**2*x(5) + 95*x(1)*x(3)**2*x(4) + 100*x(1)*x(3)**2*x(5) + 190*x(2)*x(3)*x(4)**2 + 190*x(2)*x(3)**2*x(4) - 190*x(2)**2*x(3)*x(4) + 100*x(1)*x(4)**2*x(5) + 105*x(2)*x(3)**2*x(5) - 5*x(2)**2*x(3)*x(5) + 105*x(2)*x(4)**2*x(5) - 5*x(2)**2*x(4)*x(5) + 95*x(3)*x(4)**2*x(5) + 95*x(3)**2*x(4)*x(5) - 190*x(1)*x(2)*x(3)*x(4) - 10*x(1)*x(2)*x(3)*x(5) - 10*x(1)*x(2)*x(4)*x(5) - 190*x(1)*x(3)*x(4)*x(5) - 190*x(2)*x(3)*x(4)*x(5)))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(1,5) = x(3)*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(4)*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + x(3)*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(4)*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(3)**3*(((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (24*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))))*(x(3) - x(4)))/3 - (x(4)**3*(((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (24*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))))*(x(3) - x(4)))/3 - x(3)**2*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(4)**2*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(4) + 2*x(1)*x(5) + 2*x(2)*x(4) + 2*x(2)*x(5) + 2*x(4)*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + (576*(x(3) - x(4))**6*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) + (192*x(3)**3*(x(3) - x(4))**3*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (192*x(4)**3*(x(3) - x(4))**3*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/((x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) + (144*x(3)**5*(x(3) - x(4))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (144*x(4)**5*(x(3) - x(4))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2)**2)/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2) - (6*x(3)**4*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) + (6*x(4)**4*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) - (24*x(3)**2*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5))) + (24*x(4)**2*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(4) + 6*x(5))/((x(1) - x(3))*(x(3) - x(4))*(x(3) - x(5))) - ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(5)))/((x(1) - x(4))*(x(2) - x(4))*(x(3) - x(4))*(x(4) - x(5))) + ((x(2) - x(3))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1)*x(2) - x(1)*x(3) - x(1)*x(4) - x(2)*x(3) - x(1)*x(5) - x(2)*x(4) - x(2)*x(5) + x(3)*x(4) + x(3)*x(5) + x(4)*x(5) + x(1)**2 + x(2)**2))/((x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(1,6) = (4*(x(3) - x(4))**2*(5*x(1)**5*x(2)**2*x(3) + 5*x(1)**5*x(2)**2*x(4) - 10*x(1)**5*x(2)**2*x(5) - 100*x(1)**5*x(2)*x(3)**2 + 185*x(1)**5*x(2)*x(3)*x(4) + 5*x(1)**5*x(2)*x(3)*x(5) - 105*x(1)**5*x(2)*x(4)**2 + 15*x(1)**5*x(2)*x(4)*x(5) - 100*x(1)**5*x(3)**3 + 290*x(1)**5*x(3)**2*x(4) + 110*x(1)**5*x(3)**2*x(5) - 290*x(1)**5*x(3)*x(4)**2 - 185*x(1)**5*x(3)*x(4)*x(5) - 20*x(1)**5*x(3)*x(5)**2 + 100*x(1)**5*x(4)**3 + 95*x(1)**5*x(4)**2*x(5) - 10*x(1)**5*x(4)*x(5)**2 + 10*x(1)**5*x(5)**3 + 5*x(1)**4*x(2)**3*x(3) + 5*x(1)**4*x(2)**3*x(4) - 10*x(1)**4*x(2)**3*x(5) - 210*x(1)**4*x(2)**2*x(3)**2 + 360*x(1)**4*x(2)**2*x(3)*x(4) + 20*x(1)**4*x(2)**2*x(3)*x(5) - 220*x(1)**4*x(2)**2*x(4)**2 + 40*x(1)**4*x(2)**2*x(4)*x(5) + 10*x(1)**4*x(2)**2*x(5)**2 + 100*x(1)**4*x(2)*x(3)**3 + 410*x(1)**4*x(2)*x(3)**2*x(4) + 210*x(1)**4*x(2)*x(3)**2*x(5) - 1035*x(1)**4*x(2)*x(3)*x(4)**2 - 395*x(1)**4*x(2)*x(3)*x(4)*x(5) - 45*x(1)**4*x(2)*x(3)*x(5)**2 + 615*x(1)**4*x(2)*x(4)**3 + 155*x(1)**4*x(2)*x(4)**2*x(5) - 35*x(1)**4*x(2)*x(4)*x(5)**2 + 20*x(1)**4*x(2)*x(5)**3 + 300*x(1)**4*x(3)**4 - 480*x(1)**4*x(3)**3*x(4) - 320*x(1)**4*x(3)**3*x(5) - 290*x(1)**4*x(3)**2*x(4)**2 + 160*x(1)**4*x(3)**2*x(4)*x(5) + 20*x(1)**4*x(3)**2*x(5)**2 + 870*x(1)**4*x(3)*x(4)**3 + 455*x(1)**4*x(3)*x(4)**2*x(5) + 45*x(1)**4*x(3)*x(4)*x(5)**2 + 20*x(1)**4*x(3)*x(5)**3 - 400*x(1)**4*x(4)**4 - 385*x(1)**4*x(4)**3*x(5) + 35*x(1)**4*x(4)**2*x(5)**2 - 10*x(1)**4*x(4)*x(5)**3 - 20*x(1)**4*x(5)**4 + 5*x(1)**3*x(2)**4*x(3) + 5*x(1)**3*x(2)**4*x(4) - 10*x(1)**3*x(2)**4*x(5) - 210*x(1)**3*x(2)**3*x(3)**2 + 360*x(1)**3*x(2)**3*x(3)*x(4) + 20*x(1)**3*x(2)**3*x(3)*x(5) - 220*x(1)**3*x(2)**3*x(4)**2 + 40*x(1)**3*x(2)**3*x(4)*x(5) + 10*x(1)**3*x(2)**3*x(5)**2 + 205*x(1)**3*x(2)**2*x(3)**3 + 535*x(1)**3*x(2)**2*x(3)**2*x(4) + 320*x(1)**3*x(2)**2*x(3)**2*x(5) - 1490*x(1)**3*x(2)**2*x(3)*x(4)**2 - 610*x(1)**3*x(2)**2*x(3)*x(4)*x(5) - 85*x(1)**3*x(2)**2*x(3)*x(5)**2 + 930*x(1)**3*x(2)**2*x(4)**3 + 240*x(1)**3*x(2)**2*x(4)**2*x(5) - 75*x(1)**3*x(2)**2*x(4)*x(5)**2 + 30*x(1)**3*x(2)**2*x(5)**3 - 1262*x(1)**3*x(2)*x(3)**4 + 4973*x(1)**3*x(2)*x(3)**3*x(4) - 735*x(1)**3*x(2)*x(3)**3*x(5) - 9027*x(1)**3*x(2)*x(3)**2*x(4)**2 + 425*x(1)**3*x(2)*x(3)**2*x(4)*x(5) + 150*x(1)**3*x(2)*x(3)**2*x(5)**2 + 8043*x(1)**3*x(2)*x(3)*x(4)**3 + 1045*x(1)**3*x(2)*x(3)*x(4)**2*x(5) - 70*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 40*x(1)**3*x(2)*x(3)*x(5)**3 - 2877*x(1)**3*x(2)*x(4)**4 - 855*x(1)**3*x(2)*x(4)**3*x(5) + 210*x(1)**3*x(2)*x(4)**2*x(5)**2 - 20*x(1)**3*x(2)*x(4)*x(5)**3 - 40*x(1)**3*x(2)*x(5)**4 - 1862*x(1)**3*x(3)**5 + 7700*x(1)**3*x(3)**4*x(4) + 1672*x(1)**3*x(3)**4*x(5) - 14170*x(1)**3*x(3)**3*x(4)**2 - 5513*x(1)**3*x(3)**3*x(4)*x(5) + 420*x(1)**3*x(3)**3*x(5)**2 + 14750*x(1)**3*x(3)**2*x(4)**3 + 8447*x(1)**3*x(3)**2*x(4)**2*x(5) - 710*x(1)**3*x(3)**2*x(4)*x(5)**2 - 260*x(1)**3*x(3)**2*x(5)**3 - 8580*x(1)**3*x(3)*x(4)**4 - 6703*x(1)**3*x(3)*x(4)**3*x(5) + 175*x(1)**3*x(3)*x(4)**2*x(5)**2 + 320*x(1)**3*x(3)*x(4)*x(5)**3 + 20*x(1)**3*x(3)*x(5)**4 + 2162*x(1)**3*x(4)**5 + 2247*x(1)**3*x(4)**4*x(5) + 55*x(1)**3*x(4)**3*x(5)**2 - 230*x(1)**3*x(4)**2*x(5)**3 + 50*x(1)**3*x(4)*x(5)**4 + 10*x(1)**3*x(5)**5 + 5*x(1)**2*x(2)**5*x(3) + 5*x(1)**2*x(2)**5*x(4) - 10*x(1)**2*x(2)**5*x(5) - 210*x(1)**2*x(2)**4*x(3)**2 + 360*x(1)**2*x(2)**4*x(3)*x(4) + 20*x(1)**2*x(2)**4*x(3)*x(5) - 220*x(1)**2*x(2)**4*x(4)**2 + 40*x(1)**2*x(2)**4*x(4)*x(5) + 10*x(1)**2*x(2)**4*x(5)**2 + 205*x(1)**2*x(2)**3*x(3)**3 + 535*x(1)**2*x(2)**3*x(3)**2*x(4) + 320*x(1)**2*x(2)**3*x(3)**2*x(5) - 1490*x(1)**2*x(2)**3*x(3)*x(4)**2 - 610*x(1)**2*x(2)**3*x(3)*x(4)*x(5) - 85*x(1)**2*x(2)**3*x(3)*x(5)**2 + 930*x(1)**2*x(2)**3*x(4)**3 + 240*x(1)**2*x(2)**3*x(4)**2*x(5) - 75*x(1)**2*x(2)**3*x(4)*x(5)**2 + 30*x(1)**2*x(2)**3*x(5)**3 - 1162*x(1)**2*x(2)**2*x(3)**4 + 4473*x(1)**2*x(2)**2*x(3)**3*x(4) - 1055*x(1)**2*x(2)**2*x(3)**3*x(5) - 8657*x(1)**2*x(2)**2*x(3)**2*x(4)**2 + 685*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) + 280*x(1)**2*x(2)**2*x(3)**2*x(5)**2 + 8308*x(1)**2*x(2)**2*x(3)*x(4)**3 + 1330*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) - 185*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 + 45*x(1)**2*x(2)**2*x(3)*x(5)**3 - 3182*x(1)**2*x(2)**2*x(4)**4 - 1160*x(1)**2*x(2)**2*x(4)**3*x(5) + 355*x(1)**2*x(2)**2*x(4)**2*x(5)**2 - 25*x(1)**2*x(2)**2*x(4)*x(5)**3 - 50*x(1)**2*x(2)**2*x(5)**4 - 2062*x(1)**2*x(2)*x(3)**5 + 11214*x(1)**2*x(2)*x(3)**4*x(4) + 5206*x(1)**2*x(2)*x(3)**4*x(5) - 25676*x(1)**2*x(2)*x(3)**3*x(4)**2 - 17369*x(1)**2*x(2)*x(3)**3*x(4)*x(5) + 430*x(1)**2*x(2)*x(3)**3*x(5)**2 + 31874*x(1)**2*x(2)*x(3)**2*x(4)**3 + 25801*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) - 1070*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 - 410*x(1)**2*x(2)*x(3)**2*x(5)**3 - 21231*x(1)**2*x(2)*x(3)*x(4)**4 - 19569*x(1)**2*x(2)*x(3)*x(4)**3*x(5) + 655*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 + 470*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 + 35*x(1)**2*x(2)*x(3)*x(5)**4 + 5991*x(1)**2*x(2)*x(4)**5 + 6271*x(1)**2*x(2)*x(4)**4*x(5) - 315*x(1)**2*x(2)*x(4)**3*x(5)**2 - 350*x(1)**2*x(2)*x(4)**2*x(5)**3 + 85*x(1)**2*x(2)*x(4)*x(5)**4 + 20*x(1)**2*x(2)*x(5)**5 + 1662*x(1)**2*x(3)**6 - 4186*x(1)**2*x(3)**5*x(4) + 1862*x(1)**2*x(3)**5*x(5) - 870*x(1)**2*x(3)**4*x(4)**2 - 11644*x(1)**2*x(3)**4*x(4)*x(5) - 3944*x(1)**2*x(3)**4*x(5)**2 + 14750*x(1)**2*x(3)**3*x(4)**3 + 27416*x(1)**2*x(3)**3*x(4)**2*x(5) + 12826*x(1)**2*x(3)**3*x(4)*x(5)**2 + 420*x(1)**2*x(3)**3*x(5)**3 - 21880*x(1)**2*x(3)**2*x(4)**4 - 32854*x(1)**2*x(3)**2*x(4)**3*x(5) - 17414*x(1)**2*x(3)**2*x(4)**2*x(5)**2 - 150*x(1)**2*x(3)**2*x(4)*x(5)**3 + 20*x(1)**2*x(3)**2*x(5)**4 + 14048*x(1)**2*x(3)*x(4)**5 + 20491*x(1)**2*x(3)*x(4)**4*x(5) + 11711*x(1)**2*x(3)*x(4)**3*x(5)**2 - 495*x(1)**2*x(3)*x(4)**2*x(5)**3 - 35*x(1)**2*x(3)*x(4)*x(5)**4 - 20*x(1)**2*x(3)*x(5)**5 - 3524*x(1)**2*x(4)**6 - 5381*x(1)**2*x(4)**5*x(5) - 3299*x(1)**2*x(4)**4*x(5)**2 + 545*x(1)**2*x(4)**3*x(5)**3 - 25*x(1)**2*x(4)**2*x(5)**4 - 30*x(1)**2*x(4)*x(5)**5 - 100*x(1)*x(2)**5*x(3)**2 + 185*x(1)*x(2)**5*x(3)*x(4) + 5*x(1)*x(2)**5*x(3)*x(5) - 105*x(1)*x(2)**5*x(4)**2 + 15*x(1)*x(2)**5*x(4)*x(5) + 100*x(1)*x(2)**4*x(3)**3 + 410*x(1)*x(2)**4*x(3)**2*x(4) + 210*x(1)*x(2)**4*x(3)**2*x(5) - 1035*x(1)*x(2)**4*x(3)*x(4)**2 - 395*x(1)*x(2)**4*x(3)*x(4)*x(5) - 45*x(1)*x(2)**4*x(3)*x(5)**2 + 615*x(1)*x(2)**4*x(4)**3 + 155*x(1)*x(2)**4*x(4)**2*x(5) - 35*x(1)*x(2)**4*x(4)*x(5)**2 + 20*x(1)*x(2)**4*x(5)**3 - 1262*x(1)*x(2)**3*x(3)**4 + 4973*x(1)*x(2)**3*x(3)**3*x(4) - 735*x(1)*x(2)**3*x(3)**3*x(5) - 9027*x(1)*x(2)**3*x(3)**2*x(4)**2 + 425*x(1)*x(2)**3*x(3)**2*x(4)*x(5) + 150*x(1)*x(2)**3*x(3)**2*x(5)**2 + 8043*x(1)*x(2)**3*x(3)*x(4)**3 + 1045*x(1)*x(2)**3*x(3)*x(4)**2*x(5) - 70*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 + 40*x(1)*x(2)**3*x(3)*x(5)**3 - 2877*x(1)*x(2)**3*x(4)**4 - 855*x(1)*x(2)**3*x(4)**3*x(5) + 210*x(1)*x(2)**3*x(4)**2*x(5)**2 - 20*x(1)*x(2)**3*x(4)*x(5)**3 - 40*x(1)*x(2)**3*x(5)**4 - 2062*x(1)*x(2)**2*x(3)**5 + 11214*x(1)*x(2)**2*x(3)**4*x(4) + 5206*x(1)*x(2)**2*x(3)**4*x(5) - 25676*x(1)*x(2)**2*x(3)**3*x(4)**2 - 17369*x(1)*x(2)**2*x(3)**3*x(4)*x(5) + 430*x(1)*x(2)**2*x(3)**3*x(5)**2 + 31874*x(1)*x(2)**2*x(3)**2*x(4)**3 + 25801*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) - 1070*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 - 410*x(1)*x(2)**2*x(3)**2*x(5)**3 - 21231*x(1)*x(2)**2*x(3)*x(4)**4 - 19569*x(1)*x(2)**2*x(3)*x(4)**3*x(5) + 655*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 + 470*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 + 35*x(1)*x(2)**2*x(3)*x(5)**4 + 5991*x(1)*x(2)**2*x(4)**5 + 6271*x(1)*x(2)**2*x(4)**4*x(5) - 315*x(1)*x(2)**2*x(4)**3*x(5)**2 - 350*x(1)*x(2)**2*x(4)**2*x(5)**3 + 85*x(1)*x(2)**2*x(4)*x(5)**4 + 20*x(1)*x(2)**2*x(5)**5 + 3324*x(1)*x(2)*x(3)**6 - 11896*x(1)*x(2)*x(3)**5*x(4) + 200*x(1)*x(2)*x(3)**5*x(5) + 12608*x(1)*x(2)*x(3)**4*x(4)**2 - 10592*x(1)*x(2)*x(3)**4*x(4)*x(5) - 5616*x(1)*x(2)*x(3)**4*x(5)**2 + 5173*x(1)*x(2)*x(3)**3*x(4)**3 + 36753*x(1)*x(2)*x(3)**3*x(4)**2*x(5) + 19169*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 + 525*x(1)*x(2)*x(3)**3*x(5)**3 - 22377*x(1)*x(2)*x(3)**2*x(4)**4 - 53507*x(1)*x(2)*x(3)**2*x(4)**3*x(5) - 26471*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 - 95*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 + 40*x(1)*x(2)*x(3)**2*x(5)**4 + 18424*x(1)*x(2)*x(3)*x(4)**5 + 37558*x(1)*x(2)*x(3)*x(4)**4*x(5) + 17529*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 - 755*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 - 45*x(1)*x(2)*x(3)*x(4)*x(5)**4 - 35*x(1)*x(2)*x(3)*x(5)**5 - 5286*x(1)*x(2)*x(4)**6 - 10672*x(1)*x(2)*x(4)**5*x(5) - 4641*x(1)*x(2)*x(4)**4*x(5)**2 + 765*x(1)*x(2)*x(4)**3*x(5)**3 - 35*x(1)*x(2)*x(4)**2*x(5)**4 - 45*x(1)*x(2)*x(4)*x(5)**5 - 3324*x(1)*x(3)**6*x(4) - 3324*x(1)*x(3)**6*x(5) + 13958*x(1)*x(3)**5*x(4)**2 + 12296*x(1)*x(3)**5*x(4)*x(5) + 1862*x(1)*x(3)**5*x(5)**2 - 22560*x(1)*x(3)**4*x(4)**3 - 12978*x(1)*x(3)**4*x(4)**2*x(5) - 822*x(1)*x(3)**4*x(4)*x(5)**2 + 1672*x(1)*x(3)**4*x(5)**3 + 15430*x(1)*x(3)**3*x(4)**4 - 6153*x(1)*x(3)**3*x(4)**3*x(5) - 10607*x(1)*x(3)**3*x(4)**2*x(5)**2 - 6773*x(1)*x(3)**3*x(4)*x(5)**3 - 320*x(1)*x(3)**3*x(5)**4 - 780*x(1)*x(3)**2*x(4)**5 + 23747*x(1)*x(3)**2*x(4)**4*x(5) + 21343*x(1)*x(3)**2*x(4)**3*x(5)**2 + 9697*x(1)*x(3)**2*x(4)**2*x(5)**3 + 330*x(1)*x(3)**2*x(4)*x(5)**4 + 110*x(1)*x(3)**2*x(5)**5 - 4386*x(1)*x(3)*x(4)**6 - 18644*x(1)*x(3)*x(4)**5*x(5) - 16407*x(1)*x(3)*x(4)**4*x(5)**2 - 6003*x(1)*x(3)*x(4)**3*x(5)**3 + 90*x(1)*x(3)*x(4)**2*x(5)**4 - 145*x(1)*x(3)*x(4)*x(5)**5 + 1662*x(1)*x(4)**7 + 5086*x(1)*x(4)**6*x(5) + 4781*x(1)*x(4)**5*x(5)**2 + 1247*x(1)*x(4)**4*x(5)**3 - 210*x(1)*x(4)**3*x(5)**4 + 125*x(1)*x(4)**2*x(5)**5 - 100*x(2)**5*x(3)**3 + 290*x(2)**5*x(3)**2*x(4) + 110*x(2)**5*x(3)**2*x(5) - 290*x(2)**5*x(3)*x(4)**2 - 185*x(2)**5*x(3)*x(4)*x(5) - 20*x(2)**5*x(3)*x(5)**2 + 100*x(2)**5*x(4)**3 + 95*x(2)**5*x(4)**2*x(5) - 10*x(2)**5*x(4)*x(5)**2 + 10*x(2)**5*x(5)**3 + 300*x(2)**4*x(3)**4 - 480*x(2)**4*x(3)**3*x(4) - 320*x(2)**4*x(3)**3*x(5) - 290*x(2)**4*x(3)**2*x(4)**2 + 160*x(2)**4*x(3)**2*x(4)*x(5) + 20*x(2)**4*x(3)**2*x(5)**2 + 870*x(2)**4*x(3)*x(4)**3 + 455*x(2)**4*x(3)*x(4)**2*x(5) + 45*x(2)**4*x(3)*x(4)*x(5)**2 + 20*x(2)**4*x(3)*x(5)**3 - 400*x(2)**4*x(4)**4 - 385*x(2)**4*x(4)**3*x(5) + 35*x(2)**4*x(4)**2*x(5)**2 - 10*x(2)**4*x(4)*x(5)**3 - 20*x(2)**4*x(5)**4 - 1862*x(2)**3*x(3)**5 + 7700*x(2)**3*x(3)**4*x(4) + 1672*x(2)**3*x(3)**4*x(5) - 14170*x(2)**3*x(3)**3*x(4)**2 - 5513*x(2)**3*x(3)**3*x(4)*x(5) + 420*x(2)**3*x(3)**3*x(5)**2 + 14750*x(2)**3*x(3)**2*x(4)**3 + 8447*x(2)**3*x(3)**2*x(4)**2*x(5) - 710*x(2)**3*x(3)**2*x(4)*x(5)**2 - 260*x(2)**3*x(3)**2*x(5)**3 - 8580*x(2)**3*x(3)*x(4)**4 - 6703*x(2)**3*x(3)*x(4)**3*x(5) + 175*x(2)**3*x(3)*x(4)**2*x(5)**2 + 320*x(2)**3*x(3)*x(4)*x(5)**3 + 20*x(2)**3*x(3)*x(5)**4 + 2162*x(2)**3*x(4)**5 + 2247*x(2)**3*x(4)**4*x(5) + 55*x(2)**3*x(4)**3*x(5)**2 - 230*x(2)**3*x(4)**2*x(5)**3 + 50*x(2)**3*x(4)*x(5)**4 + 10*x(2)**3*x(5)**5 + 1662*x(2)**2*x(3)**6 - 4186*x(2)**2*x(3)**5*x(4) + 1862*x(2)**2*x(3)**5*x(5) - 870*x(2)**2*x(3)**4*x(4)**2 - 11644*x(2)**2*x(3)**4*x(4)*x(5) - 3944*x(2)**2*x(3)**4*x(5)**2 + 14750*x(2)**2*x(3)**3*x(4)**3 + 27416*x(2)**2*x(3)**3*x(4)**2*x(5) + 12826*x(2)**2*x(3)**3*x(4)*x(5)**2 + 420*x(2)**2*x(3)**3*x(5)**3 - 21880*x(2)**2*x(3)**2*x(4)**4 - 32854*x(2)**2*x(3)**2*x(4)**3*x(5) - 17414*x(2)**2*x(3)**2*x(4)**2*x(5)**2 - 150*x(2)**2*x(3)**2*x(4)*x(5)**3 + 20*x(2)**2*x(3)**2*x(5)**4 + 14048*x(2)**2*x(3)*x(4)**5 + 20491*x(2)**2*x(3)*x(4)**4*x(5) + 11711*x(2)**2*x(3)*x(4)**3*x(5)**2 - 495*x(2)**2*x(3)*x(4)**2*x(5)**3 - 35*x(2)**2*x(3)*x(4)*x(5)**4 - 20*x(2)**2*x(3)*x(5)**5 - 3524*x(2)**2*x(4)**6 - 5381*x(2)**2*x(4)**5*x(5) - 3299*x(2)**2*x(4)**4*x(5)**2 + 545*x(2)**2*x(4)**3*x(5)**3 - 25*x(2)**2*x(4)**2*x(5)**4 - 30*x(2)**2*x(4)*x(5)**5 - 3324*x(2)*x(3)**6*x(4) - 3324*x(2)*x(3)**6*x(5) + 13958*x(2)*x(3)**5*x(4)**2 + 12296*x(2)*x(3)**5*x(4)*x(5) + 1862*x(2)*x(3)**5*x(5)**2 - 22560*x(2)*x(3)**4*x(4)**3 - 12978*x(2)*x(3)**4*x(4)**2*x(5) - 822*x(2)*x(3)**4*x(4)*x(5)**2 + 1672*x(2)*x(3)**4*x(5)**3 + 15430*x(2)*x(3)**3*x(4)**4 - 6153*x(2)*x(3)**3*x(4)**3*x(5) - 10607*x(2)*x(3)**3*x(4)**2*x(5)**2 - 6773*x(2)*x(3)**3*x(4)*x(5)**3 - 320*x(2)*x(3)**3*x(5)**4 - 780*x(2)*x(3)**2*x(4)**5 + 23747*x(2)*x(3)**2*x(4)**4*x(5) + 21343*x(2)*x(3)**2*x(4)**3*x(5)**2 + 9697*x(2)*x(3)**2*x(4)**2*x(5)**3 + 330*x(2)*x(3)**2*x(4)*x(5)**4 + 110*x(2)*x(3)**2*x(5)**5 - 4386*x(2)*x(3)*x(4)**6 - 18644*x(2)*x(3)*x(4)**5*x(5) - 16407*x(2)*x(3)*x(4)**4*x(5)**2 - 6003*x(2)*x(3)*x(4)**3*x(5)**3 + 90*x(2)*x(3)*x(4)**2*x(5)**4 - 145*x(2)*x(3)*x(4)*x(5)**5 + 1662*x(2)*x(4)**7 + 5086*x(2)*x(4)**6*x(5) + 4781*x(2)*x(4)**5*x(5)**2 + 1247*x(2)*x(4)**4*x(5)**3 - 210*x(2)*x(4)**3*x(5)**4 + 125*x(2)*x(4)**2*x(5)**5 + 1662*x(3)**6*x(4)**2 + 3324*x(3)**6*x(4)*x(5) + 1662*x(3)**6*x(5)**2 - 7910*x(3)**5*x(4)**3 - 14158*x(3)**5*x(4)**2*x(5) - 8110*x(3)**5*x(4)*x(5)**2 - 1862*x(3)**5*x(5)**3 + 15430*x(3)**4*x(4)**4 + 22950*x(3)**4*x(4)**3*x(5) + 13948*x(3)**4*x(4)**2*x(5)**2 + 4766*x(3)**4*x(4)*x(5)**3 + 300*x(3)**4*x(5)**4 - 15430*x(3)**3*x(4)**5 - 15430*x(3)**3*x(4)**4*x(5) - 8887*x(3)**3*x(4)**3*x(5)**2 - 2639*x(3)**3*x(4)**2*x(5)**3 - 60*x(3)**3*x(4)*x(5)**4 - 100*x(3)**3*x(5)**5 + 7910*x(3)**2*x(4)**6 + 390*x(3)**2*x(4)**5*x(5) - 1577*x(3)**2*x(4)**4*x(5)**2 - 3239*x(3)**2*x(4)**3*x(5)**3 - 440*x(3)**2*x(4)**2*x(5)**4 + 80*x(3)**2*x(4)*x(5)**5 - 1662*x(3)*x(4)**7 + 4586*x(3)*x(4)**6*x(5) + 4496*x(3)*x(4)**5*x(5)**2 + 4496*x(3)*x(4)**4*x(5)**3 + 125*x(3)*x(4)**3*x(5)**4 + 65*x(3)*x(4)**2*x(5)**5 - 1662*x(4)**7*x(5) - 1562*x(4)**6*x(5)**2 - 1562*x(4)**5*x(5)**3 + 205*x(4)**4*x(5)**4 - 105*x(4)**3*x(5)**5))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(1,7) = (4*(x(3) - x(4))**2*(- 5*x(1)**4*x(2)*x(3) - 5*x(1)**4*x(2)*x(4) + 10*x(1)**4*x(2)*x(5) + 100*x(1)**4*x(3)**2 - 190*x(1)**4*x(3)*x(4) - 5*x(1)**4*x(3)*x(5) + 100*x(1)**4*x(4)**2 - 5*x(1)**4*x(4)*x(5) - 5*x(1)**3*x(2)**2*x(3) - 5*x(1)**3*x(2)**2*x(4) + 10*x(1)**3*x(2)**2*x(5) + 210*x(1)**3*x(2)*x(3)**2 - 370*x(1)**3*x(2)*x(3)*x(4) - 20*x(1)**3*x(2)*x(3)*x(5) + 210*x(1)**3*x(2)*x(4)**2 - 20*x(1)**3*x(2)*x(4)*x(5) - 10*x(1)**3*x(2)*x(5)**2 - 300*x(1)**3*x(3)**3 + 280*x(1)**3*x(3)**2*x(4) + 10*x(1)**3*x(3)**2*x(5) + 280*x(1)**3*x(3)*x(4)**2 + 10*x(1)**3*x(3)*x(4)*x(5) + 5*x(1)**3*x(3)*x(5)**2 - 300*x(1)**3*x(4)**3 + 10*x(1)**3*x(4)**2*x(5) + 5*x(1)**3*x(4)*x(5)**2 - 5*x(1)**2*x(2)**3*x(3) - 5*x(1)**2*x(2)**3*x(4) + 10*x(1)**2*x(2)**3*x(5) + 210*x(1)**2*x(2)**2*x(3)**2 - 370*x(1)**2*x(2)**2*x(3)*x(4) - 20*x(1)**2*x(2)**2*x(3)*x(5) + 210*x(1)**2*x(2)**2*x(4)**2 - 20*x(1)**2*x(2)**2*x(4)*x(5) - 10*x(1)**2*x(2)**2*x(5)**2 - 505*x(1)**2*x(2)*x(3)**3 + 455*x(1)**2*x(2)*x(3)**2*x(4) + 10*x(1)**2*x(2)*x(3)**2*x(5) + 455*x(1)**2*x(2)*x(3)*x(4)**2 + 30*x(1)**2*x(2)*x(3)*x(4)*x(5) + 25*x(1)**2*x(2)*x(3)*x(5)**2 - 505*x(1)**2*x(2)*x(4)**3 + 10*x(1)**2*x(2)*x(4)**2*x(5) + 25*x(1)**2*x(2)*x(4)*x(5)**2 + 1862*x(1)**2*x(3)**4 - 6138*x(1)**2*x(3)**3*x(4) + 95*x(1)**2*x(3)**3*x(5) + 8612*x(1)**2*x(3)**2*x(4)**2 - 105*x(1)**2*x(3)**2*x(4)*x(5) - 110*x(1)**2*x(3)**2*x(5)**2 - 6138*x(1)**2*x(3)*x(4)**3 - 105*x(1)**2*x(3)*x(4)**2*x(5) + 180*x(1)**2*x(3)*x(4)*x(5)**2 + 1862*x(1)**2*x(4)**4 + 95*x(1)**2*x(4)**3*x(5) - 110*x(1)**2*x(4)**2*x(5)**2 - 5*x(1)*x(2)**4*x(3) - 5*x(1)*x(2)**4*x(4) + 10*x(1)*x(2)**4*x(5) + 210*x(1)*x(2)**3*x(3)**2 - 370*x(1)*x(2)**3*x(3)*x(4) - 20*x(1)*x(2)**3*x(3)*x(5) + 210*x(1)*x(2)**3*x(4)**2 - 20*x(1)*x(2)**3*x(4)*x(5) - 10*x(1)*x(2)**3*x(5)**2 - 505*x(1)*x(2)**2*x(3)**3 + 455*x(1)*x(2)**2*x(3)**2*x(4) + 10*x(1)*x(2)**2*x(3)**2*x(5) + 455*x(1)*x(2)**2*x(3)*x(4)**2 + 30*x(1)*x(2)**2*x(3)*x(4)*x(5) + 25*x(1)*x(2)**2*x(3)*x(5)**2 - 505*x(1)*x(2)**2*x(4)**3 + 10*x(1)*x(2)**2*x(4)**2*x(5) + 25*x(1)*x(2)**2*x(4)*x(5)**2 + 1962*x(1)*x(2)*x(3)**4 - 6028*x(1)*x(2)*x(3)**3*x(4) + 200*x(1)*x(2)*x(3)**3*x(5) + 8232*x(1)*x(2)*x(3)**2*x(4)**2 - 200*x(1)*x(2)*x(3)**2*x(4)*x(5) - 220*x(1)*x(2)*x(3)**2*x(5)**2 - 6028*x(1)*x(2)*x(3)*x(4)**3 - 200*x(1)*x(2)*x(3)*x(4)**2*x(5) + 340*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 1962*x(1)*x(2)*x(4)**4 + 200*x(1)*x(2)*x(4)**3*x(5) - 220*x(1)*x(2)*x(4)**2*x(5)**2 - 1662*x(1)*x(3)**5 + 4386*x(1)*x(3)**4*x(4) - 1762*x(1)*x(3)**4*x(5) - 2744*x(1)*x(3)**3*x(4)**2 + 6248*x(1)*x(3)**3*x(4)*x(5) + 205*x(1)*x(3)**3*x(5)**2 - 2744*x(1)*x(3)**2*x(4)**3 - 8992*x(1)*x(3)**2*x(4)**2*x(5) - 175*x(1)*x(3)**2*x(4)*x(5)**2 + 4386*x(1)*x(3)*x(4)**4 + 6248*x(1)*x(3)*x(4)**3*x(5) - 175*x(1)*x(3)*x(4)**2*x(5)**2 - 1662*x(1)*x(4)**5 - 1762*x(1)*x(4)**4*x(5) + 205*x(1)*x(4)**3*x(5)**2 + 100*x(2)**4*x(3)**2 - 190*x(2)**4*x(3)*x(4) - 5*x(2)**4*x(3)*x(5) + 100*x(2)**4*x(4)**2 - 5*x(2)**4*x(4)*x(5) - 300*x(2)**3*x(3)**3 + 280*x(2)**3*x(3)**2*x(4) + 10*x(2)**3*x(3)**2*x(5) + 280*x(2)**3*x(3)*x(4)**2 + 10*x(2)**3*x(3)*x(4)*x(5) + 5*x(2)**3*x(3)*x(5)**2 - 300*x(2)**3*x(4)**3 + 10*x(2)**3*x(4)**2*x(5) + 5*x(2)**3*x(4)*x(5)**2 + 1862*x(2)**2*x(3)**4 - 6138*x(2)**2*x(3)**3*x(4) + 95*x(2)**2*x(3)**3*x(5) + 8612*x(2)**2*x(3)**2*x(4)**2 - 105*x(2)**2*x(3)**2*x(4)*x(5) - 110*x(2)**2*x(3)**2*x(5)**2 - 6138*x(2)**2*x(3)*x(4)**3 - 105*x(2)**2*x(3)*x(4)**2*x(5) + 180*x(2)**2*x(3)*x(4)*x(5)**2 + 1862*x(2)**2*x(4)**4 + 95*x(2)**2*x(4)**3*x(5) - 110*x(2)**2*x(4)**2*x(5)**2 - 1662*x(2)*x(3)**5 + 4386*x(2)*x(3)**4*x(4) - 1762*x(2)*x(3)**4*x(5) - 2744*x(2)*x(3)**3*x(4)**2 + 6248*x(2)*x(3)**3*x(4)*x(5) + 205*x(2)*x(3)**3*x(5)**2 - 2744*x(2)*x(3)**2*x(4)**3 - 8992*x(2)*x(3)**2*x(4)**2*x(5) - 175*x(2)*x(3)**2*x(4)*x(5)**2 + 4386*x(2)*x(3)*x(4)**4 + 6248*x(2)*x(3)*x(4)**3*x(5) - 175*x(2)*x(3)*x(4)**2*x(5)**2 - 1662*x(2)*x(4)**5 - 1762*x(2)*x(4)**4*x(5) + 205*x(2)*x(4)**3*x(5)**2 + 1662*x(3)**5*x(4) + 1662*x(3)**5*x(5) - 6248*x(3)**4*x(4)**2 - 4586*x(3)**4*x(4)*x(5) - 100*x(3)**4*x(5)**2 + 9182*x(3)**3*x(4)**3 + 2934*x(3)**3*x(4)**2*x(5) - 10*x(3)**3*x(4)*x(5)**2 - 6248*x(3)**2*x(4)**4 + 2934*x(3)**2*x(4)**3*x(5) + 190*x(3)**2*x(4)**2*x(5)**2 + 1662*x(3)*x(4)**5 - 4586*x(3)*x(4)**4*x(5) - 10*x(3)*x(4)**3*x(5)**2 + 1662*x(4)**5*x(5) - 100*x(4)**4*x(5)**2))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(1,8) = (24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(3) - x(4))**6 + (x(3)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(3) - x(4)))/5 - (x(4)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(3) - x(4)))/5 + x(3)*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(4)*(x(3) - x(4))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(3)**3*(x(3) - x(4))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 - (x(4)**3*(x(3) - x(4))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 + (x(3)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(3) - x(4))**3)/3 - (x(4)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(3) - x(4))**3)/3 + x(3)*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(4)*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(3)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(4)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(3) - x(4))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - x(3)**2*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(4)**2*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - (x(3)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2 + (x(4)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(3) - x(4))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2

                            Cb(1,9) = -(4*(x(3) - x(4))**2*(10*x(1)**3*x(2)**3 - 10*x(1)**3*x(2)**2*x(3) - 20*x(1)**3*x(2)**2*x(4) + 95*x(1)**3*x(2)*x(3)**2 - 185*x(1)**3*x(2)*x(3)*x(4) + 15*x(1)**3*x(2)*x(3)*x(5) + 110*x(1)**3*x(2)*x(4)**2 + 5*x(1)**3*x(2)*x(4)*x(5) - 10*x(1)**3*x(2)*x(5)**2 + 100*x(1)**3*x(3)**3 - 290*x(1)**3*x(3)**2*x(4) - 105*x(1)**3*x(3)**2*x(5) + 290*x(1)**3*x(3)*x(4)**2 + 185*x(1)**3*x(3)*x(4)*x(5) + 5*x(1)**3*x(3)*x(5)**2 - 100*x(1)**3*x(4)**3 - 100*x(1)**3*x(4)**2*x(5) + 5*x(1)**3*x(4)*x(5)**2 - 10*x(1)**2*x(2)**3*x(3) - 20*x(1)**2*x(2)**3*x(4) + 205*x(1)**2*x(2)**2*x(3)**2 - 365*x(1)**2*x(2)**2*x(3)*x(4) + 15*x(1)**2*x(2)**2*x(3)*x(5) + 240*x(1)**2*x(2)**2*x(4)**2 + 5*x(1)**2*x(2)**2*x(4)*x(5) - 10*x(1)**2*x(2)**2*x(5)**2 + 5*x(1)**2*x(2)*x(3)**3 - 495*x(1)**2*x(2)*x(3)**2*x(4) - 215*x(1)**2*x(2)*x(3)**2*x(5) + 960*x(1)**2*x(2)*x(3)*x(4)**2 + 355*x(1)**2*x(2)*x(3)*x(4)*x(5) - 520*x(1)**2*x(2)*x(4)**3 - 210*x(1)**2*x(2)*x(4)**2*x(5) + 20*x(1)**2*x(2)*x(4)*x(5)**2 + 10*x(1)**2*x(2)*x(5)**3 - 200*x(1)**2*x(3)**4 + 290*x(1)**2*x(3)**3*x(4) + 205*x(1)**2*x(3)**3*x(5) + 290*x(1)**2*x(3)**2*x(4)**2 - 85*x(1)**2*x(3)**2*x(4)*x(5) - 680*x(1)**2*x(3)*x(4)**3 - 370*x(1)**2*x(3)*x(4)**2*x(5) - 5*x(1)**2*x(3)*x(5)**3 + 300*x(1)**2*x(4)**4 + 300*x(1)**2*x(4)**3*x(5) - 10*x(1)**2*x(4)**2*x(5)**2 - 5*x(1)**2*x(4)*x(5)**3 + 95*x(1)*x(2)**3*x(3)**2 - 185*x(1)*x(2)**3*x(3)*x(4) + 15*x(1)*x(2)**3*x(3)*x(5) + 110*x(1)*x(2)**3*x(4)**2 + 5*x(1)*x(2)**3*x(4)*x(5) - 10*x(1)*x(2)**3*x(5)**2 + 5*x(1)*x(2)**2*x(3)**3 - 495*x(1)*x(2)**2*x(3)**2*x(4) - 215*x(1)*x(2)**2*x(3)**2*x(5) + 960*x(1)*x(2)**2*x(3)*x(4)**2 + 355*x(1)*x(2)**2*x(3)*x(4)*x(5) - 520*x(1)*x(2)**2*x(4)**3 - 210*x(1)*x(2)**2*x(4)**2*x(5) + 20*x(1)*x(2)**2*x(4)*x(5)**2 + 10*x(1)*x(2)**2*x(5)**3 + 1362*x(1)*x(2)*x(3)**4 - 5768*x(1)*x(2)*x(3)**3*x(4) + 300*x(1)*x(2)*x(3)**3*x(5) + 9677*x(1)*x(2)*x(3)**2*x(4)**2 - 70*x(1)*x(2)*x(3)**2*x(4)*x(5) + 15*x(1)*x(2)*x(3)**2*x(5)**2 - 7513*x(1)*x(2)*x(3)*x(4)**3 - 655*x(1)*x(2)*x(3)*x(4)**2*x(5) + 15*x(1)*x(2)*x(3)*x(4)*x(5)**2 - 15*x(1)*x(2)*x(3)*x(5)**3 + 2272*x(1)*x(2)*x(4)**4 + 505*x(1)*x(2)*x(4)**3*x(5) - 10*x(1)*x(2)*x(4)**2*x(5)**2 - 25*x(1)*x(2)*x(4)*x(5)**3 + 1662*x(1)*x(3)**5 - 7710*x(1)*x(3)**4*x(4) - 1562*x(1)*x(3)**4*x(5) + 15040*x(1)*x(3)**3*x(4)**2 + 5948*x(1)*x(3)**3*x(4)*x(5) - 205*x(1)*x(3)**3*x(5)**2 - 15430*x(1)*x(3)**2*x(4)**3 - 9087*x(1)*x(3)**2*x(4)**2*x(5) + 285*x(1)*x(3)**2*x(4)*x(5)**2 + 105*x(1)*x(3)**2*x(5)**3 + 8300*x(1)*x(3)*x(4)**4 + 6533*x(1)*x(3)*x(4)**3*x(5) - 15*x(1)*x(3)*x(4)**2*x(5)**2 - 185*x(1)*x(3)*x(4)*x(5)**3 - 1862*x(1)*x(4)**5 - 1862*x(1)*x(4)**4*x(5) - 95*x(1)*x(4)**3*x(5)**2 + 110*x(1)*x(4)**2*x(5)**3 + 100*x(2)**3*x(3)**3 - 290*x(2)**3*x(3)**2*x(4) - 105*x(2)**3*x(3)**2*x(5) + 290*x(2)**3*x(3)*x(4)**2 + 185*x(2)**3*x(3)*x(4)*x(5) + 5*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 - 100*x(2)**3*x(4)**2*x(5) + 5*x(2)**3*x(4)*x(5)**2 - 200*x(2)**2*x(3)**4 + 290*x(2)**2*x(3)**3*x(4) + 205*x(2)**2*x(3)**3*x(5) + 290*x(2)**2*x(3)**2*x(4)**2 - 85*x(2)**2*x(3)**2*x(4)*x(5) - 680*x(2)**2*x(3)*x(4)**3 - 370*x(2)**2*x(3)*x(4)**2*x(5) - 5*x(2)**2*x(3)*x(5)**3 + 300*x(2)**2*x(4)**4 + 300*x(2)**2*x(4)**3*x(5) - 10*x(2)**2*x(4)**2*x(5)**2 - 5*x(2)**2*x(4)*x(5)**3 + 1662*x(2)*x(3)**5 - 7710*x(2)*x(3)**4*x(4) - 1562*x(2)*x(3)**4*x(5) + 15040*x(2)*x(3)**3*x(4)**2 + 5948*x(2)*x(3)**3*x(4)*x(5) - 205*x(2)*x(3)**3*x(5)**2 - 15430*x(2)*x(3)**2*x(4)**3 - 9087*x(2)*x(3)**2*x(4)**2*x(5) + 285*x(2)*x(3)**2*x(4)*x(5)**2 + 105*x(2)*x(3)**2*x(5)**3 + 8300*x(2)*x(3)*x(4)**4 + 6533*x(2)*x(3)*x(4)**3*x(5) - 15*x(2)*x(3)*x(4)**2*x(5)**2 - 185*x(2)*x(3)*x(4)*x(5)**3 - 1862*x(2)*x(4)**5 - 1862*x(2)*x(4)**4*x(5) - 95*x(2)*x(4)**3*x(5)**2 + 110*x(2)*x(4)**2*x(5)**3 - 1662*x(3)**5*x(4) - 1662*x(3)**5*x(5) + 7910*x(3)**4*x(4)**2 + 7910*x(3)**4*x(4)*x(5) + 1762*x(3)**4*x(5)**2 - 15430*x(3)**3*x(4)**3 - 15430*x(3)**3*x(4)**2*x(5) - 6338*x(3)**3*x(4)*x(5)**2 - 100*x(3)**3*x(5)**3 + 15430*x(3)**2*x(4)**4 + 15430*x(3)**2*x(4)**3*x(5) + 9087*x(3)**2*x(4)**2*x(5)**2 + 90*x(3)**2*x(4)*x(5)**3 - 7910*x(3)*x(4)**5 - 7910*x(3)*x(4)**4*x(5) - 6143*x(3)*x(4)**3*x(5)**2 + 95*x(3)*x(4)**2*x(5)**3 + 1662*x(4)**6 + 1662*x(4)**5*x(5) + 1662*x(4)**4*x(5)**2 - 105*x(4)**3*x(5)**3))/(5*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            Cb(1,10) = (4*(x(3) - x(4))**2*(5*x(1)**2*x(2)**2 - 5*x(1)**2*x(2)*x(3) - 5*x(1)**2*x(2)*x(4) + 50*x(1)**2*x(3)**2 - 95*x(1)**2*x(3)*x(4) + 50*x(1)**2*x(4)**2 - 5*x(1)*x(2)**2*x(3) - 5*x(1)*x(2)**2*x(4) + 105*x(1)*x(2)*x(3)**2 - 190*x(1)*x(2)*x(3)*x(4) + 105*x(1)*x(2)*x(4)**2 - 100*x(1)*x(3)**3 + 95*x(1)*x(3)**2*x(4) + 95*x(1)*x(3)*x(4)**2 - 100*x(1)*x(4)**3 + 50*x(2)**2*x(3)**2 - 95*x(2)**2*x(3)*x(4) + 50*x(2)**2*x(4)**2 - 100*x(2)*x(3)**3 + 95*x(2)*x(3)**2*x(4) + 95*x(2)*x(3)*x(4)**2 - 100*x(2)*x(4)**3 + 831*x(3)**4 - 3124*x(3)**3*x(4) + 4591*x(3)**2*x(4)**2 - 3124*x(3)*x(4)**3 + 831*x(4)**4))/(5*(x(1) - x(5))**2*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            x = x_cb(i-1:i+3)
                            x = x - x(1)

                            Cb(2,1) = (4*(x(2) - x(3))**2*(831*x(2)**4 - 3124*x(2)**3*x(3) - 100*x(2)**3*x(4) - 100*x(2)**3*x(5) + 4591*x(2)**2*x(3)**2 + 95*x(2)**2*x(3)*x(4) + 95*x(2)**2*x(3)*x(5) + 50*x(2)**2*x(4)**2 + 105*x(2)**2*x(4)*x(5) + 50*x(2)**2*x(5)**2 - 3124*x(2)*x(3)**3 + 95*x(2)*x(3)**2*x(4) + 95*x(2)*x(3)**2*x(5) - 95*x(2)*x(3)*x(4)**2 - 190*x(2)*x(3)*x(4)*x(5) - 95*x(2)*x(3)*x(5)**2 - 5*x(2)*x(4)**2*x(5) - 5*x(2)*x(4)*x(5)**2 + 831*x(3)**4 - 100*x(3)**3*x(4) - 100*x(3)**3*x(5) + 50*x(3)**2*x(4)**2 + 105*x(3)**2*x(4)*x(5) + 50*x(3)**2*x(5)**2 - 5*x(3)*x(4)**2*x(5) - 5*x(3)*x(4)*x(5)**2 + 5*x(4)**2*x(5)**2))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2)

                            Cb(2,2) = -(4*(x(2) - x(3))**2*(- 105*x(1)**3*x(2)**3 + 95*x(1)**3*x(2)**2*x(3) + 110*x(1)**3*x(2)**2*x(4) + 110*x(1)**3*x(2)**2*x(5) + 90*x(1)**3*x(2)*x(3)**2 - 185*x(1)**3*x(2)*x(3)*x(4) - 185*x(1)**3*x(2)*x(3)*x(5) - 5*x(1)**3*x(2)*x(4)**2 - 25*x(1)**3*x(2)*x(4)*x(5) - 5*x(1)**3*x(2)*x(5)**2 - 100*x(1)**3*x(3)**3 + 105*x(1)**3*x(3)**2*x(4) + 105*x(1)**3*x(3)**2*x(5) - 5*x(1)**3*x(3)*x(4)**2 - 15*x(1)**3*x(3)*x(4)*x(5) - 5*x(1)**3*x(3)*x(5)**2 + 10*x(1)**3*x(4)**2*x(5) + 10*x(1)**3*x(4)*x(5)**2 + 1662*x(1)**2*x(2)**4 - 6143*x(1)**2*x(2)**3*x(3) - 95*x(1)**2*x(2)**3*x(4) - 95*x(1)**2*x(2)**3*x(5) + 9087*x(1)**2*x(2)**2*x(3)**2 - 15*x(1)**2*x(2)**2*x(3)*x(4) - 15*x(1)**2*x(2)**2*x(3)*x(5) - 10*x(1)**2*x(2)**2*x(4)**2 - 10*x(1)**2*x(2)**2*x(4)*x(5) - 10*x(1)**2*x(2)**2*x(5)**2 - 6338*x(1)**2*x(2)*x(3)**3 + 285*x(1)**2*x(2)*x(3)**2*x(4) + 285*x(1)**2*x(2)*x(3)**2*x(5) + 15*x(1)**2*x(2)*x(3)*x(4)*x(5) + 5*x(1)**2*x(2)*x(4)**3 + 20*x(1)**2*x(2)*x(4)**2*x(5) + 20*x(1)**2*x(2)*x(4)*x(5)**2 + 5*x(1)**2*x(2)*x(5)**3 + 1762*x(1)**2*x(3)**4 - 205*x(1)**2*x(3)**3*x(4) - 205*x(1)**2*x(3)**3*x(5) + 15*x(1)**2*x(3)**2*x(4)*x(5) + 5*x(1)**2*x(3)*x(4)**3 + 5*x(1)**2*x(3)*x(5)**3 - 10*x(1)**2*x(4)**3*x(5) - 10*x(1)**2*x(4)**2*x(5)**2 - 10*x(1)**2*x(4)*x(5)**3 + 1662*x(1)*x(2)**5 - 7910*x(1)*x(2)**4*x(3) - 1862*x(1)*x(2)**4*x(4) - 1862*x(1)*x(2)**4*x(5) + 15430*x(1)*x(2)**3*x(3)**2 + 6533*x(1)*x(2)**3*x(3)*x(4) + 6533*x(1)*x(2)**3*x(3)*x(5) + 300*x(1)*x(2)**3*x(4)**2 + 505*x(1)*x(2)**3*x(4)*x(5) + 300*x(1)*x(2)**3*x(5)**2 - 15430*x(1)*x(2)**2*x(3)**3 - 9087*x(1)*x(2)**2*x(3)**2*x(4) - 9087*x(1)*x(2)**2*x(3)**2*x(5) - 370*x(1)*x(2)**2*x(3)*x(4)**2 - 655*x(1)*x(2)**2*x(3)*x(4)*x(5) - 370*x(1)*x(2)**2*x(3)*x(5)**2 - 100*x(1)*x(2)**2*x(4)**3 - 210*x(1)*x(2)**2*x(4)**2*x(5) - 210*x(1)*x(2)**2*x(4)*x(5)**2 - 100*x(1)*x(2)**2*x(5)**3 + 7910*x(1)*x(2)*x(3)**4 + 5948*x(1)*x(2)*x(3)**3*x(4) + 5948*x(1)*x(2)*x(3)**3*x(5) - 85*x(1)*x(2)*x(3)**2*x(4)**2 - 70*x(1)*x(2)*x(3)**2*x(4)*x(5) - 85*x(1)*x(2)*x(3)**2*x(5)**2 + 185*x(1)*x(2)*x(3)*x(4)**3 + 355*x(1)*x(2)*x(3)*x(4)**2*x(5) + 355*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 185*x(1)*x(2)*x(3)*x(5)**3 + 5*x(1)*x(2)*x(4)**3*x(5) + 5*x(1)*x(2)*x(4)**2*x(5)**2 + 5*x(1)*x(2)*x(4)*x(5)**3 - 1662*x(1)*x(3)**5 - 1562*x(1)*x(3)**4*x(4) - 1562*x(1)*x(3)**4*x(5) + 205*x(1)*x(3)**3*x(4)**2 + 300*x(1)*x(3)**3*x(4)*x(5) + 205*x(1)*x(3)**3*x(5)**2 - 105*x(1)*x(3)**2*x(4)**3 - 215*x(1)*x(3)**2*x(4)**2*x(5) - 215*x(1)*x(3)**2*x(4)*x(5)**2 - 105*x(1)*x(3)**2*x(5)**3 + 15*x(1)*x(3)*x(4)**3*x(5) + 15*x(1)*x(3)*x(4)**2*x(5)**2 + 15*x(1)*x(3)*x(4)*x(5)**3 + 1662*x(2)**6 - 7910*x(2)**5*x(3) - 1862*x(2)**5*x(4) - 1862*x(2)**5*x(5) + 15430*x(2)**4*x(3)**2 + 8300*x(2)**4*x(3)*x(4) + 8300*x(2)**4*x(3)*x(5) + 300*x(2)**4*x(4)**2 + 2272*x(2)**4*x(4)*x(5) + 300*x(2)**4*x(5)**2 - 15430*x(2)**3*x(3)**3 - 15430*x(2)**3*x(3)**2*x(4) - 15430*x(2)**3*x(3)**2*x(5) - 680*x(2)**3*x(3)*x(4)**2 - 7513*x(2)**3*x(3)*x(4)*x(5) - 680*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 - 520*x(2)**3*x(4)**2*x(5) - 520*x(2)**3*x(4)*x(5)**2 - 100*x(2)**3*x(5)**3 + 7910*x(2)**2*x(3)**4 + 15040*x(2)**2*x(3)**3*x(4) + 15040*x(2)**2*x(3)**3*x(5) + 290*x(2)**2*x(3)**2*x(4)**2 + 9677*x(2)**2*x(3)**2*x(4)*x(5) + 290*x(2)**2*x(3)**2*x(5)**2 + 290*x(2)**2*x(3)*x(4)**3 + 960*x(2)**2*x(3)*x(4)**2*x(5) + 960*x(2)**2*x(3)*x(4)*x(5)**2 + 290*x(2)**2*x(3)*x(5)**3 + 110*x(2)**2*x(4)**3*x(5) + 240*x(2)**2*x(4)**2*x(5)**2 + 110*x(2)**2*x(4)*x(5)**3 - 1662*x(2)*x(3)**5 - 7710*x(2)*x(3)**4*x(4) - 7710*x(2)*x(3)**4*x(5) + 290*x(2)*x(3)**3*x(4)**2 - 5768*x(2)*x(3)**3*x(4)*x(5) + 290*x(2)*x(3)**3*x(5)**2 - 290*x(2)*x(3)**2*x(4)**3 - 495*x(2)*x(3)**2*x(4)**2*x(5) - 495*x(2)*x(3)**2*x(4)*x(5)**2 - 290*x(2)*x(3)**2*x(5)**3 - 185*x(2)*x(3)*x(4)**3*x(5) - 365*x(2)*x(3)*x(4)**2*x(5)**2 - 185*x(2)*x(3)*x(4)*x(5)**3 - 20*x(2)*x(4)**3*x(5)**2 - 20*x(2)*x(4)**2*x(5)**3 + 1662*x(3)**5*x(4) + 1662*x(3)**5*x(5) - 200*x(3)**4*x(4)**2 + 1362*x(3)**4*x(4)*x(5) - 200*x(3)**4*x(5)**2 + 100*x(3)**3*x(4)**3 + 5*x(3)**3*x(4)**2*x(5) + 5*x(3)**3*x(4)*x(5)**2 + 100*x(3)**3*x(5)**3 + 95*x(3)**2*x(4)**3*x(5) + 205*x(3)**2*x(4)**2*x(5)**2 + 95*x(3)**2*x(4)*x(5)**3 - 10*x(3)*x(4)**3*x(5)**2 - 10*x(3)*x(4)**2*x(5)**3 + 10*x(4)**3*x(5)**3))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(2,3) = (4*(x(2) - x(3))**2*(- 100*x(1)**2*x(2)**4 - 10*x(1)**2*x(2)**3*x(3) + 205*x(1)**2*x(2)**3*x(4) + 205*x(1)**2*x(2)**3*x(5) + 190*x(1)**2*x(2)**2*x(3)**2 - 175*x(1)**2*x(2)**2*x(3)*x(4) - 175*x(1)**2*x(2)**2*x(3)*x(5) - 110*x(1)**2*x(2)**2*x(4)**2 - 220*x(1)**2*x(2)**2*x(4)*x(5) - 110*x(1)**2*x(2)**2*x(5)**2 - 10*x(1)**2*x(2)*x(3)**3 - 175*x(1)**2*x(2)*x(3)**2*x(4) - 175*x(1)**2*x(2)*x(3)**2*x(5) + 180*x(1)**2*x(2)*x(3)*x(4)**2 + 340*x(1)**2*x(2)*x(3)*x(4)*x(5) + 180*x(1)**2*x(2)*x(3)*x(5)**2 + 5*x(1)**2*x(2)*x(4)**3 + 25*x(1)**2*x(2)*x(4)**2*x(5) + 25*x(1)**2*x(2)*x(4)*x(5)**2 + 5*x(1)**2*x(2)*x(5)**3 - 100*x(1)**2*x(3)**4 + 205*x(1)**2*x(3)**3*x(4) + 205*x(1)**2*x(3)**3*x(5) - 110*x(1)**2*x(3)**2*x(4)**2 - 220*x(1)**2*x(3)**2*x(4)*x(5) - 110*x(1)**2*x(3)**2*x(5)**2 + 5*x(1)**2*x(3)*x(4)**3 + 25*x(1)**2*x(3)*x(4)**2*x(5) + 25*x(1)**2*x(3)*x(4)*x(5)**2 + 5*x(1)**2*x(3)*x(5)**3 - 10*x(1)**2*x(4)**3*x(5) - 10*x(1)**2*x(4)**2*x(5)**2 - 10*x(1)**2*x(4)*x(5)**3 + 1662*x(1)*x(2)**5 - 4586*x(1)*x(2)**4*x(3) - 1762*x(1)*x(2)**4*x(4) - 1762*x(1)*x(2)**4*x(5) + 2934*x(1)*x(2)**3*x(3)**2 + 6248*x(1)*x(2)**3*x(3)*x(4) + 6248*x(1)*x(2)**3*x(3)*x(5) + 95*x(1)*x(2)**3*x(4)**2 + 200*x(1)*x(2)**3*x(4)*x(5) + 95*x(1)*x(2)**3*x(5)**2 + 2934*x(1)*x(2)**2*x(3)**3 - 8992*x(1)*x(2)**2*x(3)**2*x(4) - 8992*x(1)*x(2)**2*x(3)**2*x(5) - 105*x(1)*x(2)**2*x(3)*x(4)**2 - 200*x(1)*x(2)**2*x(3)*x(4)*x(5) - 105*x(1)*x(2)**2*x(3)*x(5)**2 + 10*x(1)*x(2)**2*x(4)**3 + 10*x(1)*x(2)**2*x(4)**2*x(5) + 10*x(1)*x(2)**2*x(4)*x(5)**2 + 10*x(1)*x(2)**2*x(5)**3 - 4586*x(1)*x(2)*x(3)**4 + 6248*x(1)*x(2)*x(3)**3*x(4) + 6248*x(1)*x(2)*x(3)**3*x(5) - 105*x(1)*x(2)*x(3)**2*x(4)**2 - 200*x(1)*x(2)*x(3)**2*x(4)*x(5) - 105*x(1)*x(2)*x(3)**2*x(5)**2 + 10*x(1)*x(2)*x(3)*x(4)**3 + 30*x(1)*x(2)*x(3)*x(4)**2*x(5) + 30*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 10*x(1)*x(2)*x(3)*x(5)**3 - 5*x(1)*x(2)*x(4)**4 - 20*x(1)*x(2)*x(4)**3*x(5) - 20*x(1)*x(2)*x(4)**2*x(5)**2 - 20*x(1)*x(2)*x(4)*x(5)**3 - 5*x(1)*x(2)*x(5)**4 + 1662*x(1)*x(3)**5 - 1762*x(1)*x(3)**4*x(4) - 1762*x(1)*x(3)**4*x(5) + 95*x(1)*x(3)**3*x(4)**2 + 200*x(1)*x(3)**3*x(4)*x(5) + 95*x(1)*x(3)**3*x(5)**2 + 10*x(1)*x(3)**2*x(4)**3 + 10*x(1)*x(3)**2*x(4)**2*x(5) + 10*x(1)*x(3)**2*x(4)*x(5)**2 + 10*x(1)*x(3)**2*x(5)**3 - 5*x(1)*x(3)*x(4)**4 - 20*x(1)*x(3)*x(4)**3*x(5) - 20*x(1)*x(3)*x(4)**2*x(5)**2 - 20*x(1)*x(3)*x(4)*x(5)**3 - 5*x(1)*x(3)*x(5)**4 + 10*x(1)*x(4)**4*x(5) + 10*x(1)*x(4)**3*x(5)**2 + 10*x(1)*x(4)**2*x(5)**3 + 10*x(1)*x(4)*x(5)**4 + 1662*x(2)**5*x(3) - 1662*x(2)**5*x(4) - 1662*x(2)**5*x(5) - 6248*x(2)**4*x(3)**2 + 4386*x(2)**4*x(3)*x(4) + 4386*x(2)**4*x(3)*x(5) + 1862*x(2)**4*x(4)**2 + 1962*x(2)**4*x(4)*x(5) + 1862*x(2)**4*x(5)**2 + 9182*x(2)**3*x(3)**3 - 2744*x(2)**3*x(3)**2*x(4) - 2744*x(2)**3*x(3)**2*x(5) - 6138*x(2)**3*x(3)*x(4)**2 - 6028*x(2)**3*x(3)*x(4)*x(5) - 6138*x(2)**3*x(3)*x(5)**2 - 300*x(2)**3*x(4)**3 - 505*x(2)**3*x(4)**2*x(5) - 505*x(2)**3*x(4)*x(5)**2 - 300*x(2)**3*x(5)**3 - 6248*x(2)**2*x(3)**4 - 2744*x(2)**2*x(3)**3*x(4) - 2744*x(2)**2*x(3)**3*x(5) + 8612*x(2)**2*x(3)**2*x(4)**2 + 8232*x(2)**2*x(3)**2*x(4)*x(5) + 8612*x(2)**2*x(3)**2*x(5)**2 + 280*x(2)**2*x(3)*x(4)**3 + 455*x(2)**2*x(3)*x(4)**2*x(5) + 455*x(2)**2*x(3)*x(4)*x(5)**2 + 280*x(2)**2*x(3)*x(5)**3 + 100*x(2)**2*x(4)**4 + 210*x(2)**2*x(4)**3*x(5) + 210*x(2)**2*x(4)**2*x(5)**2 + 210*x(2)**2*x(4)*x(5)**3 + 100*x(2)**2*x(5)**4 + 1662*x(2)*x(3)**5 + 4386*x(2)*x(3)**4*x(4) + 4386*x(2)*x(3)**4*x(5) - 6138*x(2)*x(3)**3*x(4)**2 - 6028*x(2)*x(3)**3*x(4)*x(5) - 6138*x(2)*x(3)**3*x(5)**2 + 280*x(2)*x(3)**2*x(4)**3 + 455*x(2)*x(3)**2*x(4)**2*x(5) + 455*x(2)*x(3)**2*x(4)*x(5)**2 + 280*x(2)*x(3)**2*x(5)**3 - 190*x(2)*x(3)*x(4)**4 - 370*x(2)*x(3)*x(4)**3*x(5) - 370*x(2)*x(3)*x(4)**2*x(5)**2 - 370*x(2)*x(3)*x(4)*x(5)**3 - 190*x(2)*x(3)*x(5)**4 - 5*x(2)*x(4)**4*x(5) - 5*x(2)*x(4)**3*x(5)**2 - 5*x(2)*x(4)**2*x(5)**3 - 5*x(2)*x(4)*x(5)**4 - 1662*x(3)**5*x(4) - 1662*x(3)**5*x(5) + 1862*x(3)**4*x(4)**2 + 1962*x(3)**4*x(4)*x(5) + 1862*x(3)**4*x(5)**2 - 300*x(3)**3*x(4)**3 - 505*x(3)**3*x(4)**2*x(5) - 505*x(3)**3*x(4)*x(5)**2 - 300*x(3)**3*x(5)**3 + 100*x(3)**2*x(4)**4 + 210*x(3)**2*x(4)**3*x(5) + 210*x(3)**2*x(4)**2*x(5)**2 + 210*x(3)**2*x(4)*x(5)**3 + 100*x(3)**2*x(5)**4 - 5*x(3)*x(4)**4*x(5) - 5*x(3)*x(4)**3*x(5)**2 - 5*x(3)*x(4)**2*x(5)**3 - 5*x(3)*x(4)*x(5)**4))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(2,4) = -(4*(x(2) - x(3))**2*(9182*x(2)**2*x(3)**2 + 100*x(2)**2*x(4)**2 + 100*x(3)**2*x(4)**2 - 100*x(1)*x(2)**3 - 100*x(1)*x(3)**3 - 6248*x(2)*x(3)**3 - 6248*x(2)**3*x(3) - 200*x(2)**3*x(4) - 100*x(2)**3*x(5) - 200*x(3)**3*x(4) - 100*x(3)**3*x(5) + 1662*x(2)**4 + 1662*x(3)**4 + 95*x(1)*x(2)*x(3)**2 + 95*x(1)*x(2)**2*x(3) - 5*x(1)*x(2)*x(4)**2 + 105*x(1)*x(2)**2*x(4) - 5*x(1)*x(3)*x(4)**2 + 100*x(1)*x(2)**2*x(5) + 105*x(1)*x(3)**2*x(4) + 100*x(1)*x(3)**2*x(5) - 190*x(2)*x(3)*x(4)**2 + 190*x(2)*x(3)**2*x(4) + 190*x(2)**2*x(3)*x(4) + 10*x(1)*x(4)**2*x(5) + 95*x(2)*x(3)**2*x(5) + 95*x(2)**2*x(3)*x(5) - 5*x(2)*x(4)**2*x(5) + 105*x(2)**2*x(4)*x(5) - 5*x(3)*x(4)**2*x(5) + 105*x(3)**2*x(4)*x(5) - 190*x(1)*x(2)*x(3)*x(4) - 190*x(1)*x(2)*x(3)*x(5) - 10*x(1)*x(2)*x(4)*x(5) - 10*x(1)*x(3)*x(4)*x(5) - 190*x(2)*x(3)*x(4)*x(5)))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(2,5) = (4*(x(2) - x(3))**2*(60*x(1)**6*x(2)**2 - 90*x(1)**6*x(2)*x(3) - 15*x(1)**6*x(2)*x(4) - 15*x(1)**6*x(2)*x(5) + 50*x(1)**6*x(3)**2 - 5*x(1)**6*x(3)*x(4) - 5*x(1)**6*x(3)*x(5) + 5*x(1)**6*x(4)**2 + 10*x(1)**6*x(4)*x(5) + 5*x(1)**6*x(5)**2 - 105*x(1)**5*x(2)**3 - 25*x(1)**5*x(2)**2*x(3) - 10*x(1)**5*x(2)**2*x(4) - 10*x(1)**5*x(2)**2*x(5) + 270*x(1)**5*x(2)*x(3)**2 + 25*x(1)**5*x(2)*x(3)*x(4) + 25*x(1)**5*x(2)*x(3)*x(5) + 25*x(1)**5*x(2)*x(4)**2 + 35*x(1)**5*x(2)*x(4)*x(5) + 25*x(1)**5*x(2)*x(5)**2 - 200*x(1)**5*x(3)**3 + 15*x(1)**5*x(3)**2*x(4) + 15*x(1)**5*x(3)**2*x(5) - 5*x(1)**5*x(3)*x(4)**2 - 15*x(1)**5*x(3)*x(4)*x(5) - 5*x(1)**5*x(3)*x(5)**2 - 10*x(1)**5*x(4)**3 - 20*x(1)**5*x(4)**2*x(5) - 20*x(1)**5*x(4)*x(5)**2 - 10*x(1)**5*x(5)**3 + 726*x(1)**4*x(2)**4 - 2819*x(1)**4*x(2)**3*x(3) + 220*x(1)**4*x(2)**3*x(4) + 220*x(1)**4*x(2)**3*x(5) + 4551*x(1)**4*x(2)**2*x(3)**2 - 260*x(1)**4*x(2)**2*x(3)*x(4) - 260*x(1)**4*x(2)**2*x(3)*x(5) - 115*x(1)**4*x(2)**2*x(4)**2 - 120*x(1)**4*x(2)**2*x(4)*x(5) - 115*x(1)**4*x(2)**2*x(5)**2 - 3494*x(1)**4*x(2)*x(3)**3 + 15*x(1)**4*x(2)*x(3)**2*x(4) + 15*x(1)**4*x(2)*x(3)**2*x(5) + 130*x(1)**4*x(2)*x(3)*x(4)**2 + 105*x(1)**4*x(2)*x(3)*x(4)*x(5) + 130*x(1)**4*x(2)*x(3)*x(5)**2 - 5*x(1)**4*x(2)*x(4)**3 - 10*x(1)**4*x(2)*x(4)**2*x(5) - 10*x(1)**4*x(2)*x(4)*x(5)**2 - 5*x(1)**4*x(2)*x(5)**3 + 1081*x(1)**4*x(3)**4 + 85*x(1)**4*x(3)**3*x(4) + 85*x(1)**4*x(3)**3*x(5) - 115*x(1)**4*x(3)**2*x(4)**2 - 115*x(1)**4*x(3)**2*x(4)*x(5) - 115*x(1)**4*x(3)**2*x(5)**2 + 25*x(1)**4*x(3)*x(4)**3 + 50*x(1)**4*x(3)*x(4)**2*x(5) + 50*x(1)**4*x(3)*x(4)*x(5)**2 + 25*x(1)**4*x(3)*x(5)**3 + 5*x(1)**4*x(4)**4 + 10*x(1)**4*x(4)**3*x(5) + 15*x(1)**4*x(4)**2*x(5)**2 + 10*x(1)**4*x(4)*x(5)**3 + 5*x(1)**4*x(5)**4 + 1557*x(1)**3*x(2)**5 - 7605*x(1)**3*x(2)**4*x(3) - 1542*x(1)**3*x(2)**4*x(4) - 1542*x(1)**3*x(2)**4*x(5) + 15225*x(1)**3*x(2)**3*x(3)**2 + 5623*x(1)**3*x(2)**3*x(3)*x(4) + 5623*x(1)**3*x(2)**3*x(3)*x(5) - 30*x(1)**3*x(2)**3*x(4)**2 - 275*x(1)**3*x(2)**3*x(4)*x(5) - 30*x(1)**3*x(2)**3*x(5)**2 - 15615*x(1)**3*x(2)**2*x(3)**3 - 8517*x(1)**3*x(2)**2*x(3)**2*x(4) - 8517*x(1)**3*x(2)**2*x(3)**2*x(5) + 310*x(1)**3*x(2)**2*x(3)*x(4)**2 + 585*x(1)**3*x(2)**2*x(3)*x(4)*x(5) + 310*x(1)**3*x(2)**2*x(3)*x(5)**2 + 20*x(1)**3*x(2)**2*x(4)**3 + 180*x(1)**3*x(2)**2*x(4)**2*x(5) + 180*x(1)**3*x(2)**2*x(4)*x(5)**2 + 20*x(1)**3*x(2)**2*x(5)**3 + 8200*x(1)**3*x(2)*x(3)**4 + 6203*x(1)**3*x(2)*x(3)**3*x(4) + 6203*x(1)**3*x(2)*x(3)**3*x(5) - 545*x(1)**3*x(2)*x(3)**2*x(4)**2 - 545*x(1)**3*x(2)*x(3)**2*x(4)*x(5) - 545*x(1)**3*x(2)*x(3)**2*x(5)**2 + 25*x(1)**3*x(2)*x(3)*x(4)**3 - 125*x(1)**3*x(2)*x(3)*x(4)**2*x(5) - 125*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 25*x(1)**3*x(2)*x(3)*x(5)**3 - 5*x(1)**3*x(2)*x(4)**4 - 25*x(1)**3*x(2)*x(4)**3*x(5) - 60*x(1)**3*x(2)*x(4)**2*x(5)**2 - 25*x(1)**3*x(2)*x(4)*x(5)**3 - 5*x(1)**3*x(2)*x(5)**4 - 1762*x(1)**3*x(3)**5 - 1857*x(1)**3*x(3)**4*x(4) - 1857*x(1)**3*x(3)**4*x(5) + 325*x(1)**3*x(3)**3*x(4)**2 + 235*x(1)**3*x(3)**3*x(4)*x(5) + 325*x(1)**3*x(3)**3*x(5)**2 - 15*x(1)**3*x(3)**2*x(4)**3 + 75*x(1)**3*x(3)**2*x(4)**2*x(5) + 75*x(1)**3*x(3)**2*x(4)*x(5)**2 - 15*x(1)**3*x(3)**2*x(5)**3 - 15*x(1)**3*x(3)*x(4)**4 - 35*x(1)**3*x(3)*x(4)**3*x(5) - 60*x(1)**3*x(3)*x(4)**2*x(5)**2 - 35*x(1)**3*x(3)*x(4)*x(5)**3 - 15*x(1)**3*x(3)*x(5)**4 + 10*x(1)**3*x(4)**3*x(5)**2 + 10*x(1)**3*x(4)**2*x(5)**3 + 2493*x(1)**2*x(2)**6 - 12591*x(1)**2*x(2)**5*x(3) - 3519*x(1)**2*x(2)**5*x(4) - 3519*x(1)**2*x(2)**5*x(5) + 26900*x(1)**2*x(2)**4*x(3)**2 + 15985*x(1)**2*x(2)**4*x(3)*x(4) + 15985*x(1)**2*x(2)**4*x(3)*x(5) + 1166*x(1)**2*x(2)**4*x(4)**2 + 3904*x(1)**2*x(2)**4*x(4)*x(5) + 1166*x(1)**2*x(2)**4*x(5)**2 - 30855*x(1)**2*x(2)**3*x(3)**3 - 30355*x(1)**2*x(2)**3*x(3)**2*x(4) - 30355*x(1)**2*x(2)**3*x(3)**2*x(5) - 3564*x(1)**2*x(2)**3*x(3)*x(4)**2 - 12971*x(1)**2*x(2)**3*x(3)*x(4)*x(5) - 3564*x(1)**2*x(2)**3*x(3)*x(5)**2 - 185*x(1)**2*x(2)**3*x(4)**3 - 455*x(1)**2*x(2)**3*x(4)**2*x(5) - 455*x(1)**2*x(2)**3*x(4)*x(5)**2 - 185*x(1)**2*x(2)**3*x(5)**3 + 19770*x(1)**2*x(2)**2*x(3)**4 + 30165*x(1)**2*x(2)**2*x(3)**3*x(4) + 30165*x(1)**2*x(2)**2*x(3)**3*x(5) + 4251*x(1)**2*x(2)**2*x(3)**2*x(4)**2 + 17619*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) + 4251*x(1)**2*x(2)**2*x(3)**2*x(5)**2 + 260*x(1)**2*x(2)**2*x(3)*x(4)**3 + 480*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) + 480*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 + 260*x(1)**2*x(2)**2*x(3)*x(5)**3 + 45*x(1)**2*x(2)**2*x(4)**4 + 55*x(1)**2*x(2)**2*x(4)**3*x(5) + 90*x(1)**2*x(2)**2*x(4)**2*x(5)**2 + 55*x(1)**2*x(2)**2*x(4)*x(5)**3 + 45*x(1)**2*x(2)**2*x(5)**4 - 6548*x(1)**2*x(2)*x(3)**5 - 15700*x(1)**2*x(2)*x(3)**4*x(4) - 15700*x(1)**2*x(2)*x(3)**4*x(5) - 2329*x(1)**2*x(2)*x(3)**3*x(4)**2 - 11481*x(1)**2*x(2)*x(3)**3*x(4)*x(5) - 2329*x(1)**2*x(2)*x(3)**3*x(5)**2 - 30*x(1)**2*x(2)*x(3)**2*x(4)**3 + 210*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) + 210*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 - 30*x(1)**2*x(2)*x(3)**2*x(5)**3 - 90*x(1)**2*x(2)*x(3)*x(4)**4 - 175*x(1)**2*x(2)*x(3)*x(4)**3*x(5) - 240*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 - 175*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 - 90*x(1)**2*x(2)*x(3)*x(5)**4 + 15*x(1)**2*x(2)*x(4)**4*x(5) + 40*x(1)**2*x(2)*x(4)**3*x(5)**2 + 40*x(1)**2*x(2)*x(4)**2*x(5)**3 + 15*x(1)**2*x(2)*x(4)*x(5)**4 + 831*x(1)**2*x(3)**6 + 3424*x(1)**2*x(3)**5*x(4) + 3424*x(1)**2*x(3)**5*x(5) + 521*x(1)**2*x(3)**4*x(4)**2 + 3109*x(1)**2*x(3)**4*x(4)*x(5) + 521*x(1)**2*x(3)**4*x(5)**2 - 105*x(1)**2*x(3)**3*x(4)**3 - 415*x(1)**2*x(3)**3*x(4)**2*x(5) - 415*x(1)**2*x(3)**3*x(4)*x(5)**2 - 105*x(1)**2*x(3)**3*x(5)**3 + 65*x(1)**2*x(3)**2*x(4)**4 + 130*x(1)**2*x(3)**2*x(4)**3*x(5) + 210*x(1)**2*x(3)**2*x(4)**2*x(5)**2 + 130*x(1)**2*x(3)**2*x(4)*x(5)**3 + 65*x(1)**2*x(3)**2*x(5)**4 + 5*x(1)**2*x(3)*x(4)**4*x(5) + 5*x(1)**2*x(3)*x(4)*x(5)**4 - 10*x(1)**2*x(4)**4*x(5)**2 - 10*x(1)**2*x(4)**3*x(5)**3 - 10*x(1)**2*x(4)**2*x(5)**4 + 1662*x(1)*x(2)**7 - 9572*x(1)*x(2)**6*x(3) - 3524*x(1)*x(2)**6*x(4) - 3524*x(1)*x(2)**6*x(5) + 23340*x(1)*x(2)**5*x(3)**2 + 17967*x(1)*x(2)**5*x(3)*x(4) + 17967*x(1)*x(2)**5*x(3)*x(5) + 2162*x(1)*x(2)**5*x(4)**2 + 5891*x(1)*x(2)**5*x(4)*x(5) + 2162*x(1)*x(2)**5*x(5)**2 - 30860*x(1)*x(2)**4*x(3)**3 - 38960*x(1)*x(2)**4*x(3)**2*x(4) - 38960*x(1)*x(2)**4*x(3)**2*x(5) - 9065*x(1)*x(2)**4*x(3)*x(4)**2 - 25755*x(1)*x(2)**4*x(3)*x(4)*x(5) - 9065*x(1)*x(2)**4*x(3)*x(5)**2 - 400*x(1)*x(2)**4*x(4)**3 - 2877*x(1)*x(2)**4*x(4)**2*x(5) - 2877*x(1)*x(2)**4*x(4)*x(5)**2 - 400*x(1)*x(2)**4*x(5)**3 + 23340*x(1)*x(2)**3*x(3)**4 + 45895*x(1)*x(2)**3*x(3)**3*x(4) + 45895*x(1)*x(2)**3*x(3)**3*x(5) + 15905*x(1)*x(2)**3*x(3)**2*x(4)**2 + 47055*x(1)*x(2)**3*x(3)**2*x(4)*x(5) + 15905*x(1)*x(2)**3*x(3)**2*x(5)**2 + 955*x(1)*x(2)**3*x(3)*x(4)**3 + 8713*x(1)*x(2)**3*x(3)*x(4)**2*x(5) + 8713*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 + 955*x(1)*x(2)**3*x(3)*x(5)**3 + 100*x(1)*x(2)**3*x(4)**4 + 615*x(1)*x(2)**3*x(4)**3*x(5) + 930*x(1)*x(2)**3*x(4)**2*x(5)**2 + 615*x(1)*x(2)**3*x(4)*x(5)**3 + 100*x(1)*x(2)**3*x(5)**4 - 9572*x(1)*x(2)**2*x(3)**5 - 30850*x(1)*x(2)**2*x(3)**4*x(4) - 30850*x(1)*x(2)**2*x(3)**4*x(5) - 14745*x(1)*x(2)**2*x(3)**3*x(4)**2 - 45125*x(1)*x(2)**2*x(3)**3*x(4)*x(5) - 14745*x(1)*x(2)**2*x(3)**3*x(5)**2 - 575*x(1)*x(2)**2*x(3)**2*x(4)**3 - 10257*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) - 10257*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 - 575*x(1)*x(2)**2*x(3)**2*x(5)**3 - 285*x(1)*x(2)**2*x(3)*x(4)**4 - 1095*x(1)*x(2)**2*x(3)*x(4)**3*x(5) - 1650*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 - 1095*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 - 285*x(1)*x(2)**2*x(3)*x(5)**4 - 105*x(1)*x(2)**2*x(4)**4*x(5) - 220*x(1)*x(2)**2*x(4)**3*x(5)**2 - 220*x(1)*x(2)**2*x(4)**2*x(5)**3 - 105*x(1)*x(2)**2*x(4)*x(5)**4 + 1662*x(1)*x(2)*x(3)**6 + 11134*x(1)*x(2)*x(3)**5*x(4) + 11134*x(1)*x(2)*x(3)**5*x(5) + 7305*x(1)*x(2)*x(3)**4*x(4)**2 + 22820*x(1)*x(2)*x(3)**4*x(4)*x(5) + 7305*x(1)*x(2)*x(3)**4*x(5)**2 - 185*x(1)*x(2)*x(3)**3*x(4)**3 + 5483*x(1)*x(2)*x(3)**3*x(4)**2*x(5) + 5483*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 - 185*x(1)*x(2)*x(3)**3*x(5)**3 + 290*x(1)*x(2)*x(3)**2*x(4)**4 + 605*x(1)*x(2)*x(3)**2*x(4)**3*x(5) + 915*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 + 605*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 + 290*x(1)*x(2)*x(3)**2*x(5)**4 + 170*x(1)*x(2)*x(3)*x(4)**4*x(5) + 325*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 + 325*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 + 170*x(1)*x(2)*x(3)*x(4)*x(5)**4 + 5*x(1)*x(2)*x(4)**4*x(5)**2 + 5*x(1)*x(2)*x(4)**3*x(5)**3 + 5*x(1)*x(2)*x(4)**2*x(5)**4 - 1662*x(1)*x(3)**6*x(4) - 1662*x(1)*x(3)**6*x(5) - 1562*x(1)*x(3)**5*x(4)**2 - 4886*x(1)*x(3)**5*x(4)*x(5) - 1562*x(1)*x(3)**5*x(5)**2 + 205*x(1)*x(3)**4*x(4)**3 - 1152*x(1)*x(3)**4*x(4)**2*x(5) - 1152*x(1)*x(3)**4*x(4)*x(5)**2 + 205*x(1)*x(3)**4*x(5)**3 - 105*x(1)*x(3)**3*x(4)**4 - 5*x(1)*x(3)**3*x(4)**3*x(5) - 15*x(1)*x(3)**3*x(4)**2*x(5)**2 - 5*x(1)*x(3)**3*x(4)*x(5)**3 - 105*x(1)*x(3)**3*x(5)**4 - 105*x(1)*x(3)**2*x(4)**4*x(5) - 215*x(1)*x(3)**2*x(4)**3*x(5)**2 - 215*x(1)*x(3)**2*x(4)**2*x(5)**3 - 105*x(1)*x(3)**2*x(4)*x(5)**4 + 15*x(1)*x(3)*x(4)**4*x(5)**2 + 15*x(1)*x(3)*x(4)**3*x(5)**3 + 15*x(1)*x(3)*x(4)**2*x(5)**4 + 831*x(2)**8 - 4786*x(2)**7*x(3) - 1762*x(2)**7*x(4) - 1762*x(2)**7*x(5) + 11670*x(2)**6*x(3)**2 + 9867*x(2)**6*x(3)*x(4) + 9867*x(2)**6*x(3)*x(5) + 1081*x(2)**6*x(4)**2 + 3829*x(2)**6*x(4)*x(5) + 1081*x(2)**6*x(5)**2 - 15430*x(2)**5*x(3)**3 - 23535*x(2)**5*x(3)**2*x(4) - 23535*x(2)**5*x(3)**2*x(5) - 5571*x(2)**5*x(3)*x(4)**2 - 18957*x(2)**5*x(3)*x(4)*x(5) - 5571*x(2)**5*x(3)*x(5)**2 - 200*x(2)**5*x(4)**3 - 2477*x(2)**5*x(4)**2*x(5) - 2477*x(2)**5*x(4)*x(5)**2 - 200*x(2)**5*x(5)**3 + 11670*x(2)**4*x(3)**4 + 30665*x(2)**4*x(3)**3*x(4) + 30665*x(2)**4*x(3)**3*x(5) + 12350*x(2)**4*x(3)**2*x(4)**2 + 39940*x(2)**4*x(3)**2*x(4)*x(5) + 12350*x(2)**4*x(3)**2*x(5)**2 + 685*x(2)**4*x(3)*x(4)**3 + 10165*x(2)**4*x(3)*x(4)**2*x(5) + 10165*x(2)**4*x(3)*x(4)*x(5)**2 + 685*x(2)**4*x(3)*x(5)**3 + 50*x(2)**4*x(4)**4 + 515*x(2)**4*x(4)**3*x(5) + 1776*x(2)**4*x(4)**2*x(5)**2 + 515*x(2)**4*x(4)*x(5)**3 + 50*x(2)**4*x(5)**4 - 4786*x(2)**3*x(3)**5 - 23045*x(2)**3*x(3)**4*x(4) - 23045*x(2)**3*x(3)**4*x(5) - 15235*x(2)**3*x(3)**3*x(4)**2 - 45905*x(2)**3*x(3)**3*x(4)*x(5) - 15235*x(2)**3*x(3)**3*x(5)**2 - 775*x(2)**3*x(3)**2*x(4)**3 - 17275*x(2)**3*x(3)**2*x(4)**2*x(5) - 17275*x(2)**3*x(3)**2*x(4)*x(5)**2 - 775*x(2)**3*x(3)**2*x(5)**3 - 195*x(2)**3*x(3)*x(4)**4 - 1365*x(2)**3*x(3)*x(4)**3*x(5) - 5364*x(2)**3*x(3)*x(4)**2*x(5)**2 - 1365*x(2)**3*x(3)*x(4)*x(5)**3 - 195*x(2)**3*x(3)*x(5)**4 - 105*x(2)**3*x(4)**4*x(5) - 445*x(2)**3*x(4)**3*x(5)**2 - 445*x(2)**3*x(4)**2*x(5)**3 - 105*x(2)**3*x(4)*x(5)**4 + 831*x(2)**2*x(3)**6 + 9472*x(2)**2*x(3)**5*x(4) + 9472*x(2)**2*x(3)**5*x(5) + 11130*x(2)**2*x(3)**4*x(4)**2 + 30365*x(2)**2*x(3)**4*x(4)*x(5) + 11130*x(2)**2*x(3)**4*x(5)**2 + 195*x(2)**2*x(3)**3*x(4)**3 + 15345*x(2)**2*x(3)**3*x(4)**2*x(5) + 15345*x(2)**2*x(3)**3*x(4)*x(5)**2 + 195*x(2)**2*x(3)**3*x(5)**3 + 290*x(2)**2*x(3)**2*x(4)**4 + 1155*x(2)**2*x(3)**2*x(4)**3*x(5) + 6291*x(2)**2*x(3)**2*x(4)**2*x(5)**2 + 1155*x(2)**2*x(3)**2*x(4)*x(5)**3 + 290*x(2)**2*x(3)**2*x(5)**4 + 290*x(2)**2*x(3)*x(4)**4*x(5) + 860*x(2)**2*x(3)*x(4)**3*x(5)**2 + 860*x(2)**2*x(3)*x(4)**2*x(5)**3 + 290*x(2)**2*x(3)*x(4)*x(5)**4 + 65*x(2)**2*x(4)**4*x(5)**2 + 145*x(2)**2*x(4)**3*x(5)**3 + 65*x(2)**2*x(4)**2*x(5)**4 - 1662*x(2)*x(3)**6*x(4) - 1662*x(2)*x(3)**6*x(5) - 4586*x(2)*x(3)**5*x(4)**2 - 10934*x(2)*x(3)**5*x(4)*x(5) - 4586*x(2)*x(3)**5*x(5)**2 + 195*x(2)*x(3)**4*x(4)**3 - 7220*x(2)*x(3)**4*x(4)**2*x(5) - 7220*x(2)*x(3)**4*x(4)*x(5)**2 + 195*x(2)*x(3)**4*x(5)**3 - 195*x(2)*x(3)**3*x(4)**4 - 205*x(2)*x(3)**3*x(4)**3*x(5) - 3339*x(2)*x(3)**3*x(4)**2*x(5)**2 - 205*x(2)*x(3)**3*x(4)*x(5)**3 - 195*x(2)*x(3)**3*x(5)**4 - 285*x(2)*x(3)**2*x(4)**4*x(5) - 580*x(2)*x(3)**2*x(4)**3*x(5)**2 - 580*x(2)*x(3)**2*x(4)**2*x(5)**3 - 285*x(2)*x(3)**2*x(4)*x(5)**4 - 90*x(2)*x(3)*x(4)**4*x(5)**2 - 175*x(2)*x(3)*x(4)**3*x(5)**3 - 90*x(2)*x(3)*x(4)**2*x(5)**4 - 15*x(2)*x(4)**4*x(5)**3 - 15*x(2)*x(4)**3*x(5)**4 + 831*x(3)**6*x(4)**2 + 1662*x(3)**6*x(4)*x(5) + 831*x(3)**6*x(5)**2 - 100*x(3)**5*x(4)**3 + 1462*x(3)**5*x(4)**2*x(5) + 1462*x(3)**5*x(4)*x(5)**2 - 100*x(3)**5*x(5)**3 + 50*x(3)**4*x(4)**4 - 100*x(3)**4*x(4)**3*x(5) + 681*x(3)**4*x(4)**2*x(5)**2 - 100*x(3)**4*x(4)*x(5)**3 + 50*x(3)**4*x(5)**4 + 100*x(3)**3*x(4)**4*x(5) + 105*x(3)**3*x(4)**3*x(5)**2 + 105*x(3)**3*x(4)**2*x(5)**3 + 100*x(3)**3*x(4)*x(5)**4 + 45*x(3)**2*x(4)**4*x(5)**2 + 100*x(3)**2*x(4)**3*x(5)**3 + 45*x(3)**2*x(4)**2*x(5)**4 - 5*x(3)*x(4)**4*x(5)**3 - 5*x(3)*x(4)**3*x(5)**4 + 5*x(4)**4*x(5)**4))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2)

                            Cb(2,6) = (4*(x(2) - x(3))**2*(- 105*x(1)**5*x(2)**3 + 65*x(1)**5*x(2)**2*x(3) + 125*x(1)**5*x(2)**2*x(4) + 125*x(1)**5*x(2)**2*x(5) + 80*x(1)**5*x(2)*x(3)**2 - 145*x(1)**5*x(2)*x(3)*x(4) - 145*x(1)**5*x(2)*x(3)*x(5) - 30*x(1)**5*x(2)*x(4)**2 - 45*x(1)**5*x(2)*x(4)*x(5) - 30*x(1)**5*x(2)*x(5)**2 - 100*x(1)**5*x(3)**3 + 110*x(1)**5*x(3)**2*x(4) + 110*x(1)**5*x(3)**2*x(5) - 20*x(1)**5*x(3)*x(4)**2 - 35*x(1)**5*x(3)*x(4)*x(5) - 20*x(1)**5*x(3)*x(5)**2 + 10*x(1)**5*x(4)**3 + 20*x(1)**5*x(4)**2*x(5) + 20*x(1)**5*x(4)*x(5)**2 + 10*x(1)**5*x(5)**3 + 205*x(1)**4*x(2)**4 + 125*x(1)**4*x(2)**3*x(3) - 210*x(1)**4*x(2)**3*x(4) - 210*x(1)**4*x(2)**3*x(5) - 440*x(1)**4*x(2)**2*x(3)**2 + 90*x(1)**4*x(2)**2*x(3)*x(4) + 90*x(1)**4*x(2)**2*x(3)*x(5) - 25*x(1)**4*x(2)**2*x(4)**2 - 35*x(1)**4*x(2)**2*x(4)*x(5) - 25*x(1)**4*x(2)**2*x(5)**2 - 60*x(1)**4*x(2)*x(3)**3 + 330*x(1)**4*x(2)*x(3)**2*x(4) + 330*x(1)**4*x(2)*x(3)**2*x(5) - 35*x(1)**4*x(2)*x(3)*x(4)**2 - 45*x(1)**4*x(2)*x(3)*x(4)*x(5) - 35*x(1)**4*x(2)*x(3)*x(5)**2 + 50*x(1)**4*x(2)*x(4)**3 + 85*x(1)**4*x(2)*x(4)**2*x(5) + 85*x(1)**4*x(2)*x(4)*x(5)**2 + 50*x(1)**4*x(2)*x(5)**3 + 300*x(1)**4*x(3)**4 - 320*x(1)**4*x(3)**3*x(4) - 320*x(1)**4*x(3)**3*x(5) + 20*x(1)**4*x(3)**2*x(4)**2 + 40*x(1)**4*x(3)**2*x(4)*x(5) + 20*x(1)**4*x(3)**2*x(5)**2 + 20*x(1)**4*x(3)*x(4)**3 + 35*x(1)**4*x(3)*x(4)**2*x(5) + 35*x(1)**4*x(3)*x(4)*x(5)**2 + 20*x(1)**4*x(3)*x(5)**3 - 20*x(1)**4*x(4)**4 - 40*x(1)**4*x(4)**3*x(5) - 50*x(1)**4*x(4)**2*x(5)**2 - 40*x(1)**4*x(4)*x(5)**3 - 20*x(1)**4*x(5)**4 - 1562*x(1)**3*x(2)**5 + 4496*x(1)**3*x(2)**4*x(3) + 1247*x(1)**3*x(2)**4*x(4) + 1247*x(1)**3*x(2)**4*x(5) - 3239*x(1)**3*x(2)**3*x(3)**2 - 6003*x(1)**3*x(2)**3*x(3)*x(4) - 6003*x(1)**3*x(2)**3*x(3)*x(5) + 545*x(1)**3*x(2)**3*x(4)**2 + 765*x(1)**3*x(2)**3*x(4)*x(5) + 545*x(1)**3*x(2)**3*x(5)**2 - 2639*x(1)**3*x(2)**2*x(3)**3 + 9697*x(1)**3*x(2)**2*x(3)**2*x(4) + 9697*x(1)**3*x(2)**2*x(3)**2*x(5) - 495*x(1)**3*x(2)**2*x(3)*x(4)**2 - 755*x(1)**3*x(2)**2*x(3)*x(4)*x(5) - 495*x(1)**3*x(2)**2*x(3)*x(5)**2 - 230*x(1)**3*x(2)**2*x(4)**3 - 350*x(1)**3*x(2)**2*x(4)**2*x(5) - 350*x(1)**3*x(2)**2*x(4)*x(5)**2 - 230*x(1)**3*x(2)**2*x(5)**3 + 4766*x(1)**3*x(2)*x(3)**4 - 6773*x(1)**3*x(2)*x(3)**3*x(4) - 6773*x(1)**3*x(2)*x(3)**3*x(5) - 150*x(1)**3*x(2)*x(3)**2*x(4)**2 - 95*x(1)**3*x(2)*x(3)**2*x(4)*x(5) - 150*x(1)**3*x(2)*x(3)**2*x(5)**2 + 320*x(1)**3*x(2)*x(3)*x(4)**3 + 470*x(1)**3*x(2)*x(3)*x(4)**2*x(5) + 470*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 320*x(1)**3*x(2)*x(3)*x(5)**3 - 10*x(1)**3*x(2)*x(4)**4 - 20*x(1)**3*x(2)*x(4)**3*x(5) - 25*x(1)**3*x(2)*x(4)**2*x(5)**2 - 20*x(1)**3*x(2)*x(4)*x(5)**3 - 10*x(1)**3*x(2)*x(5)**4 - 1862*x(1)**3*x(3)**5 + 1672*x(1)**3*x(3)**4*x(4) + 1672*x(1)**3*x(3)**4*x(5) + 420*x(1)**3*x(3)**3*x(4)**2 + 525*x(1)**3*x(3)**3*x(4)*x(5) + 420*x(1)**3*x(3)**3*x(5)**2 - 260*x(1)**3*x(3)**2*x(4)**3 - 410*x(1)**3*x(3)**2*x(4)**2*x(5) - 410*x(1)**3*x(3)**2*x(4)*x(5)**2 - 260*x(1)**3*x(3)**2*x(5)**3 + 20*x(1)**3*x(3)*x(4)**4 + 40*x(1)**3*x(3)*x(4)**3*x(5) + 45*x(1)**3*x(3)*x(4)**2*x(5)**2 + 40*x(1)**3*x(3)*x(4)*x(5)**3 + 20*x(1)**3*x(3)*x(5)**4 + 10*x(1)**3*x(4)**5 + 20*x(1)**3*x(4)**4*x(5) + 30*x(1)**3*x(4)**3*x(5)**2 + 30*x(1)**3*x(4)**2*x(5)**3 + 20*x(1)**3*x(4)*x(5)**4 + 10*x(1)**3*x(5)**5 - 1562*x(1)**2*x(2)**6 + 4496*x(1)**2*x(2)**5*x(3) + 4781*x(1)**2*x(2)**5*x(4) + 4781*x(1)**2*x(2)**5*x(5) - 1577*x(1)**2*x(2)**4*x(3)**2 - 16407*x(1)**2*x(2)**4*x(3)*x(4) - 16407*x(1)**2*x(2)**4*x(3)*x(5) - 3299*x(1)**2*x(2)**4*x(4)**2 - 4641*x(1)**2*x(2)**4*x(4)*x(5) - 3299*x(1)**2*x(2)**4*x(5)**2 - 8887*x(1)**2*x(2)**3*x(3)**3 + 21343*x(1)**2*x(2)**3*x(3)**2*x(4) + 21343*x(1)**2*x(2)**3*x(3)**2*x(5) + 11711*x(1)**2*x(2)**3*x(3)*x(4)**2 + 17529*x(1)**2*x(2)**3*x(3)*x(4)*x(5) + 11711*x(1)**2*x(2)**3*x(3)*x(5)**2 + 55*x(1)**2*x(2)**3*x(4)**3 - 315*x(1)**2*x(2)**3*x(4)**2*x(5) - 315*x(1)**2*x(2)**3*x(4)*x(5)**2 + 55*x(1)**2*x(2)**3*x(5)**3 + 13948*x(1)**2*x(2)**2*x(3)**4 - 10607*x(1)**2*x(2)**2*x(3)**3*x(4) - 10607*x(1)**2*x(2)**2*x(3)**3*x(5) - 17414*x(1)**2*x(2)**2*x(3)**2*x(4)**2 - 26471*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) - 17414*x(1)**2*x(2)**2*x(3)**2*x(5)**2 + 175*x(1)**2*x(2)**2*x(3)*x(4)**3 + 655*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) + 655*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 + 175*x(1)**2*x(2)**2*x(3)*x(5)**3 + 35*x(1)**2*x(2)**2*x(4)**4 + 210*x(1)**2*x(2)**2*x(4)**3*x(5) + 355*x(1)**2*x(2)**2*x(4)**2*x(5)**2 + 210*x(1)**2*x(2)**2*x(4)*x(5)**3 + 35*x(1)**2*x(2)**2*x(5)**4 - 8110*x(1)**2*x(2)*x(3)**5 - 822*x(1)**2*x(2)*x(3)**4*x(4) - 822*x(1)**2*x(2)*x(3)**4*x(5) + 12826*x(1)**2*x(2)*x(3)**3*x(4)**2 + 19169*x(1)**2*x(2)*x(3)**3*x(4)*x(5) + 12826*x(1)**2*x(2)*x(3)**3*x(5)**2 - 710*x(1)**2*x(2)*x(3)**2*x(4)**3 - 1070*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) - 1070*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 - 710*x(1)**2*x(2)*x(3)**2*x(5)**3 + 45*x(1)**2*x(2)*x(3)*x(4)**4 - 70*x(1)**2*x(2)*x(3)*x(4)**3*x(5) - 185*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 - 70*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 + 45*x(1)**2*x(2)*x(3)*x(5)**4 - 10*x(1)**2*x(2)*x(4)**5 - 35*x(1)**2*x(2)*x(4)**4*x(5) - 75*x(1)**2*x(2)*x(4)**3*x(5)**2 - 75*x(1)**2*x(2)*x(4)**2*x(5)**3 - 35*x(1)**2*x(2)*x(4)*x(5)**4 - 10*x(1)**2*x(2)*x(5)**5 + 1662*x(1)**2*x(3)**6 + 1862*x(1)**2*x(3)**5*x(4) + 1862*x(1)**2*x(3)**5*x(5) - 3944*x(1)**2*x(3)**4*x(4)**2 - 5616*x(1)**2*x(3)**4*x(4)*x(5) - 3944*x(1)**2*x(3)**4*x(5)**2 + 420*x(1)**2*x(3)**3*x(4)**3 + 430*x(1)**2*x(3)**3*x(4)**2*x(5) + 430*x(1)**2*x(3)**3*x(4)*x(5)**2 + 420*x(1)**2*x(3)**3*x(5)**3 + 20*x(1)**2*x(3)**2*x(4)**4 + 150*x(1)**2*x(3)**2*x(4)**3*x(5) + 280*x(1)**2*x(3)**2*x(4)**2*x(5)**2 + 150*x(1)**2*x(3)**2*x(4)*x(5)**3 + 20*x(1)**2*x(3)**2*x(5)**4 - 20*x(1)**2*x(3)*x(4)**5 - 45*x(1)**2*x(3)*x(4)**4*x(5) - 85*x(1)**2*x(3)*x(4)**3*x(5)**2 - 85*x(1)**2*x(3)*x(4)**2*x(5)**3 - 45*x(1)**2*x(3)*x(4)*x(5)**4 - 20*x(1)**2*x(3)*x(5)**5 + 10*x(1)**2*x(4)**4*x(5)**2 + 10*x(1)**2*x(4)**3*x(5)**3 + 10*x(1)**2*x(4)**2*x(5)**4 - 1662*x(1)*x(2)**7 + 4586*x(1)*x(2)**6*x(3) + 5086*x(1)*x(2)**6*x(4) + 5086*x(1)*x(2)**6*x(5) + 390*x(1)*x(2)**5*x(3)**2 - 18644*x(1)*x(2)**5*x(3)*x(4) - 18644*x(1)*x(2)**5*x(3)*x(5) - 5381*x(1)*x(2)**5*x(4)**2 - 10672*x(1)*x(2)**5*x(4)*x(5) - 5381*x(1)*x(2)**5*x(5)**2 - 15430*x(1)*x(2)**4*x(3)**3 + 23747*x(1)*x(2)**4*x(3)**2*x(4) + 23747*x(1)*x(2)**4*x(3)**2*x(5) + 20491*x(1)*x(2)**4*x(3)*x(4)**2 + 37558*x(1)*x(2)**4*x(3)*x(4)*x(5) + 20491*x(1)*x(2)**4*x(3)*x(5)**2 + 2247*x(1)*x(2)**4*x(4)**3 + 6271*x(1)*x(2)**4*x(4)**2*x(5) + 6271*x(1)*x(2)**4*x(4)*x(5)**2 + 2247*x(1)*x(2)**4*x(5)**3 + 22950*x(1)*x(2)**3*x(3)**4 - 6153*x(1)*x(2)**3*x(3)**3*x(4) - 6153*x(1)*x(2)**3*x(3)**3*x(5) - 32854*x(1)*x(2)**3*x(3)**2*x(4)**2 - 53507*x(1)*x(2)**3*x(3)**2*x(4)*x(5) - 32854*x(1)*x(2)**3*x(3)**2*x(5)**2 - 6703*x(1)*x(2)**3*x(3)*x(4)**3 - 19569*x(1)*x(2)**3*x(3)*x(4)**2*x(5) - 19569*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 - 6703*x(1)*x(2)**3*x(3)*x(5)**3 - 385*x(1)*x(2)**3*x(4)**4 - 855*x(1)*x(2)**3*x(4)**3*x(5) - 1160*x(1)*x(2)**3*x(4)**2*x(5)**2 - 855*x(1)*x(2)**3*x(4)*x(5)**3 - 385*x(1)*x(2)**3*x(5)**4 - 14158*x(1)*x(2)**2*x(3)**5 - 12978*x(1)*x(2)**2*x(3)**4*x(4) - 12978*x(1)*x(2)**2*x(3)**4*x(5) + 27416*x(1)*x(2)**2*x(3)**3*x(4)**2 + 36753*x(1)*x(2)**2*x(3)**3*x(4)*x(5) + 27416*x(1)*x(2)**2*x(3)**3*x(5)**2 + 8447*x(1)*x(2)**2*x(3)**2*x(4)**3 + 25801*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) + 25801*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 + 8447*x(1)*x(2)**2*x(3)**2*x(5)**3 + 455*x(1)*x(2)**2*x(3)*x(4)**4 + 1045*x(1)*x(2)**2*x(3)*x(4)**3*x(5) + 1330*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 + 1045*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 + 455*x(1)*x(2)**2*x(3)*x(5)**4 + 95*x(1)*x(2)**2*x(4)**5 + 155*x(1)*x(2)**2*x(4)**4*x(5) + 240*x(1)*x(2)**2*x(4)**3*x(5)**2 + 240*x(1)*x(2)**2*x(4)**2*x(5)**3 + 155*x(1)*x(2)**2*x(4)*x(5)**4 + 95*x(1)*x(2)**2*x(5)**5 + 3324*x(1)*x(2)*x(3)**6 + 12296*x(1)*x(2)*x(3)**5*x(4) + 12296*x(1)*x(2)*x(3)**5*x(5) - 11644*x(1)*x(2)*x(3)**4*x(4)**2 - 10592*x(1)*x(2)*x(3)**4*x(4)*x(5) - 11644*x(1)*x(2)*x(3)**4*x(5)**2 - 5513*x(1)*x(2)*x(3)**3*x(4)**3 - 17369*x(1)*x(2)*x(3)**3*x(4)**2*x(5) - 17369*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 - 5513*x(1)*x(2)*x(3)**3*x(5)**3 + 160*x(1)*x(2)*x(3)**2*x(4)**4 + 425*x(1)*x(2)*x(3)**2*x(4)**3*x(5) + 685*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 + 425*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 + 160*x(1)*x(2)*x(3)**2*x(5)**4 - 185*x(1)*x(2)*x(3)*x(4)**5 - 395*x(1)*x(2)*x(3)*x(4)**4*x(5) - 610*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 - 610*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 - 395*x(1)*x(2)*x(3)*x(4)*x(5)**4 - 185*x(1)*x(2)*x(3)*x(5)**5 + 15*x(1)*x(2)*x(4)**5*x(5) + 40*x(1)*x(2)*x(4)**4*x(5)**2 + 40*x(1)*x(2)*x(4)**3*x(5)**3 + 40*x(1)*x(2)*x(4)**2*x(5)**4 + 15*x(1)*x(2)*x(4)*x(5)**5 - 3324*x(1)*x(3)**6*x(4) - 3324*x(1)*x(3)**6*x(5) + 1862*x(1)*x(3)**5*x(4)**2 + 200*x(1)*x(3)**5*x(4)*x(5) + 1862*x(1)*x(3)**5*x(5)**2 + 1672*x(1)*x(3)**4*x(4)**3 + 5206*x(1)*x(3)**4*x(4)**2*x(5) + 5206*x(1)*x(3)**4*x(4)*x(5)**2 + 1672*x(1)*x(3)**4*x(5)**3 - 320*x(1)*x(3)**3*x(4)**4 - 735*x(1)*x(3)**3*x(4)**3*x(5) - 1055*x(1)*x(3)**3*x(4)**2*x(5)**2 - 735*x(1)*x(3)**3*x(4)*x(5)**3 - 320*x(1)*x(3)**3*x(5)**4 + 110*x(1)*x(3)**2*x(4)**5 + 210*x(1)*x(3)**2*x(4)**4*x(5) + 320*x(1)*x(3)**2*x(4)**3*x(5)**2 + 320*x(1)*x(3)**2*x(4)**2*x(5)**3 + 210*x(1)*x(3)**2*x(4)*x(5)**4 + 110*x(1)*x(3)**2*x(5)**5 + 5*x(1)*x(3)*x(4)**5*x(5) + 20*x(1)*x(3)*x(4)**4*x(5)**2 + 20*x(1)*x(3)*x(4)**3*x(5)**3 + 20*x(1)*x(3)*x(4)**2*x(5)**4 + 5*x(1)*x(3)*x(4)*x(5)**5 - 10*x(1)*x(4)**5*x(5)**2 - 10*x(1)*x(4)**4*x(5)**3 - 10*x(1)*x(4)**3*x(5)**4 - 10*x(1)*x(4)**2*x(5)**5 - 1662*x(2)**7*x(3) + 1662*x(2)**7*x(4) + 1662*x(2)**7*x(5) + 7910*x(2)**6*x(3)**2 - 4386*x(2)**6*x(3)*x(4) - 4386*x(2)**6*x(3)*x(5) - 3524*x(2)**6*x(4)**2 - 5286*x(2)**6*x(4)*x(5) - 3524*x(2)**6*x(5)**2 - 15430*x(2)**5*x(3)**3 - 780*x(2)**5*x(3)**2*x(4) - 780*x(2)**5*x(3)**2*x(5) + 14048*x(2)**5*x(3)*x(4)**2 + 18424*x(2)**5*x(3)*x(4)*x(5) + 14048*x(2)**5*x(3)*x(5)**2 + 2162*x(2)**5*x(4)**3 + 5991*x(2)**5*x(4)**2*x(5) + 5991*x(2)**5*x(4)*x(5)**2 + 2162*x(2)**5*x(5)**3 + 15430*x(2)**4*x(3)**4 + 15430*x(2)**4*x(3)**3*x(4) + 15430*x(2)**4*x(3)**3*x(5) - 21880*x(2)**4*x(3)**2*x(4)**2 - 22377*x(2)**4*x(3)**2*x(4)*x(5) - 21880*x(2)**4*x(3)**2*x(5)**2 - 8580*x(2)**4*x(3)*x(4)**3 - 21231*x(2)**4*x(3)*x(4)**2*x(5) - 21231*x(2)**4*x(3)*x(4)*x(5)**2 - 8580*x(2)**4*x(3)*x(5)**3 - 400*x(2)**4*x(4)**4 - 2877*x(2)**4*x(4)**3*x(5) - 3182*x(2)**4*x(4)**2*x(5)**2 - 2877*x(2)**4*x(4)*x(5)**3 - 400*x(2)**4*x(5)**4 - 7910*x(2)**3*x(3)**5 - 22560*x(2)**3*x(3)**4*x(4) - 22560*x(2)**3*x(3)**4*x(5) + 14750*x(2)**3*x(3)**3*x(4)**2 + 5173*x(2)**3*x(3)**3*x(4)*x(5) + 14750*x(2)**3*x(3)**3*x(5)**2 + 14750*x(2)**3*x(3)**2*x(4)**3 + 31874*x(2)**3*x(3)**2*x(4)**2*x(5) + 31874*x(2)**3*x(3)**2*x(4)*x(5)**2 + 14750*x(2)**3*x(3)**2*x(5)**3 + 870*x(2)**3*x(3)*x(4)**4 + 8043*x(2)**3*x(3)*x(4)**3*x(5) + 8308*x(2)**3*x(3)*x(4)**2*x(5)**2 + 8043*x(2)**3*x(3)*x(4)*x(5)**3 + 870*x(2)**3*x(3)*x(5)**4 + 100*x(2)**3*x(4)**5 + 615*x(2)**3*x(4)**4*x(5) + 930*x(2)**3*x(4)**3*x(5)**2 + 930*x(2)**3*x(4)**2*x(5)**3 + 615*x(2)**3*x(4)*x(5)**4 + 100*x(2)**3*x(5)**5 + 1662*x(2)**2*x(3)**6 + 13958*x(2)**2*x(3)**5*x(4) + 13958*x(2)**2*x(3)**5*x(5) - 870*x(2)**2*x(3)**4*x(4)**2 + 12608*x(2)**2*x(3)**4*x(4)*x(5) - 870*x(2)**2*x(3)**4*x(5)**2 - 14170*x(2)**2*x(3)**3*x(4)**3 - 25676*x(2)**2*x(3)**3*x(4)**2*x(5) - 25676*x(2)**2*x(3)**3*x(4)*x(5)**2 - 14170*x(2)**2*x(3)**3*x(5)**3 - 290*x(2)**2*x(3)**2*x(4)**4 - 9027*x(2)**2*x(3)**2*x(4)**3*x(5) - 8657*x(2)**2*x(3)**2*x(4)**2*x(5)**2 - 9027*x(2)**2*x(3)**2*x(4)*x(5)**3 - 290*x(2)**2*x(3)**2*x(5)**4 - 290*x(2)**2*x(3)*x(4)**5 - 1035*x(2)**2*x(3)*x(4)**4*x(5) - 1490*x(2)**2*x(3)*x(4)**3*x(5)**2 - 1490*x(2)**2*x(3)*x(4)**2*x(5)**3 - 1035*x(2)**2*x(3)*x(4)*x(5)**4 - 290*x(2)**2*x(3)*x(5)**5 - 105*x(2)**2*x(4)**5*x(5) - 220*x(2)**2*x(4)**4*x(5)**2 - 220*x(2)**2*x(4)**3*x(5)**3 - 220*x(2)**2*x(4)**2*x(5)**4 - 105*x(2)**2*x(4)*x(5)**5 - 3324*x(2)*x(3)**6*x(4) - 3324*x(2)*x(3)**6*x(5) - 4186*x(2)*x(3)**5*x(4)**2 - 11896*x(2)*x(3)**5*x(4)*x(5) - 4186*x(2)*x(3)**5*x(5)**2 + 7700*x(2)*x(3)**4*x(4)**3 + 11214*x(2)*x(3)**4*x(4)**2*x(5) + 11214*x(2)*x(3)**4*x(4)*x(5)**2 + 7700*x(2)*x(3)**4*x(5)**3 - 480*x(2)*x(3)**3*x(4)**4 + 4973*x(2)*x(3)**3*x(4)**3*x(5) + 4473*x(2)*x(3)**3*x(4)**2*x(5)**2 + 4973*x(2)*x(3)**3*x(4)*x(5)**3 - 480*x(2)*x(3)**3*x(5)**4 + 290*x(2)*x(3)**2*x(4)**5 + 410*x(2)*x(3)**2*x(4)**4*x(5) + 535*x(2)*x(3)**2*x(4)**3*x(5)**2 + 535*x(2)*x(3)**2*x(4)**2*x(5)**3 + 410*x(2)*x(3)**2*x(4)*x(5)**4 + 290*x(2)*x(3)**2*x(5)**5 + 185*x(2)*x(3)*x(4)**5*x(5) + 360*x(2)*x(3)*x(4)**4*x(5)**2 + 360*x(2)*x(3)*x(4)**3*x(5)**3 + 360*x(2)*x(3)*x(4)**2*x(5)**4 + 185*x(2)*x(3)*x(4)*x(5)**5 + 5*x(2)*x(4)**5*x(5)**2 + 5*x(2)*x(4)**4*x(5)**3 + 5*x(2)*x(4)**3*x(5)**4 + 5*x(2)*x(4)**2*x(5)**5 + 1662*x(3)**6*x(4)**2 + 3324*x(3)**6*x(4)*x(5) + 1662*x(3)**6*x(5)**2 - 1862*x(3)**5*x(4)**3 - 2062*x(3)**5*x(4)**2*x(5) - 2062*x(3)**5*x(4)*x(5)**2 - 1862*x(3)**5*x(5)**3 + 300*x(3)**4*x(4)**4 - 1262*x(3)**4*x(4)**3*x(5) - 1162*x(3)**4*x(4)**2*x(5)**2 - 1262*x(3)**4*x(4)*x(5)**3 + 300*x(3)**4*x(5)**4 - 100*x(3)**3*x(4)**5 + 100*x(3)**3*x(4)**4*x(5) + 205*x(3)**3*x(4)**3*x(5)**2 + 205*x(3)**3*x(4)**2*x(5)**3 + 100*x(3)**3*x(4)*x(5)**4 - 100*x(3)**3*x(5)**5 - 100*x(3)**2*x(4)**5*x(5) - 210*x(3)**2*x(4)**4*x(5)**2 - 210*x(3)**2*x(4)**3*x(5)**3 - 210*x(3)**2*x(4)**2*x(5)**4 - 100*x(3)**2*x(4)*x(5)**5 + 5*x(3)*x(4)**5*x(5)**2 + 5*x(3)*x(4)**4*x(5)**3 + 5*x(3)*x(4)**3*x(5)**4 + 5*x(3)*x(4)**2*x(5)**5))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2*(x(3) - x(5)))
                        
                            Cb(2,7) = (4*(x(2) - x(3))**2*(105*x(1)**4*x(2)**2 - 185*x(1)**4*x(2)*x(3) - 20*x(1)**4*x(2)*x(4) - 5*x(1)**4*x(2)*x(5) + 100*x(1)**4*x(3)**2 - 10*x(1)**4*x(3)*x(4) - 5*x(1)**4*x(3)*x(5) + 10*x(1)**4*x(4)**2 + 10*x(1)**4*x(4)*x(5) - 205*x(1)**3*x(2)**3 + 85*x(1)**3*x(2)**2*x(3) + 110*x(1)**3*x(2)**2*x(4) + 370*x(1)**3*x(2)*x(3)**2 - 170*x(1)**3*x(2)*x(3)*x(4) + 10*x(1)**3*x(2)*x(4)**2 + 10*x(1)**3*x(2)*x(4)*x(5) + 5*x(1)**3*x(2)*x(5)**2 - 300*x(1)**3*x(3)**3 + 120*x(1)**3*x(3)**2*x(4) + 10*x(1)**3*x(3)**2*x(5) - 10*x(1)**3*x(3)*x(4)**2 - 10*x(1)**3*x(3)*x(4)*x(5) + 5*x(1)**3*x(3)*x(5)**2 - 10*x(1)**3*x(4)**3 - 10*x(1)**3*x(4)**2*x(5) - 10*x(1)**3*x(4)*x(5)**2 + 1562*x(1)**2*x(2)**4 - 5948*x(1)**2*x(2)**3*x(3) + 110*x(1)**2*x(2)**3*x(4) + 205*x(1)**2*x(2)**3*x(5) + 9087*x(1)**2*x(2)**2*x(3)**2 - 300*x(1)**2*x(2)**2*x(3)*x(4) - 285*x(1)**2*x(2)**2*x(3)*x(5) - 120*x(1)**2*x(2)**2*x(4)**2 - 120*x(1)**2*x(2)**2*x(4)*x(5) - 105*x(1)**2*x(2)**2*x(5)**2 - 6533*x(1)**2*x(2)*x(3)**3 + 300*x(1)**2*x(2)*x(3)**2*x(4) + 15*x(1)**2*x(2)*x(3)**2*x(5) + 170*x(1)**2*x(2)*x(3)*x(4)**2 + 170*x(1)**2*x(2)*x(3)*x(4)*x(5) + 185*x(1)**2*x(2)*x(3)*x(5)**2 + 10*x(1)**2*x(2)*x(4)**3 + 10*x(1)**2*x(2)*x(4)**2*x(5) + 10*x(1)**2*x(2)*x(4)*x(5)**2 + 1862*x(1)**2*x(3)**4 - 110*x(1)**2*x(3)**3*x(4) + 95*x(1)**2*x(3)**3*x(5) - 110*x(1)**2*x(3)**2*x(4)**2 - 110*x(1)**2*x(3)**2*x(4)*x(5) - 110*x(1)**2*x(3)**2*x(5)**2 + 20*x(1)**2*x(3)*x(4)**3 + 20*x(1)**2*x(3)*x(4)**2*x(5) + 20*x(1)**2*x(3)*x(4)*x(5)**2 + 1562*x(1)*x(2)**5 - 7715*x(1)*x(2)**4*x(3) - 1657*x(1)*x(2)**4*x(4) - 1562*x(1)*x(2)**4*x(5) + 15430*x(1)*x(2)**3*x(3)**2 + 6043*x(1)*x(2)**3*x(3)*x(4) + 5853*x(1)*x(2)**3*x(3)*x(5) + 190*x(1)*x(2)**3*x(4)**2 - 15*x(1)*x(2)**3*x(4)*x(5) - 15625*x(1)*x(2)**2*x(3)**3 - 8792*x(1)*x(2)**2*x(3)**2*x(4) - 8797*x(1)*x(2)**2*x(3)**2*x(5) - 75*x(1)*x(2)**2*x(3)*x(4)**2 + 205*x(1)*x(2)**2*x(3)*x(4)*x(5) + 200*x(1)*x(2)**2*x(3)*x(5)**2 - 95*x(1)*x(2)**2*x(4)**3 + 30*x(1)*x(2)**2*x(4)**2*x(5) + 10*x(1)*x(2)**2*x(4)*x(5)**2 + 8010*x(1)*x(2)*x(3)**4 + 6038*x(1)*x(2)*x(3)**3*x(4) + 6238*x(1)*x(2)*x(3)**3*x(5) - 375*x(1)*x(2)*x(3)**2*x(4)**2 - 380*x(1)*x(2)*x(3)**2*x(4)*x(5) - 385*x(1)*x(2)*x(3)**2*x(5)**2 + 185*x(1)*x(2)*x(3)*x(4)**3 + 5*x(1)*x(2)*x(3)*x(4)**2*x(5) - 15*x(1)*x(2)*x(4)**3*x(5) - 20*x(1)*x(2)*x(4)**2*x(5)**2 - 1662*x(1)*x(3)**5 - 1662*x(1)*x(3)**4*x(4) - 1762*x(1)*x(3)**4*x(5) + 310*x(1)*x(3)**3*x(4)**2 + 210*x(1)*x(3)**3*x(4)*x(5) + 205*x(1)*x(3)**3*x(5)**2 - 110*x(1)*x(3)**2*x(4)**3 - 5*x(1)*x(3)**2*x(4)**2*x(5) - 10*x(1)*x(3)**2*x(4)*x(5)**2 - 5*x(1)*x(3)*x(4)**3*x(5) - 10*x(1)*x(3)*x(4)**2*x(5)**2 + 10*x(1)*x(4)**3*x(5)**2 + 1662*x(2)**6 - 7910*x(2)**5*x(3) - 1862*x(2)**5*x(4) - 1762*x(2)**5*x(5) + 15430*x(2)**4*x(3)**2 + 8300*x(2)**4*x(3)*x(4) + 8105*x(2)**4*x(3)*x(5) + 300*x(2)**4*x(4)**2 + 2067*x(2)**4*x(4)*x(5) + 100*x(2)**4*x(5)**2 - 15430*x(2)**3*x(3)**3 - 15430*x(2)**3*x(3)**2*x(4) - 15430*x(2)**3*x(3)**2*x(5) - 680*x(2)**3*x(3)*x(4)**2 - 7023*x(2)**3*x(3)*x(4)*x(5) - 195*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 - 410*x(2)**3*x(4)**2*x(5) - 205*x(2)**3*x(4)*x(5)**2 + 7910*x(2)**2*x(3)**4 + 15040*x(2)**2*x(3)**3*x(4) + 15235*x(2)**2*x(3)**3*x(5) + 290*x(2)**2*x(3)**2*x(4)**2 + 9382*x(2)**2*x(3)**2*x(4)*x(5) + 290*x(2)**2*x(3)*x(4)**3 + 665*x(2)**2*x(3)*x(4)**2*x(5) + 385*x(2)**2*x(3)*x(4)*x(5)**2 + 105*x(2)**2*x(4)**3*x(5) + 110*x(2)**2*x(4)**2*x(5)**2 - 1662*x(2)*x(3)**5 - 7710*x(2)*x(3)**4*x(4) - 7810*x(2)*x(3)**4*x(5) + 290*x(2)*x(3)**3*x(4)**2 - 5858*x(2)*x(3)**3*x(4)*x(5) + 195*x(2)*x(3)**3*x(5)**2 - 290*x(2)*x(3)**2*x(4)**3 - 205*x(2)*x(3)**2*x(4)**2*x(5) - 200*x(2)*x(3)**2*x(4)*x(5)**2 - 185*x(2)*x(3)*x(4)**3*x(5) - 185*x(2)*x(3)*x(4)**2*x(5)**2 - 5*x(2)*x(4)**3*x(5)**2 + 1662*x(3)**5*x(4) + 1662*x(3)**5*x(5) - 200*x(3)**4*x(4)**2 + 1462*x(3)**4*x(4)*x(5) - 100*x(3)**4*x(5)**2 + 100*x(3)**3*x(4)**3 - 100*x(3)**3*x(4)**2*x(5) + 100*x(3)**2*x(4)**3*x(5) + 105*x(3)**2*x(4)**2*x(5)**2 - 5*x(3)*x(4)**3*x(5)**2))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(2,8) = (24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(2) - x(3))**6 + (x(2)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(2) - x(3)))/5 - (x(3)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(2) - x(3)))/5 + x(2)*(x(2) - x(3))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(3)*(x(2) - x(3))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(2)**3*(x(2) - x(3))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 - (x(3)**3*(x(2) - x(3))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 + (x(2)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(2) - x(3))**3)/3 - (x(3)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(2) - x(3))**3)/3 + x(2)*(x(2) - x(3))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(3)*(x(2) - x(3))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(2)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(2) - x(3))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(3)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(2) - x(3))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - x(2)**2*(x(2) - x(3))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(3)**2*(x(2) - x(3))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - (x(2)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(2) - x(3))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2 + (x(3)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(2) - x(3))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2

                            Cb(2,9) = (4*(x(2) - x(3))**2*(- 100*x(1)**3*x(2)**3 + 85*x(1)**3*x(2)**2*x(3) + 110*x(1)**3*x(2)**2*x(4) + 105*x(1)**3*x(2)**2*x(5) + 85*x(1)**3*x(2)*x(3)**2 - 160*x(1)**3*x(2)*x(3)*x(4) - 180*x(1)**3*x(2)*x(3)*x(5) - 20*x(1)**3*x(2)*x(4)**2 - 20*x(1)**3*x(2)*x(4)*x(5) - 5*x(1)**3*x(2)*x(5)**2 - 100*x(1)**3*x(3)**3 + 110*x(1)**3*x(3)**2*x(4) + 105*x(1)**3*x(3)**2*x(5) - 20*x(1)**3*x(3)*x(4)**2 - 20*x(1)**3*x(3)*x(4)*x(5) - 5*x(1)**3*x(3)*x(5)**2 + 10*x(1)**3*x(4)**3 + 10*x(1)**3*x(4)**2*x(5) + 10*x(1)**3*x(4)*x(5)**2 + 200*x(1)**2*x(2)**4 + 15*x(1)**2*x(2)**3*x(3) - 310*x(1)**2*x(2)**3*x(4) - 205*x(1)**2*x(2)**3*x(5) - 380*x(1)**2*x(2)**2*x(3)**2 + 270*x(1)**2*x(2)**2*x(3)*x(4) + 190*x(1)**2*x(2)**2*x(3)*x(5) + 110*x(1)**2*x(2)**2*x(4)**2 + 110*x(1)**2*x(2)**2*x(4)*x(5) + 15*x(1)**2*x(2)*x(3)**3 + 270*x(1)**2*x(2)*x(3)**2*x(4) + 190*x(1)**2*x(2)*x(3)**2*x(5) - 200*x(1)**2*x(2)*x(3)*x(4)**2 - 200*x(1)**2*x(2)*x(3)*x(4)*x(5) - 10*x(1)**2*x(2)*x(3)*x(5)**2 + 10*x(1)**2*x(2)*x(4)**3 + 10*x(1)**2*x(2)*x(4)**2*x(5) + 10*x(1)**2*x(2)*x(4)*x(5)**2 + 5*x(1)**2*x(2)*x(5)**3 + 200*x(1)**2*x(3)**4 - 310*x(1)**2*x(3)**3*x(4) - 205*x(1)**2*x(3)**3*x(5) + 110*x(1)**2*x(3)**2*x(4)**2 + 110*x(1)**2*x(3)**2*x(4)*x(5) + 10*x(1)**2*x(3)*x(4)**3 + 10*x(1)**2*x(3)*x(4)**2*x(5) + 10*x(1)**2*x(3)*x(4)*x(5)**2 + 5*x(1)**2*x(3)*x(5)**3 - 10*x(1)**2*x(4)**4 - 10*x(1)**2*x(4)**3*x(5) - 10*x(1)**2*x(4)**2*x(5)**2 - 10*x(1)**2*x(4)*x(5)**3 - 1662*x(1)*x(2)**5 + 4686*x(1)*x(2)**4*x(3) + 1662*x(1)*x(2)**4*x(4) + 1562*x(1)*x(2)**4*x(5) - 3029*x(1)*x(2)**3*x(3)**2 - 6358*x(1)*x(2)**3*x(3)*x(4) - 6358*x(1)*x(2)**3*x(3)*x(5) + 110*x(1)*x(2)**3*x(4)**2 + 110*x(1)*x(2)**3*x(4)*x(5) + 205*x(1)*x(2)**3*x(5)**2 - 3029*x(1)*x(2)**2*x(3)**3 + 9372*x(1)*x(2)**2*x(3)**2*x(4) + 9562*x(1)*x(2)**2*x(3)**2*x(5) - 70*x(1)*x(2)**2*x(3)*x(4)**2 - 70*x(1)*x(2)**2*x(3)*x(4)*x(5) - 180*x(1)*x(2)**2*x(3)*x(5)**2 - 120*x(1)*x(2)**2*x(4)**3 - 120*x(1)*x(2)**2*x(4)**2*x(5) - 120*x(1)*x(2)**2*x(4)*x(5)**2 - 105*x(1)*x(2)**2*x(5)**3 + 4686*x(1)*x(2)*x(3)**4 - 6358*x(1)*x(2)*x(3)**3*x(4) - 6358*x(1)*x(2)*x(3)**3*x(5) - 70*x(1)*x(2)*x(3)**2*x(4)**2 - 70*x(1)*x(2)*x(3)**2*x(4)*x(5) - 180*x(1)*x(2)*x(3)**2*x(5)**2 + 170*x(1)*x(2)*x(3)*x(4)**3 + 170*x(1)*x(2)*x(3)*x(4)**2*x(5) + 170*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 190*x(1)*x(2)*x(3)*x(5)**3 + 10*x(1)*x(2)*x(4)**4 + 10*x(1)*x(2)*x(4)**3*x(5) + 10*x(1)*x(2)*x(4)**2*x(5)**2 + 10*x(1)*x(2)*x(4)*x(5)**3 - 1662*x(1)*x(3)**5 + 1662*x(1)*x(3)**4*x(4) + 1562*x(1)*x(3)**4*x(5) + 110*x(1)*x(3)**3*x(4)**2 + 110*x(1)*x(3)**3*x(4)*x(5) + 205*x(1)*x(3)**3*x(5)**2 - 120*x(1)*x(3)**2*x(4)**3 - 120*x(1)*x(3)**2*x(4)**2*x(5) - 120*x(1)*x(3)**2*x(4)*x(5)**2 - 105*x(1)*x(3)**2*x(5)**3 + 10*x(1)*x(3)*x(4)**4 + 10*x(1)*x(3)*x(4)**3*x(5) + 10*x(1)*x(3)*x(4)**2*x(5)**2 + 10*x(1)*x(3)*x(4)*x(5)**3 - 1662*x(2)**5*x(3) + 1662*x(2)**5*x(4) + 1662*x(2)**5*x(5) + 6248*x(2)**4*x(3)**2 - 4386*x(2)**4*x(3)*x(4) - 4486*x(2)**4*x(3)*x(5) - 1862*x(2)**4*x(4)**2 - 1862*x(2)**4*x(4)*x(5) - 1762*x(2)**4*x(5)**2 - 9182*x(2)**3*x(3)**3 + 2744*x(2)**3*x(3)**2*x(4) + 2839*x(2)**3*x(3)**2*x(5) + 6138*x(2)**3*x(3)*x(4)**2 + 6138*x(2)**3*x(3)*x(4)*x(5) + 6243*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 + 300*x(2)**3*x(4)**2*x(5) + 300*x(2)**3*x(4)*x(5)**2 + 100*x(2)**3*x(5)**3 + 6248*x(2)**2*x(3)**4 + 2744*x(2)**2*x(3)**3*x(4) + 2839*x(2)**2*x(3)**3*x(5) - 8612*x(2)**2*x(3)**2*x(4)**2 - 8612*x(2)**2*x(3)**2*x(4)*x(5) - 8992*x(2)**2*x(3)**2*x(5)**2 - 280*x(2)**2*x(3)*x(4)**3 - 280*x(2)**2*x(3)*x(4)**2*x(5) - 280*x(2)**2*x(3)*x(4)*x(5)**2 - 95*x(2)**2*x(3)*x(5)**3 - 100*x(2)**2*x(4)**4 - 100*x(2)**2*x(4)**3*x(5) - 100*x(2)**2*x(4)**2*x(5)**2 - 100*x(2)**2*x(4)*x(5)**3 - 1662*x(2)*x(3)**5 - 4386*x(2)*x(3)**4*x(4) - 4486*x(2)*x(3)**4*x(5) + 6138*x(2)*x(3)**3*x(4)**2 + 6138*x(2)*x(3)**3*x(4)*x(5) + 6243*x(2)*x(3)**3*x(5)**2 - 280*x(2)*x(3)**2*x(4)**3 - 280*x(2)*x(3)**2*x(4)**2*x(5) - 280*x(2)*x(3)**2*x(4)*x(5)**2 - 95*x(2)*x(3)**2*x(5)**3 + 190*x(2)*x(3)*x(4)**4 + 190*x(2)*x(3)*x(4)**3*x(5) + 190*x(2)*x(3)*x(4)**2*x(5)**2 + 190*x(2)*x(3)*x(4)*x(5)**3 + 1662*x(3)**5*x(4) + 1662*x(3)**5*x(5) - 1862*x(3)**4*x(4)**2 - 1862*x(3)**4*x(4)*x(5) - 1762*x(3)**4*x(5)**2 + 300*x(3)**3*x(4)**3 + 300*x(3)**3*x(4)**2*x(5) + 300*x(3)**3*x(4)*x(5)**2 + 100*x(3)**3*x(5)**3 - 100*x(3)**2*x(4)**4 - 100*x(3)**2*x(4)**3*x(5) - 100*x(3)**2*x(4)**2*x(5)**2 - 100*x(3)**2*x(4)*x(5)**3))/(5*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            Cb(2,10) = (4*(x(2) - x(3))**2*(50*x(1)**2*x(2)**2 - 95*x(1)**2*x(2)*x(3) - 5*x(1)**2*x(2)*x(4) + 50*x(1)**2*x(3)**2 - 5*x(1)**2*x(3)*x(4) + 5*x(1)**2*x(4)**2 - 100*x(1)*x(2)**3 + 95*x(1)*x(2)**2*x(3) + 105*x(1)*x(2)**2*x(4) + 95*x(1)*x(2)*x(3)**2 - 190*x(1)*x(2)*x(3)*x(4) - 5*x(1)*x(2)*x(4)**2 - 100*x(1)*x(3)**3 + 105*x(1)*x(3)**2*x(4) - 5*x(1)*x(3)*x(4)**2 + 831*x(2)**4 - 3124*x(2)**3*x(3) - 100*x(2)**3*x(4) + 4591*x(2)**2*x(3)**2 + 95*x(2)**2*x(3)*x(4) + 50*x(2)**2*x(4)**2 - 3124*x(2)*x(3)**3 + 95*x(2)*x(3)**2*x(4) - 95*x(2)*x(3)*x(4)**2 + 831*x(3)**4 - 100*x(3)**3*x(4) + 50*x(3)**2*x(4)**2))/(5*(x(1) - x(5))**2*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            x = x_cb(i:i+4)
                            x = x - x(1)

                            Cb(3,1) = -(4*(x(1) - x(2))**2*(- 996*x(1)**4 + 3309*x(1)**3*x(2) + 225*x(1)**3*x(3) + 225*x(1)**3*x(4) + 225*x(1)**3*x(5) - 4551*x(1)**2*x(2)**2 - 275*x(1)**2*x(2)*x(3) - 275*x(1)**2*x(2)*x(4) - 275*x(1)**2*x(2)*x(5) - 60*x(1)**2*x(3)**2 - 140*x(1)**2*x(3)*x(4) - 140*x(1)**2*x(3)*x(5) - 60*x(1)**2*x(4)**2 - 140*x(1)**2*x(4)*x(5) - 60*x(1)**2*x(5)**2 + 3024*x(1)*x(2)**3 + 10*x(1)*x(2)**2*x(3) + 10*x(1)*x(2)**2*x(4) + 10*x(1)*x(2)**2*x(5) + 90*x(1)*x(2)*x(3)**2 + 175*x(1)*x(2)*x(3)*x(4) + 175*x(1)*x(2)*x(3)*x(5) + 90*x(1)*x(2)*x(4)**2 + 175*x(1)*x(2)*x(4)*x(5) + 90*x(1)*x(2)*x(5)**2 + 15*x(1)*x(3)**2*x(4) + 15*x(1)*x(3)**2*x(5) + 15*x(1)*x(3)*x(4)**2 + 45*x(1)*x(3)*x(4)*x(5) + 15*x(1)*x(3)*x(5)**2 + 15*x(1)*x(4)**2*x(5) + 15*x(1)*x(4)*x(5)**2 - 831*x(2)**4 + 100*x(2)**3*x(3) + 100*x(2)**3*x(4) + 100*x(2)**3*x(5) - 50*x(2)**2*x(3)**2 - 105*x(2)**2*x(3)*x(4) - 105*x(2)**2*x(3)*x(5) - 50*x(2)**2*x(4)**2 - 105*x(2)**2*x(4)*x(5) - 50*x(2)**2*x(5)**2 + 5*x(2)*x(3)**2*x(4) + 5*x(2)*x(3)**2*x(5) + 5*x(2)*x(3)*x(4)**2 + 15*x(2)*x(3)*x(4)*x(5) + 5*x(2)*x(3)*x(5)**2 + 5*x(2)*x(4)**2*x(5) + 5*x(2)*x(4)*x(5)**2 - 5*x(3)**2*x(4)**2 - 10*x(3)**2*x(4)*x(5) - 5*x(3)**2*x(5)**2 - 10*x(3)*x(4)**2*x(5) - 10*x(3)*x(4)*x(5)**2 - 5*x(4)**2*x(5)**2))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2)

                            Cb(3,2) = (4*(x(1) - x(2))**2*(- 1767*x(1)**6 + 4371*x(1)**5*x(2) + 2077*x(1)**5*x(3) + 2077*x(1)**5*x(4) + 2077*x(1)**5*x(5) - 4461*x(1)**4*x(2)**2 - 4311*x(1)**4*x(2)*x(3) - 4311*x(1)**4*x(2)*x(4) - 4311*x(1)**4*x(2)*x(5) - 415*x(1)**4*x(3)**2 - 2622*x(1)**4*x(3)*x(4) - 2622*x(1)**4*x(3)*x(5) - 415*x(1)**4*x(4)**2 - 2622*x(1)**4*x(4)*x(5) - 415*x(1)**4*x(5)**2 + 3669*x(1)**3*x(2)**3 + 2279*x(1)**3*x(2)**2*x(3) + 2279*x(1)**3*x(2)**2*x(4) + 2279*x(1)**3*x(2)**2*x(5) + 5*x(1)**3*x(2)*x(3)**2 + 6338*x(1)**3*x(2)*x(3)*x(4) + 6338*x(1)**3*x(2)*x(3)*x(5) + 5*x(1)**3*x(2)*x(4)**2 + 6338*x(1)**3*x(2)*x(4)*x(5) + 5*x(1)**3*x(2)*x(5)**2 + 105*x(1)**3*x(3)**3 + 670*x(1)**3*x(3)**2*x(4) + 670*x(1)**3*x(3)**2*x(5) + 670*x(1)**3*x(3)*x(4)**2 + 1470*x(1)**3*x(3)*x(4)*x(5) + 670*x(1)**3*x(3)*x(5)**2 + 105*x(1)**3*x(4)**3 + 670*x(1)**3*x(4)**2*x(5) + 670*x(1)**3*x(4)*x(5)**2 + 105*x(1)**3*x(5)**3 - 4716*x(1)**2*x(2)**4 + 2619*x(1)**2*x(2)**3*x(3) + 2619*x(1)**2*x(2)**3*x(4) + 2619*x(1)**2*x(2)**3*x(5) + 600*x(1)**2*x(2)**2*x(3)**2 - 7947*x(1)**2*x(2)**2*x(3)*x(4) - 7947*x(1)**2*x(2)**2*x(3)*x(5) + 600*x(1)**2*x(2)**2*x(4)**2 - 7947*x(1)**2*x(2)**2*x(4)*x(5) + 600*x(1)**2*x(2)**2*x(5)**2 - 65*x(1)**2*x(2)*x(3)**3 - 510*x(1)**2*x(2)*x(3)**2*x(4) - 510*x(1)**2*x(2)*x(3)**2*x(5) - 510*x(1)**2*x(2)*x(3)*x(4)**2 - 1080*x(1)**2*x(2)*x(3)*x(4)*x(5) - 510*x(1)**2*x(2)*x(3)*x(5)**2 - 65*x(1)**2*x(2)*x(4)**3 - 510*x(1)**2*x(2)*x(4)**2*x(5) - 510*x(1)**2*x(2)*x(4)*x(5)**2 - 65*x(1)**2*x(2)*x(5)**3 - 125*x(1)**2*x(3)**3*x(4) - 125*x(1)**2*x(3)**3*x(5) - 285*x(1)**2*x(3)**2*x(4)**2 - 555*x(1)**2*x(3)**2*x(4)*x(5) - 285*x(1)**2*x(3)**2*x(5)**2 - 125*x(1)**2*x(3)*x(4)**3 - 555*x(1)**2*x(3)*x(4)**2*x(5) - 555*x(1)**2*x(3)*x(4)*x(5)**2 - 125*x(1)**2*x(3)*x(5)**3 - 125*x(1)**2*x(4)**3*x(5) - 285*x(1)**2*x(4)**2*x(5)**2 - 125*x(1)**2*x(4)*x(5)**3 + 4386*x(1)*x(2)**5 - 4166*x(1)*x(2)**4*x(3) - 4166*x(1)*x(2)**4*x(4) - 4166*x(1)*x(2)**4*x(5) - 140*x(1)*x(2)**3*x(3)**2 + 5853*x(1)*x(2)**3*x(3)*x(4) + 5853*x(1)*x(2)**3*x(3)*x(5) - 140*x(1)*x(2)**3*x(4)**2 + 5853*x(1)*x(2)**3*x(4)*x(5) - 140*x(1)*x(2)**3*x(5)**2 - 80*x(1)*x(2)**2*x(3)**3 - 270*x(1)*x(2)**2*x(3)**2*x(4) - 270*x(1)*x(2)**2*x(3)**2*x(5) - 270*x(1)*x(2)**2*x(3)*x(4)**2 - 585*x(1)*x(2)**2*x(3)*x(4)*x(5) - 270*x(1)*x(2)**2*x(3)*x(5)**2 - 80*x(1)*x(2)**2*x(4)**3 - 270*x(1)*x(2)**2*x(4)**2*x(5) - 270*x(1)*x(2)**2*x(4)*x(5)**2 - 80*x(1)*x(2)**2*x(5)**3 + 145*x(1)*x(2)*x(3)**3*x(4) + 145*x(1)*x(2)*x(3)**3*x(5) + 285*x(1)*x(2)*x(3)**2*x(4)**2 + 555*x(1)*x(2)*x(3)**2*x(4)*x(5) + 285*x(1)*x(2)*x(3)**2*x(5)**2 + 145*x(1)*x(2)*x(3)*x(4)**3 + 555*x(1)*x(2)*x(3)*x(4)**2*x(5) + 555*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 145*x(1)*x(2)*x(3)*x(5)**3 + 145*x(1)*x(2)*x(4)**3*x(5) + 285*x(1)*x(2)*x(4)**2*x(5)**2 + 145*x(1)*x(2)*x(4)*x(5)**3 + 30*x(1)*x(3)**3*x(4)**2 + 45*x(1)*x(3)**3*x(4)*x(5) + 30*x(1)*x(3)**3*x(5)**2 + 30*x(1)*x(3)**2*x(4)**3 + 105*x(1)*x(3)**2*x(4)**2*x(5) + 105*x(1)*x(3)**2*x(4)*x(5)**2 + 30*x(1)*x(3)**2*x(5)**3 + 45*x(1)*x(3)*x(4)**3*x(5) + 105*x(1)*x(3)*x(4)**2*x(5)**2 + 45*x(1)*x(3)*x(4)*x(5)**3 + 30*x(1)*x(4)**3*x(5)**2 + 30*x(1)*x(4)**2*x(5)**3 - 1662*x(2)**6 + 1862*x(2)**5*x(3) + 1862*x(2)**5*x(4) + 1862*x(2)**5*x(5) - 300*x(2)**4*x(3)**2 - 2272*x(2)**4*x(3)*x(4) - 2272*x(2)**4*x(3)*x(5) - 300*x(2)**4*x(4)**2 - 2272*x(2)**4*x(4)*x(5) - 300*x(2)**4*x(5)**2 + 100*x(2)**3*x(3)**3 + 520*x(2)**3*x(3)**2*x(4) + 520*x(2)**3*x(3)**2*x(5) + 520*x(2)**3*x(3)*x(4)**2 + 1155*x(2)**3*x(3)*x(4)*x(5) + 520*x(2)**3*x(3)*x(5)**2 + 100*x(2)**3*x(4)**3 + 520*x(2)**3*x(4)**2*x(5) + 520*x(2)**3*x(4)*x(5)**2 + 100*x(2)**3*x(5)**3 - 110*x(2)**2*x(3)**3*x(4) - 110*x(2)**2*x(3)**3*x(5) - 240*x(2)**2*x(3)**2*x(4)**2 - 480*x(2)**2*x(3)**2*x(4)*x(5) - 240*x(2)**2*x(3)**2*x(5)**2 - 110*x(2)**2*x(3)*x(4)**3 - 480*x(2)**2*x(3)*x(4)**2*x(5) - 480*x(2)**2*x(3)*x(4)*x(5)**2 - 110*x(2)**2*x(3)*x(5)**3 - 110*x(2)**2*x(4)**3*x(5) - 240*x(2)**2*x(4)**2*x(5)**2 - 110*x(2)**2*x(4)*x(5)**3 + 20*x(2)*x(3)**3*x(4)**2 + 35*x(2)*x(3)**3*x(4)*x(5) + 20*x(2)*x(3)**3*x(5)**2 + 20*x(2)*x(3)**2*x(4)**3 + 75*x(2)*x(3)**2*x(4)**2*x(5) + 75*x(2)*x(3)**2*x(4)*x(5)**2 + 20*x(2)*x(3)**2*x(5)**3 + 35*x(2)*x(3)*x(4)**3*x(5) + 75*x(2)*x(3)*x(4)**2*x(5)**2 + 35*x(2)*x(3)*x(4)*x(5)**3 + 20*x(2)*x(4)**3*x(5)**2 + 20*x(2)*x(4)**2*x(5)**3 - 10*x(3)**3*x(4)**3 - 20*x(3)**3*x(4)**2*x(5) - 20*x(3)**3*x(4)*x(5)**2 - 10*x(3)**3*x(5)**3 - 20*x(3)**2*x(4)**3*x(5) - 30*x(3)**2*x(4)**2*x(5)**2 - 20*x(3)**2*x(4)*x(5)**3 - 20*x(3)*x(4)**3*x(5)**2 - 20*x(3)*x(4)**2*x(5)**3 - 10*x(4)**3*x(5)**3))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5)))

                            Cb(3,3) = -(4*(x(1) - x(2))**2*(- 1767*x(1)**5*x(2) - 1767*x(1)**5*x(3) + 1767*x(1)**5*x(4) + 1767*x(1)**5*x(5) + 6343*x(1)**4*x(2)**2 + 4681*x(1)**4*x(2)*x(3) - 4266*x(1)**4*x(2)*x(4) - 4266*x(1)**4*x(2)*x(5) + 310*x(1)**4*x(3)**2 + 1767*x(1)**4*x(3)*x(4) + 1767*x(1)**4*x(3)*x(5) - 2077*x(1)**4*x(4)**2 - 2182*x(1)**4*x(4)*x(5) - 2077*x(1)**4*x(5)**2 - 9092*x(1)**3*x(2)**3 - 2844*x(1)**3*x(2)**2*x(3) + 2374*x(1)**3*x(2)**2*x(4) + 2374*x(1)**3*x(2)**2*x(5) - 45*x(1)**3*x(2)*x(3)**2 - 6473*x(1)**3*x(2)*x(3)*x(4) - 6473*x(1)**3*x(2)*x(3)*x(5) + 6303*x(1)**3*x(2)*x(4)**2 + 6183*x(1)**3*x(2)*x(4)*x(5) + 6303*x(1)**3*x(2)*x(5)**2 - 105*x(1)**3*x(3)**3 - 440*x(1)**3*x(3)**2*x(4) - 440*x(1)**3*x(3)**2*x(5) + 130*x(1)**3*x(3)*x(4)**2 + 25*x(1)**3*x(3)*x(4)*x(5) + 130*x(1)**3*x(3)*x(5)**2 + 415*x(1)**3*x(4)**3 + 630*x(1)**3*x(4)**2*x(5) + 630*x(1)**3*x(4)*x(5)**2 + 415*x(1)**3*x(5)**3 + 6148*x(1)**2*x(2)**4 - 3034*x(1)**2*x(2)**3*x(3) + 2859*x(1)**2*x(2)**3*x(4) + 2859*x(1)**2*x(2)**3*x(5) - 455*x(1)**2*x(2)**2*x(3)**2 + 9272*x(1)**2*x(2)**2*x(3)*x(4) + 9272*x(1)**2*x(2)**2*x(3)*x(5) - 8447*x(1)**2*x(2)**2*x(4)**2 - 8077*x(1)**2*x(2)**2*x(4)*x(5) - 8447*x(1)**2*x(2)**2*x(5)**2 + 65*x(1)**2*x(2)*x(3)**3 + 425*x(1)**2*x(2)*x(3)**2*x(4) + 425*x(1)**2*x(2)*x(3)**2*x(5) - 30*x(1)**2*x(2)*x(3)*x(4)**2 + 85*x(1)**2*x(2)*x(3)*x(4)*x(5) - 30*x(1)**2*x(2)*x(3)*x(5)**2 - 455*x(1)**2*x(2)*x(4)**3 - 620*x(1)**2*x(2)*x(4)**2*x(5) - 620*x(1)**2*x(2)*x(4)*x(5)**2 - 455*x(1)**2*x(2)*x(5)**3 + 125*x(1)**2*x(3)**3*x(4) + 125*x(1)**2*x(3)**3*x(5) + 130*x(1)**2*x(3)**2*x(4)**2 + 260*x(1)**2*x(3)**2*x(4)*x(5) + 130*x(1)**2*x(3)**2*x(5)**2 - 150*x(1)**2*x(3)*x(4)**3 - 170*x(1)**2*x(3)*x(4)**2*x(5) - 170*x(1)**2*x(3)*x(4)*x(5)**2 - 150*x(1)**2*x(3)*x(5)**3 - 105*x(1)**2*x(4)**4 - 220*x(1)**2*x(4)**3*x(5) - 220*x(1)**2*x(4)**2*x(5)**2 - 220*x(1)**2*x(4)*x(5)**3 - 105*x(1)**2*x(5)**4 - 1662*x(1)*x(2)**5 + 4586*x(1)*x(2)**4*x(3) - 4286*x(1)*x(2)**4*x(4) - 4286*x(1)*x(2)**4*x(5) + 120*x(1)*x(2)**3*x(3)**2 - 6258*x(1)*x(2)**3*x(3)*x(4) - 6258*x(1)*x(2)**3*x(3)*x(5) + 5933*x(1)*x(2)**3*x(4)**2 + 5818*x(1)*x(2)**3*x(4)*x(5) + 5933*x(1)*x(2)**3*x(5)**2 + 80*x(1)*x(2)**2*x(3)**3 + 155*x(1)*x(2)**2*x(3)**2*x(4) + 155*x(1)*x(2)**2*x(3)**2*x(5) - 60*x(1)*x(2)**2*x(3)*x(4)**2 + 40*x(1)*x(2)**2*x(3)*x(4)*x(5) - 60*x(1)*x(2)**2*x(3)*x(5)**2 - 170*x(1)*x(2)**2*x(4)**3 - 335*x(1)*x(2)**2*x(4)**2*x(5) - 335*x(1)*x(2)**2*x(4)*x(5)**2 - 170*x(1)*x(2)**2*x(5)**3 - 145*x(1)*x(2)*x(3)**3*x(4) - 145*x(1)*x(2)*x(3)**3*x(5) - 190*x(1)*x(2)*x(3)**2*x(4)**2 - 345*x(1)*x(2)*x(3)**2*x(4)*x(5) - 190*x(1)*x(2)*x(3)**2*x(5)**2 + 150*x(1)*x(2)*x(3)*x(4)**3 + 110*x(1)*x(2)*x(3)*x(4)**2*x(5) + 110*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 150*x(1)*x(2)*x(3)*x(5)**3 + 185*x(1)*x(2)*x(4)**4 + 360*x(1)*x(2)*x(4)**3*x(5) + 360*x(1)*x(2)*x(4)**2*x(5)**2 + 360*x(1)*x(2)*x(4)*x(5)**3 + 185*x(1)*x(2)*x(5)**4 - 30*x(1)*x(3)**3*x(4)**2 - 45*x(1)*x(3)**3*x(4)*x(5) - 30*x(1)*x(3)**3*x(5)**2 + 10*x(1)*x(3)**2*x(4)**3 - 10*x(1)*x(3)**2*x(4)**2*x(5) - 10*x(1)*x(3)**2*x(4)*x(5)**2 + 10*x(1)*x(3)**2*x(5)**3 + 20*x(1)*x(3)*x(4)**4 + 50*x(1)*x(3)*x(4)**3*x(5) + 50*x(1)*x(3)*x(4)**2*x(5)**2 + 50*x(1)*x(3)*x(4)*x(5)**3 + 20*x(1)*x(3)*x(5)**4 + 5*x(1)*x(4)**4*x(5) + 5*x(1)*x(4)**3*x(5)**2 + 5*x(1)*x(4)**2*x(5)**3 + 5*x(1)*x(4)*x(5)**4 - 1662*x(2)**5*x(3) + 1662*x(2)**5*x(4) + 1662*x(2)**5*x(5) + 200*x(2)**4*x(3)**2 + 1662*x(2)**4*x(3)*x(4) + 1662*x(2)**4*x(3)*x(5) - 1862*x(2)**4*x(4)**2 - 1962*x(2)**4*x(4)*x(5) - 1862*x(2)**4*x(5)**2 - 100*x(2)**3*x(3)**3 - 310*x(2)**3*x(3)**2*x(4) - 310*x(2)**3*x(3)**2*x(5) + 110*x(2)**3*x(3)*x(4)**2 + 10*x(2)**3*x(3)*x(4)*x(5) + 110*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 + 505*x(2)**3*x(4)**2*x(5) + 505*x(2)**3*x(4)*x(5)**2 + 300*x(2)**3*x(5)**3 + 110*x(2)**2*x(3)**3*x(4) + 110*x(2)**2*x(3)**3*x(5) + 110*x(2)**2*x(3)**2*x(4)**2 + 225*x(2)**2*x(3)**2*x(4)*x(5) + 110*x(2)**2*x(3)**2*x(5)**2 - 120*x(2)**2*x(3)*x(4)**3 - 130*x(2)**2*x(3)*x(4)**2*x(5) - 130*x(2)**2*x(3)*x(4)*x(5)**2 - 120*x(2)**2*x(3)*x(5)**3 - 100*x(2)**2*x(4)**4 - 210*x(2)**2*x(4)**3*x(5) - 210*x(2)**2*x(4)**2*x(5)**2 - 210*x(2)**2*x(4)*x(5)**3 - 100*x(2)**2*x(5)**4 - 20*x(2)*x(3)**3*x(4)**2 - 35*x(2)*x(3)**3*x(4)*x(5) - 20*x(2)*x(3)**3*x(5)**2 + 10*x(2)*x(3)**2*x(4)**3 + 10*x(2)*x(3)**2*x(5)**3 + 10*x(2)*x(3)*x(4)**4 + 30*x(2)*x(3)*x(4)**3*x(5) + 30*x(2)*x(3)*x(4)**2*x(5)**2 + 30*x(2)*x(3)*x(4)*x(5)**3 + 10*x(2)*x(3)*x(5)**4 + 5*x(2)*x(4)**4*x(5) + 5*x(2)*x(4)**3*x(5)**2 + 5*x(2)*x(4)**2*x(5)**3 + 5*x(2)*x(4)*x(5)**4 + 10*x(3)**3*x(4)**3 + 20*x(3)**3*x(4)**2*x(5) + 20*x(3)**3*x(4)*x(5)**2 + 10*x(3)**3*x(5)**3 - 10*x(3)**2*x(4)**4 - 10*x(3)**2*x(4)**3*x(5) - 10*x(3)**2*x(4)**2*x(5)**2 - 10*x(3)**2*x(4)*x(5)**3 - 10*x(3)**2*x(5)**4 - 10*x(3)*x(4)**4*x(5) - 10*x(3)*x(4)**3*x(5)**2 - 10*x(3)*x(4)**2*x(5)**3 - 10*x(3)*x(4)*x(5)**4))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(3,4) = (4*(x(1) - x(2))**2*(- 1767*x(1)**4 + 6343*x(1)**3*x(2) + 310*x(1)**3*x(3) + 310*x(1)**3*x(4) + 105*x(5)*x(1)**3 - 9092*x(1)**2*x(2)**2 - 375*x(1)**2*x(2)*x(3) - 375*x(1)**2*x(2)*x(4) - 95*x(5)*x(1)**2*x(2) - 105*x(1)**2*x(3)**2 - 235*x(1)**2*x(3)*x(4) - 110*x(5)*x(1)**2*x(3) - 105*x(1)**2*x(4)**2 - 110*x(5)*x(1)**2*x(4) + 6148*x(1)*x(2)**3 - 85*x(1)*x(2)**2*x(3) - 85*x(1)*x(2)**2*x(4) - 90*x(5)*x(1)*x(2)**2 + 185*x(1)*x(2)*x(3)**2 + 365*x(1)*x(2)*x(3)*x(4) + 185*x(5)*x(1)*x(2)*x(3) + 185*x(1)*x(2)*x(4)**2 + 185*x(5)*x(1)*x(2)*x(4) + 20*x(1)*x(3)**2*x(4) + 5*x(5)*x(1)*x(3)**2 + 20*x(1)*x(3)*x(4)**2 + 25*x(5)*x(1)*x(3)*x(4) + 5*x(5)*x(1)*x(4)**2 - 1662*x(2)**4 + 200*x(2)**3*x(3) + 200*x(2)**3*x(4) + 100*x(5)*x(2)**3 - 100*x(2)**2*x(3)**2 - 210*x(2)**2*x(3)*x(4) - 105*x(5)*x(2)**2*x(3) - 100*x(2)**2*x(4)**2 - 105*x(5)*x(2)**2*x(4) + 10*x(2)*x(3)**2*x(4) + 5*x(5)*x(2)*x(3)**2 + 10*x(2)*x(3)*x(4)**2 + 15*x(5)*x(2)*x(3)*x(4) + 5*x(5)*x(2)*x(4)**2 - 10*x(3)**2*x(4)**2 - 10*x(5)*x(3)**2*x(4) - 10*x(5)*x(3)*x(4)**2))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(5))*(x(3) - x(5)))

                            Cb(3,5) = (4*(x(1) - x(2))**2*(831*x(1)**8 - 1362*x(1)**7*x(2) - 1762*x(1)**7*x(3) - 1762*x(1)**7*x(4) - 1762*x(1)**7*x(5) + 996*x(1)**6*x(2)**2 + 2514*x(1)**6*x(2)*x(3) + 2514*x(1)**6*x(2)*x(4) + 2514*x(1)**6*x(2)*x(5) + 1081*x(1)**6*x(3)**2 + 3829*x(1)**6*x(3)*x(4) + 3829*x(1)**6*x(3)*x(5) + 1081*x(1)**6*x(4)**2 + 3829*x(1)**6*x(4)*x(5) + 1081*x(1)**6*x(5)**2 - 1722*x(1)**5*x(2)**3 - 270*x(1)**5*x(2)**2*x(3) - 270*x(1)**5*x(2)**2*x(4) - 270*x(1)**5*x(2)**2*x(5) - 932*x(1)**5*x(2)*x(3)**2 - 6340*x(1)**5*x(2)*x(3)*x(4) - 6340*x(1)**5*x(2)*x(3)*x(5) - 932*x(1)**5*x(2)*x(4)**2 - 6340*x(1)**5*x(2)*x(4)*x(5) - 932*x(1)**5*x(2)*x(5)**2 - 200*x(1)**5*x(3)**3 - 2477*x(1)**5*x(3)**2*x(4) - 2477*x(1)**5*x(3)**2*x(5) - 2477*x(1)**5*x(3)*x(4)**2 - 6726*x(1)**5*x(3)*x(4)*x(5) - 2477*x(1)**5*x(3)*x(5)**2 - 200*x(1)**5*x(4)**3 - 2477*x(1)**5*x(4)**2*x(5) - 2477*x(1)**5*x(4)*x(5)**2 - 200*x(1)**5*x(5)**3 + 2694*x(1)**4*x(2)**4 - 722*x(1)**4*x(2)**3*x(3) - 722*x(1)**4*x(2)**3*x(4) - 722*x(1)**4*x(2)**3*x(5) - 961*x(1)**4*x(2)**2*x(3)**2 + 2719*x(1)**4*x(2)**2*x(3)*x(4) + 2719*x(1)**4*x(2)**2*x(3)*x(5) - 961*x(1)**4*x(2)**2*x(4)**2 + 2719*x(1)**4*x(2)**2*x(4)*x(5) - 961*x(1)**4*x(2)**2*x(5)**2 - 230*x(1)**4*x(2)*x(3)**3 + 3636*x(1)**4*x(2)*x(3)**2*x(4) + 3636*x(1)**4*x(2)*x(3)**2*x(5) + 3636*x(1)**4*x(2)*x(3)*x(4)**2 + 11718*x(1)**4*x(2)*x(3)*x(4)*x(5) + 3636*x(1)**4*x(2)*x(3)*x(5)**2 - 230*x(1)**4*x(2)*x(4)**3 + 3636*x(1)**4*x(2)*x(4)**2*x(5) + 3636*x(1)**4*x(2)*x(4)*x(5)**2 - 230*x(1)**4*x(2)*x(5)**3 + 50*x(1)**4*x(3)**4 + 515*x(1)**4*x(3)**3*x(4) + 515*x(1)**4*x(3)**3*x(5) + 1776*x(1)**4*x(3)**2*x(4)**2 + 3652*x(1)**4*x(3)**2*x(4)*x(5) + 1776*x(1)**4*x(3)**2*x(5)**2 + 515*x(1)**4*x(3)*x(4)**3 + 3652*x(1)**4*x(3)*x(4)**2*x(5) + 3652*x(1)**4*x(3)*x(4)*x(5)**2 + 515*x(1)**4*x(3)*x(5)**3 + 50*x(1)**4*x(4)**4 + 515*x(1)**4*x(4)**3*x(5) + 1776*x(1)**4*x(4)**2*x(5)**2 + 515*x(1)**4*x(4)*x(5)**3 + 50*x(1)**4*x(5)**4 - 1722*x(1)**3*x(2)**5 - 722*x(1)**3*x(2)**4*x(3) - 722*x(1)**3*x(2)**4*x(4) - 722*x(1)**3*x(2)**4*x(5) + 2124*x(1)**3*x(2)**3*x(3)**2 + 764*x(1)**3*x(2)**3*x(3)*x(4) + 764*x(1)**3*x(2)**3*x(3)*x(5) + 2124*x(1)**3*x(2)**3*x(4)**2 + 764*x(1)**3*x(2)**3*x(4)*x(5) + 2124*x(1)**3*x(2)**3*x(5)**2 + 310*x(1)**3*x(2)**2*x(3)**3 - 1729*x(1)**3*x(2)**2*x(3)**2*x(4) - 1729*x(1)**3*x(2)**2*x(3)**2*x(5) - 1729*x(1)**3*x(2)**2*x(3)*x(4)**2 - 6252*x(1)**3*x(2)**2*x(3)*x(4)*x(5) - 1729*x(1)**3*x(2)**2*x(3)*x(5)**2 + 310*x(1)**3*x(2)**2*x(4)**3 - 1729*x(1)**3*x(2)**2*x(4)**2*x(5) - 1729*x(1)**3*x(2)**2*x(4)*x(5)**2 + 310*x(1)**3*x(2)**2*x(5)**3 + 10*x(1)**3*x(2)*x(3)**4 + 130*x(1)**3*x(2)*x(3)**3*x(4) + 130*x(1)**3*x(2)*x(3)**3*x(5) - 2874*x(1)**3*x(2)*x(3)**2*x(4)**2 - 5728*x(1)**3*x(2)*x(3)**2*x(4)*x(5) - 2874*x(1)**3*x(2)*x(3)**2*x(5)**2 + 130*x(1)**3*x(2)*x(3)*x(4)**3 - 5728*x(1)**3*x(2)*x(3)*x(4)**2*x(5) - 5728*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 130*x(1)**3*x(2)*x(3)*x(5)**3 + 10*x(1)**3*x(2)*x(4)**4 + 130*x(1)**3*x(2)*x(4)**3*x(5) - 2874*x(1)**3*x(2)*x(4)**2*x(5)**2 + 130*x(1)**3*x(2)*x(4)*x(5)**3 + 10*x(1)**3*x(2)*x(5)**4 - 105*x(1)**3*x(3)**4*x(4) - 105*x(1)**3*x(3)**4*x(5) - 445*x(1)**3*x(3)**3*x(4)**2 - 880*x(1)**3*x(3)**3*x(4)*x(5) - 445*x(1)**3*x(3)**3*x(5)**2 - 445*x(1)**3*x(3)**2*x(4)**3 - 1560*x(1)**3*x(3)**2*x(4)**2*x(5) - 1560*x(1)**3*x(3)**2*x(4)*x(5)**2 - 445*x(1)**3*x(3)**2*x(5)**3 - 105*x(1)**3*x(3)*x(4)**4 - 880*x(1)**3*x(3)*x(4)**3*x(5) - 1560*x(1)**3*x(3)*x(4)**2*x(5)**2 - 880*x(1)**3*x(3)*x(4)*x(5)**3 - 105*x(1)**3*x(3)*x(5)**4 - 105*x(1)**3*x(4)**4*x(5) - 445*x(1)**3*x(4)**3*x(5)**2 - 445*x(1)**3*x(4)**2*x(5)**3 - 105*x(1)**3*x(4)*x(5)**4 + 996*x(1)**2*x(2)**6 - 270*x(1)**2*x(2)**5*x(3) - 270*x(1)**2*x(2)**5*x(4) - 270*x(1)**2*x(2)**5*x(5) - 961*x(1)**2*x(2)**4*x(3)**2 + 2719*x(1)**2*x(2)**4*x(3)*x(4) + 2719*x(1)**2*x(2)**4*x(3)*x(5) - 961*x(1)**2*x(2)**4*x(4)**2 + 2719*x(1)**2*x(2)**4*x(4)*x(5) - 961*x(1)**2*x(2)**4*x(5)**2 + 310*x(1)**2*x(2)**3*x(3)**3 - 1729*x(1)**2*x(2)**3*x(3)**2*x(4) - 1729*x(1)**2*x(2)**3*x(3)**2*x(5) - 1729*x(1)**2*x(2)**3*x(3)*x(4)**2 - 6252*x(1)**2*x(2)**3*x(3)*x(4)*x(5) - 1729*x(1)**2*x(2)**3*x(3)*x(5)**2 + 310*x(1)**2*x(2)**3*x(4)**3 - 1729*x(1)**2*x(2)**3*x(4)**2*x(5) - 1729*x(1)**2*x(2)**3*x(4)*x(5)**2 + 310*x(1)**2*x(2)**3*x(5)**3 - 75*x(1)**2*x(2)**2*x(3)**4 - 780*x(1)**2*x(2)**2*x(3)**3*x(4) - 780*x(1)**2*x(2)**2*x(3)**3*x(5) + 3231*x(1)**2*x(2)**2*x(3)**2*x(4)**2 + 6252*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) + 3231*x(1)**2*x(2)**2*x(3)**2*x(5)**2 - 780*x(1)**2*x(2)**2*x(3)*x(4)**3 + 6252*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) + 6252*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 - 780*x(1)**2*x(2)**2*x(3)*x(5)**3 - 75*x(1)**2*x(2)**2*x(4)**4 - 780*x(1)**2*x(2)**2*x(4)**3*x(5) + 3231*x(1)**2*x(2)**2*x(4)**2*x(5)**2 - 780*x(1)**2*x(2)**2*x(4)*x(5)**3 - 75*x(1)**2*x(2)**2*x(5)**4 + 60*x(1)**2*x(2)*x(3)**4*x(4) + 60*x(1)**2*x(2)*x(3)**4*x(5) + 225*x(1)**2*x(2)*x(3)**3*x(4)**2 + 480*x(1)**2*x(2)*x(3)**3*x(4)*x(5) + 225*x(1)**2*x(2)*x(3)**3*x(5)**2 + 225*x(1)**2*x(2)*x(3)**2*x(4)**3 + 810*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) + 810*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 + 225*x(1)**2*x(2)*x(3)**2*x(5)**3 + 60*x(1)**2*x(2)*x(3)*x(4)**4 + 480*x(1)**2*x(2)*x(3)*x(4)**3*x(5) + 810*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 + 480*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 + 60*x(1)**2*x(2)*x(3)*x(5)**4 + 60*x(1)**2*x(2)*x(4)**4*x(5) + 225*x(1)**2*x(2)*x(4)**3*x(5)**2 + 225*x(1)**2*x(2)*x(4)**2*x(5)**3 + 60*x(1)**2*x(2)*x(4)*x(5)**4 + 65*x(1)**2*x(3)**4*x(4)**2 + 125*x(1)**2*x(3)**4*x(4)*x(5) + 65*x(1)**2*x(3)**4*x(5)**2 + 145*x(1)**2*x(3)**3*x(4)**3 + 415*x(1)**2*x(3)**3*x(4)**2*x(5) + 415*x(1)**2*x(3)**3*x(4)*x(5)**2 + 145*x(1)**2*x(3)**3*x(5)**3 + 65*x(1)**2*x(3)**2*x(4)**4 + 415*x(1)**2*x(3)**2*x(4)**3*x(5) + 690*x(1)**2*x(3)**2*x(4)**2*x(5)**2 + 415*x(1)**2*x(3)**2*x(4)*x(5)**3 + 65*x(1)**2*x(3)**2*x(5)**4 + 125*x(1)**2*x(3)*x(4)**4*x(5) + 415*x(1)**2*x(3)*x(4)**3*x(5)**2 + 415*x(1)**2*x(3)*x(4)**2*x(5)**3 + 125*x(1)**2*x(3)*x(4)*x(5)**4 + 65*x(1)**2*x(4)**4*x(5)**2 + 145*x(1)**2*x(4)**3*x(5)**3 + 65*x(1)**2*x(4)**2*x(5)**4 - 1362*x(1)*x(2)**7 + 2514*x(1)*x(2)**6*x(3) + 2514*x(1)*x(2)**6*x(4) + 2514*x(1)*x(2)**6*x(5) - 932*x(1)*x(2)**5*x(3)**2 - 6340*x(1)*x(2)**5*x(3)*x(4) - 6340*x(1)*x(2)**5*x(3)*x(5) - 932*x(1)*x(2)**5*x(4)**2 - 6340*x(1)*x(2)**5*x(4)*x(5) - 932*x(1)*x(2)**5*x(5)**2 - 230*x(1)*x(2)**4*x(3)**3 + 3636*x(1)*x(2)**4*x(3)**2*x(4) + 3636*x(1)*x(2)**4*x(3)**2*x(5) + 3636*x(1)*x(2)**4*x(3)*x(4)**2 + 11718*x(1)*x(2)**4*x(3)*x(4)*x(5) + 3636*x(1)*x(2)**4*x(3)*x(5)**2 - 230*x(1)*x(2)**4*x(4)**3 + 3636*x(1)*x(2)**4*x(4)**2*x(5) + 3636*x(1)*x(2)**4*x(4)*x(5)**2 - 230*x(1)*x(2)**4*x(5)**3 + 10*x(1)*x(2)**3*x(3)**4 + 130*x(1)*x(2)**3*x(3)**3*x(4) + 130*x(1)*x(2)**3*x(3)**3*x(5) - 2874*x(1)*x(2)**3*x(3)**2*x(4)**2 - 5728*x(1)*x(2)**3*x(3)**2*x(4)*x(5) - 2874*x(1)*x(2)**3*x(3)**2*x(5)**2 + 130*x(1)*x(2)**3*x(3)*x(4)**3 - 5728*x(1)*x(2)**3*x(3)*x(4)**2*x(5) - 5728*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 + 130*x(1)*x(2)**3*x(3)*x(5)**3 + 10*x(1)*x(2)**3*x(4)**4 + 130*x(1)*x(2)**3*x(4)**3*x(5) - 2874*x(1)*x(2)**3*x(4)**2*x(5)**2 + 130*x(1)*x(2)**3*x(4)*x(5)**3 + 10*x(1)*x(2)**3*x(5)**4 + 60*x(1)*x(2)**2*x(3)**4*x(4) + 60*x(1)*x(2)**2*x(3)**4*x(5) + 225*x(1)*x(2)**2*x(3)**3*x(4)**2 + 480*x(1)*x(2)**2*x(3)**3*x(4)*x(5) + 225*x(1)*x(2)**2*x(3)**3*x(5)**2 + 225*x(1)*x(2)**2*x(3)**2*x(4)**3 + 810*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) + 810*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 + 225*x(1)*x(2)**2*x(3)**2*x(5)**3 + 60*x(1)*x(2)**2*x(3)*x(4)**4 + 480*x(1)*x(2)**2*x(3)*x(4)**3*x(5) + 810*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 + 480*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 + 60*x(1)*x(2)**2*x(3)*x(5)**4 + 60*x(1)*x(2)**2*x(4)**4*x(5) + 225*x(1)*x(2)**2*x(4)**3*x(5)**2 + 225*x(1)*x(2)**2*x(4)**2*x(5)**3 + 60*x(1)*x(2)**2*x(4)*x(5)**4 - 55*x(1)*x(2)*x(3)**4*x(4)**2 - 130*x(1)*x(2)*x(3)**4*x(4)*x(5) - 55*x(1)*x(2)*x(3)**4*x(5)**2 - 110*x(1)*x(2)*x(3)**3*x(4)**3 - 350*x(1)*x(2)*x(3)**3*x(4)**2*x(5) - 350*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 - 110*x(1)*x(2)*x(3)**3*x(5)**3 - 55*x(1)*x(2)*x(3)**2*x(4)**4 - 350*x(1)*x(2)*x(3)**2*x(4)**3*x(5) - 570*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 - 350*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 - 55*x(1)*x(2)*x(3)**2*x(5)**4 - 130*x(1)*x(2)*x(3)*x(4)**4*x(5) - 350*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 - 350*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 - 130*x(1)*x(2)*x(3)*x(4)*x(5)**4 - 55*x(1)*x(2)*x(4)**4*x(5)**2 - 110*x(1)*x(2)*x(4)**3*x(5)**3 - 55*x(1)*x(2)*x(4)**2*x(5)**4 - 15*x(1)*x(3)**4*x(4)**3 - 30*x(1)*x(3)**4*x(4)**2*x(5) - 30*x(1)*x(3)**4*x(4)*x(5)**2 - 15*x(1)*x(3)**4*x(5)**3 - 15*x(1)*x(3)**3*x(4)**4 - 60*x(1)*x(3)**3*x(4)**3*x(5) - 90*x(1)*x(3)**3*x(4)**2*x(5)**2 - 60*x(1)*x(3)**3*x(4)*x(5)**3 - 15*x(1)*x(3)**3*x(5)**4 - 30*x(1)*x(3)**2*x(4)**4*x(5) - 90*x(1)*x(3)**2*x(4)**3*x(5)**2 - 90*x(1)*x(3)**2*x(4)**2*x(5)**3 - 30*x(1)*x(3)**2*x(4)*x(5)**4 - 30*x(1)*x(3)*x(4)**4*x(5)**2 - 60*x(1)*x(3)*x(4)**3*x(5)**3 - 30*x(1)*x(3)*x(4)**2*x(5)**4 - 15*x(1)*x(4)**4*x(5)**3 - 15*x(1)*x(4)**3*x(5)**4 + 831*x(2)**8 - 1762*x(2)**7*x(3) - 1762*x(2)**7*x(4) - 1762*x(2)**7*x(5) + 1081*x(2)**6*x(3)**2 + 3829*x(2)**6*x(3)*x(4) + 3829*x(2)**6*x(3)*x(5) + 1081*x(2)**6*x(4)**2 + 3829*x(2)**6*x(4)*x(5) + 1081*x(2)**6*x(5)**2 - 200*x(2)**5*x(3)**3 - 2477*x(2)**5*x(3)**2*x(4) - 2477*x(2)**5*x(3)**2*x(5) - 2477*x(2)**5*x(3)*x(4)**2 - 6726*x(2)**5*x(3)*x(4)*x(5) - 2477*x(2)**5*x(3)*x(5)**2 - 200*x(2)**5*x(4)**3 - 2477*x(2)**5*x(4)**2*x(5) - 2477*x(2)**5*x(4)*x(5)**2 - 200*x(2)**5*x(5)**3 + 50*x(2)**4*x(3)**4 + 515*x(2)**4*x(3)**3*x(4) + 515*x(2)**4*x(3)**3*x(5) + 1776*x(2)**4*x(3)**2*x(4)**2 + 3652*x(2)**4*x(3)**2*x(4)*x(5) + 1776*x(2)**4*x(3)**2*x(5)**2 + 515*x(2)**4*x(3)*x(4)**3 + 3652*x(2)**4*x(3)*x(4)**2*x(5) + 3652*x(2)**4*x(3)*x(4)*x(5)**2 + 515*x(2)**4*x(3)*x(5)**3 + 50*x(2)**4*x(4)**4 + 515*x(2)**4*x(4)**3*x(5) + 1776*x(2)**4*x(4)**2*x(5)**2 + 515*x(2)**4*x(4)*x(5)**3 + 50*x(2)**4*x(5)**4 - 105*x(2)**3*x(3)**4*x(4) - 105*x(2)**3*x(3)**4*x(5) - 445*x(2)**3*x(3)**3*x(4)**2 - 880*x(2)**3*x(3)**3*x(4)*x(5) - 445*x(2)**3*x(3)**3*x(5)**2 - 445*x(2)**3*x(3)**2*x(4)**3 - 1560*x(2)**3*x(3)**2*x(4)**2*x(5) - 1560*x(2)**3*x(3)**2*x(4)*x(5)**2 - 445*x(2)**3*x(3)**2*x(5)**3 - 105*x(2)**3*x(3)*x(4)**4 - 880*x(2)**3*x(3)*x(4)**3*x(5) - 1560*x(2)**3*x(3)*x(4)**2*x(5)**2 - 880*x(2)**3*x(3)*x(4)*x(5)**3 - 105*x(2)**3*x(3)*x(5)**4 - 105*x(2)**3*x(4)**4*x(5) - 445*x(2)**3*x(4)**3*x(5)**2 - 445*x(2)**3*x(4)**2*x(5)**3 - 105*x(2)**3*x(4)*x(5)**4 + 65*x(2)**2*x(3)**4*x(4)**2 + 125*x(2)**2*x(3)**4*x(4)*x(5) + 65*x(2)**2*x(3)**4*x(5)**2 + 145*x(2)**2*x(3)**3*x(4)**3 + 415*x(2)**2*x(3)**3*x(4)**2*x(5) + 415*x(2)**2*x(3)**3*x(4)*x(5)**2 + 145*x(2)**2*x(3)**3*x(5)**3 + 65*x(2)**2*x(3)**2*x(4)**4 + 415*x(2)**2*x(3)**2*x(4)**3*x(5) + 690*x(2)**2*x(3)**2*x(4)**2*x(5)**2 + 415*x(2)**2*x(3)**2*x(4)*x(5)**3 + 65*x(2)**2*x(3)**2*x(5)**4 + 125*x(2)**2*x(3)*x(4)**4*x(5) + 415*x(2)**2*x(3)*x(4)**3*x(5)**2 + 415*x(2)**2*x(3)*x(4)**2*x(5)**3 + 125*x(2)**2*x(3)*x(4)*x(5)**4 + 65*x(2)**2*x(4)**4*x(5)**2 + 145*x(2)**2*x(4)**3*x(5)**3 + 65*x(2)**2*x(4)**2*x(5)**4 - 15*x(2)*x(3)**4*x(4)**3 - 30*x(2)*x(3)**4*x(4)**2*x(5) - 30*x(2)*x(3)**4*x(4)*x(5)**2 - 15*x(2)*x(3)**4*x(5)**3 - 15*x(2)*x(3)**3*x(4)**4 - 60*x(2)*x(3)**3*x(4)**3*x(5) - 90*x(2)*x(3)**3*x(4)**2*x(5)**2 - 60*x(2)*x(3)**3*x(4)*x(5)**3 - 15*x(2)*x(3)**3*x(5)**4 - 30*x(2)*x(3)**2*x(4)**4*x(5) - 90*x(2)*x(3)**2*x(4)**3*x(5)**2 - 90*x(2)*x(3)**2*x(4)**2*x(5)**3 - 30*x(2)*x(3)**2*x(4)*x(5)**4 - 30*x(2)*x(3)*x(4)**4*x(5)**2 - 60*x(2)*x(3)*x(4)**3*x(5)**3 - 30*x(2)*x(3)*x(4)**2*x(5)**4 - 15*x(2)*x(4)**4*x(5)**3 - 15*x(2)*x(4)**3*x(5)**4 + 5*x(3)**4*x(4)**4 + 10*x(3)**4*x(4)**3*x(5) + 15*x(3)**4*x(4)**2*x(5)**2 + 10*x(3)**4*x(4)*x(5)**3 + 5*x(3)**4*x(5)**4 + 10*x(3)**3*x(4)**4*x(5) + 20*x(3)**3*x(4)**3*x(5)**2 + 20*x(3)**3*x(4)**2*x(5)**3 + 10*x(3)**3*x(4)*x(5)**4 + 15*x(3)**2*x(4)**4*x(5)**2 + 20*x(3)**2*x(4)**3*x(5)**3 + 15*x(3)**2*x(4)**2*x(5)**4 + 10*x(3)*x(4)**4*x(5)**3 + 10*x(3)*x(4)**3*x(5)**4 + 5*x(4)**4*x(5)**4))/(5*(x(1) - x(3))**2*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2)

                            Cb(3,6) = (4*(x(1) - x(2))**2*(- 1662*x(1)**7*x(2) - 1662*x(1)**7*x(3) + 1662*x(1)**7*x(4) + 1662*x(1)**7*x(5) + 4486*x(1)**6*x(2)**2 + 4586*x(1)**6*x(2)*x(3) - 962*x(1)**6*x(2)*x(4) - 962*x(1)**6*x(2)*x(5) + 1862*x(1)**6*x(3)**2 + 1662*x(1)**6*x(3)*x(4) + 1662*x(1)**6*x(3)*x(5) - 3524*x(1)**6*x(4)**2 - 5286*x(1)**6*x(4)*x(5) - 3524*x(1)**6*x(5)**2 - 4606*x(1)**5*x(2)**3 - 4596*x(1)**5*x(2)**2*x(3) - 4251*x(1)**5*x(2)**2*x(4) - 4251*x(1)**5*x(2)**2*x(5) - 2604*x(1)**5*x(2)*x(3)**2 - 6558*x(1)**5*x(2)*x(3)*x(4) - 6558*x(1)**5*x(2)*x(3)*x(5) + 6695*x(1)**5*x(2)*x(4)**2 + 7442*x(1)**5*x(2)*x(4)*x(5) + 6695*x(1)**5*x(2)*x(5)**2 - 300*x(1)**5*x(3)**3 - 3834*x(1)**5*x(3)**2*x(4) - 3834*x(1)**5*x(3)**2*x(5) + 1972*x(1)**5*x(3)*x(4)**2 + 310*x(1)**5*x(3)*x(4)*x(5) + 1972*x(1)**5*x(3)*x(5)**2 + 2162*x(1)**5*x(4)**3 + 5991*x(1)**5*x(4)**2*x(5) + 5991*x(1)**5*x(4)*x(5)**2 + 2162*x(1)**5*x(5)**3 + 3504*x(1)**4*x(2)**4 + 1652*x(1)**4*x(2)**3*x(3) + 3681*x(1)**4*x(2)**3*x(4) + 3681*x(1)**4*x(2)**3*x(5) - 1792*x(1)**4*x(2)**2*x(3)**2 + 10804*x(1)**4*x(2)**2*x(3)*x(4) + 10804*x(1)**4*x(2)**2*x(3)*x(5) - 2844*x(1)**4*x(2)**2*x(4)**2 + 5096*x(1)**4*x(2)**2*x(4)*x(5) - 2844*x(1)**4*x(2)**2*x(5)**2 - 340*x(1)**4*x(2)*x(3)**3 + 8812*x(1)**4*x(2)*x(3)**2*x(4) + 8812*x(1)**4*x(2)*x(3)**2*x(5) - 3816*x(1)**4*x(2)*x(3)*x(4)**2 + 1190*x(1)**4*x(2)*x(3)*x(4)*x(5) - 3816*x(1)**4*x(2)*x(3)*x(5)**2 - 3941*x(1)**4*x(2)*x(4)**3 - 12148*x(1)**4*x(2)*x(4)**2*x(5) - 12148*x(1)**4*x(2)*x(4)*x(5)**2 - 3941*x(1)**4*x(2)*x(5)**3 + 100*x(1)**4*x(3)**4 + 720*x(1)**4*x(3)**3*x(4) + 720*x(1)**4*x(3)**3*x(5) + 1972*x(1)**4*x(3)**2*x(4)**2 + 4254*x(1)**4*x(3)**2*x(4)*x(5) + 1972*x(1)**4*x(3)**2*x(5)**2 - 2392*x(1)**4*x(3)*x(4)**3 - 2812*x(1)**4*x(3)*x(4)**2*x(5) - 2812*x(1)**4*x(3)*x(4)*x(5)**2 - 2392*x(1)**4*x(3)*x(5)**3 - 400*x(1)**4*x(4)**4 - 2877*x(1)**4*x(4)**3*x(5) - 3182*x(1)**4*x(4)**2*x(5)**2 - 2877*x(1)**4*x(4)*x(5)**3 - 400*x(1)**4*x(5)**4 - 4606*x(1)**3*x(2)**5 + 1652*x(1)**3*x(2)**4*x(3) + 3681*x(1)**3*x(2)**4*x(4) + 3681*x(1)**3*x(2)**4*x(5) + 5378*x(1)**3*x(2)**3*x(3)**2 - 11986*x(1)**3*x(2)**3*x(3)*x(4) - 11986*x(1)**3*x(2)**3*x(3)*x(5) - 1104*x(1)**3*x(2)**3*x(4)**2 - 15254*x(1)**3*x(2)**3*x(4)*x(5) - 1104*x(1)**3*x(2)**3*x(5)**2 + 490*x(1)**3*x(2)**2*x(3)**3 - 5218*x(1)**3*x(2)**2*x(3)**2*x(4) - 5218*x(1)**3*x(2)**2*x(3)**2*x(5) + 2144*x(1)**3*x(2)**2*x(3)*x(4)**2 - 1110*x(1)**3*x(2)**2*x(3)*x(4)*x(5) + 2144*x(1)**3*x(2)**2*x(3)*x(5)**2 + 1974*x(1)**3*x(2)**2*x(4)**3 + 6622*x(1)**3*x(2)**2*x(4)**2*x(5) + 6622*x(1)**3*x(2)**2*x(4)*x(5)**2 + 1974*x(1)**3*x(2)**2*x(5)**3 + 20*x(1)**3*x(2)*x(3)**4 + 150*x(1)**3*x(2)*x(3)**3*x(4) + 150*x(1)**3*x(2)*x(3)**3*x(5) - 6338*x(1)**3*x(2)*x(3)**2*x(4)**2 - 12586*x(1)**3*x(2)*x(3)**2*x(4)*x(5) - 6338*x(1)**3*x(2)*x(3)**2*x(5)**2 + 5998*x(1)**3*x(2)*x(3)*x(4)**3 + 5658*x(1)**3*x(2)*x(3)*x(4)**2*x(5) + 5658*x(1)**3*x(2)*x(3)*x(4)*x(5)**2 + 5998*x(1)**3*x(2)*x(3)*x(5)**3 - 45*x(1)**3*x(2)*x(4)**4 + 5998*x(1)**3*x(2)*x(4)**3*x(5) + 5848*x(1)**3*x(2)*x(4)**2*x(5)**2 + 5998*x(1)**3*x(2)*x(4)*x(5)**3 - 45*x(1)**3*x(2)*x(5)**4 - 210*x(1)**3*x(3)**4*x(4) - 210*x(1)**3*x(3)**4*x(5) - 550*x(1)**3*x(3)**3*x(4)**2 - 1090*x(1)**3*x(3)**3*x(4)*x(5) - 550*x(1)**3*x(3)**3*x(5)**2 + 130*x(1)**3*x(3)**2*x(4)**3 - 290*x(1)**3*x(3)**2*x(4)**2*x(5) - 290*x(1)**3*x(3)**2*x(4)*x(5)**2 + 130*x(1)**3*x(3)**2*x(5)**3 + 530*x(1)**3*x(3)*x(4)**4 + 1190*x(1)**3*x(3)*x(4)**3*x(5) + 1300*x(1)**3*x(3)*x(4)**2*x(5)**2 + 1190*x(1)**3*x(3)*x(4)*x(5)**3 + 530*x(1)**3*x(3)*x(5)**4 + 100*x(1)**3*x(4)**5 + 615*x(1)**3*x(4)**4*x(5) + 930*x(1)**3*x(4)**3*x(5)**2 + 930*x(1)**3*x(4)**2*x(5)**3 + 615*x(1)**3*x(4)*x(5)**4 + 100*x(1)**3*x(5)**5 + 4486*x(1)**2*x(2)**6 - 4596*x(1)**2*x(2)**5*x(3) - 4251*x(1)**2*x(2)**5*x(4) - 4251*x(1)**2*x(2)**5*x(5) - 1792*x(1)**2*x(2)**4*x(3)**2 + 10804*x(1)**2*x(2)**4*x(3)*x(4) + 10804*x(1)**2*x(2)**4*x(3)*x(5) - 2844*x(1)**2*x(2)**4*x(4)**2 + 5096*x(1)**2*x(2)**4*x(4)*x(5) - 2844*x(1)**2*x(2)**4*x(5)**2 + 490*x(1)**2*x(2)**3*x(3)**3 - 5218*x(1)**2*x(2)**3*x(3)**2*x(4) - 5218*x(1)**2*x(2)**3*x(3)**2*x(5) + 2144*x(1)**2*x(2)**3*x(3)*x(4)**2 - 1110*x(1)**2*x(2)**3*x(3)*x(4)*x(5) + 2144*x(1)**2*x(2)**3*x(3)*x(5)**2 + 1974*x(1)**2*x(2)**3*x(4)**3 + 6622*x(1)**2*x(2)**3*x(4)**2*x(5) + 6622*x(1)**2*x(2)**3*x(4)*x(5)**2 + 1974*x(1)**2*x(2)**3*x(5)**3 - 150*x(1)**2*x(2)**2*x(3)**4 - 1170*x(1)**2*x(2)**2*x(3)**3*x(4) - 1170*x(1)**2*x(2)**2*x(3)**3*x(5) + 8832*x(1)**2*x(2)**2*x(3)**2*x(4)**2 + 17154*x(1)**2*x(2)**2*x(3)**2*x(4)*x(5) + 8832*x(1)**2*x(2)**2*x(3)**2*x(5)**2 - 7872*x(1)**2*x(2)**2*x(3)*x(4)**3 - 6912*x(1)**2*x(2)**2*x(3)*x(4)**2*x(5) - 6912*x(1)**2*x(2)**2*x(3)*x(4)*x(5)**2 - 7872*x(1)**2*x(2)**2*x(3)*x(5)**3 + 720*x(1)**2*x(2)**2*x(4)**4 - 6852*x(1)**2*x(2)**2*x(4)**3*x(5) - 6132*x(1)**2*x(2)**2*x(4)**2*x(5)**2 - 6852*x(1)**2*x(2)**2*x(4)*x(5)**3 + 720*x(1)**2*x(2)**2*x(5)**4 + 120*x(1)**2*x(2)*x(3)**4*x(4) + 120*x(1)**2*x(2)*x(3)**4*x(5) + 350*x(1)**2*x(2)*x(3)**3*x(4)**2 + 710*x(1)**2*x(2)*x(3)**3*x(4)*x(5) + 350*x(1)**2*x(2)*x(3)**3*x(5)**2 - 10*x(1)**2*x(2)*x(3)**2*x(4)**3 + 330*x(1)**2*x(2)*x(3)**2*x(4)**2*x(5) + 330*x(1)**2*x(2)*x(3)**2*x(4)*x(5)**2 - 10*x(1)**2*x(2)*x(3)**2*x(5)**3 - 370*x(1)**2*x(2)*x(3)*x(4)**4 - 750*x(1)**2*x(2)*x(3)*x(4)**3*x(5) - 780*x(1)**2*x(2)*x(3)*x(4)**2*x(5)**2 - 750*x(1)**2*x(2)*x(3)*x(4)*x(5)**3 - 370*x(1)**2*x(2)*x(3)*x(5)**4 - 85*x(1)**2*x(2)*x(4)**5 - 510*x(1)**2*x(2)*x(4)**4*x(5) - 750*x(1)**2*x(2)*x(4)**3*x(5)**2 - 750*x(1)**2*x(2)*x(4)**2*x(5)**3 - 510*x(1)**2*x(2)*x(4)*x(5)**4 - 85*x(1)**2*x(2)*x(5)**5 + 130*x(1)**2*x(3)**4*x(4)**2 + 250*x(1)**2*x(3)**4*x(4)*x(5) + 130*x(1)**2*x(3)**4*x(5)**2 + 130*x(1)**2*x(3)**3*x(4)**3 + 390*x(1)**2*x(3)**3*x(4)**2*x(5) + 390*x(1)**2*x(3)**3*x(4)*x(5)**2 + 130*x(1)**2*x(3)**3*x(5)**3 - 150*x(1)**2*x(3)**2*x(4)**4 - 170*x(1)**2*x(3)**2*x(4)**3*x(5) - 60*x(1)**2*x(3)**2*x(4)**2*x(5)**2 - 170*x(1)**2*x(3)**2*x(4)*x(5)**3 - 150*x(1)**2*x(3)**2*x(5)**4 - 110*x(1)**2*x(3)*x(4)**5 - 370*x(1)**2*x(3)*x(4)**4*x(5) - 500*x(1)**2*x(3)*x(4)**3*x(5)**2 - 500*x(1)**2*x(3)*x(4)**2*x(5)**3 - 370*x(1)**2*x(3)*x(4)*x(5)**4 - 110*x(1)**2*x(3)*x(5)**5 - 105*x(1)**2*x(4)**5*x(5) - 220*x(1)**2*x(4)**4*x(5)**2 - 220*x(1)**2*x(4)**3*x(5)**3 - 220*x(1)**2*x(4)**2*x(5)**4 - 105*x(1)**2*x(4)*x(5)**5 - 1662*x(1)*x(2)**7 + 4586*x(1)*x(2)**6*x(3) - 962*x(1)*x(2)**6*x(4) - 962*x(1)*x(2)**6*x(5) - 2604*x(1)*x(2)**5*x(3)**2 - 6558*x(1)*x(2)**5*x(3)*x(4) - 6558*x(1)*x(2)**5*x(3)*x(5) + 6695*x(1)*x(2)**5*x(4)**2 + 7442*x(1)*x(2)**5*x(4)*x(5) + 6695*x(1)*x(2)**5*x(5)**2 - 340*x(1)*x(2)**4*x(3)**3 + 8812*x(1)*x(2)**4*x(3)**2*x(4) + 8812*x(1)*x(2)**4*x(3)**2*x(5) - 3816*x(1)*x(2)**4*x(3)*x(4)**2 + 1190*x(1)*x(2)**4*x(3)*x(4)*x(5) - 3816*x(1)*x(2)**4*x(3)*x(5)**2 - 3941*x(1)*x(2)**4*x(4)**3 - 12148*x(1)*x(2)**4*x(4)**2*x(5) - 12148*x(1)*x(2)**4*x(4)*x(5)**2 - 3941*x(1)*x(2)**4*x(5)**3 + 20*x(1)*x(2)**3*x(3)**4 + 150*x(1)*x(2)**3*x(3)**3*x(4) + 150*x(1)*x(2)**3*x(3)**3*x(5) - 6338*x(1)*x(2)**3*x(3)**2*x(4)**2 - 12586*x(1)*x(2)**3*x(3)**2*x(4)*x(5) - 6338*x(1)*x(2)**3*x(3)**2*x(5)**2 + 5998*x(1)*x(2)**3*x(3)*x(4)**3 + 5658*x(1)*x(2)**3*x(3)*x(4)**2*x(5) + 5658*x(1)*x(2)**3*x(3)*x(4)*x(5)**2 + 5998*x(1)*x(2)**3*x(3)*x(5)**3 - 45*x(1)*x(2)**3*x(4)**4 + 5998*x(1)*x(2)**3*x(4)**3*x(5) + 5848*x(1)*x(2)**3*x(4)**2*x(5)**2 + 5998*x(1)*x(2)**3*x(4)*x(5)**3 - 45*x(1)*x(2)**3*x(5)**4 + 120*x(1)*x(2)**2*x(3)**4*x(4) + 120*x(1)*x(2)**2*x(3)**4*x(5) + 350*x(1)*x(2)**2*x(3)**3*x(4)**2 + 710*x(1)*x(2)**2*x(3)**3*x(4)*x(5) + 350*x(1)*x(2)**2*x(3)**3*x(5)**2 - 10*x(1)*x(2)**2*x(3)**2*x(4)**3 + 330*x(1)*x(2)**2*x(3)**2*x(4)**2*x(5) + 330*x(1)*x(2)**2*x(3)**2*x(4)*x(5)**2 - 10*x(1)*x(2)**2*x(3)**2*x(5)**3 - 370*x(1)*x(2)**2*x(3)*x(4)**4 - 750*x(1)*x(2)**2*x(3)*x(4)**3*x(5) - 780*x(1)*x(2)**2*x(3)*x(4)**2*x(5)**2 - 750*x(1)*x(2)**2*x(3)*x(4)*x(5)**3 - 370*x(1)*x(2)**2*x(3)*x(5)**4 - 85*x(1)*x(2)**2*x(4)**5 - 510*x(1)*x(2)**2*x(4)**4*x(5) - 750*x(1)*x(2)**2*x(4)**3*x(5)**2 - 750*x(1)*x(2)**2*x(4)**2*x(5)**3 - 510*x(1)*x(2)**2*x(4)*x(5)**4 - 85*x(1)*x(2)**2*x(5)**5 - 110*x(1)*x(2)*x(3)**4*x(4)**2 - 260*x(1)*x(2)*x(3)**4*x(4)*x(5) - 110*x(1)*x(2)*x(3)**4*x(5)**2 - 170*x(1)*x(2)*x(3)**3*x(4)**3 - 450*x(1)*x(2)*x(3)**3*x(4)**2*x(5) - 450*x(1)*x(2)*x(3)**3*x(4)*x(5)**2 - 170*x(1)*x(2)*x(3)**3*x(5)**3 + 120*x(1)*x(2)*x(3)**2*x(4)**4 + 70*x(1)*x(2)*x(3)**2*x(4)**3*x(5) - 90*x(1)*x(2)*x(3)**2*x(4)**2*x(5)**2 + 70*x(1)*x(2)*x(3)**2*x(4)*x(5)**3 + 120*x(1)*x(2)*x(3)**2*x(5)**4 + 160*x(1)*x(2)*x(3)*x(4)**5 + 440*x(1)*x(2)*x(3)*x(4)**4*x(5) + 550*x(1)*x(2)*x(3)*x(4)**3*x(5)**2 + 550*x(1)*x(2)*x(3)*x(4)**2*x(5)**3 + 440*x(1)*x(2)*x(3)*x(4)*x(5)**4 + 160*x(1)*x(2)*x(3)*x(5)**5 + 180*x(1)*x(2)*x(4)**5*x(5) + 350*x(1)*x(2)*x(4)**4*x(5)**2 + 350*x(1)*x(2)*x(4)**3*x(5)**3 + 350*x(1)*x(2)*x(4)**2*x(5)**4 + 180*x(1)*x(2)*x(4)*x(5)**5 - 30*x(1)*x(3)**4*x(4)**3 - 60*x(1)*x(3)**4*x(4)**2*x(5) - 60*x(1)*x(3)**4*x(4)*x(5)**2 - 30*x(1)*x(3)**4*x(5)**3 + 10*x(1)*x(3)**3*x(4)**4 - 10*x(1)*x(3)**3*x(4)**3*x(5) - 30*x(1)*x(3)**3*x(4)**2*x(5)**2 - 10*x(1)*x(3)**3*x(4)*x(5)**3 + 10*x(1)*x(3)**3*x(5)**4 + 20*x(1)*x(3)**2*x(4)**5 + 50*x(1)*x(3)**2*x(4)**4*x(5) + 50*x(1)*x(3)**2*x(4)**3*x(5)**2 + 50*x(1)*x(3)**2*x(4)**2*x(5)**3 + 50*x(1)*x(3)**2*x(4)*x(5)**4 + 20*x(1)*x(3)**2*x(5)**5 + 20*x(1)*x(3)*x(4)**5*x(5) + 50*x(1)*x(3)*x(4)**4*x(5)**2 + 50*x(1)*x(3)*x(4)**3*x(5)**3 + 50*x(1)*x(3)*x(4)**2*x(5)**4 + 20*x(1)*x(3)*x(4)*x(5)**5 + 5*x(1)*x(4)**5*x(5)**2 + 5*x(1)*x(4)**4*x(5)**3 + 5*x(1)*x(4)**3*x(5)**4 + 5*x(1)*x(4)**2*x(5)**5 - 1662*x(2)**7*x(3) + 1662*x(2)**7*x(4) + 1662*x(2)**7*x(5) + 1862*x(2)**6*x(3)**2 + 1662*x(2)**6*x(3)*x(4) + 1662*x(2)**6*x(3)*x(5) - 3524*x(2)**6*x(4)**2 - 5286*x(2)**6*x(4)*x(5) - 3524*x(2)**6*x(5)**2 - 300*x(2)**5*x(3)**3 - 3834*x(2)**5*x(3)**2*x(4) - 3834*x(2)**5*x(3)**2*x(5) + 1972*x(2)**5*x(3)*x(4)**2 + 310*x(2)**5*x(3)*x(4)*x(5) + 1972*x(2)**5*x(3)*x(5)**2 + 2162*x(2)**5*x(4)**3 + 5991*x(2)**5*x(4)**2*x(5) + 5991*x(2)**5*x(4)*x(5)**2 + 2162*x(2)**5*x(5)**3 + 100*x(2)**4*x(3)**4 + 720*x(2)**4*x(3)**3*x(4) + 720*x(2)**4*x(3)**3*x(5) + 1972*x(2)**4*x(3)**2*x(4)**2 + 4254*x(2)**4*x(3)**2*x(4)*x(5) + 1972*x(2)**4*x(3)**2*x(5)**2 - 2392*x(2)**4*x(3)*x(4)**3 - 2812*x(2)**4*x(3)*x(4)**2*x(5) - 2812*x(2)**4*x(3)*x(4)*x(5)**2 - 2392*x(2)**4*x(3)*x(5)**3 - 400*x(2)**4*x(4)**4 - 2877*x(2)**4*x(4)**3*x(5) - 3182*x(2)**4*x(4)**2*x(5)**2 - 2877*x(2)**4*x(4)*x(5)**3 - 400*x(2)**4*x(5)**4 - 210*x(2)**3*x(3)**4*x(4) - 210*x(2)**3*x(3)**4*x(5) - 550*x(2)**3*x(3)**3*x(4)**2 - 1090*x(2)**3*x(3)**3*x(4)*x(5) - 550*x(2)**3*x(3)**3*x(5)**2 + 130*x(2)**3*x(3)**2*x(4)**3 - 290*x(2)**3*x(3)**2*x(4)**2*x(5) - 290*x(2)**3*x(3)**2*x(4)*x(5)**2 + 130*x(2)**3*x(3)**2*x(5)**3 + 530*x(2)**3*x(3)*x(4)**4 + 1190*x(2)**3*x(3)*x(4)**3*x(5) + 1300*x(2)**3*x(3)*x(4)**2*x(5)**2 + 1190*x(2)**3*x(3)*x(4)*x(5)**3 + 530*x(2)**3*x(3)*x(5)**4 + 100*x(2)**3*x(4)**5 + 615*x(2)**3*x(4)**4*x(5) + 930*x(2)**3*x(4)**3*x(5)**2 + 930*x(2)**3*x(4)**2*x(5)**3 + 615*x(2)**3*x(4)*x(5)**4 + 100*x(2)**3*x(5)**5 + 130*x(2)**2*x(3)**4*x(4)**2 + 250*x(2)**2*x(3)**4*x(4)*x(5) + 130*x(2)**2*x(3)**4*x(5)**2 + 130*x(2)**2*x(3)**3*x(4)**3 + 390*x(2)**2*x(3)**3*x(4)**2*x(5) + 390*x(2)**2*x(3)**3*x(4)*x(5)**2 + 130*x(2)**2*x(3)**3*x(5)**3 - 150*x(2)**2*x(3)**2*x(4)**4 - 170*x(2)**2*x(3)**2*x(4)**3*x(5) - 60*x(2)**2*x(3)**2*x(4)**2*x(5)**2 - 170*x(2)**2*x(3)**2*x(4)*x(5)**3 - 150*x(2)**2*x(3)**2*x(5)**4 - 110*x(2)**2*x(3)*x(4)**5 - 370*x(2)**2*x(3)*x(4)**4*x(5) - 500*x(2)**2*x(3)*x(4)**3*x(5)**2 - 500*x(2)**2*x(3)*x(4)**2*x(5)**3 - 370*x(2)**2*x(3)*x(4)*x(5)**4 - 110*x(2)**2*x(3)*x(5)**5 - 105*x(2)**2*x(4)**5*x(5) - 220*x(2)**2*x(4)**4*x(5)**2 - 220*x(2)**2*x(4)**3*x(5)**3 - 220*x(2)**2*x(4)**2*x(5)**4 - 105*x(2)**2*x(4)*x(5)**5 - 30*x(2)*x(3)**4*x(4)**3 - 60*x(2)*x(3)**4*x(4)**2*x(5) - 60*x(2)*x(3)**4*x(4)*x(5)**2 - 30*x(2)*x(3)**4*x(5)**3 + 10*x(2)*x(3)**3*x(4)**4 - 10*x(2)*x(3)**3*x(4)**3*x(5) - 30*x(2)*x(3)**3*x(4)**2*x(5)**2 - 10*x(2)*x(3)**3*x(4)*x(5)**3 + 10*x(2)*x(3)**3*x(5)**4 + 20*x(2)*x(3)**2*x(4)**5 + 50*x(2)*x(3)**2*x(4)**4*x(5) + 50*x(2)*x(3)**2*x(4)**3*x(5)**2 + 50*x(2)*x(3)**2*x(4)**2*x(5)**3 + 50*x(2)*x(3)**2*x(4)*x(5)**4 + 20*x(2)*x(3)**2*x(5)**5 + 20*x(2)*x(3)*x(4)**5*x(5) + 50*x(2)*x(3)*x(4)**4*x(5)**2 + 50*x(2)*x(3)*x(4)**3*x(5)**3 + 50*x(2)*x(3)*x(4)**2*x(5)**4 + 20*x(2)*x(3)*x(4)*x(5)**5 + 5*x(2)*x(4)**5*x(5)**2 + 5*x(2)*x(4)**4*x(5)**3 + 5*x(2)*x(4)**3*x(5)**4 + 5*x(2)*x(4)**2*x(5)**5 + 10*x(3)**4*x(4)**4 + 20*x(3)**4*x(4)**3*x(5) + 30*x(3)**4*x(4)**2*x(5)**2 + 20*x(3)**4*x(4)*x(5)**3 + 10*x(3)**4*x(5)**4 - 10*x(3)**3*x(4)**5 - 10*x(3)**3*x(4)**4*x(5) - 10*x(3)**3*x(4)**3*x(5)**2 - 10*x(3)**3*x(4)**2*x(5)**3 - 10*x(3)**3*x(4)*x(5)**4 - 10*x(3)**3*x(5)**5 - 10*x(3)**2*x(4)**5*x(5) - 10*x(3)**2*x(4)**4*x(5)**2 - 10*x(3)**2*x(4)**3*x(5)**3 - 10*x(3)**2*x(4)**2*x(5)**4 - 10*x(3)**2*x(4)*x(5)**5 - 10*x(3)*x(4)**5*x(5)**2 - 10*x(3)*x(4)**4*x(5)**3 - 10*x(3)*x(4)**3*x(5)**4 - 10*x(3)*x(4)**2*x(5)**5))/(5*(x(1) - x(3))*(x(1) - x(4))**2*(x(1) - x(5))**2*(x(2) - x(4))**2*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(3,7) = (4*(x(1) - x(2))**2*(1662*x(1)**6 - 4486*x(1)**5*x(2) - 1862*x(1)**5*x(3) - 1862*x(1)**5*x(4) - 1762*x(1)**5*x(5) + 4606*x(1)**4*x(2)**2 + 4371*x(1)**4*x(2)*x(3) + 4371*x(1)**4*x(2)*x(4) + 4476*x(1)**4*x(2)*x(5) + 300*x(1)**4*x(3)**2 + 2272*x(1)**4*x(3)*x(4) + 2067*x(1)**4*x(3)*x(5) + 300*x(1)**4*x(4)**2 + 2067*x(1)**4*x(4)*x(5) + 100*x(1)**4*x(5)**2 - 3504*x(1)**3*x(2)**3 - 2579*x(1)**3*x(2)**2*x(3) - 2579*x(1)**3*x(2)**2*x(4) - 2754*x(1)**3*x(2)**2*x(5) + 30*x(1)**3*x(2)*x(3)**2 - 6178*x(1)**3*x(2)*x(3)*x(4) - 6208*x(1)**3*x(2)*x(3)*x(5) + 30*x(1)**3*x(2)*x(4)**2 - 6208*x(1)**3*x(2)*x(4)*x(5) + 10*x(1)**3*x(2)*x(5)**2 - 100*x(1)**3*x(3)**3 - 520*x(1)**3*x(3)**2*x(4) - 410*x(1)**3*x(3)**2*x(5) - 520*x(1)**3*x(3)*x(4)**2 - 830*x(1)**3*x(3)*x(4)*x(5) - 205*x(1)**3*x(3)*x(5)**2 - 100*x(1)**3*x(4)**3 - 410*x(1)**3*x(4)**2*x(5) - 205*x(1)**3*x(4)*x(5)**2 + 4606*x(1)**2*x(2)**4 - 2579*x(1)**2*x(2)**3*x(3) - 2579*x(1)**2*x(2)**3*x(4) - 2754*x(1)**2*x(2)**3*x(5) - 550*x(1)**2*x(2)**2*x(3)**2 + 8122*x(1)**2*x(2)**2*x(3)*x(4) + 8452*x(1)**2*x(2)**2*x(3)*x(5) - 550*x(1)**2*x(2)**2*x(4)**2 + 8452*x(1)**2*x(2)**2*x(4)*x(5) - 190*x(1)**2*x(2)**2*x(5)**2 + 85*x(1)**2*x(2)*x(3)**3 + 405*x(1)**2*x(2)*x(3)**2*x(4) + 350*x(1)**2*x(2)*x(3)**2*x(5) + 405*x(1)**2*x(2)*x(3)*x(4)**2 + 670*x(1)**2*x(2)*x(3)*x(4)*x(5) + 175*x(1)**2*x(2)*x(3)*x(5)**2 + 85*x(1)**2*x(2)*x(4)**3 + 350*x(1)**2*x(2)*x(4)**2*x(5) + 175*x(1)**2*x(2)*x(4)*x(5)**2 + 110*x(1)**2*x(3)**3*x(4) + 105*x(1)**2*x(3)**3*x(5) + 240*x(1)**2*x(3)**2*x(4)**2 + 345*x(1)**2*x(3)**2*x(4)*x(5) + 110*x(1)**2*x(3)**2*x(5)**2 + 110*x(1)**2*x(3)*x(4)**3 + 345*x(1)**2*x(3)*x(4)**2*x(5) + 220*x(1)**2*x(3)*x(4)*x(5)**2 + 105*x(1)**2*x(4)**3*x(5) + 110*x(1)**2*x(4)**2*x(5)**2 - 4486*x(1)*x(2)**5 + 4371*x(1)*x(2)**4*x(3) + 4371*x(1)*x(2)**4*x(4) + 4476*x(1)*x(2)**4*x(5) + 30*x(1)*x(2)**3*x(3)**2 - 6178*x(1)*x(2)**3*x(3)*x(4) - 6208*x(1)*x(2)**3*x(3)*x(5) + 30*x(1)*x(2)**3*x(4)**2 - 6208*x(1)*x(2)**3*x(4)*x(5) + 10*x(1)*x(2)**3*x(5)**2 + 85*x(1)*x(2)**2*x(3)**3 + 405*x(1)*x(2)**2*x(3)**2*x(4) + 350*x(1)*x(2)**2*x(3)**2*x(5) + 405*x(1)*x(2)**2*x(3)*x(4)**2 + 670*x(1)*x(2)**2*x(3)*x(4)*x(5) + 175*x(1)*x(2)**2*x(3)*x(5)**2 + 85*x(1)*x(2)**2*x(4)**3 + 350*x(1)*x(2)**2*x(4)**2*x(5) + 175*x(1)*x(2)**2*x(4)*x(5)**2 - 160*x(1)*x(2)*x(3)**3*x(4) - 180*x(1)*x(2)*x(3)**3*x(5) - 320*x(1)*x(2)*x(3)**2*x(4)**2 - 500*x(1)*x(2)*x(3)**2*x(4)*x(5) - 180*x(1)*x(2)*x(3)**2*x(5)**2 - 160*x(1)*x(2)*x(3)*x(4)**3 - 500*x(1)*x(2)*x(3)*x(4)**2*x(5) - 340*x(1)*x(2)*x(3)*x(4)*x(5)**2 - 180*x(1)*x(2)*x(4)**3*x(5) - 180*x(1)*x(2)*x(4)**2*x(5)**2 - 20*x(1)*x(3)**3*x(4)**2 - 20*x(1)*x(3)**3*x(4)*x(5) - 5*x(1)*x(3)**3*x(5)**2 - 20*x(1)*x(3)**2*x(4)**3 - 40*x(1)*x(3)**2*x(4)**2*x(5) - 25*x(1)*x(3)**2*x(4)*x(5)**2 - 20*x(1)*x(3)*x(4)**3*x(5) - 25*x(1)*x(3)*x(4)**2*x(5)**2 - 5*x(1)*x(4)**3*x(5)**2 + 1662*x(2)**6 - 1862*x(2)**5*x(3) - 1862*x(2)**5*x(4) - 1762*x(2)**5*x(5) + 300*x(2)**4*x(3)**2 + 2272*x(2)**4*x(3)*x(4) + 2067*x(2)**4*x(3)*x(5) + 300*x(2)**4*x(4)**2 + 2067*x(2)**4*x(4)*x(5) + 100*x(2)**4*x(5)**2 - 100*x(2)**3*x(3)**3 - 520*x(2)**3*x(3)**2*x(4) - 410*x(2)**3*x(3)**2*x(5) - 520*x(2)**3*x(3)*x(4)**2 - 830*x(2)**3*x(3)*x(4)*x(5) - 205*x(2)**3*x(3)*x(5)**2 - 100*x(2)**3*x(4)**3 - 410*x(2)**3*x(4)**2*x(5) - 205*x(2)**3*x(4)*x(5)**2 + 110*x(2)**2*x(3)**3*x(4) + 105*x(2)**2*x(3)**3*x(5) + 240*x(2)**2*x(3)**2*x(4)**2 + 345*x(2)**2*x(3)**2*x(4)*x(5) + 110*x(2)**2*x(3)**2*x(5)**2 + 110*x(2)**2*x(3)*x(4)**3 + 345*x(2)**2*x(3)*x(4)**2*x(5) + 220*x(2)**2*x(3)*x(4)*x(5)**2 + 105*x(2)**2*x(4)**3*x(5) + 110*x(2)**2*x(4)**2*x(5)**2 - 20*x(2)*x(3)**3*x(4)**2 - 20*x(2)*x(3)**3*x(4)*x(5) - 5*x(2)*x(3)**3*x(5)**2 - 20*x(2)*x(3)**2*x(4)**3 - 40*x(2)*x(3)**2*x(4)**2*x(5) - 25*x(2)*x(3)**2*x(4)*x(5)**2 - 20*x(2)*x(3)*x(4)**3*x(5) - 25*x(2)*x(3)*x(4)**2*x(5)**2 - 5*x(2)*x(4)**3*x(5)**2 + 10*x(3)**3*x(4)**3 + 10*x(3)**3*x(4)**2*x(5) + 10*x(3)**3*x(4)*x(5)**2 + 10*x(3)**2*x(4)**3*x(5) + 10*x(3)**2*x(4)**2*x(5)**2 + 10*x(3)*x(4)**3*x(5)**2))/(5*(x(1) - x(3))*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5)))

                            Cb(3,8) = (24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(1) - x(2))**6 + (x(1)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(1) - x(2)))/5 - (x(2)**5*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(1) - x(2)))/5 + x(1)*(x(1) - x(2))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(2)*(x(1) - x(2))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + (x(1)**3*(x(1) - x(2))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 - (x(2)**3*(x(1) - x(2))*(((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 + 2*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))))/3 + (x(1)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(1) - x(2))**3)/3 - (x(2)**3*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2*(x(1) - x(2))**3)/3 + x(1)*(x(1) - x(2))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(2)*(x(1) - x(2))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))**2 - x(1)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1) - x(2))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(2)**2*(24/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (24*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1) - x(2))**3*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - x(1)**2*(x(1) - x(2))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) + x(2)**2*(x(1) - x(2))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*((2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(2)*x(3) + 2*x(1)*x(5) + 2*x(2)*x(5) + 2*x(3)*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(2*x(1)*x(2) + 2*x(1)*x(3) + 2*x(1)*x(4) + 2*x(2)*x(3) + 2*x(2)*x(4) + 2*x(3)*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))) - (x(1)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1) - x(2))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2 + (x(2)**4*(12/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - (12*(x(3) - x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5))))*(x(1) - x(2))*((6*x(1) + 6*x(2) + 6*x(3) + 6*x(5))/((x(1) - x(4))*(x(2) - x(4))*(x(4) - x(5))) - ((x(3) - x(4))*(6*x(1) + 6*x(2) + 6*x(3) + 6*x(4)))/((x(1) - x(5))*(x(2) - x(5))*(x(3) - x(5))*(x(4) - x(5)))))/2

                            Cb(3,9) = (4*(x(1) - x(2))**2*(- 1662*x(1)**5*x(2) - 1662*x(1)**5*x(3) + 1662*x(1)**5*x(4) + 1662*x(1)**5*x(5) + 6248*x(1)**4*x(2)**2 + 4686*x(1)**4*x(2)*x(3) - 4386*x(1)**4*x(2)*x(4) - 4486*x(1)**4*x(2)*x(5) + 200*x(1)**4*x(3)**2 + 1662*x(1)**4*x(3)*x(4) + 1562*x(1)**4*x(3)*x(5) - 1862*x(1)**4*x(4)**2 - 1862*x(1)**4*x(4)*x(5) - 1762*x(1)**4*x(5)**2 - 9182*x(1)**3*x(2)**3 - 3029*x(1)**3*x(2)**2*x(3) + 2744*x(1)**3*x(2)**2*x(4) + 2839*x(1)**3*x(2)**2*x(5) + 15*x(1)**3*x(2)*x(3)**2 - 6358*x(1)**3*x(2)*x(3)*x(4) - 6358*x(1)**3*x(2)*x(3)*x(5) + 6138*x(1)**3*x(2)*x(4)**2 + 6138*x(1)**3*x(2)*x(4)*x(5) + 6243*x(1)**3*x(2)*x(5)**2 - 100*x(1)**3*x(3)**3 - 310*x(1)**3*x(3)**2*x(4) - 205*x(1)**3*x(3)**2*x(5) + 110*x(1)**3*x(3)*x(4)**2 + 110*x(1)**3*x(3)*x(4)*x(5) + 205*x(1)**3*x(3)*x(5)**2 + 300*x(1)**3*x(4)**3 + 300*x(1)**3*x(4)**2*x(5) + 300*x(1)**3*x(4)*x(5)**2 + 100*x(1)**3*x(5)**3 + 6248*x(1)**2*x(2)**4 - 3029*x(1)**2*x(2)**3*x(3) + 2744*x(1)**2*x(2)**3*x(4) + 2839*x(1)**2*x(2)**3*x(5) - 380*x(1)**2*x(2)**2*x(3)**2 + 9372*x(1)**2*x(2)**2*x(3)*x(4) + 9562*x(1)**2*x(2)**2*x(3)*x(5) - 8612*x(1)**2*x(2)**2*x(4)**2 - 8612*x(1)**2*x(2)**2*x(4)*x(5) - 8992*x(1)**2*x(2)**2*x(5)**2 + 85*x(1)**2*x(2)*x(3)**3 + 270*x(1)**2*x(2)*x(3)**2*x(4) + 190*x(1)**2*x(2)*x(3)**2*x(5) - 70*x(1)**2*x(2)*x(3)*x(4)**2 - 70*x(1)**2*x(2)*x(3)*x(4)*x(5) - 180*x(1)**2*x(2)*x(3)*x(5)**2 - 280*x(1)**2*x(2)*x(4)**3 - 280*x(1)**2*x(2)*x(4)**2*x(5) - 280*x(1)**2*x(2)*x(4)*x(5)**2 - 95*x(1)**2*x(2)*x(5)**3 + 110*x(1)**2*x(3)**3*x(4) + 105*x(1)**2*x(3)**3*x(5) + 110*x(1)**2*x(3)**2*x(4)**2 + 110*x(1)**2*x(3)**2*x(4)*x(5) - 120*x(1)**2*x(3)*x(4)**3 - 120*x(1)**2*x(3)*x(4)**2*x(5) - 120*x(1)**2*x(3)*x(4)*x(5)**2 - 105*x(1)**2*x(3)*x(5)**3 - 100*x(1)**2*x(4)**4 - 100*x(1)**2*x(4)**3*x(5) - 100*x(1)**2*x(4)**2*x(5)**2 - 100*x(1)**2*x(4)*x(5)**3 - 1662*x(1)*x(2)**5 + 4686*x(1)*x(2)**4*x(3) - 4386*x(1)*x(2)**4*x(4) - 4486*x(1)*x(2)**4*x(5) + 15*x(1)*x(2)**3*x(3)**2 - 6358*x(1)*x(2)**3*x(3)*x(4) - 6358*x(1)*x(2)**3*x(3)*x(5) + 6138*x(1)*x(2)**3*x(4)**2 + 6138*x(1)*x(2)**3*x(4)*x(5) + 6243*x(1)*x(2)**3*x(5)**2 + 85*x(1)*x(2)**2*x(3)**3 + 270*x(1)*x(2)**2*x(3)**2*x(4) + 190*x(1)*x(2)**2*x(3)**2*x(5) - 70*x(1)*x(2)**2*x(3)*x(4)**2 - 70*x(1)*x(2)**2*x(3)*x(4)*x(5) - 180*x(1)*x(2)**2*x(3)*x(5)**2 - 280*x(1)*x(2)**2*x(4)**3 - 280*x(1)*x(2)**2*x(4)**2*x(5) - 280*x(1)*x(2)**2*x(4)*x(5)**2 - 95*x(1)*x(2)**2*x(5)**3 - 160*x(1)*x(2)*x(3)**3*x(4) - 180*x(1)*x(2)*x(3)**3*x(5) - 200*x(1)*x(2)*x(3)**2*x(4)**2 - 200*x(1)*x(2)*x(3)**2*x(4)*x(5) - 10*x(1)*x(2)*x(3)**2*x(5)**2 + 170*x(1)*x(2)*x(3)*x(4)**3 + 170*x(1)*x(2)*x(3)*x(4)**2*x(5) + 170*x(1)*x(2)*x(3)*x(4)*x(5)**2 + 190*x(1)*x(2)*x(3)*x(5)**3 + 190*x(1)*x(2)*x(4)**4 + 190*x(1)*x(2)*x(4)**3*x(5) + 190*x(1)*x(2)*x(4)**2*x(5)**2 + 190*x(1)*x(2)*x(4)*x(5)**3 - 20*x(1)*x(3)**3*x(4)**2 - 20*x(1)*x(3)**3*x(4)*x(5) - 5*x(1)*x(3)**3*x(5)**2 + 10*x(1)*x(3)**2*x(4)**3 + 10*x(1)*x(3)**2*x(4)**2*x(5) + 10*x(1)*x(3)**2*x(4)*x(5)**2 + 5*x(1)*x(3)**2*x(5)**3 + 10*x(1)*x(3)*x(4)**4 + 10*x(1)*x(3)*x(4)**3*x(5) + 10*x(1)*x(3)*x(4)**2*x(5)**2 + 10*x(1)*x(3)*x(4)*x(5)**3 - 1662*x(2)**5*x(3) + 1662*x(2)**5*x(4) + 1662*x(2)**5*x(5) + 200*x(2)**4*x(3)**2 + 1662*x(2)**4*x(3)*x(4) + 1562*x(2)**4*x(3)*x(5) - 1862*x(2)**4*x(4)**2 - 1862*x(2)**4*x(4)*x(5) - 1762*x(2)**4*x(5)**2 - 100*x(2)**3*x(3)**3 - 310*x(2)**3*x(3)**2*x(4) - 205*x(2)**3*x(3)**2*x(5) + 110*x(2)**3*x(3)*x(4)**2 + 110*x(2)**3*x(3)*x(4)*x(5) + 205*x(2)**3*x(3)*x(5)**2 + 300*x(2)**3*x(4)**3 + 300*x(2)**3*x(4)**2*x(5) + 300*x(2)**3*x(4)*x(5)**2 + 100*x(2)**3*x(5)**3 + 110*x(2)**2*x(3)**3*x(4) + 105*x(2)**2*x(3)**3*x(5) + 110*x(2)**2*x(3)**2*x(4)**2 + 110*x(2)**2*x(3)**2*x(4)*x(5) - 120*x(2)**2*x(3)*x(4)**3 - 120*x(2)**2*x(3)*x(4)**2*x(5) - 120*x(2)**2*x(3)*x(4)*x(5)**2 - 105*x(2)**2*x(3)*x(5)**3 - 100*x(2)**2*x(4)**4 - 100*x(2)**2*x(4)**3*x(5) - 100*x(2)**2*x(4)**2*x(5)**2 - 100*x(2)**2*x(4)*x(5)**3 - 20*x(2)*x(3)**3*x(4)**2 - 20*x(2)*x(3)**3*x(4)*x(5) - 5*x(2)*x(3)**3*x(5)**2 + 10*x(2)*x(3)**2*x(4)**3 + 10*x(2)*x(3)**2*x(4)**2*x(5) + 10*x(2)*x(3)**2*x(4)*x(5)**2 + 5*x(2)*x(3)**2*x(5)**3 + 10*x(2)*x(3)*x(4)**4 + 10*x(2)*x(3)*x(4)**3*x(5) + 10*x(2)*x(3)*x(4)**2*x(5)**2 + 10*x(2)*x(3)*x(4)*x(5)**3 + 10*x(3)**3*x(4)**3 + 10*x(3)**3*x(4)**2*x(5) + 10*x(3)**3*x(4)*x(5)**2 - 10*x(3)**2*x(4)**4 - 10*x(3)**2*x(4)**3*x(5) - 10*x(3)**2*x(4)**2*x(5)**2 - 10*x(3)**2*x(4)*x(5)**3))/(5*(x(1) - x(4))*(x(1) - x(5))**2*(x(2) - x(4))*(x(2) - x(5))**2*(x(3) - x(5))**2)

                            Cb(3,10) = (4*(x(1) - x(2))**2*(831*x(1)**4 - 3124*x(1)**3*x(2) - 100*x(1)**3*x(3) - 100*x(1)**3*x(4) + 4591*x(1)**2*x(2)**2 + 95*x(1)**2*x(2)*x(3) + 95*x(1)**2*x(2)*x(4) + 50*x(1)**2*x(3)**2 + 105*x(1)**2*x(3)*x(4) + 50*x(1)**2*x(4)**2 - 3124*x(1)*x(2)**3 + 95*x(1)*x(2)**2*x(3) + 95*x(1)*x(2)**2*x(4) - 95*x(1)*x(2)*x(3)**2 - 190*x(1)*x(2)*x(3)*x(4) - 95*x(1)*x(2)*x(4)**2 - 5*x(1)*x(3)**2*x(4) - 5*x(1)*x(3)*x(4)**2 + 831*x(2)**4 - 100*x(2)**3*x(3) - 100*x(2)**3*x(4) + 50*x(2)**2*x(3)**2 + 105*x(2)**2*x(3)*x(4) + 50*x(2)**2*x(4)**2 - 5*x(2)*x(3)**2*x(4) - 5*x(2)*x(3)*x(4)**2 + 5*x(3)**2*x(4)**2))/(5*(x(1) - x(5))**2*(x(2) - x(5))**2*(x(3) - x(5))**2)

                        end do

                    else ! TENO
                        ! (Fu, et al., 2016) Table 2 (for right flux)
                        d_cbL_${XYZ}$ (0, :) = 18d0/35d0
                        d_cbL_${XYZ}$ (1, :) = 3d0/35d0
                        d_cbL_${XYZ}$ (2, :) = 9d0/35d0
                        d_cbL_${XYZ}$ (3, :) = 1d0/35d0
                        d_cbL_${XYZ}$ (4, :) = 4d0/35d0

                        d_cbR_${XYZ}$ (0, :) = 18d0/35d0
                        d_cbR_${XYZ}$ (1, :) = 9d0/35d0
                        d_cbR_${XYZ}$ (2, :) = 3d0/35d0
                        d_cbR_${XYZ}$ (3, :) = 4d0/35d0
                        d_cbR_${XYZ}$ (4, :) = 1d0/35d0

                    end if
                end if

            end if
        #:endfor

! END: Computing WENO Coefficients ================================
        if (weno_dir == 1) then
            !$acc update device(poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x)
        elseif (weno_dir == 2) then
            !$acc update device(poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y)
        else
            !$acc update device(poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z)
        end if

        ! Nullifying WENO coefficients and cell-boundary locations pointers

        nullify (s_cb)

    end subroutine s_compute_weno_coefficients

    subroutine s_weno(v_vf, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, &
                      norm_dir, weno_dir, &
                      is1_weno_d, is2_weno_d, is3_weno_d)

        type(scalar_field), dimension(1:), intent(in) :: v_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(inout) :: vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in) :: norm_dir
        integer, intent(in) :: weno_dir
        type(int_bounds_info), intent(in) :: is1_weno_d, is2_weno_d, is3_weno_d

        real(kind(0d0)), dimension(-weno_polyn:weno_polyn - 1) :: dvd
        real(kind(0d0)), dimension(0:weno_num_stencils) :: poly
        real(kind(0d0)), dimension(0:weno_num_stencils) :: alpha
        real(kind(0d0)), dimension(0:weno_num_stencils) :: omega
        real(kind(0d0)), dimension(0:weno_num_stencils) :: beta
        real(kind(0d0)), dimension(0:weno_num_stencils) :: delta
        real(kind(0d0)), dimension(-3:3) :: v ! temporary field value array for clarity (WENO7 only)
        real(kind(0d0)) :: tau

        integer :: i, j, k, l

        is1_weno = is1_weno_d
        is2_weno = is2_weno_d
        is3_weno = is3_weno_d

        !$acc update device(is1_weno, is2_weno, is3_weno)

        if (weno_order /= 1) then
            call s_initialize_weno(v_vf, &
                                   norm_dir, weno_dir)
        end if

        if (weno_order == 1) then
            if (weno_dir == 1) then
                !$acc parallel loop collapse(4) default(present)
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
                !$acc end parallel loop
            else if (weno_dir == 2) then
                !$acc parallel loop collapse(4) default(present)
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
                !$acc end parallel loop
            else if (weno_dir == 3) then
                !$acc parallel loop collapse(4) default(present)
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
                !$acc end parallel loop
            end if
        elseif (weno_order == 3) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    !$acc parallel loop collapse(4) gang vector default(present) private(beta,dvd,poly,omega,alpha,tau)
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
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1d0 + d_cbL_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Borges, et al. (2008)
                                        tau = abs(beta(1) - beta(0))
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1d0 + tau/beta)

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
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1d0 + d_cbR_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        alpha = d_cbR_${XYZ}$ (:, j)*(1d0 + tau/beta)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                end if
            #:endfor
        elseif (weno_order == 5) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    !$acc parallel loop vector gang collapse(3) default(present) private(dvd, poly, beta, alpha, omega, tau, delta)
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                !$acc loop seq
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
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1d0 + d_cbL_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Borges, et al. (2008)
                                        tau = abs(beta(2) - beta(0))                   ! Equation 25
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1d0 + tau/beta)  ! Equation 28 (note: weno_eps was already added to beta)

                                    elseif (teno) then
                                        ! Fu, et al. (2016)
                                        ! Fu''s code: https://dx.doi.org/10.13140/RG.2.2.36250.34247
                                        tau = abs(beta(2) - beta(0))
                                        alpha = (1d0 + tau/beta)**6d0              ! Equation 22 (reuse alpha as gamma; pick C=1 & q=6)
                                        omega = alpha/sum(alpha)                    ! Equation 25 (reuse omega as xi)
                                        delta = merge(0d0, 1d0, omega < teno_CT)    ! Equation 26
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
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1d0 + d_cbR_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        alpha = d_cbR_${XYZ}$ (:, j)*(1d0 + tau/beta)

                                    elseif (teno) then
                                        alpha = delta*d_cbR_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    if (mp_weno) then
                        call s_preserve_monotonicity(v_rs_ws_${XYZ}$, vL_rs_vf_${XYZ}$, &
                                                     vR_rs_vf_${XYZ}$)
                    end if
                end if
            #:endfor
        elseif (weno_order == 7) then
            #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (weno_dir == ${WENO_DIR}$) then
                    !$acc parallel loop vector gang collapse(3) default(present) private(poly, beta, alpha, omega, tau, delta, v)
                    ! Note: dvd is not used as the equations are not cast in terms of the differences
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                !$acc loop seq
                                do i = 1, v_size

                                    v = v_rs_ws_${XYZ}$ (j - 3:j + 3, k, l, i) ! temporary field value array for clarity

                                    if (.not. teno) then
                                        ! (Balsara & Shu, 2000) Page 11 Table I
                                        ! poly(0) = ( 1d0*v(-3) -  5d0*v(-2) + 13d0*v(-1) + 3d0*v( 0)) / 12d0 !&
                                        ! poly(1) = (-1d0*v(-2) +  7d0*v(-1) +  7d0*v( 0) - 1d0*v( 1)) / 12d0 !&
                                        ! poly(2) = ( 3d0*v(-1) + 13d0*v( 0) -  5d0*v( 1) + 1d0*v( 2)) / 12d0 !&
                                        ! poly(3) = (25d0*v( 0) - 23d0*v( 1) + 13d0*v( 2) - 3d0*v( 3)) / 12d0 !&
                                        poly(0) = ( CpL(0,1)*v(-3) + CpL(0,2)*v(-2) + CpL(0,3)*v(-1) + CpL(0,4)*v( 0)) !&
                                        poly(1) = ( CpL(1,1)*v(-2) + CpL(1,2)*v(-1) + CpL(1,3)*v( 0) + CpL(1,4)*v( 1)) !&
                                        poly(2) = ( CpL(2,1)*v(-1) + CpL(2,2)*v( 0) + CpL(2,3)*v( 1) + CpL(2,4)*v( 2)) !&
                                        poly(3) = ( CpL(3,1)*v( 0) + CpL(3,2)*v( 1) + CpL(3,3)*v( 2) + CpL(3,4)*v( 3)) !&
                                    else
                                        ! (Fu, et al., 2016) Table 1
                                        ! Note: Unlike TENO5, TENO7 stencils differ from WENO7 stencils
                                        ! See Figure 2 (right) for right-sided flux (at i+1/2)
                                        ! Here we need the left-sided flux, so we flip the weights with respect to the x=i point
                                        ! But we need to keep the stencil order to reuse the beta coefficients
                                        poly(0) = ( 2d0*v(-1) +  5d0*v( 0) -  1d0*v( 1)) / 6d0 !&
                                        poly(1) = (11d0*v( 0) -  7d0*v( 1) +  2d0*v( 2)) / 6d0 !&
                                        poly(2) = (-1d0*v(-2) +  5d0*v(-1) +  2d0*v( 0)) / 6d0 !&
                                        poly(3) = (25d0*v( 0) - 23d0*v( 1) + 13d0*v( 2) - 3d0*v( 3)) / 12d0 !&
                                        poly(4) = ( 1d0*v(-3) -  5d0*v(-2) + 13d0*v(-1) + 3d0*v( 0)) / 12d0 !&
                                    end if

                                    if (.not. teno) then
                                        ! (Balsara & Shu, 2000) Page 11 Section III.a
                                        ! Note: parentheses are needed to group the terms before '+ weno_eps' to avoid unintended floating point errors
                                        ! beta(0) = ( v(-3)*(547d0*v(-3) - 3882d0*v(-2) +  4642d0*v(-1) - 1854d0*v( 0)) & !&
                                        !           + v(-2)*(              7043d0*v(-2) - 17246d0*v(-1) + 7042d0*v( 0)) & !&
                                        !           + v(-1)*(                             11003d0*v(-1) - 9402d0*v( 0)) & !&
                                        !           + v( 0)*(                                             2107d0*v( 0)) ) & !&
                                        !           + weno_eps !&

                                        ! beta(1) = ( v(-2)*(267d0*v(-2) - 1642d0*v(-1) + 1602d0*v( 0) -  494d0*v( 1)) & !&
                                        !           + v(-1)*(              2843d0*v(-1) - 5966d0*v( 0) + 1922d0*v( 1)) & !&
                                        !           + v( 0)*(                             3443d0*v( 0) - 2522d0*v( 1)) & !&
                                        !           + v( 1)*(                                             547d0*v( 1)) ) & !&
                                        !           + weno_eps !&

                                        ! beta(2) = ( v(-1)*(547d0*v(-1) - 2522d0*v( 0) + 1922d0*v( 1) -  494d0*v( 2)) & !&
                                        !           + v( 0)*(              3443d0*v( 0) - 5966d0*v( 1) + 1602d0*v( 2)) & !&
                                        !           + v( 1)*(                             2843d0*v( 1) - 1642d0*v( 2)) & !&
                                        !           + v( 2)*(                                             267d0*v( 2)) ) & !&
                                        !           + weno_eps !&

                                        ! beta(3) = ( v( 0)*(2107d0*v( 0) - 9402d0*v( 1) +  7042d0*v( 2) - 1854d0*v( 3)) & !&
                                        !           + v( 1)*(              11003d0*v( 1) - 17246d0*v( 2) + 4642d0*v( 3)) & !&
                                        !           + v( 2)*(                               7043d0*v( 2) - 3882d0*v( 3)) & !&
                                        !           + v( 3)*(                                               547d0*v( 3)) ) & !&
                                        !           + weno_eps !&

                                        beta(0) = ( Cb(0,1)*v(-3)*v(-3) + Cb(0,2)*v(-3)*v(-2) + Cb(0,3)*v(-3)*v(-1) + Cb(0,4)*v(-3)*v( 0) & !&
                                                  + Cb(0,5)*v(-2)*v(-2) + Cb(0,6)*v(-2)*v(-1) + Cb(0,7)*v(-2)*v( 0) & !&
                                                  + Cb(0,8)*v(-1)*v(-1) + Cb(0,9)*v(-1)*v( 0) & !&
                                                  + Cb(0,10)*v( 0)*v( 0) ) & !&
                                                  + weno_eps !&

                                        beta(1) = ( Cb(1,1)*v(-2)*v(-2) + Cb(1,2)*v(-2)*v(-1) + Cb(1,3)*v(-2)*v( 0) + Cb(1,4)*v(-2)*v( 1) & !&
                                                  + Cb(1,5)*v(-1)*v(-1) + Cb(1,6)*v(-1)*v( 0) + Cb(1,7)*v(-1)*v( 1) & !&
                                                  + Cb(1,8)*v( 0)*v( 0) + Cb(1,9)*v( 0)*v( 1) & !&
                                                  + Cb(1,10)*v( 1)*v( 1) ) & !&
                                                  + weno_eps !&

                                        beta(2) = ( Cb(2,1)*v(-1)*v(-1) + Cb(2,2)*v(-1)*v( 0) + Cb(2,3)*v(-1)*v( 1) + Cb(2,4)*v(-1)*v( 2) & !&
                                                  + Cb(2,5)*v( 0)*v( 0) + Cb(2,6)*v( 0)*v( 1) + Cb(2,7)*v( 0)*v( 2) & !&
                                                  + Cb(2,8)*v( 1)*v( 1) + Cb(2,9)*v( 1)*v( 2) & !&
                                                  + Cb(2,10)*v( 2)*v( 2) ) & !&
                                                  + weno_eps !&

                                        beta(3) = ( Cb(3,1)*v( 0)*v( 0) + Cb(3,2)*v( 0)*v( 1) + Cb(3,3)*v( 0)*v( 2) + Cb(3,4)*v( 0)*v( 3) & !&
                                                  + Cb(3,5)*v( 1)*v( 1) + Cb(3,6)*v( 1)*v( 2) + Cb(3,7)*v( 1)*v( 3) & !&
                                                  + Cb(3,8)*v( 2)*v( 2) + Cb(3,9)*v( 2)*v( 3) & !&
                                                  + Cb(3,10)*v( 3)*v( 3) ) & !&
                                                  + weno_eps !&

                                    else ! TENO
                                        ! High-Order Low-Dissipation Targeted ENO Schemes for Ideal Magnetohydrodynamics (Fu & Tang, 2019) Section 3.2
                                        beta(0) = 13d0/12d0*(v(-1) - 2d0*v( 0) + v( 1))**2d0 + ((    v(-1)             -     v( 1))**2d0)/4d0 + weno_eps !&
                                        beta(1) = 13d0/12d0*(v( 0) - 2d0*v( 1) + v( 2))**2d0 + ((3d0*v( 0) - 4d0*v( 1) +     v( 2))**2d0)/4d0 + weno_eps !&
                                        beta(2) = 13d0/12d0*(v(-2) - 2d0*v(-1) + v( 0))**2d0 + ((    v(-2) - 4d0*v(-1) + 3d0*v( 0))**2d0)/4d0 + weno_eps !&

                                        beta(3) = ( v( 0)*(2107d0*v( 0) - 9402d0*v( 1) +  7042d0*v( 2) - 1854d0*v( 3)) & !&
                                                  + v( 1)*(              11003d0*v( 1) - 17246d0*v( 2) + 4642d0*v( 3)) & !&
                                                  + v( 2)*(                               7043d0*v( 2) - 3882d0*v( 3)) & !&
                                                  + v( 3)*(                                               547d0*v( 3)) ) / 240d0 & !&
                                                  + weno_eps !&

                                        beta(4) = ( v(-3)*(547d0*v(-3) - 3882d0*v(-2) +  4642d0*v(-1) - 1854d0*v( 0)) & !&
                                                  + v(-2)*(              7043d0*v(-2) - 17246d0*v(-1) + 7042d0*v( 0)) & !&
                                                  + v(-1)*(                             11003d0*v(-1) - 9402d0*v( 0)) & !&
                                                  + v( 0)*(                                             2107d0*v( 0)) ) / 240d0 & !&
                                                  + weno_eps !&
                                    end if

                                    if (wenojs) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbL_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbL_${XYZ}$ (:, j)*(1d0 + d_cbL_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbL_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbL_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        ! Castro, et al. (2010)
                                        ! Don & Borges (2013) also helps
                                        tau = abs(beta(3) - beta(0)) ! Equation 50
                                        alpha = d_cbL_${XYZ}$ (:, j)*(1d0 + (tau/beta)**wenoz_q) ! q = 2,3,4 for stability

                                    elseif (teno) then
                                        tau = abs(beta(4) - beta(3)) ! Note the reordering of stencils
                                        alpha = (1d0 + tau/beta)**6d0
                                        omega = alpha/sum(alpha)
                                        delta = merge(0d0, 1d0, omega < teno_CT)
                                        alpha = delta*d_cbL_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                    if (.not. teno) then
                                        ! poly(0) = (-3d0*v(-3) + 13d0*v(-2) - 23d0*v(-1) + 25d0*v( 0)) / 12d0 !&
                                        ! poly(1) = ( 1d0*v(-2) -  5d0*v(-1) + 13d0*v( 0) +  3d0*v( 1)) / 12d0 !&
                                        ! poly(2) = (-1d0*v(-1) +  7d0*v( 0) +  7d0*v( 1) -  1d0*v( 2)) / 12d0 !&
                                        ! poly(3) = ( 3d0*v( 0) + 13d0*v( 1) -  5d0*v( 2) +  1d0*v( 3)) / 12d0 !&
                                        poly(0) = ( CpR(0,1)*v(-3) + CpR(0,2)*v(-2) + CpR(0,3)*v(-1) + CpR(0,4)*v( 0)) !&
                                        poly(1) = ( CpR(1,1)*v(-2) + CpR(1,2)*v(-1) + CpR(1,3)*v( 0) + CpR(1,4)*v( 1)) !&
                                        poly(2) = ( CpR(2,1)*v(-1) + CpR(2,2)*v( 0) + CpR(2,3)*v( 1) + CpR(2,4)*v( 2)) !&
                                        poly(3) = ( CpR(3,1)*v( 0) + CpR(3,2)*v( 1) + CpR(3,3)*v( 2) + CpR(3,4)*v( 3)) !&
                                    else
                                        poly(0) = (-1d0*v(-1) +  5d0*v( 0) +  2d0*v( 1)) / 6d0 !&
                                        poly(1) = ( 2d0*v( 0) +  5d0*v( 1) -  1d0*v( 2)) / 6d0 !&
                                        poly(2) = ( 2d0*v(-2) -  7d0*v(-1) + 11d0*v( 0)) / 6d0 !&
                                        poly(3) = ( 3d0*v( 0) + 13d0*v( 1) -  5d0*v( 2) +  1d0*v( 3)) / 12d0 !&
                                        poly(4) = (-3d0*v(-3) + 13d0*v(-2) - 23d0*v(-1) + 25d0*v( 0)) / 12d0 !&
                                    end if

                                    if (wenojs) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)

                                    elseif (mapped_weno) then
                                        alpha = d_cbR_${XYZ}$ (:, j)/(beta*beta)
                                        omega = alpha/sum(alpha)
                                        alpha = (d_cbR_${XYZ}$ (:, j)*(1d0 + d_cbR_${XYZ}$ (:, j) - 3d0*omega) + omega**2d0) &
                                                *(omega/(d_cbR_${XYZ}$ (:, j)**2d0 + omega*(1d0 - 2d0*d_cbR_${XYZ}$ (:, j))))

                                    elseif (wenoz) then
                                        alpha = d_cbR_${XYZ}$ (:, j)*(1d0 + (tau/beta)**wenoz_q)

                                    elseif (teno) then
                                        alpha = delta*d_cbR_${XYZ}$ (:, j)

                                    end if

                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = sum(omega*poly)

                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                end if
            #:endfor
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
                                 norm_dir, weno_dir)

        type(scalar_field), dimension(:), intent(IN) :: v_vf

        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: weno_dir

        integer :: j, k, l, q

        ! Determining the number of cell-average variables which will be
        ! WENO-reconstructed and mapping their indical bounds in the x-,
        ! y- and z-directions to those in the s1-, s2- and s3-directions
        ! as to reshape the inputted data in the coordinate direction of
        ! the WENO reconstruction
        v_size = ubound(v_vf, 1)
        !$acc update device(v_size)

        if (weno_dir == 1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, v_size
                do q = is3_weno%beg, is3_weno%end
                    do l = is2_weno%beg, is2_weno%end
                        do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                            v_rs_ws_x(k, l, q, j) = v_vf(j)%sf(k, l, q)
                        end do
                    end do
                end do
            end do
            !$acc end parallel loop
        end if
        ! ==================================================================

        ! Reshaping/Projecting onto Characteristic Fields in y-direction ===
        if (n == 0) return

        if (weno_dir == 2) then
#if MFC_cuTENSOR
            if (cu_tensor) then
                if (p == 0) then
                    block
                        use CuTensorEx

                        !$acc host_data use_device(v_rs_ws_x, v_rs_ws_y)
                        v_rs_ws_y = reshape(v_rs_ws_x, shape=[n + 1 + 2*buff_size, m + 2*buff_size + 1, p + 1, sys_size], order=[2, 1, 3, 4])
                        !$acc end host_data
                    end block
                else
                    block
                        use CuTensorEx

                        !$acc host_data use_device(v_rs_ws_x, v_rs_ws_y)
                        v_rs_ws_y = reshape(v_rs_ws_x, shape=[n + 1 + 2*buff_size, m + 2*buff_size + 1, p + 1 + 2*buff_size, sys_size], order=[2, 1, 3, 4])
                        !$acc end host_data
                    end block
                end if
            else
#endif
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, v_size
                    do q = is3_weno%beg, is3_weno%end
                        do l = is2_weno%beg, is2_weno%end
                            do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                                v_rs_ws_y(k, l, q, j) = v_vf(j)%sf(l, k, q)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
#if MFC_cuTENSOR
            end if
#endif
        end if

        ! ==================================================================

        ! Reshaping/Projecting onto Characteristic Fields in z-direction ===
        if (p == 0) return
        if (weno_dir == 3) then
#if MFC_cuTENSOR
            if (cu_tensor) then
                block
                    use CuTensorEx

                    !$acc host_data use_device(v_rs_ws_x, v_rs_ws_z)
                    v_rs_ws_z = reshape(v_rs_ws_x, shape=[p + 1 + 2*buff_size, n + 2*buff_size + 1, m + 2*buff_size + 1, sys_size], order=[3, 2, 1, 4])
                    !$acc end host_data
                end block
            else
#endif
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, v_size
                    do q = is3_weno%beg, is3_weno%end
                        do l = is2_weno%beg, is2_weno%end
                            do k = is1_weno%beg - weno_polyn, is1_weno%end + weno_polyn
                                v_rs_ws_z(k, l, q, j) = v_vf(j)%sf(q, l, k)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
#if MFC_cuTENSOR
            end if
#endif
        end if

        ! ==================================================================

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
        !!  @param k Second-coordinate cell index
        !!  @param l Third-coordinate cell index
    subroutine s_preserve_monotonicity(v_rs_ws, vL_rs_vf, vR_rs_vf)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(IN) :: v_rs_ws
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: vL_rs_vf, vR_rs_vf

        integer :: i, j, k, l

        real(kind(0d0)), dimension(-1:1) :: d !< Curvature measures at the zone centers

        real(kind(0d0)) :: d_MD, d_LC !<
            !! Median (md) curvature and large curvature (LC) measures

        ! The left and right upper bounds (UL), medians, large curvatures,
        ! minima, and maxima of the WENO-reconstructed values of the cell-
        ! average variables.
        real(kind(0d0)) :: vL_UL, vR_UL
        real(kind(0d0)) :: vL_MD, vR_MD
        real(kind(0d0)) :: vL_LC, vR_LC
        real(kind(0d0)) :: vL_min, vR_min
        real(kind(0d0)) :: vL_max, vR_max

        real(kind(0d0)), parameter :: alpha = 2d0 !>
            !! Determines the maximum Courant–Friedrichs–Lewy (CFL) number that
            !! may be utilized with the scheme. In theory, for stability, a CFL
            !! number less than 1/(1+alpha) is necessary. The default value for
            !! alpha is 2.

        real(kind(0d0)), parameter :: beta = 4d0/3d0 !<
            !! Determines the amount of freedom available from utilizing a large
            !! value for the local curvature. The default value for beta is 4/3.

        real(kind(0d0)), parameter :: alpha_mp = 2d0
        real(kind(0d0)), parameter :: beta_mp = 4d0/3d0

        !$acc parallel loop gang vector collapse (4)  default(present) private(d)
        do l = is3_weno%beg, is3_weno%end
            do k = is2_weno%beg, is2_weno%end
                do j = is1_weno%beg, is1_weno%end
                    do i = 1, v_size
                        d(-1) = v_rs_ws(j, k, l, i) &
                                + v_rs_ws(j - 2, k, l, i) &
                                - v_rs_ws(j - 1, k, l, i) &
                                *2d0
                        d(0) = v_rs_ws(j + 1, k, l, i) &
                               + v_rs_ws(j - 1, k, l, i) &
                               - v_rs_ws(j, k, l, i) &
                               *2d0
                        d(1) = v_rs_ws(j + 2, k, l, i) &
                               + v_rs_ws(j, k, l, i) &
                               - v_rs_ws(j + 1, k, l, i) &
                               *2d0

                        d_MD = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
                               *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                                    *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
                               *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                                    abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

                        d_LC = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
                               *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                                    *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
                               *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                                    abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

                        vL_UL = v_rs_ws(j, k, l, i) &
                                - (v_rs_ws(j + 1, k, l, i) &
                                   - v_rs_ws(j, k, l, i))*alpha_mp

                        vL_MD = (v_rs_ws(j, k, l, i) &
                                 + v_rs_ws(j - 1, k, l, i) &
                                 - d_MD)*5d-1

                        vL_LC = v_rs_ws(j, k, l, i) &
                                - (v_rs_ws(j + 1, k, l, i) &
                                   - v_rs_ws(j, k, l, i))*5d-1 + beta_mp*d_LC

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
                                               + (sign(5d-1, vL_min - vL_rs_vf(j, k, l, i)) &
                                                  + sign(5d-1, vL_max - vL_rs_vf(j, k, l, i))) &
                                               *min(abs(vL_min - vL_rs_vf(j, k, l, i)), &
                                                    abs(vL_max - vL_rs_vf(j, k, l, i)))
                        ! END: Left Monotonicity Preserving Bound ==========================

                        ! Right Monotonicity Preserving Bound ==============================
                        d(-1) = v_rs_ws(j, k, l, i) &
                                + v_rs_ws(j - 2, k, l, i) &
                                - v_rs_ws(j - 1, k, l, i) &
                                *2d0
                        d(0) = v_rs_ws(j + 1, k, l, i) &
                               + v_rs_ws(j - 1, k, l, i) &
                               - v_rs_ws(j, k, l, i) &
                               *2d0
                        d(1) = v_rs_ws(j + 2, k, l, i) &
                               + v_rs_ws(j, k, l, i) &
                               - v_rs_ws(j + 1, k, l, i) &
                               *2d0

                        d_MD = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
                               *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                                    *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
                               *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                                    abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

                        d_LC = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
                               *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                                    *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
                               *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                                    abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

                        vR_UL = v_rs_ws(j, k, l, i) &
                                + (v_rs_ws(j, k, l, i) &
                                   - v_rs_ws(j - 1, k, l, i))*alpha_mp

                        vR_MD = (v_rs_ws(j, k, l, i) &
                                 + v_rs_ws(j + 1, k, l, i) &
                                 - d_MD)*5d-1

                        vR_LC = v_rs_ws(j, k, l, i) &
                                + (v_rs_ws(j, k, l, i) &
                                   - v_rs_ws(j - 1, k, l, i))*5d-1 + beta_mp*d_LC

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
                                               + (sign(5d-1, vR_min - vR_rs_vf(j, k, l, i)) &
                                                  + sign(5d-1, vR_max - vR_rs_vf(j, k, l, i))) &
                                               *min(abs(vR_min - vR_rs_vf(j, k, l, i)), &
                                                    abs(vR_max - vR_rs_vf(j, k, l, i)))
                        ! END: Right Monotonicity Preserving Bound =========================
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine s_preserve_monotonicity

    !>  Module deallocation and/or disassociation procedures
    subroutine s_finalize_weno_module()

        if (weno_order == 1) return

        ! Deallocating the WENO-stencil of the WENO-reconstructed variables

        !deallocate(vL_rs_vf_x, vR_rs_vf_x)
        @:DEALLOCATE(v_rs_ws_x)

        ! Deallocating WENO coefficients in x-direction ====================
        @:DEALLOCATE(poly_coef_cbL_x, poly_coef_cbR_x)
        @:DEALLOCATE(d_cbL_x, d_cbR_x)
        @:DEALLOCATE(beta_coef_x)
        ! ==================================================================

        ! Deallocating WENO coefficients in y-direction ====================
        if (n == 0) return

        !deallocate(vL_rs_vf_y, vR_rs_vf_y)
        @:DEALLOCATE(v_rs_ws_y)

        @:DEALLOCATE(poly_coef_cbL_y, poly_coef_cbR_y)
        @:DEALLOCATE(d_cbL_y, d_cbR_y)
        @:DEALLOCATE(beta_coef_y)
        ! ==================================================================

        ! Deallocating WENO coefficients in z-direction ====================
        if (p == 0) return

        !deallocate(vL_rs_vf_z, vR_rs_vf_z)
        @:DEALLOCATE(v_rs_ws_z)

        @:DEALLOCATE(poly_coef_cbL_z, poly_coef_cbR_z)
        @:DEALLOCATE(d_cbL_z, d_cbR_z)
        @:DEALLOCATE(beta_coef_z)
        ! ==================================================================

    end subroutine s_finalize_weno_module

end module m_weno
