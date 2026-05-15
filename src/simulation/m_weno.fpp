!>
!! @file
!! @brief Contains module m_weno
#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief WENO/WENO-Z/TENO reconstruction with optional monotonicity-preserving bounds and mapped weights
module m_weno

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    ! $:USE_GPU_MODULE()

    use m_mpi_proxy
    use m_thinc, only: s_thinc_compression
    use m_nvtx

    private; public :: s_initialize_weno_module, s_finalize_weno_module, s_weno, s_pack_weno_input_arr

    !> @name The cell-average variables that will be WENO-reconstructed unpacked into an array for performance
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: v_rs_weno
    !> @}
    $:GPU_DECLARE(create='[v_rs_weno]')

    ! WENO Coefficients

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at the left and right quadrature points (QP), in
    !! the x-, y- and z-directions. Note that the first dimension of the array identifies the polynomial, the second dimension
    !! identifies the position of its coefficients and the last dimension denotes the cell-location in the relevant coordinate
    !! direction.
    !> @{
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbL_x
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbL_y
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbL_z
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbR_x
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbR_y
    real(wp), target, allocatable, dimension(:,:,:) :: poly_coef_cbR_z
    !> @}
    $:GPU_DECLARE(create='[poly_coef_cbL_x, poly_coef_cbL_y, poly_coef_cbL_z]')
    $:GPU_DECLARE(create='[poly_coef_cbR_x, poly_coef_cbR_y, poly_coef_cbR_z]')

    !> @name The ideal weights at the left and the right cell-boundaries and at the left and the right quadrature points, in x-, y-
    !! and z-directions. Note that the first dimension of the array identifies the weight, while the last denotes the cell-location
    !! in the relevant coordinate direction.
    !> @{
    real(wp), target, allocatable, dimension(:,:) :: d_cbL_x
    real(wp), target, allocatable, dimension(:,:) :: d_cbL_y
    real(wp), target, allocatable, dimension(:,:) :: d_cbL_z
    real(wp), target, allocatable, dimension(:,:) :: d_cbR_x
    real(wp), target, allocatable, dimension(:,:) :: d_cbR_y
    real(wp), target, allocatable, dimension(:,:) :: d_cbR_z
    !> @}
    $:GPU_DECLARE(create='[d_cbL_x, d_cbL_y, d_cbL_z, d_cbR_x, d_cbR_y, d_cbR_z]')

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note that the first array dimension identifies the
    !! smoothness indicator, the second identifies the position of its coefficients and the last denotes the cell-location in the
    !! relevant coordinate direction.
    !> @{
    real(wp), target, allocatable, dimension(:,:,:) :: beta_coef_x
    real(wp), target, allocatable, dimension(:,:,:) :: beta_coef_y
    real(wp), target, allocatable, dimension(:,:,:) :: beta_coef_z
    !> @}
    $:GPU_DECLARE(create='[beta_coef_x, beta_coef_y, beta_coef_z]')

    ! END: WENO Coefficients

    integer :: v_size  !< Number of WENO-reconstructed cell-average variables
    $:GPU_DECLARE(create='[v_size]')

    logical :: uniform_grid(3)  !< True if grid spacing is uniform in each direction
    $:GPU_DECLARE(create='[uniform_grid]')

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1_weno, is2_weno, is3_weno
#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
    $:GPU_DECLARE(create='[is1_weno, is2_weno, is3_weno]')
#endif
    !
    !> @}

contains

    !> Initialize the WENO module
    impure subroutine s_initialize_weno_module

        if (weno_order == 1) return

        ! Allocating/Computing WENO Coefficients in x-direction
        is1_weno%beg = -buff_size; is1_weno%end = m - is1_weno%beg
        if (n == 0) then
            is2_weno%beg = 0
        else
            is2_weno%beg = -buff_size
        end if

        is2_weno%end = n - is2_weno%beg

        if (p == 0) then
            is3_weno%beg = 0
        else
            is3_weno%beg = -buff_size
        end if

        is3_weno%end = p - is3_weno%beg

        @:ALLOCATE(poly_coef_cbL_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_x(0:weno_num_stencils, is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_x(0:weno_num_stencils, is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_x(is1_weno%beg + weno_polyn:is1_weno%end - weno_polyn, 0:weno_polyn, &
                   & 0:weno_polyn*(weno_polyn + 1)/2 - 1))
        ! Number of cross terms for dvd = (k-1)(k-1+1)/2, where weno_polyn = k-1 Note: k-1 not k because we are using value
        ! differences (dvd) not the values themselves

        call s_compute_weno_coefficients(1, is1_weno)

        @:ALLOCATE(v_rs_weno(is1_weno%beg:is1_weno%end, is2_weno%beg:is2_weno%end, is3_weno%beg:is3_weno%end, 1:sys_size))

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

        @:ALLOCATE(poly_coef_cbL_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_y(0:weno_num_stencils, is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_y(0:weno_num_stencils, is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_y(is2_weno%beg + weno_polyn:is2_weno%end - weno_polyn, 0:weno_polyn, &
                   & 0:weno_polyn*(weno_polyn + 1)/2 - 1))

        call s_compute_weno_coefficients(2, is2_weno)

        ! Allocating/Computing WENO Coefficients in z-direction
        if (p == 0) return

        is2_weno%beg = -buff_size; is2_weno%end = n - is2_weno%beg
        is1_weno%beg = -buff_size; is1_weno%end = m - is1_weno%beg
        is3_weno%beg = -buff_size; is3_weno%end = p - is3_weno%beg

        @:ALLOCATE(poly_coef_cbL_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))
        @:ALLOCATE(poly_coef_cbR_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, 0:weno_polyn - 1))

        @:ALLOCATE(d_cbL_z(0:weno_num_stencils, is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn))
        @:ALLOCATE(d_cbR_z(0:weno_num_stencils, is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn))

        @:ALLOCATE(beta_coef_z(is3_weno%beg + weno_polyn:is3_weno%end - weno_polyn, 0:weno_polyn, &
                   & 0:weno_polyn*(weno_polyn + 1)/2 - 1))

        call s_compute_weno_coefficients(3, is3_weno)

    end subroutine s_initialize_weno_module

    !> Compute WENO polynomial coefficients, ideal weights, and smoothness indicators for a given direction
    subroutine s_compute_weno_coefficients(weno_dir, is)

        ! Compute WENO coefficients for a given coordinate direction. Shu (1997)
        integer, intent(in)               :: weno_dir
        type(int_bounds_info), intent(in) :: is
        integer                           :: s
        real(wp), pointer, dimension(:)   :: s_cb => null()  !< Cell-boundary locations in the s-direction
        type(bc_dir_t)                    :: bc_s            !< Boundary conditions (BC) in the s-direction
        integer                           :: i               !< Generic loop iterator
        real(wp)                          :: w(1:8)          !< Intermediate var for ideal weights: s_cb across overall stencil
        real(wp)                          :: ys(1:4)         !< Intermediate var for poly & beta: diff(s_cb) across sub-stencil
        real(wp)                          :: h0              !< Reference spacing for uniform-grid detection

        ! Determine cell count, boundary locations, and BCs for selected WENO direction

        if (weno_dir == 1) then
            s = m; s_cb => x%cb; bc_s = bc%x
        else if (weno_dir == 2) then
            s = n; s_cb => y%cb; bc_s = bc%y
        else
            s = p; s_cb => z%cb; bc_s = bc%z
        end if

        #:for WENO_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            ! Computing WENO3 Coefficients
            if (weno_dir == ${WENO_DIR}$) then
                if (weno_order == 3) then
                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn
                        ! Polynomial reconstruction coefficients
                        poly_coef_cbR_${XYZ}$ (i + 1, 0, 0) = (s_cb(i) - s_cb(i + 1))/(s_cb(i) - s_cb(i + 2))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, 0) = (s_cb(i) - s_cb(i + 1))/(s_cb(i - 1) - s_cb(i + 1))

                        poly_coef_cbL_${XYZ}$ (i + 1, 0, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 0, 0)
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, 0) = -poly_coef_cbR_${XYZ}$ (i + 1, 1, 0)

                        ! Ideal (linear) weights
                        d_cbR_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i + 1))/(s_cb(i - 1) - s_cb(i + 2))
                        d_cbL_${XYZ}$ (0, i + 1) = (s_cb(i - 1) - s_cb(i))/(s_cb(i - 1) - s_cb(i + 2))

                        d_cbR_${XYZ}$ (1, i + 1) = 1._wp - d_cbR_${XYZ}$ (0, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1._wp - d_cbL_${XYZ}$ (0, i + 1)

                        ! Smoothness indicator coefficients
                        beta_coef_${XYZ}$ (i + 1, 0, 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp/(s_cb(i) - s_cb(i + 2))**2._wp
                        beta_coef_${XYZ}$ (i + 1, 1, 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp/(s_cb(i - 1) - s_cb(i + 1))**2._wp
                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood of beginning and end Riemann state extrapolation
                    ! BC to avoid any contributions from outside of the physical domain during the WENO reconstruction
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
                else if (weno_order == 5) then
                    do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn
                        ! Polynomial reconstruction coefficients
                        poly_coef_cbR_${XYZ}$ (i + 1, 0, &
                                               & 0) = ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/((s_cb(i) - s_cb(i &
                                               & + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, &
                                               & 0) = ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) &
                                               & - s_cb(i + 2))*(s_cb(i + 2) - s_cb(i)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 1, &
                                               & 1) = ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/((s_cb(i - 1) &
                                               & - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, &
                                               & 1) = ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/((s_cb(i - 2) &
                                               & - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, &
                                               & 0) = ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/((s_cb(i) - s_cb(i + 3)) &
                                               & *(s_cb(i + 3) - s_cb(i + 1)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, &
                                               & 0) = ((s_cb(i) - s_cb(i - 1))*(s_cb(i) - s_cb(i + 1)))/((s_cb(i - 1) - s_cb(i &
                                               & + 2))*(s_cb(i) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 1, &
                                               & 1) = ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/((s_cb(i - 1) - s_cb(i &
                                               & + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, &
                                               & 1) = ((s_cb(i - 1) - s_cb(i))*(s_cb(i) - s_cb(i + 1)))/((s_cb(i - 2) - s_cb(i)) &
                                               & *(s_cb(i - 2) - s_cb(i + 1)))

                        poly_coef_cbR_${XYZ}$ (i + 1, 0, &
                                               & 1) = ((s_cb(i) - s_cb(i + 2)) + (s_cb(i + 1) - s_cb(i + 3)))/((s_cb(i) - s_cb(i &
                                               & + 2))*(s_cb(i) - s_cb(i + 3)))*((s_cb(i) - s_cb(i + 1)))
                        poly_coef_cbR_${XYZ}$ (i + 1, 2, &
                                               & 0) = ((s_cb(i - 2) - s_cb(i + 1)) + (s_cb(i - 1) - s_cb(i + 1)))/((s_cb(i - 1) &
                                               & - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 2)))*((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 0, &
                                               & 1) = ((s_cb(i) - s_cb(i + 2)) + (s_cb(i) - s_cb(i + 3)))/((s_cb(i) - s_cb(i + 2)) &
                                               & *(s_cb(i) - s_cb(i + 3)))*((s_cb(i + 1) - s_cb(i)))
                        poly_coef_cbL_${XYZ}$ (i + 1, 2, &
                                               & 0) = ((s_cb(i - 2) - s_cb(i)) + (s_cb(i - 1) - s_cb(i + 1)))/((s_cb(i - 2) &
                                               & - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))*((s_cb(i) - s_cb(i + 1)))

                        ! Ideal (linear) weights
                        d_cbR_${XYZ}$ (0, &
                                       & i + 1) = ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/((s_cb(i - 2) &
                                       & - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbR_${XYZ}$ (2, &
                                       & i + 1) = ((s_cb(i + 1) - s_cb(i + 2))*(s_cb(i + 1) - s_cb(i + 3)))/((s_cb(i - 2) &
                                       & - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))
                        d_cbL_${XYZ}$ (0, &
                                       & i + 1) = ((s_cb(i - 2) - s_cb(i))*(s_cb(i) - s_cb(i - 1)))/((s_cb(i - 2) - s_cb(i + 3)) &
                                       & *(s_cb(i + 3) - s_cb(i - 1)))
                        d_cbL_${XYZ}$ (2, &
                                       & i + 1) = ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))/((s_cb(i - 2) - s_cb(i + 2)) &
                                       & *(s_cb(i - 2) - s_cb(i + 3)))

                        d_cbR_${XYZ}$ (1, i + 1) = 1._wp - d_cbR_${XYZ}$ (0, i + 1) - d_cbR_${XYZ}$ (2, i + 1)
                        d_cbL_${XYZ}$ (1, i + 1) = 1._wp - d_cbL_${XYZ}$ (0, i + 1) - d_cbL_${XYZ}$ (2, i + 1)

                        ! Smoothness indicator coefficients
                        beta_coef_${XYZ}$ (i + 1, 0, &
                                           & 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1)) &
                                           & **2._wp)/((s_cb(i) - s_cb(i + 3))**2._wp*(s_cb(i + 1) - s_cb(i + 3))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 0, &
                                           & 1) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(19._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - s_cb(i + 1)) + 2._wp*(s_cb(i + 2) - s_cb(i)) &
                                           & *((s_cb(i + 2) - s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - s_cb(i + 2)) &
                                           & *(s_cb(i) - s_cb(i + 3))**2._wp*(s_cb(i + 3) - s_cb(i + 1)))

                        beta_coef_${XYZ}$ (i + 1, 0, &
                                           & 2) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) &
                                           & + ((s_cb(i + 2) - s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2._wp)/((s_cb(i) - s_cb(i &
                                           & + 2))**2._wp*(s_cb(i) - s_cb(i + 3))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 1, &
                                           & 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + (s_cb(i) - s_cb(i - 1))**2._wp + (s_cb(i) - s_cb(i - 1))*(s_cb(i + 1) - s_cb(i))) &
                                           & /((s_cb(i - 1) - s_cb(i + 2))**2._wp*(s_cb(i) - s_cb(i + 2))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 1, &
                                           & 1) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*((s_cb(i) - s_cb(i + 1))*((s_cb(i) &
                                           & - s_cb(i - 1)) + 20._wp*(s_cb(i + 1) - s_cb(i))) + (2._wp*(s_cb(i) - s_cb(i - 1)) &
                                           & + (s_cb(i + 1) - s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - s_cb(i - 1)) &
                                           & *(s_cb(i - 1) - s_cb(i + 2))**2._wp*(s_cb(i + 2) - s_cb(i)))

                        beta_coef_${XYZ}$ (i + 1, 1, &
                                           & 2) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1)) &
                                           & **2._wp)/((s_cb(i - 1) - s_cb(i + 1))**2._wp*(s_cb(i - 1) - s_cb(i + 2))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 2, &
                                           & 0) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(12._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - s_cb(i - 1)))**2._wp + 3._wp*((s_cb(i) &
                                           & - s_cb(i - 2)) + (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) &
                                           & - s_cb(i + 1))**2._wp*(s_cb(i - 1) - s_cb(i + 1))**2._wp)

                        beta_coef_${XYZ}$ (i + 1, 2, &
                                           & 1) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(19._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - s_cb(i + 1))) + 2._wp*(s_cb(i + 1) - s_cb(i &
                                           & - 1))*((s_cb(i) - s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) &
                                           & - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2._wp*(s_cb(i + 1) - s_cb(i - 1)))

                        beta_coef_${XYZ}$ (i + 1, 2, &
                                           & 2) = 4._wp*(s_cb(i) - s_cb(i + 1))**2._wp*(10._wp*(s_cb(i + 1) - s_cb(i))**2._wp &
                                           & + (s_cb(i) - s_cb(i - 1))**2._wp + (s_cb(i) - s_cb(i - 1))*(s_cb(i + 1) - s_cb(i))) &
                                           & /((s_cb(i - 2) - s_cb(i))**2._wp*(s_cb(i - 2) - s_cb(i + 1))**2._wp)
                    end do

                    ! Modifying the ideal weights coefficients in the neighborhood of beginning and end Riemann state extrapolation
                    ! BC to avoid any contributions from outside of the physical domain during the WENO reconstruction
                    if (null_weights) then
                        if (bc_s%beg == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (1:2,0) = 0._wp; d_cbR_${XYZ}$ (0, 0) = 1._wp
                            d_cbL_${XYZ}$ (1:2,0) = 0._wp; d_cbL_${XYZ}$ (0, 0) = 1._wp
                            d_cbR_${XYZ}$ (2, 1) = 0._wp; d_cbR_${XYZ}$ (:,1) = d_cbR_${XYZ}$ (:,1)/sum(d_cbR_${XYZ}$ (:,1))
                            d_cbL_${XYZ}$ (2, 1) = 0._wp; d_cbL_${XYZ}$ (:,1) = d_cbL_${XYZ}$ (:,1)/sum(d_cbL_${XYZ}$ (:,1))
                        end if

                        if (bc_s%end == BC_RIEMANN_EXTRAP) then
                            d_cbR_${XYZ}$ (0, s - 1) = 0._wp; d_cbR_${XYZ}$ (:,s - 1) = d_cbR_${XYZ}$ (:, &
                                           & s - 1)/sum(d_cbR_${XYZ}$ (:,s - 1))
                            d_cbL_${XYZ}$ (0, s - 1) = 0._wp; d_cbL_${XYZ}$ (:,s - 1) = d_cbL_${XYZ}$ (:, &
                                           & s - 1)/sum(d_cbL_${XYZ}$ (:,s - 1))
                            d_cbR_${XYZ}$ (0:1,s) = 0._wp; d_cbR_${XYZ}$ (2, s) = 1._wp
                            d_cbL_${XYZ}$ (0:1,s) = 0._wp; d_cbL_${XYZ}$ (2, s) = 1._wp
                        end if
                    end if
                else
                    if (.not. teno) then
                        do i = is%beg - 1 + weno_polyn, is%end - 1 - weno_polyn
                            ! Reference: Shu (1997) "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes
                            ! for Hyperbolic Conservation Laws" Equation 2.20: Polynomial Coefficients (poly_coef_cb) Equation 2.61:
                            ! Smoothness Indicators (beta_coef) To reduce computational cost, we leverage the fact that all
                            ! polynomial coefficients in a stencil sum to 1 and compute the polynomial coefficients (poly_coef_cb)
                            ! for the cell value differences (dvd) instead of the values themselves. The computation of coefficients
                            ! is further simplified by using grid spacing (ys or w) rather than the grid locations (s_cb) directly.
                            ! Ideal weights (d_cb) are obtained by comparing the grid location coefficients of the polynomial
                            ! coefficients. The smoothness indicators (beta_coef) are calculated through numerical differentiation
                            ! and integration of each cross term of the polynomial coefficients, using the cell value differences
                            ! (dvd) instead of the values themselves. While the polynomial coefficients sum to 1, the derivative of
                            ! 1 is 0, which means it does not create additional cross terms in the smoothness indicators.

                            w = s_cb(i - 3:i + 4) - s_cb(i)  ! Offset using s_cb(i) to reduce floating point error
                            d_cbR_${XYZ}$ (0, &
                                           & i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7)) &
                                           & *(w(1) - w(8)))
                            d_cbR_${XYZ}$ (1, &
                                           & i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1) &
                                           & *w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7) &
                                           & *w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7)) &
                                           & *(w(2) - w(8)))
                            d_cbR_${XYZ}$ (2, &
                                           & i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2) &
                                           & *w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7) &
                                           & *w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8)) &
                                           & *(w(3) - w(8)))
                            d_cbR_${XYZ}$ (3, &
                                           & i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8)) &
                                           & *(w(3) - w(8)))

                            w = s_cb(i + 4:i - 3:-1) - s_cb(i)
                            d_cbL_${XYZ}$ (0, &
                                           & i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(3) - w(5)))/((w(1) - w(8))*(w(2) - w(8)) &
                                           & *(w(3) - w(8)))
                            d_cbL_${XYZ}$ (1, &
                                           & i + 1) = ((w(1) - w(5))*(w(2) - w(5))*(w(5) - w(8))*(w(1)*w(2) + w(1)*w(3) + w(2) &
                                           & *w(3) - w(1)*w(7) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) - w(3)*w(7) - w(3)*w(8) + w(7) &
                                           & *w(8) + w(7)**2 + w(8)**2))/((w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7))*(w(2) - w(8)) &
                                           & *(w(3) - w(8)))
                            d_cbL_${XYZ}$ (2, &
                                           & i + 1) = ((w(1) - w(5))*(w(5) - w(7))*(w(5) - w(8))*(w(1)*w(2) - w(1)*w(6) - w(1) &
                                           & *w(7) - w(2)*w(6) - w(1)*w(8) - w(2)*w(7) - w(2)*w(8) + w(6)*w(7) + w(6)*w(8) + w(7) &
                                           & *w(8) + w(1)**2 + w(2)**2))/((w(1) - w(6))*(w(1) - w(7))*(w(1) - w(8))*(w(2) - w(7)) &
                                           & *(w(2) - w(8)))
                            d_cbL_${XYZ}$ (3, &
                                           & i + 1) = ((w(5) - w(6))*(w(5) - w(7))*(w(5) - w(8)))/((w(1) - w(6))*(w(1) - w(7)) &
                                           & *(w(1) - w(8)))
                            ! Note: Left has the reversed order of both points and coefficients compared to the right

                            ys = s_cb(i + 1:i + 4) - s_cb(i:i + 3)
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, &
                                                   & 0) = (ys(1)*ys(2)*(ys(2) + ys(3)))/((ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, &
                                                   & 1) = -(ys(1)*ys(2)*(3*ys(2)**2 + 6*ys(2)*ys(3) + 3*ys(2)*ys(4) + 2*ys(1) &
                                                   & *ys(2) + 3*ys(3)**2 + 3*ys(3)*ys(4) + 2*ys(1)*ys(3) + ys(4)**2 + ys(1)*ys(4)) &
                                                   & )/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 0, &
                                                   & 2) = (ys(1)*(ys(1)**2 + 3*ys(1)*ys(2) + 2*ys(1)*ys(3) + ys(4)*ys(1) + 3*ys(2) &
                                                   & **2 + 4*ys(2)*ys(3) + 2*ys(4)*ys(2) + ys(3)**2 + ys(4)*ys(3)))/((ys(1) &
                                                   & + ys(2))*(ys(1) + ys(2) + ys(3))*(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i:i + 3) - s_cb(i - 1:i + 2)
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, &
                                                   & 0) = -(ys(2)*ys(3)*(ys(1) + ys(2)))/((ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, &
                                                   & 1) = (ys(2)*(ys(1) + ys(2))*(ys(2)**2 + 4*ys(2)*ys(3) + 2*ys(2)*ys(4) + ys(1) &
                                                   & *ys(2) + 3*ys(3)**2 + 3*ys(3)*ys(4) + 2*ys(1)*ys(3) + ys(4)**2 + ys(1)*ys(4)) &
                                                   & )/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 1, &
                                                   & 2) = (ys(2)*ys(3)*(ys(3) + ys(4)))/((ys(1) + ys(2))*(ys(1) + ys(2) + ys(3)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i - 1:i + 2) - s_cb(i - 2:i + 1)
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, &
                                                   & 0) = (ys(3)*(ys(2) + ys(3))*(ys(1) + ys(2) + ys(3)))/((ys(3) + ys(4))*(ys(2) &
                                                   & + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, &
                                                   & 1) = (ys(3)*ys(4)*(ys(1)**2 + 3*ys(1)*ys(2) + 3*ys(1)*ys(3) + ys(4)*ys(1) &
                                                   & + 3*ys(2)**2 + 6*ys(2)*ys(3) + 2*ys(4)*ys(2) + 3*ys(3)**2 + 2*ys(4)*ys(3))) &
                                                   & /((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 2, &
                                                   & 2) = -(ys(3)*ys(4)*(ys(2) + ys(3)))/((ys(1) + ys(2))*(ys(1) + ys(2) + ys(3)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i - 2:i + 1) - s_cb(i - 3:i)
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, &
                                                   & 0) = (ys(4)*(ys(2)**2 + 4*ys(2)*ys(3) + 4*ys(2)*ys(4) + ys(1)*ys(2) + 3*ys(3) &
                                                   & **2 + 6*ys(3)*ys(4) + 2*ys(1)*ys(3) + 3*ys(4)**2 + 2*ys(1)*ys(4)))/((ys(3) &
                                                   & + ys(4))*(ys(2) + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, &
                                                   & 1) = -(ys(4)*(ys(3) + ys(4))*(ys(1)**2 + 3*ys(1)*ys(2) + 3*ys(1)*ys(3) &
                                                   & + 2*ys(1)*ys(4) + 3*ys(2)**2 + 6*ys(2)*ys(3) + 4*ys(2)*ys(4) + 3*ys(3)**2 &
                                                   & + 4*ys(3)*ys(4) + ys(4)**2))/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) &
                                                   & + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbR_${XYZ}$ (i + 1, 3, &
                                                   & 2) = (ys(4)*(ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)))/((ys(1) + ys(2))*(ys(1) &
                                                   & + ys(2) + ys(3))*(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i + 1:i - 2:-1) - s_cb(i:i - 3:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, &
                                                   & 2) = (ys(1)*ys(2)*(ys(2) + ys(3)))/((ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, &
                                                   & 1) = -(ys(1)*ys(2)*(3*ys(2)**2 + 6*ys(2)*ys(3) + 3*ys(2)*ys(4) + 2*ys(1) &
                                                   & *ys(2) + 3*ys(3)**2 + 3*ys(3)*ys(4) + 2*ys(1)*ys(3) + ys(4)**2 + ys(1)*ys(4)) &
                                                   & )/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 3, &
                                                   & 0) = (ys(1)*(ys(1)**2 + 3*ys(1)*ys(2) + 2*ys(1)*ys(3) + ys(4)*ys(1) + 3*ys(2) &
                                                   & **2 + 4*ys(2)*ys(3) + 2*ys(4)*ys(2) + ys(3)**2 + ys(4)*ys(3)))/((ys(1) &
                                                   & + ys(2))*(ys(1) + ys(2) + ys(3))*(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i + 2:i - 1:-1) - s_cb(i + 1:i - 2:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, &
                                                   & 2) = -(ys(2)*ys(3)*(ys(1) + ys(2)))/((ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, &
                                                   & 1) = (ys(2)*(ys(1) + ys(2))*(ys(2)**2 + 4*ys(2)*ys(3) + 2*ys(2)*ys(4) + ys(1) &
                                                   & *ys(2) + 3*ys(3)**2 + 3*ys(3)*ys(4) + 2*ys(1)*ys(3) + ys(4)**2 + ys(1)*ys(4)) &
                                                   & )/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 2, &
                                                   & 0) = (ys(2)*ys(3)*(ys(3) + ys(4)))/((ys(1) + ys(2))*(ys(1) + ys(2) + ys(3)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i + 3:i:-1) - s_cb(i + 2:i - 1:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, &
                                                   & 2) = (ys(3)*(ys(2) + ys(3))*(ys(1) + ys(2) + ys(3)))/((ys(3) + ys(4))*(ys(2) &
                                                   & + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, &
                                                   & 1) = (ys(3)*ys(4)*(ys(1)**2 + 3*ys(1)*ys(2) + 3*ys(1)*ys(3) + ys(4)*ys(1) &
                                                   & + 3*ys(2)**2 + 6*ys(2)*ys(3) + 2*ys(4)*ys(2) + 3*ys(3)**2 + 2*ys(4)*ys(3))) &
                                                   & /((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))*(ys(1) &
                                                   & + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 1, &
                                                   & 0) = -(ys(3)*ys(4)*(ys(2) + ys(3)))/((ys(1) + ys(2))*(ys(1) + ys(2) + ys(3)) &
                                                   & *(ys(1) + ys(2) + ys(3) + ys(4)))

                            ys = s_cb(i + 4:i + 1:-1) - s_cb(i + 3:i:-1)
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, &
                                                   & 2) = (ys(4)*(ys(2)**2 + 4*ys(2)*ys(3) + 4*ys(2)*ys(4) + ys(1)*ys(2) + 3*ys(3) &
                                                   & **2 + 6*ys(3)*ys(4) + 2*ys(1)*ys(3) + 3*ys(4)**2 + 2*ys(1)*ys(4)))/((ys(3) &
                                                   & + ys(4))*(ys(2) + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, &
                                                   & 1) = -(ys(4)*(ys(3) + ys(4))*(ys(1)**2 + 3*ys(1)*ys(2) + 3*ys(1)*ys(3) &
                                                   & + 2*ys(1)*ys(4) + 3*ys(2)**2 + 6*ys(2)*ys(3) + 4*ys(2)*ys(4) + 3*ys(3)**2 &
                                                   & + 4*ys(3)*ys(4) + ys(4)**2))/((ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))*(ys(2) &
                                                   & + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4)))
                            poly_coef_cbL_${XYZ}$ (i + 1, 0, &
                                                   & 0) = (ys(4)*(ys(3) + ys(4))*(ys(2) + ys(3) + ys(4)))/((ys(1) + ys(2))*(ys(1) &
                                                   & + ys(2) + ys(3))*(ys(1) + ys(2) + ys(3) + ys(4)))

                            poly_coef_cbL_${XYZ}$ (i + 1,:,:) = -poly_coef_cbL_${XYZ}$ (i + 1,:,:)
                            ! Note: negative sign as the direction of taking the difference (dvd) is reversed

                            ys = s_cb(i - 2:i + 1) - s_cb(i - 3:i)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 0) = (4*ys(4)**2*(5*ys(1)**2*ys(2)**2 + 20*ys(1)**2*ys(2)*ys(3) + 15*ys(1) &
                                               & **2*ys(2)*ys(4) + 20*ys(1)**2*ys(3)**2 + 30*ys(1)**2*ys(3)*ys(4) + 60*ys(1) &
                                               & **2*ys(4)**2 + 10*ys(1)*ys(2)**3 + 60*ys(1)*ys(2)**2*ys(3) + 45*ys(1)*ys(2) &
                                               & **2*ys(4) + 110*ys(1)*ys(2)*ys(3)**2 + 165*ys(1)*ys(2)*ys(3)*ys(4) + 260*ys(1) &
                                               & *ys(2)*ys(4)**2 + 60*ys(1)*ys(3)**3 + 135*ys(1)*ys(3)**2*ys(4) + 400*ys(1)*ys(3) &
                                               & *ys(4)**2 + 225*ys(1)*ys(4)**3 + 5*ys(2)**4 + 40*ys(2)**3*ys(3) + 30*ys(2) &
                                               & **3*ys(4) + 110*ys(2)**2*ys(3)**2 + 165*ys(2)**2*ys(3)*ys(4) + 260*ys(2)**2*ys(4) &
                                               & **2 + 120*ys(2)*ys(3)**3 + 270*ys(2)*ys(3)**2*ys(4) + 800*ys(2)*ys(3)*ys(4)**2 &
                                               & + 450*ys(2)*ys(4)**3 + 45*ys(3)**4 + 135*ys(3)**3*ys(4) + 600*ys(3)**2*ys(4)**2 &
                                               & + 675*ys(3)*ys(4)**3 + 996*ys(4)**4))/(5*(ys(3) + ys(4))**2*(ys(2) + ys(3) &
                                               & + ys(4))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 1) = -(4*ys(4)**2*(10*ys(1)**3*ys(2)*ys(3) + 5*ys(1)**3*ys(2)*ys(4) + 20*ys(1) &
                                               & **3*ys(3)**2 + 25*ys(1)**3*ys(3)*ys(4) + 105*ys(1)**3*ys(4)**2 + 40*ys(1) &
                                               & **2*ys(2)**2*ys(3) + 20*ys(1)**2*ys(2)**2*ys(4) + 130*ys(1)**2*ys(2)*ys(3)**2 &
                                               & + 155*ys(1)**2*ys(2)*ys(3)*ys(4) + 535*ys(1)**2*ys(2)*ys(4)**2 + 90*ys(1) &
                                               & **2*ys(3)**3 + 165*ys(1)**2*ys(3)**2*ys(4) + 790*ys(1)**2*ys(3)*ys(4)**2 &
                                               & + 415*ys(1)**2*ys(4)**3 + 60*ys(1)*ys(2)**3*ys(3) + 30*ys(1)*ys(2)**3*ys(4) &
                                               & + 270*ys(1)*ys(2)**2*ys(3)**2 + 315*ys(1)*ys(2)**2*ys(3)*ys(4) + 975*ys(1)*ys(2) &
                                               & **2*ys(4)**2 + 360*ys(1)*ys(2)*ys(3)**3 + 645*ys(1)*ys(2)*ys(3)**2*ys(4) &
                                               & + 2850*ys(1)*ys(2)*ys(3)*ys(4)**2 + 1460*ys(1)*ys(2)*ys(4)**3 + 150*ys(1)*ys(3) &
                                               & **4 + 360*ys(1)*ys(3)**3*ys(4) + 2000*ys(1)*ys(3)**2*ys(4)**2 + 2005*ys(1)*ys(3) &
                                               & *ys(4)**3 + 2077*ys(1)*ys(4)**4 + 30*ys(2)**4*ys(3) + 15*ys(2)**4*ys(4) &
                                               & + 180*ys(2)**3*ys(3)**2 + 210*ys(2)**3*ys(3)*ys(4) + 650*ys(2)**3*ys(4)**2 &
                                               & + 360*ys(2)**2*ys(3)**3 + 645*ys(2)**2*ys(3)**2*ys(4) + 2850*ys(2)**2*ys(3)*ys(4) &
                                               & **2 + 1460*ys(2)**2*ys(4)**3 + 300*ys(2)*ys(3)**4 + 720*ys(2)*ys(3)**3*ys(4) &
                                               & + 4000*ys(2)*ys(3)**2*ys(4)**2 + 4010*ys(2)*ys(3)*ys(4)**3 + 4154*ys(2)*ys(4)**4 &
                                               & + 90*ys(3)**5 + 270*ys(3)**4*ys(4) + 1800*ys(3)**3*ys(4)**2 + 2655*ys(3)**2*ys(4) &
                                               & **3 + 4464*ys(3)*ys(4)**4 + 1767*ys(4)**5))/(5*(ys(2) + ys(3))*(ys(3) + ys(4)) &
                                               & *(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))**2*(ys(1) + ys(2) + ys(3) &
                                               & + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 2) = (4*ys(4)**2*(10*ys(2)**3*ys(3) + 5*ys(2)**3*ys(4) + 50*ys(2)**2*ys(3)**2 &
                                               & + 60*ys(2)**2*ys(3)*ys(4) + 10*ys(1)*ys(2)**2*ys(3) + 215*ys(2)**2*ys(4)**2 &
                                               & + 5*ys(1)*ys(2)**2*ys(4) + 70*ys(2)*ys(3)**3 + 130*ys(2)*ys(3)**2*ys(4) &
                                               & + 30*ys(1)*ys(2)*ys(3)**2 + 775*ys(2)*ys(3)*ys(4)**2 + 35*ys(1)*ys(2)*ys(3)*ys(4) &
                                               & + 415*ys(2)*ys(4)**3 + 110*ys(1)*ys(2)*ys(4)**2 + 30*ys(3)**4 + 75*ys(3)**3*ys(4) &
                                               & + 20*ys(1)*ys(3)**3 + 665*ys(3)**2*ys(4)**2 + 35*ys(1)*ys(3)**2*ys(4) + 725*ys(3) &
                                               & *ys(4)**3 + 220*ys(1)*ys(3)*ys(4)**2 + 1767*ys(4)**4 + 105*ys(1)*ys(4)**3)) &
                                               & /(5*(ys(1) + ys(2))*(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) &
                                               & + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 3) = (4*ys(4)**2*(5*ys(1)**4*ys(3)**2 + 5*ys(1)**4*ys(3)*ys(4) + 50*ys(1) &
                                               & **4*ys(4)**2 + 30*ys(1)**3*ys(2)*ys(3)**2 + 30*ys(1)**3*ys(2)*ys(3)*ys(4) &
                                               & + 300*ys(1)**3*ys(2)*ys(4)**2 + 30*ys(1)**3*ys(3)**3 + 45*ys(1)**3*ys(3)**2*ys(4) &
                                               & + 415*ys(1)**3*ys(3)*ys(4)**2 + 200*ys(1)**3*ys(4)**3 + 75*ys(1)**2*ys(2) &
                                               & **2*ys(3)**2 + 75*ys(1)**2*ys(2)**2*ys(3)*ys(4) + 750*ys(1)**2*ys(2)**2*ys(4)**2 &
                                               & + 150*ys(1)**2*ys(2)*ys(3)**3 + 225*ys(1)**2*ys(2)*ys(3)**2*ys(4) + 2075*ys(1) &
                                               & **2*ys(2)*ys(3)*ys(4)**2 + 1000*ys(1)**2*ys(2)*ys(4)**3 + 75*ys(1)**2*ys(3)**4 &
                                               & + 150*ys(1)**2*ys(3)**3*ys(4) + 1390*ys(1)**2*ys(3)**2*ys(4)**2 + 1315*ys(1) &
                                               & **2*ys(3)*ys(4)**3 + 1081*ys(1)**2*ys(4)**4 + 90*ys(1)*ys(2)**3*ys(3)**2 &
                                               & + 90*ys(1)*ys(2)**3*ys(3)*ys(4) + 900*ys(1)*ys(2)**3*ys(4)**2 + 270*ys(1)*ys(2) &
                                               & **2*ys(3)**3 + 405*ys(1)*ys(2)**2*ys(3)**2*ys(4) + 3735*ys(1)*ys(2)**2*ys(3) &
                                               & *ys(4)**2 + 1800*ys(1)*ys(2)**2*ys(4)**3 + 270*ys(1)*ys(2)*ys(3)**4 + 540*ys(1) &
                                               & *ys(2)*ys(3)**3*ys(4) + 5025*ys(1)*ys(2)*ys(3)**2*ys(4)**2 + 4755*ys(1)*ys(2) &
                                               & *ys(3)*ys(4)**3 + 4224*ys(1)*ys(2)*ys(4)**4 + 90*ys(1)*ys(3)**5 + 225*ys(1)*ys(3) &
                                               & **4*ys(4) + 2190*ys(1)*ys(3)**3*ys(4)**2 + 3060*ys(1)*ys(3)**2*ys(4)**3 &
                                               & + 4529*ys(1)*ys(3)*ys(4)**4 + 1762*ys(1)*ys(4)**5 + 45*ys(2)**4*ys(3)**2 &
                                               & + 45*ys(2)**4*ys(3)*ys(4) + 450*ys(2)**4*ys(4)**2 + 180*ys(2)**3*ys(3)**3 &
                                               & + 270*ys(2)**3*ys(3)**2*ys(4) + 2490*ys(2)**3*ys(3)*ys(4)**2 + 1200*ys(2) &
                                               & **3*ys(4)**3 + 270*ys(2)**2*ys(3)**4 + 540*ys(2)**2*ys(3)**3*ys(4) + 5025*ys(2) &
                                               & **2*ys(3)**2*ys(4)**2 + 4755*ys(2)**2*ys(3)*ys(4)**3 + 4224*ys(2)**2*ys(4)**4 &
                                               & + 180*ys(2)*ys(3)**5 + 450*ys(2)*ys(3)**4*ys(4) + 4380*ys(2)*ys(3)**3*ys(4)**2 &
                                               & + 6120*ys(2)*ys(3)**2*ys(4)**3 + 9058*ys(2)*ys(3)*ys(4)**4 + 3524*ys(2)*ys(4)**5 &
                                               & + 45*ys(3)**6 + 135*ys(3)**5*ys(4) + 1395*ys(3)**4*ys(4)**2 + 2565*ys(3)**3*ys(4) &
                                               & **3 + 4884*ys(3)**2*ys(4)**4 + 3624*ys(3)*ys(4)**5 + 831*ys(4)**6))/(5*(ys(2) &
                                               & + ys(3))**2*(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) + ys(4))**2*(ys(1) + ys(2) &
                                               & + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 4) = -(4*ys(4)**2*(10*ys(1)**2*ys(2)*ys(3)**2 + 10*ys(1)**2*ys(2)*ys(3)*ys(4) &
                                               & + 100*ys(1)**2*ys(2)*ys(4)**2 + 10*ys(1)**2*ys(3)**3 + 15*ys(1)**2*ys(3)**2*ys(4) &
                                               & + 205*ys(1)**2*ys(3)*ys(4)**2 + 100*ys(1)**2*ys(4)**3 + 30*ys(1)*ys(2)**2*ys(3) &
                                               & **2 + 30*ys(1)*ys(2)**2*ys(3)*ys(4) + 300*ys(1)*ys(2)**2*ys(4)**2 + 60*ys(1) &
                                               & *ys(2)*ys(3)**3 + 90*ys(1)*ys(2)*ys(3)**2*ys(4) + 1030*ys(1)*ys(2)*ys(3)*ys(4) &
                                               & **2 + 500*ys(1)*ys(2)*ys(4)**3 + 30*ys(1)*ys(3)**4 + 60*ys(1)*ys(3)**3*ys(4) &
                                               & + 835*ys(1)*ys(3)**2*ys(4)**2 + 805*ys(1)*ys(3)*ys(4)**3 + 1762*ys(1)*ys(4)**4 &
                                               & + 30*ys(2)**3*ys(3)**2 + 30*ys(2)**3*ys(3)*ys(4) + 300*ys(2)**3*ys(4)**2 &
                                               & + 90*ys(2)**2*ys(3)**3 + 135*ys(2)**2*ys(3)**2*ys(4) + 1445*ys(2)**2*ys(3)*ys(4) &
                                               & **2 + 700*ys(2)**2*ys(4)**3 + 90*ys(2)*ys(3)**4 + 180*ys(2)*ys(3)**3*ys(4) &
                                               & + 2205*ys(2)*ys(3)**2*ys(4)**2 + 2115*ys(2)*ys(3)*ys(4)**3 + 3624*ys(2)*ys(4)**4 &
                                               & + 30*ys(3)**5 + 75*ys(3)**4*ys(4) + 1060*ys(3)**3*ys(4)**2 + 1515*ys(3)**2*ys(4) &
                                               & **3 + 3824*ys(3)*ys(4)**4 + 1662*ys(4)**5))/(5*(ys(1) + ys(2))*(ys(2) + ys(3)) &
                                               & *(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) &
                                               & + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 3, &
                                               & 5) = (4*ys(4)**2*(5*ys(2)**2*ys(3)**2 + 5*ys(2)**2*ys(3)*ys(4) + 50*ys(2) &
                                               & **2*ys(4)**2 + 10*ys(2)*ys(3)**3 + 15*ys(2)*ys(3)**2*ys(4) + 205*ys(2)*ys(3) &
                                               & *ys(4)**2 + 100*ys(2)*ys(4)**3 + 5*ys(3)**4 + 10*ys(3)**3*ys(4) + 205*ys(3) &
                                               & **2*ys(4)**2 + 200*ys(3)*ys(4)**3 + 831*ys(4)**4))/(5*(ys(1) + ys(2))**2*(ys(1) &
                                               & + ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)

                            ys = s_cb(i - 1:i + 2) - s_cb(i - 2:i + 1)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 0) = (4*ys(3)**2*(5*ys(1)**2*ys(2)**2 + 5*ys(1)**2*ys(2)*ys(3) + 50*ys(1) &
                                               & **2*ys(3)**2 + 10*ys(1)*ys(2)**3 + 15*ys(1)*ys(2)**2*ys(3) + 205*ys(1)*ys(2) &
                                               & *ys(3)**2 + 100*ys(1)*ys(3)**3 + 5*ys(2)**4 + 10*ys(2)**3*ys(3) + 205*ys(2) &
                                               & **2*ys(3)**2 + 200*ys(2)*ys(3)**3 + 831*ys(3)**4))/(5*(ys(3) + ys(4))**2*(ys(2) &
                                               & + ys(3) + ys(4))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 1) = (4*ys(3)**2*(5*ys(1)**3*ys(2)*ys(3) + 10*ys(1)**3*ys(2)*ys(4) - 95*ys(1) &
                                               & **3*ys(3)**2 + 5*ys(1)**3*ys(3)*ys(4) + 20*ys(1)**2*ys(2)**2*ys(3) + 40*ys(1) &
                                               & **2*ys(2)**2*ys(4) - 465*ys(1)**2*ys(2)*ys(3)**2 + 55*ys(1)**2*ys(2)*ys(3)*ys(4) &
                                               & + 10*ys(1)**2*ys(2)*ys(4)**2 - 285*ys(1)**2*ys(3)**3 + 20*ys(1)**2*ys(3)**2*ys(4) &
                                               & + 5*ys(1)**2*ys(3)*ys(4)**2 + 30*ys(1)*ys(2)**3*ys(3) + 60*ys(1)*ys(2)**3*ys(4) &
                                               & - 825*ys(1)*ys(2)**2*ys(3)**2 + 135*ys(1)*ys(2)**2*ys(3)*ys(4) + 30*ys(1)*ys(2) &
                                               & **2*ys(4)**2 - 1040*ys(1)*ys(2)*ys(3)**3 + 100*ys(1)*ys(2)*ys(3)**2*ys(4) &
                                               & + 35*ys(1)*ys(2)*ys(3)*ys(4)**2 - 1847*ys(1)*ys(3)**4 + 125*ys(1)*ys(3)**3*ys(4) &
                                               & + 110*ys(1)*ys(3)**2*ys(4)**2 + 15*ys(2)**4*ys(3) + 30*ys(2)**4*ys(4) - 550*ys(2) &
                                               & **3*ys(3)**2 + 90*ys(2)**3*ys(3)*ys(4) + 20*ys(2)**3*ys(4)**2 - 1040*ys(2) &
                                               & **2*ys(3)**3 + 100*ys(2)**2*ys(3)**2*ys(4) + 35*ys(2)**2*ys(3)*ys(4)**2 &
                                               & - 3694*ys(2)*ys(3)**4 + 250*ys(2)*ys(3)**3*ys(4) + 220*ys(2)*ys(3)**2*ys(4)**2 &
                                               & - 3219*ys(3)**5 - 1452*ys(3)**4*ys(4) + 105*ys(3)**3*ys(4)**2))/(5*(ys(2) + ys(3) &
                                               & )*(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))**2*(ys(1) &
                                               & + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 2) = -(4*ys(3)**2*(5*ys(2)**3*ys(3) - 95*ys(2)*ys(3)**3 - 190*ys(2)**2*ys(3)**2 &
                                               & + 10*ys(2)**3*ys(4) + 100*ys(3)**3*ys(4) - 1562*ys(3)**4 - 95*ys(1)*ys(2)*ys(3) &
                                               & **2 + 5*ys(1)*ys(2)**2*ys(3) + 10*ys(1)*ys(2)**2*ys(4) + 100*ys(1)*ys(3)**2*ys(4) &
                                               & + 205*ys(2)*ys(3)**2*ys(4) + 15*ys(2)**2*ys(3)*ys(4) + 10*ys(1)*ys(2)*ys(3)*ys(4) &
                                               & ))/(5*(ys(1) + ys(2))*(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) &
                                               & + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 3) = (4*ys(3)**2*(50*ys(1)**4*ys(3)**2 + 5*ys(1)**4*ys(3)*ys(4) + 5*ys(1) &
                                               & **4*ys(4)**2 + 300*ys(1)**3*ys(2)*ys(3)**2 + 30*ys(1)**3*ys(2)*ys(3)*ys(4) &
                                               & + 30*ys(1)**3*ys(2)*ys(4)**2 + 200*ys(1)**3*ys(3)**3 + 25*ys(1)**3*ys(3)**2*ys(4) &
                                               & + 35*ys(1)**3*ys(3)*ys(4)**2 + 10*ys(1)**3*ys(4)**3 + 750*ys(1)**2*ys(2)**2*ys(3) &
                                               & **2 + 75*ys(1)**2*ys(2)**2*ys(3)*ys(4) + 75*ys(1)**2*ys(2)**2*ys(4)**2 &
                                               & + 1000*ys(1)**2*ys(2)*ys(3)**3 + 125*ys(1)**2*ys(2)*ys(3)**2*ys(4) + 175*ys(1) &
                                               & **2*ys(2)*ys(3)*ys(4)**2 + 50*ys(1)**2*ys(2)*ys(4)**3 + 1081*ys(1)**2*ys(3)**4 &
                                               & - 50*ys(1)**2*ys(3)**3*ys(4) - 10*ys(1)**2*ys(3)**2*ys(4)**2 + 45*ys(1)**2*ys(3) &
                                               & *ys(4)**3 + 5*ys(1)**2*ys(4)**4 + 900*ys(1)*ys(2)**3*ys(3)**2 + 90*ys(1)*ys(2) &
                                               & **3*ys(3)*ys(4) + 90*ys(1)*ys(2)**3*ys(4)**2 + 1800*ys(1)*ys(2)**2*ys(3)**3 &
                                               & + 225*ys(1)*ys(2)**2*ys(3)**2*ys(4) + 315*ys(1)*ys(2)**2*ys(3)*ys(4)**2 &
                                               & + 90*ys(1)*ys(2)**2*ys(4)**3 + 4224*ys(1)*ys(2)*ys(3)**4 - 120*ys(1)*ys(2)*ys(3) &
                                               & **3*ys(4) + 25*ys(1)*ys(2)*ys(3)**2*ys(4)**2 + 165*ys(1)*ys(2)*ys(3)*ys(4)**3 &
                                               & + 20*ys(1)*ys(2)*ys(4)**4 + 3324*ys(1)*ys(3)**5 + 1407*ys(1)*ys(3)**4*ys(4) &
                                               & - 100*ys(1)*ys(3)**3*ys(4)**2 + 70*ys(1)*ys(3)**2*ys(4)**3 + 15*ys(1)*ys(3)*ys(4) &
                                               & **4 + 450*ys(2)**4*ys(3)**2 + 45*ys(2)**4*ys(3)*ys(4) + 45*ys(2)**4*ys(4)**2 &
                                               & + 1200*ys(2)**3*ys(3)**3 + 150*ys(2)**3*ys(3)**2*ys(4) + 210*ys(2)**3*ys(3)*ys(4) &
                                               & **2 + 60*ys(2)**3*ys(4)**3 + 4224*ys(2)**2*ys(3)**4 - 120*ys(2)**2*ys(3)**3*ys(4) &
                                               & + 25*ys(2)**2*ys(3)**2*ys(4)**2 + 165*ys(2)**2*ys(3)*ys(4)**3 + 20*ys(2)**2*ys(4) &
                                               & **4 + 6648*ys(2)*ys(3)**5 + 2814*ys(2)*ys(3)**4*ys(4) - 200*ys(2)*ys(3)**3*ys(4) &
                                               & **2 + 140*ys(2)*ys(3)**2*ys(4)**3 + 30*ys(2)*ys(3)*ys(4)**4 + 3174*ys(3)**6 &
                                               & + 3039*ys(3)**5*ys(4) + 771*ys(3)**4*ys(4)**2 + 135*ys(3)**3*ys(4)**3 + 60*ys(3) &
                                               & **2*ys(4)**4))/(5*(ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) &
                                               & + ys(4))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 4) = -(4*ys(3)**2*(100*ys(1)**2*ys(2)*ys(3)**2 + 10*ys(1)**2*ys(2)*ys(3)*ys(4) &
                                               & + 10*ys(1)**2*ys(2)*ys(4)**2 - 95*ys(1)**2*ys(3)**2*ys(4) + 5*ys(1)**2*ys(3) &
                                               & *ys(4)**2 + 300*ys(1)*ys(2)**2*ys(3)**2 + 30*ys(1)*ys(2)**2*ys(3)*ys(4) &
                                               & + 30*ys(1)*ys(2)**2*ys(4)**2 + 200*ys(1)*ys(2)*ys(3)**3 - 260*ys(1)*ys(2)*ys(3) &
                                               & **2*ys(4) + 50*ys(1)*ys(2)*ys(3)*ys(4)**2 + 10*ys(1)*ys(2)*ys(4)**3 + 1562*ys(1) &
                                               & *ys(3)**4 - 190*ys(1)*ys(3)**3*ys(4) + 15*ys(1)*ys(3)**2*ys(4)**2 + 5*ys(1)*ys(3) &
                                               & *ys(4)**3 + 300*ys(2)**3*ys(3)**2 + 30*ys(2)**3*ys(3)*ys(4) + 30*ys(2)**3*ys(4) &
                                               & **2 + 400*ys(2)**2*ys(3)**3 - 235*ys(2)**2*ys(3)**2*ys(4) + 85*ys(2)**2*ys(3) &
                                               & *ys(4)**2 + 20*ys(2)**2*ys(4)**3 + 3224*ys(2)*ys(3)**4 - 460*ys(2)*ys(3)**3*ys(4) &
                                               & - 35*ys(2)*ys(3)**2*ys(4)**2 + 25*ys(2)*ys(3)*ys(4)**3 + 3124*ys(3)**5 &
                                               & + 1467*ys(3)**4*ys(4) + 110*ys(3)**3*ys(4)**2 + 105*ys(3)**2*ys(4)**3))/(5*(ys(1) &
                                               & + ys(2))*(ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) + ys(4)) &
                                               & *(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 2, &
                                               & 5) = (4*ys(3)**2*(50*ys(2)**2*ys(3)**2 + 5*ys(2)**2*ys(3)*ys(4) + 5*ys(2) &
                                               & **2*ys(4)**2 - 95*ys(2)*ys(3)**2*ys(4) + 5*ys(2)*ys(3)*ys(4)**2 + 781*ys(3)**4 &
                                               & + 50*ys(3)**2*ys(4)**2))/(5*(ys(1) + ys(2))**2*(ys(1) + ys(2) + ys(3))**2*(ys(1) &
                                               & + ys(2) + ys(3) + ys(4))**2)

                            ys = s_cb(i:i + 3) - s_cb(i - 1:i + 2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 0) = (4*ys(2)**2*(50*ys(1)**2*ys(2)**2 + 5*ys(1)**2*ys(2)*ys(3) + 5*ys(1) &
                                               & **2*ys(3)**2 - 95*ys(1)*ys(2)**2*ys(3) + 5*ys(1)*ys(2)*ys(3)**2 + 781*ys(2)**4 &
                                               & + 50*ys(2)**2*ys(3)**2))/(5*(ys(3) + ys(4))**2*(ys(2) + ys(3) + ys(4))**2*(ys(1) &
                                               & + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 1) = -(4*ys(2)**2*(105*ys(1)**3*ys(2)**2 + 25*ys(1)**3*ys(2)*ys(3) + 5*ys(1) &
                                               & **3*ys(2)*ys(4) + 20*ys(1)**3*ys(3)**2 + 10*ys(1)**3*ys(3)*ys(4) + 110*ys(1) &
                                               & **2*ys(2)**3 - 35*ys(1)**2*ys(2)**2*ys(3) + 15*ys(1)**2*ys(2)**2*ys(4) + 85*ys(1) &
                                               & **2*ys(2)*ys(3)**2 + 50*ys(1)**2*ys(2)*ys(3)*ys(4) + 5*ys(1)**2*ys(2)*ys(4)**2 &
                                               & + 30*ys(1)**2*ys(3)**3 + 30*ys(1)**2*ys(3)**2*ys(4) + 10*ys(1)**2*ys(3)*ys(4)**2 &
                                               & + 1467*ys(1)*ys(2)**4 - 460*ys(1)*ys(2)**3*ys(3) - 190*ys(1)*ys(2)**3*ys(4) &
                                               & - 235*ys(1)*ys(2)**2*ys(3)**2 - 260*ys(1)*ys(2)**2*ys(3)*ys(4) - 95*ys(1)*ys(2) &
                                               & **2*ys(4)**2 + 30*ys(1)*ys(2)*ys(3)**3 + 30*ys(1)*ys(2)*ys(3)**2*ys(4) + 10*ys(1) &
                                               & *ys(2)*ys(3)*ys(4)**2 + 3124*ys(2)**5 + 3224*ys(2)**4*ys(3) + 1562*ys(2)**4*ys(4) &
                                               & + 400*ys(2)**3*ys(3)**2 + 200*ys(2)**3*ys(3)*ys(4) + 300*ys(2)**2*ys(3)**3 &
                                               & + 300*ys(2)**2*ys(3)**2*ys(4) + 100*ys(2)**2*ys(3)*ys(4)**2))/(5*(ys(2) + ys(3)) &
                                               & *(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))**2*(ys(1) &
                                               & + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 2) = -(4*ys(2)**2*(100*ys(1)*ys(2)**3 - 190*ys(2)**2*ys(3)**2 + 10*ys(1)*ys(3) &
                                               & **3 + 5*ys(2)*ys(3)**3 - 95*ys(2)**3*ys(3) - 1562*ys(2)**4 + 15*ys(1)*ys(2)*ys(3) &
                                               & **2 + 205*ys(1)*ys(2)**2*ys(3) + 100*ys(1)*ys(2)**2*ys(4) + 10*ys(1)*ys(3) &
                                               & **2*ys(4) + 5*ys(2)*ys(3)**2*ys(4) - 95*ys(2)**2*ys(3)*ys(4) + 10*ys(1)*ys(2) &
                                               & *ys(3)*ys(4)))/(5*(ys(1) + ys(2))*(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) &
                                               & + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 3) = (4*ys(2)**2*(60*ys(1)**4*ys(2)**2 + 30*ys(1)**4*ys(2)*ys(3) + 15*ys(1) &
                                               & **4*ys(2)*ys(4) + 20*ys(1)**4*ys(3)**2 + 20*ys(1)**4*ys(3)*ys(4) + 5*ys(1) &
                                               & **4*ys(4)**2 + 135*ys(1)**3*ys(2)**3 + 140*ys(1)**3*ys(2)**2*ys(3) + 70*ys(1) &
                                               & **3*ys(2)**2*ys(4) + 165*ys(1)**3*ys(2)*ys(3)**2 + 165*ys(1)**3*ys(2)*ys(3)*ys(4) &
                                               & + 45*ys(1)**3*ys(2)*ys(4)**2 + 60*ys(1)**3*ys(3)**3 + 90*ys(1)**3*ys(3)**2*ys(4) &
                                               & + 50*ys(1)**3*ys(3)*ys(4)**2 + 10*ys(1)**3*ys(4)**3 + 771*ys(1)**2*ys(2)**4 &
                                               & - 200*ys(1)**2*ys(2)**3*ys(3) - 100*ys(1)**2*ys(2)**3*ys(4) + 25*ys(1)**2*ys(2) &
                                               & **2*ys(3)**2 + 25*ys(1)**2*ys(2)**2*ys(3)*ys(4) - 10*ys(1)**2*ys(2)**2*ys(4)**2 &
                                               & + 210*ys(1)**2*ys(2)*ys(3)**3 + 315*ys(1)**2*ys(2)*ys(3)**2*ys(4) + 175*ys(1) &
                                               & **2*ys(2)*ys(3)*ys(4)**2 + 35*ys(1)**2*ys(2)*ys(4)**3 + 45*ys(1)**2*ys(3)**4 &
                                               & + 90*ys(1)**2*ys(3)**3*ys(4) + 75*ys(1)**2*ys(3)**2*ys(4)**2 + 30*ys(1)**2*ys(3) &
                                               & *ys(4)**3 + 5*ys(1)**2*ys(4)**4 + 3039*ys(1)*ys(2)**5 + 2814*ys(1)*ys(2)**4*ys(3) &
                                               & + 1407*ys(1)*ys(2)**4*ys(4) - 120*ys(1)*ys(2)**3*ys(3)**2 - 120*ys(1)*ys(2) &
                                               & **3*ys(3)*ys(4) - 50*ys(1)*ys(2)**3*ys(4)**2 + 150*ys(1)*ys(2)**2*ys(3)**3 &
                                               & + 225*ys(1)*ys(2)**2*ys(3)**2*ys(4) + 125*ys(1)*ys(2)**2*ys(3)*ys(4)**2 &
                                               & + 25*ys(1)*ys(2)**2*ys(4)**3 + 45*ys(1)*ys(2)*ys(3)**4 + 90*ys(1)*ys(2)*ys(3) &
                                               & **3*ys(4) + 75*ys(1)*ys(2)*ys(3)**2*ys(4)**2 + 30*ys(1)*ys(2)*ys(3)*ys(4)**3 &
                                               & + 5*ys(1)*ys(2)*ys(4)**4 + 3174*ys(2)**6 + 6648*ys(2)**5*ys(3) + 3324*ys(2) &
                                               & **5*ys(4) + 4224*ys(2)**4*ys(3)**2 + 4224*ys(2)**4*ys(3)*ys(4) + 1081*ys(2) &
                                               & **4*ys(4)**2 + 1200*ys(2)**3*ys(3)**3 + 1800*ys(2)**3*ys(3)**2*ys(4) + 1000*ys(2) &
                                               & **3*ys(3)*ys(4)**2 + 200*ys(2)**3*ys(4)**3 + 450*ys(2)**2*ys(3)**4 + 900*ys(2) &
                                               & **2*ys(3)**3*ys(4) + 750*ys(2)**2*ys(3)**2*ys(4)**2 + 300*ys(2)**2*ys(3)*ys(4) &
                                               & **3 + 50*ys(2)**2*ys(4)**4))/(5*(ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3)) &
                                               & **2*(ys(2) + ys(3) + ys(4))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 4) = (4*ys(2)**2*(105*ys(1)**2*ys(2)**3 + 220*ys(1)**2*ys(2)**2*ys(3) + 110*ys(1) &
                                               & **2*ys(2)**2*ys(4) + 35*ys(1)**2*ys(2)*ys(3)**2 + 35*ys(1)**2*ys(2)*ys(3)*ys(4) &
                                               & + 5*ys(1)**2*ys(2)*ys(4)**2 + 20*ys(1)**2*ys(3)**3 + 30*ys(1)**2*ys(3)**2*ys(4) &
                                               & + 10*ys(1)**2*ys(3)*ys(4)**2 - 1452*ys(1)*ys(2)**4 + 250*ys(1)*ys(2)**3*ys(3) &
                                               & + 125*ys(1)*ys(2)**3*ys(4) + 100*ys(1)*ys(2)**2*ys(3)**2 + 100*ys(1)*ys(2) &
                                               & **2*ys(3)*ys(4) + 20*ys(1)*ys(2)**2*ys(4)**2 + 90*ys(1)*ys(2)*ys(3)**3 &
                                               & + 135*ys(1)*ys(2)*ys(3)**2*ys(4) + 55*ys(1)*ys(2)*ys(3)*ys(4)**2 + 5*ys(1)*ys(2) &
                                               & *ys(4)**3 + 30*ys(1)*ys(3)**4 + 60*ys(1)*ys(3)**3*ys(4) + 40*ys(1)*ys(3)**2*ys(4) &
                                               & **2 + 10*ys(1)*ys(3)*ys(4)**3 - 3219*ys(2)**5 - 3694*ys(2)**4*ys(3) - 1847*ys(2) &
                                               & **4*ys(4) - 1040*ys(2)**3*ys(3)**2 - 1040*ys(2)**3*ys(3)*ys(4) - 285*ys(2) &
                                               & **3*ys(4)**2 - 550*ys(2)**2*ys(3)**3 - 825*ys(2)**2*ys(3)**2*ys(4) - 465*ys(2) &
                                               & **2*ys(3)*ys(4)**2 - 95*ys(2)**2*ys(4)**3 + 15*ys(2)*ys(3)**4 + 30*ys(2)*ys(3) &
                                               & **3*ys(4) + 20*ys(2)*ys(3)**2*ys(4)**2 + 5*ys(2)*ys(3)*ys(4)**3))/(5*(ys(1) &
                                               & + ys(2))*(ys(2) + ys(3))*(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) + ys(4)) &
                                               & *(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 1, &
                                               & 5) = (4*ys(2)**2*(831*ys(2)**4 + 200*ys(2)**3*ys(3) + 100*ys(2)**3*ys(4) &
                                               & + 205*ys(2)**2*ys(3)**2 + 205*ys(2)**2*ys(3)*ys(4) + 50*ys(2)**2*ys(4)**2 &
                                               & + 10*ys(2)*ys(3)**3 + 15*ys(2)*ys(3)**2*ys(4) + 5*ys(2)*ys(3)*ys(4)**2 + 5*ys(3) &
                                               & **4 + 10*ys(3)**3*ys(4) + 5*ys(3)**2*ys(4)**2))/(5*(ys(1) + ys(2))**2*(ys(1) &
                                               & + ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)

                            ys = s_cb(i + 1:i + 4) - s_cb(i:i + 3)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 0) = (4*ys(1)**2*(831*ys(1)**4 + 200*ys(1)**3*ys(2) + 100*ys(1)**3*ys(3) &
                                               & + 205*ys(1)**2*ys(2)**2 + 205*ys(1)**2*ys(2)*ys(3) + 50*ys(1)**2*ys(3)**2 &
                                               & + 10*ys(1)*ys(2)**3 + 15*ys(1)*ys(2)**2*ys(3) + 5*ys(1)*ys(2)*ys(3)**2 + 5*ys(2) &
                                               & **4 + 10*ys(2)**3*ys(3) + 5*ys(2)**2*ys(3)**2))/(5*(ys(3) + ys(4))**2*(ys(2) &
                                               & + ys(3) + ys(4))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 1) = -(4*ys(1)**2*(1662*ys(1)**5 + 3824*ys(1)**4*ys(2) + 3624*ys(1)**4*ys(3) &
                                               & + 1762*ys(1)**4*ys(4) + 1515*ys(1)**3*ys(2)**2 + 2115*ys(1)**3*ys(2)*ys(3) &
                                               & + 805*ys(1)**3*ys(2)*ys(4) + 700*ys(1)**3*ys(3)**2 + 500*ys(1)**3*ys(3)*ys(4) &
                                               & + 100*ys(1)**3*ys(4)**2 + 1060*ys(1)**2*ys(2)**3 + 2205*ys(1)**2*ys(2)**2*ys(3) &
                                               & + 835*ys(1)**2*ys(2)**2*ys(4) + 1445*ys(1)**2*ys(2)*ys(3)**2 + 1030*ys(1) &
                                               & **2*ys(2)*ys(3)*ys(4) + 205*ys(1)**2*ys(2)*ys(4)**2 + 300*ys(1)**2*ys(3)**3 &
                                               & + 300*ys(1)**2*ys(3)**2*ys(4) + 100*ys(1)**2*ys(3)*ys(4)**2 + 75*ys(1)*ys(2)**4 &
                                               & + 180*ys(1)*ys(2)**3*ys(3) + 60*ys(1)*ys(2)**3*ys(4) + 135*ys(1)*ys(2)**2*ys(3) &
                                               & **2 + 90*ys(1)*ys(2)**2*ys(3)*ys(4) + 15*ys(1)*ys(2)**2*ys(4)**2 + 30*ys(1)*ys(2) &
                                               & *ys(3)**3 + 30*ys(1)*ys(2)*ys(3)**2*ys(4) + 10*ys(1)*ys(2)*ys(3)*ys(4)**2 &
                                               & + 30*ys(2)**5 + 90*ys(2)**4*ys(3) + 30*ys(2)**4*ys(4) + 90*ys(2)**3*ys(3)**2 &
                                               & + 60*ys(2)**3*ys(3)*ys(4) + 10*ys(2)**3*ys(4)**2 + 30*ys(2)**2*ys(3)**3 &
                                               & + 30*ys(2)**2*ys(3)**2*ys(4) + 10*ys(2)**2*ys(3)*ys(4)**2))/(5*(ys(2) + ys(3)) &
                                               & *(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3))*(ys(2) + ys(3) + ys(4))**2*(ys(1) &
                                               & + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 2) = (4*ys(1)**2*(1767*ys(1)**4 + 725*ys(1)**3*ys(2) + 415*ys(1)**3*ys(3) &
                                               & + 105*ys(4)*ys(1)**3 + 665*ys(1)**2*ys(2)**2 + 775*ys(1)**2*ys(2)*ys(3) &
                                               & + 220*ys(4)*ys(1)**2*ys(2) + 215*ys(1)**2*ys(3)**2 + 110*ys(4)*ys(1)**2*ys(3) &
                                               & + 75*ys(1)*ys(2)**3 + 130*ys(1)*ys(2)**2*ys(3) + 35*ys(4)*ys(1)*ys(2)**2 &
                                               & + 60*ys(1)*ys(2)*ys(3)**2 + 35*ys(4)*ys(1)*ys(2)*ys(3) + 5*ys(1)*ys(3)**3 &
                                               & + 5*ys(4)*ys(1)*ys(3)**2 + 30*ys(2)**4 + 70*ys(2)**3*ys(3) + 20*ys(4)*ys(2)**3 &
                                               & + 50*ys(2)**2*ys(3)**2 + 30*ys(4)*ys(2)**2*ys(3) + 10*ys(2)*ys(3)**3 + 10*ys(4) &
                                               & *ys(2)*ys(3)**2))/(5*(ys(1) + ys(2))*(ys(3) + ys(4))*(ys(1) + ys(2) + ys(3)) &
                                               & *(ys(2) + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 3) = (4*ys(1)**2*(831*ys(1)**6 + 3624*ys(1)**5*ys(2) + 3524*ys(1)**5*ys(3) &
                                               & + 1762*ys(1)**5*ys(4) + 4884*ys(1)**4*ys(2)**2 + 9058*ys(1)**4*ys(2)*ys(3) &
                                               & + 4529*ys(1)**4*ys(2)*ys(4) + 4224*ys(1)**4*ys(3)**2 + 4224*ys(1)**4*ys(3)*ys(4) &
                                               & + 1081*ys(1)**4*ys(4)**2 + 2565*ys(1)**3*ys(2)**3 + 6120*ys(1)**3*ys(2)**2*ys(3) &
                                               & + 3060*ys(1)**3*ys(2)**2*ys(4) + 4755*ys(1)**3*ys(2)*ys(3)**2 + 4755*ys(1) &
                                               & **3*ys(2)*ys(3)*ys(4) + 1315*ys(1)**3*ys(2)*ys(4)**2 + 1200*ys(1)**3*ys(3)**3 &
                                               & + 1800*ys(1)**3*ys(3)**2*ys(4) + 1000*ys(1)**3*ys(3)*ys(4)**2 + 200*ys(1) &
                                               & **3*ys(4)**3 + 1395*ys(1)**2*ys(2)**4 + 4380*ys(1)**2*ys(2)**3*ys(3) + 2190*ys(1) &
                                               & **2*ys(2)**3*ys(4) + 5025*ys(1)**2*ys(2)**2*ys(3)**2 + 5025*ys(1)**2*ys(2) &
                                               & **2*ys(3)*ys(4) + 1390*ys(1)**2*ys(2)**2*ys(4)**2 + 2490*ys(1)**2*ys(2)*ys(3)**3 &
                                               & + 3735*ys(1)**2*ys(2)*ys(3)**2*ys(4) + 2075*ys(1)**2*ys(2)*ys(3)*ys(4)**2 &
                                               & + 415*ys(1)**2*ys(2)*ys(4)**3 + 450*ys(1)**2*ys(3)**4 + 900*ys(1)**2*ys(3) &
                                               & **3*ys(4) + 750*ys(1)**2*ys(3)**2*ys(4)**2 + 300*ys(1)**2*ys(3)*ys(4)**3 &
                                               & + 50*ys(1)**2*ys(4)**4 + 135*ys(1)*ys(2)**5 + 450*ys(1)*ys(2)**4*ys(3) &
                                               & + 225*ys(1)*ys(2)**4*ys(4) + 540*ys(1)*ys(2)**3*ys(3)**2 + 540*ys(1)*ys(2) &
                                               & **3*ys(3)*ys(4) + 150*ys(1)*ys(2)**3*ys(4)**2 + 270*ys(1)*ys(2)**2*ys(3)**3 &
                                               & + 405*ys(1)*ys(2)**2*ys(3)**2*ys(4) + 225*ys(1)*ys(2)**2*ys(3)*ys(4)**2 &
                                               & + 45*ys(1)*ys(2)**2*ys(4)**3 + 45*ys(1)*ys(2)*ys(3)**4 + 90*ys(1)*ys(2)*ys(3) &
                                               & **3*ys(4) + 75*ys(1)*ys(2)*ys(3)**2*ys(4)**2 + 30*ys(1)*ys(2)*ys(3)*ys(4)**3 &
                                               & + 5*ys(1)*ys(2)*ys(4)**4 + 45*ys(2)**6 + 180*ys(2)**5*ys(3) + 90*ys(2)**5*ys(4) &
                                               & + 270*ys(2)**4*ys(3)**2 + 270*ys(2)**4*ys(3)*ys(4) + 75*ys(2)**4*ys(4)**2 &
                                               & + 180*ys(2)**3*ys(3)**3 + 270*ys(2)**3*ys(3)**2*ys(4) + 150*ys(2)**3*ys(3)*ys(4) &
                                               & **2 + 30*ys(2)**3*ys(4)**3 + 45*ys(2)**2*ys(3)**4 + 90*ys(2)**2*ys(3)**3*ys(4) &
                                               & + 75*ys(2)**2*ys(3)**2*ys(4)**2 + 30*ys(2)**2*ys(3)*ys(4)**3 + 5*ys(2)**2*ys(4) &
                                               & **4))/(5*(ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3))**2*(ys(2) + ys(3) + ys(4)) &
                                               & **2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 4) = -(4*ys(1)**2*(1767*ys(1)**5 + 4464*ys(1)**4*ys(2) + 4154*ys(1)**4*ys(3) &
                                               & + 2077*ys(1)**4*ys(4) + 2655*ys(1)**3*ys(2)**2 + 4010*ys(1)**3*ys(2)*ys(3) &
                                               & + 2005*ys(1)**3*ys(2)*ys(4) + 1460*ys(1)**3*ys(3)**2 + 1460*ys(1)**3*ys(3)*ys(4) &
                                               & + 415*ys(1)**3*ys(4)**2 + 1800*ys(1)**2*ys(2)**3 + 4000*ys(1)**2*ys(2)**2*ys(3) &
                                               & + 2000*ys(1)**2*ys(2)**2*ys(4) + 2850*ys(1)**2*ys(2)*ys(3)**2 + 2850*ys(1) &
                                               & **2*ys(2)*ys(3)*ys(4) + 790*ys(1)**2*ys(2)*ys(4)**2 + 650*ys(1)**2*ys(3)**3 &
                                               & + 975*ys(1)**2*ys(3)**2*ys(4) + 535*ys(1)**2*ys(3)*ys(4)**2 + 105*ys(1)**2*ys(4) &
                                               & **3 + 270*ys(1)*ys(2)**4 + 720*ys(1)*ys(2)**3*ys(3) + 360*ys(1)*ys(2)**3*ys(4) &
                                               & + 645*ys(1)*ys(2)**2*ys(3)**2 + 645*ys(1)*ys(2)**2*ys(3)*ys(4) + 165*ys(1)*ys(2) &
                                               & **2*ys(4)**2 + 210*ys(1)*ys(2)*ys(3)**3 + 315*ys(1)*ys(2)*ys(3)**2*ys(4) &
                                               & + 155*ys(1)*ys(2)*ys(3)*ys(4)**2 + 25*ys(1)*ys(2)*ys(4)**3 + 15*ys(1)*ys(3)**4 &
                                               & + 30*ys(1)*ys(3)**3*ys(4) + 20*ys(1)*ys(3)**2*ys(4)**2 + 5*ys(1)*ys(3)*ys(4)**3 &
                                               & + 90*ys(2)**5 + 300*ys(2)**4*ys(3) + 150*ys(2)**4*ys(4) + 360*ys(2)**3*ys(3)**2 &
                                               & + 360*ys(2)**3*ys(3)*ys(4) + 90*ys(2)**3*ys(4)**2 + 180*ys(2)**2*ys(3)**3 &
                                               & + 270*ys(2)**2*ys(3)**2*ys(4) + 130*ys(2)**2*ys(3)*ys(4)**2 + 20*ys(2)**2*ys(4) &
                                               & **3 + 30*ys(2)*ys(3)**4 + 60*ys(2)*ys(3)**3*ys(4) + 40*ys(2)*ys(3)**2*ys(4)**2 &
                                               & + 10*ys(2)*ys(3)*ys(4)**3))/(5*(ys(1) + ys(2))*(ys(2) + ys(3))*(ys(1) + ys(2) &
                                               & + ys(3))**2*(ys(2) + ys(3) + ys(4))*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                            beta_coef_${XYZ}$ (i + 1, 0, &
                                               & 5) = (4*ys(1)**2*(996*ys(1)**4 + 675*ys(1)**3*ys(2) + 450*ys(1)**3*ys(3) &
                                               & + 225*ys(1)**3*ys(4) + 600*ys(1)**2*ys(2)**2 + 800*ys(1)**2*ys(2)*ys(3) &
                                               & + 400*ys(1)**2*ys(2)*ys(4) + 260*ys(1)**2*ys(3)**2 + 260*ys(1)**2*ys(3)*ys(4) &
                                               & + 60*ys(1)**2*ys(4)**2 + 135*ys(1)*ys(2)**3 + 270*ys(1)*ys(2)**2*ys(3) &
                                               & + 135*ys(1)*ys(2)**2*ys(4) + 165*ys(1)*ys(2)*ys(3)**2 + 165*ys(1)*ys(2)*ys(3) &
                                               & *ys(4) + 30*ys(1)*ys(2)*ys(4)**2 + 30*ys(1)*ys(3)**3 + 45*ys(1)*ys(3)**2*ys(4) &
                                               & + 15*ys(1)*ys(3)*ys(4)**2 + 45*ys(2)**4 + 120*ys(2)**3*ys(3) + 60*ys(2)**3*ys(4) &
                                               & + 110*ys(2)**2*ys(3)**2 + 110*ys(2)**2*ys(3)*ys(4) + 20*ys(2)**2*ys(4)**2 &
                                               & + 40*ys(2)*ys(3)**3 + 60*ys(2)*ys(3)**2*ys(4) + 20*ys(2)*ys(3)*ys(4)**2 + 5*ys(3) &
                                               & **4 + 10*ys(3)**3*ys(4) + 5*ys(3)**2*ys(4)**2))/(5*(ys(1) + ys(2))**2*(ys(1) &
                                               & + ys(2) + ys(3))**2*(ys(1) + ys(2) + ys(3) + ys(4))**2)
                        end do
                    else
                        ! (Fu, et al., 2016) Table 2 (for right flux)
                        d_cbL_${XYZ}$ (0,:) = 18._wp/35._wp
                        d_cbL_${XYZ}$ (1,:) = 3._wp/35._wp
                        d_cbL_${XYZ}$ (2,:) = 9._wp/35._wp
                        d_cbL_${XYZ}$ (3,:) = 1._wp/35._wp
                        d_cbL_${XYZ}$ (4,:) = 4._wp/35._wp

                        d_cbR_${XYZ}$ (0,:) = 18._wp/35._wp
                        d_cbR_${XYZ}$ (1,:) = 9._wp/35._wp
                        d_cbR_${XYZ}$ (2,:) = 3._wp/35._wp
                        d_cbR_${XYZ}$ (3,:) = 4._wp/35._wp
                        d_cbR_${XYZ}$ (4,:) = 1._wp/35._wp
                    end if
                end if
            end if
        #:endfor

        ! Detect whether grid spacing is uniform (enables cancellation-free sum-of-squares beta). Tolerance uses sqrt(epsilon) so it
        ! works in both double and single precision: ~1.5e-8 relative in double, ~3.5e-4 in single - above FP noise, below real
        ! stretching.
        uniform_grid(weno_dir) = .true.
        h0 = (s_cb(s) - s_cb(0))/real(s, wp)
        do i = 0, s - 1
            if (abs((s_cb(i + 1) - s_cb(i)) - h0) > sqrt(epsilon(h0))*abs(h0)) then
                uniform_grid(weno_dir) = .false.
                exit
            end if
        end do

        if (weno_dir == 1) then
            $:GPU_UPDATE(device='[poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, uniform_grid]')
        else if (weno_dir == 2) then
            $:GPU_UPDATE(device='[poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, uniform_grid]')
        else
            $:GPU_UPDATE(device='[poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, uniform_grid]')
        end if

        ! Nullifying WENO coefficients and cell-boundary locations pointers

        nullify (s_cb)

    end subroutine s_compute_weno_coefficients

    subroutine s_pack_weno_input_arr(v_vf)

        type(scalar_field), dimension(1:), intent(in) :: v_vf
        integer                                       :: i, j, k, l, n_vars

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, v_size
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        v_rs_weno(j, k, l, i) = v_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pack_weno_input_arr

    !> Perform WENO reconstruction of left and right cell-boundary values from cell-averaged variables
    subroutine s_weno(v_vf, vL_rs_vf_x, vR_rs_vf_x, weno_dir, is1_weno_d, &

        & is2_weno_d, is3_weno_d)

        type(scalar_field), dimension(1:), intent(in)                                          :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf_x
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vR_rs_vf_x
        integer, intent(in)                                                                    :: weno_dir
        type(int_bounds_info), intent(in)                                                      :: is1_weno_d, is2_weno_d, is3_weno_d

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(-3:2) :: dvd
            real(wp), dimension(0:4)  :: poly
            real(wp), dimension(0:4)  :: alpha
            real(wp), dimension(0:4)  :: omega
            real(wp), dimension(0:4)  :: beta
            real(wp), dimension(0:4)  :: delta
        #:else
            real(wp), dimension(-weno_polyn:weno_polyn - 1) :: dvd
            real(wp), dimension(0:weno_num_stencils)        :: poly
            real(wp), dimension(0:weno_num_stencils)        :: alpha
            real(wp), dimension(0:weno_num_stencils)        :: omega
            real(wp), dimension(0:weno_num_stencils)        :: beta
            real(wp), dimension(0:weno_num_stencils)        :: delta
        #:endif
        real(wp), dimension(-3:3) :: v  !< temporary field value array for clarity (WENO7 only)
        real(wp)                  :: tau
        integer                   :: i, j, k, l, q
        real(wp)                  :: vp0, vp1, vp2, vp3, vm1, vm2, vm3

        is1_weno = is1_weno_d
        is2_weno = is2_weno_d
        is3_weno = is3_weno_d

        $:GPU_UPDATE(device='[is1_weno, is2_weno, is3_weno]')

        v_size = ubound(v_vf, 1)
        $:GPU_UPDATE(device='[v_size]')

        if (weno_order == 1) then
            if (weno_dir == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, v_size
                    do l = is3_weno%beg, is3_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do j = is1_weno%beg, is1_weno%end
                                vL_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                vR_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (weno_dir == 2) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, v_size
                    do l = is3_weno%beg, is3_weno%end
                        do j = is1_weno%beg, is1_weno%end
                            do k = is2_weno%beg, is2_weno%end
                                vL_rs_vf_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                                vR_rs_vf_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (weno_dir == 3) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, v_size
                    do j = is1_weno%beg, is1_weno%end
                        do k = is2_weno%beg, is2_weno%end
                            do l = is3_weno%beg, is3_weno%end
                                vL_rs_vf_x(l, k, j, i) = v_vf(i)%sf(l, k, j)
                                vR_rs_vf_x(l, k, j, i) = v_vf(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        if (weno_order /= 1) then
            call s_pack_weno_input_arr(v_vf)
        end if

        if (weno_order == 3) then
            #:for WENO_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1_weno', 'is2_weno', 'is3_weno'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2_weno', 'is1_weno', 'is3_weno'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3_weno', 'is2_weno', 'is1_weno')]
                #:set SV = STENCIL_VAR
                #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
                if (weno_dir == ${WENO_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[beta, dvd, poly, omega, alpha, tau, q, vp0, vp1, vm1]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                do i = 1, v_size
                                    ! reconstruct from left side

                                    alpha(:) = 0._wp

                                    vp0 = v_rs_weno(${SF('')}$, i)
                                    vm1 = v_rs_weno(${SF(' - 1')}$, i)
                                    vp1 = v_rs_weno(${SF(' + 1')}$, i)

                                    dvd(0) = vp1 - vp0
                                    dvd(-1) = vp0 - vm1

                                    poly(0) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 0, 0)*dvd(0)
                                    poly(1) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 1, 0)*dvd(-1)

                                    beta(0) = beta_coef_${XYZ}$ (${SV}$, 0, 0)*dvd(0)*dvd(0) + weno_eps
                                    beta(1) = beta_coef_${XYZ}$ (${SV}$, 1, 0)*dvd(-1)*dvd(-1) + weno_eps

                                    if (wenojs) then
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                        end do
                                    else if (mapped_weno) then
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                        end do
                                        omega = alpha/sum(alpha)
                                        do q = 0, weno_num_stencils
                                            alpha(q) = (d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbL_${XYZ}$ (q, &
                                                  & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbL_${XYZ}$ (q, &
                                                  & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbL_${XYZ}$ (q, ${SV}$))))
                                        end do
                                    else if (wenoz) then
                                        ! Borges, et al. (2008)
                                        tau = abs(beta(1) - beta(0))
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + tau/beta(q))
                                        end do
                                    end if
                                    omega = alpha/sum(alpha)

                                    vL_rs_vf_x(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)

                                    ! reconstruct from right side

                                    poly(0) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 0, 0)*dvd(0)
                                    poly(1) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 1, 0)*dvd(-1)

                                    if (wenojs) then
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                        end do
                                    else if (mapped_weno) then
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                        end do
                                        omega = alpha/sum(alpha)
                                        do q = 0, weno_num_stencils
                                            alpha(q) = (d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbR_${XYZ}$ (q, &
                                                  & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbR_${XYZ}$ (q, &
                                                  & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbR_${XYZ}$ (q, ${SV}$))))
                                        end do
                                    else if (wenoz) then
                                        do q = 0, weno_num_stencils
                                            alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + tau/beta(q))
                                        end do
                                    end if
                                    omega = alpha/sum(alpha)

                                    vR_rs_vf_x(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endfor
        end if
        if (weno_order == 5) then
            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 1
                #:for WENO_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1_weno', 'is2_weno', 'is3_weno'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2_weno', 'is1_weno', 'is3_weno'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3_weno', 'is2_weno', 'is1_weno')]
                    #:set SV = STENCIL_VAR
                    #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
                    if (weno_dir == ${WENO_DIR}$) then
                        $:GPU_PARALLEL_LOOP(collapse=3,private='[dvd, poly, beta, alpha, omega, tau, delta, q, vp0, vm1, vm2, &
                                            & vp1, vp2]')
                        do l = ${Z_BND}$%beg, ${Z_BND}$%end
                            do k = ${Y_BND}$%beg, ${Y_BND}$%end
                                do j = ${X_BND}$%beg, ${X_BND}$%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, v_size
                                        ! reconstruct from left side

                                        alpha(:) = 0._wp

                                        vp0 = v_rs_weno(${SF('')}$, i)
                                        vm1 = v_rs_weno(${SF(' - 1')}$, i)
                                        vm2 = v_rs_weno(${SF(' - 2')}$, i)
                                        vp1 = v_rs_weno(${SF(' + 1')}$, i)
                                        vp2 = v_rs_weno(${SF(' + 2')}$, i)

                                        dvd(1) = vp2 - vp1
                                        dvd(0) = vp1 - vp0
                                        dvd(-1) = vp0 - vm1
                                        dvd(-2) = vm1 - vm2

                                        poly(0) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 0, &
                                             & 0)*dvd(1) + poly_coef_cbL_${XYZ}$ (${SV}$, 0, 1)*dvd(0)
                                        poly(1) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 1, &
                                             & 0)*dvd(0) + poly_coef_cbL_${XYZ}$ (${SV}$, 1, 1)*dvd(-1)
                                        poly(2) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 2, &
                                             & 0)*dvd(-1) + poly_coef_cbL_${XYZ}$ (${SV}$, 2, 1)*dvd(-2)

                                        if (uniform_grid(${WENO_DIR}$)) then
                                            beta(0) = 13._wp/12._wp*(dvd(1) - dvd(0))**2 + 0.25_wp*(dvd(1) - 3._wp*dvd(0))**2 &
                                                 & + weno_eps
                                            beta(1) = 13._wp/12._wp*(dvd(0) - dvd(-1))**2 + 0.25_wp*(dvd(0) + dvd(-1))**2 + weno_eps
                                            beta(2) = 13._wp/12._wp*(dvd(-1) - dvd(-2))**2 + 0.25_wp*(3._wp*dvd(-1) - dvd(-2))**2 &
                                                 & + weno_eps
                                        else
                                            beta(0) = beta_coef_${XYZ}$ (${SV}$, 0, 0)*dvd(1)*dvd(1) + beta_coef_${XYZ}$ (${SV}$, &
                                                 & 0, 1)*dvd(1)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, 0, 2)*dvd(0)*dvd(0) + weno_eps
                                            beta(1) = beta_coef_${XYZ}$ (${SV}$, 1, 0)*dvd(0)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, &
                                                 & 1, 1)*dvd(0)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 1, &
                                                 & 2)*dvd(-1)*dvd(-1) + weno_eps
                                            beta(2) = beta_coef_${XYZ}$ (${SV}$, 2, &
                                                 & 0)*dvd(-1)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 2, &
                                                 & 1)*dvd(-1)*dvd(-2) + beta_coef_${XYZ}$ (${SV}$, 2, 2)*dvd(-2)*dvd(-2) + weno_eps
                                        end if

                                        if (wenojs) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                        else if (mapped_weno) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                            omega = alpha/sum(alpha)
                                            do q = 0, weno_num_stencils
                                                alpha(q) = (d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbL_${XYZ}$ (q, &
                                                      & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbL_${XYZ}$ (q, &
                                                      & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbL_${XYZ}$ (q, ${SV}$))))
                                            end do
                                        else if (wenoz) then
                                            ! Borges, et al. (2008)

                                            tau = abs(beta(2) - beta(0))  ! Equation 25
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + (tau/beta(q)))
                                                ! Equation 28 (note: weno_eps was already added to beta)
                                            end do
                                        else if (teno) then
                                            ! Fu, et al. (2016) Fu''s code: https://dx.doi.org/10.13140/RG.2.2.36250.34247
                                            tau = abs(beta(2) - beta(0))
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                alpha(q) = 1._wp + tau/beta(q)  ! Equation 22 (reuse alpha as gamma; pick C=1 & q=6)
                                                ! Equation 22 cont. (some CPU compilers cannot optimize x**6.0)
                                                alpha(q) = (alpha(q)**3._wp)**2._wp
                                            end do
                                            omega = alpha/sum(alpha)  ! Equation 25 (reuse omega as xi)

                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                if (omega(q) < teno_CT) then  ! Equation 26
                                                    delta(q) = 0._wp
                                                else
                                                    delta(q) = 1._wp
                                                end if
                                                alpha(q) = delta(q)*d_cbL_${XYZ}$ (q, ${SV}$)  ! Equation 27
                                            end do
                                        end if

                                        omega(0) = alpha(0)/(alpha(0) + alpha(1) + alpha(2))
                                        omega(1) = alpha(1)/(alpha(0) + alpha(1) + alpha(2))
                                        omega(2) = alpha(2)/(alpha(0) + alpha(1) + alpha(2))

                                        vL_rs_vf_x(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2)

                                        ! reconstruct from right side

                                        poly(0) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 0, &
                                             & 0)*dvd(1) + poly_coef_cbR_${XYZ}$ (${SV}$, 0, 1)*dvd(0)
                                        poly(1) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 1, &
                                             & 0)*dvd(0) + poly_coef_cbR_${XYZ}$ (${SV}$, 1, 1)*dvd(-1)
                                        poly(2) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 2, &
                                             & 0)*dvd(-1) + poly_coef_cbR_${XYZ}$ (${SV}$, 2, 1)*dvd(-2)

                                        if (wenojs) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                        else if (mapped_weno) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                            omega = alpha/sum(alpha)
                                            do q = 0, weno_num_stencils
                                                alpha(q) = (d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbR_${XYZ}$ (q, &
                                                      & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbR_${XYZ}$ (q, &
                                                      & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbR_${XYZ}$ (q, ${SV}$))))
                                            end do
                                        else if (wenoz) then
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + (tau/beta(q)))
                                            end do
                                        else if (teno) then
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                alpha(q) = delta(q)*d_cbR_${XYZ}$ (q, ${SV}$)
                                            end do
                                        end if

                                        omega(0) = alpha(0)/(alpha(0) + alpha(1) + alpha(2))
                                        omega(1) = alpha(1)/(alpha(0) + alpha(1) + alpha(2))
                                        omega(2) = alpha(2)/(alpha(0) + alpha(1) + alpha(2))

                                        vR_rs_vf_x(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        if (mp_weno) then
                            call s_preserve_monotonicity(v_rs_weno, vL_rs_vf_x, vR_rs_vf_x, weno_dir)
                        end if
                    end if
                #:endfor
            #:endif
        end if
        if (weno_order == 7) then
            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 2
                #:for WENO_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1_weno', 'is2_weno', 'is3_weno'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2_weno', 'is1_weno', 'is3_weno'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3_weno', 'is2_weno', 'is1_weno')]
                    #:set SV = STENCIL_VAR
                    #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
                    if (weno_dir == ${WENO_DIR}$) then
                        $:GPU_PARALLEL_LOOP(collapse=3,private='[poly, beta, alpha, omega, tau, delta, dvd, v, q, vp0, vp1, vp2, &
                                            & vp3, vm1, vm2, vm3]')
                        do l = ${Z_BND}$%beg, ${Z_BND}$%end
                            do k = ${Y_BND}$%beg, ${Y_BND}$%end
                                do j = ${X_BND}$%beg, ${X_BND}$%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, v_size
                                        alpha(:) = 0._wp

                                        vp0 = v_rs_weno(${SF('')}$, i)
                                        vm1 = v_rs_weno(${SF(' - 1')}$, i)
                                        vm2 = v_rs_weno(${SF(' - 2')}$, i)
                                        vm3 = v_rs_weno(${SF(' - 3')}$, i)
                                        vp1 = v_rs_weno(${SF(' + 1')}$, i)
                                        vp2 = v_rs_weno(${SF(' + 2')}$, i)
                                        vp3 = v_rs_weno(${SF(' + 3')}$, i)

                                        if (teno) then
                                            v(-3) = vm3
                                            v(-2) = vm2
                                            v(-1) = vm1
                                            v(0) = vp0
                                            v(1) = vp1
                                            v(2) = vp2
                                            v(3) = vp3
                                        end if

                                        if (.not. teno) then
                                            dvd(2) = vp3 - vp2
                                            dvd(1) = vp2 - vp1
                                            dvd(0) = vp1 - vp0
                                            dvd(-1) = vp0 - vm1
                                            dvd(-2) = vm1 - vm2
                                            dvd(-3) = vm2 - vm3

                                            poly(3) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 0, &
                                                 & 0)*dvd(2) + poly_coef_cbL_${XYZ}$ (${SV}$, 0, &
                                                 & 1)*dvd(1) + poly_coef_cbL_${XYZ}$ (${SV}$, 0, 2)*dvd(0)
                                            poly(2) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 1, &
                                                 & 0)*dvd(1) + poly_coef_cbL_${XYZ}$ (${SV}$, 1, &
                                                 & 1)*dvd(0) + poly_coef_cbL_${XYZ}$ (${SV}$, 1, 2)*dvd(-1)
                                            poly(1) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 2, &
                                                 & 0)*dvd(0) + poly_coef_cbL_${XYZ}$ (${SV}$, 2, &
                                                 & 1)*dvd(-1) + poly_coef_cbL_${XYZ}$ (${SV}$, 2, 2)*dvd(-2)
                                            poly(0) = vp0 + poly_coef_cbL_${XYZ}$ (${SV}$, 3, &
                                                 & 0)*dvd(-1) + poly_coef_cbL_${XYZ}$ (${SV}$, 3, &
                                                 & 1)*dvd(-2) + poly_coef_cbL_${XYZ}$ (${SV}$, 3, 2)*dvd(-3)
                                        else
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                ! (Fu, et al., 2016) Table 1 Note: Unlike TENO5, TENO7 stencils differ from WENO7
                                                ! stencils See Figure 2 (right) for right-sided flux (at i+1/2) Here we need the
                                                ! left-sided flux, so we flip the weights with respect to the x=i point But we need
                                                ! to keep the stencil order to reuse the beta coefficients
                                                poly(0) = (2._wp*v(-1) + 5._wp*v(0) - 1._wp*v(1))/6._wp
                                                poly(1) = (11._wp*v(0) - 7._wp*v(1) + 2._wp*v(2))/6._wp
                                                poly(2) = (-1._wp*v(-2) + 5._wp*v(-1) + 2._wp*v(0))/6._wp
                                                poly(3) = (25._wp*v(0) - 23._wp*v(1) + 13._wp*v(2) - 3._wp*v(3))/12._wp
                                                poly(4) = (1._wp*v(-3) - 5._wp*v(-2) + 13._wp*v(-1) + 3._wp*v(0))/12._wp
                                            #:endif
                                        end if

                                        if (.not. teno) then
                                            beta(3) = beta_coef_${XYZ}$ (${SV}$, 0, 0)*dvd(2)*dvd(2) + beta_coef_${XYZ}$ (${SV}$, &
                                                 & 0, 1)*dvd(2)*dvd(1) + beta_coef_${XYZ}$ (${SV}$, 0, &
                                                 & 2)*dvd(2)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, 0, &
                                                 & 3)*dvd(1)*dvd(1) + beta_coef_${XYZ}$ (${SV}$, 0, &
                                                 & 4)*dvd(1)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, 0, 5)*dvd(0)*dvd(0) + weno_eps

                                            beta(2) = beta_coef_${XYZ}$ (${SV}$, 1, 0)*dvd(1)*dvd(1) + beta_coef_${XYZ}$ (${SV}$, &
                                                 & 1, 1)*dvd(1)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, 1, &
                                                 & 2)*dvd(1)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 1, &
                                                 & 3)*dvd(0)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, 1, &
                                                 & 4)*dvd(0)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 1, 5)*dvd(-1)*dvd(-1) + weno_eps

                                            beta(1) = beta_coef_${XYZ}$ (${SV}$, 2, 0)*dvd(0)*dvd(0) + beta_coef_${XYZ}$ (${SV}$, &
                                                 & 2, 1)*dvd(0)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 2, &
                                                 & 2)*dvd(0)*dvd(-2) + beta_coef_${XYZ}$ (${SV}$, 2, &
                                                 & 3)*dvd(-1)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 2, &
                                                 & 4)*dvd(-1)*dvd(-2) + beta_coef_${XYZ}$ (${SV}$, 2, 5)*dvd(-2)*dvd(-2) + weno_eps

                                            beta(0) = beta_coef_${XYZ}$ (${SV}$, 3, &
                                                 & 0)*dvd(-1)*dvd(-1) + beta_coef_${XYZ}$ (${SV}$, 3, &
                                                 & 1)*dvd(-1)*dvd(-2) + beta_coef_${XYZ}$ (${SV}$, 3, &
                                                 & 2)*dvd(-1)*dvd(-3) + beta_coef_${XYZ}$ (${SV}$, 3, &
                                                 & 3)*dvd(-2)*dvd(-2) + beta_coef_${XYZ}$ (${SV}$, 3, &
                                                 & 4)*dvd(-2)*dvd(-3) + beta_coef_${XYZ}$ (${SV}$, 3, 5)*dvd(-3)*dvd(-3) + weno_eps
                                        else
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                ! High-Order Low-Dissipation Targeted ENO Schemes for Ideal Magnetohydrodynamics (Fu
                                                ! & Tang, 2019) Section 3.2
                                                beta(0) = 13._wp/12._wp*(v(-1) - 2._wp*v(0) + v(1))**2._wp + ((v(-1) - v(1)) &
                                                     & **2._wp)/4._wp + weno_eps
                                                beta(1) = 13._wp/12._wp*(v(0) - 2._wp*v(1) + v(2))**2._wp + ((3._wp*v(0) &
                                                     & - 4._wp*v(1) + v(2))**2._wp)/4._wp + weno_eps
                                                beta(2) = 13._wp/12._wp*(v(-2) - 2._wp*v(-1) + v(0))**2._wp + ((v(-2) &
                                                     & - 4._wp*v(-1) + 3._wp*v(0))**2._wp)/4._wp + weno_eps

                                                beta(3) = (v(0)*(2107._wp*v(0) - 9402._wp*v(1) + 7042._wp*v(2) - 1854._wp*v(3)) &
                                                     & + v(1)*(11003._wp*v(1) - 17246._wp*v(2) + 4642._wp*v(3)) + v(2) &
                                                     & *(7043._wp*v(2) - 3882._wp*v(3)) + v(3)*(547._wp*v(3)))/240._wp + weno_eps

                                                beta(4) = (v(-3)*(547._wp*v(-3) - 3882._wp*v(-2) + 4642._wp*v(-1) - 1854._wp*v(0)) &
                                                     & + v(-2)*(7043._wp*v(-2) - 17246._wp*v(-1) + 7042._wp*v(0)) + v(-1) &
                                                     & *(11003._wp*v(-1) - 9402._wp*v(0)) + v(0)*(2107._wp*v(0)))/240._wp + weno_eps
                                            #:endif
                                        end if

                                        if (wenojs) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                        else if (mapped_weno) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                            omega = alpha/sum(alpha)
                                            do q = 0, weno_num_stencils
                                                alpha(q) = (d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbL_${XYZ}$ (q, &
                                                      & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbL_${XYZ}$ (q, &
                                                      & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbL_${XYZ}$ (q, ${SV}$))))
                                            end do
                                        else if (wenoz) then
                                            ! Castro, et al. (2010) Don & Borges (2013) also helps
                                            tau = abs(beta(3) - beta(0))  ! Equation 50
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                ! wenoz_q = 2,3,4 for stability
                                                alpha(q) = d_cbL_${XYZ}$ (q, ${SV}$)*(1._wp + (tau/beta(q))**wenoz_q)
                                            end do
                                        else if (teno) then
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                tau = abs(beta(4) - beta(3))  ! Note the reordering of stencils
                                                alpha = 1._wp + tau/beta
                                                alpha = (alpha**3._wp)**2._wp  ! some CPU compilers cannot optimize x**6.0
                                                omega = alpha/sum(alpha)

                                                $:GPU_LOOP(parallelism='[seq]')
                                                do q = 0, weno_num_stencils
                                                    if (omega(q) < teno_CT) then  ! Equation 26
                                                        delta(q) = 0._wp
                                                    else
                                                        delta(q) = 1._wp
                                                    end if
                                                    alpha(q) = delta(q)*d_cbL_${XYZ}$ (q, ${SV}$)  ! Equation 27
                                                end do
                                            #:endif
                                        end if

                                        omega = alpha/sum(alpha)

                                        vL_rs_vf_x(j, k, l, &
                                                   & i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2) + omega(3)*poly(3)

                                        if (teno) then
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                vL_rs_vf_x(j, k, l, i) = vL_rs_vf_x(j, k, l, i) + omega(4)*poly(4)
                                            #:endif
                                        end if

                                        if (.not. teno) then
                                            poly(3) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 0, &
                                                 & 0)*dvd(2) + poly_coef_cbR_${XYZ}$ (${SV}$, 0, &
                                                 & 1)*dvd(1) + poly_coef_cbR_${XYZ}$ (${SV}$, 0, 2)*dvd(0)
                                            poly(2) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 1, &
                                                 & 0)*dvd(1) + poly_coef_cbR_${XYZ}$ (${SV}$, 1, &
                                                 & 1)*dvd(0) + poly_coef_cbR_${XYZ}$ (${SV}$, 1, 2)*dvd(-1)
                                            poly(1) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 2, &
                                                 & 0)*dvd(0) + poly_coef_cbR_${XYZ}$ (${SV}$, 2, &
                                                 & 1)*dvd(-1) + poly_coef_cbR_${XYZ}$ (${SV}$, 2, 2)*dvd(-2)
                                            poly(0) = vp0 + poly_coef_cbR_${XYZ}$ (${SV}$, 3, &
                                                 & 0)*dvd(-1) + poly_coef_cbR_${XYZ}$ (${SV}$, 3, &
                                                 & 1)*dvd(-2) + poly_coef_cbR_${XYZ}$ (${SV}$, 3, 2)*dvd(-3)
                                        else
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                poly(0) = (-1._wp*v(-1) + 5._wp*v(0) + 2._wp*v(1))/6._wp
                                                poly(1) = (2._wp*v(0) + 5._wp*v(1) - 1._wp*v(2))/6._wp
                                                poly(2) = (2._wp*v(-2) - 7._wp*v(-1) + 11._wp*v(0))/6._wp
                                                poly(3) = (3._wp*v(0) + 13._wp*v(1) - 5._wp*v(2) + 1._wp*v(3))/12._wp
                                                poly(4) = (-3._wp*v(-3) + 13._wp*v(-2) - 23._wp*v(-1) + 25._wp*v(0))/12._wp
                                            #:endif
                                        end if

                                        if (wenojs) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                        else if (mapped_weno) then
                                            do q = 0, weno_num_stencils
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)/(beta(q)**2._wp)
                                            end do
                                            omega = alpha/sum(alpha)
                                            do q = 0, weno_num_stencils
                                                alpha(q) = (d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + d_cbR_${XYZ}$ (q, &
                                                      & ${SV}$) - 3._wp*omega(q)) + omega(q)**2._wp)*(omega(q)/(d_cbR_${XYZ}$ (q, &
                                                      & ${SV}$)**2._wp + omega(q)*(1._wp - 2._wp*d_cbR_${XYZ}$ (q, ${SV}$))))
                                            end do
                                        else if (wenoz) then
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                ! wenoz_q = 2,3,4 for stability
                                                alpha(q) = d_cbR_${XYZ}$ (q, ${SV}$)*(1._wp + (tau/beta(q))**wenoz_q)
                                            end do
                                        else if (teno) then
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 0, weno_num_stencils
                                                alpha(q) = delta(q)*d_cbR_${XYZ}$ (q, ${SV}$)
                                            end do
                                        end if

                                        omega = alpha/sum(alpha)

                                        vR_rs_vf_x(j, k, l, &
                                                   & i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2) + omega(3)*poly(3)

                                        if (teno) then
                                            #:if not MFC_CASE_OPTIMIZATION or weno_num_stencils > 3
                                                vR_rs_vf_x(j, k, l, i) = vR_rs_vf_x(j, k, l, i) + omega(4)*poly(4)
                                            #:endif
                                        end if
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                #:endfor
            #:endif
        end if

        if (int_comp > 0 .and. v_size >= eqn_idx%adv%end) then
            call nvtxStartRange("WENO-INTCOMP")
            call s_thinc_compression(v_rs_weno, vL_rs_vf_x, vR_rs_vf_x, weno_dir, is1_weno, is2_weno, is3_weno)
            call nvtxEndRange()
        end if

    end subroutine s_weno

    !> Enforce monotonicity-preserving bounds on the WENO reconstruction
    subroutine s_preserve_monotonicity(v_rs_ws, vL_rs_vf, vR_rs_vf, weno_dir)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: v_rs_ws
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf, vR_rs_vf
        integer, intent(in) :: weno_dir
        integer :: i, j, k, l
        real(wp), dimension(-1:1) :: d  !< Curvature measures at the zone centers
        real(wp) :: d_MD, d_LC          !< Median (md) curvature and large curvature (LC) measures
        ! The left and right upper bounds (UL), medians, large curvatures, minima, and maxima of the WENO-reconstructed values of
        ! the cell- average variables.
        real(wp)            :: vL_UL, vR_UL
        real(wp)            :: vL_MD, vR_MD
        real(wp)            :: vL_LC, vR_LC
        real(wp)            :: vL_min, vR_min
        real(wp)            :: vL_max, vR_max
        real(wp), parameter :: alpha = 2._wp       !< Max CFL stability parameter (CFL < 1/(1+alpha))
        real(wp), parameter :: beta = 4._wp/3._wp  !< Local curvature freedom parameter
        real(wp), parameter :: alpha_mp = 2._wp
        real(wp), parameter :: beta_mp = 4._wp/3._wp
        real(wp)            :: vp0, vp1, vp2, vm1, vm2

        #:for WENO_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1_weno', 'is2_weno', 'is3_weno'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2_weno', 'is1_weno', 'is3_weno'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3_weno', 'is2_weno', 'is1_weno')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (weno_dir == ${WENO_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[d, vp0, vp1, vp2, vm1, vm2]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            do i = 1, v_size
                                ! Second-order undivided differences for curvature estimation

                                vp0 = v_rs_ws(${SF('')}$, i)
                                vm1 = v_rs_ws(${SF(' - 1')}$, i)
                                vm2 = v_rs_ws(${SF(' - 2')}$, i)
                                vp1 = v_rs_ws(${SF(' + 1')}$, i)
                                vp2 = v_rs_ws(${SF(' + 2')}$, i)

                                d(-1) = vp0 + vm2 - vm1*2._wp
                                d(0) = vp1 + vm1 - vp0*2._wp
                                d(1) = vp2 + vp0 - vp1*2._wp

                                ! Median function for oscillation detection
                                d_MD = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1)))*abs((sign(1._wp, &
                                        & 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1)))*(sign(1._wp, &
                                        & 4._wp*d(-1) - d(0)) + sign(1._wp, d(0))))*min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), &
                                        & abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                                d_LC = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0)))*abs((sign(1._wp, &
                                        & 4._wp*d(0) - d(1)) + sign(1._wp, d(0)))*(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, &
                                        & d(1))))*min(abs(4._wp*d(0) - d(1)), abs(d(0)), abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                                vL_UL = vp0 - (vp1 - vp0)*alpha_mp

                                vL_MD = (vp0 + vm1 - d_MD)*5.e-1_wp

                                vL_LC = vp0 - (vp1 - vp0)*5.e-1_wp + beta_mp*d_LC

                                vL_min = max(min(vp0, vm1, vL_MD), min(vp0, vL_UL, vL_LC))

                                vL_max = min(max(vp0, vm1, vL_MD), max(vp0, vL_UL, vL_LC))

                                vL_rs_vf(j, k, l, i) = vL_rs_vf(j, k, l, i) + (sign(5.e-1_wp, vL_min - vL_rs_vf(j, k, l, &
                                         & i)) + sign(5.e-1_wp, vL_max - vL_rs_vf(j, k, l, i)))*min(abs(vL_min - vL_rs_vf(j, k, &
                                         & l, i)), abs(vL_max - vL_rs_vf(j, k, l, i)))
                                ! END: Left Monotonicity Preserving Bound

                                ! Right Monotonicity Preserving Bound
                                d(-1) = vp0 + vm2 - vm1*2._wp
                                d(0) = vp1 + vm1 - vp0*2._wp
                                d(1) = vp2 + vp0 - vp1*2._wp

                                d_MD = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0)))*abs((sign(1._wp, &
                                        & 4._wp*d(0) - d(1)) + sign(1._wp, d(0)))*(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, &
                                        & d(1))))*min(abs(4._wp*d(0) - d(1)), abs(d(0)), abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                                d_LC = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1)))*abs((sign(1._wp, &
                                        & 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1)))*(sign(1._wp, &
                                        & 4._wp*d(-1) - d(0)) + sign(1._wp, d(0))))*min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), &
                                        & abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                                vR_UL = vp0 + (vp0 - vm1)*alpha_mp

                                vR_MD = (vp0 + vp1 - d_MD)*5.e-1_wp

                                vR_LC = vp0 + (vp0 - vm1)*5.e-1_wp + beta_mp*d_LC

                                vR_min = max(min(vp0, vp1, vR_MD), min(vp0, vR_UL, vR_LC))

                                vR_max = min(max(vp0, vp1, vR_MD), max(vp0, vR_UL, vR_LC))

                                vR_rs_vf(j, k, l, i) = vR_rs_vf(j, k, l, i) + (sign(5.e-1_wp, vR_min - vR_rs_vf(j, k, l, &
                                         & i)) + sign(5.e-1_wp, vR_max - vR_rs_vf(j, k, l, i)))*min(abs(vR_min - vR_rs_vf(j, k, &
                                         & l, i)), abs(vR_max - vR_rs_vf(j, k, l, i)))
                                ! END: Right Monotonicity Preserving Bound
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_preserve_monotonicity

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_weno_module()

        if (weno_order == 1) return

        ! Deallocating the WENO-stencil of the WENO-reconstructed variables

        @:DEALLOCATE(v_rs_weno)

        ! Deallocating WENO coefficients in x-direction
        @:DEALLOCATE(poly_coef_cbL_x, poly_coef_cbR_x)
        @:DEALLOCATE(d_cbL_x, d_cbR_x)
        @:DEALLOCATE(beta_coef_x)

        ! Deallocating WENO coefficients in y-direction
        if (n == 0) return

        @:DEALLOCATE(poly_coef_cbL_y, poly_coef_cbR_y)
        @:DEALLOCATE(d_cbL_y, d_cbR_y)
        @:DEALLOCATE(beta_coef_y)

        ! Deallocating WENO coefficients in z-direction
        if (p == 0) return

        @:DEALLOCATE(poly_coef_cbL_z, poly_coef_cbR_z)
        @:DEALLOCATE(d_cbL_z, d_cbR_z)
        @:DEALLOCATE(beta_coef_z)

    end subroutine s_finalize_weno_module

end module m_weno
