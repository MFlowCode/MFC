!>
!! @file m_cbc.f90
!! @brief Contains module m_cbc

!> @brief The module features a large database of characteristic boundary
!!              conditions (CBC) for the Euler system of equations. This system
!!              is augmented by the appropriate advection equations utilized to
!!              capture the material interfaces. The closure is achieved by the
!!              stiffened equation of state and mixture relations. At this time,
!!              the following CBC are available:
!!                           1) Slip Wall
!!                           2) Nonreflecting Subsonic Buffer
!!                           3) Nonreflecting Subsonic Inflow
!!                           4) Nonreflecting Subsonic Outflow
!!                           5) Force-Free Subsonic Outflow
!!                           6) Constant Pressure Subsonic Outflow
!!                           7) Supersonic Inflow
!!                           8) Supersonic Outflow
!!              Please refer to Thompson (1987, 1990) for detailed descriptions.

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'

module m_cbc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_compute_cbc
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_cbc_module, s_cbc, s_finalize_cbc_module

    !! The cell-average primitive variables. They are obtained by reshaping (RS)
    !! q_prim_vf in the coordinate direction normal to the domain boundary along
    !! which the CBC is applied.
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsz_vf)
    !$acc declare link(q_prim_rsx_vf, q_prim_rsy_vf, q_prim_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsz_vf
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), F_rs_vf, F_src_rs_vf)
    !$acc declare link(F_rs_vf, F_src_rs_vf)
#else
    type(scalar_field), allocatable, dimension(:) :: F_rs_vf, F_src_rs_vf !<
#endif
    !! Cell-average fluxes (src - source). These are directly determined from the
    !! cell-average primitive variables, q_prims_rs_vf, and not a Riemann solver.

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), F_rsx_vf, F_src_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), F_rsy_vf, F_src_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), F_rsz_vf, F_src_rsz_vf)
    !$acc declare link(F_rsx_vf, F_src_rsx_vf, F_rsy_vf, F_src_rsy_vf, F_rsz_vf, F_src_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: F_rsx_vf, F_src_rsx_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: F_rsy_vf, F_src_rsy_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: F_rsz_vf, F_src_rsz_vf !<
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsx_vf, flux_src_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsy_vf, flux_src_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsz_vf, flux_src_rsz_vf)
    !$acc declare link(flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf, flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsx_vf, flux_src_rsx_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsy_vf, flux_src_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsz_vf, flux_src_rsz_vf
#endif

    real(kind(0d0)) :: c           !< Cell averaged speed of sound
    real(kind(0d0)), dimension(2) :: Re          !< Cell averaged Reynolds numbers
    !$acc declare create(c, Re)

    real(kind(0d0)) :: dpres_ds !< Spatial derivatives in s-dir of pressure
!$acc declare create(dpres_ds)
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), ds)
    !$acc declare link(ds)
#else
    real(kind(0d0)), allocatable, dimension(:) :: ds !< Cell-width distribution in the s-direction
#endif

    ! CBC Coefficients =========================================================
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), fd_coef_x)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), fd_coef_y)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), fd_coef_z)
    !$acc declare link(fd_coef_x, fd_coef_y, fd_coef_z)
#else
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coef_x !< Finite diff. coefficients x-dir
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coef_y !< Finite diff. coefficients y-dir
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coef_z !< Finite diff. coefficients z-dir
#endif
    !! The first dimension identifies the location of a coefficient in the FD
    !! formula, while the last dimension denotes the location of the CBC.

    ! Bug with NVHPC when using nullified pointers in a declare create
    !    real(kind(0d0)), pointer, dimension(:, :) :: fd_coef => null()
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), pi_coef_x)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), pi_coef_y)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), pi_coef_z)
    !$acc declare link(pi_coef_x, pi_coef_y, pi_coef_z)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: pi_coef_x !< Polynomial interpolant coefficients in x-dir
    real(kind(0d0)), allocatable, dimension(:, :, :) :: pi_coef_y !< Polynomial interpolant coefficients in y-dir
    real(kind(0d0)), allocatable, dimension(:, :, :) :: pi_coef_z !< Polynomial interpolant coefficients in z-dir
#endif
    !! The first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the location of the CBC.

    ! ==========================================================================

    type(int_bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions
    !$acc declare create(is1, is2, is3)

    integer :: dj
    integer :: bcxb, bcxe, bcyb, bcye, bczb, bcze
    integer :: cbc_dir, cbc_loc
!$acc declare create(dj, bcxb, bcxe, bcyb, bcye, bczb, bcze, cbc_dir, cbc_loc)

#ifndef CRAY_ACC_WAR
!$acc declare create(q_prim_rsx_vf, q_prim_rsy_vf, q_prim_rsz_vf,  F_rsx_vf, F_src_rsx_vf,flux_rsx_vf, flux_src_rsx_vf, &
!$acc                 F_rsy_vf, F_src_rsy_vf,flux_rsy_vf, flux_src_rsy_vf, F_rsz_vf, F_src_rsz_vf,flux_rsz_vf, flux_src_rsz_vf,Re, &
!$acc                 ds,fd_coef_x,fd_coef_y,fd_coef_z,      &
!$acc                 pi_coef_x,pi_coef_y,pi_coef_z)
#endif

contains

    @:s_compute_speed_of_sound()

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_cbc_module

        integer :: i
        logical :: is_cbc

        call s_any_cbc_boundaries(is_cbc)

        if (is_cbc .eqv. .false.) return

        if (n == 0) then
            is2%beg = 0

        else
            is2%beg = -buff_size
        end if

        is2%end = n - is2%beg

        if (p == 0) then
            is3%beg = 0

        else
            is3%beg = -buff_size
        end if
        is3%end = p - is3%beg

        @:ALLOCATE_GLOBAL(q_prim_rsx_vf(0:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))

        if (weno_order > 1) then

            @:ALLOCATE_GLOBAL(F_rsx_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:adv_idx%end))

            @:ALLOCATE_GLOBAL(F_src_rsx_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        @:ALLOCATE_GLOBAL(flux_rsx_vf(-1:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:adv_idx%end))

        @:ALLOCATE_GLOBAL(flux_src_rsx_vf(-1:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        if (n > 0) then

            if (m == 0) then
                is2%beg = 0

            else
                is2%beg = -buff_size
            end if

            is2%end = m - is2%beg

            if (p == 0) then
                is3%beg = 0

            else
                is3%beg = -buff_size
            end if
            is3%end = p - is3%beg

            @:ALLOCATE_GLOBAL(q_prim_rsy_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1) then

                @:ALLOCATE_GLOBAL(F_rsy_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, 1:adv_idx%end))

                @:ALLOCATE_GLOBAL(F_src_rsy_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, adv_idx%beg:adv_idx%end))

            end if

            @:ALLOCATE_GLOBAL(flux_rsy_vf(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:adv_idx%end))

            @:ALLOCATE_GLOBAL(flux_src_rsy_vf(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        if (p > 0) then

            if (n == 0) then
                is2%beg = 0

            else
                is2%beg = -buff_size
            end if

            is2%end = n - is2%beg

            if (m == 0) then
                is3%beg = 0

            else
                is3%beg = -buff_size
            end if
            is3%end = m - is3%beg

            @:ALLOCATE_GLOBAL(q_prim_rsz_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1) then

                @:ALLOCATE_GLOBAL(F_rsz_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, 1:adv_idx%end))

                @:ALLOCATE_GLOBAL(F_src_rsz_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, adv_idx%beg:adv_idx%end))

            end if

            @:ALLOCATE_GLOBAL(flux_rsz_vf(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:adv_idx%end))

            @:ALLOCATE_GLOBAL(flux_src_rsz_vf(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        ! Allocating the cell-width distribution in the s-direction
        @:ALLOCATE_GLOBAL(ds(0:buff_size))

        ! Allocating/Computing CBC Coefficients in x-direction =============
        if (all((/bc_x%beg, bc_x%end/) <= -5) .and. all((/bc_x%beg, bc_x%end/) >= -13)) then

            @:ALLOCATE_GLOBAL(fd_coef_x(0:buff_size, -1:1))

            if (weno_order > 1) then
                @:ALLOCATE_GLOBAL(pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
            end if

            call s_compute_cbc_coefficients(1, -1)
            call s_compute_cbc_coefficients(1, 1)

        elseif (bc_x%beg <= -5 .and. bc_x%beg >= -13) then

            @:ALLOCATE_GLOBAL(fd_coef_x(0:buff_size, -1:-1))

            if (weno_order > 1) then
                @:ALLOCATE_GLOBAL(pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
            end if

            call s_compute_cbc_coefficients(1, -1)

        elseif (bc_x%end <= -5 .and. bc_x%end >= -13) then

            @:ALLOCATE_GLOBAL(fd_coef_x(0:buff_size, 1:1))

            if (weno_order > 1) then
                @:ALLOCATE_GLOBAL(pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
            end if

            call s_compute_cbc_coefficients(1, 1)

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in y-direction =============
        if (n > 0) then

            if (all((/bc_y%beg, bc_y%end/) <= -5) .and. all((/bc_y%beg, bc_y%end/) >= -13)) then

                @:ALLOCATE_GLOBAL(fd_coef_y(0:buff_size, -1:1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(2, -1)
                call s_compute_cbc_coefficients(2, 1)

            elseif (bc_y%beg <= -5 .and. bc_y%beg >= -13) then

                @:ALLOCATE_GLOBAL(fd_coef_y(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(2, -1)

            elseif (bc_y%end <= -5 .and. bc_y%end >= -13) then

                @:ALLOCATE_GLOBAL(fd_coef_y(0:buff_size, 1:1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(2, 1)

            end if

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in z-direction =============
        if (p > 0) then

            if (all((/bc_z%beg, bc_z%end/) <= -5) .and. all((/bc_z%beg, bc_z%end/) >= -13)) then

                @:ALLOCATE_GLOBAL(fd_coef_z(0:buff_size, -1:1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(3, -1)
                call s_compute_cbc_coefficients(3, 1)

            elseif (bc_z%beg <= -5 .and. bc_z%beg >= -13) then

                @:ALLOCATE_GLOBAL(fd_coef_z(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(3, -1)

            elseif (bc_z%end <= -5 .and. bc_z%end >= -13) then

                @:ALLOCATE_GLOBAL(fd_coef_z(0:buff_size, 1:1))

                if (weno_order > 1) then
                    @:ALLOCATE_GLOBAL(pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(3, 1)

            end if

        end if
        ! ==================================================================

        !$acc update device(fd_coef_x, fd_coef_y, fd_coef_z, pi_coef_x, pi_coef_y, pi_coef_z)

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        bcxb = bc_x%beg
        bcxe = bc_x%end

        !$acc update device(bcxb, bcxe)

        if (n > 0) then
            bcyb = bc_y%beg
            bcye = bc_y%end

            !$acc update device(bcyb, bcye)
        end if

        if (p > 0) then
            bczb = bc_z%beg
            bcze = bc_z%end

            !$acc update device(bczb, bcze)
        end if

    end subroutine s_initialize_cbc_module

    !>  Compute CBC coefficients
        !!  @param cbc_dir_in CBC coordinate direction
        !!  @param cbc_loc_in CBC coordinate location
    subroutine s_compute_cbc_coefficients(cbc_dir_in, cbc_loc_in)
        ! Description: The purpose of this subroutine is to compute the grid
        !              dependent FD and PI coefficients, or CBC coefficients,
        !              provided the CBC coordinate direction and location.

        ! CBC coordinate direction and location
        integer, intent(in) :: cbc_dir_in, cbc_loc_in

        ! Cell-boundary locations in the s-direction
        real(kind(0d0)), dimension(0:buff_size + 1) :: s_cb

        ! Generic loop iterator
        integer :: i

        ! Associating CBC coefficients pointers
        call s_associate_cbc_coefficients_pointers(cbc_dir_in, cbc_loc_in)

        ! Determining the cell-boundary locations in the s-direction
        s_cb(0) = 0d0

        do i = 0, buff_size
            s_cb(i + 1) = s_cb(i) + ds(i)
        end do

        ! Computing CBC1 Coefficients ======================================
        #:for CBC_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (cbc_dir_in == ${CBC_DIR}$) then
                if (weno_order == 1) then

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0d0
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -2d0/(ds(0) + ds(1))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -fd_coef_${XYZ}$ (0, cbc_loc_in)

                    ! ==================================================================

                    ! Computing CBC2 Coefficients ======================================
                elseif (weno_order == 3) then

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0d0
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -6d0/(3d0*ds(0) + 2d0*ds(1) - ds(2))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -4d0*fd_coef_${XYZ}$ (0, cbc_loc_in)/3d0
                    fd_coef_${XYZ}$ (2, cbc_loc_in) = fd_coef_${XYZ}$ (0, cbc_loc_in)/3d0

                    pi_coef_${XYZ}$ (0, 0, cbc_loc_in) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

                    ! ==================================================================

                    ! Computing CBC4 Coefficients ======================================
                else

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0d0
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -50d0/(25d0*ds(0) + 2d0*ds(1) &
                                                             - 1d1*ds(2) + 1d1*ds(3) &
                                                             - 3d0*ds(4))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -48d0*fd_coef_${XYZ}$ (0, cbc_loc_in)/25d0
                    fd_coef_${XYZ}$ (2, cbc_loc_in) = 36d0*fd_coef_${XYZ}$ (0, cbc_loc_in)/25d0
                    fd_coef_${XYZ}$ (3, cbc_loc_in) = -16d0*fd_coef_${XYZ}$ (0, cbc_loc_in)/25d0
                    fd_coef_${XYZ}$ (4, cbc_loc_in) = 3d0*fd_coef_${XYZ}$ (0, cbc_loc_in)/25d0

                    pi_coef_${XYZ}$ (0, 0, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                         (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                               (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
                    pi_coef_${XYZ}$ (0, 1, cbc_loc_in) = &
                        ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                         ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                          (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                               (s_cb(4) - s_cb(1)))))/ &
                        ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                         (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                    pi_coef_${XYZ}$ (0, 2, cbc_loc_in) = &
                        (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                             (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                                    (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                        ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                         (s_cb(0) - s_cb(4)))
                    pi_coef_${XYZ}$ (1, 0, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                         (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                               (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
                    pi_coef_${XYZ}$ (1, 1, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                         ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                          (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                               (s_cb(2) - s_cb(4)))))/ &
                        ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                         (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                    pi_coef_${XYZ}$ (1, 2, cbc_loc_in) = &
                        ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                         (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                               (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

                end if
            end if
        #:endfor

        ! END: Computing CBC4 Coefficients =================================

        ! Nullifying CBC coefficients

    end subroutine s_compute_cbc_coefficients

    !!  The goal of the procedure is to associate the FD and PI
    !!      coefficients, or CBC coefficients, with the appropriate
    !!      targets, based on the coordinate direction and location
    !!      of the CBC.
    !!  @param cbc_dir_in CBC coordinate direction
    !!  @param cbc_loc_in CBC coordinate location
    subroutine s_associate_cbc_coefficients_pointers(cbc_dir_in, cbc_loc_in)

        integer, intent(in) :: cbc_dir_in, cbc_loc_in

        integer :: i !< Generic loop iterator

        ! Associating CBC Coefficients in x-direction ======================
        if (cbc_dir_in == 1) then

            !fd_coef => fd_coef_x; if (weno_order > 1) pi_coef => pi_coef_x

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dx(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dx(m - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in y-direction ======================
        elseif (cbc_dir_in == 2) then

            !fd_coef => fd_coef_y; if (weno_order > 1) pi_coef => pi_coef_y

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dy(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dy(n - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in z-direction ======================
        else

            !fd_coef => fd_coef_z; if (weno_order > 1) pi_coef => pi_coef_z

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dz(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dz(p - i)
                end do
            end if

        end if

        !$acc update device(ds)

        ! ==================================================================

    end subroutine s_associate_cbc_coefficients_pointers

    !>  The following is the implementation of the CBC based on
        !!      the work of Thompson (1987, 1990) on hyperbolic systems.
        !!      The CBC is indirectly applied in the computation of the
        !!      right-hand-side (RHS) near the relevant domain boundary
        !!      through the modification of the fluxes.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir_norm CBC coordinate direction
        !!  @param cbc_loc_norm CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                     cbc_dir_norm, cbc_loc_norm, &
                     ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf

        integer, intent(in) :: cbc_dir_norm, cbc_loc_norm

        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! First-order time derivatives of the partial densities, density,
        ! velocity, pressure, advection variables, and the specific heat
        ! ratio and liquid stiffness functions
        real(kind(0d0)), dimension(num_fluids) :: dalpha_rho_dt
        real(kind(0d0)) :: drho_dt
        real(kind(0d0)), dimension(num_dims) :: dvel_dt
        real(kind(0d0)) :: dpres_dt
        real(kind(0d0)), dimension(num_fluids) :: dadv_dt
        real(kind(0d0)) :: dgamma_dt
        real(kind(0d0)) :: dpi_inf_dt
        real(kind(0d0)) :: dqv_dt
        real(kind(0d0)), dimension(contxe) :: alpha_rho, dalpha_rho_ds, mf
        real(kind(0d0)), dimension(2) :: Re_cbc
        real(kind(0d0)), dimension(num_dims) :: vel, dvel_ds
        real(kind(0d0)), dimension(num_fluids) :: adv, dadv_ds
        real(kind(0d0)), dimension(sys_size) :: L
        real(kind(0d0)), dimension(3) :: lambda

        real(kind(0d0)) :: rho         !< Cell averaged density
        real(kind(0d0)) :: pres        !< Cell averaged pressure
        real(kind(0d0)) :: E           !< Cell averaged energy
        real(kind(0d0)) :: H           !< Cell averaged enthalpy
        real(kind(0d0)) :: gamma       !< Cell averaged specific heat ratio
        real(kind(0d0)) :: pi_inf      !< Cell averaged liquid stiffness
        real(kind(0d0)) :: qv          !< Cell averaged fluid reference energy
        real(kind(0d0)) :: c

        real(kind(0d0)) :: vel_K_sum, vel_dv_dt_sum

        integer :: i, j, k, r, q !< Generic loop iterators

        real(kind(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for Wood mixture sound speed

        ! Reshaping of inputted data and association of the FD and PI
        ! coefficients, or CBC coefficients, respectively, hinging on
        ! selected CBC coordinate direction

        cbc_dir = cbc_dir_norm
        cbc_loc = cbc_loc_norm

        !$acc update device(cbc_dir, cbc_loc)

        call s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                              ix, iy, iz)

        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)

        #:for CBC_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (cbc_dir == ${CBC_DIR}$) then

                ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
                if (weno_order == 3) then

                    call s_convert_primitive_to_flux_variables(q_prim_rs${XYZ}$_vf, &
                                                               F_rs${XYZ}$_vf, &
                                                               F_src_rs${XYZ}$_vf, &
                                                               is1, is2, is3, starty, startz)

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = 1, advxe
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                flux_rs${XYZ}$_vf(0, k, r, i) = F_rs${XYZ}$_vf(0, k, r, i) &
                                                                + pi_coef_${XYZ}$ (0, 0, cbc_loc)* &
                                                                (F_rs${XYZ}$_vf(1, k, r, i) - &
                                                                 F_rs${XYZ}$_vf(0, k, r, i))
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = advxb, advxe
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                flux_src_rs${XYZ}$_vf(0, k, r, i) = F_src_rs${XYZ}$_vf(0, k, r, i) + &
                                                                    (F_src_rs${XYZ}$_vf(1, k, r, i) - &
                                                                     F_src_rs${XYZ}$_vf(0, k, r, i)) &
                                                                    *pi_coef_${XYZ}$ (0, 0, cbc_loc)
                            end do
                        end do
                    end do
                    ! ==================================================================

                    ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
                elseif (weno_order == 5) then
                    call s_convert_primitive_to_flux_variables(q_prim_rs${XYZ}$_vf, &
                                                               F_rs${XYZ}$_vf, &
                                                               F_src_rs${XYZ}$_vf, &
                                                               is1, is2, is3, starty, startz)

                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, advxe
                        do j = 0, 1
                            do r = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    flux_rs${XYZ}$_vf(j, k, r, i) = F_rs${XYZ}$_vf(j, k, r, i) &
                                                                    + pi_coef_${XYZ}$ (j, 0, cbc_loc)* &
                                                                    (F_rs${XYZ}$_vf(3, k, r, i) - &
                                                                     F_rs${XYZ}$_vf(2, k, r, i)) &
                                                                    + pi_coef_${XYZ}$ (j, 1, cbc_loc)* &
                                                                    (F_rs${XYZ}$_vf(2, k, r, i) - &
                                                                     F_rs${XYZ}$_vf(1, k, r, i)) &
                                                                    + pi_coef_${XYZ}$ (j, 2, cbc_loc)* &
                                                                    (F_rs${XYZ}$_vf(1, k, r, i) - &
                                                                     F_rs${XYZ}$_vf(0, k, r, i))
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = advxb, advxe
                        do j = 0, 1
                            do r = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    flux_src_rs${XYZ}$_vf(j, k, r, i) = F_src_rs${XYZ}$_vf(j, k, r, i) + &
                                                                        (F_src_rs${XYZ}$_vf(3, k, r, i) - &
                                                                         F_src_rs${XYZ}$_vf(2, k, r, i)) &
                                                                        *pi_coef_${XYZ}$ (j, 0, cbc_loc) + &
                                                                        (F_src_rs${XYZ}$_vf(2, k, r, i) - &
                                                                         F_src_rs${XYZ}$_vf(1, k, r, i)) &
                                                                        *pi_coef_${XYZ}$ (j, 1, cbc_loc) + &
                                                                        (F_src_rs${XYZ}$_vf(1, k, r, i) - &
                                                                         F_src_rs${XYZ}$_vf(0, k, r, i)) &
                                                                        *pi_coef_${XYZ}$ (j, 2, cbc_loc)
                                end do
                            end do
                        end do
                    end do

                end if
                ! ==================================================================

                ! FD2 or FD4 of RHS at j = 0 =======================================
                !$acc parallel loop collapse(2) gang vector default(present) private(alpha_rho, vel, adv, mf, dvel_ds, dadv_ds, Re_cbc, dalpha_rho_ds,dvel_dt, dadv_dt, dalpha_rho_dt,L, lambda)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end

                        ! Transferring the Primitive Variables =======================
                        !$acc loop seq
                        do i = 1, contxe
                            alpha_rho(i) = q_prim_rs${XYZ}$_vf(0, k, r, i)
                        end do

                        !$acc loop seq
                        do i = 1, num_dims
                            vel(i) = q_prim_rs${XYZ}$_vf(0, k, r, contxe + i)
                        end do

                        vel_K_sum = 0d0
                        !$acc loop seq
                        do i = 1, num_dims
                            vel_K_sum = vel_K_sum + vel(i)**2d0
                        end do

                        pres = q_prim_rs${XYZ}$_vf(0, k, r, E_idx)

                        !$acc loop seq
                        do i = 1, advxe - E_idx
                            adv(i) = q_prim_rs${XYZ}$_vf(0, k, r, E_idx + i)
                        end do

                        if (bubbles) then
                            call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, adv, alpha_rho, Re_cbc, 0, k, r)
                        else
                            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, adv, alpha_rho, Re_cbc, 0, k, r)
                        end if

                        !$acc loop seq
                        do i = 1, contxe
                            mf(i) = alpha_rho(i)/rho
                        end do

                        E = gamma*pres + pi_inf + 5d-1*rho*vel_K_sum

                        H = (E + pres)/rho

                        ! Compute mixture sound speed
                        call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_K_sum, c)
                        ! ============================================================

                        ! First-Order Spatial Derivatives of Primitive Variables =====

                        !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = 0d0
                        end do

                        !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = 0d0
                        end do

                        dpres_ds = 0d0
                        !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = 0d0
                        end do

                        !$acc loop seq
                        do j = 0, buff_size

                            !$acc loop seq
                            do i = 1, contxe
                                dalpha_rho_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, i)* &
                                                   fd_coef_${XYZ}$ (j, cbc_loc) + &
                                                   dalpha_rho_ds(i)
                            end do
                            !$acc loop seq
                            do i = 1, num_dims
                                dvel_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, contxe + i)* &
                                             fd_coef_${XYZ}$ (j, cbc_loc) + &
                                             dvel_ds(i)
                            end do

                            dpres_ds = q_prim_rs${XYZ}$_vf(j, k, r, E_idx)* &
                                       fd_coef_${XYZ}$ (j, cbc_loc) + &
                                       dpres_ds
                            !$acc loop seq
                            do i = 1, advxe - E_idx
                                dadv_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, E_idx + i)* &
                                             fd_coef_${XYZ}$ (j, cbc_loc) + &
                                             dadv_ds(i)
                            end do
                        end do
                        ! ============================================================

                        ! First-Order Temporal Derivatives of Primitive Variables ====
                        lambda(1) = vel(dir_idx(1)) - c
                        lambda(2) = vel(dir_idx(1))
                        lambda(3) = vel(dir_idx(1)) + c

                        if ((cbc_loc == -1 .and. bc${XYZ}$b == -5) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -5)) then
                            call s_compute_slip_wall_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -6) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -6)) then
                            call s_compute_nonreflecting_subsonic_buffer_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -7) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -7)) then
                            call s_compute_nonreflecting_subsonic_inflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -8) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -8)) then
                            call s_compute_nonreflecting_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -9) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -9)) then
                            call s_compute_force_free_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -10) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -10)) then
                            call s_compute_constant_pressure_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == -11) .or. (cbc_loc == 1 .and. bc${XYZ}$e == -11)) then
                            call s_compute_supersonic_inflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else
                            call s_compute_supersonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        end if

                        ! Be careful about the cylindrical coordinate!
                        if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                            dpres_dt = -5d-1*(L(advxe) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                       /y_cc(n)
                        else
                            dpres_dt = -5d-1*(L(advxe) + L(1))
                        end if

                        !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_dt(i) = &
                                -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                        end do

                        !$acc loop seq
                        do i = 1, num_dims
                            dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                                  (L(1) - L(advxe))/(2d0*rho*c) + &
                                                  (dir_flg(dir_idx(i)) - 1d0)* &
                                                  L(momxb + i - 1)
                        end do

                        vel_dv_dt_sum = 0d0
                        !$acc loop seq
                        do i = 1, num_dims
                            vel_dv_dt_sum = vel_dv_dt_sum + vel(i)*dvel_dt(i)
                        end do

                        ! The treatment of void fraction source is unclear
                        if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                            !$acc loop seq
                            do i = 1, advxe - E_idx
                                dadv_dt(i) = -L(momxe + i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                            end do
                        else
                            !$acc loop seq
                            do i = 1, advxe - E_idx
                                dadv_dt(i) = -L(momxe + i)
                            end do
                        end if

                        drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0; dqv_dt = 0d0

                        if (model_eqns == 1) then
                            drho_dt = dalpha_rho_dt(1)
                            dgamma_dt = dadv_dt(1)
                            dpi_inf_dt = dadv_dt(2)
                        else
                            !$acc loop seq
                            do i = 1, num_fluids
                                drho_dt = drho_dt + dalpha_rho_dt(i)
                                dgamma_dt = dgamma_dt + dadv_dt(i)*gammas(i)
                                dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_infs(i)
                                dqv_dt = dqv_dt + dalpha_rho_dt(i)*qvs(i)
                            end do
                        end if
                        ! ============================================================

                        ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
                        !$acc loop seq
                        do i = 1, contxe
                            flux_rs${XYZ}$_vf(-1, k, r, i) = flux_rs${XYZ}$_vf(0, k, r, i) &
                                                             + ds(0)*dalpha_rho_dt(i)
                        end do

                        !$acc loop seq
                        do i = momxb, momxe
                            flux_rs${XYZ}$_vf(-1, k, r, i) = flux_rs${XYZ}$_vf(0, k, r, i) &
                                                             + ds(0)*(vel(i - contxe)*drho_dt &
                                                                      + rho*dvel_dt(i - contxe))
                        end do

                        flux_rs${XYZ}$_vf(-1, k, r, E_idx) = flux_rs${XYZ}$_vf(0, k, r, E_idx) &
                                                             + ds(0)*(pres*dgamma_dt &
                                                                      + gamma*dpres_dt &
                                                                      + dpi_inf_dt &
                                                                      + dqv_dt &
                                                                      + rho*vel_dv_dt_sum &
                                                                      + 5d-1*drho_dt*vel_K_sum)

                        if (riemann_solver == 1) then
                            !$acc loop seq
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(-1, k, r, i) = 0d0
                            end do

                            !$acc loop seq
                            do i = advxb, advxe
                                flux_src_rs${XYZ}$_vf(-1, k, r, i) = &
                                    1d0/max(abs(vel(dir_idx(1))), sgm_eps) &
                                    *sign(1d0, vel(dir_idx(1))) &
                                    *(flux_rs${XYZ}$_vf(0, k, r, i) &
                                      + vel(dir_idx(1)) &
                                      *flux_src_rs${XYZ}$_vf(0, k, r, i) &
                                      + ds(0)*dadv_dt(i - E_idx))
                            end do

                        else

                            !$acc loop seq
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(-1, k, r, i) = flux_rs${XYZ}$_vf(0, k, r, i) + &
                                                                 ds(0)*dadv_dt(i - E_idx)
                            end do

                            !$acc loop seq
                            do i = advxb, advxe
                                flux_src_rs${XYZ}$_vf(-1, k, r, i) = flux_src_rs${XYZ}$_vf(0, k, r, i)
                            end do

                        end if
                        ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============

                    end do
                end do
            end if
        #:endfor

        ! END: FD2 or FD4 of RHS at j = 0 ==================================

        ! The reshaping of outputted data and disssociation of the FD and PI
        ! coefficients, or CBC coefficients, respectively, based on selected
        ! CBC coordinate direction.
        call s_finalize_cbc(flux_vf, flux_src_vf, &
                            ix, iy, iz)
    end subroutine s_cbc

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      selected CBC.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                                ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: flux_vf, flux_src_vf

        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, r !< Generic loop iterators

        ! Configuring the coordinate direction indexes and flags

        ! Determining the indicial shift based on CBC location

        ! END: Allocation/Association of Primitive and Flux Variables ======

        if (cbc_dir == 1) then
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        elseif (cbc_dir == 2) then
            is1%beg = 0; is1%end = buff_size; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0d0, 1d0, 0d0/)
        else
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0d0, 0d0, 1d0/)
        end if

        dj = max(0, cbc_loc)
        !$acc update device(is1, is2, is3, dj)
        !$acc update device( dir_idx, dir_flg)

        ! Reshaping Inputted Data in x-direction ===========================
        if (cbc_dir == 1) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsx_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(dj*(m - 2*j) + j, k, r)
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsx_vf(j, k, r, momxb) = &
                            q_prim_vf(momxb)%sf(dj*(m - 2*j) + j, k, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsx_vf(j, k, r, i) = &
                                flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsx_vf(j, k, r, momxb) = &
                            flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsx_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsx_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Inputted Data in x-direction ======================

            ! Reshaping Inputted Data in y-direction ===========================
        elseif (cbc_dir == 2) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsy_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(k, dj*(n - 2*j) + j, r)
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsy_vf(j, k, r, momxb + 1) = &
                            q_prim_vf(momxb + 1)%sf(k, dj*(n - 2*j) + j, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsy_vf(j, k, r, i) = &
                                flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsy_vf(j, k, r, momxb + 1) = &
                            flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsy_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsy_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Inputted Data in y-direction ======================

            ! Reshaping Inputted Data in z-direction ===========================
        else

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsz_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(r, k, dj*(p - 2*j) + j)
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsz_vf(j, k, r, momxe) = &
                            q_prim_vf(momxe)%sf(r, k, dj*(p - 2*j) + j)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsz_vf(j, k, r, i) = &
                                flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsz_vf(j, k, r, momxe) = &
                            flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsz_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsz_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Inputted Data in z-direction ======================

        ! Association of the procedural pointer to the appropriate procedure
        ! that will be utilized in the evaluation of L variables for the CBC

        ! ==================================================================

        ! ==================================================================

    end subroutine s_initialize_cbc

    !>  Deallocation and/or the disassociation procedures that
        !!      are necessary in order to finalize the CBC application
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_finalize_cbc(flux_vf, flux_src_vf, &
                              ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf

        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, r !< Generic loop iterators

        ! Determining the indicial shift based on CBC location
        dj = max(0, cbc_loc)
        !$acc update device(dj)

        ! Reshaping Outputted Data in x-direction ==========================
        if (cbc_dir == 1) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_rsx_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do
            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_rsx_vf(j, k, r, momxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                    flux_src_rsx_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_src_rsx_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in x-direction =====================

            ! Reshaping Outputted Data in y-direction ==========================
        elseif (cbc_dir == 2) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_rsy_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_rsy_vf(j, k, r, momxb + 1)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                    flux_src_rsy_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_src_rsy_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Outputted Data in y-direction =====================

            ! Reshaping Outputted Data in z-direction ==========================
        else

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_rsz_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_rsz_vf(j, k, r, momxe)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                    flux_src_rsz_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_src_rsz_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Outputted Data in z-direction =====================

        ! Deallocation/Disassociation of Primitive and Flux Variables ======

        ! ==================================================================

        ! Nullifying procedural pointer used in evaluation of L for the CBC

    end subroutine s_finalize_cbc

    ! Detext if the problem has any characteristic boundary conditions
    subroutine s_any_cbc_boundaries(toggle)

        logical, intent(inout) :: toggle

        toggle = .false.

        #:for BC in {-5, -6, -7, -8, -9, -10, -11, -12, -13}
            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == ${BC}$)) then
                toggle = .true.
            end if
        #:endfor

    end subroutine

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_cbc_module

        logical :: is_cbc

        call s_any_cbc_boundaries(is_cbc)

        if (is_cbc .eqv. .false.) return

        ! Deallocating the cell-average primitive variables
        @:DEALLOCATE_GLOBAL(q_prim_rsx_vf)
        if (weno_order > 1) then
            @:DEALLOCATE_GLOBAL(F_rsx_vf, F_src_rsx_vf)
        end if
        @:DEALLOCATE_GLOBAL(flux_rsx_vf, flux_src_rsx_vf)

        if (n > 0) then
            @:DEALLOCATE_GLOBAL(q_prim_rsy_vf)
            if (weno_order > 1) then
                @:DEALLOCATE_GLOBAL(F_rsy_vf, F_src_rsy_vf)
            end if
            @:DEALLOCATE_GLOBAL(flux_rsy_vf, flux_src_rsy_vf)
        end if
        if (p > 0) then
            @:DEALLOCATE_GLOBAL(q_prim_rsz_vf)
            if (weno_order > 1) then
                @:DEALLOCATE_GLOBAL(F_rsz_vf, F_src_rsz_vf)
            end if
            @:DEALLOCATE_GLOBAL(flux_rsz_vf, flux_src_rsz_vf)
        end if

        ! Deallocating the cell-width distribution in the s-direction
        @:DEALLOCATE_GLOBAL(ds)

        ! Deallocating CBC Coefficients in x-direction =====================
        if (any((/bc_x%beg, bc_x%end/) <= -5) .and. any((/bc_x%beg, bc_x%end/) >= -13)) then
            @:DEALLOCATE_GLOBAL(fd_coef_x)
            if (weno_order > 1) then
                @:DEALLOCATE_GLOBAL(pi_coef_x)
            end if
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in y-direction =====================
        if (n > 0 .and. any((/bc_y%beg, bc_y%end/) <= -5) .and. &
            any((/bc_y%beg, bc_y%end/) >= -13 .and. bc_y%beg /= -14)) then
            @:DEALLOCATE_GLOBAL(fd_coef_y)
            if (weno_order > 1) then
                @:DEALLOCATE_GLOBAL(pi_coef_y)
            end if
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in z-direction =====================
        if (p > 0 .and. any((/bc_z%beg, bc_z%end/) <= -5) .and. any((/bc_z%beg, bc_z%end/) >= -13)) then
            @:DEALLOCATE_GLOBAL(fd_coef_z)
            if (weno_order > 1) then
                @:DEALLOCATE_GLOBAL(pi_coef_z)
            end if
        end if
        ! ==================================================================

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables

    end subroutine s_finalize_cbc_module

end module m_cbc
