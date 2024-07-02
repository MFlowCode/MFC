!>
!! @file m_qbmm.f90
!! @brief Contains module m_qbmm

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief This module is used to compute moment inversion via qbmm
module m_qbmm

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic           !< Functions to compare floating point numbers

    use m_helper

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_qbmm_module, s_mom_inv, s_coeff, s_compute_qbmm_rhs

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), momrhs)
    !$acc declare link(momrhs)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: momrhs
    !$acc declare create(momrhs)
#endif
    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: nterms = ${nterms}$
    #:else
        integer :: nterms
        !$acc declare create(nterms)
    #:endif

    type(int_bounds_info) :: is1_qbmm, is2_qbmm, is3_qbmm
!$acc declare create(is1_qbmm, is2_qbmm, is3_qbmm)

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), bubrs)
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:, :), bubmoms)
    !$acc declare link(bubrs, bubmoms)
#else
    integer, allocatable, dimension(:) :: bubrs
    integer, allocatable, dimension(:, :) :: bubmoms
    !$acc declare create(bubrs, bubmoms)
#endif

contains

    subroutine s_initialize_qbmm_module

        integer :: i1, i2, q, i, j

        #:if not MFC_CASE_OPTIMIZATION

            if (bubble_model == 2) then
                ! Keller-Miksis without viscosity/surface tension
                nterms = 32
            else if (bubble_model == 3) then
                ! Rayleigh-Plesset with viscosity/surface tension
                nterms = 7
            end if

            !$acc enter data copyin(nterms)
            !$acc update device(nterms)

        #:endif

        @:ALLOCATE_GLOBAL(momrhs(3, 0:2, 0:2, nterms, nb))
        momrhs = 0d0

        ! Assigns the required RHS moments for moment transport equations
        ! The rhs%(:,3) is only to be used for R0 quadrature, not for computing X/Y indices
        ! Accounts for different governing equations in polytropic and non-polytropic models
        if (.not. polytropic) then
            do q = 1, nb
                do i1 = 0, 2; do i2 = 0, 2
                        if ((i1 + i2) <= 2) then
                            if (bubble_model == 3) then
                                momrhs(1, i1, i2, 1, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 1, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 1, q) = 0d0

                                momrhs(1, i1, i2, 2, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 2, q) = 1.d0 + i2
                                momrhs(3, i1, i2, 2, q) = 0d0

                                momrhs(1, i1, i2, 3, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 3, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 3, q) = 0d0

                                momrhs(1, i1, i2, 4, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 4, q) = 1.d0 + i2
                                momrhs(3, i1, i2, 4, q) = 0d0

                                if (.not. f_is_default(Re_inv)) then
                                    ! add viscosity
                                    momrhs(1, i1, i2, 5, q) = -2.d0 + i1
                                    momrhs(2, i1, i2, 5, q) = i2
                                    momrhs(3, i1, i2, 5, q) = 0d0
                                end if

                                if (.not. f_is_default(Web)) then
                                    ! add surface tension
                                    momrhs(1, i1, i2, 6, q) = -2.d0 + i1
                                    momrhs(2, i1, i2, 6, q) = -1.d0 + i2
                                    momrhs(3, i1, i2, 6, q) = 0d0
                                end if

                                momrhs(1, i1, i2, 7, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 7, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 7, q) = 0d0

                            else if (bubble_model == 2) then
                                ! KM with approximation of 1/(1-V/C) = 1+V/C
                                momrhs(1, i1, i2, 1, q) = -1d0 + i1
                                momrhs(2, i1, i2, 1, q) = 1d0 + i2
                                momrhs(3, i1, i2, 1, q) = 0d0

                                momrhs(1, i1, i2, 2, q) = -1d0 + i1
                                momrhs(2, i1, i2, 2, q) = 2d0 + i2
                                momrhs(3, i1, i2, 2, q) = 0d0

                                momrhs(1, i1, i2, 3, q) = -1d0 + i1
                                momrhs(2, i1, i2, 3, q) = 3d0 + i2
                                momrhs(3, i1, i2, 3, q) = 0d0

                                momrhs(1, i1, i2, 4, q) = -1d0 + i1
                                momrhs(2, i1, i2, 4, q) = -1d0 + i2
                                momrhs(3, i1, i2, 4, q) = 0d0

                                momrhs(1, i1, i2, 5, q) = -1d0 + i1
                                momrhs(2, i1, i2, 5, q) = i2
                                momrhs(3, i1, i2, 5, q) = 0d0

                                momrhs(1, i1, i2, 6, q) = -1d0 + i1
                                momrhs(2, i1, i2, 6, q) = 1d0 + i2
                                momrhs(3, i1, i2, 6, q) = 0d0

                                momrhs(1, i1, i2, 7, q) = -1d0 + i1
                                momrhs(2, i1, i2, 7, q) = -1d0 + i2
                                momrhs(3, i1, i2, 7, q) = 0d0

                                momrhs(1, i1, i2, 8, q) = -1d0 + i1
                                momrhs(2, i1, i2, 8, q) = i2
                                momrhs(3, i1, i2, 8, q) = 0d0

                                momrhs(1, i1, i2, 9, q) = -1d0 + i1
                                momrhs(2, i1, i2, 9, q) = 1d0 + i2
                                momrhs(3, i1, i2, 9, q) = 0d0

                                momrhs(1, i1, i2, 10, q) = -1d0 + i1
                                momrhs(2, i1, i2, 10, q) = i2
                                momrhs(3, i1, i2, 10, q) = 0d0

                                momrhs(1, i1, i2, 11, q) = -1d0 + i1
                                momrhs(2, i1, i2, 11, q) = 1d0 + i2
                                momrhs(3, i1, i2, 11, q) = 0d0

                                momrhs(1, i1, i2, 12, q) = -1d0 + i1
                                momrhs(2, i1, i2, 12, q) = 1d0 + i2
                                momrhs(3, i1, i2, 12, q) = 0d0

                                momrhs(1, i1, i2, 13, q) = -1d0 + i1
                                momrhs(2, i1, i2, 13, q) = -1d0 + i2
                                momrhs(3, i1, i2, 13, q) = 0d0

                                momrhs(1, i1, i2, 14, q) = -1d0 + i1
                                momrhs(2, i1, i2, 14, q) = i2
                                momrhs(3, i1, i2, 14, q) = 0d0

                                momrhs(1, i1, i2, 15, q) = -1d0 + i1
                                momrhs(2, i1, i2, 15, q) = 1d0 + i2
                                momrhs(3, i1, i2, 15, q) = 0d0

                                momrhs(1, i1, i2, 16, q) = -2d0 + i1
                                momrhs(2, i1, i2, 16, q) = i2
                                momrhs(3, i1, i2, 16, q) = 0d0

                                momrhs(1, i1, i2, 17, q) = -2d0 + i1
                                momrhs(2, i1, i2, 17, q) = -1d0 + i2
                                momrhs(3, i1, i2, 17, q) = 0d0

                                momrhs(1, i1, i2, 18, q) = -2d0 + i1
                                momrhs(2, i1, i2, 18, q) = 1d0 + i2
                                momrhs(3, i1, i2, 18, q) = 0d0

                                momrhs(1, i1, i2, 19, q) = -2d0 + i1
                                momrhs(2, i1, i2, 19, q) = 2d0 + i2
                                momrhs(3, i1, i2, 19, q) = 0d0

                                momrhs(1, i1, i2, 20, q) = -2d0 + i1
                                momrhs(2, i1, i2, 20, q) = -1d0 + i2
                                momrhs(3, i1, i2, 20, q) = 0d0

                                momrhs(1, i1, i2, 21, q) = -2d0 + i1
                                momrhs(2, i1, i2, 21, q) = i2
                                momrhs(3, i1, i2, 21, q) = 0d0

                                momrhs(1, i1, i2, 22, q) = -2d0 + i1
                                momrhs(2, i1, i2, 22, q) = -1d0 + i2
                                momrhs(3, i1, i2, 22, q) = 0d0

                                momrhs(1, i1, i2, 23, q) = -2d0 + i1
                                momrhs(2, i1, i2, 23, q) = i2
                                momrhs(3, i1, i2, 23, q) = 0d0

                                momrhs(1, i1, i2, 24, q) = -3d0 + i1
                                momrhs(2, i1, i2, 24, q) = i2
                                momrhs(3, i1, i2, 24, q) = 0d0

                                momrhs(1, i1, i2, 25, q) = -3d0 + i1
                                momrhs(2, i1, i2, 25, q) = -1d0 + i2
                                momrhs(3, i1, i2, 25, q) = 0d0

                                momrhs(1, i1, i2, 26, q) = -2d0 + i1
                                momrhs(2, i1, i2, 26, q) = i2
                                momrhs(3, i1, i2, 26, q) = 0d0

                                momrhs(1, i1, i2, 27, q) = -1d0 + i1
                                momrhs(2, i1, i2, 27, q) = -1d0 + i2
                                momrhs(3, i1, i2, 27, q) = 0d0

                                momrhs(1, i1, i2, 28, q) = -1d0 + i1
                                momrhs(2, i1, i2, 28, q) = i2
                                momrhs(3, i1, i2, 28, q) = 0d0

                                momrhs(1, i1, i2, 29, q) = -2d0 + i1
                                momrhs(2, i1, i2, 29, q) = i2
                                momrhs(3, i1, i2, 29, q) = 0d0

                                momrhs(1, i1, i2, 30, q) = -1d0 + i1
                                momrhs(2, i1, i2, 30, q) = -1d0 + i2
                                momrhs(3, i1, i2, 30, q) = 0d0

                                momrhs(1, i1, i2, 31, q) = -1d0 + i1
                                momrhs(2, i1, i2, 31, q) = i2
                                momrhs(3, i1, i2, 31, q) = 0d0

                                momrhs(1, i1, i2, 32, q) = -2d0 + i1
                                momrhs(2, i1, i2, 32, q) = i2
                                momrhs(3, i1, i2, 32, q) = 0d0
                            end if
                        end if
                    end do; end do
            end do

        else
            do q = 1, nb
                do i1 = 0, 2; do i2 = 0, 2
                        if ((i1 + i2) <= 2) then
                            if (bubble_model == 3) then
                                momrhs(1, i1, i2, 1, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 1, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 1, q) = 0d0

                                momrhs(1, i1, i2, 2, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 2, q) = 1.d0 + i2
                                momrhs(3, i1, i2, 2, q) = 0d0

                                momrhs(1, i1, i2, 3, q) = -1.d0 + i1 - 3.d0*gam
                                momrhs(2, i1, i2, 3, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 3, q) = 3.d0*gam

                                momrhs(1, i1, i2, 4, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 4, q) = 1.d0 + i2
                                momrhs(3, i1, i2, 4, q) = 0d0

                                if (.not. f_is_default(Re_inv)) then
                                    ! add viscosity
                                    momrhs(1, i1, i2, 5, q) = -2.d0 + i1
                                    momrhs(2, i1, i2, 5, q) = i2
                                    momrhs(3, i1, i2, 5, q) = 0d0
                                end if

                                if (.not. f_is_default(Web)) then
                                    ! add surface tension
                                    momrhs(1, i1, i2, 6, q) = -2.d0 + i1
                                    momrhs(2, i1, i2, 6, q) = -1.d0 + i2
                                    momrhs(3, i1, i2, 6, q) = 0d0
                                end if

                                momrhs(1, i1, i2, 7, q) = -1.d0 + i1
                                momrhs(2, i1, i2, 7, q) = -1.d0 + i2
                                momrhs(3, i1, i2, 7, q) = 0d0

                            else if (bubble_model == 2) then
                                ! KM with approximation of 1/(1-V/C) = 1+V/C
                                momrhs(1, i1, i2, 1, q) = -1d0 + i1
                                momrhs(2, i1, i2, 1, q) = 1d0 + i2
                                momrhs(3, i1, i2, 1, q) = 0d0

                                momrhs(1, i1, i2, 2, q) = -1d0 + i1
                                momrhs(2, i1, i2, 2, q) = 2d0 + i2
                                momrhs(3, i1, i2, 2, q) = 0d0

                                momrhs(1, i1, i2, 3, q) = -1d0 + i1
                                momrhs(2, i1, i2, 3, q) = 3d0 + i2
                                momrhs(3, i1, i2, 3, q) = 0d0

                                momrhs(1, i1, i2, 4, q) = -1d0 + i1
                                momrhs(2, i1, i2, 4, q) = -1d0 + i2
                                momrhs(3, i1, i2, 4, q) = 0d0

                                momrhs(1, i1, i2, 5, q) = -1d0 + i1
                                momrhs(2, i1, i2, 5, q) = i2
                                momrhs(3, i1, i2, 5, q) = 0d0

                                momrhs(1, i1, i2, 6, q) = -1d0 + i1
                                momrhs(2, i1, i2, 6, q) = 1d0 + i2
                                momrhs(3, i1, i2, 6, q) = 0d0

                                momrhs(1, i1, i2, 7, q) = -1d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 7, q) = -1d0 + i2
                                momrhs(3, i1, i2, 7, q) = 3d0*gam

                                momrhs(1, i1, i2, 8, q) = -1d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 8, q) = i2
                                momrhs(3, i1, i2, 8, q) = 3d0*gam

                                momrhs(1, i1, i2, 9, q) = -1d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 9, q) = 1d0 + i2
                                momrhs(3, i1, i2, 9, q) = 3d0*gam

                                momrhs(1, i1, i2, 10, q) = -1d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 10, q) = i2
                                momrhs(3, i1, i2, 10, q) = 3d0*gam

                                momrhs(1, i1, i2, 11, q) = -1d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 11, q) = 1d0 + i2
                                momrhs(3, i1, i2, 11, q) = 3d0*gam

                                momrhs(1, i1, i2, 12, q) = -1d0 + i1
                                momrhs(2, i1, i2, 12, q) = 1d0 + i2
                                momrhs(3, i1, i2, 12, q) = 0d0

                                momrhs(1, i1, i2, 13, q) = -1d0 + i1
                                momrhs(2, i1, i2, 13, q) = -1d0 + i2
                                momrhs(3, i1, i2, 13, q) = 0d0

                                momrhs(1, i1, i2, 14, q) = -1d0 + i1
                                momrhs(2, i1, i2, 14, q) = i2
                                momrhs(3, i1, i2, 14, q) = 0d0

                                momrhs(1, i1, i2, 15, q) = -1d0 + i1
                                momrhs(2, i1, i2, 15, q) = 1d0 + i2
                                momrhs(3, i1, i2, 15, q) = 0d0

                                momrhs(1, i1, i2, 16, q) = -2d0 + i1
                                momrhs(2, i1, i2, 16, q) = i2
                                momrhs(3, i1, i2, 16, q) = 0d0

                                momrhs(1, i1, i2, 17, q) = -2d0 + i1
                                momrhs(2, i1, i2, 17, q) = -1d0 + i2
                                momrhs(3, i1, i2, 17, q) = 0d0

                                momrhs(1, i1, i2, 18, q) = -2d0 + i1
                                momrhs(2, i1, i2, 18, q) = 1d0 + i2
                                momrhs(3, i1, i2, 18, q) = 0d0

                                momrhs(1, i1, i2, 19, q) = -2d0 + i1
                                momrhs(2, i1, i2, 19, q) = 2d0 + i2
                                momrhs(3, i1, i2, 19, q) = 0d0

                                momrhs(1, i1, i2, 20, q) = -2d0 + i1
                                momrhs(2, i1, i2, 20, q) = -1d0 + i2
                                momrhs(3, i1, i2, 20, q) = 0d0

                                momrhs(1, i1, i2, 21, q) = -2d0 + i1
                                momrhs(2, i1, i2, 21, q) = i2
                                momrhs(3, i1, i2, 21, q) = 0d0

                                momrhs(1, i1, i2, 22, q) = -2d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 22, q) = -1d0 + i2
                                momrhs(3, i1, i2, 22, q) = 3d0*gam

                                momrhs(1, i1, i2, 23, q) = -2d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 23, q) = i2
                                momrhs(3, i1, i2, 23, q) = 3d0*gam

                                momrhs(1, i1, i2, 24, q) = -3d0 + i1
                                momrhs(2, i1, i2, 24, q) = i2
                                momrhs(3, i1, i2, 24, q) = 0d0

                                momrhs(1, i1, i2, 25, q) = -3d0 + i1
                                momrhs(2, i1, i2, 25, q) = -1d0 + i2
                                momrhs(3, i1, i2, 25, q) = 0d0

                                momrhs(1, i1, i2, 26, q) = -2d0 + i1 - 3d0*gam
                                momrhs(2, i1, i2, 26, q) = i2
                                momrhs(3, i1, i2, 26, q) = 3d0*gam

                            end if
                        end if
                    end do; end do
            end do
        end if

        !$acc update device(momrhs)

        @:ALLOCATE_GLOBAL(bubrs(1:nb))
        @:ALLOCATE_GLOBAL(bubmoms(1:nb, 1:nmom))

        do i = 1, nb
            bubrs(i) = bub_idx%rs(i)
        end do
        !$acc update device(bubrs)

        do j = 1, nmom
            do i = 1, nb
                bubmoms(i, j) = bub_idx%moms(i, j)
            end do
        end do
        !$acc update device(bubmoms)

    end subroutine s_initialize_qbmm_module

    subroutine s_compute_qbmm_rhs(idir, q_cons_vf, q_prim_vf, rhs_vf, flux_n_vf, pb, rhs_pb, mv, rhs_mv)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in) :: flux_n_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, rhs_pb
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: mv, rhs_mv

        integer :: i, j, k, l, q

        real(kind(0d0)) :: nb_q, nb_dot, R, R2, nR, nR2, nR_dot, nR2_dot, var, AX

        if (idir == 1) then

            !Non-polytropic qbmm needs to account for change in bubble radius due to a change in nb
            if (.not. polytropic) then
                !$acc parallel loop collapse(5) gang vector default(present) private(nb_q, nR, nR2, R, R2, nb_dot, nR_dot, nR2_dot, var, AX)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    nb_q = q_cons_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR = q_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2 = q_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    R = q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    R2 = q_prim_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    if (R2 - R**2d0 > 0d0) then
                                        var = R2 - R**2d0
                                    else
                                        var = verysmall
                                    end if

                                    if (q <= 2) then
                                        AX = R - dsqrt(var)
                                    else
                                        AX = R + dsqrt(var)
                                    end if

                                    nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dx(j)*AX*nb_q**2)* &
                                                            (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))

                                    if (q <= 2) then
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dx(j)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dx(j)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))

                                    else
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dx(j)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dx(j)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    end if

                                end do
                            end do
                        end do
                    end do
                end do
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do q = 0, n
                    do i = 0, m

                        rhs_vf(alf_idx)%sf(i, q, l) = rhs_vf(alf_idx)%sf(i, q, l) + mom_sp(2)%sf(i, q, l)
                        j = bubxb
                        !$acc loop seq
                        do k = 1, nb
                            rhs_vf(j)%sf(i, q, l) = &
                                rhs_vf(j)%sf(i, q, l) + mom_3d(0, 0, k)%sf(i, q, l)
                            rhs_vf(j + 1)%sf(i, q, l) = &
                                rhs_vf(j + 1)%sf(i, q, l) + mom_3d(1, 0, k)%sf(i, q, l)
                            rhs_vf(j + 2)%sf(i, q, l) = &
                                rhs_vf(j + 2)%sf(i, q, l) + mom_3d(0, 1, k)%sf(i, q, l)
                            rhs_vf(j + 3)%sf(i, q, l) = &
                                rhs_vf(j + 3)%sf(i, q, l) + mom_3d(2, 0, k)%sf(i, q, l)
                            rhs_vf(j + 4)%sf(i, q, l) = &
                                rhs_vf(j + 4)%sf(i, q, l) + mom_3d(1, 1, k)%sf(i, q, l)
                            rhs_vf(j + 5)%sf(i, q, l) = &
                                rhs_vf(j + 5)%sf(i, q, l) + mom_3d(0, 2, k)%sf(i, q, l)
                            j = j + 6
                        end do

                    end do
                end do
            end do

        elseif (idir == 2) then

            !Non-polytropic qbmm needs to account for change in bubble radius due to a change in nb
            if (.not. polytropic) then
                !$acc parallel loop collapse(5) gang vector default(present) private(nb_q, nR, nR2, R, R2, nb_dot, nR_dot, nR2_dot, var, AX)
                do i = 1, nb
                    do q = 1, nnode
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    nb_q = q_cons_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR = q_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2 = q_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    R = q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    R2 = q_prim_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    if (R2 - R**2d0 > 0d0) then
                                        var = R2 - R**2d0
                                    else
                                        var = verysmall
                                    end if

                                    if (q <= 2) then
                                        AX = R - dsqrt(var)
                                    else
                                        AX = R + dsqrt(var)
                                    end if

                                    nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                    rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dy(k)*AX*nb_q**2)* &
                                                            (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))

                                    if (q <= 2) then
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dy(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dy(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))

                                    else
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dy(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dy(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    end if

                                end do
                            end do
                        end do
                    end do
                end do
            end if

        elseif (idir == 3) then

            if (.not. polytropic) then
                if (grid_geometry == 3) then
                    !Non-polytropic qbmm needs to account for change in bubble radius due to a change in nb
                    !$acc parallel loop collapse(5) gang vector default(present) private(nb_q, nR, nR2, R, R2, nb_dot, nR_dot, nR2_dot, var, AX)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 0, m
                                        nb_q = q_cons_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                        nR = q_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        nR2 = q_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                        R = q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        R2 = q_prim_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                        if (R2 - R**2d0 > 0d0) then
                                            var = R2 - R**2d0
                                        else
                                            var = verysmall
                                        end if

                                        if (q <= 2) then
                                            AX = R - dsqrt(var)
                                        else
                                            AX = R + dsqrt(var)
                                        end if

                                        nb_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l))
                                        nR_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l))
                                        nR2_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l))

                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*y_cc(k)*AX*nb_q**2)* &
                                                                (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))

                                        if (q <= 2) then
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dz(l)*y_cc(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dz(l)*y_cc(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))

                                        else
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*y_cc(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*y_cc(k)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                else
                    !Non-polytropic qbmm needs to account for change in bubble radius due to a change in nb
                    !$acc parallel loop collapse(5) gang vector default(present) private(nb_q, nR, nR2, R, R2, nb_dot, nR_dot, nR2_dot, var, AX)
                    do i = 1, nb
                        do q = 1, nnode
                            do l = 0, p
                                do k = 0, n
                                    do j = 0, m
                                        nb_q = q_cons_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                        nR = q_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        nR2 = q_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                        R = q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        R2 = q_prim_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                        if (R2 - R**2d0 > 0d0) then
                                            var = R2 - R**2d0
                                        else
                                            var = verysmall
                                        end if

                                        if (q <= 2) then
                                            AX = R - dsqrt(var)
                                        else
                                            AX = R + dsqrt(var)
                                        end if

                                        nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                        nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)

                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*AX*nb_q**2)* &
                                                                (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))

                                        if (q <= 2) then
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dz(l)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3d0*gam/(dz(l)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))

                                        else
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3d0*gam/(dz(l)*AX*nb_q**2*dsqrt(var)*2d0)* &
                                                                    (-2d0*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        end if

                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end if

    end subroutine

!Coefficient array for non-polytropic model (pb and mv values are accounted in wght_pb and wght_mv)

    subroutine s_coeff_nonpoly(pres, rho, c, coeffs)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_coeff_nonpoly
#else
        !$acc routine seq
#endif
        real(kind(0.d0)), intent(in) :: pres, rho, c
        real(kind(0.d0)), dimension(nterms, 0:2, 0:2), intent(out) :: coeffs

        integer :: i1, i2, q

        coeffs = 0d0

        do i2 = 0, 2; do i1 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1d0*i2*pres/rho
                        coeffs(2, i1, i2) = -3d0*i2/2d0
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (.not. f_is_default(Re_inv)) coeffs(5, i1, i2) = -4d0*i2*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(6, i1, i2) = -2d0*i2/Web/rho
                        coeffs(7, i1, i2) = 0d0
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3d0*i2/2d0
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2d0*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2d0*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2d0*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3d0*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3d0*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                        coeffs(13, i1, i2) = 0d0
                        coeffs(14, i1, i2) = 0d0
                        coeffs(15, i1, i2) = 0d0
                        if (.not. f_is_default(Re_inv)) coeffs(16, i1, i2) = -i2*4d0*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(17, i1, i2) = -i2*2d0/Web/rho
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(18, i1, i2) = i2*6d0*Re_inv/(rho*c)
                            coeffs(19, i1, i2) = -i2*2d0*Re_inv/(rho*c*c)
                            coeffs(20, i1, i2) = i2*4d0*pres*Re_inv/(rho*rho*c)
                            coeffs(21, i1, i2) = i2*4d0*pres*Re_inv/(rho*rho*c*c)
                            coeffs(22, i1, i2) = -i2*4d0/(rho*rho*c)
                            coeffs(23, i1, i2) = -i2*4d0/(rho*rho*c*c)
                            coeffs(24, i1, i2) = i2*16d0*Re_inv*Re_inv/(rho*rho*c)
                            if (.not. f_is_default(Web)) then
                                coeffs(25, i1, i2) = i2*8d0*Re_inv/Web/(rho*rho*c)
                            end if
                            coeffs(26, i1, i2) = -12d0*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                        coeffs(27, i1, i2) = 3d0*i2*gam*R_v*Tw/(c*rho)
                        coeffs(28, i1, i2) = 3d0*i2*gam*R_v*Tw/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(29, i1, i2) = 12d0*i2*gam*R_v*Tw*Re_inv/(rho*rho*c*c)
                        end if
                        coeffs(30, i1, i2) = 3d0*i2*gam/(c*rho)
                        coeffs(31, i1, i2) = 3d0*i2*gam/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(32, i1, i2) = 12d0*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                    end if
                end if
            end do; end do

    end subroutine s_coeff_nonpoly

!Coefficient array for polytropic model (pb for each R0 bin accounted for in wght_pb)
    subroutine s_coeff(pres, rho, c, coeffs)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_coeff
#else
        !$acc routine seq
#endif

        real(kind(0.d0)), intent(inout) :: pres, rho, c
        real(kind(0.d0)), dimension(nterms, 0:2, 0:2), intent(out) :: coeffs

        integer :: i1, i2, q

        coeffs = 0d0

        do i2 = 0, 2; do i1 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1d0*i2*pres/rho
                        coeffs(2, i1, i2) = -3d0*i2/2d0
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (.not. f_is_default(Re_inv)) coeffs(5, i1, i2) = -4d0*i2*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(6, i1, i2) = -2d0*i2/Web/rho
                        coeffs(7, i1, i2) = i2*pv/rho
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3d0*i2/2d0
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2d0*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2d0*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2d0*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3d0*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3d0*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                        coeffs(13, i1, i2) = i2*(pv)/rho
                        coeffs(14, i1, i2) = 2d0*i2*(pv)/(c*rho)
                        coeffs(15, i1, i2) = i2*(pv)/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) coeffs(16, i1, i2) = -i2*4d0*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(17, i1, i2) = -i2*2d0/Web/rho
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(18, i1, i2) = i2*6d0*Re_inv/(rho*c)
                            coeffs(19, i1, i2) = -i2*2d0*Re_inv/(rho*c*c)
                            coeffs(20, i1, i2) = i2*4d0*pres*Re_inv/(rho*rho*c)
                            coeffs(21, i1, i2) = i2*4d0*pres*Re_inv/(rho*rho*c*c)
                            coeffs(22, i1, i2) = -i2*4d0/(rho*rho*c)
                            coeffs(23, i1, i2) = -i2*4d0/(rho*rho*c*c)
                            coeffs(24, i1, i2) = i2*16d0*Re_inv*Re_inv/(rho*rho*c)
                            if (.not. f_is_default(Web)) then
                                coeffs(25, i1, i2) = i2*8d0*Re_inv/Web/(rho*rho*c)
                            end if
                            coeffs(26, i1, i2) = -12d0*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                    end if
                end if
            end do; end do

    end subroutine s_coeff

    subroutine s_mom_inv(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix, iy, iz, nbub_sc)

        type(scalar_field), dimension(:), intent(inout) :: q_cons_vf, q_prim_vf
        type(scalar_field), dimension(:), intent(inout) :: momsp
        type(scalar_field), dimension(0:, 0:, :), intent(inout) :: moms3d
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, rhs_pb
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: mv, rhs_mv
        type(int_bounds_info), intent(in) :: ix, iy, iz
        real(kind(0d0)), dimension(startx:, starty:, startz:) :: nbub_sc !> Unused Variable not sure what to put as intent

        real(kind(0d0)), dimension(nmom) :: moms, msum
        real(kind(0d0)), dimension(nnode, nb) :: wght, abscX, abscY, wght_pb, wght_mv, wght_ht, ht
        real(kind(0d0)), dimension(nterms, 0:2, 0:2) :: mom3d_terms, coeff
        real(kind(0d0)) :: pres, rho, nbub, c, alf, R3, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, T_bar, grad_T
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: n_tait, B_tait

        integer :: j, k, l, q, r, s !< Loop variables
        integer :: id1, id2, id3
        integer :: i1, i2

        is1_qbmm = ix; is2_qbmm = iy; is3_qbmm = iz

        !$acc update device(is1_qbmm, is2_qbmm, is3_qbmm)

        !$acc parallel loop collapse(3) gang vector default(present) private(moms, msum, wght, abscX, abscY, wght_pb, wght_mv, wght_ht, coeff, ht, r, q, n_tait, B_tait, pres, rho, nbub, c, alf, R3, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, T_bar, grad_T)
        do id3 = is3_qbmm%beg, is3_qbmm%end
            do id2 = is2_qbmm%beg, is2_qbmm%end
                do id1 = is1_qbmm%beg, is1_qbmm%end

                    alf = q_prim_vf(alf_idx)%sf(id1, id2, id3)
                    pres = q_prim_vf(E_idx)%sf(id1, id2, id3)
                    rho = q_prim_vf(contxb)%sf(id1, id2, id3)
                    if (bubble_model == 2) then
                        n_tait = gammas(1)
                        n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'
                        B_tait = pi_infs(1)*(n_tait - 1)/n_tait
                        c = n_tait*(pres + B_tait)*(1d0 - alf)/(rho)

                        if (c > 0.d0) then
                            c = DSQRT(c)
                        else
                            c = sgm_eps
                        end if
                    end if

                    if (polytropic) then
                        call s_coeff(pres, rho, c, coeff)
                    else
                        call s_coeff_nonpoly(pres, rho, c, coeff)
                    end if

                    ! SHB: Manually adjusted pressure here for no-coupling case
                    ! pres = 1d0/0.3d0
                    if (alf > small_alf) then

                        nbub = q_cons_vf(bubxb)%sf(id1, id2, id3)

                        !$acc loop seq
                        do q = 1, nb
                            !Initialize moment set for each R0 bin
                            !$acc loop seq
                            do r = 2, nmom
                                moms(r) = q_prim_vf(bubmoms(q, r))%sf(id1, id2, id3)
                            end do

                            moms(1) = 1d0

                            call s_chyqmom(moms, wght(:, q), abscX(:, q), abscY(:, q))

                            if (polytropic) then
                                !Account for bubble pressure pb0 at each R0 bin
                                !$acc loop seq
                                do j = 1, nnode
                                    wght_pb(j, q) = wght(j, q)*(pb0(q) - pv)
                                end do
                            else
                                !Account for bubble pressure, mass transfer rate and heat transfer rate in wght_pb, wght_mv and wght_ht using Preston model
                                !$acc loop seq
                                do j = 1, nnode
                                    chi_vw = 1.d0/(1.d0 + R_v/R_n*(pb(id1, id2, id3, j, q)/pv - 1.d0))
                                    x_vw = M_n*chi_vw/(M_v + (M_n - M_v)*chi_vw)
                                    k_mw = x_vw*k_v(q)/(x_vw + (1.d0 - x_vw)*phi_vn) &
                                           + (1.d0 - x_vw)*k_n(q)/(x_vw*phi_nv + 1.d0 - x_vw)
                                    rho_mw = pv/(chi_vw*R_v*Tw)
                                    rhs_mv(id1, id2, id3, j, q) = -Re_trans_c(q)*((mv(id1, id2, id3, j, q)/(mv(id1, id2, id3, j, q) + mass_n0(q))) - chi_vw)
                                    rhs_mv(id1, id2, id3, j, q) = rho_mw*rhs_mv(id1, id2, id3, j, q)/Pe_c/(1.d0 - chi_vw)/abscX(j, q)

                                    T_bar = Tw*(pb(id1, id2, id3, j, q)/pb0(q))*(abscX(j, q)/R0(q))**3 &
                                            *(mass_n0(q) + mass_v0(q))/(mass_n0(q) + mv(id1, id2, id3, j, q))
                                    grad_T = -Re_trans_T(q)*(T_bar - Tw)
                                    ht(j, q) = pb0(q)*k_mw*grad_T/Pe_T(q)/abscX(j, q)

                                    wght_pb(j, q) = wght(j, q)*(pb(id1, id2, id3, j, q))
                                    wght_mv(j, q) = wght(j, q)*(rhs_mv(id1, id2, id3, j, q))
                                    wght_ht(j, q) = wght(j, q)*ht(j, q)
                                end do
                            end if

                            !Compute change in moments due to bubble dynamics
                            r = 1
                            !$acc loop seq
                            do i2 = 0, 2
                                !$acc loop seq
                                do i1 = 0, 2
                                    if ((i1 + i2) <= 2) then
                                        momsum = 0d0
                                        !$acc loop seq
                                        do j = 1, nterms
                                            ! Account for term with pb in Rayleigh Plesset equation
                                            if (bubble_model == 3 .and. j == 3) then
                                                momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                         *f_quad2D(abscX(:, q), abscY(:, q), wght_pb(:, q), momrhs(:, i1, i2, j, q))
                                                ! Account for terms with pb in Keller-Miksis equation
                                            else if (bubble_model == 2 .and. ((j >= 7 .and. j <= 9) .or. (j >= 22 .and. j <= 23) .or. (j >= 10 .and. j <= 11) .or. (j == 26))) then
                                                momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                         *f_quad2D(abscX(:, q), abscY(:, q), wght_pb(:, q), momrhs(:, i1, i2, j, q))
                                                ! Account for terms with mass transfer rate in Keller-Miksis equation
                                            else if (bubble_model == 2 .and. (j >= 27 .and. j <= 29) .and. (.not. polytropic)) then
                                                momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                         *f_quad2D(abscX(:, q), abscY(:, q), wght_mv(:, q), momrhs(:, i1, i2, j, q))
                                                ! Account for terms with heat transfer rate in Keller-Miksis equation
                                            else if (bubble_model == 2 .and. (j >= 30 .and. j <= 32) .and. (.not. polytropic)) then
                                                momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                         *f_quad2D(abscX(:, q), abscY(:, q), wght_ht(:, q), momrhs(:, i1, i2, j, q))
                                            else
                                                momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q)) &
                                                         *f_quad2D(abscX(:, q), abscY(:, q), wght(:, q), momrhs(:, i1, i2, j, q))
                                            end if

                                        end do

                                        moms3d(i1, i2, q)%sf(id1, id2, id3) = nbub*momsum
                                        msum(r) = momsum
                                        r = r + 1

                                    end if
                                end do
                            end do

                            ! Compute change in pb and mv for non-polytroic model
                            if (.not. polytropic) then
                                !$acc loop seq
                                do j = 1, nnode
                                    ! Compute Rdot (drdt) at quadrature node in the ODE for pb (note this is not the same as bubble variable Rdot)
                                    drdt = msum(2)
                                    if (moms(4) - moms(2)**2d0 > 0d0) then
                                        if (j == 1 .or. j == 2) then
                                            drdt2 = -1d0/(2d0*dsqrt(moms(4) - moms(2)**2d0))
                                        else
                                            drdt2 = 1d0/(2d0*dsqrt(moms(4) - moms(2)**2d0))
                                        end if
                                    else
                                        ! Edge case where variance < 0
                                        if (j == 1 .or. j == 2) then
                                            drdt2 = -1d0/(2d0*dsqrt(verysmall))
                                        else
                                            drdt2 = 1d0/(2d0*dsqrt(verysmall))
                                        end if
                                    end if

                                    drdt2 = drdt2*(msum(3) - 2d0*moms(2)*msum(2))
                                    drdt = drdt + drdt2

                                    rhs_pb(id1, id2, id3, j, q) = (-3d0*gam*drdt/abscX(j, q))*(pb(id1, id2, id3, j, q))
                                    rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, q) + (3d0*gam/abscX(j, q))*rhs_mv(id1, id2, id3, j, q)*R_v*Tw
                                    rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, q) + (3d0*gam/abscX(j, q))*ht(j, q)
                                    rhs_mv(id1, id2, id3, j, q) = rhs_mv(id1, id2, id3, j, q)*(4d0*pi*abscX(j, q)**2d0)
                                end do

                            end if
                        end do

                        ! Compute special high-order moments
                        momsp(1)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3d0, 0d0, 0d0)
                        momsp(2)%sf(id1, id2, id3) = 4.d0*pi*nbub*f_quad(abscX, abscY, wght, 2d0, 1d0, 0d0)
                        momsp(3)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3d0, 2d0, 0d0)
                        if (abs(gam - 1.d0) <= 1.d-4) then
                            ! Gam \approx 1, don't risk imaginary quadrature
                            momsp(4)%sf(id1, id2, id3) = 1.d0
                        else
                            !Special moment with bubble pressure pb
                            if (polytropic) then
                                momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght_pb, 3d0*(1d0 - gam), 0d0, 3d0*gam) + pv*f_quad(abscX, abscY, wght, 3d0, 0d0, 0d0) &
                                                             - 4d0*Re_inv*f_quad(abscX, abscY, wght, 2d0, 1d0, 0d0) - (2d0/Web)*f_quad(abscX, abscY, wght, 2d0, 0d0, 0d0)
                            else
                                momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght_pb, 3d0, 0d0, 0d0) &
                                                             - 4d0*Re_inv*f_quad(abscX, abscY, wght, 2d0, 1d0, 0d0) - (2d0/Web)*f_quad(abscX, abscY, wght, 2d0, 0d0, 0d0)
                            end if
                        end if

                    else
                        !$acc loop seq
                        do q = 1, nb
                            !$acc loop seq
                            do i1 = 0, 2
                                !$acc loop seq
                                do i2 = 0, 2
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = 0d0
                                end do
                            end do
                        end do

                        momsp(1)%sf(id1, id2, id3) = 0d0
                        momsp(2)%sf(id1, id2, id3) = 0d0
                        momsp(3)%sf(id1, id2, id3) = 0d0
                        momsp(4)%sf(id1, id2, id3) = 0d0

                    end if

                end do
            end do
        end do

    end subroutine s_mom_inv

    subroutine s_chyqmom(momin, wght, abscX, abscY)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_chyqmom
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(nmom), intent(in) :: momin
        real(kind(0d0)), dimension(nnode), intent(inout) :: wght, abscX, abscY

        real(kind(0d0)), dimension(0:2, 0:2) :: moms
        real(kind(0d0)), dimension(3) :: M1, M3
        real(kind(0d0)), dimension(2) :: myrho, myrho3, up, up3, Vf
        real(kind(0d0)) :: bu, bv, d20, d11, d02, c20, c11, c02
        real(kind(0d0)) :: mu2avg, mu2, vp21, vp22, rho21, rho22

        moms(0, 0) = momin(1)
        moms(1, 0) = momin(2)
        moms(0, 1) = momin(3)
        moms(2, 0) = momin(4)
        moms(1, 1) = momin(5)
        moms(0, 2) = momin(6)

        bu = moms(1, 0)/moms(0, 0)
        bv = moms(0, 1)/moms(0, 0)
        d20 = moms(2, 0)/moms(0, 0)
        d11 = moms(1, 1)/moms(0, 0)
        d02 = moms(0, 2)/moms(0, 0)

        c20 = d20 - bu**2d0; 
        c11 = d11 - bu*bv; 
        c02 = d02 - bv**2d0; 
        M1 = (/1d0, 0d0, c20/)
        call s_hyqmom(myrho, up, M1)
        Vf = c11*up/c20

        mu2avg = c02 - sum(myrho(:)*(Vf(:)**2d0))

        mu2avg = maxval((/mu2avg, 0d0/))
        mu2 = mu2avg
        M3 = (/1d0, 0d0, mu2/)
        call s_hyqmom(myrho3, up3, M3)

        vp21 = up3(1)
        vp22 = up3(2)
        rho21 = myrho3(1)
        rho22 = myrho3(2)

        wght(1) = myrho(1)*rho21
        wght(2) = myrho(1)*rho22
        wght(3) = myrho(2)*rho21
        wght(4) = myrho(2)*rho22
        wght = moms(0, 0)*wght

        abscX(1) = up(1)
        abscX(2) = up(1)
        abscX(3) = up(2)
        abscX(4) = up(2)
        abscX = bu + abscX

        abscY(1) = Vf(1) + vp21
        abscY(2) = Vf(1) + vp22
        abscY(3) = Vf(2) + vp21
        abscY(4) = Vf(2) + vp22
        abscY = bv + abscY

    end subroutine s_chyqmom

    subroutine s_hyqmom(frho, fup, fmom)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_hyqmom
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(2), intent(inout) :: frho, fup
        real(kind(0d0)), dimension(3), intent(in) :: fmom

        real(kind(0d0)) :: bu, d2, c2

        bu = fmom(2)/fmom(1)
        d2 = fmom(3)/fmom(1)
        c2 = d2 - bu**2d0
        frho(1) = fmom(1)/2d0; 
        frho(2) = fmom(1)/2d0; 
        c2 = maxval((/c2, verysmall/))
        fup(1) = bu - DSQRT(c2)
        fup(2) = bu + DSQRT(c2)

    end subroutine s_hyqmom

    function f_quad(abscX, abscY, wght_in, q, r, s)
        !$acc routine seq
        real(kind(0.d0)), dimension(nnode, nb), intent(in) :: abscX, abscY, wght_in
        real(kind(0.d0)), intent(in) :: q, r, s

        real(kind(0.d0)) :: f_quad_RV, f_quad
        integer :: i

        f_quad = 0d0
        do i = 1, nb
            f_quad_RV = sum(wght_in(:, i)*(abscX(:, i)**q)*(abscY(:, i)**r))
            f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
        end do

    end function f_quad

    function f_quad2D(abscX, abscY, wght_in, pow)
        !$acc routine seq
        real(kind(0.d0)), dimension(nnode), intent(in) :: abscX, abscY, wght_in
        real(kind(0.d0)), dimension(3), intent(in) :: pow

        real(kind(0.d0)) :: f_quad2D

        f_quad2D = sum(wght_in(:)*(abscX(:)**pow(1))*(abscY(:)**pow(2)))
    end function f_quad2D

end module m_qbmm
