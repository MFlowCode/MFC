!>
!! @file m_qbmm.f90
!! @brief Contains module m_qbmm

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief This module is used to compute moment inversion via qbmm
module m_qbmm

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic           !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_initialize_qbmm_module, s_mom_inv, s_coeff, s_compute_qbmm_rhs

    real(wp), allocatable, dimension(:, :, :, :, :) :: momrhs
    $:GPU_DECLARE(create='[momrhs]')

    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: nterms = ${nterms}$
    #:else
        integer :: nterms
        $:GPU_DECLARE(create='[nterms]')
    #:endif

    type(int_bounds_info) :: is1_qbmm, is2_qbmm, is3_qbmm
    $:GPU_DECLARE(create='[is1_qbmm,is2_qbmm,is3_qbmm]')

    integer, allocatable, dimension(:) :: bubrs
    integer, allocatable, dimension(:, :) :: bubmoms
    $:GPU_DECLARE(create='[bubrs,bubmoms]')

contains

    impure subroutine s_initialize_qbmm_module

        integer :: i1, i2, q, i, j

        #:if not MFC_CASE_OPTIMIZATION

            if (bubble_model == 2) then
                ! Keller-Miksis without viscosity/surface tension
                nterms = 32
            else if (bubble_model == 3) then
                ! Rayleigh-Plesset with viscosity/surface tension
                nterms = 7
            end if

            $:GPU_ENTER_DATA(copyin='[nterms]')
            $:GPU_UPDATE(device='[nterms]')

        #:endif

        @:ALLOCATE(momrhs(1:3, 0:2, 0:2, 1:nterms, 1:nb))
        momrhs = 0._wp

        ! Assigns the required RHS moments for moment transport equations
        ! The rhs%(:,3) is only to be used for R0 quadrature, not for computing X/Y indices
        ! Accounts for different governing equations in polytropic and non-polytropic models
        if (.not. polytropic) then
            do q = 1, nb
                do i1 = 0, 2; do i2 = 0, 2
                        if ((i1 + i2) <= 2) then
                            if (bubble_model == 3) then
                                momrhs(1, i1, i2, 1, q) = -1._wp + i1
                                momrhs(2, i1, i2, 1, q) = -1._wp + i2
                                momrhs(3, i1, i2, 1, q) = 0._wp

                                momrhs(1, i1, i2, 2, q) = -1._wp + i1
                                momrhs(2, i1, i2, 2, q) = 1._wp + i2
                                momrhs(3, i1, i2, 2, q) = 0._wp

                                momrhs(1, i1, i2, 3, q) = -1._wp + i1
                                momrhs(2, i1, i2, 3, q) = -1._wp + i2
                                momrhs(3, i1, i2, 3, q) = 0._wp

                                momrhs(1, i1, i2, 4, q) = -1._wp + i1
                                momrhs(2, i1, i2, 4, q) = 1._wp + i2
                                momrhs(3, i1, i2, 4, q) = 0._wp

                                if (.not. f_is_default(Re_inv)) then
                                    ! add viscosity
                                    momrhs(1, i1, i2, 5, q) = -2._wp + i1
                                    momrhs(2, i1, i2, 5, q) = i2
                                    momrhs(3, i1, i2, 5, q) = 0._wp
                                end if

                                if (.not. f_is_default(Web)) then
                                    ! add surface tension
                                    momrhs(1, i1, i2, 6, q) = -2._wp + i1
                                    momrhs(2, i1, i2, 6, q) = -1._wp + i2
                                    momrhs(3, i1, i2, 6, q) = 0._wp
                                end if

                                momrhs(1, i1, i2, 7, q) = -1._wp + i1
                                momrhs(2, i1, i2, 7, q) = -1._wp + i2
                                momrhs(3, i1, i2, 7, q) = 0._wp

                            else if (bubble_model == 2) then
                                ! KM with approximation of 1/(1-V/C) = 1+V/C
                                momrhs(1, i1, i2, 1, q) = -1._wp + i1
                                momrhs(2, i1, i2, 1, q) = 1._wp + i2
                                momrhs(3, i1, i2, 1, q) = 0._wp

                                momrhs(1, i1, i2, 2, q) = -1._wp + i1
                                momrhs(2, i1, i2, 2, q) = 2._wp + i2
                                momrhs(3, i1, i2, 2, q) = 0._wp

                                momrhs(1, i1, i2, 3, q) = -1._wp + i1
                                momrhs(2, i1, i2, 3, q) = 3._wp + i2
                                momrhs(3, i1, i2, 3, q) = 0._wp

                                momrhs(1, i1, i2, 4, q) = -1._wp + i1
                                momrhs(2, i1, i2, 4, q) = -1._wp + i2
                                momrhs(3, i1, i2, 4, q) = 0._wp

                                momrhs(1, i1, i2, 5, q) = -1._wp + i1
                                momrhs(2, i1, i2, 5, q) = i2
                                momrhs(3, i1, i2, 5, q) = 0._wp

                                momrhs(1, i1, i2, 6, q) = -1._wp + i1
                                momrhs(2, i1, i2, 6, q) = 1._wp + i2
                                momrhs(3, i1, i2, 6, q) = 0._wp

                                momrhs(1, i1, i2, 7, q) = -1._wp + i1
                                momrhs(2, i1, i2, 7, q) = -1._wp + i2
                                momrhs(3, i1, i2, 7, q) = 0._wp

                                momrhs(1, i1, i2, 8, q) = -1._wp + i1
                                momrhs(2, i1, i2, 8, q) = i2
                                momrhs(3, i1, i2, 8, q) = 0._wp

                                momrhs(1, i1, i2, 9, q) = -1._wp + i1
                                momrhs(2, i1, i2, 9, q) = 1._wp + i2
                                momrhs(3, i1, i2, 9, q) = 0._wp

                                momrhs(1, i1, i2, 10, q) = -1._wp + i1
                                momrhs(2, i1, i2, 10, q) = i2
                                momrhs(3, i1, i2, 10, q) = 0._wp

                                momrhs(1, i1, i2, 11, q) = -1._wp + i1
                                momrhs(2, i1, i2, 11, q) = 1._wp + i2
                                momrhs(3, i1, i2, 11, q) = 0._wp

                                momrhs(1, i1, i2, 12, q) = -1._wp + i1
                                momrhs(2, i1, i2, 12, q) = 1._wp + i2
                                momrhs(3, i1, i2, 12, q) = 0._wp

                                momrhs(1, i1, i2, 13, q) = -1._wp + i1
                                momrhs(2, i1, i2, 13, q) = -1._wp + i2
                                momrhs(3, i1, i2, 13, q) = 0._wp

                                momrhs(1, i1, i2, 14, q) = -1._wp + i1
                                momrhs(2, i1, i2, 14, q) = i2
                                momrhs(3, i1, i2, 14, q) = 0._wp

                                momrhs(1, i1, i2, 15, q) = -1._wp + i1
                                momrhs(2, i1, i2, 15, q) = 1._wp + i2
                                momrhs(3, i1, i2, 15, q) = 0._wp

                                momrhs(1, i1, i2, 16, q) = -2._wp + i1
                                momrhs(2, i1, i2, 16, q) = i2
                                momrhs(3, i1, i2, 16, q) = 0._wp

                                momrhs(1, i1, i2, 17, q) = -2._wp + i1
                                momrhs(2, i1, i2, 17, q) = -1._wp + i2
                                momrhs(3, i1, i2, 17, q) = 0._wp

                                momrhs(1, i1, i2, 18, q) = -2._wp + i1
                                momrhs(2, i1, i2, 18, q) = 1._wp + i2
                                momrhs(3, i1, i2, 18, q) = 0._wp

                                momrhs(1, i1, i2, 19, q) = -2._wp + i1
                                momrhs(2, i1, i2, 19, q) = 2._wp + i2
                                momrhs(3, i1, i2, 19, q) = 0._wp

                                momrhs(1, i1, i2, 20, q) = -2._wp + i1
                                momrhs(2, i1, i2, 20, q) = -1._wp + i2
                                momrhs(3, i1, i2, 20, q) = 0._wp

                                momrhs(1, i1, i2, 21, q) = -2._wp + i1
                                momrhs(2, i1, i2, 21, q) = i2
                                momrhs(3, i1, i2, 21, q) = 0._wp

                                momrhs(1, i1, i2, 22, q) = -2._wp + i1
                                momrhs(2, i1, i2, 22, q) = -1._wp + i2
                                momrhs(3, i1, i2, 22, q) = 0._wp

                                momrhs(1, i1, i2, 23, q) = -2._wp + i1
                                momrhs(2, i1, i2, 23, q) = i2
                                momrhs(3, i1, i2, 23, q) = 0._wp

                                momrhs(1, i1, i2, 24, q) = -3._wp + i1
                                momrhs(2, i1, i2, 24, q) = i2
                                momrhs(3, i1, i2, 24, q) = 0._wp

                                momrhs(1, i1, i2, 25, q) = -3._wp + i1
                                momrhs(2, i1, i2, 25, q) = -1._wp + i2
                                momrhs(3, i1, i2, 25, q) = 0._wp

                                momrhs(1, i1, i2, 26, q) = -2._wp + i1
                                momrhs(2, i1, i2, 26, q) = i2
                                momrhs(3, i1, i2, 26, q) = 0._wp

                                momrhs(1, i1, i2, 27, q) = -1._wp + i1
                                momrhs(2, i1, i2, 27, q) = -1._wp + i2
                                momrhs(3, i1, i2, 27, q) = 0._wp

                                momrhs(1, i1, i2, 28, q) = -1._wp + i1
                                momrhs(2, i1, i2, 28, q) = i2
                                momrhs(3, i1, i2, 28, q) = 0._wp

                                momrhs(1, i1, i2, 29, q) = -2._wp + i1
                                momrhs(2, i1, i2, 29, q) = i2
                                momrhs(3, i1, i2, 29, q) = 0._wp

                                momrhs(1, i1, i2, 30, q) = -1._wp + i1
                                momrhs(2, i1, i2, 30, q) = -1._wp + i2
                                momrhs(3, i1, i2, 30, q) = 0._wp

                                momrhs(1, i1, i2, 31, q) = -1._wp + i1
                                momrhs(2, i1, i2, 31, q) = i2
                                momrhs(3, i1, i2, 31, q) = 0._wp

                                momrhs(1, i1, i2, 32, q) = -2._wp + i1
                                momrhs(2, i1, i2, 32, q) = i2
                                momrhs(3, i1, i2, 32, q) = 0._wp
                            end if
                        end if
                    end do; end do
            end do

        else
            do q = 1, nb
                do i1 = 0, 2; do i2 = 0, 2
                        if ((i1 + i2) <= 2) then
                            if (bubble_model == 3) then
                                momrhs(1, i1, i2, 1, q) = -1._wp + i1
                                momrhs(2, i1, i2, 1, q) = -1._wp + i2
                                momrhs(3, i1, i2, 1, q) = 0._wp

                                momrhs(1, i1, i2, 2, q) = -1._wp + i1
                                momrhs(2, i1, i2, 2, q) = 1._wp + i2
                                momrhs(3, i1, i2, 2, q) = 0._wp

                                momrhs(1, i1, i2, 3, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 3, q) = -1._wp + i2
                                momrhs(3, i1, i2, 3, q) = 3._wp*gam

                                momrhs(1, i1, i2, 4, q) = -1._wp + i1
                                momrhs(2, i1, i2, 4, q) = 1._wp + i2
                                momrhs(3, i1, i2, 4, q) = 0._wp

                                if (.not. f_is_default(Re_inv)) then
                                    ! add viscosity
                                    momrhs(1, i1, i2, 5, q) = -2._wp + i1
                                    momrhs(2, i1, i2, 5, q) = i2
                                    momrhs(3, i1, i2, 5, q) = 0._wp
                                end if

                                if (.not. f_is_default(Web)) then
                                    ! add surface tension
                                    momrhs(1, i1, i2, 6, q) = -2._wp + i1
                                    momrhs(2, i1, i2, 6, q) = -1._wp + i2
                                    momrhs(3, i1, i2, 6, q) = 0._wp
                                end if

                                momrhs(1, i1, i2, 7, q) = -1._wp + i1
                                momrhs(2, i1, i2, 7, q) = -1._wp + i2
                                momrhs(3, i1, i2, 7, q) = 0._wp

                            else if (bubble_model == 2) then
                                ! KM with approximation of 1/(1-V/C) = 1+V/C
                                momrhs(1, i1, i2, 1, q) = -1._wp + i1
                                momrhs(2, i1, i2, 1, q) = 1._wp + i2
                                momrhs(3, i1, i2, 1, q) = 0._wp

                                momrhs(1, i1, i2, 2, q) = -1._wp + i1
                                momrhs(2, i1, i2, 2, q) = 2._wp + i2
                                momrhs(3, i1, i2, 2, q) = 0._wp

                                momrhs(1, i1, i2, 3, q) = -1._wp + i1
                                momrhs(2, i1, i2, 3, q) = 3._wp + i2
                                momrhs(3, i1, i2, 3, q) = 0._wp

                                momrhs(1, i1, i2, 4, q) = -1._wp + i1
                                momrhs(2, i1, i2, 4, q) = -1._wp + i2
                                momrhs(3, i1, i2, 4, q) = 0._wp

                                momrhs(1, i1, i2, 5, q) = -1._wp + i1
                                momrhs(2, i1, i2, 5, q) = i2
                                momrhs(3, i1, i2, 5, q) = 0._wp

                                momrhs(1, i1, i2, 6, q) = -1._wp + i1
                                momrhs(2, i1, i2, 6, q) = 1._wp + i2
                                momrhs(3, i1, i2, 6, q) = 0._wp

                                momrhs(1, i1, i2, 7, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 7, q) = -1._wp + i2
                                momrhs(3, i1, i2, 7, q) = 3._wp*gam

                                momrhs(1, i1, i2, 8, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 8, q) = i2
                                momrhs(3, i1, i2, 8, q) = 3._wp*gam

                                momrhs(1, i1, i2, 9, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 9, q) = 1._wp + i2
                                momrhs(3, i1, i2, 9, q) = 3._wp*gam

                                momrhs(1, i1, i2, 10, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 10, q) = i2
                                momrhs(3, i1, i2, 10, q) = 3._wp*gam

                                momrhs(1, i1, i2, 11, q) = -1._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 11, q) = 1._wp + i2
                                momrhs(3, i1, i2, 11, q) = 3._wp*gam

                                momrhs(1, i1, i2, 12, q) = -1._wp + i1
                                momrhs(2, i1, i2, 12, q) = 1._wp + i2
                                momrhs(3, i1, i2, 12, q) = 0._wp

                                momrhs(1, i1, i2, 13, q) = -1._wp + i1
                                momrhs(2, i1, i2, 13, q) = -1._wp + i2
                                momrhs(3, i1, i2, 13, q) = 0._wp

                                momrhs(1, i1, i2, 14, q) = -1._wp + i1
                                momrhs(2, i1, i2, 14, q) = i2
                                momrhs(3, i1, i2, 14, q) = 0._wp

                                momrhs(1, i1, i2, 15, q) = -1._wp + i1
                                momrhs(2, i1, i2, 15, q) = 1._wp + i2
                                momrhs(3, i1, i2, 15, q) = 0._wp

                                momrhs(1, i1, i2, 16, q) = -2._wp + i1
                                momrhs(2, i1, i2, 16, q) = i2
                                momrhs(3, i1, i2, 16, q) = 0._wp

                                momrhs(1, i1, i2, 17, q) = -2._wp + i1
                                momrhs(2, i1, i2, 17, q) = -1._wp + i2
                                momrhs(3, i1, i2, 17, q) = 0._wp

                                momrhs(1, i1, i2, 18, q) = -2._wp + i1
                                momrhs(2, i1, i2, 18, q) = 1._wp + i2
                                momrhs(3, i1, i2, 18, q) = 0._wp

                                momrhs(1, i1, i2, 19, q) = -2._wp + i1
                                momrhs(2, i1, i2, 19, q) = 2._wp + i2
                                momrhs(3, i1, i2, 19, q) = 0._wp

                                momrhs(1, i1, i2, 20, q) = -2._wp + i1
                                momrhs(2, i1, i2, 20, q) = -1._wp + i2
                                momrhs(3, i1, i2, 20, q) = 0._wp

                                momrhs(1, i1, i2, 21, q) = -2._wp + i1
                                momrhs(2, i1, i2, 21, q) = i2
                                momrhs(3, i1, i2, 21, q) = 0._wp

                                momrhs(1, i1, i2, 22, q) = -2._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 22, q) = -1._wp + i2
                                momrhs(3, i1, i2, 22, q) = 3._wp*gam

                                momrhs(1, i1, i2, 23, q) = -2._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 23, q) = i2
                                momrhs(3, i1, i2, 23, q) = 3._wp*gam

                                momrhs(1, i1, i2, 24, q) = -3._wp + i1
                                momrhs(2, i1, i2, 24, q) = i2
                                momrhs(3, i1, i2, 24, q) = 0._wp

                                momrhs(1, i1, i2, 25, q) = -3._wp + i1
                                momrhs(2, i1, i2, 25, q) = -1._wp + i2
                                momrhs(3, i1, i2, 25, q) = 0._wp

                                momrhs(1, i1, i2, 26, q) = -2._wp + i1 - 3._wp*gam
                                momrhs(2, i1, i2, 26, q) = i2
                                momrhs(3, i1, i2, 26, q) = 3._wp*gam

                            end if
                        end if
                    end do; end do
            end do
        end if

        $:GPU_UPDATE(device='[momrhs]')

        @:ALLOCATE(bubrs(1:nb))
        @:ALLOCATE(bubmoms(1:nb, 1:nmom))

        do i = 1, nb
            bubrs(i) = bub_idx%rs(i)
        end do
        $:GPU_UPDATE(device='[bubrs]')

        do j = 1, nmom
            do i = 1, nb
                bubmoms(i, j) = bub_idx%moms(i, j)
            end do
        end do
        $:GPU_UPDATE(device='[bubmoms]')

    end subroutine s_initialize_qbmm_module

    pure subroutine s_compute_qbmm_rhs(idir, q_cons_vf, q_prim_vf, rhs_vf, flux_n_vf, pb, rhs_pb)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in) :: flux_n_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, rhs_pb

        integer :: i, j, k, l, q
        real(wp) :: nb_q, nb_dot, R, R2, nR, nR2, nR_dot, nR2_dot, var, AX
        logical :: is_axisym

        select case (idir)
        case (1)
            is_axisym = .false.
        case (2)
            is_axisym = .false.
        case (3)
            is_axisym = (grid_geometry == 3)
        end select

        if (.not. polytropic) then
            $:GPU_PARALLEL_LOOP(collapse=5,private='[nb_q,nR,nR2,R,R2,nb_dot,nR_dot,nR2_dot,var,AX]')
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
                                var = max(R2 - R**2._wp, verysmall)
                                if (q <= 2) then
                                    AX = R - sqrt(var)
                                else
                                    AX = R + sqrt(var)
                                end if
                                select case (idir)
                                case (1)
                                    nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j - 1, k, l) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)
                                    rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dx(j)*AX*nb_q**2)* &
                                                            (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))
                                case (2)
                                    nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                    nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                    nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k - 1, l) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)
                                    rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dy(k)*AX*nb_q**2)* &
                                                            (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))
                                case (3)
                                    if (is_axisym) then
                                        nb_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l))
                                        nR_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l))
                                        nR2_dot = q_prim_vf(contxe + idir)%sf(j, k, l)*(flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*y_cc(k)*AX*nb_q**2)* &
                                                                (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))
                                    else
                                        nb_dot = flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + (i - 1)*nmom)%sf(j, k, l)
                                        nR_dot = flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)
                                        nR2_dot = flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l - 1) - flux_n_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*AX*nb_q**2)* &
                                                                (nR_dot*nb_q - nR*nb_dot)*(pb(j, k, l, q, i))
                                    end if
                                end select
                                if (q <= 2) then
                                    select case (idir)
                                    case (1)
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dx(j)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dx(j)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    case (2)
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dy(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dy(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    case (3)
                                        if (is_axisym) then
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dz(l)*y_cc(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dz(l)*y_cc(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        else
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dz(l)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) + 3._wp*gam/(dz(l)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        end if
                                    end select
                                else
                                    select case (idir)
                                    case (1)
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dx(j)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dx(j)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    case (2)
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dy(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                        rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dy(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                    case (3)
                                        if (is_axisym) then
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*y_cc(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*y_cc(k)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        else
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (nR2_dot*nb_q - nR2*nb_dot)*(pb(j, k, l, q, i))
                                            rhs_pb(j, k, l, q, i) = rhs_pb(j, k, l, q, i) - 3._wp*gam/(dz(l)*AX*nb_q**2*sqrt(var)*2._wp)* &
                                                                    (-2._wp*(nR/nb_q)*(nR_dot*nb_q - nR*nb_dot))*(pb(j, k, l, q, i))
                                        end if
                                    end select
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! The following block is not repeated and is left as is
        if (idir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do q = 0, n
                    do i = 0, m
                        rhs_vf(alf_idx)%sf(i, q, l) = rhs_vf(alf_idx)%sf(i, q, l) + mom_sp(2)%sf(i, q, l)
                        j = bubxb
                        $:GPU_LOOP(parallelism='[seq]')
                        do k = 1, nb
                            rhs_vf(j)%sf(i, q, l) = rhs_vf(j)%sf(i, q, l) + mom_3d(0, 0, k)%sf(i, q, l)
                            rhs_vf(j + 1)%sf(i, q, l) = rhs_vf(j + 1)%sf(i, q, l) + mom_3d(1, 0, k)%sf(i, q, l)
                            rhs_vf(j + 2)%sf(i, q, l) = rhs_vf(j + 2)%sf(i, q, l) + mom_3d(0, 1, k)%sf(i, q, l)
                            rhs_vf(j + 3)%sf(i, q, l) = rhs_vf(j + 3)%sf(i, q, l) + mom_3d(2, 0, k)%sf(i, q, l)
                            rhs_vf(j + 4)%sf(i, q, l) = rhs_vf(j + 4)%sf(i, q, l) + mom_3d(1, 1, k)%sf(i, q, l)
                            rhs_vf(j + 5)%sf(i, q, l) = rhs_vf(j + 5)%sf(i, q, l) + mom_3d(0, 2, k)%sf(i, q, l)
                            j = j + 6
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_compute_qbmm_rhs

    !Coefficient array for non-polytropic model (pb and mv values are accounted in wght_pb and wght_mv)
    pure subroutine s_coeff_nonpoly(pres, rho, c, coeffs)
        $:GPU_ROUTINE(function_name='s_coeff_nonpoly',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(in) :: pres, rho, c
        real(wp), dimension(nterms, 0:2, 0:2), intent(out) :: coeffs

        integer :: i1, i2

        coeffs = 0._wp

        do i2 = 0, 2; do i1 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1._wp*i2*pres/rho
                        coeffs(2, i1, i2) = -3._wp*i2/2._wp
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (.not. f_is_default(Re_inv)) coeffs(5, i1, i2) = -4._wp*i2*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(6, i1, i2) = -2._wp*i2/Web/rho
                        coeffs(7, i1, i2) = 0._wp
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3._wp*i2/2._wp
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2._wp*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2._wp*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2._wp*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3._wp*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3._wp*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                        coeffs(13, i1, i2) = 0._wp
                        coeffs(14, i1, i2) = 0._wp
                        coeffs(15, i1, i2) = 0._wp
                        if (.not. f_is_default(Re_inv)) coeffs(16, i1, i2) = -i2*4._wp*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(17, i1, i2) = -i2*2._wp/Web/rho
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(18, i1, i2) = i2*6._wp*Re_inv/(rho*c)
                            coeffs(19, i1, i2) = -i2*2._wp*Re_inv/(rho*c*c)
                            coeffs(20, i1, i2) = i2*4._wp*pres*Re_inv/(rho*rho*c)
                            coeffs(21, i1, i2) = i2*4._wp*pres*Re_inv/(rho*rho*c*c)
                            coeffs(22, i1, i2) = -i2*4._wp*Re_inv/(rho*rho*c)
                            coeffs(23, i1, i2) = -i2*4._wp*Re_inv/(rho*rho*c*c)
                            coeffs(24, i1, i2) = i2*16._wp*Re_inv*Re_inv/(rho*rho*c)
                            if (.not. f_is_default(Web)) then
                                coeffs(25, i1, i2) = i2*8._wp*Re_inv/Web/(rho*rho*c)
                            end if
                            coeffs(26, i1, i2) = -12._wp*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                        coeffs(27, i1, i2) = 3._wp*i2*gam*R_v*Tw/(c*rho)
                        coeffs(28, i1, i2) = 3._wp*i2*gam*R_v*Tw/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(29, i1, i2) = 12._wp*i2*gam*R_v*Tw*Re_inv/(rho*rho*c*c)
                        end if
                        coeffs(30, i1, i2) = 3._wp*i2*gam/(c*rho)
                        coeffs(31, i1, i2) = 3._wp*i2*gam/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(32, i1, i2) = 12._wp*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                    end if
                end if
            end do; end do

    end subroutine s_coeff_nonpoly

!Coefficient array for polytropic model (pb for each R0 bin accounted for in wght_pb)
    pure subroutine s_coeff(pres, rho, c, coeffs)
        $:GPU_ROUTINE(function_name='s_coeff',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(in) :: pres, rho, c
        real(wp), dimension(nterms, 0:2, 0:2), intent(out) :: coeffs

        integer :: i1, i2

        coeffs = 0._wp

        do i2 = 0, 2; do i1 = 0, 2
                if ((i1 + i2) <= 2) then
                    if (bubble_model == 3) then
                        ! RPE
                        coeffs(1, i1, i2) = -1._wp*i2*pres/rho
                        coeffs(2, i1, i2) = -3._wp*i2/2._wp
                        coeffs(3, i1, i2) = i2/rho
                        coeffs(4, i1, i2) = i1
                        if (.not. f_is_default(Re_inv)) coeffs(5, i1, i2) = -4._wp*i2*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(6, i1, i2) = -2._wp*i2/Web/rho
                        coeffs(7, i1, i2) = i2*pv/rho
                    else if (bubble_model == 2) then
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1, i1, i2) = -3._wp*i2/2._wp
                        coeffs(2, i1, i2) = -i2/c
                        coeffs(3, i1, i2) = i2/(2._wp*c*c)
                        coeffs(4, i1, i2) = -i2*pres/rho
                        coeffs(5, i1, i2) = -2._wp*i2*pres/(c*rho)
                        coeffs(6, i1, i2) = -i2*pres/(c*c*rho)
                        coeffs(7, i1, i2) = i2/rho
                        coeffs(8, i1, i2) = 2._wp*i2/(c*rho)
                        coeffs(9, i1, i2) = i2/(c*c*rho)
                        coeffs(10, i1, i2) = -3._wp*i2*gam/(c*rho)
                        coeffs(11, i1, i2) = -3._wp*i2*gam/(c*c*rho)
                        coeffs(12, i1, i2) = i1
                        coeffs(13, i1, i2) = i2*(pv)/rho
                        coeffs(14, i1, i2) = 2._wp*i2*(pv)/(c*rho)
                        coeffs(15, i1, i2) = i2*(pv)/(c*c*rho)
                        if (.not. f_is_default(Re_inv)) coeffs(16, i1, i2) = -i2*4._wp*Re_inv/rho
                        if (.not. f_is_default(Web)) coeffs(17, i1, i2) = -i2*2._wp/Web/rho
                        if (.not. f_is_default(Re_inv)) then
                            coeffs(18, i1, i2) = i2*6._wp*Re_inv/(rho*c)
                            coeffs(19, i1, i2) = -i2*2._wp*Re_inv/(rho*c*c)
                            coeffs(20, i1, i2) = i2*4._wp*pres*Re_inv/(rho*rho*c)
                            coeffs(21, i1, i2) = i2*4._wp*pres*Re_inv/(rho*rho*c*c)
                            coeffs(22, i1, i2) = -i2*4._wp*Re_inv/(rho*rho*c)
                            coeffs(23, i1, i2) = -i2*4._wp*Re_inv/(rho*rho*c*c)
                            coeffs(24, i1, i2) = i2*16._wp*Re_inv*Re_inv/(rho*rho*c)
                            if (.not. f_is_default(Web)) then
                                coeffs(25, i1, i2) = i2*8._wp*Re_inv/Web/(rho*rho*c)
                            end if
                            coeffs(26, i1, i2) = -12._wp*i2*gam*Re_inv/(rho*rho*c*c)
                        end if
                    end if
                end if
            end do; end do

    end subroutine s_coeff

    subroutine s_mom_inv(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix, iy, iz)

        type(scalar_field), dimension(:), intent(inout) :: q_cons_vf, q_prim_vf
        type(scalar_field), dimension(:), intent(inout) :: momsp
        type(scalar_field), dimension(0:, 0:, :), intent(inout) :: moms3d
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: pb, rhs_pb
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:, 1:), intent(inout) :: mv, rhs_mv
        type(int_bounds_info), intent(in) :: ix, iy, iz

        real(wp), dimension(nmom) :: moms, msum
        real(wp), dimension(nnode, nb) :: wght, abscX, abscY, wght_pb, wght_mv, wght_ht, ht
        real(wp), dimension(nterms, 0:2, 0:2) :: coeff
        real(wp) :: pres, rho, nbub, c, alf, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, T_bar, grad_T
        real(wp) :: n_tait, B_tait
        integer :: id1, id2, id3, i1, i2, j, q, r

        is1_qbmm = ix; is2_qbmm = iy; is3_qbmm = iz
        $:GPU_UPDATE(device='[is1_qbmm,is2_qbmm,is3_qbmm]')

        $:GPU_PARALLEL_LOOP(collapse=3, private='[moms, msum, wght, abscX, &
            & abscY, wght_pb, wght_mv, wght_ht, coeff, ht, r, q, &
            & n_tait, B_tait, pres, rho, nbub, c, alf, momsum, &
            & drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, T_bar, grad_T]')
        do id3 = is3_qbmm%beg, is3_qbmm%end
            do id2 = is2_qbmm%beg, is2_qbmm%end
                do id1 = is1_qbmm%beg, is1_qbmm%end

                    alf = q_prim_vf(alf_idx)%sf(id1, id2, id3)
                    pres = q_prim_vf(E_idx)%sf(id1, id2, id3)
                    rho = q_prim_vf(contxb)%sf(id1, id2, id3)

                    if (bubble_model == 2) then
                        n_tait = 1._wp/gammas(1) + 1._wp
                        B_tait = pi_infs(1)*(n_tait - 1)/n_tait
                        c = n_tait*(pres + B_tait)*(1._wp - alf)/(rho)
                        c = merge(sqrt(c), sgm_eps, c > 0._wp)
                    end if

                    call s_coeff_selector(pres, rho, c, coeff, polytropic)

                    if (alf > small_alf) then
                        nbub = q_cons_vf(bubxb)%sf(id1, id2, id3)
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            ! Gather moments for this bubble bin
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = 2, nmom
                                moms(r) = q_prim_vf(bubmoms(q, r))%sf(id1, id2, id3)
                            end do
                            moms(1) = 1._wp
                            call s_chyqmom(moms, wght(:, q), abscX(:, q), abscY(:, q))

                            if (polytropic) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, nnode
                                    wght_pb(j, q) = wght(j, q)*(pb0(q) - pv)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, nnode
                                    chi_vw = 1._wp/(1._wp + R_v/R_n*(pb(id1, id2, id3, j, q)/pv - 1._wp))
                                    x_vw = M_n*chi_vw/(M_v + (M_n - M_v)*chi_vw)
                                    k_mw = x_vw*k_v(q)/(x_vw + (1._wp - x_vw)*phi_vn) + (1._wp - x_vw)*k_n(q)/(x_vw*phi_nv + 1._wp - x_vw)
                                    rho_mw = pv/(chi_vw*R_v*Tw)
                                    rhs_mv(id1, id2, id3, j, q) = -Re_trans_c(q)*((mv(id1, id2, id3, j, q)/(mv(id1, id2, id3, j, q) + mass_n0(q))) - chi_vw)
                                    rhs_mv(id1, id2, id3, j, q) = rho_mw*rhs_mv(id1, id2, id3, j, q)/Pe_c/(1._wp - chi_vw)/abscX(j, q)
                                    T_bar = Tw*(pb(id1, id2, id3, j, q)/pb0(q))*(abscX(j, q)/R0(q))**3*(mass_n0(q) + mass_v0(q))/(mass_n0(q) + mv(id1, id2, id3, j, q))
                                    grad_T = -Re_trans_T(q)*(T_bar - Tw)
                                    ht(j, q) = pb0(q)*k_mw*grad_T/Pe_T(q)/abscX(j, q)
                                    wght_pb(j, q) = wght(j, q)*(pb(id1, id2, id3, j, q))
                                    wght_mv(j, q) = wght(j, q)*(rhs_mv(id1, id2, id3, j, q))
                                    wght_ht(j, q) = wght(j, q)*ht(j, q)
                                end do
                            end if

                            ! Compute change in moments due to bubble dynamics
                            r = 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do i2 = 0, 2
                                $:GPU_LOOP(parallelism='[seq]')
                                do i1 = 0, 2
                                    if ((i1 + i2) <= 2) then
                                        momsum = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do j = 1, nterms
                                            select case (bubble_model)
                                            case (3)
                                                if (j == 3) then
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght_pb(:, q), momrhs(:, i1, i2, j, q))
                                                else
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght(:, q), momrhs(:, i1, i2, j, q))
                                                end if
                                            case (2)
                                                if ((j >= 7 .and. j <= 9) .or. (j >= 22 .and. j <= 23) .or. (j >= 10 .and. j <= 11) .or. (j == 26)) then
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght_pb(:, q), momrhs(:, i1, i2, j, q))
                                                else if ((j >= 27 .and. j <= 29) .and. (.not. polytropic)) then
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght_mv(:, q), momrhs(:, i1, i2, j, q))
                                                else if ((j >= 30 .and. j <= 32) .and. (.not. polytropic)) then
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght_ht(:, q), momrhs(:, i1, i2, j, q))
                                                else
                                                    momsum = momsum + coeff(j, i1, i2)*(R0(q)**momrhs(3, i1, i2, j, q))*f_quad2D(abscX(:, q), abscY(:, q), wght(:, q), momrhs(:, i1, i2, j, q))
                                                end if
                                            end select
                                        end do
                                        moms3d(i1, i2, q)%sf(id1, id2, id3) = nbub*momsum
                                        msum(r) = momsum
                                        r = r + 1
                                    end if
                                end do
                            end do

                            ! Compute change in pb and mv for non-polytropic model
                            if (.not. polytropic) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, nnode
                                    drdt = msum(2)
                                    drdt2 = merge(-1._wp, 1._wp, j == 1 .or. j == 2)/(2._wp*sqrt(merge(moms(4) - moms(2)**2._wp, verysmall, moms(4) - moms(2)**2._wp > 0._wp)))
                                    drdt2 = drdt2*(msum(3) - 2._wp*moms(2)*msum(2))
                                    drdt = drdt + drdt2
                                    rhs_pb(id1, id2, id3, j, q) = (-3._wp*gam*drdt/abscX(j, q))*(pb(id1, id2, id3, j, q))
                                    rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, q) + (3._wp*gam/abscX(j, q))*rhs_mv(id1, id2, id3, j, q)*R_v*Tw
                                    rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, q) + (3._wp*gam/abscX(j, q))*ht(j, q)
                                    rhs_mv(id1, id2, id3, j, q) = rhs_mv(id1, id2, id3, j, q)*(4._wp*pi*abscX(j, q)**2._wp)
                                end do
                            end if
                        end do

                        ! Compute special high-order moments
                        momsp(1)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3._wp, 0._wp, 0._wp)
                        momsp(2)%sf(id1, id2, id3) = 4._wp*pi*nbub*f_quad(abscX, abscY, wght, 2._wp, 1._wp, 0._wp)
                        momsp(3)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght, 3._wp, 2._wp, 0._wp)
                        if (abs(gam - 1._wp) <= 1.e-4_wp) then
                            momsp(4)%sf(id1, id2, id3) = 1._wp
                        else
                            if (polytropic) then
                                momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght_pb, 3._wp*(1._wp - gam), 0._wp, 3._wp*gam) + pv*f_quad(abscX, abscY, wght, 3._wp, 0._wp, 0._wp) - 4._wp*Re_inv*f_quad(abscX, abscY, wght, 2._wp, 1._wp, 0._wp) - (2._wp/Web)*f_quad(abscX, abscY, wght, 2._wp, 0._wp, 0._wp)
                            else
                                momsp(4)%sf(id1, id2, id3) = f_quad(abscX, abscY, wght_pb, 3._wp, 0._wp, 0._wp) - 4._wp*Re_inv*f_quad(abscX, abscY, wght, 2._wp, 1._wp, 0._wp) - (2._wp/Web)*f_quad(abscX, abscY, wght, 2._wp, 0._wp, 0._wp)
                            end if
                        end if
                    else
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            $:GPU_LOOP(parallelism='[seq]')
                            do i1 = 0, 2
                                $:GPU_LOOP(parallelism='[seq]')
                                do i2 = 0, 2
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = 0._wp
                                end do
                            end do
                        end do
                        momsp(1)%sf(id1, id2, id3) = 0._wp
                        momsp(2)%sf(id1, id2, id3) = 0._wp
                        momsp(3)%sf(id1, id2, id3) = 0._wp
                        momsp(4)%sf(id1, id2, id3) = 0._wp
                    end if
                end do
            end do
        end do

    contains
        ! Helper to select the correct coefficient routine
        subroutine s_coeff_selector(pres, rho, c, coeff, polytropic)
            $:GPU_ROUTINE(function_name='s_coeff_selector',parallelism='[seq]', &
                & cray_inline=True)
            real(wp), intent(in) :: pres, rho, c
            real(wp), dimension(nterms, 0:2, 0:2), intent(out) :: coeff
            logical, intent(in) :: polytropic
            if (polytropic) then
                call s_coeff(pres, rho, c, coeff)
            else
                call s_coeff_nonpoly(pres, rho, c, coeff)
            end if
        end subroutine s_coeff_selector

        pure subroutine s_chyqmom(momin, wght, abscX, abscY)
            $:GPU_ROUTINE(function_name='s_chyqmom',parallelism='[seq]', &
                & cray_inline=True)

            real(wp), dimension(nmom), intent(in) :: momin
            real(wp), dimension(nnode), intent(inout) :: wght, abscX, abscY

            ! Local variables
            real(wp), dimension(0:2, 0:2) :: moms
            real(wp), dimension(3) :: M1, M3
            real(wp), dimension(2) :: myrho, myrho3, up, up3, Vf
            real(wp) :: bu, bv, d20, d11, d_02, c20, c11, c02
            real(wp) :: mu2, vp21, vp22, rho21, rho22

            ! Assign moments to 2D array for clarity
            moms(0, 0) = momin(1)
            moms(1, 0) = momin(2)
            moms(0, 1) = momin(3)
            moms(2, 0) = momin(4)
            moms(1, 1) = momin(5)
            moms(0, 2) = momin(6)

            ! Compute means and central moments
            bu = moms(1, 0)/moms(0, 0)
            bv = moms(0, 1)/moms(0, 0)
            d20 = moms(2, 0)/moms(0, 0)
            d11 = moms(1, 1)/moms(0, 0)
            d_02 = moms(0, 2)/moms(0, 0)

            c20 = d20 - bu**2._wp
            c11 = d11 - bu*bv
            c02 = d_02 - bv**2._wp

            ! First 1D quadrature (X direction)
            M1 = (/1._wp, 0._wp, c20/)
            call s_hyqmom(myrho, up, M1)
            Vf = c11*up/c20

            ! Second 1D quadrature (Y direction, conditional on X)
            mu2 = max(0._wp, c02 - sum(myrho*(Vf**2._wp)))
            M3 = (/1._wp, 0._wp, mu2/)
            call s_hyqmom(myrho3, up3, M3)

            ! Assign roots and weights for 2D quadrature
            vp21 = up3(1)
            vp22 = up3(2)
            rho21 = myrho3(1)
            rho22 = myrho3(2)

            ! Compute weights (vectorized)
            wght = moms(0, 0)*[myrho(1)*rho21, myrho(1)*rho22, myrho(2)*rho21, myrho(2)*rho22]

            ! Compute abscissas (vectorized)
            abscX = bu + [up(1), up(1), up(2), up(2)]
            abscY = bv + [Vf(1) + vp21, Vf(1) + vp22, Vf(2) + vp21, Vf(2) + vp22]

        end subroutine s_chyqmom

        pure subroutine s_hyqmom(frho, fup, fmom)
            $:GPU_ROUTINE(function_name='s_hyqmom',parallelism='[seq]', &
                & cray_inline=True)

            real(wp), dimension(2), intent(inout) :: frho, fup
            real(wp), dimension(3), intent(in) :: fmom

            real(wp) :: bu, d2, c2

            bu = fmom(2)/fmom(1)
            d2 = fmom(3)/fmom(1)
            c2 = d2 - bu**2._wp
            frho(1) = fmom(1)/2._wp; 
            frho(2) = fmom(1)/2._wp; 
            c2 = maxval((/c2, verysmall/))
            fup(1) = bu - sqrt(c2)
            fup(2) = bu + sqrt(c2)

        end subroutine s_hyqmom

        pure function f_quad(abscX, abscY, wght_in, q, r, s)
            $:GPU_ROUTINE(parallelism='[seq]')
            real(wp), dimension(nnode, nb), intent(in) :: abscX, abscY, wght_in
            real(wp), intent(in) :: q, r, s

            real(wp) :: f_quad_RV, f_quad
            integer :: i

            f_quad = 0._wp
            do i = 1, nb
                f_quad_RV = sum(wght_in(:, i)*(abscX(:, i)**q)*(abscY(:, i)**r))
                f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
            end do

        end function f_quad

        pure function f_quad2D(abscX, abscY, wght_in, pow)
            $:GPU_ROUTINE(parallelism='[seq]')
            real(wp), dimension(nnode), intent(in) :: abscX, abscY, wght_in
            real(wp), dimension(3), intent(in) :: pow

            real(wp) :: f_quad2D

            f_quad2D = sum(wght_in(:)*(abscX(:)**pow(1))*(abscY(:)**pow(2)))
        end function f_quad2D

    end subroutine s_mom_inv

end module m_qbmm
