# 1 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

# 1 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/common/include/macros.fpp" 1
# 6 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/common/include/macros.fpp"


# 25 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/common/include/macros.fpp"



# 45 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/common/include/macros.fpp"
# 6 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp" 2

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_fftw

    use m_nvtx
    ! ==========================================================================

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

!$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(int_bounds_info) :: ix, iy, iz !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Setting the indical bounds in the x-, y- and z-directions
        ix%beg = -buff_size; ix%end = m + buff_size

        if (n > 0) then

            iy%beg = -buff_size; iy%end = n + buff_size

            if (p > 0) then
                iz%beg = -buff_size; iz%end = p + buff_size
            else
                iz%beg = 0; iz%end = 0
            end if

        else

            iy%beg = 0; iy%end = 0
            iz%beg = 0; iz%end = 0

        end if

        ! Allocating the cell-average conservative variables

# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_cons_ts(1:num_ts))
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_cons_ts(1:num_ts))
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 92 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        do i = 1, num_ts

# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_cons_ts(i)%vf(1:sys_size))
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_cons_ts(i)%vf(1:sys_size))
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 95 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        end do

        do i = 1, num_ts
            do j = 1, sys_size

# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_cons_ts(i)%vf(j)%sf(ix%beg:ix%end,                                                 iy%beg:iy%end,                      &
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    &                           iz%beg:iz%end))
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_cons_ts(i)%vf(j)%sf(ix%beg:ix%end,                                                 iy%beg:iy%end,                                                 iz%beg:iz%end))
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 100 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 103 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
            end do
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then

# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_ts(0:3))
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_ts(0:3))
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 108 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


            do i = 0, 3

# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_ts(i)%vf(1:sys_size))
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_ts(i)%vf(1:sys_size))
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 111 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do

            do i = 0, 3
                do j = 1, sys_size

# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_ts(i)%vf(j)%sf(ix%beg:ix%end,                                                     iy%beg:iy%end,                  &
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    &                                   iz%beg:iz%end))
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_ts(i)%vf(j)%sf(ix%beg:ix%end,                                                     iy%beg:iy%end,                                                     iz%beg:iz%end))
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 116 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 119 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
                end do
            end do
        end if

        ! Allocating the cell-average primitive variables

# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_vf(1:sys_size))
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_vf(1:sys_size))
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 124 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        
        do i = 1, adv_idx%end

# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_vf(i)%sf(ix%beg:ix%end,                                       iy%beg:iy%end,                                      &
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    & iz%beg:iz%end))
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,                                       iy%beg:iy%end,                                       iz%beg:iz%end))
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 127 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 130 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end

# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                  &
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    &         iz%beg:iz%end))
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                           iz%beg:iz%end))
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 134 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 137 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
            end do
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end

# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                  &
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    &         iz%beg:iz%end))
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                           iz%beg:iz%end))
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 143 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 146 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end

# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                  &
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
    &         iz%beg:iz%end))
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,                                           iy%beg:iy%end,                                           iz%beg:iz%end))
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 151 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 154 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
            end do
        end if

        ! Allocating the cell-average RHS variables

# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(rhs_vf(1:sys_size))
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(rhs_vf(1:sys_size))
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 158 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        do i = 1, sys_size

# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
allocate(rhs_vf(i)%sf(0:m, 0:n, 0:p))
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc enter data create(rhs_vf(i)%sf(0:m, 0:n, 0:p))
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  ALLOCATE (src/common/macros.fpp) ====
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ================================================
# 161 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        end do

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !print *, q_cons_ts(1)%vf(cont_idx%beg)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(E_idx)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(adv_idx%end)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(mom_idx%beg)%sf(102,0,0)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 3 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 3 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2d0*dt*rhs_vf(i)%sf(j, k, l))/3d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
!$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size

# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_cons_ts(i)%vf(j)%sf)
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_cons_ts(i)%vf(j)%sf)
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 485 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do


# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_cons_ts(i)%vf)
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_cons_ts(i)%vf)
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 488 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        end do


# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_cons_ts)
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_cons_ts)
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 492 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size

# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_ts(i)%vf(j)%sf)
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_ts(i)%vf(j)%sf)
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 498 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

                end do

# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_ts(i)%vf)
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_ts(i)%vf)
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 500 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do

# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_ts)
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_ts)
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 502 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end

# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_vf(i)%sf)
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_vf(i)%sf)
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 507 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end

# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_vf(i)%sf)
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_vf(i)%sf)
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 512 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end

# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_vf(i)%sf)
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_vf(i)%sf)
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 518 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end

# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_vf(i)%sf)
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_vf(i)%sf)
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 524 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

            end do
        end if


# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(q_prim_vf)
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(q_prim_vf)
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 528 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size

# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(rhs_vf(i)%sf)
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(rhs_vf(i)%sf)
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 532 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

        end do


# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
deallocate(rhs_vf)
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
!$acc exit data delete(rhs_vf)
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"

# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"
! ==================================================
# 535 "/Users/benwilfong/Documents/software/MFC-Wilfong/src/simulation/m_time_steppers.fpp"


        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
