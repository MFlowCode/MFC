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
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!                   4) 4th Order RK
!!                   5) 5th Order RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), private, allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(bounds_info) :: ix, iy, iz !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        elseif (time_stepper == 4) then
            num_ts = 3
        else
            num_ts = 6
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
        allocate (q_cons_ts(1:num_ts))

        do i = 1, num_ts
            allocate (q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                allocate (q_cons_ts(i)%vf(j)%sf(ix%beg:ix%end, &
                                                iy%beg:iy%end, &
                                                iz%beg:iz%end))
            end do
        end do

        ! Allocating the cell-average primitive ts variables
        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            allocate (q_prim_ts(0:3))

            do i = 0, 3
                allocate (q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    allocate (q_prim_ts(i)%vf(j)%sf(ix%beg:ix%end, &
                                                    iy%beg:iy%end, &
                                                    iz%beg:iz%end))
                end do
            end do
        end if

        ! Allocating the cell-average primitive variables
        allocate (q_prim_vf(1:sys_size))

        do i = mom_idx%beg, E_idx
            allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
        end do

        if (bubbles) then
            do i = bub_idx%beg, sys_size
                allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
            end do
        end if


        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
            end do
        end if

        ! Allocating the cell-average RHS variables
        allocate (rhs_vf(1:sys_size))

        do i = 1, sys_size
            allocate (rhs_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 1 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        if (adv_alphan) then
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            end do
        else
            do i = adv_idx%beg, sys_size
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            end do
        end if

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
        if (DEBUG) print *, 'got rhs'

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if
        if (DEBUG) print *, 'wrote runtime info'

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        if (adv_alphan) then
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => null()
            end do
        else
            do i = adv_idx%beg, sys_size ! adv_idx%end
                q_prim_vf(i)%sf => null()
            end do
        end if
        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 2 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        ! ==================================================================

        ! Stage 2 of 2 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + dt*rhs_vf(i)%sf)/2d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i, j !< Generic loop iterator

        ! Stage 1 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (3d0*q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + dt*rhs_vf(i)%sf)/4d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*dt*rhs_vf(i)%sf)/3d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> Adaptive SSP RK23 time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_23_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)) :: relerr, absval, tmp
        real(kind(0d0)) :: dtmin, dtmax

        integer :: i, j !< Generic loop iterator

        ! Stage 1 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) call s_write_run_time_information(q_prim_vf, t_step)
        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) &
            call s_time_step_cycling(t_step)

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        ! RK2 estimate
        do i = 1, sys_size
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = ( &
                                                   q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                                                   + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                                                   + dt*rhs_vf(i)%sf)/2d0
        end do

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)

        ! Stage 2 of RK3
        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (3d0*q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + dt*rhs_vf(i)%sf)/4d0
        end do

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*dt*rhs_vf(i)%sf)/3d0
        end do

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        ! ==================================================================

        ! Approximate error =================================================
        ! err = (q_cons_ts(1)%vf(i)%sf - q_cons_ts(3)%vf(i)%sf) / &
        !     q_cons_ts(1)%vf(i)

        ! PRINT*, '          '
        ! DO i = 1,sys_size
        !     PRINT*, 'MAXVAL', i,  MAXVAL( ABS( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p)  ), 1 )
        ! END DO

        ! DO i = 1,sys_size
        !     PRINT*, 'ABSERR', i,  MAXVAL( ABS( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) - &
        !         q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p)  ), 1 )
        ! END DO

        relerr = 0d0
        do i = 1, sys_size
            absval = maxval(abs(q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p)))
            if (absval >= 1d-10) then
                relerr = max(relerr, maxval(abs( &
                                            (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) - &
                                             q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p))/ &
                                            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p))) &
                             )
            end if
        end do

        if (num_procs > 1) then
            tmp = relerr
            call s_mpi_allreduce_max(tmp, relerr)
        end if

        dtmin = 0.002d0/2d0
        dtmax = 0.002d0*2d0

        dt = dt*min(max(sqrt(t_tol/(2d0*relerr)), 0.3d0), 2d0)
        dt = max(min(dtmax, dt), dtmin)
        ! dt = 0.0015d0

        if (proc_rank == 0) print *, 'RELERR:', relerr
        if (proc_rank == 0) print *, 'dt/dt0:', dt/dt0
        ! IF (proc_rank==0) PRINT*, '---t/T:', mytime/finaltime

        ! ==================================================================

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_23_order_tvd_rk ! ------------------------------------

    !> 4th order RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_4th_order_rk(t_step) ! ------------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 4 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/2d0
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/6d0
        end do

        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
        end if
        ! ==================================================================

        ! Stage 2 of 4 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/2d0
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/3d0
        end do

        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
        end if
        ! ==================================================================

        ! Stage 3 of 4 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/3d0
        end do

        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
        end if
        ! ==================================================================

        ! Stage 4 of 4 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf/6d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_4th_order_rk ! ----------------------------------------

    !> 5th order RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_5th_order_rk(t_step) ! ------------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 2d-1*dt*rhs_vf(i)%sf
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 75d-3*dt*rhs_vf(i)%sf
            q_cons_ts(4)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 3d-1*dt*rhs_vf(i)%sf
            q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                - (11d0/54d0)*dt*rhs_vf(i)%sf
            q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (1631d0/55296d0)*dt*rhs_vf(i)%sf
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (37d0/378d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if
        ! ==================================================================

        ! Stage 2 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(3)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 225d-3*dt*rhs_vf(i)%sf
            q_cons_ts(4)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(4)%vf(i)%sf(0:m, 0:n, 0:p) &
                - 9d-1*dt*rhs_vf(i)%sf
            q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 25d-1*dt*rhs_vf(i)%sf
            q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (175d0/512d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
        end if
        ! ==================================================================

        ! Stage 3 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(3)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(3)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(3)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(4)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(4)%vf(i)%sf(0:m, 0:n, 0:p) &
                + 12d-1*dt*rhs_vf(i)%sf
            q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) &
                - (7d1/27d0)*dt*rhs_vf(i)%sf
            q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (575d0/13824d0)*dt*rhs_vf(i)%sf
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (25d1/621d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if
        ! ==================================================================

        ! Stage 4 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(4)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(4)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(4)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(5)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (35d0/27d0)*dt*rhs_vf(i)%sf
            q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (44275d0/110592d0)*dt*rhs_vf(i)%sf
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (125d0/594d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) then
            call s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if
        ! ==================================================================

        ! Stage 5 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(5)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(5)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(5)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(6)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (253d0/4096d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
        ! ==================================================================

        ! Stage 6 of 6 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(6)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(6)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(6)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + (512d0/1771d0)*dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_5th_order_rk ! ----------------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

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
                deallocate (q_cons_ts(i)%vf(j)%sf)
            end do

            deallocate (q_cons_ts(i)%vf)

        end do

        deallocate (q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    deallocate (q_prim_ts(i)%vf(j)%sf)
                end do
                deallocate (q_prim_ts(i)%vf)
            end do
            deallocate (q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = mom_idx%beg, E_idx
            deallocate (q_prim_vf(i)%sf)
        end do
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                deallocate (q_prim_vf(i)%sf)
            end do
        end if

        deallocate (q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            deallocate (rhs_vf(i)%sf)
        end do

        deallocate (rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
