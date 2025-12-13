!>
!! @file m_bubbles_EL.fpp
!! @brief Contains module m_bubbles_EL

#:include 'macros.fpp'

!> @brief This module is used to to compute the volume-averaged bubble model
module m_bubbles_EL

    use m_global_parameters             !< Definitions of the global parameters

    use m_mpi_proxy                     !< Message passing interface (MPI) module proxy

    use m_bubbles_EL_kernels            !< Definitions of the kernel functions

    use m_bubbles                       !< General bubble dynamics procedures

    use m_variables_conversion          !< State variables type conversion procedures

    use m_compile_specific

    use m_boundary_common

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_helper

    implicit none

    !(nBub)
    integer, allocatable, dimension(:, :) :: lag_id                 !< Global and local IDs
    real(wp), allocatable, dimension(:) :: bub_R0            !< Initial bubble radius
    real(wp), allocatable, dimension(:) :: Rmax_stats        !< Maximum radius
    real(wp), allocatable, dimension(:) :: Rmin_stats        !< Minimum radius
    $:GPU_DECLARE(create='[lag_id, bub_R0, Rmax_stats, Rmin_stats]')

    real(wp), allocatable, dimension(:) :: gas_mg            !< Bubble's gas mass
    real(wp), allocatable, dimension(:) :: gas_betaT         !< heatflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: gas_betaC         !< massflux model (Preston et al., 2007)
    real(wp), allocatable, dimension(:) :: bub_dphidt        !< subgrid velocity potential (Maeda & Colonius, 2018)
    $:GPU_DECLARE(create='[gas_mg, gas_betaT, gas_betaC, bub_dphidt]')

    !(nBub, 1 -> actual val or 2 -> temp val)
    real(wp), allocatable, dimension(:, :) :: gas_p          !< Pressure in the bubble
    real(wp), allocatable, dimension(:, :) :: gas_mv         !< Vapor mass in the bubble
    real(wp), allocatable, dimension(:, :) :: intfc_rad      !< Bubble radius
    real(wp), allocatable, dimension(:, :) :: intfc_vel      !< Velocity of the bubble interface
    $:GPU_DECLARE(create='[gas_p, gas_mv, intfc_rad, intfc_vel]')
    !(nBub, 1-> x or 2->y or 3 ->z, 1 -> actual or 2 -> temporal val)
    real(wp), allocatable, dimension(:, :, :) :: mtn_pos     !< Bubble's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_posPrev !< Bubble's previous position
    real(wp), allocatable, dimension(:, :, :) :: mtn_vel     !< Bubble's velocity
    real(wp), allocatable, dimension(:, :, :) :: mtn_s       !< Bubble's computational cell position in real format
    $:GPU_DECLARE(create='[mtn_pos, mtn_posPrev, mtn_vel, mtn_s]')
    !(nBub, 1-> x or 2->y or 3 ->z, time-stage)
    real(wp), allocatable, dimension(:, :) :: intfc_draddt   !< Time derivative of bubble's radius
    real(wp), allocatable, dimension(:, :) :: intfc_dveldt   !< Time derivative of bubble's interface velocity
    real(wp), allocatable, dimension(:, :) :: gas_dpdt       !< Time derivative of gas pressure
    real(wp), allocatable, dimension(:, :) :: gas_dmvdt      !< Time derivative of the vapor mass in the bubble
    real(wp), allocatable, dimension(:, :, :) :: mtn_dposdt  !< Time derivative of the bubble's position
    real(wp), allocatable, dimension(:, :, :) :: mtn_dveldt  !< Time derivative of the bubble's velocity
    $:GPU_DECLARE(create='[intfc_draddt, intfc_dveldt, gas_dpdt, gas_dmvdt, mtn_dposdt, mtn_dveldt]')

    integer, private :: lag_num_ts                                  !<  Number of time stages in the time-stepping scheme

    $:GPU_DECLARE(create='[lag_num_ts]')

    integer :: nBubs                            !< Number of bubbles in the local domain
    real(wp) :: Rmax_glb, Rmin_glb       !< Maximum and minimum bubbe size in the local domain
    !< Projection of the lagrangian particles in the Eulerian framework
    type(scalar_field), dimension(:), allocatable :: q_beta
    integer :: q_beta_idx                       !< Size of the q_beta vector field

    $:GPU_DECLARE(create='[nBubs,Rmax_glb,Rmin_glb,q_beta,q_beta_idx]')

contains

    !> Initializes the lagrangian subgrid bubble solver
        !! @param q_cons_vf Initial conservative variables
    impure subroutine s_initialize_bubbles_EL_module(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        integer :: nBubs_glb, i

        ! Setting number of time-stages for selected time-stepping scheme
        lag_num_ts = time_stepper

        ! Allocate space for the Eulerian fields needed to map the effect of the bubbles
        if (lag_params%solver_approach == 1) then
            ! One-way coupling
            q_beta_idx = 3
        elseif (lag_params%solver_approach == 2) then
            ! Two-way coupling
            q_beta_idx = 4
            if (p == 0) then
                !Subgrid noise model for 2D approximation
                q_beta_idx = 6
            end if
        else
            call s_mpi_abort('Please check the lag_params%solver_approach input')
        end if

        $:GPU_UPDATE(device='[lag_num_ts, q_beta_idx]')

        @:ALLOCATE(q_beta(1:q_beta_idx))

        do i = 1, q_beta_idx
            @:ALLOCATE(q_beta(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
        end do

        do i = 1, q_beta_idx
            @:ACC_SETUP_SFs(q_beta(i))
        end do

        ! Allocating space for lagrangian variables
        nBubs_glb = lag_params%nBubs_glb

        @:ALLOCATE(lag_id(1:nBubs_glb, 1:2))
        @:ALLOCATE(bub_R0(1:nBubs_glb))
        @:ALLOCATE(Rmax_stats(1:nBubs_glb))
        @:ALLOCATE(Rmin_stats(1:nBubs_glb))
        @:ALLOCATE(gas_mg(1:nBubs_glb))
        @:ALLOCATE(gas_betaT(1:nBubs_glb))
        @:ALLOCATE(gas_betaC(1:nBubs_glb))
        @:ALLOCATE(bub_dphidt(1:nBubs_glb))
        @:ALLOCATE(gas_p(1:nBubs_glb, 1:2))
        @:ALLOCATE(gas_mv(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_rad(1:nBubs_glb, 1:2))
        @:ALLOCATE(intfc_vel(1:nBubs_glb, 1:2))
        @:ALLOCATE(mtn_pos(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_posPrev(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_vel(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(mtn_s(1:nBubs_glb, 1:3, 1:2))
        @:ALLOCATE(intfc_draddt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(intfc_dveldt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dpdt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(gas_dmvdt(1:nBubs_glb, 1:lag_num_ts))
        @:ALLOCATE(mtn_dposdt(1:nBubs_glb, 1:3, 1:lag_num_ts))
        @:ALLOCATE(mtn_dveldt(1:nBubs_glb, 1:3, 1:lag_num_ts))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol

        ! Starting bubbles
        call s_read_input_bubbles(q_cons_vf)

    end subroutine s_initialize_bubbles_EL_module

    !> The purpose of this procedure is to obtain the initial bubbles' information
        !! @param q_cons_vf Conservative variables
    impure subroutine s_read_input_bubbles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        real(wp), dimension(8) :: inputBubble
        real(wp) :: qtime
        integer :: id, bub_id, save_count
        integer :: i, ios
        logical :: file_exist, indomain

        character(LEN=path_len + 2*name_len) :: path_D_dir !<

        ! Initialize number of particles
        bub_id = 0
        id = 0

        ! Read the input lag_bubble file or restart point
        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

        if (save_count == 0) then
            if (proc_rank == 0) print *, 'Reading lagrange bubbles input file.'
            inquire (file='input/lag_bubbles.dat', exist=file_exist)
            if (file_exist) then
                open (94, file='input/lag_bubbles.dat', form='formatted', iostat=ios)
                do while (ios == 0)
                    read (94, *, iostat=ios) (inputBubble(i), i=1, 8)
                    if (ios /= 0) cycle
                    indomain = particle_in_domain_physical(inputBubble(1:3))
                    id = id + 1
                    if (id > lag_params%nBubs_glb .and. proc_rank == 0) then
                        call s_mpi_abort("Current number of bubbles is larger than nBubs_glb")
                    end if
                    if (indomain) then
                        bub_id = bub_id + 1
                        call s_add_bubbles(inputBubble, q_cons_vf, bub_id)
                        lag_id(bub_id, 1) = id      !global ID
                        lag_id(bub_id, 2) = bub_id  !local ID
                        nBubs = bub_id              ! local number of bubbles
                    end if
                end do
                close (94)
            else
                call s_mpi_abort("Initialize the lagrange bubbles in input/lag_bubbles.dat")
            end if
        else
            if (proc_rank == 0) print *, 'Restarting lagrange bubbles at save_count: ', save_count
            call s_restart_bubbles(bub_id, save_count)
        end if

        print *, " Lagrange bubbles running, in proc", proc_rank, "number:", bub_id, "/", id

        $:GPU_UPDATE(device='[bubbles_lagrange, lag_params]')

        $:GPU_UPDATE(device='[lag_id,bub_R0,Rmax_stats,Rmin_stats,gas_mg, &
            & gas_betaT,gas_betaC,bub_dphidt,gas_p,gas_mv, &
            & intfc_rad,intfc_vel,mtn_pos,mtn_posPrev,mtn_vel, &
            & mtn_s,intfc_draddt,intfc_dveldt,gas_dpdt,gas_dmvdt, &
            & mtn_dposdt,mtn_dveldt,nBubs]')

        Rmax_glb = min(dflt_real, -dflt_real)
        Rmin_glb = max(dflt_real, -dflt_real)
        $:GPU_UPDATE(device='[Rmax_glb, Rmin_glb]')

        $:GPU_UPDATE(device='[dx,dy,dz,x_cb,x_cc,y_cb,y_cc,z_cb,z_cc]')

        !Populate temporal variables
        call s_transfer_data_to_tmp()
        call s_smear_voidfraction()

        if (lag_params%write_bubbles) call s_write_lag_particles(qtime)

        if (save_count == 0) then
            ! Create ./D directory
            write (path_D_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
            call my_inquire(path_D_dir, file_exist)
            if (.not. file_exist) call s_create_directory(trim(path_D_dir))
            call s_write_restart_lag_bubbles(save_count) ! Needed for post_processing
            call s_write_void_evol(qtime)
        end if

    end subroutine s_read_input_bubbles

    !> The purpose of this procedure is to obtain the information of the bubbles when starting fresh
        !! @param inputBubble Bubble information
        !! @param q_cons_vf Conservative variables
        !! @param bub_id Local id of the bubble
    impure subroutine s_add_bubbles(inputBubble, q_cons_vf, bub_id)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(8), intent(in) :: inputBubble
        integer, intent(in) :: bub_id
        integer :: i

        real(wp) :: pliq, volparticle, concvap, totalmass, kparticle, cpparticle
        real(wp) :: omegaN_local, PeG, PeT, rhol, pcrit, qv, gamma, pi_inf, dynP
        integer, dimension(3) :: cell
        real(wp), dimension(2) :: Re
        real(wp) :: massflag, heatflag, Re_trans, Im_trans

        massflag = 0._wp
        heatflag = 0._wp
        if (lag_params%massTransfer_model) massflag = 1._wp
        if (lag_params%heatTransfer_model) heatflag = 1._wp

        bub_R0(bub_id) = inputBubble(7)
        Rmax_stats(bub_id) = min(dflt_real, -dflt_real)
        Rmin_stats(bub_id) = max(dflt_real, -dflt_real)
        bub_dphidt(bub_id) = 0._wp
        intfc_rad(bub_id, 1) = inputBubble(7)
        intfc_vel(bub_id, 1) = inputBubble(8)
        mtn_pos(bub_id, 1:3, 1) = inputBubble(1:3)
        mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1)
        mtn_vel(bub_id, 1:3, 1) = inputBubble(4:6)

        if (cyl_coord .and. p == 0) then
            mtn_pos(bub_id, 2, 1) = sqrt(mtn_pos(bub_id, 2, 1)**2._wp + &
                                         mtn_pos(bub_id, 3, 1)**2._wp)
            !Storing azimuthal angle (-Pi to Pi)) into the third coordinate variable
            mtn_pos(bub_id, 3, 1) = atan2(inputBubble(3), inputBubble(2))
            mtn_posPrev(bub_id, 1:3, 1) = mtn_pos(bub_id, 1:3, 1)
        end if

        cell = -buff_size
        call s_locate_cell(mtn_pos(bub_id, 1:3, 1), cell, mtn_s(bub_id, 1:3, 1))

        ! Check if the bubble is located in the ghost cell of a symmetric, or wall boundary
        if ((any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) < 0) .or. &
            (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(1) > m) .or. &
            (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) < 0) .or. &
            (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(2) > n)) then
            call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric or wall boundary.")
        end if

        if (p > 0) then
            if ((any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) < 0) .or. &
                (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. cell(3) > p)) then
                call s_mpi_abort("Lagrange bubble is in the ghost cells of a symmetric or wall boundary.")
            end if
        end if

        call s_convert_to_mixture_variables(q_cons_vf, cell(1), cell(2), cell(3), &
                                            rhol, gamma, pi_inf, qv, Re)
        dynP = 0._wp
        do i = 1, num_dims
            dynP = dynP + 0.5_wp*q_cons_vf(contxe + i)%sf(cell(1), cell(2), cell(3))**2/rhol
        end do
        pliq = (q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - dynP - pi_inf)/gamma
        if (pliq < 0) print *, "Negative pressure", proc_rank, &
            q_cons_vf(E_idx)%sf(cell(1), cell(2), cell(3)), pi_inf, gamma, pliq, cell, dynP

        ! Initial particle pressure
        gas_p(bub_id, 1) = pliq + 2._wp*(1._wp/Web)/bub_R0(bub_id)
        if (.not. f_approx_equal((1._wp/Web), 0._wp)) then
            pcrit = pv - 4._wp*(1._wp/Web)/(3._wp*sqrt(3._wp*gas_p(bub_id, 1)*bub_R0(bub_id)**3._wp/(2._wp*(1._wp/Web))))
            pref = gas_p(bub_id, 1)
        else
            pcrit = 0._wp
        end if

        ! Initial particle mass
        volparticle = 4._wp/3._wp*pi*bub_R0(bub_id)**3._wp ! volume
        gas_mv(bub_id, 1) = pv*volparticle*(1._wp/(R_v*Tw))*(massflag) ! vapermass
        gas_mg(bub_id) = (gas_p(bub_id, 1) - pv*(massflag))*volparticle*(1._wp/(R_g*Tw)) ! gasmass
        if (gas_mg(bub_id) <= 0._wp) then
            call s_mpi_abort("The initial mass of gas inside the bubble is negative. Check the initial conditions.")
        end if
        totalmass = gas_mg(bub_id) + gas_mv(bub_id, 1) ! totalmass

        ! Bubble natural frequency
        concvap = gas_mv(bub_id, 1)/(gas_mv(bub_id, 1) + gas_mg(bub_id))
        omegaN_local = (3._wp*(gas_p(bub_id, 1) - pv*(massflag)) + 4._wp*(1._wp/Web)/bub_R0(bub_id))/rhol
        if (pv*(massflag) > gas_p(bub_id, 1)) then
            call s_mpi_abort("Lagrange bubble initially located in a region with pressure below the vapor pressure.")
        end if
        omegaN_local = sqrt(omegaN_local/bub_R0(bub_id)**2._wp)

        cpparticle = concvap*cp_v + (1._wp - concvap)*cp_g
        kparticle = concvap*k_vl + (1._wp - concvap)*k_gl

        ! Mass and heat transfer coefficients (based on Preston 2007)
        PeT = totalmass/volparticle*cpparticle*bub_R0(bub_id)**2._wp*omegaN_local/kparticle
        call s_transcoeff(1._wp, PeT, Re_trans, Im_trans)
        gas_betaT(bub_id) = Re_trans*(heatflag)*kparticle

        PeG = bub_R0(bub_id)**2._wp*omegaN_local/vd
        call s_transcoeff(1._wp, PeG, Re_trans, Im_trans)
        gas_betaC(bub_id) = Re_trans*(massflag)*vd

        if (gas_mg(bub_id) <= 0._wp) then
            call s_mpi_abort("Negative gas mass in the bubble, check if the bubble is in the domain.")
        end if

    end subroutine s_add_bubbles

    !> The purpose of this procedure is to obtain the information of the bubbles from a restart point.
        !! @param bub_id Local ID of the particle
        !! @param save_count File identifier
    impure subroutine s_restart_bubbles(bub_id, save_count)

        integer, intent(inout) :: bub_id, save_count

        character(LEN=path_len + 2*name_len) :: file_loc
        real(wp) :: file_time, file_dt
        integer :: file_num_procs, file_tot_part, tot_part

#ifdef MFC_MPI
        real(wp), dimension(20) :: inputvals
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: view

        integer, dimension(3) :: cell
        logical :: indomain, particle_file, file_exist

        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer :: ifile, ierr, tot_data, id
        integer :: i

        integer, dimension(:), allocatable :: proc_bubble_counts
        real(wp), dimension(1:1, 1:lag_io_vars) :: dummy
        dummy = 0._wp

        ! Construct file path
        write (file_loc, '(A,I0,A)') 'lag_bubbles_', save_count, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)

        ! Check if file exists
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort('Restart file '//trim(file_loc)//' does not exist!')
        end if

        if (.not. parallel_io) return

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            call MPI_FILE_READ(ifile, file_tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_READ(ifile, file_time, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_dt, 1, mpi_p, status, ierr)
            call MPI_FILE_READ(ifile, file_num_procs, 1, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(file_tot_part, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_time, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_dt, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(file_num_procs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        allocate (proc_bubble_counts(file_num_procs))

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip to processor counts position
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs), &
                       MPI_OFFSET_KIND)
            call MPI_FILE_SEEK(ifile, disp, MPI_SEEK_SET, ierr)
            call MPI_FILE_READ(ifile, proc_bubble_counts, file_num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BCAST(proc_bubble_counts, file_num_procs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! Set time variables from file
        mytime = file_time
        dt = file_dt

        bub_id = proc_bubble_counts(proc_rank + 1)

        start_idx_part(1) = 0
        do i = 1, proc_rank
            start_idx_part(1) = start_idx_part(1) + proc_bubble_counts(i)
        end do

        start_idx_part(2) = 0
        lsizes(1) = bub_id
        lsizes(2) = lag_io_vars

        gsizes(1) = file_tot_part
        gsizes(2) = lag_io_vars

        if (bub_id > 0) then

            allocate (MPI_IO_DATA_lag_bubbles(bub_id, 1:lag_io_vars))

            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                          MPI_ORDER_FORTRAN, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) + &
                       file_num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA_lag_bubbles, &
                                   lag_io_vars*bub_id, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)

            nBubs = bub_id

            do i = 1, bub_id
                lag_id(i, 1) = int(MPI_IO_DATA_lag_bubbles(i, 1))
                mtn_pos(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 2:4)
                mtn_posPrev(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 5:7)
                mtn_vel(i, 1:3, 1) = MPI_IO_DATA_lag_bubbles(i, 8:10)
                intfc_rad(i, 1) = MPI_IO_DATA_lag_bubbles(i, 11)
                intfc_vel(i, 1) = MPI_IO_DATA_lag_bubbles(i, 12)
                bub_R0(i) = MPI_IO_DATA_lag_bubbles(i, 13)
                Rmax_stats(i) = MPI_IO_DATA_lag_bubbles(i, 14)
                Rmin_stats(i) = MPI_IO_DATA_lag_bubbles(i, 15)
                bub_dphidt(i) = MPI_IO_DATA_lag_bubbles(i, 16)
                gas_p(i, 1) = MPI_IO_DATA_lag_bubbles(i, 17)
                gas_mv(i, 1) = MPI_IO_DATA_lag_bubbles(i, 18)
                gas_mg(i) = MPI_IO_DATA_lag_bubbles(i, 19)
                gas_betaT(i) = MPI_IO_DATA_lag_bubbles(i, 20)
                gas_betaC(i) = MPI_IO_DATA_lag_bubbles(i, 21)
                cell = -buff_size
                call s_locate_cell(mtn_pos(i, 1:3, 1), cell, mtn_s(i, 1:3, 1))
            end do

            deallocate (MPI_IO_DATA_lag_bubbles)

        else
            nBubs = 0

            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, &
                               mpi_info_int, ifile, ierr)

            ! Skip extended header
            disp = int(sizeof(file_tot_part) + 2*sizeof(file_time) + sizeof(file_num_procs) + &
                       file_num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ_ALL(ifile, dummy, 0, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
            call MPI_TYPE_FREE(view, ierr)
        end if

        if (proc_rank == 0) then
            write (*, '(A,I0,A,I0)') 'Read ', file_tot_part, ' particles from restart file at t_step = ', save_count
            write (*, '(A,E15.7,A,E15.7)') 'Restart time = ', mytime, ', dt = ', dt
        end if

        deallocate (proc_bubble_counts)
#endif

    end subroutine s_restart_bubbles

    !>  Contains the bubble dynamics subroutines.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param rhs_vf Calculated change of conservative variables
        !! @param t_step Current time step
        !! @param stage Current stage in the time-stepper algorithm
    subroutine s_compute_bubble_EL_dynamics(q_prim_vf, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: stage

        real(wp) :: myVapFlux
        real(wp) :: preterm1, term2, paux, pint, Romega, term1_fac
        real(wp) :: myR_m, mygamma_m, myPb, myMass_n, myMass_v
        real(wp) :: myR, myV, myBeta_c, myBeta_t, myR0, myPbdot, myMvdot
        real(wp) :: myPinf, aux1, aux2, myCson, myRho
        real(wp) :: gamma, pi_inf, qv
        real(wp), dimension(contxe) :: myalpha_rho, myalpha
        real(wp), dimension(2) :: Re
        integer, dimension(3) :: cell

        integer :: adap_dt_stop_max, adap_dt_stop !< Fail-safe exit if max iteration count reached
        real(wp) :: dmalf, dmntait, dmBtait, dm_bub_adv_src, dm_divu !< Dummy variables for unified subgrid bubble subroutines

        integer :: i, k, l

        call nvtxStartRange("LAGRANGE-BUBBLE-DYNAMICS")

        ! Subgrid p_inf model based on Maeda and Colonius (2018).
        if (lag_params%pressure_corrector) then
            ! Calculate velocity potentials (valid for one bubble per cell)
            $:GPU_PARALLEL_LOOP(private='[k,cell]')
            do k = 1, nBubs
                call s_get_pinf(k, q_prim_vf, 2, paux, cell, preterm1, term2, Romega)
                myR0 = bub_R0(k)
                myR = intfc_rad(k, 2)
                myV = intfc_vel(k, 2)
                myPb = gas_p(k, 2)
                pint = f_cpbw_KM(myR0, myR, myV, myPb)
                pint = pint + 0.5_wp*myV**2._wp
                if (lag_params%cluster_type == 2) then
                    bub_dphidt(k) = (paux - pint) + term2
                    ! Accounting for the potential induced by the bubble averaged over the control volume
                    ! Note that this is based on the incompressible flow assumption near the bubble.
                    term1_fac = 3._wp/2._wp*(myR*(Romega**2._wp - myR**2._wp))/(Romega**3._wp - myR**3._wp)
                    bub_dphidt(k) = bub_dphidt(k)/(1._wp - term1_fac)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Radial motion model
        adap_dt_stop_max = 0
        $:GPU_PARALLEL_LOOP(private='[k,i,myalpha_rho,myalpha,Re,cell]', &
            & reduction='[[adap_dt_stop_max]]',reductionOp='[MAX]', &
            & copy='[adap_dt_stop_max]',copyin='[stage]')
        do k = 1, nBubs
            ! Keller-Miksis model

            ! Current bubble state
            myPb = gas_p(k, 2)
            myMass_n = gas_mg(k)
            myMass_v = gas_mv(k, 2)
            myR = intfc_rad(k, 2)
            myV = intfc_vel(k, 2)
            myBeta_c = gas_betaC(k)
            myBeta_t = gas_betaT(k)
            myR0 = bub_R0(k)

            ! Vapor and heat fluxes
            call s_vflux(myR, myV, myPb, myMass_v, k, myVapFlux, myMass_n, myBeta_c, myR_m, mygamma_m)
            myPbdot = f_bpres_dot(myVapFlux, myR, myV, myPb, myMass_v, k, myBeta_t, myR_m, mygamma_m)
            myMvdot = 4._wp*pi*myR**2._wp*myVapFlux

            ! Obtaining driving pressure
            call s_get_pinf(k, q_prim_vf, 1, myPinf, cell, aux1, aux2)

            ! Obtain liquid density and computing speed of sound from pinf
            call s_compute_species_fraction(q_prim_vf, cell(1), cell(2), cell(3), myalpha_rho, myalpha)
            call s_convert_species_to_mixture_variables_acc(myRho, gamma, pi_inf, qv, myalpha, &
                                                            myalpha_rho, Re)
            call s_compute_cson_from_pinf(q_prim_vf, myPinf, cell, myRho, gamma, pi_inf, myCson)

            ! Adaptive time stepping
            adap_dt_stop = 0

            if (adap_dt) then

                call s_advance_step(myRho, myPinf, myR, myV, myR0, myPb, myPbdot, dmalf, &
                                    dmntait, dmBtait, dm_bub_adv_src, dm_divu, &
                                    k, myMass_v, myMass_n, myBeta_c, &
                                    myBeta_t, myCson, adap_dt_stop)

                ! Update bubble state
                intfc_rad(k, 1) = myR
                intfc_vel(k, 1) = myV
                gas_p(k, 1) = myPb
                gas_mv(k, 1) = myMass_v

            else

                ! Radial acceleration from bubble models
                intfc_dveldt(k, stage) = f_rddot(myRho, myPinf, myR, myV, myR0, &
                                                 myPb, myPbdot, dmalf, dmntait, dmBtait, &
                                                 dm_bub_adv_src, dm_divu, &
                                                 myCson)
                intfc_draddt(k, stage) = myV
                gas_dmvdt(k, stage) = myMvdot
                gas_dpdt(k, stage) = myPbdot

            end if

            adap_dt_stop_max = max(adap_dt_stop_max, adap_dt_stop)

        end do
        $:END_GPU_PARALLEL_LOOP()

        if (adap_dt .and. adap_dt_stop_max > 0) call s_mpi_abort("Adaptive time stepping failed to converge.")

        ! Bubbles remain in a fixed position
        $:GPU_PARALLEL_LOOP(collapse=2, private='[k,l]', copyin='[stage]')
        do k = 1, nBubs
            do l = 1, 3
                mtn_dposdt(k, l, stage) = 0._wp
                mtn_dveldt(k, l, stage) = 0._wp
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_compute_bubble_EL_dynamics

    !>  The purpose of this subroutine is to obtain the bubble source terms based on Maeda and Colonius (2018)
        !!      and add them to the RHS scalar field.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Conservative variables
        !! @param rhs_vf Time derivative of the conservative variables
    subroutine s_compute_bubbles_EL_source(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l

        if (.not. adap_dt) call s_smear_voidfraction()

        if (lag_params%solver_approach == 2) then

            ! (q / (1 - beta)) * d(beta)/dt source
            if (p == 0) then
                $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            do l = 1, E_idx
                                if (q_beta(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                    rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                            q_cons_vf(l)%sf(i, j, k)*(q_beta(2)%sf(i, j, k) + &
                                                                                      q_beta(5)%sf(i, j, k))

                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            do l = 1, E_idx
                                if (q_beta(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                    rhs_vf(l)%sf(i, j, k) = rhs_vf(l)%sf(i, j, k) + &
                                                            q_cons_vf(l)%sf(i, j, k)/q_beta(1)%sf(i, j, k)* &
                                                            q_beta(2)%sf(i, j, k)
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            do l = 1, num_dims

                call s_gradient_dir(q_prim_vf(E_idx)%sf, q_beta(3)%sf, l)

                ! (q / (1 - beta)) * d(beta)/dt source
                $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            if (q_beta(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(contxe + l)%sf(i, j, k) = rhs_vf(contxe + l)%sf(i, j, k) - &
                                                                 (1._wp - q_beta(1)%sf(i, j, k))/ &
                                                                 q_beta(1)%sf(i, j, k)* &
                                                                 q_beta(3)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                !source in energy
                $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
                do k = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(2)%beg, idwbuff(2)%end
                        do i = idwbuff(1)%beg, idwbuff(1)%end
                            q_beta(3)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)*q_prim_vf(contxe + l)%sf(i, j, k)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_gradient_dir(q_beta(3)%sf, q_beta(4)%sf, l)

                ! (beta / (1 - beta)) * d(Pu)/dl source
                $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
                do k = 0, p
                    do j = 0, n
                        do i = 0, m
                            if (q_beta(1)%sf(i, j, k) > (1._wp - lag_params%valmaxvoid)) then
                                rhs_vf(E_idx)%sf(i, j, k) = rhs_vf(E_idx)%sf(i, j, k) - &
                                                            q_beta(4)%sf(i, j, k)*(1._wp - q_beta(1)%sf(i, j, k))/ &
                                                            q_beta(1)%sf(i, j, k)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do

        end if

    end subroutine s_compute_bubbles_EL_source

    !>  This procedure computes the speed of sound from a given driving pressure
        !! @param bub_id Bubble id
        !! @param q_prim_vf Primitive variables
        !! @param pinf Driving pressure
        !! @param cell Bubble cell
        !! @param rhol Liquid density
        !! @param gamma Liquid specific heat ratio
        !! @param pi_inf Liquid stiffness
        !! @param cson Calculated speed of sound
    subroutine s_compute_cson_from_pinf(q_prim_vf, pinf, cell, rhol, gamma, pi_inf, cson)
        $:GPU_ROUTINE(function_name='s_compute_cson_from_pinf', &
            & parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(in) :: pinf, rhol, gamma, pi_inf
        integer, dimension(3), intent(in) :: cell
        real(wp), intent(out) :: cson

        real(wp) :: E, H
        real(wp), dimension(num_dims) :: vel
        integer :: i

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            vel(i) = q_prim_vf(i + contxe)%sf(cell(1), cell(2), cell(3))
        end do
        E = gamma*pinf + pi_inf + 0.5_wp*rhol*dot_product(vel, vel)
        H = (E + pinf)/rhol
        cson = sqrt((H - 0.5_wp*dot_product(vel, vel))/gamma)

    end subroutine s_compute_cson_from_pinf

    !>  The purpose of this subroutine is to smear the effect of the bubbles in the Eulerian framework
    subroutine s_smear_voidfraction()

        integer :: i, j, k, l

        call nvtxStartRange("BUBBLES-LAGRANGE-KERNELS")

        $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
        do i = 1, q_beta_idx
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_beta(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_smoothfunction(nBubs, intfc_rad, intfc_vel, &
                              mtn_s, mtn_pos, q_beta)

        !Store 1-beta
        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    q_beta(1)%sf(j, k, l) = 1._wp - q_beta(1)%sf(j, k, l)
                    ! Limiting void fraction given max value
                    q_beta(1)%sf(j, k, l) = max(q_beta(1)%sf(j, k, l), &
                                                1._wp - lag_params%valmaxvoid)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_smear_voidfraction

    !> The purpose of this procedure is obtain the bubble driving pressure p_inf
        !! @param bub_id Particle identifier
        !! @param q_prim_vf  Primitive variables
        !! @param ptype 1: p at infinity, 2: averaged P at the bubble location
        !! @param f_pinfl Driving pressure
        !! @param cell Bubble cell
        !! @param Romega Control volume radius
    subroutine s_get_pinf(bub_id, q_prim_vf, ptype, f_pinfl, cell, preterm1, term2, Romega)
        $:GPU_ROUTINE(function_name='s_get_pinf',parallelism='[seq]', &
            & cray_inline=True)

        integer, intent(in) :: bub_id, ptype
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), intent(out) :: f_pinfl
        integer, dimension(3), intent(out) :: cell
        real(wp), intent(out), optional :: preterm1, term2, Romega

        real(wp), dimension(3) :: scoord, psi
        real(wp) :: dc, vol, aux
        real(wp) :: volgas, term1, Rbeq, denom
        real(wp) :: charvol, charpres, charvol2, charpres2
        integer, dimension(3) :: cellaux
        integer :: i, j, k
        integer :: smearGrid, smearGridz
        logical :: celloutside

        scoord = mtn_s(bub_id, 1:3, 2)
        f_pinfl = 0._wp

        !< Find current bubble cell
        cell(:) = int(scoord(:))
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

        if ((lag_params%cluster_type == 1)) then
            !< Getting p_cell in terms of only the current cell by interpolation

            !< Getting the cell volulme as Omega
            if (p > 0) then
                vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            else
                if (cyl_coord) then
                    vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
                else
                    vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                end if
            end if

            !< Obtain bilinear interpolation coefficients, based on the current location of the bubble.
            psi(1) = (scoord(1) - real(cell(1)))*dx(cell(1)) + x_cb(cell(1) - 1)
            if (cell(1) == (m + buff_size)) then
                cell(1) = cell(1) - 1
                psi(1) = 1._wp
            else if (cell(1) == (-buff_size)) then
                psi(1) = 0._wp
            else
                if (psi(1) < x_cc(cell(1))) cell(1) = cell(1) - 1
                psi(1) = abs((psi(1) - x_cc(cell(1)))/(x_cc(cell(1) + 1) - x_cc(cell(1))))
            end if

            psi(2) = (scoord(2) - real(cell(2)))*dy(cell(2)) + y_cb(cell(2) - 1)
            if (cell(2) == (n + buff_size)) then
                cell(2) = cell(2) - 1
                psi(2) = 1._wp
            else if (cell(2) == (-buff_size)) then
                psi(2) = 0._wp
            else
                if (psi(2) < y_cc(cell(2))) cell(2) = cell(2) - 1
                psi(2) = abs((psi(2) - y_cc(cell(2)))/(y_cc(cell(2) + 1) - y_cc(cell(2))))
            end if

            if (p > 0) then
                psi(3) = (scoord(3) - real(cell(3)))*dz(cell(3)) + z_cb(cell(3) - 1)
                if (cell(3) == (p + buff_size)) then
                    cell(3) = cell(3) - 1
                    psi(3) = 1._wp
                else if (cell(3) == (-buff_size)) then
                    psi(3) = 0._wp
                else
                    if (psi(3) < z_cc(cell(3))) cell(3) = cell(3) - 1
                    psi(3) = abs((psi(3) - z_cc(cell(3)))/(z_cc(cell(3) + 1) - z_cc(cell(3))))
                end if
            else
                psi(3) = 0._wp
            end if

            !< Perform bilinear interpolation
            if (p == 0) then  !2D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)
            else              !3D
                f_pinfl = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))*(1._wp - psi(1))*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3))*psi(1)*(1._wp - psi(2))*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3))*psi(1)*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3))*(1._wp - psi(1))*psi(2)*(1._wp - psi(3))
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 1)*(1._wp - psi(1))*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3) + 1)*psi(1)*(1._wp - psi(2))*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2) + 1, cell(3) + 1)*psi(1)*psi(2)*psi(3)
                f_pinfl = f_pinfl + q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3) + 1)*(1._wp - psi(1))*psi(2)*psi(3)
            end if

            !R_Omega
            dc = (3._wp*vol/(4._wp*pi))**(1._wp/3._wp)

        else if (lag_params%cluster_type >= 2) then
            ! Bubble dynamic closure from Maeda and Colonius (2018)

            ! Include the cell that contains the bubble (mapCells+1+mapCells)
            smearGrid = mapCells - (-mapCells) + 1
            smearGridz = smearGrid
            if (p == 0) smearGridz = 1

            charvol = 0._wp
            charpres = 0._wp
            charvol2 = 0._wp
            charpres2 = 0._wp
            vol = 0._wp

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, smearGrid
                $:GPU_LOOP(parallelism='[seq]')
                do j = 1, smearGrid
                    $:GPU_LOOP(parallelism='[seq]')
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !< check if the current cell is outside the computational domain or not (including ghost cells)
                        celloutside = .false.
                        if (num_dims == 2) then
                            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if
                            if (cyl_coord .and. y_cc(cellaux(2)) < 0._wp) then
                                celloutside = .true.
                            end if
                            if ((cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        else
                            if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                                celloutside = .true.
                            end if

                            if ((cellaux(3) > p + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                                celloutside = .true.
                            end if
                        end if
                        if (.not. celloutside) then
                            if (cyl_coord .and. (p == 0) .and. (y_cc(cellaux(2)) < 0._wp)) then
                                celloutside = .true.
                            end if
                        end if

                        if (.not. celloutside) then
                            !< Obtaining the cell volulme
                            if (p > 0) then
                                vol = dx(cellaux(1))*dy(cellaux(2))*dz(cellaux(3))
                            else
                                if (cyl_coord) then
                                    vol = dx(cellaux(1))*dy(cellaux(2))*y_cc(cellaux(2))*2._wp*pi
                                else
                                    vol = dx(cellaux(1))*dy(cellaux(2))*lag_params%charwidth
                                end if
                            end if
                            !< Update values
                            charvol = charvol + vol
                            charpres = charpres + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3))*vol
                            charvol2 = charvol2 + vol*q_beta(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                            charpres2 = charpres2 + q_prim_vf(E_idx)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        *vol*q_beta(1)%sf(cellaux(1), cellaux(2), cellaux(3))
                        end if

                    end do
                end do
            end do

            f_pinfl = charpres2/charvol2
            vol = charvol
            dc = (3._wp*abs(vol)/(4._wp*pi))**(1._wp/3._wp)

        end if

        if (lag_params%pressure_corrector) then

            !Valid if only one bubble exists per cell
            volgas = intfc_rad(bub_id, 2)**3._wp
            denom = intfc_rad(bub_id, 2)**2._wp
            term1 = bub_dphidt(bub_id)*intfc_rad(bub_id, 2)**2._wp
            term2 = intfc_vel(bub_id, 2)*intfc_rad(bub_id, 2)**2._wp

            Rbeq = volgas**(1._wp/3._wp) !surrogate bubble radius
            aux = dc**3._wp - Rbeq**3._wp
            term2 = term2/denom
            term2 = 3._wp/2._wp*term2**2._wp*Rbeq**3._wp*(1._wp - Rbeq/dc)/aux
            preterm1 = 3._wp/2._wp*Rbeq*(dc**2._wp - Rbeq**2._wp)/(aux*denom)

            !Control volume radius
            if (ptype == 2) Romega = dc

            ! Getting p_inf
            if (ptype == 1) then
                f_pinfl = f_pinfl + preterm1*term1 + term2
            end if

        end if

    end subroutine s_get_pinf

    !>  This subroutine updates the Lagrange variables using the tvd RK time steppers.
        !!      The time derivative of the bubble variables must be stored at every stage to avoid precision errors.
        !! @param stage Current tvd RK stage
    impure subroutine s_update_lagrange_tdv_rk(stage)

        integer, intent(in) :: stage

        integer :: k

        if (time_stepper == 1) then ! 1st order TVD RK
            $:GPU_PARALLEL_LOOP(private='[k]')
            do k = 1, nBubs
                !u{1} = u{n} +  dt * RHS{n}
                intfc_rad(k, 1) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                intfc_vel(k, 1) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                gas_p(k, 1) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                gas_mv(k, 1) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_transfer_data_to_tmp()
            call s_write_void_evol(mytime)
            if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

            if (lag_params%write_bubbles) then
                $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                call s_write_lag_particles(mytime)
            end if

        elseif (time_stepper == 2) then ! 2nd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                end do
                $:END_GPU_PARALLEL_LOOP()

            elseif (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} + (1/2) * dt * (RHS{n} + RHS{1})
                    intfc_rad(k, 1) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/2._wp
                    intfc_vel(k, 1) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/2._wp
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/2._wp
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/2._wp
                    gas_p(k, 1) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/2._wp
                    gas_mv(k, 1) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/2._wp
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp()
                call s_write_void_evol(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime)
                end if

            end if

        elseif (time_stepper == 3) then ! 3rd order TVD RK
            if (stage == 1) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{1} = u{n} +  dt * RHS{n}
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*intfc_draddt(k, 1)
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*intfc_dveldt(k, 1)
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*mtn_dposdt(k, 1:3, 1)
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*mtn_dveldt(k, 1:3, 1)
                    gas_p(k, 2) = gas_p(k, 1) + dt*gas_dpdt(k, 1)
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*gas_dmvdt(k, 1)
                end do
                $:END_GPU_PARALLEL_LOOP()

            elseif (stage == 2) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{2} = u{n} + (1/4) * dt * [RHS{n} + RHS{1}]
                    intfc_rad(k, 2) = intfc_rad(k, 1) + dt*(intfc_draddt(k, 1) + intfc_draddt(k, 2))/4._wp
                    intfc_vel(k, 2) = intfc_vel(k, 1) + dt*(intfc_dveldt(k, 1) + intfc_dveldt(k, 2))/4._wp
                    mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1) + dt*(mtn_dposdt(k, 1:3, 1) + mtn_dposdt(k, 1:3, 2))/4._wp
                    mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1) + dt*(mtn_dveldt(k, 1:3, 1) + mtn_dveldt(k, 1:3, 2))/4._wp
                    gas_p(k, 2) = gas_p(k, 1) + dt*(gas_dpdt(k, 1) + gas_dpdt(k, 2))/4._wp
                    gas_mv(k, 2) = gas_mv(k, 1) + dt*(gas_dmvdt(k, 1) + gas_dmvdt(k, 2))/4._wp
                end do
                $:END_GPU_PARALLEL_LOOP()
            elseif (stage == 3) then
                $:GPU_PARALLEL_LOOP(private='[k]')
                do k = 1, nBubs
                    !u{n+1} = u{n} + (2/3) * dt * [(1/4)* RHS{n} + (1/4)* RHS{1} + RHS{2}]
                    intfc_rad(k, 1) = intfc_rad(k, 1) + (2._wp/3._wp)*dt*(intfc_draddt(k, 1)/4._wp + intfc_draddt(k, 2)/4._wp + intfc_draddt(k, 3))
                    intfc_vel(k, 1) = intfc_vel(k, 1) + (2._wp/3._wp)*dt*(intfc_dveldt(k, 1)/4._wp + intfc_dveldt(k, 2)/4._wp + intfc_dveldt(k, 3))
                    mtn_pos(k, 1:3, 1) = mtn_pos(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dposdt(k, 1:3, 1)/4._wp + mtn_dposdt(k, 1:3, 2)/4._wp + mtn_dposdt(k, 1:3, 3))
                    mtn_vel(k, 1:3, 1) = mtn_vel(k, 1:3, 1) + (2._wp/3._wp)*dt*(mtn_dveldt(k, 1:3, 1)/4._wp + mtn_dveldt(k, 1:3, 2)/4._wp + mtn_dveldt(k, 1:3, 3))
                    gas_p(k, 1) = gas_p(k, 1) + (2._wp/3._wp)*dt*(gas_dpdt(k, 1)/4._wp + gas_dpdt(k, 2)/4._wp + gas_dpdt(k, 3))
                    gas_mv(k, 1) = gas_mv(k, 1) + (2._wp/3._wp)*dt*(gas_dmvdt(k, 1)/4._wp + gas_dmvdt(k, 2)/4._wp + gas_dmvdt(k, 3))
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_transfer_data_to_tmp()
                call s_write_void_evol(mytime)
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,gas_mg,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime)
                end if

            end if

        end if

    end subroutine s_update_lagrange_tdv_rk

    !> This subroutine returns the computational coordinate of the cell for the given position.
          !! @param pos Input coordinates
          !! @param cell Computational coordinate of the cell
          !! @param scoord Calculated particle coordinates
    subroutine s_locate_cell(pos, cell, scoord)

        real(wp), dimension(3), intent(in) :: pos
        real(wp), dimension(3), intent(out) :: scoord
        integer, dimension(3), intent(inout) :: cell

        integer :: i

        do while (pos(1) < x_cb(cell(1) - 1))
            cell(1) = cell(1) - 1
        end do

        do while (pos(1) > x_cb(cell(1)))
            cell(1) = cell(1) + 1
        end do

        do while (pos(2) < y_cb(cell(2) - 1))
            cell(2) = cell(2) - 1
        end do

        do while (pos(2) > y_cb(cell(2)))
            cell(2) = cell(2) + 1
        end do

        if (p > 0) then
            do while (pos(3) < z_cb(cell(3) - 1))
                cell(3) = cell(3) - 1
            end do
            do while (pos(3) > z_cb(cell(3)))
                cell(3) = cell(3) + 1
            end do
        end if

        ! The numbering of the cell of which left boundary is the domain boundary is 0.
        ! if comp.coord of the pos is s, the real coordinate of s is
        ! (the coordinate of the left boundary of the Floor(s)-th cell)
        ! + (s-(int(s))*(cell-width).
        ! In other words,  the coordinate of the center of the cell is x_cc(cell).

        !coordinates in computational space
        scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
        scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
        scoord(3) = 0._wp
        if (p > 0) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))
        cell(:) = int(scoord(:))
        do i = 1, num_dims
            if (scoord(i) < 0._wp) cell(i) = cell(i) - 1
        end do

    end subroutine s_locate_cell

    !> This subroutine transfer data into the temporal variables.
    impure subroutine s_transfer_data_to_tmp()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]')
        do k = 1, nBubs
            gas_p(k, 2) = gas_p(k, 1)
            gas_mv(k, 2) = gas_mv(k, 1)
            intfc_rad(k, 2) = intfc_rad(k, 1)
            intfc_vel(k, 2) = intfc_vel(k, 1)
            mtn_pos(k, 1:3, 2) = mtn_pos(k, 1:3, 1)
            mtn_posPrev(k, 1:3, 2) = mtn_posPrev(k, 1:3, 1)
            mtn_vel(k, 1:3, 2) = mtn_vel(k, 1:3, 1)
            mtn_s(k, 1:3, 2) = mtn_s(k, 1:3, 1)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transfer_data_to_tmp

    !> The purpose of this procedure is to determine if the global coordinates of the bubbles
        !!      are present in the current MPI processor (including ghost cells).
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain(pos_part)

        logical :: particle_in_domain
        real(wp), dimension(3), intent(in) :: pos_part

        ! 2D
        if (p == 0 .and. cyl_coord .neqv. .true.) then
            ! Defining a virtual z-axis that has the same dimensions as y-axis
            ! defined in the input file
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < lag_params%charwidth/2._wp) .and. (pos_part(3) >= -lag_params%charwidth/2._wp))
        else
            ! cyl_coord
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (abs(pos_part(2)) < y_cb(n + buff_size)) .and. (abs(pos_part(2)) >= max(y_cb(-buff_size - 1), 0._wp)))
        end if

        ! 3D
        if (p > 0) then
            particle_in_domain = ((pos_part(1) < x_cb(m + buff_size)) .and. (pos_part(1) >= x_cb(-buff_size - 1)) .and. &
                                  (pos_part(2) < y_cb(n + buff_size)) .and. (pos_part(2) >= y_cb(-buff_size - 1)) .and. &
                                  (pos_part(3) < z_cb(p + buff_size)) .and. (pos_part(3) >= z_cb(-buff_size - 1)))
        end if

        ! For symmetric and wall boundary condition
        if (any(bc_x%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) >= x_cb(-1)))
        end if
        if (any(bc_x%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
            particle_in_domain = (particle_in_domain .and. (pos_part(1) < x_cb(m)))
        end if
        if (any(bc_y%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) >= y_cb(-1)))
        end if
        if (any(bc_y%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/)) .and. (.not. cyl_coord)) then
            particle_in_domain = (particle_in_domain .and. (pos_part(2) < y_cb(n)))
        end if
        if (p > 0) then
            if (any(bc_z%beg == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) >= z_cb(-1)))
            end if
            if (any(bc_z%end == (/BC_REFLECTIVE, BC_CHAR_SLIP_WALL, BC_SLIP_WALL, BC_NO_SLIP_WALL/))) then
                particle_in_domain = (particle_in_domain .and. (pos_part(3) < z_cb(p)))
            end if
        end if

    end function particle_in_domain

    !> The purpose of this procedure is to determine if the lagrangian bubble is located in the
        !!       physical domain. The ghost cells are not part of the physical domain.
        !! @param pos_part Spatial coordinates of the bubble
    function particle_in_domain_physical(pos_part)

        logical :: particle_in_domain_physical
        real(wp), dimension(3), intent(in) :: pos_part

        particle_in_domain_physical = ((pos_part(1) < x_cb(m)) .and. (pos_part(1) >= x_cb(-1)) .and. &
                                       (pos_part(2) < y_cb(n)) .and. (pos_part(2) >= y_cb(-1)))

        if (p > 0) then
            particle_in_domain_physical = (particle_in_domain_physical .and. (pos_part(3) < z_cb(p)) .and. (pos_part(3) >= z_cb(-1)))
        end if

    end function particle_in_domain_physical

    !> The purpose of this procedure is to calculate the gradient of a scalar field along the x, y and z directions
        !!      following a second-order central difference considering uneven widths
        !! @param q Input scalar field
        !! @param dq Output gradient of q
        !! @param dir Gradient spatial direction
    subroutine s_gradient_dir(q, dq, dir)

        real(stp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:), intent(inout) :: q, dq
        integer, intent(in) :: dir

        integer :: i, j, k

        if (dir == 1) then
            ! Gradient in x dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dx(i + 1) - dx(i - 1)) &
                                      + q(i + 1, j, k)*(dx(i) + dx(i - 1)) &
                                      - q(i - 1, j, k)*(dx(i) + dx(i + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dx(i) + dx(i - 1))*(dx(i) + dx(i + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 2) then
            ! Gradient in y dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dy(j + 1) - dy(j - 1)) &
                                      + q(i, j + 1, k)*(dy(j) + dy(j - 1)) &
                                      - q(i, j - 1, k)*(dy(j) + dy(j + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dy(j) + dy(j - 1))*(dy(j) + dy(j + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (dir == 3) then
            ! Gradient in z dir.
            $:GPU_PARALLEL_LOOP(private='[i,j,k]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        dq(i, j, k) = q(i, j, k)*(dz(k + 1) - dz(k - 1)) &
                                      + q(i, j, k + 1)*(dz(k) + dz(k - 1)) &
                                      - q(i, j, k - 1)*(dz(k) + dz(k + 1))
                        dq(i, j, k) = dq(i, j, k)/ &
                                      ((dz(k) + dz(k - 1))*(dz(k) + dz(k + 1)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_gradient_dir

    !> Subroutine that writes on each time step the changes of the lagrangian bubbles.
        !!  @param q_time Current time
    impure subroutine s_write_lag_particles(qtime)

        real(wp), intent(in) :: qtime
        integer :: k

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        character(LEN=25) :: FMT

        write (file_loc, '(A,I0,A)') 'lag_bubble_evol_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (precision == 1) then
            FMT = "(A16,A14,8A16)"
        else
            FMT = "(A24,A14,8A24)"
        end if

        if (.not. file_exist) then
            open (11, FILE=trim(file_loc), FORM='formatted', position='rewind')
            write (11, FMT) 'currentTime', 'particleID', 'x', 'y', 'z', &
                'coreVaporMass', 'coreVaporConcentration', 'radius', 'interfaceVelocity', &
                'corePressure'
        else
            open (11, FILE=trim(file_loc), FORM='formatted', position='append')
        end if

        if (precision == 1) then
            FMT = "(F16.8,I14,8F16.8)"
        else
            FMT = "(F24.16,I14,8F24.16)"
        end if

        ! Cycle through list
        do k = 1, nBubs
            write (11, FMT) &
                qtime, &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                gas_mv(k, 1), &
                gas_mv(k, 1)/(gas_mv(k, 1) + gas_mg(k)), &
                intfc_rad(k, 1), &
                intfc_vel(k, 1), &
                gas_p(k, 1)
        end do

        close (11)

    end subroutine s_write_lag_particles

    !>  Subroutine that writes some useful statistics related to the volume fraction
            !!       of the particles (void fraction) in the computatioational domain
            !!       on each time step.
            !!  @param q_time Current time
    impure subroutine s_write_void_evol(qtime)

        real(wp), intent(in) :: qtime
        real(wp) :: volcell, voltot
        real(wp) :: lag_void_max, lag_void_avg, lag_vol
        real(wp) :: void_max_glb, void_avg_glb, vol_glb

        integer :: i, j, k

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        if (proc_rank == 0) then
            write (file_loc, '(A)') 'voidfraction.dat'
            file_loc = trim(case_dir)//'/D/'//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (.not. file_exist) then
                open (12, FILE=trim(file_loc), FORM='formatted', position='rewind')
                !write (12, *) 'currentTime, averageVoidFraction, ', &
                !    'maximumVoidFraction, totalParticlesVolume'
                !write (12, *) 'The averageVoidFraction value does ', &
                !    'not reflect the real void fraction in the cloud since the ', &
                !    'cells which do not have bubbles are not accounted'
            else
                open (12, FILE=trim(file_loc), FORM='formatted', position='append')
            end if
        end if

        lag_void_max = 0._wp
        lag_void_avg = 0._wp
        lag_vol = 0._wp
        $:GPU_PARALLEL_LOOP(private='[volcell]', collapse=3, reduction='[[lag_vol, lag_void_avg], [lag_void_max]]', reductionOp='[+, MAX]', copy='[lag_vol, lag_void_avg, lag_void_max]')
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    lag_void_max = max(lag_void_max, 1._wp - q_beta(1)%sf(i, j, k))
                    call s_get_char_vol(i, j, k, volcell)
                    if ((1._wp - q_beta(1)%sf(i, j, k)) > 5.0d-11) then
                        lag_void_avg = lag_void_avg + (1._wp - q_beta(1)%sf(i, j, k))*volcell
                        lag_vol = lag_vol + volcell
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

#ifdef MFC_MPI
        if (num_procs > 1) then
            call s_mpi_allreduce_max(lag_void_max, void_max_glb)
            lag_void_max = void_max_glb
            call s_mpi_allreduce_sum(lag_vol, vol_glb)
            lag_vol = vol_glb
            call s_mpi_allreduce_sum(lag_void_avg, void_avg_glb)
            lag_void_avg = void_avg_glb
        end if
#endif
        voltot = lag_void_avg
        ! This voidavg value does not reflect the real void fraction in the cloud
        ! since the cell which does not have bubbles are not accounted
        if (lag_vol > 0._wp) lag_void_avg = lag_void_avg/lag_vol

        if (proc_rank == 0) then
            write (12, '(6X,4e24.8)') &
                qtime, &
                lag_void_avg, &
                lag_void_max, &
                voltot
            close (12)
        end if

    end subroutine s_write_void_evol

    !>  Subroutine that writes the restarting files for the particles in the lagrangian solver.
        !!  @param t_step Current time step
    impure subroutine s_write_restart_lag_bubbles(t_step)

        ! Generic string used to store the address of a particular file
        integer, intent(in) :: t_step

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        integer :: bub_id, tot_part
        integer :: i, k

#ifdef MFC_MPI
        ! For Parallel I/O
        integer :: ifile, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer :: view
        integer, dimension(2) :: gsizes, lsizes, start_idx_part
        integer, dimension(num_procs) :: part_order, part_ord_mpi
        integer, dimension(num_procs) :: proc_bubble_counts
        real(wp), dimension(1:1, 1:lag_io_vars) :: dummy
        dummy = 0._wp

        bub_id = 0._wp
        if (nBubs /= 0) then
            do k = 1, nBubs
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    bub_id = bub_id + 1
                end if
            end do
        end if

        if (.not. parallel_io) return

        lsizes(1) = bub_id
        lsizes(2) = lag_io_vars

        ! Total number of particles
        call MPI_ALLREDUCE(bub_id, tot_part, 1, MPI_integer, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        call MPI_ALLGATHER(bub_id, 1, MPI_INTEGER, proc_bubble_counts, 1, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierr)

        ! Calculate starting index for this processor's particles
        call MPI_EXSCAN(lsizes(1), start_idx_part(1), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (proc_rank == 0) start_idx_part(1) = 0
        start_idx_part(2) = 0

        gsizes(1) = tot_part
        gsizes(2) = lag_io_vars

        write (file_loc, '(A,I0,A)') 'lag_bubbles_', t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)

        ! Clean up existing file
        if (proc_rank == 0) then
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (proc_rank == 0) then
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Write header using MPI I/O for consistency
            call MPI_FILE_WRITE(ifile, tot_part, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, mytime, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, dt, 1, mpi_p, status, ierr)
            call MPI_FILE_WRITE(ifile, num_procs, 1, MPI_INTEGER, status, ierr)
            call MPI_FILE_WRITE(ifile, proc_bubble_counts, num_procs, MPI_INTEGER, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (bub_id > 0) then
            allocate (MPI_IO_DATA_lag_bubbles(max(1, bub_id), 1:lag_io_vars))

            i = 1
            do k = 1, nBubs
                if (particle_in_domain_physical(mtn_pos(k, 1:3, 1))) then
                    MPI_IO_DATA_lag_bubbles(i, 1) = real(lag_id(k, 1))
                    MPI_IO_DATA_lag_bubbles(i, 2:4) = mtn_pos(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 5:7) = mtn_posPrev(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 8:10) = mtn_vel(k, 1:3, 1)
                    MPI_IO_DATA_lag_bubbles(i, 11) = intfc_rad(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 12) = intfc_vel(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 13) = bub_R0(k)
                    MPI_IO_DATA_lag_bubbles(i, 14) = Rmax_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 15) = Rmin_stats(k)
                    MPI_IO_DATA_lag_bubbles(i, 16) = bub_dphidt(k)
                    MPI_IO_DATA_lag_bubbles(i, 17) = gas_p(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 18) = gas_mv(k, 1)
                    MPI_IO_DATA_lag_bubbles(i, 19) = gas_mg(k)
                    MPI_IO_DATA_lag_bubbles(i, 20) = gas_betaT(k)
                    MPI_IO_DATA_lag_bubbles(i, 21) = gas_betaC(k)
                    i = i + 1
                end if
            end do

            call MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, start_idx_part, &
                                          MPI_ORDER_FORTRAN, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Skip header (written by rank 0)
            disp = int(sizeof(tot_part) + 2*sizeof(mytime) + sizeof(num_procs) + &
                       num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA_lag_bubbles, &
                                    lag_io_vars*bub_id, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)

            deallocate (MPI_IO_DATA_lag_bubbles)

        else
            call MPI_TYPE_CONTIGUOUS(0, mpi_p, view, ierr)
            call MPI_TYPE_COMMIT(view, ierr)

            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Skip header (written by rank 0)
            disp = int(sizeof(tot_part) + 2*sizeof(mytime) + sizeof(num_procs) + &
                       num_procs*sizeof(proc_bubble_counts(1)), MPI_OFFSET_KIND)
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)

            call MPI_FILE_WRITE_ALL(ifile, dummy, 0, mpi_p, status, ierr)

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

#endif

    end subroutine s_write_restart_lag_bubbles

    !>  This procedure calculates the maximum and minimum radius of each bubble.
    subroutine s_calculate_lag_bubble_stats()

        integer :: k

        $:GPU_PARALLEL_LOOP(private='[k]', reduction='[[Rmax_glb], [Rmin_glb]]', &
            & reductionOp='[MAX, MIN]', copy='[Rmax_glb,Rmin_glb]')
        do k = 1, nBubs
            Rmax_glb = max(Rmax_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmin_glb = min(Rmin_glb, intfc_rad(k, 1)/bub_R0(k))
            Rmax_stats(k) = max(Rmax_stats(k), intfc_rad(k, 1)/bub_R0(k))
            Rmin_stats(k) = min(Rmin_stats(k), intfc_rad(k, 1)/bub_R0(k))
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_calculate_lag_bubble_stats

    !>  Subroutine that writes the maximum and minimum radius of each bubble.
    impure subroutine s_write_lag_bubble_stats()

        integer :: k
        character(LEN=path_len + 2*name_len) :: file_loc

        character(len=20) :: FMT

        write (file_loc, '(A,I0,A)') 'stats_lag_bubbles_', proc_rank, '.dat'
        file_loc = trim(case_dir)//'/D/'//trim(file_loc)

        $:GPU_UPDATE(host='[Rmax_glb,Rmin_glb]')

        if (precision == 1) then
            FMT = "(A10,A14,5A16)"
        else
            FMT = "(A10,A14,5A24)"
        end if

        open (13, FILE=trim(file_loc), FORM='formatted', position='rewind')
        write (13, FMT) 'proc_rank', 'particleID', 'x', 'y', 'z', 'Rmax_glb', 'Rmin_glb'

        if (precision == 1) then
            FMT = "(I10,I14,5F16.8)"
        else
            FMT = "(I10,I14,5F24.16)"
        end if

        do k = 1, nBubs
            write (13, FMT) &
                proc_rank, &
                lag_id(k, 1), &
                mtn_pos(k, 1, 1), &
                mtn_pos(k, 2, 1), &
                mtn_pos(k, 3, 1), &
                Rmax_stats(k), &
                Rmin_stats(k)
        end do

        close (13)

    end subroutine s_write_lag_bubble_stats

    !> The purpose of this subroutine is to remove one specific particle if dt is too small.
          !! @param bub_id Particle id
    impure subroutine s_remove_lag_bubble(bub_id)

        integer, intent(in) :: bub_id

        integer :: i

        $:GPU_LOOP(parallelism='[seq]')
        do i = bub_id, nBubs - 1
            lag_id(i, 1) = lag_id(i + 1, 1)
            bub_R0(i) = bub_R0(i + 1)
            Rmax_stats(i) = Rmax_stats(i + 1)
            Rmin_stats(i) = Rmin_stats(i + 1)
            gas_mg(i) = gas_mg(i + 1)
            gas_betaT(i) = gas_betaT(i + 1)
            gas_betaC(i) = gas_betaC(i + 1)
            bub_dphidt(i) = bub_dphidt(i + 1)
            gas_p(i, 1:2) = gas_p(i + 1, 1:2)
            gas_mv(i, 1:2) = gas_mv(i + 1, 1:2)
            intfc_rad(i, 1:2) = intfc_rad(i + 1, 1:2)
            intfc_vel(i, 1:2) = intfc_vel(i + 1, 1:2)
            mtn_pos(i, 1:3, 1:2) = mtn_pos(i + 1, 1:3, 1:2)
            mtn_posPrev(i, 1:3, 1:2) = mtn_posPrev(i + 1, 1:3, 1:2)
            mtn_vel(i, 1:3, 1:2) = mtn_vel(i + 1, 1:3, 1:2)
            mtn_s(i, 1:3, 1:2) = mtn_s(i + 1, 1:3, 1:2)
            intfc_draddt(i, 1:lag_num_ts) = intfc_draddt(i + 1, 1:lag_num_ts)
            intfc_dveldt(i, 1:lag_num_ts) = intfc_dveldt(i + 1, 1:lag_num_ts)
            gas_dpdt(i, 1:lag_num_ts) = gas_dpdt(i + 1, 1:lag_num_ts)
            gas_dmvdt(i, 1:lag_num_ts) = gas_dmvdt(i + 1, 1:lag_num_ts)
        end do

        nBubs = nBubs - 1
        $:GPU_UPDATE(device='[nBubs]')

    end subroutine s_remove_lag_bubble

    !> The purpose of this subroutine is to deallocate variables
    impure subroutine s_finalize_lagrangian_solver()

        integer :: i

        do i = 1, q_beta_idx
            @:DEALLOCATE(q_beta(i)%sf)
        end do
        @:DEALLOCATE(q_beta)

        !Deallocating space
        @:DEALLOCATE(lag_id)
        @:DEALLOCATE(bub_R0)
        @:DEALLOCATE(Rmax_stats)
        @:DEALLOCATE(Rmin_stats)
        @:DEALLOCATE(gas_mg)
        @:DEALLOCATE(gas_betaT)
        @:DEALLOCATE(gas_betaC)
        @:DEALLOCATE(bub_dphidt)
        @:DEALLOCATE(gas_p)
        @:DEALLOCATE(gas_mv)
        @:DEALLOCATE(intfc_rad)
        @:DEALLOCATE(intfc_vel)
        @:DEALLOCATE(mtn_pos)
        @:DEALLOCATE(mtn_posPrev)
        @:DEALLOCATE(mtn_vel)
        @:DEALLOCATE(mtn_s)
        @:DEALLOCATE(intfc_draddt)
        @:DEALLOCATE(intfc_dveldt)
        @:DEALLOCATE(gas_dpdt)
        @:DEALLOCATE(gas_dmvdt)
        @:DEALLOCATE(mtn_dposdt)
        @:DEALLOCATE(mtn_dveldt)

    end subroutine s_finalize_lagrangian_solver

end module m_bubbles_EL
