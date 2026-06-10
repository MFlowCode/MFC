!>
!! @file
!! @brief Contains module m_mpi_proxy

!> @brief MPI gather and scatter operations for distributing post-process grid and flow-variable data
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use ieee_arithmetic
    use m_constants, only: format_silo

    implicit none

    !> @name Receive counts and displacement vector variables, respectively, used in enabling MPI to gather varying amounts of data
    !! from all processes to the root process
    !> @{
    integer, allocatable, dimension(:) :: recvcounts
    integer, allocatable, dimension(:) :: displs
    !> @}

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_mpi_proxy_module

#ifdef MFC_MPI
        integer :: i     !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Allocating and configuring the receive counts and the displacement vector variables used in variable-gather communication
        ! procedures. Note that these are only needed for either multidimensional runs that utilize the Silo database file format or
        ! for 1D simulations.

        if ((format == format_silo .and. n > 0) .or. n == 0) then
            allocate (recvcounts(0:num_procs - 1))
            allocate (displs(0:num_procs - 1))

            if (n == 0) then
                call MPI_GATHER(m + 1, 1, MPI_INTEGER, recvcounts(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            else if (proc_rank == 0) then
                recvcounts = 1
            end if

            if (proc_rank == 0) then
                displs(0) = 0

                do i = 1, num_procs - 1
                    displs(i) = displs(i - 1) + recvcounts(i - 1)
                end do
            end if
        end if
#endif

    end subroutine s_initialize_mpi_proxy_module

    !> Since only processor with rank 0 is in charge of reading and checking the consistency of the user provided inputs, these are
    !! not available to the remaining processors. This subroutine is then in charge of broadcasting the required information.
    impure subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI
        integer :: i     !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        ! Generated: case_dir, namelist scalars (INT/LOG/REAL), array dims (including
        !            chem_wrt_Y - latent bug fix), fluid_pp loop, bub_pp guard
        #:include 'generated_bcast.fpp'

        ! manual: m_glb, n_glb, p_glb (computed in s_read_input_file, not namelist-bound)
        call MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x/y/z member broadcasts (struct members not in NAMELIST_VARS)
        #:for VAR in [ 'bc_x%beg', 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'bc_x%isothermal_in', 'bc_y%isothermal_in', 'bc_z%isothermal_in',   &
            & 'bc_x%isothermal_out', 'bc_y%isothermal_out', 'bc_z%isothermal_out']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        ! manual: cfl_dt (runtime-computed logical), bc_io (BC-file existence)
        call MPI_BCAST(cfl_dt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bc_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: output domain and Twall bc members (struct members not in NAMELIST_VARS)
        #:for VAR in [ 'x_output%beg', 'x_output%end', 'y_output%beg',    &
            & 'y_output%end', 'z_output%beg', 'z_output%end',              &
            & 'bc_x%Twall_in', 'bc_x%Twall_out', 'bc_y%Twall_in',         &
            & 'bc_y%Twall_out','bc_z%Twall_in', 'bc_z%Twall_out']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
#endif

    end subroutine s_mpi_bcast_user_inputs

    !> Gather spatial extents from all ranks for Silo database metadata
    impure subroutine s_mpi_gather_spatial_extents(spatial_extents)

        real(wp), dimension(1:,0:), intent(inout) :: spatial_extents

#ifdef MFC_MPI
        integer  :: ierr  !< Generic flag used to identify and report MPI errors
        real(wp) :: ext_temp(0:num_procs - 1)

        ! Simulation is 3D

        if (p > 0) then
            if (grid_geometry == 3) then
                ! Minimum spatial extent in the r-direction
                call MPI_GATHERV(minval(y_cb), 1, mpi_p, spatial_extents(1, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Minimum spatial extent in the theta-direction
                call MPI_GATHERV(minval(z_cb), 1, mpi_p, spatial_extents(2, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Minimum spatial extent in the z-direction
                call MPI_GATHERV(minval(x_cb), 1, mpi_p, spatial_extents(3, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the r-direction
                call MPI_GATHERV(maxval(y_cb), 1, mpi_p, spatial_extents(4, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the theta-direction
                call MPI_GATHERV(maxval(z_cb), 1, mpi_p, spatial_extents(5, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the z-direction
                call MPI_GATHERV(maxval(x_cb), 1, mpi_p, spatial_extents(6, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)
            else
                ! Minimum spatial extent in the x-direction
                call MPI_GATHERV(minval(x_cb), 1, mpi_p, spatial_extents(1, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Minimum spatial extent in the y-direction
                call MPI_GATHERV(minval(y_cb), 1, mpi_p, spatial_extents(2, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Minimum spatial extent in the z-direction
                call MPI_GATHERV(minval(z_cb), 1, mpi_p, spatial_extents(3, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the x-direction
                call MPI_GATHERV(maxval(x_cb), 1, mpi_p, spatial_extents(4, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the y-direction
                call MPI_GATHERV(maxval(y_cb), 1, mpi_p, spatial_extents(5, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)

                ! Maximum spatial extent in the z-direction
                call MPI_GATHERV(maxval(z_cb), 1, mpi_p, spatial_extents(6, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, &
                                 & ierr)
            end if
            ! Simulation is 2D
        else if (n > 0) then
            ! Minimum spatial extent in the x-direction
            call MPI_GATHERV(minval(x_cb), 1, mpi_p, spatial_extents(1, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Minimum spatial extent in the y-direction
            call MPI_GATHERV(minval(y_cb), 1, mpi_p, spatial_extents(2, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the x-direction
            call MPI_GATHERV(maxval(x_cb), 1, mpi_p, spatial_extents(3, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the y-direction
            call MPI_GATHERV(maxval(y_cb), 1, mpi_p, spatial_extents(4, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
            ! Simulation is 1D
        else
            ! For 1D, recvcounts/displs are sized for grid defragmentation (m+1 per rank), not for scalar gathers. Use MPI_GATHER
            ! instead.

            ! Minimum spatial extent in the x-direction
            call MPI_GATHER(minval(x_cb), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) spatial_extents(1,:) = ext_temp

            ! Maximum spatial extent in the x-direction
            call MPI_GATHER(maxval(x_cb), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) spatial_extents(2,:) = ext_temp
        end if
#endif

    end subroutine s_mpi_gather_spatial_extents

    !> Collect the sub-domain cell-boundary or cell-center location data from all processors and put back together the grid of the
    !! entire computational domain on the rank 0 processor. This is only done for 1D simulations.
    impure subroutine s_mpi_defragment_1d_grid_variable

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Silo-HDF5 database format

        if (format == format_silo) then
            call MPI_GATHERV(x_cc(0), m + 1, mpi_p, x_root_cc(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Binary database format
        else
            call MPI_GATHERV(x_cb(0), m + 1, mpi_p, x_root_cb(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            if (proc_rank == 0) x_root_cb(-1) = x_cb(-1)
        end if
#endif

    end subroutine s_mpi_defragment_1d_grid_variable

    !> Gather the Silo database metadata for the flow variable's extents to boost performance of the multidimensional visualization.
    !! @param q_sf Flow variable on a single computational sub-domain
    impure subroutine s_mpi_gather_data_extents(q_sf, data_extents)

        real(wp), dimension(:,:,:), intent(in)                  :: q_sf
        real(wp), dimension(1:2,0:num_procs - 1), intent(inout) :: data_extents

#ifdef MFC_MPI
        integer  :: ierr  !< Generic flag used to identify and report MPI errors
        real(wp) :: ext_temp(0:num_procs - 1)

        if (n > 0) then
            ! Multi-D: recvcounts = 1, so strided MPI_GATHERV works correctly Minimum flow variable extent
            call MPI_GATHERV(minval(q_sf), 1, mpi_p, data_extents(1, 0), recvcounts, 2*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum flow variable extent
            call MPI_GATHERV(maxval(q_sf), 1, mpi_p, data_extents(2, 0), recvcounts, 2*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
        else
            ! 1D: recvcounts/displs are sized for grid defragmentation (m+1 per rank), not for scalar gathers. Use MPI_GATHER
            ! instead.

            ! Minimum flow variable extent
            call MPI_GATHER(minval(q_sf), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) data_extents(1,:) = ext_temp

            ! Maximum flow variable extent
            call MPI_GATHER(maxval(q_sf), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) data_extents(2,:) = ext_temp
        end if
#endif

    end subroutine s_mpi_gather_data_extents

    !> Gather the sub-domain flow variable data from all processors and reassemble it for the entire computational domain on the
    !! rank 0 processor. This is only done for 1D simulations.
    !! @param q_sf Flow variable on a single computational sub-domain
    !! @param q_root_sf Flow variable on the entire computational domain
    impure subroutine s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)

        real(wp), dimension(0:m), intent(in)    :: q_sf
        real(wp), dimension(0:m), intent(inout) :: q_root_sf

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Gathering the sub-domain flow variable data from all the processes and putting it back together for the entire
        ! computational domain on the process with rank 0

        call MPI_GATHERV(q_sf(0), m + 1, mpi_p, q_root_sf(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_defragment_1d_flow_variable

    !> Deallocation procedures for the module
    impure subroutine s_finalize_mpi_proxy_module

#ifdef MFC_MPI
        ! Deallocating the receive counts and the displacement vector variables used in variable-gather communication procedures
        if ((format == format_silo .and. n > 0) .or. n == 0) then
            deallocate (recvcounts)
            deallocate (displs)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
