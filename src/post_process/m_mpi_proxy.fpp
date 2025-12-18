!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

!> @brief  This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required communication goals
!!              for the post-process.
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_common

    use ieee_arithmetic

    implicit none

    !> @name Receive counts and displacement vector variables, respectively, used in
    !! enabling MPI to gather varying amounts of data from all processes to the
    !! root process
    !> @{
    integer, allocatable, dimension(:) :: recvcounts
    integer, allocatable, dimension(:) :: displs
    !> @}

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    impure subroutine s_initialize_mpi_proxy_module

#ifdef MFC_MPI

        integer :: i !< Generic loop iterator
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Allocating and configuring the receive counts and the displacement
        ! vector variables used in variable-gather communication procedures.
        ! Note that these are only needed for either multidimensional runs
        ! that utilize the Silo database file format or for 1D simulations.
        if ((format == 1 .and. n > 0) .or. n == 0) then

            allocate (recvcounts(0:num_procs - 1))
            allocate (displs(0:num_procs - 1))

            if (n == 0) then
                call MPI_GATHER(m + 1, 1, MPI_INTEGER, recvcounts(0), 1, &
                                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            elseif (proc_rank == 0) then
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

    !>  Since only processor with rank 0 is in charge of reading
        !!      and checking the consistency of the user provided inputs,
        !!      these are not available to the remaining processors. This
        !!      subroutine is then in charge of broadcasting the required
        !!      information.
    impure subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI
        integer :: i !< Generic loop iterator
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in [ 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',               &
            & 't_step_start', 't_step_stop', 't_step_save', 'weno_order',      &
            & 'model_eqns', 'num_fluids', 'bc_x%beg', 'bc_x%end', 'bc_y%beg',  &
            & 'bc_y%end', 'bc_z%beg', 'bc_z%end', 'flux_lim', 'format',        &
            & 'precision', 'fd_order', 'thermal', 'nb', 'relax_model',         &
            & 'n_start', 'num_ibs', 'muscl_order' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'cyl_coord', 'mpp_lim', 'mixture_err',                  &
            & 'alt_soundspeed', 'hypoelasticity', 'mhd', 'parallel_io',        &
            & 'rho_wrt', 'E_wrt', 'pres_wrt', 'gamma_wrt', 'sim_data',         &
            & 'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt', 'cons_vars_wrt', &
            & 'prim_vars_wrt', 'c_wrt', 'qm_wrt','schlieren_wrt','chem_wrt_T', &
            & 'bubbles_euler', 'qbmm', 'polytropic', 'polydisperse',           &
            & 'file_per_process', 'relax', 'cf_wrt', 'igr', 'liutex_wrt',      &
            & 'adv_n', 'ib', 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt',          &
            & 'surface_tension', 'hyperelasticity', 'bubbles_lagrange',        &
            & 'output_partial_domain', 'relativity', 'cont_damage', 'bc_io',   &
            & 'down_sample','fft_wrt' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        if (bubbles_lagrange) then
            #:for VAR in ['lag_header', 'lag_txt_wrt', 'lag_db_wrt', 'lag_id_wrt',         &
                & 'lag_pos_wrt', 'lag_pos_prev_wrt', 'lag_vel_wrt', 'lag_rad_wrt', &
                & 'lag_rvel_wrt', 'lag_r0_wrt', 'lag_rmax_wrt', 'lag_rmin_wrt',    &
                & 'lag_dphidt_wrt', 'lag_pres_wrt', 'lag_mv_wrt', 'lag_mg_wrt',    &
                & 'lag_betaT_wrt', 'lag_betaC_wrt', 'bc_io', 'down_sample' ]
                call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        call MPI_BCAST(flux_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(omega_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mom_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(vel_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_rho_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_rho_e_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        do i = 1, num_fluids_max
            call MPI_BCAST(fluid_pp(i)%gamma, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%pi_inf, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%cv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qvp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%G, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        ! Subgrid bubble parameters
        if (bubbles_euler .or. bubbles_lagrange) then
            #:for VAR in [ 'R0ref','p0ref','rho0ref','T0ref', &
                'ss','pv','vd','mu_l','mu_v','mu_g','gam_v','gam_g', &
                'M_v','M_g','k_v','k_g','cp_v','cp_g','R_v','R_g']
                call MPI_BCAST(bub_pp%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        #:for VAR in [ 'pref', 'rhoref', 'R0ref', 'poly_sigma', 'Web', 'Ca', &
            & 'Re_inv', 'Bx0', 'sigma', 't_save', 't_stop',   &
            & 'x_output%beg', 'x_output%end', 'y_output%beg', &
            & 'y_output%end', 'z_output%beg', 'z_output%end']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(schlieren_alpha(1), num_fluids_max, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_bcast_user_inputs

    !>  This subroutine gathers the Silo database metadata for
        !!      the spatial extents in order to boost the performance of
        !!      the multidimensional visualization.
        !!  @param spatial_extents Spatial extents for each processor's sub-domain. First dimension
        !!  corresponds to the minimum and maximum values, respectively, while
        !!  the second dimension corresponds to the processor rank.
    impure subroutine s_mpi_gather_spatial_extents(spatial_extents)

        real(wp), dimension(1:, 0:), intent(INOUT) :: spatial_extents

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Simulation is 3D
        if (p > 0) then
            if (grid_geometry == 3) then
                ! Minimum spatial extent in the r-direction
                call MPI_GATHERV(minval(y_cb), 1, mpi_p, &
                                 spatial_extents(1, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Minimum spatial extent in the theta-direction
                call MPI_GATHERV(minval(z_cb), 1, mpi_p, &
                                 spatial_extents(2, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Minimum spatial extent in the z-direction
                call MPI_GATHERV(minval(x_cb), 1, mpi_p, &
                                 spatial_extents(3, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the r-direction
                call MPI_GATHERV(maxval(y_cb), 1, mpi_p, &
                                 spatial_extents(4, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the theta-direction
                call MPI_GATHERV(maxval(z_cb), 1, mpi_p, &
                                 spatial_extents(5, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the z-direction
                call MPI_GATHERV(maxval(x_cb), 1, mpi_p, &
                                 spatial_extents(6, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)
            else
                ! Minimum spatial extent in the x-direction
                call MPI_GATHERV(minval(x_cb), 1, mpi_p, &
                                 spatial_extents(1, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Minimum spatial extent in the y-direction
                call MPI_GATHERV(minval(y_cb), 1, mpi_p, &
                                 spatial_extents(2, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Minimum spatial extent in the z-direction
                call MPI_GATHERV(minval(z_cb), 1, mpi_p, &
                                 spatial_extents(3, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the x-direction
                call MPI_GATHERV(maxval(x_cb), 1, mpi_p, &
                                 spatial_extents(4, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the y-direction
                call MPI_GATHERV(maxval(y_cb), 1, mpi_p, &
                                 spatial_extents(5, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)

                ! Maximum spatial extent in the z-direction
                call MPI_GATHERV(maxval(z_cb), 1, mpi_p, &
                                 spatial_extents(6, 0), recvcounts, 6*displs, &
                                 mpi_p, 0, MPI_COMM_WORLD, &
                                 ierr)
            end if
            ! Simulation is 2D
        elseif (n > 0) then

            ! Minimum spatial extent in the x-direction
            call MPI_GATHERV(minval(x_cb), 1, mpi_p, &
                             spatial_extents(1, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            ! Minimum spatial extent in the y-direction
            call MPI_GATHERV(minval(y_cb), 1, mpi_p, &
                             spatial_extents(2, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            ! Maximum spatial extent in the x-direction
            call MPI_GATHERV(maxval(x_cb), 1, mpi_p, &
                             spatial_extents(3, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            ! Maximum spatial extent in the y-direction
            call MPI_GATHERV(maxval(y_cb), 1, mpi_p, &
                             spatial_extents(4, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)
            ! Simulation is 1D
        else

            ! Minimum spatial extent in the x-direction
            call MPI_GATHERV(minval(x_cb), 1, mpi_p, &
                             spatial_extents(1, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            ! Maximum spatial extent in the x-direction
            call MPI_GATHERV(maxval(x_cb), 1, mpi_p, &
                             spatial_extents(2, 0), recvcounts, 4*displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)
        end if

#endif

    end subroutine s_mpi_gather_spatial_extents

    !>  This subroutine collects the sub-domain cell-boundary or
        !!      cell-center locations data from all of the processors and
        !!      puts back together the grid of the entire computational
        !!      domain on the rank 0 processor. This is only done for 1D
        !!      simulations.
    impure subroutine s_mpi_defragment_1d_grid_variable

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Silo-HDF5 database format
        if (format == 1) then

            call MPI_GATHERV(x_cc(0), m + 1, mpi_p, &
                             x_root_cc(0), recvcounts, displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            ! Binary database format
        else

            call MPI_GATHERV(x_cb(0), m + 1, mpi_p, &
                             x_root_cb(0), recvcounts, displs, &
                             mpi_p, 0, MPI_COMM_WORLD, &
                             ierr)

            if (proc_rank == 0) x_root_cb(-1) = x_cb(-1)

        end if

#endif

    end subroutine s_mpi_defragment_1d_grid_variable

    !>  This subroutine gathers the Silo database metadata for
        !!      the flow variable's extents as to boost performance of
        !!      the multidimensional visualization.
        !!  @param q_sf Flow variable defined on a single computational sub-domain
        !!  @param data_extents The flow variable extents on each of the processor's sub-domain.
        !!   First dimension of array corresponds to the former's minimum and
        !!  maximum values, respectively, while second dimension corresponds
        !!  to each processor's rank.
    impure subroutine s_mpi_gather_data_extents(q_sf, data_extents)

        real(wp), dimension(:, :, :), intent(in) :: q_sf
        real(wp), dimension(1:2, 0:num_procs - 1), intent(inout) :: data_extents

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Minimum flow variable extent
        call MPI_GATHERV(minval(q_sf), 1, mpi_p, &
                         data_extents(1, 0), recvcounts, 2*displs, &
                         mpi_p, 0, MPI_COMM_WORLD, ierr)

        ! Maximum flow variable extent
        call MPI_GATHERV(maxval(q_sf), 1, mpi_p, &
                         data_extents(2, 0), recvcounts, 2*displs, &
                         mpi_p, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_gather_data_extents

    !>  This subroutine gathers the sub-domain flow variable data
        !!      from all of the processors and puts it back together for
        !!      the entire computational domain on the rank 0 processor.
        !!      This is only done for 1D simulations.
        !!  @param q_sf Flow variable defined on a single computational sub-domain
        !!  @param q_root_sf Flow variable defined on the entire computational domain
    impure subroutine s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)

        real(wp), dimension(0:m), intent(in) :: q_sf
        real(wp), dimension(0:m), intent(inout) :: q_root_sf

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        ! Gathering the sub-domain flow variable data from all the processes
        ! and putting it back together for the entire computational domain
        ! on the process with rank 0
        call MPI_GATHERV(q_sf(0), m + 1, mpi_p, &
                         q_root_sf(0), recvcounts, displs, &
                         mpi_p, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_defragment_1d_flow_variable

    !> Deallocation procedures for the module
    impure subroutine s_finalize_mpi_proxy_module

#ifdef MFC_MPI

        ! Deallocating the receive counts and the displacement vector
        ! variables used in variable-gather communication procedures
        if ((format == 1 .and. n > 0) .or. n == 0) then
            deallocate (recvcounts)
            deallocate (displs)
        end if

#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
