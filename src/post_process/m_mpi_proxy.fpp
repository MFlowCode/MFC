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

    use m_helper_basic, only: s_update_cell_bounds

    implicit none

    !> @name Buffers of the conservative variables received/sent from/to neighboring
    !! processors. Note that these variables are structured as vectors rather
    !! than arrays.
    !> @{
    real(wp), allocatable, dimension(:) :: q_cons_buffer_in
    real(wp), allocatable, dimension(:) :: q_cons_buffer_out
    !> @}

    !> @name Receive counts and displacement vector variables, respectively, used in
    !! enabling MPI to gather varying amounts of data from all processes to the
    !! root process
    !> @{
    integer, allocatable, dimension(:) :: recvcounts
    integer, allocatable, dimension(:) :: displs
    !> @}

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: err_code, ierr
    !> @}

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_mpi_proxy_module

#ifdef MFC_MPI

        integer :: i !< Generic loop iterator

        ! Allocating vectorized buffer regions of conservative variables.
        ! The length of buffer vectors are set according to the size of the
        ! largest buffer region in the sub-domain.
        if (buff_size > 0) then

            ! Simulation is at least 2D
            if (n > 0) then

                ! Simulation is 3D
                if (p > 0) then

                    allocate (q_cons_buffer_in(0:buff_size* &
                                               sys_size* &
                                               (m + 2*buff_size + 1)* &
                                               (n + 2*buff_size + 1)* &
                                               (p + 2*buff_size + 1)/ &
                                               (cells_bounds%mnp_min &
                                                + 2*buff_size + 1) - 1))
                    allocate (q_cons_buffer_out(0:buff_size* &
                                                sys_size* &
                                                (m + 2*buff_size + 1)* &
                                                (n + 2*buff_size + 1)* &
                                                (p + 2*buff_size + 1)/ &
                                                (cells_bounds%mnp_min &
                                                 + 2*buff_size + 1) - 1))

                    ! Simulation is 2D
                else

                    allocate (q_cons_buffer_in(0:buff_size* &
                                               sys_size* &
                                               (cells_bounds%mn_max &
                                                + 2*buff_size + 1) - 1))
                    allocate (q_cons_buffer_out(0:buff_size* &
                                                sys_size* &
                                                (cells_bounds%mn_max &
                                                 + 2*buff_size + 1) - 1))

                end if

                ! Simulation is 1D
            else

                allocate (q_cons_buffer_in(0:buff_size*sys_size - 1))
                allocate (q_cons_buffer_out(0:buff_size*sys_size - 1))

            end if

            ! Initially zeroing out the vectorized buffer region variables
            ! to avoid possible underflow from any unused allocated memory
            q_cons_buffer_in = 0._wp
            q_cons_buffer_out = 0._wp

        end if

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
    subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI
        integer :: i !< Generic loop iterator

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in [ 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',               &
            & 't_step_start', 't_step_stop', 't_step_save', 'weno_order',      &
            & 'model_eqns', 'num_fluids', 'bc_x%beg', 'bc_x%end', 'bc_y%beg',  &
            & 'bc_y%end', 'bc_z%beg', 'bc_z%end', 'flux_lim', 'format',        &
            & 'precision', 'fd_order', 'thermal', 'nb', 'relax_model',         &
            & 'n_start', 'num_ibs' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'cyl_coord', 'mpp_lim', 'mixture_err',                  &
            & 'alt_soundspeed', 'hypoelasticity', 'mhd', 'parallel_io',        &
            & 'rho_wrt', 'E_wrt', 'pres_wrt', 'gamma_wrt', 'sim_data',         &
            & 'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt', 'cons_vars_wrt', &
            & 'prim_vars_wrt', 'c_wrt', 'qm_wrt','schlieren_wrt','chem_wrt_T', &
            & 'bubbles_euler', 'qbmm', 'polytropic', 'polydisperse',           &
            & 'file_per_process', 'relax', 'cf_wrt',                           &
            & 'adv_n', 'ib', 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt',          &
            & 'surface_tension', 'hyperelasticity', 'bubbles_lagrange',        &
            & 'output_partial_domain', 'relativity', 'cont_damage' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        call MPI_BCAST(flux_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(omega_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mom_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(vel_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_rho_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        do i = 1, num_fluids_max
            call MPI_BCAST(fluid_pp(i)%gamma, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%pi_inf, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%cv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qvp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%G, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        #:for VAR in [ 'pref', 'rhoref', 'R0ref', 'poly_sigma', 'Web', 'Ca', &
            & 'Re_inv', 'Bx0', 'sigma', 't_save', 't_stop',   &
            & 'x_output%beg', 'x_output%end', 'y_output%beg', &
            & 'y_output%end', 'z_output%beg', 'z_output%end' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(schlieren_alpha(1), num_fluids_max, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_bcast_user_inputs

    !>  This subroutine takes care of efficiently distributing
        !!      the computational domain among the available processors
        !!      as well as recomputing some of the global parameters so
        !!      that they reflect the configuration of sub-domain that
        !!      is overseen by the local processor.
    subroutine s_mpi_decompose_computational_domain

#ifdef MFC_MPI

        ! # of processors in the x-, y- and z-coordinate directions
        integer :: num_procs_x, num_procs_y, num_procs_z

        ! Temporary # of processors in x-, y- and z-coordinate directions
        ! used during the processor factorization optimization procedure
        real(wp) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z

        ! Processor factorization (fct) minimization parameter
        real(wp) :: fct_min

        ! Cartesian processor topology communicator
        integer :: MPI_COMM_CART

        ! Number of remaining cells for a particular coordinate direction
        ! after the bulk has evenly been distributed among the available
        ! processors for that coordinate direction
        integer :: rem_cells

        ! Generic loop iterators
        integer :: i, j

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        ! Performing the computational domain decomposition. The procedure
        ! is optimized by ensuring that each processor contains a close to
        ! equivalent piece of the computational domain. Note that explicit
        ! type-casting is omitted here for code legibility purposes.

        ! Generating 3D Cartesian Processor Topology

        if (n > 0) then

            if (p > 0) then

                if (cyl_coord .and. p > 0) then
                    ! Implement pencil processor blocking if using cylindrical coordinates so
                    ! that all cells in azimuthal direction are stored on a single processor.
                    ! This is necessary for efficient application of Fourier filter near axis.

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = num_procs
                    num_procs_z = 1
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                         - (n + 1)/tmp_num_procs_y)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i

                            if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                               - (n + 1)/tmp_num_procs_y) &
                                .and. &
                                (n + 1)/tmp_num_procs_y &
                                >= &
                                num_stcls_min*weno_order) then

                                num_procs_x = i
                                num_procs_y = num_procs/i
                                fct_min = abs((m + 1)/tmp_num_procs_x &
                                              - (n + 1)/tmp_num_procs_y)
                                ierr = 0

                            end if

                        end if

                    end do

                else

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = 1
                    num_procs_z = num_procs
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                         - (n + 1)/tmp_num_procs_y) &
                              + 10._wp*abs((n + 1)/tmp_num_procs_y &
                                           - (p + 1)/tmp_num_procs_z)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            do j = 1, (num_procs/i)

                                if (mod(num_procs/i, j) == 0 &
                                    .and. &
                                    (n + 1)/j >= num_stcls_min*weno_order) then

                                    tmp_num_procs_x = i
                                    tmp_num_procs_y = j
                                    tmp_num_procs_z = num_procs/(i*j)

                                    if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                                       - (n + 1)/tmp_num_procs_y) &
                                        + abs((n + 1)/tmp_num_procs_y &
                                              - (p + 1)/tmp_num_procs_z) &
                                        .and. &
                                        (p + 1)/tmp_num_procs_z &
                                        >= &
                                        num_stcls_min*weno_order) &
                                        then

                                        num_procs_x = i
                                        num_procs_y = j
                                        num_procs_z = num_procs/(i*j)
                                        fct_min = abs((m + 1)/tmp_num_procs_x &
                                                      - (n + 1)/tmp_num_procs_y) &
                                                  + abs((n + 1)/tmp_num_procs_y &
                                                        - (p + 1)/tmp_num_procs_z)
                                        ierr = 0

                                    end if

                                end if

                            end do

                        end if

                    end do

                end if

                ! Checking whether the decomposition of the computational
                ! domain was successful
                if (proc_rank == 0 .and. ierr == -1) then
                    print '(A)', 'Unable to decompose computational '// &
                        'domain for selected number of '// &
                        'processors. Exiting.'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                end if

                ! Creating a new communicator using Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/num_procs_x, &
                                                          num_procs_y, num_procs_z/), &
                                     (/.true., .true., .true./), &
                                     .false., MPI_COMM_CART, ierr)

                ! Finding corresponding Cartesian coordinates of the local
                ! processor rank in newly declared cartesian communicator
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, &
                                     proc_coords, ierr)

                ! END: Generating 3D Cartesian Processor Topology

                ! Sub-domain Global Parameters in z-direction

                ! Number of remaining cells after majority is distributed
                rem_cells = mod(p + 1, num_procs_z)

                ! Optimal number of cells per processor
                p = (p + 1)/num_procs_z - 1

                ! Distributing any remaining cells
                do i = 1, rem_cells
                    if (proc_coords(3) == i - 1) then
                        p = p + 1
                        exit
                    end if
                end do

                ! Boundary condition at the beginning
                if (proc_coords(3) > 0 .or. bc_z%beg == BC_PERIODIC) then
                    proc_coords(3) = proc_coords(3) - 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%beg, ierr)
                    proc_coords(3) = proc_coords(3) + 1
                end if

                ! Ghost zone at the beginning
                if (proc_coords(3) > 0 .and. format == 1) then
                    offset_z%beg = 2
                else
                    offset_z%beg = 0
                end if

                ! Boundary condition at the end
                if (proc_coords(3) < num_procs_z - 1 .or. bc_z%end == BC_PERIODIC) then
                    proc_coords(3) = proc_coords(3) + 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%end, ierr)
                    proc_coords(3) = proc_coords(3) - 1
                end if

                ! Ghost zone at the end
                if (proc_coords(3) < num_procs_z - 1 .and. format == 1) then
                    offset_z%end = 2
                else
                    offset_z%end = 0
                end if

                if (parallel_io) then
                    if (proc_coords(3) < rem_cells) then
                        start_idx(3) = (p + 1)*proc_coords(3)
                    else
                        start_idx(3) = (p + 1)*proc_coords(3) + rem_cells
                    end if
                end if

                ! Generating 2D Cartesian Processor Topology

            else

                ! Initial values of the processor factorization optimization
                num_procs_x = 1
                num_procs_y = num_procs
                ierr = -1

                ! Computing minimization variable for these initial values
                tmp_num_procs_x = num_procs_x
                tmp_num_procs_y = num_procs_y
                fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x &
                                     - (n + 1)/tmp_num_procs_y)

                ! Searching for optimal computational domain distribution
                do i = 1, num_procs

                    if (mod(num_procs, i) == 0 &
                        .and. &
                        (m + 1)/i >= num_stcls_min*weno_order) then

                        tmp_num_procs_x = i
                        tmp_num_procs_y = num_procs/i

                        if (fct_min >= abs((m + 1)/tmp_num_procs_x &
                                           - (n + 1)/tmp_num_procs_y) &
                            .and. &
                            (n + 1)/tmp_num_procs_y &
                            >= &
                            num_stcls_min*weno_order) then

                            num_procs_x = i
                            num_procs_y = num_procs/i
                            fct_min = abs((m + 1)/tmp_num_procs_x &
                                          - (n + 1)/tmp_num_procs_y)
                            ierr = 0

                        end if

                    end if

                end do

                ! Checking whether the decomposition of the computational
                ! domain was successful
                if (proc_rank == 0 .and. ierr == -1) then
                    print '(A)', 'Unable to decompose computational '// &
                        'domain for selected number of '// &
                        'processors. Exiting.'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                end if

                ! Creating a new communicator using Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                          num_procs_y/), (/.true., &
                                                                           .true./), .false., MPI_COMM_CART, &
                                     ierr)

                ! Finding corresponding Cartesian coordinates of the local
                ! processor rank in newly declared cartesian communicator
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                     proc_coords, ierr)

            end if

            ! END: Generating 2D Cartesian Processor Topology

            ! Sub-domain Global Parameters in y-direction

            ! Number of remaining cells after majority has been distributed
            rem_cells = mod(n + 1, num_procs_y)

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing any remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1
                    exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(2) > 0 .or. bc_y%beg == BC_PERIODIC) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%beg, ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Ghost zone at the beginning
            if (proc_coords(2) > 0 .and. format == 1) then
                offset_y%beg = 2
            else
                offset_y%beg = 0
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. bc_y%end == BC_PERIODIC) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%end, ierr)
                proc_coords(2) = proc_coords(2) - 1
            end if

            ! Ghost zone at the end
            if (proc_coords(2) < num_procs_y - 1 .and. format == 1) then
                offset_y%end = 2
            else
                offset_y%end = 0
            end if

            if (parallel_io) then
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            end if

            ! Generating 1D Cartesian Processor Topology

        else

            ! Number of processors in the coordinate direction is equal to
            ! the total number of processors available
            num_procs_x = num_procs

            ! Number of cells in undecomposed computational domain needed
            ! for sub-domain reassembly during formatted data output
            m_root = m

            ! Creating a new communicator using Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the corresponding Cartesian coordinates of the local
            ! processor rank in the newly declared cartesian communicator
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if

        ! Sub-domain Global Parameters in x-direction

        ! Number of remaining cells after majority has been distributed
        rem_cells = mod(m + 1, num_procs_x)

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing any remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1
                exit
            end if
        end do

        call s_update_cell_bounds(cells_bounds, m, n, p)

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. bc_x%beg == BC_PERIODIC) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Ghost zone at the beginning
        if (proc_coords(1) > 0 .and. format == 1 .and. n > 0) then
            offset_x%beg = 2
        else
            offset_x%beg = 0
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. bc_x%end == BC_PERIODIC) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

        ! Ghost zone at the end
        if (proc_coords(1) < num_procs_x - 1 .and. format == 1 .and. n > 0) then
            offset_x%end = 2
        else
            offset_x%end = 0
        end if

        if (parallel_io) then
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if

#endif

    end subroutine s_mpi_decompose_computational_domain

    !>  Communicates the buffer regions associated with the grid
        !!      variables with processors in charge of the neighboring
        !!      sub-domains. Note that only cell-width spacings feature
        !!      buffer regions so that no information relating to the
        !!      cell-boundary locations is communicated.
        !!  @param pbc_loc Processor boundary condition (PBC) location
        !!  @param sweep_coord Coordinate direction normal to the processor boundary
    subroutine s_mpi_sendrecv_grid_vars_buffer_regions(pbc_loc, sweep_coord)

        character(LEN=3), intent(in) :: pbc_loc
        character, intent(in) :: sweep_coord

#ifdef MFC_MPI

        ! Communications in the x-direction

        if (sweep_coord == 'x') then

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_x%end >= 0) then

                    ! Sending/receiving the data to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV(dx(m - buff_size + 1), buff_size, &
                                      mpi_p, bc_x%end, 0, &
                                      dx(-buff_size), buff_size, &
                                      mpi_p, bc_x%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_x%beg/bc_x%beg
                    call MPI_SENDRECV(dx(0), buff_size, &
                                      mpi_p, bc_x%beg, 1, &
                                      dx(-buff_size), buff_size, &
                                      mpi_p, bc_x%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_x%beg >= 0) then

                    ! Sending/receiving the data to/from bc_x%beg/bc_x%end
                    call MPI_SENDRECV(dx(0), buff_size, &
                                      mpi_p, bc_x%beg, 1, &
                                      dx(m + 1), buff_size, &
                                      mpi_p, bc_x%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_x%end/bc_x%end
                    call MPI_SENDRECV(dx(m - buff_size + 1), buff_size, &
                                      mpi_p, bc_x%end, 0, &
                                      dx(m + 1), buff_size, &
                                      mpi_p, bc_x%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            end if

            ! END: Communications in the x-direction

            ! Communications in the y-direction

        elseif (sweep_coord == 'y') then

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_y%end >= 0) then

                    ! Sending/receiving the data to/from bc_y%end/bc_y%beg
                    call MPI_SENDRECV(dy(n - buff_size + 1), buff_size, &
                                      mpi_p, bc_y%end, 0, &
                                      dy(-buff_size), buff_size, &
                                      mpi_p, bc_y%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_y%beg/bc_y%beg
                    call MPI_SENDRECV(dy(0), buff_size, &
                                      mpi_p, bc_y%beg, 1, &
                                      dy(-buff_size), buff_size, &
                                      mpi_p, bc_y%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_y%beg >= 0) then

                    ! Sending/receiving the data to/from bc_y%beg/bc_y%end
                    call MPI_SENDRECV(dy(0), buff_size, &
                                      mpi_p, bc_y%beg, 1, &
                                      dy(n + 1), buff_size, &
                                      mpi_p, bc_y%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_y%end/bc_y%end
                    call MPI_SENDRECV(dy(n - buff_size + 1), buff_size, &
                                      mpi_p, bc_y%end, 0, &
                                      dy(n + 1), buff_size, &
                                      mpi_p, bc_y%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            end if

            ! END: Communications in the y-direction

            ! Communications in the z-direction

        else

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_z%end >= 0) then

                    ! Sending/receiving the data to/from bc_z%end/bc_z%beg
                    call MPI_SENDRECV(dz(p - buff_size + 1), buff_size, &
                                      mpi_p, bc_z%end, 0, &
                                      dz(-buff_size), buff_size, &
                                      mpi_p, bc_z%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_z%beg/bc_z%beg
                    call MPI_SENDRECV(dz(0), buff_size, &
                                      mpi_p, bc_z%beg, 1, &
                                      dz(-buff_size), buff_size, &
                                      mpi_p, bc_z%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_z%beg >= 0) then

                    ! Sending/receiving the data to/from bc_z%beg/bc_z%end
                    call MPI_SENDRECV(dz(0), buff_size, &
                                      mpi_p, bc_z%beg, 1, &
                                      dz(p + 1), buff_size, &
                                      mpi_p, bc_z%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Sending/receiving the data to/from bc_z%end/bc_z%end
                    call MPI_SENDRECV(dz(p - buff_size + 1), buff_size, &
                                      mpi_p, bc_z%end, 0, &
                                      dz(p + 1), buff_size, &
                                      mpi_p, bc_z%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

            end if

        end if

        ! END: Communications in the z-direction

#endif

    end subroutine s_mpi_sendrecv_grid_vars_buffer_regions

    !>  Communicates buffer regions associated with conservative
        !!      variables with processors in charge of the neighboring
        !!      sub-domains
        !!  @param q_cons_vf Conservative variables
        !!  @param pbc_loc Processor boundary condition (PBC) location
        !!  @param sweep_coord Coordinate direction normal to the processor boundary
        !!  @param q_particle Projection of the lagrangian particles in the Eulerian framework
    subroutine s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, pbc_loc, &
                                                       sweep_coord, q_particle)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf

        character(LEN=3), intent(in) :: pbc_loc

        character, intent(in) :: sweep_coord

        type(scalar_field), &
            intent(inout), optional :: q_particle

#ifdef MFC_MPI

        integer :: i, j, k, l, r !< Generic loop iterators

        ! Communications in the x-direction

        if (sweep_coord == 'x') then

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_x%end >= 0) then

                    ! Packing the data to be sent to bc_x%end
                    do l = 0, p
                        do k = 0, n
                            do j = m - buff_size + 1, m
                                do i = 1, sys_size
                                    r = sys_size*(j - m + buff_size - 1) &
                                        + sys_size*buff_size*k + (i - 1) &
                                        + sys_size*buff_size*(n + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%end, 0, &
                                      q_cons_buffer_in(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Packing the data to be sent to bc_x%beg
                    do l = 0, p
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size*j &
                                        + sys_size*buff_size*k &
                                        + sys_size*buff_size*(n + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_x%beg/bc_x%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%beg, 1, &
                                      q_cons_buffer_in(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data received from bc_x%beg
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size, -1
                            do i = 1, sys_size
                                r = sys_size*(j + buff_size) &
                                    + sys_size*buff_size*k + (i - 1) &
                                    + sys_size*buff_size*(n + 1)*l
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_x%beg >= 0) then

                    ! Packing the data to be sent to bc_x%beg
                    do l = 0, p
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, sys_size
                                    r = (i - 1) + sys_size*j &
                                        + sys_size*buff_size*k &
                                        + sys_size*buff_size*(n + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    call MPI_SENDRECV(q_cons_buffer_out(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%beg, 1, &
                                      q_cons_buffer_in(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Packing the data to be sent to bc_x%end
                    do l = 0, p
                        do k = 0, n
                            do j = m - buff_size + 1, m
                                do i = 1, sys_size
                                    r = sys_size*(j - m + buff_size - 1) &
                                        + sys_size*buff_size*k + (i - 1) &
                                        + sys_size*buff_size*(n + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    call MPI_SENDRECV(q_cons_buffer_out(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%end, 0, &
                                      q_cons_buffer_in(0), &
                                      buff_size*sys_size*(n + 1)*(p + 1), &
                                      mpi_p, bc_x%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data received from bc_x%end
                do l = 0, p
                    do k = 0, n
                        do j = m + 1, m + buff_size
                            do i = 1, sys_size
                                r = (i - 1) + sys_size*(j - m - 1) &
                                    + sys_size*buff_size*k &
                                    + sys_size*buff_size*(n + 1)*l
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            end if

            ! END: Communications in the x-direction

            ! Communications in the y-direction

        elseif (sweep_coord == 'y') then

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_y%end >= 0) then

                    ! Packing the data to be sent to bc_y%end
                    do l = 0, p
                        do k = n - buff_size + 1, n
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k - n + buff_size - 1) + (i - 1) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        buff_size*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_y%end/bc_y%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (p + 1), mpi_p, &
                                      bc_y%end, 0, q_cons_buffer_in(0), &
                                      buff_size*sys_size* &
                                      (m + 2*buff_size + 1)*(p + 1), &
                                      mpi_p, bc_y%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Packing the data to be sent to bc_y%beg
                    do l = 0, p
                        do k = 0, buff_size - 1
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)*k &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        buff_size*l + (i - 1)
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_y%beg/bc_y%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (p + 1), mpi_p, &
                                      bc_y%beg, 1, q_cons_buffer_in(0), &
                                      buff_size*sys_size* &
                                      (m + 2*buff_size + 1)*(p + 1), &
                                      mpi_p, bc_y%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data received from bc_y%beg
                do l = 0, p
                    do k = -buff_size, -1
                        do j = -buff_size, m + buff_size
                            do i = 1, sys_size
                                r = (i - 1) + sys_size*(j + buff_size) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (k + buff_size) + sys_size* &
                                    (m + 2*buff_size + 1)*buff_size*l
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_y%beg >= 0) then

                    ! Packing the data to be sent to bc_y%beg
                    do l = 0, p
                        do k = 0, buff_size - 1
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)*k &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        buff_size*l + (i - 1)
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_y%beg/bc_y%end
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (p + 1), mpi_p, &
                                      bc_y%beg, 1, q_cons_buffer_in(0), &
                                      buff_size*sys_size* &
                                      (m + 2*buff_size + 1)*(p + 1), &
                                      mpi_p, bc_y%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Packing the data to be sent to bc_y%end
                    do l = 0, p
                        do k = n - buff_size + 1, n
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k - n + buff_size - 1) + (i - 1) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        buff_size*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_y%end/bc_y%end
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (p + 1), mpi_p, &
                                      bc_y%end, 0, q_cons_buffer_in(0), &
                                      buff_size*sys_size* &
                                      (m + 2*buff_size + 1)*(p + 1), &
                                      mpi_p, bc_y%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data received form bc_y%end
                do l = 0, p
                    do k = n + 1, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, sys_size
                                r = (i - 1) + sys_size*(j + buff_size) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (k - n - 1) + sys_size* &
                                    (m + 2*buff_size + 1)*buff_size*l
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            end if

            ! END: Communications in the y-direction

            ! Communications in the z-direction

        else

            if (pbc_loc == 'beg') then    ! Buffer region at the beginning

                ! PBC at both ends of the sub-domain
                if (bc_z%end >= 0) then

                    ! Packing the data to be sent to bc_z%end
                    do l = p - buff_size + 1, p
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k + buff_size) + sys_size* &
                                        (m + 2*buff_size + 1)* &
                                        (n + 2*buff_size + 1)* &
                                        (l - p + buff_size - 1) + (i - 1)
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_z%end/bc_z%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%end, 0, &
                                      q_cons_buffer_in(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at beginning of the sub-domain
                else

                    ! Packing the data to be sent to bc_z%beg
                    do l = 0, buff_size - 1
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k + buff_size) + (i - 1) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (n + 2*buff_size + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_z%beg/bc_z%beg
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%beg, 1, &
                                      q_cons_buffer_in(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%beg, 0, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data from bc_z%beg
                do l = -buff_size, -1
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, sys_size
                                r = sys_size*(j + buff_size) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (k + buff_size) + (i - 1) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (n + 2*buff_size + 1)*(l + buff_size)
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            else                         ! Buffer region at the end

                ! PBC at both ends of the sub-domain
                if (bc_z%beg >= 0) then

                    ! Packing the data to be sent to bc_z%beg
                    do l = 0, buff_size - 1
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k + buff_size) + (i - 1) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (n + 2*buff_size + 1)*l
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_z%beg/bc_z%end
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%beg, 1, &
                                      q_cons_buffer_in(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                    ! PBC only at end of the sub-domain
                else

                    ! Packing the data to be sent to bc_z%end
                    do l = p - buff_size + 1, p
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                do i = 1, sys_size
                                    r = sys_size*(j + buff_size) &
                                        + sys_size*(m + 2*buff_size + 1)* &
                                        (k + buff_size) + sys_size* &
                                        (m + 2*buff_size + 1)* &
                                        (n + 2*buff_size + 1)* &
                                        (l - p + buff_size - 1) + (i - 1)
                                    if (present(q_particle)) then
                                        q_cons_buffer_out(r) = &
                                            q_particle%sf(j, k, l)
                                    else
                                        q_cons_buffer_out(r) = &
                                            q_cons_vf(i)%sf(j, k, l)
                                    end if
                                end do
                            end do
                        end do
                    end do

                    ! Sending/receiving the data to/from bc_z%end/bc_z%end
                    call MPI_SENDRECV(q_cons_buffer_out(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%end, 0, &
                                      q_cons_buffer_in(0), buff_size* &
                                      sys_size*(m + 2*buff_size + 1)* &
                                      (n + 2*buff_size + 1), &
                                      mpi_p, bc_z%end, 1, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                      ierr)

                end if

                ! Unpacking the data received from bc_z%end
                do l = p + 1, p + buff_size
                    do k = -buff_size, n + buff_size
                        do j = -buff_size, m + buff_size
                            do i = 1, sys_size
                                r = sys_size*(j + buff_size) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (k + buff_size) + (i - 1) &
                                    + sys_size*(m + 2*buff_size + 1)* &
                                    (n + 2*buff_size + 1)*(l - p - 1)
                                if (present(q_particle)) then
                                    q_particle%sf(j, k, l) = q_cons_buffer_in(r)
                                else
                                    q_cons_vf(i)%sf(j, k, l) = q_cons_buffer_in(r)
                                end if
#if defined(__INTEL_COMPILER)
                                if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                                    print *, "Error", j, k, l, i
                                    error stop "NaN(s) in recv"
                                end if
#endif
                            end do
                        end do
                    end do
                end do

            end if

        end if

        ! END: Communications in the z-direction

#endif

    end subroutine s_mpi_sendrecv_cons_vars_buffer_regions

    !>  This subroutine gathers the Silo database metadata for
        !!      the spatial extents in order to boost the performance of
        !!      the multidimensional visualization.
        !!  @param spatial_extents Spatial extents for each processor's sub-domain. First dimension
        !!  corresponds to the minimum and maximum values, respectively, while
        !!  the second dimension corresponds to the processor rank.
    subroutine s_mpi_gather_spatial_extents(spatial_extents)

        real(wp), dimension(1:, 0:), intent(INOUT) :: spatial_extents

#ifdef MFC_MPI

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
        else

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

        end if

#endif

    end subroutine s_mpi_gather_spatial_extents

    !>  This subroutine collects the sub-domain cell-boundary or
        !!      cell-center locations data from all of the processors and
        !!      puts back together the grid of the entire computational
        !!      domain on the rank 0 processor. This is only done for 1D
        !!      simulations.
    subroutine s_mpi_defragment_1d_grid_variable

#ifdef MFC_MPI

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
    subroutine s_mpi_gather_data_extents(q_sf, data_extents)

        real(wp), dimension(:, :, :), intent(in) :: q_sf

        real(wp), &
            dimension(1:2, 0:num_procs - 1), &
            intent(inout) :: data_extents

#ifdef MFC_MPI

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
    subroutine s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)

        real(wp), &
            dimension(0:m), &
            intent(in) :: q_sf

        real(wp), &
            dimension(0:m), &
            intent(inout) :: q_root_sf

#ifdef MFC_MPI

        ! Gathering the sub-domain flow variable data from all the processes
        ! and putting it back together for the entire computational domain
        ! on the process with rank 0
        call MPI_GATHERV(q_sf(0), m + 1, mpi_p, &
                         q_root_sf(0), recvcounts, displs, &
                         mpi_p, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_defragment_1d_flow_variable

    !> Deallocation procedures for the module
    subroutine s_finalize_mpi_proxy_module

#ifdef MFC_MPI

        ! Deallocating the conservative variables buffer vectors
        if (buff_size > 0) then
            deallocate (q_cons_buffer_in)
            deallocate (q_cons_buffer_out)
        end if

        ! Deallocating the receive counts and the displacement vector
        ! variables used in variable-gather communication procedures
        if ((format == 1 .and. n > 0) .or. n == 0) then
            deallocate (recvcounts)
            deallocate (displs)
        end if

#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
