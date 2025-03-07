!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

!> @brief This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required pre-processing
!!              communication goals.
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_common

    use m_helper

    implicit none

    integer, private :: err_code, ierr, v_size !<
    !! Generic flags used to identify and report MPI errors

    real(wp), private, allocatable, dimension(:), target :: q_prims_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! primitive variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(wp), private, allocatable, dimension(:), target :: q_prims_buff_recv !<
    !! q_prims_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average primitive variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    integer :: halo_size

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_mpi_proxy_module

#ifdef MFC_MPI

        ! Allocating q_prims_buff_send/recv. Please note that
        ! for the sake of simplicity, both variables are provided sufficient
        ! storage to hold the largest buffer in the computational domain.


        if (n > 0) then
            if (p > 0) then

                halo_size = -1 + buff_size*sys_size* &
                                         & (m + 2*buff_size + 1)* &
                                         & (n + 2*buff_size + 1)* &
                                         & (p + 2*buff_size + 1)/ &
                                         & (min(m, n, p) + 2*buff_size + 1)
            else
                halo_size = -1 + buff_size*sys_size* &
                                         & (max(m, n) + 2*buff_size + 1)
            end if
        else
            halo_size = -1 + buff_size*sys_size
        end if

        v_size = sys_size

        allocate(q_prims_buff_send(0:halo_size))

        allocate(q_prims_buff_recv(0:ubound(q_prims_buff_send, 1)))

#endif

    end subroutine s_initialize_mpi_proxy_module


    !> Since only processor with rank 0 is in charge of reading
        !!       and checking the consistency of the user provided inputs,
        !!       these are not available to the remaining processors. This
        !!       subroutine is then in charge of broadcasting the required
        !!       information.
    subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI

        ! Generic loop iterator
        integer :: i

        ! Logistics
        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 't_step_start', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 'loops_x', 'loops_y', 'loops_z', 'model_eqns', 'num_fluids',     &
            & 'weno_order', 'precision', 'perturb_flow_fluid', &
            & 'perturb_sph_fluid', 'num_patches', 'thermal', 'nb', 'dist_type',&
            & 'R0_type', 'relax_model', 'num_ibs', 'n_start', 'elliptic_smoothing_iters' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'old_grid','old_ic','stretch_x','stretch_y','stretch_z',&
            & 'cyl_coord','mpp_lim','hypoelasticity', 'relax', 'parallel_io',  &
            & 'perturb_flow', 'perturb_sph', 'mixlayer_vel_profile',           &
            & 'mixlayer_perturb', 'bubbles_euler', 'polytropic', 'polydisperse',&
            & 'qbmm', 'file_per_process', 'adv_n', 'ib' , 'cfl_adap_dt',       &
            & 'cfl_const_dt', 'cfl_dt', 'surface_tension',                     &
            & 'hyperelasticity', 'pre_stress', 'viscous', 'bubbles_lagrange',  &
            & 'elliptic_smoothing' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(fluid_rho(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in [ 'x_domain%beg', 'x_domain%end', 'y_domain%beg',         &
            & 'y_domain%end', 'z_domain%beg', 'z_domain%end', 'a_x', 'a_y',    &
            & 'a_z', 'x_a', 'x_b', 'y_a', 'y_b', 'z_a', 'z_b', 'bc_x%beg',     &
            & 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',      &
            & 'perturb_flow_mag', 'pref', 'rhoref', 'poly_sigma', 'R0ref',     &
            & 'Web', 'Ca', 'Re_inv', 'sigR', 'sigV', 'rhoRV', 'palpha_eps',    &
            & 'ptgalpha_eps', 'sigma', 'pi_fac', 'mixlayer_vel_coef',          &
            & 'mixlayer_domain' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, num_patches_max
            #:for VAR in [ 'geometry', 'smooth_patch_id']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%non_axis_sym, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'epsilon',     &
                & 'beta', 'smooth_coeff', 'rho', 'p0', 'm0', 'r0', 'v0',       &
                & 'pres', 'gamma', 'pi_inf', 'hcid', 'cv', 'qv', 'qvp',        &
                & 'model_threshold', 'cf_val']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ '2', '3', '4', '5', '6', '7', '8', '9']
                call MPI_BCAST(patch_icpp(i)%a(${VAR}$), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%model_filepath, len(patch_icpp(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'model_translate', 'model_scale', 'model_rotate', &
                'normal', 'radii', 'vel', 'tau_e', 'alpha_rho', 'alpha' ]
                call MPI_BCAST(patch_icpp(i)%${VAR}$, size(patch_icpp(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%model_spc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            ! Broadcast IB variables
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_filepath, len(patch_ib(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_threshold, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_spc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'c', 'p', 't', 'm', 'theta']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%slip, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'model_translate', 'model_scale', 'model_rotate']
                call MPI_BCAST(patch_ib(i)%${VAR}$, size(patch_ib(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v', 'G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do
#endif

    end subroutine s_mpi_bcast_user_inputs

    !> Description: This subroutine takes care of efficiently distributing
        !!              the computational domain among the available processors
        !!             as well as recomputing some of the global parameters so
        !!              that they reflect the configuration of sub-domain that is
        !!              overseen by the local processor.
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

                ! Preliminary uniform cell-width spacing
                if (old_grid .neqv. .true.) then
                    dz = (z_domain%end - z_domain%beg)/real(p + 1, wp)
                end if

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
                if (proc_coords(3) > 0 .or. (bc_z%beg == -1 .and. num_procs_z > 1)) then
                    proc_coords(3) = proc_coords(3) - 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%beg, ierr)
                    proc_coords(3) = proc_coords(3) + 1
                end if

                ! Boundary condition at the end
                if (proc_coords(3) < num_procs_z - 1 .or. (bc_z%end == -1 .and. num_procs_z > 1)) then
                    proc_coords(3) = proc_coords(3) + 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%end, ierr)
                    proc_coords(3) = proc_coords(3) - 1
                end if

                ! Beginning and end sub-domain boundary locations
                if (parallel_io .neqv. .true.) then
                    if (old_grid .neqv. .true.) then
                        if (proc_coords(3) < rem_cells) then
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3))
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1) &
                                                                  - (num_procs_z - rem_cells))
                        else
                            z_domain%beg = z_domain%beg + dz*real((p + 1)* &
                                                                  proc_coords(3) + rem_cells)
                            z_domain%end = z_domain%end - dz*real((p + 1)* &
                                                                  (num_procs_z - proc_coords(3) - 1))
                        end if
                    end if
                else
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

            ! Preliminary uniform cell-width spacing
            if (old_grid .neqv. .true.) then
                dy = (y_domain%end - y_domain%beg)/real(n + 1, wp)
            end if

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
            if (proc_coords(2) > 0 .or. (bc_y%beg == -1 .and. num_procs_y > 1)) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%beg, ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. (bc_y%end == -1 .and. num_procs_y > 1)) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%end, ierr)
                proc_coords(2) = proc_coords(2) - 1
            end if

            ! Beginning and end sub-domain boundary locations
            if (parallel_io .neqv. .true.) then
                if (old_grid .neqv. .true.) then
                    if (proc_coords(2) < rem_cells) then
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2))
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1) &
                                                              - (num_procs_y - rem_cells))
                    else
                        y_domain%beg = y_domain%beg + dy*real((n + 1)* &
                                                              proc_coords(2) + rem_cells)
                        y_domain%end = y_domain%end - dy*real((n + 1)* &
                                                              (num_procs_y - proc_coords(2) - 1))
                    end if
                end if
            else
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

        ! Preliminary uniform cell-width spacing
        if (old_grid .neqv. .true.) then
            dx = (x_domain%end - x_domain%beg)/real(m + 1, wp)
        end if

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing any remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1
                exit
            end if
        end do

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. (bc_x%beg == -1 .and. num_procs_x > 1)) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. (bc_x%end == -1 .and. num_procs_x > 1)) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

        ! Beginning and end sub-domain boundary locations
        if (parallel_io .neqv. .true.) then
            if (old_grid .neqv. .true.) then
                if (proc_coords(1) < rem_cells) then
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1))
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1) &
                                                          - (num_procs_x - rem_cells))
                else
                    x_domain%beg = x_domain%beg + dx*real((m + 1)* &
                                                          proc_coords(1) + rem_cells)
                    x_domain%end = x_domain%end - dx*real((m + 1)* &
                                                          (num_procs_x - proc_coords(1) - 1))
                end if
            end if
        else
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if

#endif

    end subroutine s_mpi_decompose_computational_domain

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average primitive variables by communicating
        !!      with the neighboring processors.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_variables_buffers(q_prim_vf, &
                                                mpi_dir, &
                                                pbc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0

        integer :: pack_offset, unpack_offset

        real(wp), pointer :: p_send, p_recv

#ifdef MFC_MPI

        buffer_counts = (/ &
                buff_size*sys_size*(n + 1)*(p + 1), &
                buff_size*sys_size*(m + 2*buff_size + 1)*(p + 1), &
                buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                /)

        buffer_count = buffer_counts(mpi_dir)
        boundary_conditions = (/bc_x, bc_y, bc_z/)
        beg_end = (/boundary_conditions(mpi_dir)%beg, boundary_conditions(mpi_dir)%end/)
        beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0

        ! Implements:
        ! pbc_loc  bc_x >= 0 -> [send/recv]_tag  [dst/src]_proc
        ! -1 (=0)      0            ->     [1,0]       [0,0]      | 0 0 [1,0] [beg,beg]
        ! -1 (=0)      1            ->     [0,0]       [1,0]      | 0 1 [0,0] [end,beg]
        ! +1 (=1)      0            ->     [0,1]       [1,1]      | 1 0 [0,1] [end,end]
        ! +1 (=1)      1            ->     [1,1]       [0,1]      | 1 1 [1,1] [beg,end]

        send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
        recv_tag = f_logical_to_int(pbc_loc == 1)

        dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
        src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))

        grid_dims = (/m, n, p/)

        pack_offset = 0
        if (f_xor(pbc_loc == 1, beg_end_geq_0)) then
            pack_offset = grid_dims(mpi_dir) - buff_size + 1
        end if

        unpack_offset = 0
        if (pbc_loc == 1) then
            unpack_offset = grid_dims(mpi_dir) + buff_size + 1
        end if

        ! Pack Buffer to Send
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, buff_size - 1
                                do i = 1, sys_size
                                    r = (i - 1) + v_size*(j + buff_size*(k + (n + 1)*l))
                                    q_prims_buff_send(r) = q_prim_vf(i)%sf(j + pack_offset, k, l)
                                end do
                            end do
                        end do
                    end do
                #:elif mpi_dir == 2
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, p
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k + buff_size*l))
                                    q_prims_buff_send(r) = q_prim_vf(i)%sf(j, k + pack_offset, l)
                                end do
                            end do
                        end do
                    end do
                #:else
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, buff_size - 1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                    q_prims_buff_send(r) = q_prim_vf(i)%sf(j, k, l + pack_offset)
                                end do
                            end do
                        end do
                    end do
                #:endif
            end if
        #:endfor

        p_send => q_prims_buff_send(0)
        p_recv => q_prims_buff_recv(0)

        call MPI_SENDRECV( &
        p_send, buffer_count, mpi_p, dst_proc, send_tag, &
        p_recv, buffer_count, mpi_p, src_proc, recv_tag, &
        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        
        ! Unpack Received Buffer
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do l = 0, p
                        do k = 0, n
                            do j = -buff_size, -1
                                do i = 1, sys_size
                                    r = (i - 1) + v_size* &
                                        (j + buff_size*((k + 1) + (n + 1)*l))
                                        q_prim_vf(i)%sf(j + unpack_offset, k, l) = q_prims_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_prim_vf(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do

               #:elif mpi_dir == 2
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, p
                            do k = -buff_size, -1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + buff_size*l))
                                         q_prim_vf(i)%sf(j, k + unpack_offset, l) = q_prims_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_prim_vf(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do
                #:else
                    ! Unpacking buffer from bc_z%beg
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = -buff_size, -1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)* &
                                          (l + buff_size)))
                                          q_prim_vf(i)%sf(j, k, l + unpack_offset) = q_prims_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(q_prim_vf(i)%sf(j, k, l))) then
                                        print *, "Error", j, k, l, i
                                        error stop "NaN(s) in recv"
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do
                #:endif
            end if
        #:endfor

#endif

    end subroutine s_mpi_sendrecv_variables_buffers

end module m_mpi_proxy
