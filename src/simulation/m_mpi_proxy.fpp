!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
module m_mpi_proxy

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_common

    use ieee_arithmetic
    ! ==========================================================================

    implicit none

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), q_cons_buff_send)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), q_cons_buff_recv)
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), ib_buff_send)
    @:CRAY_DECLARE_GLOBAL(integer, dimension(:), ib_buff_recv)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), c_divs_buff_send)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), c_divs_buff_recv)
    !$acc declare link(q_cons_buff_recv, q_cons_buff_send)
    !$acc declare link(ib_buff_send, ib_buff_recv)
    !$acc declare link(c_divs_buff_send, c_divs_buff_recv)
#else
    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! conservative variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average conservative variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    real(kind(0d0)), private, allocatable, dimension(:) :: c_divs_buff_send !<
    !! c_divs_buff_send is utilized to send and unpack the buffer of the cell-
    !! centered color function derivatives, for a single computational domain
    !! boundary at the time, to the the relevant neighboring processor

    real(kind(0d0)), private, allocatable, dimension(:) :: c_divs_buff_recv
    !! c_divs_buff_recv is utilized to receiver and unpack the buffer of the cell-
    !! centered color function derivatives, for a single computational domain
    !! boundary at the time, from the relevant neighboring processor

    integer, private, allocatable, dimension(:) :: ib_buff_send !<
    !! This variable is utilized to pack and send the buffer of the immersed
    !! boundary markers, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    integer, private, allocatable, dimension(:) :: ib_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the
    !! immersed boundary markers, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !$acc declare create(q_cons_buff_send, q_cons_buff_recv)
    !$acc declare create( ib_buff_send, ib_buff_recv)
    !$acc declare create(c_divs_buff_send, c_divs_buff_recv)
#endif
    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: err_code, ierr, v_size
    !> @}
    !$acc declare create(v_size)

    !real :: s_time, e_time
    !real :: compress_time, mpi_time, decompress_time
    !integer :: nCalls_time = 0

    integer :: nVars !< nVars for surface tension communication
    !$acc declare create(nVars)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_mpi_proxy_module

#ifdef MFC_MPI

        ! Allocating q_cons_buff_send/recv and ib_buff_send/recv. Please note that
        ! for the sake of simplicity, both variables are provided sufficient
        ! storage to hold the largest buffer in the computational domain.

        if (qbmm .and. .not. polytropic) then
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*(sys_size + 2*nb*4)))
            end if

            @:ALLOCATE_GLOBAL(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size + 2*nb*4
        else
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE_GLOBAL(q_cons_buff_send(0:-1 + buff_size*sys_size))
            end if

            @:ALLOCATE_GLOBAL(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

            v_size = sys_size
        end if

        if (.not. f_is_default(sigma)) then
            nVars = num_dims + 1
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE_GLOBAL(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE_GLOBAL(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE_GLOBAL(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)))
            end if

            @:ALLOCATE_GLOBAL(c_divs_buff_recv(0:ubound(c_divs_buff_send, 1)))
        end if
        !$acc update device(v_size, nVars)

#endif

    end subroutine s_initialize_mpi_proxy_module

    !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.
    subroutine s_mpi_bcast_user_inputs

#ifdef MFC_MPI

        integer :: i, j !< Generic loop iterator

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['k_x', 'k_y', 'k_z', 'w_x', 'w_y', 'w_z', 'p_x', 'p_y', &
            & 'p_z', 'g_x', 'g_y', 'g_z']
            call MPI_BCAST(${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in ['t_step_old', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 't_step_start','t_step_stop','t_step_save','t_step_print',       &
            & 'model_eqns','time_stepper', 'riemann_solver', 'low_Mach',       &
            & 'wave_speeds', 'avg_state', 'precision', 'bc_x%beg', 'bc_x%end', &
            & 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',  'fd_order',     &
            & 'num_probes', 'num_integrals', 'bubble_model', 'thermal',        &
            & 'R0_type', 'num_source', 'relax_model', 'num_ibs']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'adv_alphan', 'mpp_lim',    &
            & 'mp_weno', 'rdma_mpi', 'weno_flat', 'riemann_flat',                &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles', 'polytropic',          &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref','R0ref','Web','Ca', 'sigma', &
            & 'Re_inv', 'poly_sigma', 'palpha_eps', 'ptgalpha_eps', 'pi_fac',    &
            & 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve2', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3', &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3', &
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',    &
            & 'z_domain%beg', 'z_domain%end']
            call MPI_BCAST(${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(mapped_weno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(teno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v','G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', &
                & 'x_centroid', 'y_centroid', 'c', 'm', 'p', 't', 'theta', 'slip' ]
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        do j = 1, num_probes_max
            do i = 1, 3
                call MPI_BCAST(acoustic(j)%loc(i), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            end do

            call MPI_BCAST(acoustic(j)%dipole, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'pulse', 'support', 'num_elements', 'element_on' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'mag', 'length', 'height', &
                'wavelength', 'frequency', 'gauss_sigma_dist', 'gauss_sigma_time', &
                'npulse', 'dir', 'delay', 'foc_length', 'aperture', &
                'element_spacing_angle', 'element_polygon_ratio', 'rotate_angle' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'x','y','z' ]
                call MPI_BCAST(probe(j)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax' ]
                call MPI_BCAST(integral(j)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

#endif

    end subroutine s_mpi_bcast_user_inputs

    !>  The purpose of this procedure is to optimally decompose
        !!      the computational domain among the available processors.
        !!      This is performed by attempting to award each processor,
        !!      in each of the coordinate directions, approximately the
        !!      same number of cells, and then recomputing the affected
        !!      global parameters.
    subroutine s_mpi_decompose_computational_domain

#ifdef MFC_MPI

        integer :: num_procs_x, num_procs_y, num_procs_z !<
            !! Optimal number of processors in the x-, y- and z-directions

        real(kind(0d0)) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z !<
            !! Non-optimal number of processors in the x-, y- and z-directions

        real(kind(0d0)) :: fct_min !<
            !! Processor factorization (fct) minimization parameter

        integer :: MPI_COMM_CART !<
            !! Cartesian processor topology communicator

        integer :: rem_cells !<
            !! Remaining number of cells, in a particular coordinate direction,
            !! after the majority is divided up among the available processors

        integer :: i, j !< Generic loop iterators

        if (num_procs == 1 .and. parallel_io) then
            do i = 1, num_dims
                start_idx(i) = 0
            end do
            return
        end if

        ! 3D Cartesian Processor Topology ==================================
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
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
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

                    ! Initial estimate of optimal processor topology
                    num_procs_x = 1
                    num_procs_y = 1
                    num_procs_z = num_procs
                    ierr = -1

                    ! Benchmarking the quality of this initial guess
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                       - (n + 1)/tmp_num_procs_y) &
                              + 10d0*abs((n + 1)/tmp_num_procs_y &
                                         - (p + 1)/tmp_num_procs_z)

                    ! Optimization of the initial processor topology
                    do i = 1, num_procs

                        if (mod(num_procs, i) == 0 &
                            .and. &
                            (m + 1)/i >= num_stcls_min*weno_order) then

                            do j = 1, num_procs/i

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

                ! Verifying that a valid decomposition of the computational
                ! domain has been established. If not, the simulation exits.
                if (proc_rank == 0 .and. ierr == -1) then
                    call s_mpi_abort('Unsupported combination of values '// &
                                     'of num_procs, m, n, p and '// &
                                     'weno_order. Exiting ...')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/num_procs_x, &
                                                          num_procs_y, num_procs_z/), &
                                     (/.true., .true., .true./), &
                                     .false., MPI_COMM_CART, ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, &
                                     proc_coords, ierr)
                ! END: 3D Cartesian Processor Topology =============================

                ! Global Parameters for z-direction ================================

                ! Number of remaining cells
                rem_cells = mod(p + 1, num_procs_z)

                ! Optimal number of cells per processor
                p = (p + 1)/num_procs_z - 1

                ! Distributing the remaining cells
                do i = 1, rem_cells
                    if (proc_coords(3) == i - 1) then
                        p = p + 1; exit
                    end if
                end do

                ! Boundary condition at the beginning
                if (proc_coords(3) > 0 .or. bc_z%beg == -1) then
                    proc_coords(3) = proc_coords(3) - 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%beg, ierr)
                    proc_coords(3) = proc_coords(3) + 1
                end if

                ! Boundary condition at the end
                if (proc_coords(3) < num_procs_z - 1 .or. bc_z%end == -1) then
                    proc_coords(3) = proc_coords(3) + 1
                    call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                       bc_z%end, ierr)
                    proc_coords(3) = proc_coords(3) - 1
                end if

                if (parallel_io) then
                    if (proc_coords(3) < rem_cells) then
                        start_idx(3) = (p + 1)*proc_coords(3)
                    else
                        start_idx(3) = (p + 1)*proc_coords(3) + rem_cells
                    end if
                end if
                ! ==================================================================

                ! 2D Cartesian Processor Topology ==================================
            else

                ! Initial estimate of optimal processor topology
                num_procs_x = 1
                num_procs_y = num_procs
                ierr = -1

                ! Benchmarking the quality of this initial guess
                tmp_num_procs_x = num_procs_x
                tmp_num_procs_y = num_procs_y
                fct_min = 10d0*abs((m + 1)/tmp_num_procs_x &
                                   - (n + 1)/tmp_num_procs_y)

                ! Optimization of the initial processor topology
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

                ! Verifying that a valid decomposition of the computational
                ! domain has been established. If not, the simulation exits.
                if (proc_rank == 0 .and. ierr == -1) then
                    call s_mpi_abort('Unsupported combination of values '// &
                                     'of num_procs, m, n and '// &
                                     'weno_order. Exiting ...')
                end if

                ! Creating new communicator using the Cartesian topology
                call MPI_CART_CREATE(MPI_COMM_WORLD, 2, (/num_procs_x, &
                                                          num_procs_y/), (/.true., &
                                                                           .true./), .false., MPI_COMM_CART, &
                                     ierr)

                ! Finding the Cartesian coordinates of the local process
                call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, &
                                     proc_coords, ierr)

            end if
            ! END: 2D Cartesian Processor Topology =============================

            ! Global Parameters for y-direction ================================

            ! Number of remaining cells
            rem_cells = mod(n + 1, num_procs_y)

            ! Optimal number of cells per processor
            n = (n + 1)/num_procs_y - 1

            ! Distributing the remaining cells
            do i = 1, rem_cells
                if (proc_coords(2) == i - 1) then
                    n = n + 1; exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(2) > 0 .or. bc_y%beg == -1) then
                proc_coords(2) = proc_coords(2) - 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%beg, ierr)
                proc_coords(2) = proc_coords(2) + 1
            end if

            ! Boundary condition at the end
            if (proc_coords(2) < num_procs_y - 1 .or. bc_y%end == -1) then
                proc_coords(2) = proc_coords(2) + 1
                call MPI_CART_RANK(MPI_COMM_CART, proc_coords, &
                                   bc_y%end, ierr)
                proc_coords(2) = proc_coords(2) - 1
            end if

            if (parallel_io) then
                if (proc_coords(2) < rem_cells) then
                    start_idx(2) = (n + 1)*proc_coords(2)
                else
                    start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
                end if
            end if

            ! ==================================================================

            ! 1D Cartesian Processor Topology ==================================
        else

            ! Optimal processor topology
            num_procs_x = num_procs

            ! Creating new communicator using the Cartesian topology
            call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), &
                                 (/.true./), .false., MPI_COMM_CART, &
                                 ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, &
                                 proc_coords, ierr)

        end if
        ! ==================================================================

        ! Global Parameters for x-direction ================================

        ! Number of remaining cells
        rem_cells = mod(m + 1, num_procs_x)

        ! Optimal number of cells per processor
        m = (m + 1)/num_procs_x - 1

        ! Distributing the remaining cells
        do i = 1, rem_cells
            if (proc_coords(1) == i - 1) then
                m = m + 1; exit
            end if
        end do

        ! Boundary condition at the beginning
        if (proc_coords(1) > 0 .or. bc_x%beg == -1) then
            proc_coords(1) = proc_coords(1) - 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
            proc_coords(1) = proc_coords(1) + 1
        end if

        ! Boundary condition at the end
        if (proc_coords(1) < num_procs_x - 1 .or. bc_x%end == -1) then
            proc_coords(1) = proc_coords(1) + 1
            call MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
            proc_coords(1) = proc_coords(1) - 1
        end if

        if (parallel_io) then
            if (proc_coords(1) < rem_cells) then
                start_idx(1) = (m + 1)*proc_coords(1)
            else
                start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
            end if
        end if
        ! ==================================================================

#endif

    end subroutine s_mpi_decompose_computational_domain

    !>  The goal of this procedure is to populate the buffers of
        !!      the grid variables by communicating with the neighboring
        !!      processors. Note that only the buffers of the cell-width
        !!      distributions are handled in such a way. This is because
        !!      the buffers of cell-boundary locations may be calculated
        !!      directly from those of the cell-width distributions.
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_grid_variables_buffers(mpi_dir, pbc_loc)

        integer, intent(in) :: mpi_dir
        integer, intent(in) :: pbc_loc

        integer :: dst_proc(1:3)

#ifdef MFC_MPI

        ! MPI Communication in x-direction =================================
        if (mpi_dir == 1) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_x%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                        dx(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                        dx(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_x%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_x%beg/bc_x%end
                    call MPI_SENDRECV( &
                        dx(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                        dx(m + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_x%end/bc_x%end
                    call MPI_SENDRECV( &
                        dx(m - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                        dx(m + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in x-direction ============================

            ! MPI Communication in y-direction =================================
        elseif (mpi_dir == 2) then

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_y%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_y%end/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                        dy(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_y%beg/bc_y%beg
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                        dy(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_y%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_y%beg/bc_y%end
                    call MPI_SENDRECV( &
                        dy(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                        dy(n + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_y%end/bc_y%end
                    call MPI_SENDRECV( &
                        dy(n - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                        dy(n + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if
            ! END: MPI Communication in y-direction ============================

            ! MPI Communication in z-direction =================================
        else

            if (pbc_loc == -1) then      ! PBC at the beginning

                if (bc_z%end >= 0) then      ! PBC at the beginning and end

                    ! Send/receive buffer to/from bc_z%end/bc_z%beg
                    call MPI_SENDRECV( &
                        dz(p - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                        dz(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the beginning only

                    ! Send/receive buffer to/from bc_z%beg/bc_z%beg
                    call MPI_SENDRECV( &
                        dz(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                        dz(-buff_size), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            else                        ! PBC at the end

                if (bc_z%beg >= 0) then      ! PBC at the end and beginning

                    ! Send/receive buffer to/from bc_z%beg/bc_z%end
                    call MPI_SENDRECV( &
                        dz(0), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                        dz(p + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                else                        ! PBC at the end only

                    ! Send/receive buffer to/from bc_z%end/bc_z%end
                    call MPI_SENDRECV( &
                        dz(p - buff_size + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                        dz(p + 1), buff_size, &
                        MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                end if

            end if

        end if
        ! END: MPI Communication in z-direction ============================

#endif

    end subroutine s_mpi_sendrecv_grid_variables_buffers

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
    subroutine s_mpi_sendrecv_variables_buffers(q_cons_vf, &
                                                pb, mv, &
                                                mpi_dir, &
                                                pbc_loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0

        integer :: pack_offsets(1:3), unpack_offsets(1:3)
        integer :: pack_offset, unpack_offset

#ifdef MFC_MPI

        !$acc update device(v_size)

        if (qbmm .and. .not. polytropic) then
            buffer_counts = (/ &
                            buff_size*(sys_size + 2*nb*4)*(n + 1)*(p + 1), &
                            buff_size*(sys_size + 2*nb*4)*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        else
            buffer_counts = (/ &
                            buff_size*sys_size*(n + 1)*(p + 1), &
                            buff_size*sys_size*(m + 2*buff_size + 1)*(p + 1), &
                            buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                            /)
        end if

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
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j + pack_offset, k, l)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(4) gang vector default(present) private(r)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, buff_size - 1
                                    do i = sys_size + 1, sys_size + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                (j + buff_size*(k + (n + 1)*l))
                                            q_cons_buff_send(r) = pb(j + pack_offset, k, l, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, buff_size - 1
                                    do i = sys_size + 1, sys_size + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                (j + buff_size*(k + (n + 1)*l))
                                            q_cons_buff_send(r) = mv(j + pack_offset, k, l, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:elif mpi_dir == 2
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, p
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k + buff_size*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k + pack_offset, l)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, p
                                do k = 0, buff_size - 1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 (k + buff_size*l))
                                            q_cons_buff_send(r) = pb(j, k + pack_offset, l, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, p
                                do k = 0, buff_size - 1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 (k + buff_size*l))
                                            q_cons_buff_send(r) = mv(j, k + pack_offset, l, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:else
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, buff_size - 1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j, k, l + pack_offset)
                                end do
                            end do
                        end do
                    end do

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, buff_size - 1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                            q_cons_buff_send(r) = pb(j, k, l + pack_offset, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, buff_size - 1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                            q_cons_buff_send(r) = mv(j, k, l + pack_offset, i - sys_size, q)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:endif
            end if
        #:endfor

        ! Send/Recv
        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    !$acc host_data use_device(q_cons_buff_recv, q_cons_buff_send, ib_buff_recv, ib_buff_send)
                #:else
                    !$acc update host(q_cons_buff_send, ib_buff_send)
                #:endif

                call MPI_SENDRECV( &
                    q_cons_buff_send(0), buffer_count, MPI_DOUBLE_PRECISION, dst_proc, send_tag, &
                    q_cons_buff_recv(0), buffer_count, MPI_DOUBLE_PRECISION, src_proc, recv_tag, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                #:if rdma_mpi
                    !$acc end host_data
                    !$acc wait
                #:else
                    !$acc update device(q_cons_buff_recv)
                #:endif
            end if
        #:endfor

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
                                    q_cons_vf(i)%sf(j + unpack_offset, k, l) = q_cons_buff_recv(r)
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

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do l = 0, p
                            do k = 0, n
                                do j = -buff_size, -1
                                    do i = sys_size + 1, sys_size + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                (j + buff_size*((k + 1) + (n + 1)*l))
                                            pb(j + unpack_offset, k, l, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do l = 0, p
                            do k = 0, n
                                do j = -buff_size, -1
                                    do i = sys_size + 1, sys_size + 4
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                (j + buff_size*((k + 1) + (n + 1)*l))
                                            mv(j + unpack_offset, k, l, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:elif mpi_dir == 2
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, sys_size
                        do l = 0, p
                            do k = -buff_size, -1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + v_size* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + buff_size*l))
                                    q_cons_vf(i)%sf(j, k + unpack_offset, l) = q_cons_buff_recv(r)
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

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, p
                                do k = -buff_size, -1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + buff_size*l))
                                            pb(j, k + unpack_offset, l, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = 0, p
                                do k = -buff_size, -1
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + buff_size*l))
                                            mv(j, k + unpack_offset, l, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
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
                                    q_cons_vf(i)%sf(j, k, l + unpack_offset) = q_cons_buff_recv(r)
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

                    if (qbmm .and. .not. polytropic) then
                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = -buff_size, -1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)* &
                                                  (l + buff_size)))
                                            pb(j, k, l + unpack_offset, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(5) gang vector default(present) private(r)
                        do i = sys_size + 1, sys_size + 4
                            do l = -buff_size, -1
                                do k = -buff_size, n + buff_size
                                    do j = -buff_size, m + buff_size
                                        do q = 1, nb
                                            r = (i - 1) + (q - 1)*4 + nb*4 + v_size* &
                                                ((j + buff_size) + (m + 2*buff_size + 1)* &
                                                 ((k + buff_size) + (n + 2*buff_size + 1)* &
                                                  (l + buff_size)))
                                            mv(j, k, l + unpack_offset, i - sys_size, q) = q_cons_buff_recv(r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end if
                #:endif
            end if
        #:endfor

#endif

    end subroutine s_mpi_sendrecv_variables_buffers

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
    subroutine s_mpi_sendrecv_ib_buffers(ib_markers, gp_layers)

        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in) :: gp_layers

        integer :: i, j, k, l, r !< Generic loop iterators

#ifdef MFC_MPI

        if (n > 0) then
            if (p > 0) then
                @:ALLOCATE_GLOBAL(ib_buff_send(0:-1 + gp_layers * &
                                        & (m + 2*gp_layers + 1)* &
                                        & (n + 2*gp_layers + 1)* &
                                        & (p + 2*gp_layers + 1)/ &
                                        & (min(m, n, p) + 2*gp_layers + 1)))
            else
                @:ALLOCATE_GLOBAL(ib_buff_send(0:-1 + gp_layers* &
                                        & (max(m, n) + 2*gp_layers + 1)))
            end if
        else
            @:ALLOCATE_GLOBAL(ib_buff_send(0:-1 + gp_layers))
        end if
        @:ALLOCATE_GLOBAL(ib_buff_recv(0:ubound(ib_buff_send, 1)))

        !nCalls_time = nCalls_time + 1

        ! MPI Communication in x-direction =================================
        if (bc_x%beg >= 0) then      ! PBC at the beginning

            if (bc_x%end >= 0) then      ! PBC at the beginning and end

                ! Packing buffer to be sent to bc_x%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = 0, n
                        do j = m - gp_layers + 1, m
                            r = ((j - m - 1) + gp_layers*((k + 1) + (n + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send, ib_buff_recv, ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send, ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the beginning only

                ! Packing buffer to be sent to bc_x%beg
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = 0, n
                        do j = 0, gp_layers - 1
                            r = (j + gp_layers*(k + (n + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif
                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

#if defined(MFC_OpenACC)
            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if
#endif

            ! Unpacking buffer received from bc_x%beg
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = 0, p
                do k = 0, n
                    do j = -gp_layers, -1
                        r = (j + gp_layers*((k + 1) + (n + 1)*l))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if

        if (bc_x%end >= 0) then                ! PBC at the end

            if (bc_x%beg >= 0) then      ! PBC at the end and beginning

                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                ! Packing buffer to be sent to bc_x%beg
                do l = 0, p
                    do k = 0, n
                        do j = 0, gp_layers - 1
                            r = (j + gp_layers*(k + (n + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the end only

                ! Packing buffer to be sent to bc_x%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = 0, n
                        do j = m - gp_layers + 1, m
                            r = ((j - m - 1) + gp_layers*((k + 1) + (n + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if

            ! Unpacking buffer received from bc_x%end
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = 0, p
                do k = 0, n
                    do j = m + 1, m + gp_layers
                        r = ((j - m - 1) + gp_layers*(k + (n + 1)*l))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if
        ! END: MPI Communication in x-direction ============================

        ! MPI Communication in y-direction =================================

        if (bc_y%beg >= 0) then      ! PBC at the beginning

            if (bc_y%end >= 0) then      ! PBC at the beginning and end

                ! Packing buffer to be sent to bc_y%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = n - gp_layers + 1, n
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k - n + gp_layers - 1) + gp_layers*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the beginning only

                ! Packing buffer to be sent to bc_y%beg
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = 0, gp_layers - 1
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 (k + gp_layers*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

#if defined(MFC_OpenACC)
            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if
#endif

            ! Unpacking buffer received from bc_y%beg
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = 0, p
                do k = -gp_layers, -1
                    do j = -gp_layers, m + gp_layers
                        r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                             ((k + gp_layers) + gp_layers*l))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if

        if (bc_y%end >= 0) then                   ! PBC at the end

            if (bc_y%beg >= 0) then      ! PBC at the end and beginning

                ! Packing buffer to be sent to bc_y%beg
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = 0, gp_layers - 1
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 (k + gp_layers*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the end only

                ! Packing buffer to be sent to bc_y%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, p
                    do k = n - gp_layers + 1, n
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k - n + gp_layers - 1) + gp_layers*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

#if defined(MFC_OpenACC)
            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if
#endif

            ! Unpacking buffer received form bc_y%end
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = 0, p
                do k = n + 1, n + gp_layers
                    do j = -gp_layers, m + gp_layers
                        r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                             ((k - n - 1) + gp_layers*l))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if
        ! END: MPI Communication in y-direction ============================

        ! MPI Communication in z-direction =================================
        if (bc_z%beg >= 0) then      ! PBC at the beginning

            if (bc_z%end >= 0) then      ! PBC at the beginning and end

                ! Packing buffer to be sent to bc_z%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = p - gp_layers + 1, p
                    do k = -gp_layers, n + gp_layers
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k + gp_layers) + (n + 2*gp_layers + 1)* &
                                  (l - p + gp_layers - 1)))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the beginning only

                ! Packing buffer to be sent to bc_z%beg
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, gp_layers - 1
                    do k = -gp_layers, n + gp_layers
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k + gp_layers) + (n + 2*gp_layers + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif

                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

#if defined(MFC_OpenACC)
            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if
#endif

            ! Unpacking buffer from bc_z%beg
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = -gp_layers, -1
                do k = -gp_layers, n + gp_layers
                    do j = -gp_layers, m + gp_layers
                        r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                             ((k + gp_layers) + (n + 2*gp_layers + 1)* &
                              (l + gp_layers)))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if

        if (bc_z%end >= 0) then                       ! PBC at the end

            if (bc_z%beg >= 0) then      ! PBC at the end and beginning

                ! Packing buffer to be sent to bc_z%beg
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = 0, gp_layers - 1
                    do k = -gp_layers, n + gp_layers
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k + gp_layers) + (n + 2*gp_layers + 1)*l))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif
                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            else                        ! PBC at the end only

                ! Packing buffer to be sent to bc_z%end
                !$acc parallel loop collapse(3) gang vector default(present) private(r)
                do l = p - gp_layers + 1, p
                    do k = -gp_layers, n + gp_layers
                        do j = -gp_layers, m + gp_layers
                            r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                                 ((k + gp_layers) + (n + 2*gp_layers + 1)* &
                                  (l - p + gp_layers - 1)))
                            ib_buff_send(r) = ib_markers%sf(j, k, l)
                        end do
                    end do
                end do

                !call MPI_Barrier(MPI_COMM_WORLD, ierr)

#if defined(MFC_OpenACC)
                if (rdma_mpi) then
                    !$acc host_data use_device( ib_buff_recv, ib_buff_send )

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc wait
                else
#endif
                    !$acc update host(ib_buff_send)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        ib_buff_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        ib_buff_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

#if defined(MFC_OpenACC)
                end if
#endif

            end if

#if defined(MFC_OpenACC)
            if (rdma_mpi .eqv. .false.) then
                !$acc update device(ib_buff_recv)
            end if
#endif

            ! Unpacking buffer received from bc_z%end
            !$acc parallel loop collapse(3) gang vector default(present) private(r)
            do l = p + 1, p + gp_layers
                do k = -gp_layers, n + gp_layers
                    do j = -gp_layers, m + gp_layers
                        r = ((j + gp_layers) + (m + 2*gp_layers + 1)* &
                             ((k + gp_layers) + (n + 2*gp_layers + 1)* &
                              (l - p - 1)))
                        ib_markers%sf(j, k, l) = ib_buff_recv(r)
                    end do
                end do
            end do

        end if

        ! END: MPI Communication in z-direction ============================

#endif

    end subroutine s_mpi_sendrecv_ib_buffers

    subroutine s_mpi_sendrecv_capilary_variables_buffers(c_divs_vf, mpi_dir, pbc_loc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs_vf
        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0

        integer :: pack_offsets(1:3), unpack_offsets(1:3)
        integer :: pack_offset, unpack_offset

#ifdef MFC_MPI

        nVars = num_dims + 1
        !$acc update device(nVars)

        buffer_counts = (/ &
                        buff_size*nVars*(n + 1)*(p + 1), &
                        buff_size*nVars*(m + 2*buff_size + 1)*(p + 1), &
                        buff_size*nVars*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
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
                                do i = 1, nVars
                                    r = (i - 1) + nVars*(j + buff_size*(k + (n + 1)*l))
                                    c_divs_buff_send(r) = c_divs_vf(i)%sf(j + pack_offset, k, l)
                                end do
                            end do
                        end do
                    end do

                #:elif mpi_dir == 2
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, nVars
                        do l = 0, p
                            do k = 0, buff_size - 1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + nVars* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         (k + buff_size*l))
                                    c_divs_buff_send(r) = c_divs_vf(i)%sf(j, k + pack_offset, l)
                                end do
                            end do
                        end do
                    end do

                #:else
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do i = 1, nVars
                        do l = 0, buff_size - 1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + nVars* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                    c_divs_buff_send(r) = c_divs_vf(i)%sf(j, k, l + pack_offset)
                                end do
                            end do
                        end do
                    end do
                #:endif
            end if
        #:endfor

        ! Send/Recv
        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    !$acc host_data use_device(c_divs_buff_recv, c_divs_buff_send)
                #:else
                    !$acc update host(c_divs_buff_send)
                #:endif

                call MPI_SENDRECV( &
                    c_divs_buff_send(0), buffer_count, MPI_DOUBLE_PRECISION, dst_proc, send_tag, &
                    c_divs_buff_recv(0), buffer_count, MPI_DOUBLE_PRECISION, src_proc, recv_tag, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                #:if rdma_mpi
                    !$acc end host_data
                    !$acc wait
                #:else
                    !$acc update device(c_divs_buff_recv)
                #:endif
            end if
        #:endfor

        ! Unpack Received Buffer
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    !$acc parallel loop collapse(4) gang vector default(present) private(r)
                    do l = 0, p
                        do k = 0, n
                            do j = -buff_size, -1
                                do i = 1, nVars
                                    r = (i - 1) + nVars* &
                                        (j + buff_size*((k + 1) + (n + 1)*l))
                                    c_divs_vf(i)%sf(j + unpack_offset, k, l) = c_divs_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(c_divs_vf(i)%sf(j, k, l))) then
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
                    do i = 1, nVars
                        do l = 0, p
                            do k = -buff_size, -1
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + nVars* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + buff_size*l))
                                    c_divs_vf(i)%sf(j, k + unpack_offset, l) = c_divs_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(c_divs_vf(i)%sf(j, k, l))) then
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
                    do i = 1, nVars
                        do l = -buff_size, -1
                            do k = -buff_size, n + buff_size
                                do j = -buff_size, m + buff_size
                                    r = (i - 1) + nVars* &
                                        ((j + buff_size) + (m + 2*buff_size + 1)* &
                                         ((k + buff_size) + (n + 2*buff_size + 1)* &
                                          (l + buff_size)))
                                    c_divs_vf(i)%sf(j, k, l + unpack_offset) = c_divs_buff_recv(r)
#if defined(__INTEL_COMPILER)
                                    if (ieee_is_nan(c_divs_vf(i)%sf(j, k, l))) then
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

    end subroutine s_mpi_sendrecv_capilary_variables_buffers

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_mpi_proxy_module

#ifdef MFC_MPI

        ! Deallocating q_cons_buff_send and q_cons_buff_recv
        @:DEALLOCATE_GLOBAL(q_cons_buff_send, q_cons_buff_recv)
        if (ib) then
            @:DEALLOCATE_GLOBAL(ib_buff_send, ib_buff_recv)
        end if

        if (.not. f_is_default(sigma)) then
            @:DEALLOCATE_GLOBAL(c_divs_buff_send, c_divs_buff_recv)
        end if

#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
