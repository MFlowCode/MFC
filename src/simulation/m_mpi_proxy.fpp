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

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_common

    use m_nvtx

    use ieee_arithmetic

    implicit none

    integer, private, allocatable, dimension(:), target :: ib_buff_send !<
    !! This variable is utilized to pack and send the buffer of the immersed
    !! boundary markers, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    integer, private, allocatable, dimension(:), target :: ib_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the
    !! immersed boundary markers, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: err_code, ierr, v_size
    !> @}
    !$acc declare create(v_size)

contains

    subroutine s_initialize_mpi_proxy_module()

        if (ib) then
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE(ib_buff_send(0:-1 + gp_layers * &
                                            & (m + 2*gp_layers + 1)* &
                                            & (n + 2*gp_layers + 1)* &
                                            & (p + 2*gp_layers + 1)/ &
                                            & (min(m, n, p) + 2*gp_layers + 1)))
                else
                    @:ALLOCATE(ib_buff_send(0:-1 + gp_layers* &
                                            & (max(m, n) + 2*gp_layers + 1)))
                end if
            else
                @:ALLOCATE(ib_buff_send(0:-1 + gp_layers))
            end if
            @:ALLOCATE(ib_buff_recv(0:ubound(ib_buff_send, 1)))
        end if

    end subroutine s_initialize_mpi_proxy_module

contains

    !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.
    subroutine s_mpi_bcast_user_inputs()

#ifdef MFC_MPI

        integer :: i, j !< Generic loop iterator

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['k_x', 'k_y', 'k_z', 'w_x', 'w_y', 'w_z', 'p_x', 'p_y', &
            & 'p_z', 'g_x', 'g_y', 'g_z']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in ['t_step_old', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 't_step_start','t_step_stop','t_step_save','t_step_print',       &
            & 'model_eqns','time_stepper', 'riemann_solver', 'low_Mach',       &
            & 'wave_speeds', 'avg_state', 'precision', 'bc_x%beg', 'bc_x%end', &
            & 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',  'fd_order',     &
            & 'num_probes', 'num_integrals', 'bubble_model', 'thermal',        &
            & 'R0_type', 'num_source', 'relax_model', 'num_ibs', 'n_start',    &
            & 'num_bc_patches']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'mpp_lim',     &
            &  'mp_weno', 'rdma_mpi', 'weno_flat', 'riemann_flat', &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles_euler', 'polytropic',    &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z',    &
            & 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',          &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',          &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',          &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt', 'surface_tension',        &
            & 'viscous', 'shear_stress', 'bulk_stress', 'bubbles_lagrange',     &
            & 'hyperelasticity', 'rkck_adap_dt', 'bc_io', 'powell', 'cont_damage' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        if (chemistry) then
            #:for VAR in [ 'diffusion', 'reactions' ]
                call MPI_BCAST(chem_params%${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'gamma_method' ]
                call MPI_BCAST(chem_params%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        if (bubbles_lagrange) then
            #:for VAR in [ 'heatTransfer_model', 'massTransfer_model', 'pressure_corrector', &
                & 'write_bubbles', 'write_bubbles_stats']
                call MPI_BCAST(lag_params%${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in ['solver_approach', 'cluster_type', 'smooth_type', 'nBubs_glb']
                call MPI_BCAST(lag_params%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'c0', 'rho0', 'T0', 'x0', 'diffcoefvap', 'epsilonb','charwidth', &
                & 'valmaxvoid', 'Thost']
                call MPI_BCAST(lag_params%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end if

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref','R0ref','Web','Ca', 'sigma', &
            & 'Re_inv', 'poly_sigma', 'palpha_eps', 'ptgalpha_eps', 'pi_fac',    &
            & 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve2', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3', &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3', &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out', 'bc_z%pres_in','bc_z%pres_out', &
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',    &
            & 'z_domain%beg', 'z_domain%end', 'x_a', 'x_b', 'y_a', 'y_b', 'z_a', &
            & 'z_b', 't_stop', 't_save', 'cfl_target', 'rkck_tolerance', 'Bx0',  &
            & 'tau_star', 'cont_damage_s', 'alpha_bar' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(mapped_weno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(teno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz_q, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(mhd, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(relativity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v', 'cp_v','G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in','bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', &
                & 'x_centroid', 'y_centroid', 'c', 'm', 'p', 't', 'theta', 'slip' ]
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        do j = 1, num_probes_max
            do i = 1, 3
                call MPI_BCAST(acoustic(j)%loc(i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            end do

            call MPI_BCAST(acoustic(j)%dipole, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'pulse', 'support', 'num_elements', 'element_on', 'bb_num_freq' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'mag', 'length', 'height', &
                'wavelength', 'frequency', 'gauss_sigma_dist', 'gauss_sigma_time', &
                'npulse', 'dir', 'delay', 'foc_length', 'aperture', &
                'element_spacing_angle', 'element_polygon_ratio', 'rotate_angle', &
                'bb_bandwidth', 'bb_lowest_freq' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'x','y','z' ]
                call MPI_BCAST(probe(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax' ]
                call MPI_BCAST(integral(j)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

#endif

    end subroutine s_mpi_bcast_user_inputs

    !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
    subroutine s_mpi_sendrecv_ib_buffers(ib_markers, gp_layers)

        type(integer_field), intent(inout) :: ib_markers
        integer, intent(in) :: gp_layers

        integer :: i, j, k, l, r !< Generic loop iterators
        integer, pointer, dimension(:) :: p_i_send, p_i_recv

#ifdef MFC_MPI
        !nCalls_time = nCalls_time + 1

        ! MPI Communication in x-direction
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(n + 1)*(p + 1), &
                        MPI_INTEGER, bc_x%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
        ! END: MPI Communication in x-direction

        ! MPI Communication in y-direction

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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        MPI_INTEGER, bc_y%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
        ! END: MPI Communication in y-direction

        ! MPI Communication in z-direction
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%beg, 1, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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
                    p_i_send => ib_buff_send
                    p_i_recv => ib_buff_recv

                    !$acc data attach(p_i_send, p_i_recv)
                    !$acc host_data use_device(p_i_send, p_i_recv)

                    ! Send/receive buffer to/from bc_x%end/bc_x%beg
                    call MPI_SENDRECV( &
                        p_i_send(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 0, &
                        p_i_recv(0), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1), &
                        MPI_INTEGER, bc_z%end, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    !$acc end host_data
                    !$acc end data
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

        ! END: MPI Communication in z-direction

#endif

    end subroutine s_mpi_sendrecv_ib_buffers

    subroutine s_mpi_send_random_number(phi_rn, num_freq)
        integer, intent(in) :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn
#ifdef MFC_MPI
        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif
    end subroutine s_mpi_send_random_number

    subroutine s_finalize_mpi_proxy_module()

        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
