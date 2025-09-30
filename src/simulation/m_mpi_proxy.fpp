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

    integer, private, allocatable, dimension(:) :: ib_buff_send !<
    !! This variable is utilized to pack and send the buffer of the immersed
    !! boundary markers, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    integer, private, allocatable, dimension(:) :: ib_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the
    !! immersed boundary markers, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    integer :: i_halo_size
    $:GPU_DECLARE(create='[i_halo_size]')

contains

    subroutine s_initialize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            if (n > 0) then
                if (p > 0) then
                    i_halo_size = -1 + buff_size* &
                                            & (m + 2*buff_size + 1)* &
                                            & (n + 2*buff_size + 1)* &
                                            & (p + 2*buff_size + 1)/ &
                                            & (cells_bounds%mnp_min + 2*buff_size + 1)
                else
                    i_halo_size = -1 + buff_size* &
                                            & (cells_bounds%mn_max + 2*buff_size + 1)
                end if
            else
                i_halo_size = -1 + buff_size
            end if

            $:GPU_UPDATE(device='[i_halo_size]')
            @:ALLOCATE(ib_buff_send(0:i_halo_size), ib_buff_recv(0:i_halo_size))
        end if
#endif

    end subroutine s_initialize_mpi_proxy_module

    !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.
    impure subroutine s_mpi_bcast_user_inputs()

#ifdef MFC_MPI

        integer :: i, j !< Generic loop iterator
        integer :: ierr !< Generic flag used to identify and report MPI errors

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
            & 'num_source', 'relax_model', 'num_ibs', 'n_start',    &
            & 'num_bc_patches', 'num_igr_iters', 'num_igr_warm_start_iters', &
            & 'adap_dt_max_iters' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'mpp_lim',     &
            &  'mp_weno', 'rdma_mpi', 'powell', 'cont_damage', 'bc_io', &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles_euler', 'polytropic',    &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z',    &
            & 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',          &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',          &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',          &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt', 'surface_tension',       &
            & 'shear_stress', 'bulk_stress', 'bubbles_lagrange',                &
            & 'hyperelasticity', 'down_sample', 'int_comp' ]
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
            & 'z_b', 't_stop', 't_save', 'cfl_target', 'Bx0', 'alf_factor',  &
            & 'tau_star', 'cont_damage_s', 'alpha_bar', 'adap_dt_tol', &
            & 'ic_eps', 'ic_beta' ]
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
            call MPI_BCAST(igr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_pres_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_iter_solver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(viscous, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(recon_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_lim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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

        ! NVIDIA UVM variables
        call MPI_BCAST(nv_uvm_out_of_core, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_igr_temps_on_gpu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_pref_gpu, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

#endif

    end subroutine s_mpi_bcast_user_inputs

    subroutine s_mpi_sendrecv_ib_buffers(ib_markers, mpi_dir, pbc_loc)

        type(integer_field), intent(inout) :: ib_markers

        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: i, j, k, l, r, q !< Generic loop iterators

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        logical :: beg_end_geq_0, qbmm_comm

        integer :: pack_offset, unpack_offset

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors

        call nvtxStartRange("IB-MARKER-COMM-PACKBUF")

        buffer_counts = (/ &
                        buff_size*(n + 1)*(p + 1), &
                        buff_size*(m + 2*buff_size + 1)*(p + 1), &
                        buff_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
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
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, buff_size - 1
                                r = (j + buff_size*(k + (n + 1)*l))
                                ib_buff_send(r) = ib_markers%sf(j + pack_offset, k, l)
                            end do
                        end do
                    end do
                #:elif mpi_dir == 2
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = 0, p
                        do k = 0, buff_size - 1
                            do j = -buff_size, m + buff_size
                                r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                     (k + buff_size*l))
                                ib_buff_send(r) = ib_markers%sf(j, k + pack_offset, l)
                            end do
                        end do
                    end do
                #:else
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = 0, buff_size - 1
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                     ((k + buff_size) + (n + 2*buff_size + 1)*l))
                                ib_buff_send(r) = ib_markers%sf(j, k, l + pack_offset)
                            end do
                        end do
                    end do
                #:endif
            end if
        #:endfor
        call nvtxEndRange ! Packbuf

        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    #:call GPU_HOST_DATA(use_device='[ib_buff_send, ib_buff_recv]')

                        call nvtxStartRange("IB-MARKER-SENDRECV-RDMA")
                        call MPI_SENDRECV( &
                            ib_buff_send, buffer_count, MPI_INTEGER, dst_proc, send_tag, &
                            ib_buff_recv, buffer_count, MPI_INTEGER, src_proc, recv_tag, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                        call nvtxEndRange

                    #:endcall GPU_HOST_DATA
                    $:GPU_WAIT()
                #:else
                    call nvtxStartRange("IB-MARKER-DEV2HOST")
                    $:GPU_UPDATE(host='[ib_buff_send]')
                    call nvtxEndRange

                    call nvtxStartRange("IB-MARKER-SENDRECV-NO-RMDA")
                    call MPI_SENDRECV( &
                        ib_buff_send, buffer_count, MPI_INTEGER, dst_proc, send_tag, &
                        ib_buff_recv, buffer_count, MPI_INTEGER, src_proc, recv_tag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                    call nvtxEndRange

                    call nvtxStartRange("IB-MARKER-HOST2DEV")
                    $:GPU_UPDATE(device='[ib_buff_recv]')
                    call nvtxEndRange
                #:endif
            end if
        #:endfor

        ! Unpack Received Buffer
        call nvtxStartRange("IB-MARKER-COMM-UNPACKBUF")
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = 0, p
                        do k = 0, n
                            do j = -buff_size, -1
                                r = (j + buff_size*((k + 1) + (n + 1)*l))
                                ib_markers%sf(j + unpack_offset, k, l) = ib_buff_recv(r)
                            end do
                        end do
                    end do
                #:elif mpi_dir == 2
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = 0, p
                        do k = -buff_size, -1
                            do j = -buff_size, m + buff_size
                                r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                     ((k + buff_size) + buff_size*l))
                                ib_markers%sf(j, k + unpack_offset, l) = ib_buff_recv(r)
                            end do
                        end do
                    end do
                #:else
                    ! Unpacking buffer from bc_z%beg
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[r]')
                    do l = -buff_size, -1
                        do k = -buff_size, n + buff_size
                            do j = -buff_size, m + buff_size
                                r = ((j + buff_size) + (m + 2*buff_size + 1)* &
                                     ((k + buff_size) + (n + 2*buff_size + 1)* &
                                      (l + buff_size)))
                                ib_markers%sf(j, k, l + unpack_offset) = ib_buff_recv(r)
                            end do
                        end do
                    end do
                #:endif
            end if
        #:endfor
        call nvtxEndRange
#endif

    end subroutine s_mpi_sendrecv_ib_buffers

    impure subroutine s_mpi_send_random_number(phi_rn, num_freq)
        integer, intent(in) :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors
        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_send_random_number

    subroutine s_finalize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
