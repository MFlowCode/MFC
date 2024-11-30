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

    use m_nvtx

    use ieee_arithmetic

    use m_constants

    ! ==========================================================================

    implicit none

    real(kind(0d0)), private, allocatable, dimension(:), target :: c_divs_buff_send !<
    !! c_divs_buff_send is utilized to send and unpack the buffer of the cell-
    !! centered color function derivatives, for a single computational domain
    !! boundary at the time, to the the relevant neighboring processor

    real(kind(0d0)), private, allocatable, dimension(:), target :: c_divs_buff_recv
    !! c_divs_buff_recv is utilized to receiver and unpack the buffer of the cell-
    !! centered color function derivatives, for a single computational domain
    !! boundary at the time, from the relevant neighboring processor

    integer, private, allocatable, dimension(:), target :: ib_buff_send !<
    !! This variable is utilized to pack and send the buffer of the immersed
    !! boundary markers, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.

    integer, private, allocatable, dimension(:), target :: ib_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the
    !! immersed boundary markers, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.

    !$acc declare create(ib_buff_send, ib_buff_recv)
    !$acc declare create(c_divs_buff_send, c_divs_buff_recv)

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    integer, private :: err_code, ierr
    !> @}

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

        if (surface_tension) then
            nVars = num_dims + 1
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)* &
                                             & (m + 2*buff_size + 1)* &
                                             & (n + 2*buff_size + 1)* &
                                             & (p + 2*buff_size + 1)/ &
                                             & (min(m, n, p) + 2*buff_size + 1)))
                else
                    @:ALLOCATE(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)* &
                                             & (max(m, n) + 2*buff_size + 1)))
                end if
            else
                @:ALLOCATE(c_divs_buff_send(0:-1 + buff_size*(num_dims+1)))
            end if

            @:ALLOCATE(c_divs_buff_recv(0:ubound(c_divs_buff_send, 1)))
        end if
        !$acc update device(nVars)

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
            & 'R0_type', 'num_source', 'relax_model', 'num_ibs', 'n_start'       ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info','cyl_coord', 'mpp_lim',     &
            &  'mp_weno', 'rdma_mpi', 'weno_flat', 'riemann_flat', &
            & 'weno_Re_flux', 'alt_soundspeed', 'null_weights', 'mixture_err',   &
            & 'parallel_io', 'hypoelasticity', 'bubbles', 'polytropic',          &
            & 'polydisperse', 'qbmm', 'acoustic_source', 'probe_wrt', 'integral_wrt',   &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process', 'relax',          &
            & 'adv_n', 'adap_dt', 'ib', 'bodyForces', 'bf_x', 'bf_y', 'bf_z',    &
            & 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',          &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',          &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',          &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt', 'surface_tension',        &
            & 'viscous', 'shear_stress', 'bulk_stress' ]
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

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref','R0ref','Web','Ca', 'sigma', &
            & 'Re_inv', 'poly_sigma', 'palpha_eps', 'ptgalpha_eps', 'pi_fac',    &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out', 'bc_z%pres_in','bc_z%pres_out', &
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',    &
            & 'z_domain%beg', 'z_domain%end', 't_stop',  't_save', 'cfl_target']
            call MPI_BCAST(${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out']
                call MPI_BCAST(${VAR}$ (i), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(mapped_weno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(teno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz_q, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','mul0','ss','pv','gamma_v','M_v',  &
                & 'mu_v','k_v','G', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in','bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            #:endfor
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

            #:for VAR in [ 'pulse', 'support', 'num_elements', 'element_on', 'bb_num_freq' ]
                call MPI_BCAST(acoustic(j)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'mag', 'length', 'height', &
                'wavelength', 'frequency', 'gauss_sigma_dist', 'gauss_sigma_time', &
                'npulse', 'dir', 'delay', 'foc_length', 'aperture', &
                'element_spacing_angle', 'element_polygon_ratio', 'rotate_angle', &
                'bb_bandwidth', 'bb_lowest_freq' ]
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

    subroutine s_mpi_sendrecv_ib_buffers(ib_markers, bc_id_sfs)

        type(integer_field), intent(inout) :: ib_markers
        type(t_bc_id_sf), dimension(1:3, -1:1), intent(in) :: bc_id_sfs

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        integer :: buffer_counts(1:3)
        integer :: iter_dir, iter_loc
        integer, pointer :: p_send, p_recv

#ifdef MFC_MPI

        buffer_counts = (/ &
                        gp_layers*(n + 1)*(p + 1), &
                        gp_layers*(m + 2*gp_layers + 1)*(p + 1), &
                        gp_layers*(m + 2*gp_layers + 1)*(n + 2*gp_layers + 1) &
                        /)

        do iter_dir = 1, num_dims

            do iter_loc = -1, 1, 2

                call nvtxStartRange("RHS-COMM-PACKBUF")

                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="-iter_loc", pack_v_size='1', thickness='gp_layers', inner_loops=[("i", 1, 1)])
                    ib_buff_send(pack_idr) = ib_markers%sf(sx, sy, sz)
                #:endblock

                call nvtxEndRange ! Packbuf

                p_send => ib_buff_send(0)
                p_recv => ib_buff_recv(0)
                #:for rdma_mpi in [False, True]
                    if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                        #:if rdma_mpi
                            !$acc data attach(p_send, p_recv)
                            !$acc host_data use_device(p_send, p_recv)

                            call nvtxStartRange("RHS-COMM-SENDRECV-RDMA")
                        #:else
                            call nvtxStartRange("RHS-COMM-DEV2HOST")
                            !$acc update host(ib_buff_send)
                            call nvtxEndRange ! Dev2Host

                            call nvtxStartRange("RHS-COMM-SENDRECV-NO-RMDA")
                        #:endif

                        call MPI_SENDRECV( &
                            p_send, buffer_counts(iter_dir), MPI_INTEGER, neighbor_procs(iter_dir, -iter_loc), 0, &
                            p_recv, buffer_counts(iter_dir), MPI_INTEGER, neighbor_procs(iter_dir, +iter_loc), 0, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                        call nvtxEndRange ! RHS-MPI-SENDRECV-(NO)-RDMA

                        #:if rdma_mpi
                            !$acc end host_data
                            !$acc end data
                            !$acc wait
                        #:else
                            call nvtxStartRange("RHS-COMM-HOST2DEV")
                            !$acc update device(ib_buff_recv)
                            call nvtxEndRange ! Host2Dev
                        #:endif
                    end if
                #:endfor

                call nvtxStartRange("RHS-COMM-UNPACKBUF")

                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc", pack_v_size='1', thickness='gp_layers', inner_loops=[("i", 1, 1)])
                    if (bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type >= 0) then
                        ib_markers%sf(x, y, z) = ib_buff_recv(pack_idr)
                    end if
                #:endblock

                call nvtxEndRange ! Unpackbuf

            end do

        end do

#endif

    end subroutine s_mpi_sendrecv_ib_buffers

    subroutine s_mpi_sendrecv_capilary_variables_buffers(c_divs_vf, bc_id_sfs, mpi_dir, pbc_loc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs_vf
        type(t_bc_id_sf), dimension(1:3, -1:1), intent(in) :: bc_id_sfs

        integer, intent(in) :: mpi_dir, pbc_loc

        integer :: buffer_counts(1:3), buffer_count

        type(int_bounds_info) :: boundary_conditions(1:3)
        integer :: beg_end(1:2), grid_dims(1:3)
        integer :: dst_proc, src_proc, recv_tag, send_tag

        integer :: pack_offset, unpack_offset
        real(kind(0d0)), pointer :: p_send, p_recv

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

#ifdef MFC_MPI

        nVars = num_dims + 1
        !$acc update device(nVars)

        buffer_counts = (/ &
                        buff_size*nVars*(n + 1)*(p + 1), &
                        buff_size*nVars*(m + 2*buff_size + 1)*(p + 1), &
                        buff_size*nVars*(m + 2*buff_size + 1)*(n + 2*buff_size + 1) &
                        /)

        #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="-pbc_loc", inner_loops=[("i", 1, "nVars")], pack_v_size="nVars")
            c_divs_buff_send(pack_idr) = c_divs_vf(i)%sf(sx, sy, sz)
        #:endblock

        p_send => c_divs_buff_send(0)
        p_recv => c_divs_buff_recv(0)
        #:for rdma_mpi in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
                #:if rdma_mpi
                    !$acc data attach(p_send, p_recv)
                    !$acc host_data use_device(p_send, p_recv)

                    call nvtxStartRange("RHS-COMM-SENDRECV-RDMA")
                #:else
                    call nvtxStartRange("RHS-COMM-DEV2HOST")
                    !$acc update host(c_divs_buff_send)
                    call nvtxEndRange

                    call nvtxStartRange("RHS-COMM-SENDRECV-NO-RMDA")
                #:endif

                call MPI_SENDRECV( &
                    p_send, buffer_counts(mpi_dir), MPI_DOUBLE_PRECISION, neighbor_procs(mpi_dir, -pbc_loc), 0, &
                    p_recv, buffer_counts(mpi_dir), MPI_DOUBLE_PRECISION, neighbor_procs(mpi_dir, +pbc_loc), 0, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                call nvtxEndRange ! RHS-MPI-SENDRECV-(NO)-RDMA

                #:if rdma_mpi
                    !$acc end host_data
                    !$acc end data
                    !$acc wait
                #:else
                    call nvtxStartRange("RHS-COMM-HOST2DEV")
                    !$acc update device(c_divs_buff_recv)
                    call nvtxEndRange
                #:endif
            end if
        #:endfor

        #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="mpi_dir", loc="pbc_loc", inner_loops=[("i", 1, "nVars")], pack_v_size="nVars")
            if (bc_id_sfs(mpi_dir, pbc_loc)%sf(exlhs, eylhs, ezlhs)%type >= 0) then
                c_divs_vf(i)%sf(x, y, z) = c_divs_buff_recv(pack_idr)
            end if
        #:endblock

#endif

    end subroutine s_mpi_sendrecv_capilary_variables_buffers

    subroutine s_mpi_send_random_number(phi_rn, num_freq)
        integer, intent(in) :: num_freq
        real(kind(0d0)), intent(inout), dimension(1:num_freq) :: phi_rn
#ifdef MFC_MPI
        call MPI_BCAST(phi_rn, num_freq, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
    end subroutine s_mpi_send_random_number

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_mpi_proxy_module

#ifdef MFC_MPI

        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if

        if (surface_tension) then
            @:DEALLOCATE(c_divs_buff_send, c_divs_buff_recv)
        end if

#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
