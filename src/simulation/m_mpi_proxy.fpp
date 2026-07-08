!>
!! @file
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief MPI halo exchange, domain decomposition, and buffer packing/unpacking for the simulation solver
module m_mpi_proxy

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_helper_basic
    use m_helper
    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_nvtx
    use ieee_arithmetic

    implicit none

    integer, private, allocatable, dimension(:) :: ib_buff_send  !< IB marker send buffer for halo exchange
    integer, private, allocatable, dimension(:) :: ib_buff_recv  !< IB marker receive buffer for halo exchange
    integer                                     :: i_halo_size
    $:GPU_DECLARE(create='[i_halo_size]')

    real(wp), private, allocatable, dimension(:) :: amr_buff_send  !< AMR fine-halo send buffer (device-resident)
    real(wp), private, allocatable, dimension(:) :: amr_buff_recv  !< AMR fine-halo receive buffer (device-resident)

contains

    !> Initialize the MPI proxy module
    subroutine s_initialize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            if (n > 0) then
                if (p > 0) then
                    i_halo_size = -1 + buff_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)*(p + 2*buff_size + 1) &
                                                  & /(cells_bounds%mnp_min + 2*buff_size + 1)
                else
                    i_halo_size = -1 + buff_size*(cells_bounds%mn_max + 2*buff_size + 1)
                end if
            else
                i_halo_size = -1 + buff_size
            end if

            $:GPU_UPDATE(device='[i_halo_size]')
            @:ALLOCATE(ib_buff_send(0:i_halo_size), ib_buff_recv(0:i_halo_size))
        end if
#endif

    end subroutine s_initialize_mpi_proxy_module

    !> Preallocate the AMR fine-halo pack buffers at this rank's max fine extents (m_amr's preallocation cap). Called from
    !! s_initialize_amr_module on all ranks; no-op without MPI or at np=1.
    impure subroutine s_initialize_amr_mpi_buffers(max_f1, max_f2, max_f3)

        integer, intent(in) :: max_f1, max_f2, max_f3

#ifdef MFC_MPI
        integer :: sz

        if (num_procs == 1) return
        sz = sys_size*buff_size*(max_f2 + 1)*(max_f3 + 1)
        if (n_glb > 0) sz = max(sz, sys_size*buff_size*(max_f1 + 2*buff_size + 1)*(max_f3 + 1))
        if (p_glb > 0) sz = max(sz, sys_size*buff_size*(max_f1 + 2*buff_size + 1)*(max_f2 + 2*buff_size + 1))
        @:ALLOCATE(amr_buff_send(0:sz - 1), amr_buff_recv(0:sz - 1))
#endif

    end subroutine s_initialize_amr_mpi_buffers

    !> AMR fine-level halo exchange: overwrite the buff_size fine ghost layers of q_fine at CONTINUATION faces (where the block
    !! extends past this rank's intersection) with the neighbor rank's true fine data - the coarse-topology neighbor (bc_x/y/z rank
    !! encoding) holds matching fine cells by construction (mirror decomposition, lockstep in time). Sequential per direction with
    !! the coarse halo's transverse ranges, so corners propagate. Pack/unpack run as device kernels; only the buffers move through
    !! the host (no-op macros on CPU builds). Reads the COARSE grid state for participation - call BEFORE s_amr_swap_to_fine. No
    !! exchange fires at np=1 or for fully-contained blocks; fm/fn/fp are the LOCAL fine extents.
    impure subroutine s_mpi_sendrecv_amr_fine_halo(q_fine, fm, fn, fp)

        type(scalar_field), dimension(1:), intent(inout) :: q_fine
        integer, intent(in)                              :: fm, fn, fp

#ifdef MFC_MPI
        integer :: i, j, k, l, r, cnt, nbr, stag, rtag, pack_off, unpack_off, d, loc, ierr
        logical :: go

        if (num_procs == 1) return
        if (.not. amr_rank_owns_block) return
        do d = 1, num_dims
            do loc = -1, 1, 2
                ! a face participates iff the block CONTINUES past this rank's intersection there; the
                ! neighbor then participates on its matching face by construction (blocking pairwise
                ! sendrecvs in a fixed (d, loc) order cannot deadlock)
                select case (d)
                case (1)
                    cnt = sys_size*buff_size*(fn + 1)*(fp + 1)
                    if (loc == -1) then
                        go = amr_isect_lo(1) > amr_region_lo(1); nbr = bc_x%beg
                        pack_off = 0; unpack_off = -buff_size
                    else
                        go = amr_isect_hi(1) < amr_region_hi(1); nbr = bc_x%end
                        pack_off = fm - buff_size + 1; unpack_off = fm + 1
                    end if
                case (2)
                    cnt = sys_size*buff_size*(fm + 2*buff_size + 1)*(fp + 1)
                    if (loc == -1) then
                        go = amr_isect_lo(2) > amr_region_lo(2); nbr = bc_y%beg
                        pack_off = 0; unpack_off = -buff_size
                    else
                        go = amr_isect_hi(2) < amr_region_hi(2); nbr = bc_y%end
                        pack_off = fn - buff_size + 1; unpack_off = fn + 1
                    end if
                case (3)
                    cnt = sys_size*buff_size*(fm + 2*buff_size + 1)*(fn + 2*buff_size + 1)
                    if (loc == -1) then
                        go = amr_isect_lo(3) > amr_region_lo(3); nbr = bc_z%beg
                        pack_off = 0; unpack_off = -buff_size
                    else
                        go = amr_isect_hi(3) < amr_region_hi(3); nbr = bc_z%end
                        pack_off = fp - buff_size + 1; unpack_off = fp + 1
                    end if
                end select
                if (.not. go) cycle
                ! tags by data direction: 2*d = moving to the lower rank, 2*d+1 = moving to the upper rank
                if (loc == -1) then
                    stag = 2*d; rtag = 2*d + 1
                else
                    stag = 2*d + 1; rtag = 2*d
                end if
                ! pack the near-face INTERIOR layers (device kernel)
                #:for DIR in [1, 2, 3]
                    if (d == ${DIR}$) then
                        #:if DIR == 1
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, fp
                                do k = 0, fn
                                    do j = 0, buff_size - 1
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*(j + buff_size*(k + (fn + 1)*l))
                                            amr_buff_send(r) = real(q_fine(i)%sf(j + pack_off, k, l), wp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:elif DIR == 2
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, fp
                                do k = 0, buff_size - 1
                                    do j = -buff_size, fm + buff_size
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*((j + buff_size) + (fm + 2*buff_size + 1)*(k + buff_size*l))
                                            amr_buff_send(r) = real(q_fine(i)%sf(j, k + pack_off, l), wp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:else
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, buff_size - 1
                                do k = -buff_size, fn + buff_size
                                    do j = -buff_size, fm + buff_size
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*((j + buff_size) + (fm + 2*buff_size + 1)*((k + buff_size) &
                                                 & + (fn + 2*buff_size + 1)*l))
                                            amr_buff_send(r) = real(q_fine(i)%sf(j, k, l + pack_off), wp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:endif
                    end if
                #:endfor
                $:GPU_UPDATE(host='[amr_buff_send]')
                call MPI_SENDRECV(amr_buff_send, cnt, mpi_p, nbr, stag, amr_buff_recv, cnt, mpi_p, nbr, rtag, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                $:GPU_UPDATE(device='[amr_buff_recv]')
                ! unpack into the near-face GHOST layers (device kernel)
                #:for DIR in [1, 2, 3]
                    if (d == ${DIR}$) then
                        #:if DIR == 1
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, fp
                                do k = 0, fn
                                    do j = 0, buff_size - 1
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*(j + buff_size*(k + (fn + 1)*l))
                                            q_fine(i)%sf(j + unpack_off, k, l) = real(amr_buff_recv(r), stp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:elif DIR == 2
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, fp
                                do k = 0, buff_size - 1
                                    do j = -buff_size, fm + buff_size
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*((j + buff_size) + (fm + 2*buff_size + 1)*(k + buff_size*l))
                                            q_fine(i)%sf(j, k + unpack_off, l) = real(amr_buff_recv(r), stp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:else
                            $:GPU_PARALLEL_LOOP(collapse=4, private='[r]')
                            do l = 0, buff_size - 1
                                do k = -buff_size, fn + buff_size
                                    do j = -buff_size, fm + buff_size
                                        do i = 1, sys_size
                                            r = (i - 1) + sys_size*((j + buff_size) + (fm + 2*buff_size + 1)*((k + buff_size) &
                                                 & + (fn + 2*buff_size + 1)*l))
                                            q_fine(i)%sf(j, k, l + unpack_off) = real(amr_buff_recv(r), stp)
                                        end do
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:endif
                    end if
                #:endfor
            end do
        end do
#endif

    end subroutine s_mpi_sendrecv_amr_fine_halo

    !> Since only the processor with rank 0 reads and verifies the consistency of user inputs, these are initially not available to
    !! the other processors. Then, the purpose of this subroutine is to distribute the user inputs to the remaining processors in
    !! the communicator.
    impure subroutine s_mpi_bcast_user_inputs()

#ifdef MFC_MPI
        integer :: i, j  !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        ! Generated: case_dir, namelist scalars (INT/LOG/REAL), CASE_OPT guard, fluid_pp loop,
        !            bub_pp guard, lag_params guard, chem_params guard
        #:include 'generated_bcast.fpp'

        ! manual: m_glb, n_glb, p_glb (computed in s_read_input_file, not namelist-bound)
        call MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x/y/z member broadcasts (struct members not in NAMELIST_VARS)
        #:for VAR in [ 'bc_x%beg', 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end']
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'bc_x%grcbc_in', 'bc_x%grcbc_out', 'bc_x%grcbc_vel_out',  &
            & 'bc_y%grcbc_in', 'bc_y%grcbc_out', 'bc_y%grcbc_vel_out',            &
            & 'bc_z%grcbc_in', 'bc_z%grcbc_out', 'bc_z%grcbc_vel_out',            &
            & 'bc_x%isothermal_in', 'bc_y%isothermal_in', 'bc_z%isothermal_in',   &
            & 'bc_x%isothermal_out', 'bc_y%isothermal_out', 'bc_z%isothermal_out']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve3', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3',           &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3',           &
            & 'bc_x%pres_in','bc_x%pres_out','bc_y%pres_in','bc_y%pres_out',               &
            & 'bc_z%pres_in','bc_z%pres_out',                                               &
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',              &
            & 'z_domain%beg', 'z_domain%end',                                               &
            & 'bc_x%Twall_in', 'bc_x%Twall_out', 'bc_y%Twall_in', 'bc_y%Twall_out',       &
            & 'bc_z%Twall_in', 'bc_z%Twall_out']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, 3
            #:for VAR in [ 'bc_x%vel_in', 'bc_x%vel_out', 'bc_y%vel_in', 'bc_y%vel_out',  &
                & 'bc_z%vel_in', 'bc_z%vel_out' ]
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: cfl_dt (runtime-computed logical), bc_io (BC-file existence)
        call MPI_BCAST(cfl_dt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bc_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: shear_stress, bulk_stress (derived from Re_size post-init on all ranks),
        !         bodyForces (derived from bf_x/y/z)
        call MPI_BCAST(shear_stress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bulk_stress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(bodyForces, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        ! manual: bc_x per-fluid inflow arrays (loop over num_fluids_max)
        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in', &
                & 'bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: patch_ib (sim member subset differs from pre; uses count=3, adds mass/moving_ibm)
        do i = 1, num_ibs
            #:for VAR in [ 'radius', 'length_x', 'length_y', 'length_z', &
                & 'x_centroid', 'y_centroid', 'z_centroid', 'slip', 'mass']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            #:for VAR in ['vel', 'angular_vel', 'angles']
                call MPI_BCAST(patch_ib(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(patch_ib(i)%geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%airfoil_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%model_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        ! manual: ib_airfoil (kept manual alongside patch_ib)
        do i = 1, num_ib_airfoils_max
            #:for VAR in ['c', 'p', 't', 'm']
                call MPI_BCAST(ib_airfoil(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: stl_models loop (num_stl_models scalar is generated; grouped array members)
        do i = 1, num_stl_models_max
            call MPI_BCAST(stl_models(i)%model_filepath, len(stl_models(i)%model_filepath), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(stl_models(i)%model_threshold, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:for VAR in ['model_translate', 'model_scale']
                call MPI_BCAST(stl_models(i)%${VAR}$, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! manual: particle_cloud (runtime loop to num_particle_clouds; irregular member subset)
        do i = 1, num_particle_clouds
            #:for VAR in ['x_centroid', 'y_centroid', 'z_centroid', 'length_x', 'length_y', 'length_z', &
                & 'radius', 'mass', 'min_spacing']
                call MPI_BCAST(particle_cloud(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(particle_cloud(i)%num_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%moving_ibm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(particle_cloud(i)%packing_method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        end do

        ! manual: acoustic/probe/integral (combined loop; complex acoustic member set)
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

    !> Broadcast random phase numbers from rank 0 to all MPI processes
    impure subroutine s_mpi_send_random_number(phi_rn, num_freq)

        integer, intent(in)                            :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_send_random_number

    !> Finalize the MPI proxy module
    subroutine s_finalize_mpi_proxy_module()

#ifdef MFC_MPI
        if (ib) then
            @:DEALLOCATE(ib_buff_send, ib_buff_recv)
        end if
        if (allocated(amr_buff_send)) then
            @:DEALLOCATE(amr_buff_send, amr_buff_recv)
        end if
#endif

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
