#:include 'macros.fpp'

module m_boundary_conditions_common

    use m_global_parameters

    use m_mpi_proxy

    use m_helper

    use m_compile_specific

    implicit none

    private; public :: s_initialize_boundary_conditions_module, &
 s_generate_boundary_condition_patch_buffers, &
 bc_id_sfs, bc_id_has_bc, &
 s_read_boundary_condition_files, &
 s_write_boundary_condition_files

    ! Boundary condition structure. bc_id_sfs(<dir>, <loc>)%sf(<x>, <y>, <z>)
    ! holds the necessary boundary condition information for a point along the
    ! surface in the direction <dir> on the side <loc>, to apply said condition.
    type(t_bc_id_sf), dimension(1:3, -1:1) :: bc_id_sfs

#ifdef MFC_MPI
    integer :: mpi_bc_id_type
    integer, dimension(1:3, -1:1) :: mpi_bc_sf_type
#endif

    ! Optimization structure. Holds whether a boundary condition can be found
    ! on a surface. bc_id_has_bc(<dir>, <loc>, 0) is true if MPI communication
    ! is required for at least one point on the surface. In general,
    ! bc_id_has_bc(<dir>, <loc>, <id>) is true if a boundary condition with type
    ! <id> can be found on the surface on the side <loc> of the direction <dir>.
    logical, dimension(1:3, -1:1, -num_bcs_max:0) :: bc_id_has_bc

!$acc declare create(bc_id_sfs)

#ifndef MFC_PRE_PROCESS
    public :: s_populate_prim_buffers, s_populate_cons_buffers
#endif

contains

    subroutine s_initialize_boundary_conditions_module()

        integer :: iter_loc

        do iter_loc = -1, 2, 2
            @:ALLOCATE(bc_id_sfs(1, iter_loc)%sf(0:0, 0:n, 0:p))
            if (num_dims > 1) then
                @:ALLOCATE(bc_id_sfs(2, iter_loc)%sf(0:m, 0:0, 0:p))

                if (num_dims > 2) then
                    @:ALLOCATE(bc_id_sfs(3, iter_loc)%sf(0:m, 0:n, 0:0))
                end if
            end if
        end do

#ifdef MFC_MPI
        call s_create_mpi_types()
#endif

    end subroutine s_initialize_boundary_conditions_module

#ifdef MFC_MPI
    subroutine s_create_mpi_types()

        use mpi

        integer :: blocklengths(2) = (/1, 3/)
        integer(KIND=MPI_ADDRESS_KIND) :: displacements(2) = (/0, 8/)   ! Offset between fields in memory
        integer :: types(2) = (/MPI_INTEGER, MPI_DOUBLE_PRECISION/)
        integer, dimension(1:3) :: sf_extents_loc, sf_start_idx
        integer :: ierr

        integer :: iter_dir, iter_loc

        call MPI_Type_create_struct(2, blocklengths, displacements, types, mpi_bc_id_type, ierr)
        call MPI_Type_commit(mpi_bc_id_type, ierr)

        do iter_dir = 1, num_dims
            do iter_loc = -1, 1, 2

                sf_start_idx = (/0, 0, 0/)
                sf_extents_loc = shape(bc_id_sfs(iter_dir, iter_loc)%sf)

                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sf_extents_loc, sf_extents_loc, sf_start_idx, &
                                              MPI_ORDER_FORTRAN, mpi_bc_id_type, mpi_bc_sf_type(iter_dir, iter_loc), ierr)
                call MPI_TYPE_COMMIT(mpi_bc_sf_type(iter_dir, iter_loc), ierr)

            end do
        end do

    end subroutine s_create_mpi_types
#endif

    subroutine s_write_boundary_condition_files(step_dirpath)

        character(LEN=*), intent(in) :: step_dirpath

        integer :: iter_dir, iter_loc
        character(len=path_len) :: file_path

        character(len=10) :: status

#ifdef MFC_MPI
        integer :: ierr
        integer :: file_id
        integer :: offset
        character(len=7) :: proc_rank_str
#endif

#ifdef MFC_PRE_PROCESS
        if (old_grid) then
            status = 'old'
        else
            status = 'new'
        end if
#else
        status = 'unknown'
#endif

        if (parallel_io .eqv. .false.) then
            file_path = trim(step_dirpath)//'/bc.dat'
            open (1, FILE=trim(file_path), FORM='unformatted', STATUS=status)
            do iter_dir = 1, num_dims
                do iter_loc = -1, 1, 2
                    write (1) bc_id_sfs(iter_dir, iter_loc)%sf
                end do
            end do
            close (1)
        else
#ifdef MFC_MPI
            write (proc_rank_str, '(I7.7)') proc_rank
            file_path = trim(step_dirpath)//'/bc_'//trim(proc_rank_str)//'.dat'
            call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_id, ierr)

            offset = 0
            do iter_dir = 1, num_dims
                do iter_loc = -1, 1, 2
                    call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), mpi_bc_id_type, mpi_bc_sf_type(iter_dir, iter_loc), 'native', MPI_INFO_NULL, ierr)
                    call MPI_File_write_all(file_id, bc_id_sfs(iter_dir, iter_loc)%sf, 1, mpi_bc_sf_type(iter_dir, iter_loc), MPI_STATUS_IGNORE, ierr)
                    offset = offset + sizeof(bc_id_sfs(iter_dir, iter_loc)%sf)
                end do
            end do

            call MPI_File_close(file_id, ierr)
#endif
        end if

    end subroutine s_write_boundary_condition_files

    subroutine s_read_boundary_condition_files(step_dirpath)

        character(LEN=*), intent(in) :: step_dirpath

        integer :: iter_dir, iter_loc
        logical :: file_exist
        character(len=path_len) :: file_path

#ifdef MFC_MPI
        integer :: ierr
        integer :: file_id
        integer :: offset
        character(len=7) :: proc_rank_str
#endif

        if (parallel_io .eqv. .false.) then
            file_path = trim(step_dirpath)//'/bc.dat'
        else
#ifdef MFC_MPI
            write (proc_rank_str, '(I7.7)') proc_rank
            file_path = trim(step_dirpath)//'/bc_'//trim(proc_rank_str)//'.dat'
#endif
        end if

        inquire (FILE=trim(file_path), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

        if (parallel_io .eqv. .false.) then
            open (1, FILE=trim(file_path), FORM='unformatted', STATUS='unknown')
            do iter_dir = 1, num_dims
                do iter_loc = -1, 1, 2
                    read (1) bc_id_sfs(iter_dir, iter_loc)%sf
                end do
            end do
            close (1)
        else
#ifdef MFC_MPI
            call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_RDONLY, MPI_INFO_NULL, file_id, ierr)

            offset = 0
            do iter_dir = 1, num_dims
                do iter_loc = -1, 1, 2
                    call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), mpi_bc_id_type, mpi_bc_sf_type(iter_dir, iter_loc), 'native', MPI_INFO_NULL, ierr)
                    call MPI_File_read_all(file_id, bc_id_sfs(iter_dir, iter_loc)%sf, 1, mpi_bc_sf_type(iter_dir, iter_loc), MPI_STATUS_IGNORE, ierr)
                    offset = offset + sizeof(bc_id_sfs(iter_dir, iter_loc)%sf)
                end do
            end do

            call MPI_File_close(file_id, ierr)
#endif
        end if

        call s_generate_boundary_condition_lookup_buffers()

    end subroutine s_read_boundary_condition_files

    subroutine s_generate_boundary_condition_patch_buffers()

        type(bc_patch_parameters) :: bc
        integer :: iter_dir, iter_loc
        real(kind(0d0)) :: radius2

        type(int_bounds_info), dimension(1:3) :: user_input_bcs

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        user_input_bcs = [bc_x, bc_y, bc_z]

        do iter_dir = 1, num_dims
            #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="-1", thickness=1, gpu=False)
                bc_id_sfs(iter_dir, -1)%sf(exlhs, eylhs, ezlhs)%type = user_input_bcs(iter_dir)%beg
                bc_id_sfs(iter_dir, -1)%sf(exlhs, eylhs, ezlhs)%vel = user_input_bcs(iter_dir)%vel_beg
            #:endblock
            #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="+1", thickness=1, gpu=False)
                bc_id_sfs(iter_dir, +1)%sf(exlhs, eylhs, ezlhs)%type = user_input_bcs(iter_dir)%end
                bc_id_sfs(iter_dir, +1)%sf(exlhs, eylhs, ezlhs)%vel = user_input_bcs(iter_dir)%vel_end
            #:endblock
        end do

        do iter_dir = 1, num_dims
            do iter_loc = -1, 1, 2
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc", thickness=1, gpu=False)
                    do i = 1, num_bc_patches
                        bc = patch_bc(i)
                        if (bc%geometry == 1) then     ! Cuboid
                            #:for dir, name in [(1, "x"), (2, "y"), (3, "z")]
                                if (${dir}$ /= bc%dir .and. ${dir}$ <= num_dims) then
                                    if (${name}$_cc(e${name}$lhs) >= bc%centroid(${dir}$) + 0.5d0*bc%length(${dir}$)) then
                                        cycle
                                    end if
                                    if (${name}$_cc(e${name}$lhs) <= bc%centroid(${dir}$) - 0.5d0*bc%length(${dir}$)) then
                                        cycle
                                    end if
                                end if
                            #:endfor
                        elseif (bc%geometry == 2) then ! Spheroid
                            radius2 = 0d0
                            #:for dir, name in [(1, "x"), (2, "y"), (3, "z")]
                                if (${dir}$ /= bc%dir .and. ${dir}$ <= num_dims) then
                                    radius2 = radius2 + (${name}$_cc(e${name}$lhs) - bc%centroid(${dir}$))**2
                                end if
                            #:endfor
                            if (radius2 > bc%radius**2) then
                                cycle
                            end if
                        end if

                        bc_id_sfs(bc%dir, bc%loc)%sf(exlhs, eylhs, ezlhs)%type = bc%type
                        bc_id_sfs(bc%dir, bc%loc)%sf(exlhs, eylhs, ezlhs)%vel = bc%vel
                    end do
                #:endblock
            end do
        end do

        do iter_dir = 1, num_dims
            do iter_loc = -1, 1, 2
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc", thickness=1, gpu=False)
                    if (proc_nums(iter_dir) > 1 .and. ( &
                        bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type == -1 &
                        .or. (iter_loc == -1 .and. proc_coords(iter_dir) > 0) &
                        .or. (iter_loc == +1 .and. proc_coords(iter_dir) < proc_nums(iter_dir) - 1) &
                        )) then
                        bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type = neighbor_procs(iter_dir, iter_loc)
                    end if
                #:endblock

                !$acc update device(bc_id_sfs(iter_dir, iter_loc)%sf)
            end do
        end do

        call s_generate_boundary_condition_lookup_buffers()

    end subroutine s_generate_boundary_condition_patch_buffers

    subroutine s_generate_boundary_condition_lookup_buffers

        integer :: iter_dir, iter_loc

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        do iter_dir = 1, num_dims
            do iter_loc = -1, 1, 2
                bc_id_has_bc(iter_dir, iter_loc, :) = .false.
                #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc", thickness=1, gpu=False)
                    bc_id_has_bc(iter_dir, iter_loc, min(0, bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type)) = .true.
                #:endblock
            end do
        end do

    end subroutine s_generate_boundary_condition_lookup_buffers

#ifndef MFC_PRE_PROCESS
    !>  The purpose of this procedure is to populate the buffers
    !!      of the primitive variables, depending on the selected
    !!      boundary conditions.
    subroutine s_populate_prim_buffers(q_prim_vf &
#ifdef MFC_SIMULATION
                                       , pb, mv &
#endif
                                       )

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_SIMULATION
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif

        integer :: iter_dir, iter_loc

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        do iter_dir = 1, num_dims
            do iter_loc = -1, 1, 2
                if (any(bc_id_has_bc(iter_dir, iter_loc, :-1))) then
                    #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc")
                        select case (bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type)
                        case (-13:-3); ! Ghost-cell extrap.
                            !$acc loop seq
                            do i = 1, sys_size
                                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
                            end do
                        case (-2); ! Symmetry
                            !$acc loop seq
                            do i = 1, momxb + iter_dir - 2
                                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(sx, sy, sz)
                            end do

                            q_prim_vf(momxb + iter_dir - 1)%sf(x, y, z) = &
                                -q_prim_vf(momxb + iter_dir - 1)%sf(sx, sy, sz)

                            !$acc loop seq
                            do i = momxb + iter_dir, sys_size
                                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(sx, sy, sz)
                            end do
                        case (-1); ! Periodic
                            !$acc loop seq
                            do i = 1, sys_size
                                q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(px, py, pz)
                            end do
                        case (-15); ! Splip Wall
                            !$acc loop seq
                            do i = 1, sys_size
                                if (i == momxb + iter_dir - 1) then
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        -q_prim_vf(i)%sf(sx, sy, sz) + 2d0*bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%vel(iter_dir)
                                else
                                    q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
                                end if
                            end do
                        case (-16); ! No Slip Wall
                            !$acc loop seq
                            do i = 1, sys_size
                                if (i >= momxb .and. i <= momxe) then
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        -q_prim_vf(i)%sf(sx, sy, sz) + 2d0*bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%vel(i - momxb + 1)
                                else
                                    q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
                                end if
                            end do
                        case (-17); ! Velocity Inlet
                            !$acc loop seq
                            do i = 1, sys_size
                                if (i >= momxb .and. i <= momxe) then
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%vel(i - momxb + 1)
                                else
                                    q_prim_vf(i)%sf(x, y, z) = q_prim_vf(i)%sf(ex, ey, ez)
                                end if
                            end do
                        case (-14); ! Axis
                            if (z_cc(k) < pi) then
                                !$acc loop seq
                                do i = 1, momxb
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        q_prim_vf(i)%sf(sx, sy, z + ((p + 1)/2))
                                end do

                                q_prim_vf(momxb + 1)%sf(x, y, z) = &
                                    -q_prim_vf(momxb + 1)%sf(sx, sy, z + ((p + 1)/2))

                                q_prim_vf(momxe)%sf(x, y, z) = &
                                    -q_prim_vf(momxe)%sf(sx, sy, z + ((p + 1)/2))

                                !$acc loop seq
                                do i = E_idx, sys_size
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        q_prim_vf(i)%sf(sx, sy, z + ((p + 1)/2))
                                end do
                            else
                                !$acc loop seq
                                do i = 1, momxb
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        q_prim_vf(i)%sf(sx, sy, z - ((p + 1)/2))
                                end do

                                q_prim_vf(momxb + 1)%sf(x, y, z) = &
                                    -q_prim_vf(momxb + 1)%sf(sx, sy, z - ((p + 1)/2))

                                q_prim_vf(momxe)%sf(x, y, z) = &
                                    -q_prim_vf(momxe)%sf(sx, sy, z - ((p + 1)/2))

                                !$acc loop seq
                                do i = E_idx, sys_size
                                    q_prim_vf(i)%sf(x, y, z) = &
                                        q_prim_vf(i)%sf(sx, sy, z - ((p + 1)/2))
                                end do
                            end if

                        end select
                    #:endblock

#ifdef MFC_SIMULATION
                    if (qbmm .and. .not. polytropic) then
                        #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir="iter_dir", loc="iter_loc", outer_loops=[("i", 1, "nb"), ("q", 1, "nnode")])
                            select case (bc_id_sfs(iter_dir, iter_loc)%sf(exlhs, eylhs, ezlhs)%type)
                            case (-13:-3); ! Ghost-cell extrap.
                                pb(x, y, z, q, i) = pb(ex, ey, ez, q, i)
                                mv(x, y, z, q, i) = mv(ex, ey, ez, q, i)
                            case (-2); ! Symmetry
                                pb(x, y, z, q, i) = pb(sx, sy, sz, q, i)
                                mv(x, y, z, q, i) = mv(sx, sy, sz, q, i)
                            case (-1); ! Periodic
                                pb(x, y, z, q, i) = pb(px, py, pz, q, i)
                                mv(x, y, z, q, i) = mv(px, py, pz, q, i)
                            case (-14); ! Axis
                                pb(x, y, z, q, i) = pb(sx, sy, z - ((p + 1)/2), q, i)
                                mv(x, y, z, q, i) = mv(sx, sy, z - ((p + 1)/2), q, i)
                            end select
                        #:endblock
                    end if
#endif
                end if

                call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_id_sfs, &
#ifdef MFC_SIMULATION
                                                      pb, mv, &
#endif
                                                      iter_dir, iter_loc, &
                                                      bc_id_has_bc)

            end do
        end do

    end subroutine s_populate_prim_buffers

    !>  The purpose of this procedure is to populate the buffers
    !!      of the conservative variables, depending on the selected
    !!      boundary conditions.
    subroutine s_populate_cons_buffers(q_cons_vf &
#ifdef MFC_SIMULATION
                                       , pb, mv &
#endif
                                       )

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
#ifdef MFC_SIMULATION
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(inout) :: pb, mv
#endif

        call s_populate_prim_buffers(q_cons_vf &
#ifdef MFC_SIMULATION
                                     , pb, mv &
#endif
                                     )

    end subroutine s_populate_cons_buffers
#endif

end module m_boundary_conditions_common
