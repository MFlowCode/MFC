!>
!! @file m_hdf5_xdmf_output.fpp
!! @brief Contains module m_hdf5_xdmf_output

!> @brief This module provides native HDF5 + XDMF output support for MFC.
!!        XDMF is an XML-based format that describes HDF5 data for visualization
!!        tools like ParaView. This format is more portable than Silo-HDF5 and
!!        works reliably with modern ParaView versions (6.0+).
module m_hdf5_xdmf_output

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_compile_specific
    use hdf5

    implicit none

    private; public :: s_initialize_hdf5_xdmf_module, &
 s_open_hdf5_file, &
 s_write_grid_to_hdf5, &
 s_write_variable_to_hdf5, &
 s_close_hdf5_file, &
 s_write_xdmf_file, &
 s_finalize_hdf5_xdmf_module

    ! HDF5 file and group identifiers
    integer(HID_T) :: h5file_id      !< HDF5 file identifier
    integer(HID_T) :: h5dspace_id    !< HDF5 dataspace identifier
    integer(HID_T) :: h5dset_id      !< HDF5 dataset identifier
    integer(HID_T) :: h5plist_id     !< HDF5 property list identifier

    ! Directory paths for output
    character(LEN=path_len + name_len) :: hdf5_dir
    character(LEN=path_len + 2*name_len) :: proc_rank_dir_hdf5
    character(LEN=path_len + 2*name_len) :: rootdir_hdf5

    ! Error flag
    integer :: h5err

    ! Track variable names for XDMF generation
    character(LEN=name_len), allocatable, dimension(:) :: stored_varnames
    integer :: num_stored_vars

contains

    !> @brief Initialize the HDF5+XDMF output module
    subroutine s_initialize_hdf5_xdmf_module()

        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc
        logical :: dir_check
        integer :: i

        ! Initialize HDF5 library
        call h5open_f(h5err)
        if (h5err /= 0) then
            call s_mpi_abort('Failed to initialize HDF5 library. Exiting.')
        end if

        ! Create HDF5 output directory structure
        hdf5_dir = trim(case_dir)//'/hdf5_xdmf'

        write (proc_rank_dir_hdf5, '(A,I0)') '/p', proc_rank
        proc_rank_dir_hdf5 = trim(hdf5_dir)//trim(proc_rank_dir_hdf5)

        file_loc = trim(proc_rank_dir_hdf5)//'/.'
        call my_inquire(file_loc, dir_check)
        if (.not. dir_check) then
            call s_create_directory(trim(proc_rank_dir_hdf5))
        end if

        ! Root process creates root directory for XDMF files
        if (proc_rank == 0) then
            rootdir_hdf5 = trim(hdf5_dir)//'/root'
            file_loc = trim(rootdir_hdf5)//'/.'
            call my_inquire(file_loc, dir_check)
            if (.not. dir_check) then
                call s_create_directory(trim(rootdir_hdf5))
            end if
        end if

        ! Allocate variable name storage for XDMF generation
        allocate (stored_varnames(100))  ! Max 100 variables
        num_stored_vars = 0

    end subroutine s_initialize_hdf5_xdmf_module

    !> @brief Open an HDF5 file for the current time step
    !! @param t_step Current time step
    subroutine s_open_hdf5_file(t_step)

        integer, intent(in) :: t_step
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc

        ! Generate file path
        write (file_loc, '(A,I0,A)') '/', t_step, '.h5'
        file_loc = trim(proc_rank_dir_hdf5)//trim(file_loc)

        ! Create HDF5 file with default properties
        call h5fcreate_f(trim(file_loc), H5F_ACC_TRUNC_F, h5file_id, h5err)
        if (h5err /= 0) then
            call s_mpi_abort('Unable to create HDF5 file '//trim(file_loc)//'. Exiting.')
        end if

        ! Reset variable counter for this timestep
        num_stored_vars = 0

    end subroutine s_open_hdf5_file

    !> @brief Write grid coordinates to HDF5 file
    subroutine s_write_grid_to_hdf5()

        integer(HSIZE_T), dimension(1) :: dims_1d
        integer :: i

        ! Write X coordinates
        dims_1d(1) = size(x_cb)
        call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
        call h5dcreate_f(h5file_id, 'x_cb', H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
        call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, x_cb, dims_1d, h5err)
        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        ! Write Y coordinates if 2D or 3D
        if (n > 0) then
            dims_1d(1) = size(y_cb)
            call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, 'y_cb', H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, y_cb, dims_1d, h5err)
            call h5dclose_f(h5dset_id, h5err)
            call h5sclose_f(h5dspace_id, h5err)
        end if

        ! Write Z coordinates if 3D
        if (p > 0) then
            dims_1d(1) = size(z_cb)
            call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, 'z_cb', H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, z_cb, dims_1d, h5err)
            call h5dclose_f(h5dset_id, h5err)
            call h5sclose_f(h5dspace_id, h5err)
        end if

    end subroutine s_write_grid_to_hdf5

    !> @brief Write a flow variable to HDF5 file
    !! @param varname Name of the variable
    !! @param q_sf 3D array of variable data
    subroutine s_write_variable_to_hdf5(varname, q_sf)

        character(LEN=*), intent(in) :: varname
        real(wp), dimension(-offset_x%beg:, -offset_y%beg:, -offset_z%beg:), intent(in) :: q_sf

        integer(HSIZE_T), dimension(3) :: dims_3d
        integer(HSIZE_T), dimension(2) :: dims_2d
        integer(HSIZE_T), dimension(1) :: dims_1d

        ! Determine dimensions based on problem dimensionality
        if (p > 0) then
            ! 3D case
            dims_3d(1) = m + offset_x%beg + offset_x%end + 1
            dims_3d(2) = n + offset_y%beg + offset_y%end + 1
            dims_3d(3) = p + offset_z%beg + offset_z%end + 1
            call h5screate_simple_f(3, dims_3d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, q_sf, dims_3d, h5err)
        elseif (n > 0) then
            ! 2D case
            dims_2d(1) = m + offset_x%beg + offset_x%end + 1
            dims_2d(2) = n + offset_y%beg + offset_y%end + 1
            call h5screate_simple_f(2, dims_2d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, q_sf(:, :, 0), dims_2d, h5err)
        else
            ! 1D case
            dims_1d(1) = m + offset_x%beg + offset_x%end + 1
            call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), H5T_NATIVE_DOUBLE, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, H5T_NATIVE_DOUBLE, q_sf(:, 0, 0), dims_1d, h5err)
        end if

        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        ! Store variable name for XDMF generation
        num_stored_vars = num_stored_vars + 1
        stored_varnames(num_stored_vars) = trim(varname)

    end subroutine s_write_variable_to_hdf5

    !> @brief Close the current HDF5 file
    subroutine s_close_hdf5_file()

        call h5fclose_f(h5file_id, h5err)
        if (h5err /= 0) then
            call s_mpi_abort('Failed to close HDF5 file. Exiting.')
        end if

    end subroutine s_close_hdf5_file

    !> @brief Write XDMF metadata file for a time step (root process only)
    !! @param t_step Current time step
    !! @param time Physical time value
    subroutine s_write_xdmf_file(t_step, time)

        integer, intent(in) :: t_step
        real(wp), intent(in) :: time

        character(LEN=len_trim(case_dir) + 3*name_len) :: xdmf_file
        character(LEN=256) :: dims_str, h5_path
        integer :: unit_num, i, ip
        integer :: nx, ny, nz

        ! Only root process writes XDMF
        if (proc_rank /= 0) return

        ! Calculate grid dimensions (cell-centered values)
        nx = m + offset_x%beg + offset_x%end + 1
        if (n > 0) then
            ny = n + offset_y%beg + offset_y%end + 1
        else
            ny = 1
        end if
        if (p > 0) then
            nz = p + offset_z%beg + offset_z%end + 1
        else
            nz = 1
        end if

        ! Generate XDMF file path
        write (xdmf_file, '(A,I0,A)') '/timestep_', t_step, '.xdmf'
        xdmf_file = trim(rootdir_hdf5)//trim(xdmf_file)

        unit_num = 100
        open (unit=unit_num, file=trim(xdmf_file), status='replace', action='write')

        ! Write XDMF header
        write (unit_num, '(A)') '<?xml version="1.0" ?>'
        write (unit_num, '(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        write (unit_num, '(A)') '<Xdmf Version="3.0">'
        write (unit_num, '(A)') '  <Domain>'

        ! Write grid collection for all processors
        write (unit_num, '(A)') '    <Grid Name="MFC Output" GridType="Collection" CollectionType="Spatial">'
        write (unit_num, '(A,E15.8,A)') '      <Time Value="', time, '" />'

        ! Write grid for each processor
        do ip = 0, num_procs - 1
            write (h5_path, '(A,I0,A,I0,A)') '../p', ip, '/', t_step, '.h5'

            write (unit_num, '(A,I0,A)') '      <Grid Name="Processor_', ip, '" GridType="Uniform">'

            ! Topology (structured grid)
            if (p > 0) then
                write (dims_str, '(I0,A,I0,A,I0)') nz + 1, ' ', ny + 1, ' ', nx + 1
                write (unit_num, '(A,A,A)') '        <Topology TopologyType="3DRectMesh" Dimensions="', &
                    trim(dims_str), '"/>'
            elseif (n > 0) then
                write (dims_str, '(I0,A,I0)') ny + 1, ' ', nx + 1
                write (unit_num, '(A,A,A)') '        <Topology TopologyType="2DRectMesh" Dimensions="', &
                    trim(dims_str), '"/>'
            else
                write (dims_str, '(I0)') nx + 1
                write (unit_num, '(A,A,A)') '        <Topology TopologyType="2DRectMesh" Dimensions="1 ', &
                    trim(dims_str), '"/>'
            end if

            ! Geometry (rectilinear mesh coordinates)
            if (p > 0) then
                write (unit_num, '(A)') '        <Geometry GeometryType="VXVYVZ">'
                ! X coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', nx + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/x_cb'
                write (unit_num, '(A)') '          </DataItem>'
                ! Y coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', ny + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/y_cb'
                write (unit_num, '(A)') '          </DataItem>'
                ! Z coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', nz + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/z_cb'
                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Geometry>'
            elseif (n > 0) then
                write (unit_num, '(A)') '        <Geometry GeometryType="VXVY">'
                ! X coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', nx + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/x_cb'
                write (unit_num, '(A)') '          </DataItem>'
                ! Y coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', ny + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/y_cb'
                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Geometry>'
            else
                write (unit_num, '(A)') '        <Geometry GeometryType="VXVY">'
                ! X coordinates
                write (unit_num, '(A,I0,A)') '          <DataItem Dimensions="', nx + 1, &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A)') '            ', trim(h5_path), ':/x_cb'
                write (unit_num, '(A)') '          </DataItem>'
                ! Y = 0 for 1D (single value)
                write (unit_num, '(A)') '          <DataItem Dimensions="1" NumberType="Float" Precision="8" Format="XML">'
                write (unit_num, '(A)') '            0.0'
                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Geometry>'
            end if

            ! Attributes (flow variables)
            do i = 1, num_stored_vars
                write (unit_num, '(A,A,A)') '        <Attribute Name="', trim(stored_varnames(i)), &
                    '" AttributeType="Scalar" Center="Cell">'
                if (p > 0) then
                    write (dims_str, '(I0,A,I0,A,I0)') nz, ' ', ny, ' ', nx
                elseif (n > 0) then
                    write (dims_str, '(I0,A,I0)') ny, ' ', nx
                else
                    write (dims_str, '(I0)') nx
                end if
                write (unit_num, '(A,A,A)') '          <DataItem Dimensions="', trim(dims_str), &
                    '" NumberType="Float" Precision="8" Format="HDF">'
                write (unit_num, '(A,A,A,A)') '            ', trim(h5_path), ':/', trim(stored_varnames(i))
                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Attribute>'
            end do

            write (unit_num, '(A)') '      </Grid>'
        end do

        write (unit_num, '(A)') '    </Grid>'
        write (unit_num, '(A)') '  </Domain>'
        write (unit_num, '(A)') '</Xdmf>'

        close (unit_num)

    end subroutine s_write_xdmf_file

    !> @brief Finalize HDF5+XDMF module
    subroutine s_finalize_hdf5_xdmf_module()

        ! Close HDF5 library
        call h5close_f(h5err)

        ! Deallocate variable name storage
        if (allocated(stored_varnames)) deallocate (stored_varnames)

    end subroutine s_finalize_hdf5_xdmf_module

end module m_hdf5_xdmf_output
