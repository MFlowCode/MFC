!>
!! @file m_hdf5_xdmf_output.fpp
!! @brief Contains module m_hdf5_xdmf_output
!!
!! @details This module implements native HDF5 output with XDMF metadata for
!!          visualization in ParaView 5.12+ and other XDMF-compatible tools.
!!          It provides a portable alternative to Silo-HDF5 format that works
!!          reliably across different visualization software versions.
!!
!! @par Output Structure
!! The module creates the following directory structure:
!! @verbatim
!!   case_dir/hdf5_xdmf/
!!   ├── p0/                # Per-processor HDF5 files
!!   │   ├── 0.h5           # Time step 0
!!   │   └── 100.h5         # Time step 100
!!   ├── p1/                # Processor 1 files (for parallel runs)
!!   │   └── ...
!!   └── root/              # XDMF metadata files (written by rank 0)
!!       ├── timestep_0.xdmf
!!       └── timestep_100.xdmf
!! @endverbatim
!!
!! @par Supported Data Types
!! - Structured grid (rectilinear mesh) data for 1D/2D/3D simulations
!! - Particle/point mesh data for Lagrangian bubble tracking
!! - All flow variables (density, velocity, pressure, etc.)
!! - MHD magnetic field components (Bx, By, Bz)
!! - Chemistry species mass fractions
!!
!! @par Usage
!! Set format=3 in the post_process input file to enable HDF5+XDMF output.
!! Open the .xdmf files from the root/ folder in ParaView to visualize.
!!

!> @brief This module provides native HDF5 + XDMF output support for MFC.
!!
!! XDMF (eXtensible Data Model and Format) is an XML-based metadata format that
!! describes HDF5 data for visualization tools like ParaView.
!!
!! @note This module requires HDF5 with Fortran bindings (hdf5-mpi recommended).
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
 s_write_particle_mesh_to_hdf5, &
 s_write_particle_variable_to_hdf5, &
 s_close_hdf5_file, &
 s_write_xdmf_file, &
 s_write_xdmf_time_series, &
 s_finalize_hdf5_xdmf_module

    !> @name HDF5 Identifiers
    !! @{
    integer(HID_T) :: h5file_id      !< HDF5 file identifier
    integer(HID_T) :: h5dspace_id    !< HDF5 dataspace identifier
    integer(HID_T) :: h5dset_id      !< HDF5 dataset identifier
    integer(HID_T) :: h5plist_id     !< HDF5 property list identifier
    integer(HID_T) :: h5type_id      !< HDF5 datatype (single or double precision)
    !> @}

    !> @name Directory Paths
    !! @{
    character(LEN=path_len + name_len) :: hdf5_dir           !< Base HDF5 output directory
    character(LEN=path_len + 2*name_len) :: proc_rank_dir_hdf5  !< Per-processor directory
    character(LEN=path_len + 2*name_len) :: rootdir_hdf5     !< Root directory for XDMF files
    !> @}

    integer :: h5err  !< HDF5 error flag

    !> @name Variable Tracking for XDMF
    !! These arrays store variable names and particle information for XDMF generation.
    !! @{
    character(LEN=name_len), allocatable, dimension(:) :: stored_varnames  !< Grid variable names
    integer :: num_stored_vars  !< Number of stored grid variables

    character(LEN=name_len), allocatable, dimension(:) :: stored_particle_varnames  !< Particle variable names
    integer :: num_stored_particle_vars  !< Number of stored particle variables
    integer :: stored_num_particles      !< Number of particles for current timestep
    logical :: has_particle_data         !< Flag indicating particle data exists
    !> @}

contains

    !> @brief Initialize the HDF5+XDMF output module
    !!
    !! This subroutine initializes the HDF5 library and creates the directory
    !! structure for output files. It must be called before any other routines
    !! in this module.
    !!
    !! @par Directory Structure Created
    !! - hdf5_xdmf/p{rank}/ for each MPI rank
    !! - hdf5_xdmf/root/ for XDMF metadata (rank 0 only)
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
        allocate (stored_varnames(200))  ! Max 200 grid variables
        num_stored_vars = 0

        ! Allocate particle variable name storage
        allocate (stored_particle_varnames(50))  ! Max 50 particle variables
        num_stored_particle_vars = 0
        stored_num_particles = 0
        has_particle_data = .false.

        ! Set HDF5 datatype based on precision setting
        ! precision=1 -> single (32-bit), precision=2 -> double (64-bit)
        if (precision == 1) then
            h5type_id = H5T_NATIVE_REAL
        else
            h5type_id = H5T_NATIVE_DOUBLE
        end if

    end subroutine s_initialize_hdf5_xdmf_module

    !> @brief Open an HDF5 file for the current time step
    !!
    !! Creates a new HDF5 file for the specified time step in the processor's
    !! output directory. Resets variable counters for XDMF metadata generation.
    !!
    !! @param[in] t_step Current time step number
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

        ! Reset variable counters for this timestep
        num_stored_vars = 0
        num_stored_particle_vars = 0
        stored_num_particles = 0
        has_particle_data = .false.

    end subroutine s_open_hdf5_file

    !> @brief Write grid coordinates to HDF5 file
    !!
    !! Writes the cell boundary coordinates (x_cb, y_cb, z_cb) to the HDF5 file.
    !! These define the rectilinear mesh geometry for visualization.
    !! Grid coordinates are always written in double precision for accuracy.
    !! Also writes spatial extent metadata as HDF5 attributes for optimization.
    subroutine s_write_grid_to_hdf5()

        integer(HSIZE_T), dimension(1) :: dims_1d
        real(wp) :: spatial_extent(6)  ! xmin, xmax, ymin, ymax, zmin, zmax
        integer :: i

        ! Write X coordinates (always double precision for grid accuracy)
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

        ! Write spatial extent as file-level attribute for visualization optimization
        ! Format: [xmin, xmax, ymin, ymax, zmin, zmax]
        spatial_extent(1) = x_cb(-1)                              ! xmin
        spatial_extent(2) = x_cb(m)                               ! xmax
        if (n > 0) then
            spatial_extent(3) = y_cb(-1)                          ! ymin
            spatial_extent(4) = y_cb(n)                           ! ymax
        else
            spatial_extent(3) = 0._wp
            spatial_extent(4) = 0._wp
        end if
        if (p > 0) then
            spatial_extent(5) = z_cb(-1)                          ! zmin
            spatial_extent(6) = z_cb(p)                           ! zmax
        else
            spatial_extent(5) = 0._wp
            spatial_extent(6) = 0._wp
        end if
        call s_write_attribute_real_array(h5file_id, 'spatial_extent', spatial_extent, 6)

    end subroutine s_write_grid_to_hdf5

    !> @brief Write a flow variable to HDF5 file
    !!
    !! Writes cell-centered flow variable data to the HDF5 file. The data
    !! is stored in a dataset named after the variable. Respects the precision
    !! setting (precision=1 for single, precision=2 for double).
    !!
    !! @param[in] varname Name of the variable (used as HDF5 dataset name)
    !! @param[in] q_sf    3D array containing the variable data
    subroutine s_write_variable_to_hdf5(varname, q_sf)

        character(LEN=*), intent(in) :: varname
        real(wp), dimension(-offset_x%beg:, -offset_y%beg:, -offset_z%beg:), intent(in) :: q_sf

        integer(HSIZE_T), dimension(3) :: dims_3d
        integer(HSIZE_T), dimension(2) :: dims_2d
        integer(HSIZE_T), dimension(1) :: dims_1d

        real(wp) :: data_range(2)  ! min, max values for this variable

        ! Determine dimensions based on problem dimensionality
        if (p > 0) then
            ! 3D case
            dims_3d(1) = m + offset_x%beg + offset_x%end + 1
            dims_3d(2) = n + offset_y%beg + offset_y%end + 1
            dims_3d(3) = p + offset_z%beg + offset_z%end + 1
            call h5screate_simple_f(3, dims_3d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), h5type_id, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, h5type_id, q_sf, dims_3d, h5err)
            ! Compute data range for this block
            data_range(1) = minval(q_sf)
            data_range(2) = maxval(q_sf)
        elseif (n > 0) then
            ! 2D case
            dims_2d(1) = m + offset_x%beg + offset_x%end + 1
            dims_2d(2) = n + offset_y%beg + offset_y%end + 1
            call h5screate_simple_f(2, dims_2d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), h5type_id, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, h5type_id, q_sf(:, :, 0), dims_2d, h5err)
            ! Compute data range for this block
            data_range(1) = minval(q_sf(:, :, 0))
            data_range(2) = maxval(q_sf(:, :, 0))
        else
            ! 1D case
            dims_1d(1) = m + offset_x%beg + offset_x%end + 1
            call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
            call h5dcreate_f(h5file_id, trim(varname), h5type_id, h5dspace_id, h5dset_id, h5err)
            call h5dwrite_f(h5dset_id, h5type_id, q_sf(:, 0, 0), dims_1d, h5err)
            ! Compute data range for this block
            data_range(1) = minval(q_sf(:, 0, 0))
            data_range(2) = maxval(q_sf(:, 0, 0))
        end if

        ! Write data range as dataset attribute for visualization optimization
        call s_write_attribute_real_array(h5dset_id, 'data_range', data_range, 2)

        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        ! Store variable name for XDMF generation
        num_stored_vars = num_stored_vars + 1
        stored_varnames(num_stored_vars) = trim(varname)

    end subroutine s_write_variable_to_hdf5

    !> @brief Write particle/point mesh coordinates to HDF5 file
    !!
    !! Creates datasets for Lagrangian particle positions (px, py, pz).
    !! This data is used by XDMF to define a Polyvertex topology for
    !! particle visualization. Uses h5type_id for precision-aware output.
    !!
    !! @param[in] px     Array of particle x-coordinates
    !! @param[in] py     Array of particle y-coordinates
    !! @param[in] pz     Array of particle z-coordinates
    !! @param[in] nBub   Number of particles/bubbles on this processor
    subroutine s_write_particle_mesh_to_hdf5(px, py, pz, nBub)

        real(wp), dimension(:), intent(in) :: px, py, pz
        integer, intent(in) :: nBub

        integer(HSIZE_T), dimension(1) :: dims_1d
        integer(HID_T) :: grp_id

        ! Skip if no particles
        if (nBub <= 0) then
            stored_num_particles = 0
            has_particle_data = .false.
            return
        end if

        ! Create a group for particle data
        call h5gcreate_f(h5file_id, 'particles', grp_id, h5err)

        dims_1d(1) = nBub

        ! Write X coordinates
        call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
        call h5dcreate_f(grp_id, 'px', h5type_id, h5dspace_id, h5dset_id, h5err)
        call h5dwrite_f(h5dset_id, h5type_id, px, dims_1d, h5err)
        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        ! Write Y coordinates
        call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
        call h5dcreate_f(grp_id, 'py', h5type_id, h5dspace_id, h5dset_id, h5err)
        call h5dwrite_f(h5dset_id, h5type_id, py, dims_1d, h5err)
        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        ! Write Z coordinates
        call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
        call h5dcreate_f(grp_id, 'pz', h5type_id, h5dspace_id, h5dset_id, h5err)
        call h5dwrite_f(h5dset_id, h5type_id, pz, dims_1d, h5err)
        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        call h5gclose_f(grp_id, h5err)

        ! Store particle count for XDMF generation
        stored_num_particles = nBub
        has_particle_data = .true.
        num_stored_particle_vars = 0

    end subroutine s_write_particle_mesh_to_hdf5

    !> @brief Write a particle variable to HDF5 file
    !!
    !! Writes scalar data associated with particles (e.g., radius, velocity
    !! components, pressure) to the HDF5 file in the particles group.
    !! Uses h5type_id for precision-aware output.
    !!
    !! @param[in] varname Name of the variable (used as HDF5 dataset name)
    !! @param[in] data    Array of variable values (one per particle)
    !! @param[in] nBub    Number of particles/bubbles
    subroutine s_write_particle_variable_to_hdf5(varname, data, nBub)

        character(LEN=*), intent(in) :: varname
        real(wp), dimension(:), intent(in) :: data
        integer, intent(in) :: nBub

        integer(HSIZE_T), dimension(1) :: dims_1d
        integer(HID_T) :: grp_id
        character(LEN=name_len) :: dset_path

        ! Skip if no particles
        if (nBub <= 0) return

        ! Open particles group
        call h5gopen_f(h5file_id, 'particles', grp_id, h5err)

        dims_1d(1) = nBub

        ! Write variable data
        call h5screate_simple_f(1, dims_1d, h5dspace_id, h5err)
        call h5dcreate_f(grp_id, trim(varname), h5type_id, h5dspace_id, h5dset_id, h5err)
        call h5dwrite_f(h5dset_id, h5type_id, data, dims_1d, h5err)
        call h5dclose_f(h5dset_id, h5err)
        call h5sclose_f(h5dspace_id, h5err)

        call h5gclose_f(grp_id, h5err)

        ! Store variable name for XDMF generation
        num_stored_particle_vars = num_stored_particle_vars + 1
        stored_particle_varnames(num_stored_particle_vars) = trim(varname)

    end subroutine s_write_particle_variable_to_hdf5

    !> @brief Close the current HDF5 file
    subroutine s_close_hdf5_file()

        call h5fclose_f(h5file_id, h5err)
        if (h5err /= 0) then
            call s_mpi_abort('Failed to close HDF5 file. Exiting.')
        end if

    end subroutine s_close_hdf5_file

    !> @brief Write XDMF metadata file for a time step (root process only)
    !!
    !! Creates an XDMF file that describes the structure and location of data
    !! in the HDF5 files. This file can be opened directly in ParaView.
    !! The precision attribute in XDMF is set based on the global precision
    !! setting: 4 bytes for single precision, 8 bytes for double.
    !!
    !! @param[in] t_step Current time step number
    !! @param[in] time   Physical simulation time
    subroutine s_write_xdmf_file(t_step, time)

        integer, intent(in) :: t_step
        real(wp), intent(in) :: time

        character(LEN=len_trim(case_dir) + 3*name_len) :: xdmf_file
        character(LEN=256) :: dims_str, h5_path
        character(LEN=1) :: prec_str   !< Precision string for XDMF ("4" or "8")
        integer :: unit_num, i, ip
        integer :: nx, ny, nz

        ! Only root process writes XDMF
        if (proc_rank /= 0) return

        ! Set precision string based on global precision setting
        ! Grid coordinates always use double precision ("8") for accuracy
        ! Variable data uses the user-specified precision
        if (precision == 1) then
            prec_str = '4'  ! Single precision (32-bit float)
        else
            prec_str = '8'  ! Double precision (64-bit float)
        end if

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

            ! Attributes (flow variables) - use precision-aware output
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
                write (unit_num, '(A,A,A,A,A)') '          <DataItem Dimensions="', trim(dims_str), &
                    '" NumberType="Float" Precision="', prec_str, '" Format="HDF">'
                write (unit_num, '(A,A,A,A)') '            ', trim(h5_path), ':/', trim(stored_varnames(i))
                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Attribute>'
            end do

            write (unit_num, '(A)') '      </Grid>'
        end do

        ! Write particle data grids if particle data exists
        if (has_particle_data .and. stored_num_particles > 0) then
            do ip = 0, num_procs - 1
                write (h5_path, '(A,I0,A,I0,A)') '../p', ip, '/', t_step, '.h5'

                write (unit_num, '(A,I0,A)') '      <Grid Name="Particles_', ip, '" GridType="Uniform">'

                ! Topology (polyvertex for points)
                write (unit_num, '(A,I0,A)') '        <Topology TopologyType="Polyvertex" NumberOfElements="', &
                    stored_num_particles, '"/>'

                ! Geometry (XYZ coordinates) - use precision-aware output
                write (unit_num, '(A)') '        <Geometry GeometryType="XYZ">'
                write (unit_num, '(A,I0,A,A,A)') '          <DataItem Dimensions="', stored_num_particles, &
                    ' 3" NumberType="Float" Precision="', prec_str, '" Format="HDF">'

                ! Join px, py, pz into XYZ using DataItem function
                write (unit_num, '(A)') '            <DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="' &
                    //trim(adjustl(int_to_str(stored_num_particles)))//' 3">'
                write (unit_num, '(A,I0,A,A,A)') '              <DataItem Dimensions="', stored_num_particles, &
                    '" NumberType="Float" Precision="', prec_str, '" Format="HDF">'
                write (unit_num, '(A,A,A)') '                ', trim(h5_path), ':/particles/px'
                write (unit_num, '(A)') '              </DataItem>'
                write (unit_num, '(A,I0,A,A,A)') '              <DataItem Dimensions="', stored_num_particles, &
                    '" NumberType="Float" Precision="', prec_str, '" Format="HDF">'
                write (unit_num, '(A,A,A)') '                ', trim(h5_path), ':/particles/py'
                write (unit_num, '(A)') '              </DataItem>'
                write (unit_num, '(A,I0,A,A,A)') '              <DataItem Dimensions="', stored_num_particles, &
                    '" NumberType="Float" Precision="', prec_str, '" Format="HDF">'
                write (unit_num, '(A,A,A)') '                ', trim(h5_path), ':/particles/pz'
                write (unit_num, '(A)') '              </DataItem>'
                write (unit_num, '(A)') '            </DataItem>'

                write (unit_num, '(A)') '          </DataItem>'
                write (unit_num, '(A)') '        </Geometry>'

                ! Particle attributes - use precision-aware output
                do i = 1, num_stored_particle_vars
                    write (unit_num, '(A,A,A)') '        <Attribute Name="', trim(stored_particle_varnames(i)), &
                        '" AttributeType="Scalar" Center="Node">'
                    write (unit_num, '(A,I0,A,A,A)') '          <DataItem Dimensions="', stored_num_particles, &
                        '" NumberType="Float" Precision="', prec_str, '" Format="HDF">'
                    write (unit_num, '(A,A,A,A)') '            ', trim(h5_path), ':/particles/', &
                        trim(stored_particle_varnames(i))
                    write (unit_num, '(A)') '          </DataItem>'
                    write (unit_num, '(A)') '        </Attribute>'
                end do

                write (unit_num, '(A)') '      </Grid>'
            end do
        end if

        write (unit_num, '(A)') '    </Grid>'
        write (unit_num, '(A)') '  </Domain>'
        write (unit_num, '(A)') '</Xdmf>'

        close (unit_num)

    end subroutine s_write_xdmf_file

    !> @brief Write master XDMF time series file (root process only)
    !!
    !! Creates a master XDMF file that references all individual timestep
    !! XDMF files using XI:Include. This allows ParaView to load all timesteps
    !! as a single time-varying dataset using the temporal collection.
    !!
    !! @param[in] t_step_start First time step to include
    !! @param[in] t_step_stop  Last time step to include
    !! @param[in] t_step_save  Time step interval between saves
    subroutine s_write_xdmf_time_series(t_step_start, t_step_stop, t_step_save)

        integer, intent(in) :: t_step_start
        integer, intent(in) :: t_step_stop
        integer, intent(in) :: t_step_save

        character(LEN=len_trim(case_dir) + 3*name_len) :: xdmf_file
        integer :: unit_num, t_step

        ! Only root process writes XDMF
        if (proc_rank /= 0) return

        ! Generate time series XDMF file path
        xdmf_file = trim(rootdir_hdf5)//'/time_series.xdmf'

        unit_num = 101
        open (unit=unit_num, file=trim(xdmf_file), status='replace', action='write')

        ! Write XDMF header with XInclude namespace
        write (unit_num, '(A)') '<?xml version="1.0" ?>'
        write (unit_num, '(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        write (unit_num, '(A)') '<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write (unit_num, '(A)') '  <Domain>'
        write (unit_num, '(A)') '    <Grid Name="MFC Time Series" GridType="Collection" CollectionType="Temporal">'

        ! Include each timestep XDMF file
        do t_step = t_step_start, t_step_stop, t_step_save
            write (unit_num, '(A,I0,A)') '      <xi:include href="timestep_', t_step, &
                '.xdmf" xpointer="xpointer(//Xdmf/Domain/Grid)"/>'
        end do

        write (unit_num, '(A)') '    </Grid>'
        write (unit_num, '(A)') '  </Domain>'
        write (unit_num, '(A)') '</Xdmf>'

        close (unit_num)

    end subroutine s_write_xdmf_time_series

    !> @brief Convert integer to string (helper function)
    !! @param[in] i Integer to convert
    !! @return String representation
    pure function int_to_str(i) result(str)
        integer, intent(in) :: i
        character(LEN=20) :: str

        write (str, '(I0)') i
        str = adjustl(str)
    end function int_to_str

    !> @brief Write a real array as an HDF5 attribute
    !!
    !! Attaches a 1D array of real values as an attribute to an HDF5 object
    !! (file or dataset). Used for storing metadata like spatial extents
    !! and data ranges for visualization optimization.
    !!
    !! @param[in] obj_id   HDF5 object ID (file or dataset) to attach attribute to
    !! @param[in] attr_name Name of the attribute
    !! @param[in] values   Array of real values to write
    !! @param[in] n        Number of values in the array
    subroutine s_write_attribute_real_array(obj_id, attr_name, values, n)

        integer(HID_T), intent(in) :: obj_id
        character(LEN=*), intent(in) :: attr_name
        integer, intent(in) :: n
        real(wp), dimension(n), intent(in) :: values

        integer(HID_T) :: attr_space_id, attr_id
        integer(HSIZE_T), dimension(1) :: attr_dims
        integer :: attr_err

        attr_dims(1) = n

        ! Create dataspace for attribute
        call h5screate_simple_f(1, attr_dims, attr_space_id, attr_err)

        ! Create and write attribute
        call h5acreate_f(obj_id, trim(attr_name), H5T_NATIVE_DOUBLE, attr_space_id, attr_id, attr_err)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, values, attr_dims, attr_err)

        ! Close attribute and dataspace
        call h5aclose_f(attr_id, attr_err)
        call h5sclose_f(attr_space_id, attr_err)

    end subroutine s_write_attribute_real_array

    !> @brief Finalize HDF5+XDMF module
    !!
    !! Closes the HDF5 library and deallocates module storage. Should be
    !! called at the end of post-processing.
    subroutine s_finalize_hdf5_xdmf_module()

        ! Close HDF5 library
        call h5close_f(h5err)

        ! Deallocate variable name storage
        if (allocated(stored_varnames)) deallocate (stored_varnames)
        if (allocated(stored_particle_varnames)) deallocate (stored_particle_varnames)

    end subroutine s_finalize_hdf5_xdmf_module

end module m_hdf5_xdmf_output
