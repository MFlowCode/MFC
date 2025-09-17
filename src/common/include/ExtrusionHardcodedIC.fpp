!> @brief Allocate memory and read initial condition data for IC extrusion.
!>
!> @details
!>   This macro handles the complete initialization process for IC extrusion by:
!>
!>   **Memory Allocation:**
!>     - stored_values(xRows, yRows, sys_size) - stores primitive variable data from files
!>     - x_coords(nrows) - stores x-coordinates from input files
!>     - y_coords(nrows) - stores y-coordinates from input files (3D case only)
!>
!>   **File Reading Operations:**
!>     - Reads primitive variable data from multiple files with pattern:
!>       `prim.<file_number>.00.<timestep>.dat` where timestep uses `zeros_default` padding
!>     - Files are read from directory specified by `init_dir` parameter
!>     - Supports 1D, 2D, and 3D computational domains
!>
!>   **Grid Structure Detection:**
!>     - 1D/2D: Counts lines in first file to determine xRows
!>     - 3D: Analyzes coordinate patterns to determine xRows and yRows structure
!>
!>   **MPI Domain Mapping:**
!>     - Calculates global_offset_x and global_offset_y for MPI subdomain positioning
!>     - Maps file coordinates to local computational grid coordinates
!>
!>   **Data Assignment:**
!>     - Populates q_prim_vf primitive variable arrays with file data
!>     - Handles momentum component indexing with special treatment for momxe
!>     - Sets momxe component to zero for 2D/3D cases
!>
!>   **State Management:**
!>     - Uses files_loaded flag to prevent redundant file operations
!>     - Preserves data across multiple macro calls within same simulation
!>
!>   @note File pattern uses `zeros_default` parameter (default: "000000") for timestep padding
!>   @note Directory path is hardcoded in `init_dir` parameter - modify as needed
!>   @warning Aborts execution if file reading errors occur.

#:def HardcodedDimensionsExtrusion()
    integer :: xRows, yRows, nRows, iix, iiy, max_files
    integer :: f, iter, ios, ios2, unit, unit2, idx, idy, index_x, index_y, jump, line_count, ycount
    real(wp) :: x_len, x_step, y_len, y_step
    real(wp) :: dummy_x, dummy_y, dummy_z, x0, y0
    integer :: global_offset_x, global_offset_y           ! MPI subdomain offset
    real(wp) :: delta_x, delta_y
    character(len=100), dimension(sys_size) :: fileNames ! Arrays to store all data from files
    character(len=200) :: errmsg
    real(wp), allocatable :: stored_values(:, :, :)
    real(wp), allocatable :: x_coords(:), y_coords(:)
    logical :: files_loaded = .false.
    real(wp) :: domain_xstart, domain_xend, domain_ystart, domain_yend
    character(len=*), parameter :: init_dir = "/home/MFC/FilesDirectory" ! For example /home/MFC/examples/1D_Shock/D/
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "000000"  ! Default zeros (can be changed)
#:enddef

#:def HardcodedReadValues()

    if (.not. files_loaded) then
        max_files = merge(sys_size, sys_size - 1, num_dims == 1)
        do f = 1, max_files
            write (file_num_str, '(I0)') f
            fileNames(f) = trim(init_dir)//"prim."//trim(file_num_str)//".00."//zeros_default//".dat"
        end do

        ! Common file reading setup
        open (newunit=unit2, file=trim(fileNames(1)), status='old', action='read', iostat=ios2)
        if (ios2 /= 0) call s_mpi_abort("Error opening file: "//trim(fileNames(1)))

        select case (num_dims)
        case (1, 2)  ! 1D and 2D cases are similar
            ! Count lines
            line_count = 0
            do
                read (unit2, *, iostat=ios2) dummy_x, dummy_y
                if (ios2 /= 0) exit
                line_count = line_count + 1
            end do
            close (unit2)

            xRows = line_count
            yRows = 1
            index_x = 0
            if (num_dims == 2) index_x = i
            @:ALLOCATE (x_coords(xRows), stored_values(xRows, 1, sys_size))

            ! Read data from all files
            do f = 1, max_files
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) call s_mpi_abort("Error opening file: "//trim(fileNames(f)))

                do iter = 1, xRows
                    read (unit, *, iostat=ios) x_coords(iter), stored_values(iter, 1, f)
                    if (ios /= 0) call s_mpi_abort("Error reading file: "//trim(fileNames(f)))
                end do
                close (unit)
            end do

            ! Calculate offsets
            domain_xstart = x_coords(1)
            x_step = x_cc(1) - x_cc(0)
            delta_x = merge(x_cc(0) - domain_xstart + x_step/2.0, &
                            x_cc(index_x) - domain_xstart + x_step/2.0, num_dims == 1)
            global_offset_x = nint(abs(delta_x)/x_step)

        case (3)  ! 3D case - determine grid structure
            ! Find yRows by counting rows with same x
            read (unit2, *, iostat=ios2) x0, y0, dummy_z
            if (ios2 /= 0) call s_mpi_abort("Error reading first line")

            yRows = 1
            do
                read (unit2, *, iostat=ios2) dummy_x, dummy_y, dummy_z
                if (ios2 /= 0) exit
                if (dummy_x == x0 .and. dummy_y /= y0) then
                    yRows = yRows + 1
                else
                    exit
                end if
            end do
            close (unit2)

            ! Count total rows
            open (newunit=unit2, file=trim(fileNames(1)), status='old', action='read', iostat=ios2)
            nrows = 0
            do
                read (unit2, *, iostat=ios2) dummy_x, dummy_y, dummy_z
                if (ios2 /= 0) exit
                nrows = nrows + 1
            end do
            close (unit2)

            xRows = nrows/yRows
            @:ALLOCATE (x_coords(nrows), y_coords(nrows), stored_values(xRows, yRows, sys_size))
            index_x = i
            index_y = j

            ! Read all files
            do f = 1, max_files
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    if (f == 1) call s_mpi_abort("Error opening file: "//trim(fileNames(f)))
                    cycle
                end if

                iter = 0
                do iix = 1, xRows
                    do iiy = 1, yRows
                        iter = iter + 1
                        if (f == 1) then
                            read (unit, *, iostat=ios) x_coords(iter), y_coords(iter), stored_values(iix, iiy, f)
                        else
                            read (unit, *, iostat=ios) dummy_x, dummy_y, stored_values(iix, iiy, f)
                        end if
                        if (ios /= 0) call s_mpi_abort("Error reading data")
                    end do
                end do
                close (unit)
            end do

            ! Calculate offsets
            x_step = x_cc(1) - x_cc(0)
            y_step = y_cc(1) - y_cc(0)
            delta_x = x_cc(index_x) - x_coords(1) + x_step/2.0_wp
            delta_y = y_cc(index_y) - y_coords(1) + y_step/2.0_wp
            global_offset_x = nint(abs(delta_x)/x_step)
            global_offset_y = nint(abs(delta_y)/y_step)
        end select

        files_loaded = .true.
    end if

    ! Data assignment
    select case (num_dims)
    case (1)
        idx = i + 1 + global_offset_x
        do f = 1, sys_size
            q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx, 1, f)
        end do

    case (2)
        idx = i + 1 + global_offset_x - index_x
        do f = 1, sys_size - 1
            jump = merge(1, 0, f >= momxe)
            q_prim_vf(f + jump)%sf(i, j, 0) = stored_values(idx, 1, f)
        end do
        q_prim_vf(momxe)%sf(i, j, 0) = 0.0_wp

    case (3)
        idx = i + 1 + global_offset_x - index_x
        idy = j + 1 + global_offset_y - index_y
        do f = 1, sys_size - 1
            jump = merge(1, 0, f >= momxe)
            q_prim_vf(f + jump)%sf(i, j, k) = stored_values(idx, idy, f)
        end do
        q_prim_vf(momxe)%sf(i, j, k) = 0.0_wp
    end select
#:enddef

#:def HardcodedDellacation()
    if (allocated(stored_values)) then
        @:DEALLOCATE (stored_values)
        @:DEALLOCATE (x_coords)
    end if

    if (allocated(y_coords)) then
        @:DEALLOCATE (y_coords)
    end if
#:enddef
