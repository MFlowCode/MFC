#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
    integer, parameter :: nFiles = 14   ! Can be changed to any number
    integer, parameter :: nRows  = 512!
    integer :: f, iter, ios, unit, idx
    real(wp) :: x_len, x_step
    integer :: global_offset
    real(wp) :: delta
    character(len=100), dimension(nFiles) :: fileNames
    ! Arrays to store all data from files
    real(wp), dimension(nRows, nFiles) :: stored_values
    real(wp), dimension(nRows) :: x_coords
    logical :: files_loaded = .false.
    real(wp) :: domain_start, domain_end
    character(len=*), parameter :: init_dir = "/home/YourDirectory"
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "000000"  ! Default zeros (can be changed)
    
    ! Generate file names in a loop
    do f = 1, nFiles
        ! Convert file number to string with proper formatting
        if (f < 10) then
            write(file_num_str, '(I1)') f  ! Single digit
        else
            write(file_num_str, '(I2)') f  ! Double digit
        end if
            fileNames(f) = trim(init_dir) // "prim." // trim(file_num_str) // ".00." // zeros_default // ".dat"
        ! Create the filename with the pattern "prim.X.00.000000.dat"
    end do

#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100)
        ! Put your variable assignments here
        if (.not. files_loaded) then
            ! Print status message
            print *, "Loading all data files..."
                do f = 1, nFiles
                    ! Open the file for reading
                    open(newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                    if (ios /= 0) then
                        print *, "Error opening file: ", trim(fileNames(f))
                        cycle  ! Skip this file on error
                    endif
                ! Read all rows at once into memory
                    do iter = 1, nRows
                        read(unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                        if (ios /= 0) then
                            print *, "Error reading file ", trim(fileNames(f)), " at row ", iter
                            exit  ! Exit loop on error
                        endif
                    end do
                    close(unit)
                end do
                ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = x_cc(1) - x_cc(0)

            delta = x_cc(0)- domain_start-x_step/2.0
            global_offset = nint (abs(delta)/x_step)
            print *, "All data files loaded. Domain range:", domain_start, "to", domain_end
            files_loaded = .true.
        endif
        ! Simple mapping - find the closest index
        idx = i + 1 + global_offset
        do f = 1, nFiles
                q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx, f)
        end do

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
