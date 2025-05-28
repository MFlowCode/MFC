#:def Hardcoded1DVariables()

!> @file 1DHardcodedIC.fpp
!> @brief Initialize a simulation using data from a prior 1D run (for example from Cantera).
!>
!> @details
!> Reads a 1D primitive‐variable dataset and sets the
!> initial conditions for your simulation based on those 1D results.
!>
!> It provides:
!>   - Import of files named
!>     `prim.<variable>.00.<timestep>.dat`
!>     (produced when parallel I/O is disabled in `case.py`).
!>   - Configurable parameters:
!>       - @c nRows: total number of grid points in the 1D profile.
!>   - Default file directory: `examples/Case_File/D`
!>
!> @param nRows Total number of grid points in the imported 1D profile.

    ! Place any declaration of intermediate variables here
    integer, parameter :: nRows = 512   ! Number of grid points
    integer :: f, iter, ios, unit, idx
    real(wp) :: x_len, x_step
    integer :: global_offset            ! MPI subdomain offset
    real(wp) :: delta
    character(len=100), dimension(sys_size) :: fileNames ! Arrays to store all data from files
    character(len=200) :: errmsg
    real(wp), dimension(nRows, sys_size) :: stored_values  ! Imported Data
    real(wp), dimension(nRows) :: x_coords
    logical :: files_loaded = .false.
    real(wp) :: domain_start, domain_end
    character(len=*), parameter :: init_dir = "/home/YourDirectory" ! For example /home/MFC/examples/1D_Shock/D/
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "000000"  ! Default zeros (can be changed)

#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100)

        ! Generate file names in a loop
        do f = 1, sys_size
            ! Convert file number to string with proper formatting
            if (f < 10) then
                write (file_num_str, '(I1)') f  ! Single digit
            else
                write (file_num_str, '(I2)') f  ! Double digit
            end if
            fileNames(f) = trim(init_dir)//"prim."//trim(file_num_str)//".00."//zeros_default//".dat"
            ! Create the filename with the pattern "prim.X.00.000000.dat"
        end do
        ! Put your variable assignments here
        if (.not. files_loaded) then
            do f = 1, sys_size
                ! Open the file for reading
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0 .and. proc_rank == 0) then
                    write (errmsg, '(A,A)') "Error opening file: ", trim(fileNames(f))
                    call s_mpi_abort(trim(errmsg))
                end if
                ! Read all rows at once into memory
                do iter = 1, nRows
                    read (unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                    if (ios /= 0) then
                        write (errmsg, '(A,A,A,I0,A)') ' Error reading "', trim(fileNames(f)), &
                            '" at index (', iter, ')'
                        call s_mpi_abort(trim(errmsg)) ! Exit loop on error
                    end if
                end do
                close (unit)
            end do
            ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = x_cc(1) - x_cc(0)

            delta = x_cc(0) - domain_start - x_step/2.0
            global_offset = nint(abs(delta)/x_step)
            files_loaded = .true.
        end if
        ! Simple mapping - find the closest index
        idx = i + 1 + global_offset
        do f = 1, sys_size
            q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx, f)
        end do

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
