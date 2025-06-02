#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (170)

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
        
        if (.not. files_loaded) then
            !Reading the first file to calculate the number of grid points in the x direction
            line_count = 0
            open(newunit=unit2, file=trim(fileNames(1)), status='old', action='read', iostat=ios2)
                if (ios2 /= 0) then
                    write(errmsg, '(A,A)') "Error opening file: ", trim(fileNames(1))
                    call s_mpi_abort(trim(errmsg))
                end if
                do
                    read(unit2, *, iostat=ios2) dummy_x, dummy_y
                    if (ios2 /= 0) exit  ! Exit since the file has been read
                    line_count = line_count + 1
                end do
            close(unit2)
            
            xRows = line_count
            allocate(x_coords(xRows))
            allocate(stored_values(xRows, 1, sys_size))
            do f = 1, sys_size
                ! Open the file for reading
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    write (errmsg, '(A,A)') "Error opening file: ", trim(fileNames(f))
                    call s_mpi_abort(trim(errmsg))
                end if
                ! Read all rows at once into memory
                do iter = 1, xRows
                    read (unit, *, iostat=ios) x_coords(iter), stored_values(iter, 1, f)
                    if (ios /= 0) then
                        write (errmsg, '(A,A,A,I0,A)') ' Error reading "', trim(fileNames(f)), &
                            '" at index (', iter, ')'
                        call s_mpi_abort(trim(errmsg)) ! Exit loop on error
                    end if
                end do
                close (unit)
            end do
            ! Calculate domain information for mapping
            domain_xstart = x_coords(1)
            domain_xend = x_coords(xRows)
            x_step = x_cc(1) - x_cc(0)

            delta_x = x_cc(0) - domain_xstart - x_step/2.0
            global_offset_x = nint(abs(delta_x)/x_step)
            files_loaded = .true.
        end if
        ! Simple mapping - find the closest index
        idx = i + 1 + global_offset_x
        do f = 1, sys_size
            q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx, 1, f)
        end do

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select
     

#:enddef
