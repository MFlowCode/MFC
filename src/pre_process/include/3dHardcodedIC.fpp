#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps

#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.025_wp/wl

        intH = amp*(sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + sin(2._wp*pi*z_cc(k)/lam - pi/2._wp)) + h

        alph = 5e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0_wp
        lam = 1.0_wp
        amp = patch_icpp(patch_id)%a(2)
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

    case (370)
        ! Generate file names dynamically in a loop
        do f = 1, sys_size - 1
            ! Convert file number to string with proper formatting
            if (f < 10) then
                write (file_num_str, '(I1)') f  ! Single digit
            else
                write (file_num_str, '(I2)') f  ! Double digit
                ! For more than 99 files, you might need to adjust this format
            end if
            fileNames(f) = trim(init_dir)//"prim."//trim(file_num_str)//".00."//zeros_default//".dat"
        end do

        if (.not. files_loaded) then
          ! 1) Discover yRows by reading fileNames(1)  until x changes from the first line’s x.
          open(newunit=unit2, file=trim(fileNames(1)), status='old', action='read', iostat=ios2)
          if (ios2 /= 0 ) then
            write(errmsg, '(A,A)') "Error opening file for yRows detection: ", trim(fileNames(1))
            call s_mpi_abort(trim(errmsg))
          end if

          ! Read the very first line to get (x₀, y₀)
          read(unit2, *, iostat=ios2) x0, y0, dummy_z
          if (ios2 /= 0 ) then
            write(errmsg, '(A,A)') "Error reading first line of: ", trim(fileNames(1))
            call s_mpi_abort(trim(errmsg))
          end if

          ycount = 1
          do
            read(unit2, *, iostat=ios2) dummy_x, dummy_y, dummy_z
            if (ios2 /= 0) exit         ! End of File : stop counting
            if (dummy_x == x0 .and. dummy_y /= y0) then
              ! As soon as x or y changes from the first line’s (x0,y0),
              ! we know we have counted all the (x0, y1) rows.
              ycount = ycount + 1
            else
              exit
            end if
          end do
          yRows = ycount
          close(unit2)
          ! 2) Count total rows (nrows) to get xRows
          open(newunit=unit2, file=trim(fileNames(1)), status='old', action='read', iostat=ios2)
          if (ios2 /= 0 ) then
            write(errmsg, '(A,A)') "Error re‐opening file to count rows: ", trim(fileNames(1))
            call s_mpi_abort(trim(errmsg))
          end if

          nrows = 0
          do
            read(unit2, *, iostat=ios2) dummy_x, dummy_y, dummy_z
            if (ios2 /= 0) exit
            nrows = nrows + 1
          end do
          close(unit2)

          if (mod(nrows, yRows) /= 0) then
              write(errmsg, '(A,A,I0,A,I0)') &
                "File ’", trim(fileNames(1)), &
                "’ has total lines=", nrows, &
                " which is not a multiple of yRows=", yRows
              call s_mpi_abort(trim(errmsg))
          end if
            xRows = nrows / yRows
            allocate(x_coords(nRows))
            allocate(y_coords(nRows))
            allocate(stored_values(xRows, yRows, sys_size))
            index_x = i
            index_y = j
 
            open (newunit=unit, file=trim(fileNames(1)), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write (errmsg, '(A,A)') "Error opening file: ", trim(fileNames(f))
                call s_mpi_abort(trim(errmsg))
            end if

            iter = 0
            do iix = 1, xRows
                do iiy = 1, yRows
                    iter = iter + 1
                    read (unit, *, iostat=ios) x_coords(iter), y_coords(iter), stored_values(iix, iiy, 1)
                    if (ios /= 0) then
                        write (errmsg, '(A,A,A,I0,A)') 'Error reading "', trim(fileNames(1)), &
                            '" at indices (', iix, ',', iiy, ')'
                        call s_mpi_abort(trim(errmsg))
                    end if
                end do
            end do
            close (unit)

            !Now read only the values from remaining files (skip x,y columns)
            do f = 2, sys_size - 1
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle  ! Skip this file on error
                end if
                do iix = 1, xRows
                    do iiy = 1, yRows
                        ! Read and discard x,y, then read the value
                        read (unit, *, iostat=ios) dummy_x, dummy_y, stored_values(iix, iiy, f)
                        if (ios /= 0) then
                            write (errmsg, '(A,A,I0,A,I0,A)') 'Error reading "', trim(fileNames(f)), &
                                '" at indices (', iix, ',', iiy, ')'
                            call s_mpi_abort(trim(errmsg))
                        end if
                    end do
                end do
            end do
            close (unit)
           

            domain_xstart = x_coords(1)  ! First x value
            domain_ystart = y_coords(1)  ! First y value
            ! Calculate actual domain end values
            domain_xend = x_coords(nRows - yRows + 1)  ! Last x value
            domain_yend = y_coords(yRows)              ! Last y value in first x-slice
            ! Calculate simulation grid steps
            x_step = x_cc(1) - x_cc(0)
            y_step = y_cc(1) - y_cc(0)
            ! Calculate offsets
            delta_x = x_cc(index_x) - domain_xstart + x_step/2.0_wp
            delta_y = y_cc(index_y) - domain_ystart + y_step/2.0_wp
            ! Global offsets in case of multiple ranks
            global_offset_x = nint(abs(delta_x)/x_step)
            global_offset_y = nint(abs(delta_y)/y_step)

            files_loaded = .true.
        end if

        idx = i + 1 + global_offset_x - index_x  ! x-direction index in data grid
        idy = j + 1 + global_offset_y - index_y  ! y-direction index in data grid

        do f = 1, sys_size - 1
            if (f >= momxe) then
                jump = 1
            else
                jump = 0
            end if
            q_prim_vf(f + jump)%sf(i, j, k) = stored_values(idx, idy, f)
        end do
        ! Set z velocity to zero
        q_prim_vf(momxe)%sf(i, j, k) = 0.0_wp

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
