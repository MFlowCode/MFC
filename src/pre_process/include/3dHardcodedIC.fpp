#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps

    integer, parameter :: nFiles = 15   ! Can be changed to any number
    integer, parameter :: xRows  = 401
    integer, parameter :: yRows  = 403
    integer, parameter :: nRows  = xRows*yRows
    integer :: f, iter, ios, unit, idx, idy, jump, index_x, index_y
    real(wp) :: x_len, x_step, y_len, y_step
    integer :: global_offset_x, global_offset_y
    real(wp) :: delta_x, delta_y
    real(wp) :: dummy_x, dummy_y
    real(wp) :: current_x
    real(wp) :: sim_x_step, sim_y_step
    integer :: nx_data, ny_data, data_i, data_j, iix, iiy
    character(len=100), dimension(nFiles) :: fileNames
    ! Arrays to store all data from files - read once, use many times
    real(wp), dimension(xRows, yRows, nFiles) :: stored_values
    real(wp), dimension(nRows) :: x_coords, y_coords
    logical :: files_loaded = .false.
    real(wp) :: domain_xstart, domain_xend, domain_ystart, domain_yend
    character(len=*), parameter :: init_dir = "/home/YourDirectory"
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "000000"  ! Default zeros (can be changed)
    
    ! Generate file names dynamically in a loop
    do f = 1, nFiles
        ! Convert file number to string with proper formatting
        if (f < 10) then
            write(file_num_str, '(I1)') f  ! Single digit
        else
            write(file_num_str, '(I2)') f  ! Double digit
            ! For more than 99 files, you might need to adjust this format
        end if
        fileNames(f) = trim(init_dir) // "prim." // trim(file_num_str) // ".00." // zeros_default // ".dat"
    end do
    eps = 1e-9_wp
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

    case (302)
         if (.not. files_loaded) then
            ! Print status message
            print *, "Loading all data files..."
            index_x = i
            index_y = j
            open(newunit=unit, file=trim(fileNames(1)), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, "Error opening coordinate file: ", trim(fileNames(1))
                return
            endif
        
            print *, "Reading coordinates and values from first file..."
            iter = 0
            do iix = 1, xRows
                do iiy = 1, yRows
                    iter = iter + 1
                    read(unit, *, iostat=ios) x_coords(iter), y_coords(iter), stored_values(iix, iiy, 1)
                    if (ios /= 0) then
                        print *, "Error reading coordinates at row ", iter
                        exit
                    endif
                end do
            end do
            close(unit)
            !Now read only the values from remaining files (skip x,y columns)
            do f = 2, nFiles
                open(newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle
                endif
                do iix = 1, xRows
                    do iiy = 1, yRows
                        ! Read and discard x,y, then read the value
                        read(unit, *, iostat=ios) dummy_x, dummy_y, stored_values(iix, iiy, f)
                        if (ios /= 0) then
                            print *, "Error reading file ", trim(fileNames(f)), " at ix=", iix, " iy=", iiy
                            exit
                        endif
                    end do
                end do
            end do
            close(unit)

            domain_xstart = x_coords(1)  ! First x value
            domain_ystart = y_coords(1)  ! First y value
            ny_data = yRows
            nx_data = xRows 
            ! Calculate actual domain end values
            domain_xend = x_coords(nRows - yRows + 1)  ! Last x value
            domain_yend = y_coords(yRows)              ! Last y value in first x-slice
            ! Calculate simulation grid steps
            sim_x_step = x_cc(1) - x_cc(0)
            sim_y_step = y_cc(1) - y_cc(0)
            ! Calculate offsets
            delta_x = x_cc(index_x) - domain_xstart + sim_x_step/2.0_wp
            delta_y = y_cc(index_y) - domain_ystart + sim_y_step/2.0_wp
            ! Global offsets in case of multiple ranks  
            global_offset_x = nint(abs(delta_x) / sim_x_step)
            global_offset_y = nint(abs(delta_y) / sim_y_step)
            
            print *, "Domain X: ", domain_xstart, " to ", domain_xend
            print *, "Domain Y: ", domain_ystart, " to ", domain_yend
            print *, "All ", nFiles, " files loaded successfully!"
            
            files_loaded = .true.
        endif            
        
        data_i = i + 1 + global_offset_x - index_x  ! x-direction index in data grid
        data_j = j + 1 + global_offset_y - index_y  ! y-direction index in data grid

        do f = 1, nFiles 
            if (f > 3) then
                jump = 1
            else 
                jump = 0
            end if
            q_prim_vf(f + jump)%sf(i, j, k) = stored_values(data_i, data_j, f)
        enddo
        ! Set z velocity to zero
        q_prim_vf(4)%sf(i, j, k) = 0.0_wp        ! Put your variable assignments here

    case default  
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
