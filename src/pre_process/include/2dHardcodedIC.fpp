#:def Hardcoded2DVariables()

    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: factor

    integer, parameter :: nFiles = 14   ! Can be changed to any number
    integer, parameter :: nRows = 401!
    integer :: f, iter, ios, unit, idx, jump, index_1
    real(wp) :: x_len, x_step
    integer :: global_offset
    real(wp) :: delta_x
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
            write (file_num_str, '(I1)') f  ! Single digit
        else
            write (file_num_str, '(I2)') f  ! Double digit
            ! For more than 99 files, you might need to adjust this format
        end if
        fileNames(f) = trim(init_dir)//"prim."//trim(file_num_str)//".00."//zeros_default//".dat"
    end do

    eps = 1e-9_wp

#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case

    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1._wp/3._wp)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000._wp
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - eps)*1._wp
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000._wp
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2._wp))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5_wp)**2._wp + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4._wp*(1._wp - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2._wp*(-2._wp + 4*log(2._wp))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1._wp/gam)
    case (204) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.05_wp/wl

        intH = amp*sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + h

        alph = 0.5_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (205) ! 2D lung wave interaction problem
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)         !to be changed later!       !non dim amplitude

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        if (y_cc(j) > intH) then
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (206) ! 2D lung wave interaction problem - horizontal domain
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)

        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2) + h

        if (x_cc(i) > intL) then        !this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (207)
        if (.not. files_loaded) then
            ! Print status message
            index_1 = i
            do f = 1, nFiles
                ! Open the file for reading
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle  ! Skip this file on error
                end if
                ! Read all rows at once into memory
                do iter = 1, nRows
                    read (unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                    if (ios /= 0) then
                        print *, "Error reading file ", trim(fileNames(f)), " at row ", iter
                        exit  ! Exit loop on error
                    end if
                end do
                close (unit)
            end do
            ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = x_cc(1) - x_cc(0)
            delta_x = x_cc(index_1) - domain_start + x_step/2
            global_offset = nint(abs(delta_x)/x_step)
            print *, "All data files loaded. Domain range:", domain_start, "to", domain_end
            files_loaded = .true.
        end if
        ! Calculate the index in the file data array corresponding to x_cc(i)
        idx = i + 1 + global_offset - index_1
        ! Assign values from stored data for each file (with your small adjustment for f > 2)
        do f = 1, nFiles
            if (f > 2) then
                jump = 1
            else
                jump = 0
            end if
            q_prim_vf(f + jump)%sf(i, j, 0) = stored_values(idx, f)
        end do
        ! Set element 3 (perpedicular velocity v) explicitly to zero
        q_prim_vf(3)%sf(i, j, 0) = 0.0_wp

    case (250) ! MHD Orszag-Tang vortex
        ! gamma = 5/3
        !   rho = 25/(36*pi)
        !     p = 5/(12*pi)
        !     v = (-sin(2*pi*y), sin(2*pi*x), 0)
        !     B = (-sin(2*pi*y)/sqrt(4*pi), sin(4*pi*x)/sqrt(4*pi), 0)

        q_prim_vf(momxb)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = sin(2._wp*pi*x_cc(i))

        q_prim_vf(B_idx%beg)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))/sqrt(4._wp*pi)
        q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = sin(4._wp*pi*x_cc(i))/sqrt(4._wp*pi)

    case (251) ! RMHD Cylindrical Blast Wave [Mignone, 2006: Section 4.3.1]

        if (x_cc(i)**2 + y_cc(j)**2 < 0.08_wp**2) then
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0
        elseif (x_cc(i)**2 + y_cc(j)**2 <= 1._wp**2) then
            ! Linear interpolation between r=0.08 and r=1.0
            factor = (1.0_wp - sqrt(x_cc(i)**2 + y_cc(j)**2))/(1.0_wp - 0.08_wp)
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01_wp*factor + 1.e-4_wp*(1.0_wp - factor)
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0_wp*factor + 3.e-5_wp*(1.0_wp - factor)
        else
            q_prim_vf(contxb)%sf(i, j, 0) = 1.e-4_wp
            q_prim_vf(E_idx)%sf(i, j, 0) = 3.e-5_wp
        end if

    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
