#:include 'macros.fpp'

module m_sim_helpers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    use m_data_output

    use m_compile_specific
    ! ==========================================================================

    implicit none

    private; public :: s_compute_enthalpy, &
 s_read_data_files, &
 s_read_serial_data_files, &
 s_read_parallel_data_files, &
 s_compute_stability_from_dt, &
 s_compute_dt_from_cfl, &
 s_comprehensive_debug, &
 s_open_run_time_information_file, &
 s_write_run_time_information, &
 s_close_run_time_information_file, &
 s_initialize_sim_helpers_module, &
 s_finalize_sim_helpers_module

    abstract interface ! ===================================================

        !! @param q_cons_vf  Conservative variables
        subroutine s_read_abstract_data_files(q_cons_vf, ib_markers)

            import :: scalar_field, integer_field, sys_size, pres_field

            type(scalar_field), &
                dimension(sys_size), &
                intent(inout) :: q_cons_vf

            type(integer_field) :: ib_markers

        end subroutine s_read_abstract_data_files

    end interface ! ========================================================

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), icfl_sf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), vcfl_sf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), ccfl_sf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), Rc_sf)
    !$acc declare link(icfl_sf, vcfl_sf, ccfl_sf, Rc_sf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: icfl_sf  !< ICFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: vcfl_sf  !< VCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: ccfl_sf  !< CCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: Rc_sf  !< Rc stability criterion
    !$acc declare create(icfl_sf, vcfl_sf, ccfl_sf, Rc_sf)
#endif

    real(kind(0d0)) :: icfl_max_loc, icfl_max_glb !< ICFL stability extrema on local and global grids
    real(kind(0d0)) :: vcfl_max_loc, vcfl_max_glb !< VCFL stability extrema on local and global grids
    real(kind(0d0)) :: ccfl_max_loc, ccfl_max_glb !< CCFL stability extrema on local and global grids
    real(kind(0d0)) :: Rc_min_loc, Rc_min_glb !< Rc   stability extrema on local and global grids
    !$acc declare create(icfl_max_loc, icfl_max_glb, vcfl_max_loc, vcfl_max_glb, ccfl_max_loc, ccfl_max_glb, Rc_min_loc, Rc_min_glb)

    !> @name ICFL, VCFL, CCFL and Rc stability criteria extrema over all the time-steps
    !> @{
    real(kind(0d0)) :: icfl_max !< ICFL criterion maximum
    real(kind(0d0)) :: vcfl_max !< VCFL criterion maximum
    real(kind(0d0)) :: ccfl_max !< CCFL criterion maximum
    real(kind(0d0)) :: Rc_min !< Rc criterion maximum
    !> @}

contains

    subroutine s_initialize_sim_helpers_module()

        ! Allocating/initializing ICFL, VCFL, CCFL and Rc stability criteria
        @:ALLOCATE_GLOBAL(icfl_sf(0:m, 0:n, 0:p))
        icfl_max = 0d0

        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(vcfl_sf(0:m, 0:n, 0:p))
            @:ALLOCATE_GLOBAL(Rc_sf  (0:m, 0:n, 0:p))

            vcfl_max = 0d0
            Rc_min = 1d3
        end if

        ! Associate pointers for serial or parallel I/O
        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
            s_write_data_files => s_write_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
            s_write_data_files => s_write_parallel_data_files
        end if

    end subroutine s_initialize_sim_helpers_module

    !> Computes enthalpy
        !! @param q_prim_vf cell centered primitive variables
        !! @param pres mixture pressure
        !! @param rho mixture density
        !! @param gamma mixture gamma
        !! @param pi_inf mixture pi_inf
        !! @param Re mixture reynolds number
        !! @param H mixture enthalpy
        !! @param alpha component alphas
        !! @param vel directional velocities
        !! @param vel_sum squard sum of velocity components
        !! @param j x index
        !! @param k y index
        !! @param l z index
    subroutine s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)
        !$acc routine seq
        type(scalar_field), dimension(sys_size) :: q_prim_vf
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: rho, gamma, pi_inf, qv, vel_sum, E, H, pres
        real(kind(0d0)), dimension(2) :: Re
        integer :: i, j, k, l

        do i = 1, num_fluids
            alpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
            alpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
        end do

        if (bubbles) then
            call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
        else
            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
        end if

        do i = 1, num_dims
            vel(i) = q_prim_vf(contxe + i)%sf(j, k, l)
        end do

        vel_sum = 0d0
        do i = 1, num_dims
            vel_sum = vel_sum + vel(i)**2d0
        end do

        pres = q_prim_vf(E_idx)%sf(j, k, l)

        E = gamma*pres + pi_inf + 5d-1*rho*vel_sum + qv

        H = (E + pres)/rho

    end subroutine s_compute_enthalpy

    !> Computes stability criterion for a specified dt
        !! @param vel directional velocities
        !! @param c mixture speed of sound
        !! @param Re_l mixture Reynolds number
        !! @param j x index
        !! @param k y index
        !! @param l z index
        !! @param icfl_sf cell centered inviscid cfl number
        !! @param vcfl_sf (optional) cell centered viscous cfl number
        !! @param Rc_sf (optional) cell centered Rc
    subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
        !$acc routine seq
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: c, icfl_dt, vcfl_dt, rho
        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: icfl_sf
        real(kind(0d0)), dimension(0:m, 0:n, 0:p), optional :: vcfl_sf, Rc_sf
        real(kind(0d0)) :: fltr_dtheta   !<
             !! Modified dtheta accounting for Fourier filtering in azimuthal direction.
        integer :: j, k, l
        integer :: Nfq
        real(kind(0d0)), dimension(2) :: Re_l

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2d0*pi*y_cb(0)/3d0
            elseif (k <= fourier_rings) then
                Nfq = min(floor(2d0*real(k, kind(0d0))*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2d0*pi*y_cb(k - 1)/real(Nfq, kind(0d0))
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        end if

        if (p > 0) then
            !3D
            if (grid_geometry == 3) then
                icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                          dy(k)/(abs(vel(2)) + c), &
                                          fltr_dtheta/(abs(vel(3)) + c))
            else
                icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                          dy(k)/(abs(vel(2)) + c), &
                                          dz(l)/(abs(vel(3)) + c))
            end if

            if (any(Re_size > 0)) then

                if (grid_geometry == 3) then
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho) &
                                       /min(dx(j), dy(k), fltr_dtheta)**2d0

                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                         dy(k)*(abs(vel(2)) + c), &
                                         fltr_dtheta*(abs(vel(3)) + c)) &
                                     /maxval(1d0/Re_l)
                else
                    vcfl_sf(j, k, l) = maxval(dt/Re_l/rho) &
                                       /min(dx(j), dy(k), dz(l))**2d0

                    Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                         dy(k)*(abs(vel(2)) + c), &
                                         dz(l)*(abs(vel(3)) + c)) &
                                     /maxval(1d0/Re_l)
                end if

            end if

        elseif (n > 0) then
            !2D
            icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                      dy(k)/(abs(vel(2)) + c))

            if (any(Re_size > 0)) then

                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/min(dx(j), dy(k))**2d0

                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                     dy(k)*(abs(vel(2)) + c)) &
                                 /maxval(1d0/Re_l)

            end if

        else
            !1D
            icfl_sf(j, k, l) = (dt/dx(j))*(abs(vel(1)) + c)

            if (any(Re_size > 0)) then

                vcfl_sf(j, k, l) = maxval(dt/Re_l/rho)/dx(j)**2d0

                Rc_sf(j, k, l) = dx(j)*(abs(vel(1)) + c)/maxval(1d0/Re_l)

            end if

        end if

    end subroutine s_compute_stability_from_dt

    !> Computes dt for a specified CFL number
        !! @param vel directional velocities
        !! @param max_dt cell centered maximum dt
        !! @param rho cell centered density
        !! @param Re_l cell centered Reynolds number
        !! @param j x coordinate
        !! @param k y coordinate
        !! @param l z coordinate
    subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, j, k, l)
        !$acc routine seq
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: c, icfl_dt, vcfl_dt, rho
        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: max_dt
        real(kind(0d0)) :: fltr_dtheta   !<
             !! Modified dtheta accounting for Fourier filtering in azimuthal direction.
        integer :: j, k, l
        integer :: Nfq
        real(kind(0d0)), dimension(2) :: Re_l

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2d0*pi*y_cb(0)/3d0
            elseif (k <= fourier_rings) then
                Nfq = min(floor(2d0*real(k, kind(0d0))*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2d0*pi*y_cb(k - 1)/real(Nfq, kind(0d0))
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        end if

        if (p > 0) then
            !3D
            if (grid_geometry == 3) then
                icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                         dy(k)/(abs(vel(2)) + c), &
                                         fltr_dtheta/(abs(vel(3)) + c))
            else
                icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                         dy(k)/(abs(vel(2)) + c), &
                                         dz(l)/(abs(vel(3)) + c))
            end if

            if (any(Re_size > 0)) then
                if (grid_geometry == 3) then
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), fltr_dtheta)**2d0) &
                              /minval(1/(rho*Re_l))
                else
                    vcfl_dt = cfl_target*(min(dx(j), dy(k), dz(l))**2d0) &
                              /minval(1/(rho*Re_l))
                end if
            end if

        elseif (n > 0) then
            !2D
            icfl_dt = cfl_target*min(dx(j)/(abs(vel(1)) + c), &
                                     dy(k)/(abs(vel(2)) + c))

            if (any(Re_size > 0)) then
                vcfl_dt = cfl_target*(min(dx(j), dy(k))**2d0)/maxval((1/Re_l)/rho)
            end if

        else
            !1D
            icfl_dt = cfl_target*(dx(j)/(abs(vel(1)) + c))

            if (any(Re_size > 0)) then
                vcfl_dt = cfl_target*(dx(j)**2d0)/minval(1/(rho*Re_l))
            end if

        end if

        if (any(re_size > 0)) then
            max_dt(j, k, l) = min(icfl_dt, vcfl_dt)
        else
            max_dt(j, k, l) = icfl_dt
        end if

    end subroutine s_compute_dt_from_cfl

    subroutine s_comprehensive_debug(q_cons_vf, q_prim_vf, t_step, stage)

        type(scalar_field), dimension(sys_size) :: q_cons_vf, q_prim_vf
        integer, intent(in) :: t_step, stage
        integer :: j, k, l, i
        integer errors
        logical :: exists

        character(LEN=name_len) :: file_name = 'comp_debug.txt'
        character(LEN=path_len + name_len) :: file_path
        character(100) :: str_format

        ! Opening the run-time information file
        file_path = trim(case_dir)//'/'//trim(file_name)

        str_format = "(I9, A, I3, A, I4, I4, I4, A, I2, A, I5, A, I5, I5, I5)"

        open (12, FILE=trim(file_path), &
          STATUS='replace')

        errors = 0

        ! Check all variables for NaNs
        do i = 1, sys_size
            !$acc update host(q_cons_vf(i)%sf)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                            write(12, str_format) t_step, " NaN(s) in conservative variables after RK stage ", &
                                stage, " at (j,k,l) ", j, k, l, " equation", i, " proc", proc_rank, &
                                " (m, n, p)", m, n, p
                            errors = errors + 1
                        end if
                    end do
                end do
            end do
        end do

        ! Check for invalid volume fractions
        do i = advxb, advxe
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (q_cons_vf(i)%sf(j, k, l) < 0d0) then
                            write(12, str_format) t_step, " Volume fraction < 0 after RK stage ", &
                                stage, " at (j,k,l) ", j, k, l, " equation", i, " proc", proc_rank, &
                                " (m, n, p)", m, n, p
                            errors = errors + 1
                        elseif (q_cons_vf(i)%sf(j, k, l) > 1d0 + verysmall) then
                            write(12, str_format) t_step, " Volume fraction > 1 after RK stage ", &
                                stage, " at (j,k,l) ", j, k, l, " equation", i, " proc", proc_rank, &
                                " (m, n, p)", m, n, p
                            errors = errors + 1
                        end if
                    end do
                end do
            end do
        end do

        ! Check for invalid densities
        do i = contxb, contxe
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (q_cons_vf(advxb + i -1)%sf(j, k, l) < 0d0 .and. q_cons_vf(i)%sf(j, k, l) < 0d0 .or. &
                            q_cons_vf(advxb + i -1)%sf(j, k, l) > 0d0 .and. q_cons_Vf(i)%sf(j, k, l) < 0d0) then
                            print*, q_cons_vf(advxb + i - 1)%sf(j, k, l), q_cons_vf(i)%sf(j, k, l)
                            write(12, str_format) t_step, " Density is negative after RK stage ", &
                                stage, " at (j,k,l) ", j, k, l, " equation", i, " proc", proc_rank, &
                                " (m, n, p)", m, n, p
                            errors = errors + 1
                        end if
                    end do
                end do
            end do
        end do

        if (errors /= 0) then
            close(12)
            call s_write_data_files(q_cons_vf, q_prim_vf, t_step)
            call s_mpi_abort("Errors found in conservative variables")
        endif

        write(12, "(I3)") -1
        close(12)

    end subroutine s_comprehensive_debug

    !!              initial condition and grid data files. The cell-average
        !!              conservative variables constitute the former, while the
        !!              cell-boundary locations in x-, y- and z-directions make
        !!              up the latter. This procedure also calculates the cell-
        !!              width distributions from the cell-boundary locations.
        !! @param q_cons_vf Cell-averaged conservative variables
    subroutine s_read_serial_data_files(q_cons_vf, ib_markers)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(integer_field) :: ib_markers

        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the starting time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
        ! Logical used to check the existence of the data files

        integer :: i, r !< Generic loop iterator

        ! Confirming that the directory from which the initial condition and
        ! the grid data files are to be read in exists and exiting otherwise
        if (cfl_dt) then
            write (t_step_dir, '(A,I0,A,I0)') &
                trim(case_dir)//'/p_all/p', proc_rank, '/', n_start
        else
            write (t_step_dir, '(A,I0,A,I0)') &
                trim(case_dir)//'/p_all/p', proc_rank, '/', t_step_start
        end if

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

        ! Cell-boundary Locations in x-direction ===========================
        file_path = trim(t_step_dir)//'/x_cb.dat'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  ACTION='read', &
                  STATUS='old')
            read (2) x_cb(-1:m); close (2)
        else
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0))*20) + 1
                end if
            end do
        end if
        ! ==================================================================

        ! Cell-boundary Locations in y-direction ===========================
        if (n > 0) then

            file_path = trim(t_step_dir)//'/y_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) y_cb(-1:n); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
            end if

            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

        end if
        ! ==================================================================

        ! Cell-boundary Locations in z-direction ===========================
        if (p > 0) then

            file_path = trim(t_step_dir)//'/z_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) z_cb(-1:p); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
            end if

            dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
            z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0

        end if
        ! ==================================================================

        do i = 1, sys_size
            write (file_path, '(A,I0,A)') &
                trim(t_step_dir)//'/q_cons_vf', i, '.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) q_cons_vf(i)%sf(0:m, 0:n, 0:p); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
            end if
        end do

        if ((bubbles .eqv. .true.) .or. (hypoelasticity .eqv. .true.)) then
            ! Read pb and mv for non-polytropic qbmm
            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/pb', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) pb_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                        end if
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/mv', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) mv_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                        end if
                    end do
                end do
            end if
        end if
        ! ==================================================================

        ! Read IBM Data ====================================================

        if (ib) then
            do i = 1, num_ibs
                write (file_path, '(A,I0,A)') &
                    trim(t_step_dir)//'/ib.dat'
                inquire (FILE=trim(file_path), EXIST=file_exist)
                if (file_exist) then
                    open (2, FILE=trim(file_path), &
                          FORM='unformatted', &
                          ACTION='read', &
                          STATUS='old')
                    read (2) ib_markers%sf(0:m, 0:n, 0:p); close (2)
                else
                    call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                end if

                if (patch_ib(i)%c > 0) then

                    print *, "HERE Np", Np
                    allocate (airfoil_grid_u(1:Np))
                    allocate (airfoil_grid_l(1:Np))

                    write (file_path, '(A)') &
                        trim(t_step_dir)//'/airfoil_u.dat'
                    inquire (FILE=trim(file_path), EXIST=file_exist)
                    if (file_exist) then
                        open (2, FILE=trim(file_path), &
                              FORM='unformatted', &
                              ACTION='read', &
                              STATUS='old')
                        read (2) airfoil_grid_u; close (2)
                    else
                        call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                    end if

                    write (file_path, '(A)') &
                        trim(t_step_dir)//'/airfoil_l.dat'
                    inquire (FILE=trim(file_path), EXIST=file_exist)
                    if (file_exist) then
                        open (2, FILE=trim(file_path), &
                              FORM='unformatted', &
                              ACTION='read', &
                              STATUS='old')
                        read (2) airfoil_grid_l; close (2)
                    else
                        call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                    end if
                end if
            end do

        end if

    end subroutine s_read_serial_data_files

        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_data_files(q_cons_vf, ib_markers)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf
        type(integer_field) :: ib_markers

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        character(len=10) :: t_step_start_string

        integer :: i, j

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center locations
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0))*20) + 1
                    allocate (MPI_IO_airfoil_IB_DATA%var(1:2*Np))
                    print *, "HERE Np", Np
                end if
            end do
        end if

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center locations
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File '//trim(file_loc)//'is missing. Exiting...')
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing the cell width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                ! Computing the cell center locations
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0

            end if
        end if

        if (file_per_process) then
            if (cfl_dt) then
                call s_int_to_str(n_start, t_step_start_string)
                write (file_loc, '(I0,A1,I7.7,A)') n_start, '_', proc_rank, '.dat'
            else
                call s_int_to_str(t_step_start, t_step_start_string)
                write (file_loc, '(I0,A1,I7.7,A)') t_step_start, '_', proc_rank, '.dat'
            end if

            file_loc = trim(case_dir)//'/restart_data/lustre_'//trim(t_step_start_string)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                ! Initialize MPI data I/O

                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers)
                else
                    call s_initialize_mpi_data(q_cons_vf)
                end if

                ! Size of local arrays
                data_size = (m + 1)*(n + 1)*(p + 1)

                ! Resize some integers so MPI can read even the biggest file
                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                WP_MOK = int(8d0, MPI_OFFSET_KIND)
                MOK = int(1d0, MPI_OFFSET_KIND)
                str_MOK = int(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                if (bubbles .or. hypoelasticity) then

                    do i = 1, sys_size!adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                    !Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               MPI_DOUBLE_PRECISION, status, ierr)
                        end do
                    end if
                else
                    do i = 1, adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)

                if (ib) then

                    write (file_loc, '(A)') 'ib.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                           MPI_INTEGER, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
                    end if

                end if

            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
            end if
        else
            ! Open the file to read conservative variables
            if (cfl_dt) then
                 write (file_loc, '(I0,A)') n_start, '.dat'
            else
                write (file_loc, '(I0,A)') t_step_start, '.dat'
            end if
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                ! Initialize MPI data I/O

                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers)
                else

                    call s_initialize_mpi_data(q_cons_vf)

                end if


                ! Size of local arrays
                data_size = (m + 1)*(n + 1)*(p + 1)

                ! Resize some integers so MPI can read even the biggest file
                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                WP_MOK = int(8d0, MPI_OFFSET_KIND)
                MOK = int(1d0, MPI_OFFSET_KIND)
                str_MOK = int(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                if (bubbles .or. hypoelasticity) then

                    do i = 1, sys_size!adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                    !Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)
                            ! Initial displacement to skip at beginning of file
                            disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                            call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                                   'native', mpi_info_int, ierr)
                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               MPI_DOUBLE_PRECISION, status, ierr)
                        end do
                    end if
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)

                    end do
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)

                if (ib) then

                    write (file_loc, '(A)') 'ib.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                           MPI_INTEGER, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
                    end if

                end if

            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
            end if

        end if

        if (ib) then

            do j = 1, num_ibs
                if (patch_ib(j)%c > 0) then

                    print *, "HERE Np", Np

                    allocate (airfoil_grid_u(1:Np))
                    allocate (airfoil_grid_l(1:Np))

                    write (file_loc, '(A)') 'airfoil_l.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        ! Initial displacement to skip at beginning of file
                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(1), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_airfoil_IB_DATA%var(1:Np), 3*Np, &
                                           MPI_DOUBLE_PRECISION, status, ierr)

                    end if

                    write (file_loc, '(A)') 'airfoil_u.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        ! Initial displacement to skip at beginning of file
                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_airfoil_IB_DATA%view(2), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_airfoil_IB_DATA%var(Np + 1:2*Np), 3*Np, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                    end if

                    do i = 1, Np
                        airfoil_grid_l(i)%x = MPI_IO_airfoil_IB_DATA%var(i)%x
                        airfoil_grid_l(i)%y = MPI_IO_airfoil_IB_DATA%var(i)%y
                    end do

                    do i = 1, Np
                        airfoil_grid_u(i)%x = MPI_IO_airfoil_IB_DATA%var(Np + i)%x
                        airfoil_grid_u(i)%y = MPI_IO_airfoil_IB_DATA%var(Np + i)%y
                    end do

                end if
            end do
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_read_parallel_data_files

    !>  The purpose of this subroutine is to open a new or pre-
        !!          existing run-time information file and append to it the
        !!      basic header information relevant to current simulation.
        !!      In general, this requires generating a table header for
        !!      those stability criteria which will be written at every
        !!      time-step.
    subroutine s_open_run_time_information_file

        character(LEN=name_len) :: file_name = 'run_time.inf' !<
            !! Name of the run-time information file

        character(LEN=path_len + name_len) :: file_path !<
            !! Relative path to a file in the case directory

        character(LEN=8) :: file_date !<
            !! Creation date of the run-time information file

        ! Opening the run-time information file
        file_path = trim(case_dir)//'/'//trim(file_name)

        open (3, FILE=trim(file_path), &
              FORM='formatted', &
              STATUS='replace')

        write (3, '(A)') 'Description: Stability information at '// &
            'each time-step of the simulation. This'
        write (3, '(13X,A)') 'data is composed of the inviscid '// &
            'Courant–Friedrichs–Lewy (ICFL)'
        write (3, '(13X,A)') 'number, the viscous CFL (VCFL) number, '// &
            'the capillary CFL (CCFL)'
        write (3, '(13X,A)') 'number and the cell Reynolds (Rc) '// &
            'number. Please note that only'
        write (3, '(13X,A)') 'those stability conditions pertinent '// &
            'to the physics included in'
        write (3, '(13X,A)') 'the current computation are displayed.'

        call date_and_time(DATE=file_date)

        write (3, '(A)') 'Date: '//file_date(5:6)//'/'// &
            file_date(7:8)//'/'// &
            file_date(3:4)

        write (3, '(A)') ''; write (3, '(A)') ''

        ! Generating table header for the stability criteria to be outputted
        if (cfl_dt) then
            if (any(Re_size > 0)) then
                write (1, '(A)') '==== Time-steps ====== dt ===== Time ======= ICFL '// &
                    'Max ==== VCFL Max ====== Rc Min ======='
            else
                write (1, '(A)') '=========== Time-steps ============== dt ===== Time '// &
                    '============== ICFL Max ============='
            end if
        else
            if (any(Re_size > 0)) then
                write (1, '(A)') '==== Time-steps ====== Time ======= ICFL '// &
                    'Max ==== VCFL Max ====== Rc Min ======='
            else
                write (1, '(A)') '=========== Time-steps ============== Time '// &
                    '============== ICFL Max ============='
            end if
        end if

    end subroutine s_open_run_time_information_file

    !>  The goal of the procedure is to output to the run-time
        !!      information file the stability criteria extrema in the
        !!      entire computational domain and at the given time-step.
        !!      Moreover, the subroutine is also in charge of tracking
        !!      these stability criteria extrema over all time-steps.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time step
    subroutine s_write_run_time_information(q_prim_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer, intent(IN) :: t_step

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho  !< Cell-avg. partial density
        real(kind(0d0)) :: rho        !< Cell-avg. density
        real(kind(0d0)), dimension(num_dims) :: vel        !< Cell-avg. velocity
        real(kind(0d0)) :: vel_sum    !< Cell-avg. velocity sum
        real(kind(0d0)) :: pres       !< Cell-avg. pressure
        real(kind(0d0)), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(kind(0d0)) :: gamma      !< Cell-avg. sp. heat ratio
        real(kind(0d0)) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(kind(0d0)) :: qv         !< Cell-avg. fluid reference energy
        real(kind(0d0)) :: c          !< Cell-avg. sound speed
        real(kind(0d0)) :: E          !< Cell-avg. energy
        real(kind(0d0)) :: H          !< Cell-avg. enthalpy
        real(kind(0d0)), dimension(2) :: Re         !< Cell-avg. Reynolds numbers

        ! ICFL, VCFL, CCFL and Rc stability criteria extrema for the current
        ! time-step and located on both the local (loc) and the global (glb)
        ! computational domains

        real(kind(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

        integer :: i, j, k, l, q !< Generic loop iterators

        integer :: Nfq
        real(kind(0d0)) :: fltr_dtheta   !<
            !! Modified dtheta accounting for Fourier filtering in azimuthal direction.

        ! Computing Stability Criteria at Current Time-step ================
        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho, vel, alpha, Re, fltr_dtheta, Nfq)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, c)

                    if (any(Re_size > 0)) then
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
                    else
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf)
                    end if

                end do
            end do
        end do
        ! END: Computing Stability Criteria at Current Time-step ===========

        ! Determining local stability criteria extrema at current time-step

#ifdef CRAY_ACC_WAR
        !$acc update host(icfl_sf)

        if (any(Re_size > 0)) then
            !$acc update host(vcfl_sf, Rc_sf)
        end if

        icfl_max_loc = maxval(icfl_sf)

        if (any(Re_size > 0)) then
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
        end if
#else
        !$acc kernels
        icfl_max_loc = maxval(icfl_sf)
        !$acc end kernels

        if (any(Re_size > 0)) then
            !$acc kernels
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
            !$acc end kernels
        end if
#endif

        ! Determining global stability criteria extrema at current time-step
        if (num_procs > 1) then
            call s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                         vcfl_max_loc, &
                                                         ccfl_max_loc, &
                                                         Rc_min_loc, &
                                                         icfl_max_glb, &
                                                         vcfl_max_glb, &
                                                         ccfl_max_glb, &
                                                         Rc_min_glb)
        else
            icfl_max_glb = icfl_max_loc
            if (any(Re_size > 0)) vcfl_max_glb = vcfl_max_loc
            if (any(Re_size > 0)) Rc_min_glb = Rc_min_loc
        end if

        ! Determining the stability criteria extrema over all the time-steps
        if (icfl_max_glb > icfl_max) icfl_max = icfl_max_glb

        if (any(Re_size > 0)) then
            if (vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
            if (Rc_min_glb < Rc_min) Rc_min = Rc_min_glb
        end if

        ! Outputting global stability criteria extrema at current time-step
        if (proc_rank == 0) then
            if (any(Re_size > 0)) then
                write (1, '(6X,I8,F10.6,6X,6X,F10.6,6X,F9.6,6X,F9.6,6X,F10.6)') &
                    t_step, dt, t_step*dt, icfl_max_glb, &
                    vcfl_max_glb, &
                    Rc_min_glb
            else
                write (1, '(13X,I8,14X,F10.6,14X,F10.6,13X,F9.6)') &
                    t_step, dt, t_step*dt, icfl_max_glb
            end if

            if (icfl_max_glb /= icfl_max_glb) then
                call s_mpi_abort('ICFL is NaN. Exiting ...')
            elseif (icfl_max_glb > 1d0) then
                print *, 'icfl', icfl_max_glb
                call s_mpi_abort('ICFL is greater than 1.0. Exiting ...')
            end if

            do i = chemxb, chemxe
                !@:ASSERT(all(q_prim_vf(i)%sf(:,:,:) >= -1d0), "bad conc")
                !@:ASSERT(all(q_prim_vf(i)%sf(:,:,:) <=  2d0), "bad conc")
            end do

            if (vcfl_max_glb /= vcfl_max_glb) then
                call s_mpi_abort('VCFL is NaN. Exiting ...')
            elseif (vcfl_max_glb > 1d0) then
                print *, 'vcfl', vcfl_max_glb
                call s_mpi_abort('VCFL is greater than 1.0. Exiting ...')
            end if
        end if

        call s_mpi_barrier()

    end subroutine s_write_run_time_information

    !>  The goal of this subroutine is to write to the run-time
        !!      information file basic footer information applicable to
        !!      the current computation and to close the file when done.
        !!      The footer contains the stability criteria extrema over
        !!      all of the time-steps and the simulation run-time.
    subroutine s_close_run_time_information_file

        real(kind(0d0)) :: run_time !< Run-time of the simulation
        ! Writing the footer of and closing the run-time information file
        write (3, '(A)') '----------------------------------------'// &
            '----------------------------------------'
        write (3, '(A)') ''

        write (3, '(A,F9.6)') 'ICFL Max: ', icfl_max
        if (any(Re_size > 0)) write (3, '(A,F9.6)') 'VCFL Max: ', vcfl_max
        if (any(Re_size > 0)) write (3, '(A,F10.6)') 'Rc Min: ', Rc_min

        call cpu_time(run_time)

        write (3, '(A)') ''
        write (3, '(A,I0,A)') 'Run-time: ', int(anint(run_time)), 's'
        write (3, '(A)') '========================================'// &
            '========================================'
        close (3)

    end subroutine s_close_run_time_information_file

    subroutine s_finalize_sim_helpers_module()

        ! Deallocating the ICFL, VCFL, CCFL, and Rc stability criteria
        @:DEALLOCATE_GLOBAL(icfl_sf)
        if (any(Re_size > 0)) then
            @:DEALLOCATE_GLOBAL(vcfl_sf, Rc_sf)
        end if

        s_read_data_files => null()
        s_write_data_files => null()

    end subroutine s_finalize_sim_helpers_module

end module m_sim_helpers
