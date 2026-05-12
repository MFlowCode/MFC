!>
!! @file
!! @brief Contains module m_grid

!> @brief Generates uniform or stretched rectilinear grids with hyperbolic-tangent spacing
module m_grid

    use m_derived_types  ! Definitions of the derived types
    use m_global_parameters  ! Global parameters for the code
    use m_mpi_proxy  ! Message passing interface (MPI) module proxy
    use m_helper_basic
#ifdef MFC_MPI
    use mpi  ! Message passing interface (MPI) module
#endif

    implicit none

    private
    public :: s_initialize_grid_module, s_generate_grid, s_generate_serial_grid, s_generate_parallel_grid, s_finalize_grid_module

    abstract interface

        !> Abstract interface for generating a rectilinear computational grid.
        impure subroutine s_generate_abstract_grid

        end subroutine s_generate_abstract_grid
    end interface

    procedure(s_generate_abstract_grid), pointer :: s_generate_grid => null()

contains

    !> Generate a uniform or stretched rectilinear grid in serial from user parameters.
    impure subroutine s_generate_serial_grid

        ! Generic loop iterator
        integer  :: i, j    !< generic loop operators
        real(wp) :: length  !< domain lengths
        ! Uniform grid: dx = (x_end - x_beg) / (m + 1)

        dx = (x_domain%end - x_domain%beg)/real(m + 1, wp)

        do i = 0, m
            x%cc(i) = x_domain%beg + 5.e-1_wp*dx*real(2*i + 1, wp)
            x%cb(i - 1) = x_domain%beg + dx*real(i, wp)
        end do

        x%cb(m) = x_domain%end

        ! Hyperbolic tangent grid stretching
        if (stretch_x) then
            length = abs(x%cb(m) - x%cb(-1))
            x%cb = x%cb/length
            x_a = x_a/length
            x_b = x_b/length

            do j = 1, loops_x
                do i = -1, m
                    x%cb(i) = x%cb(i)/a_x*(a_x + log(cosh(a_x*(x%cb(i) - x_a))) + log(cosh(a_x*(x%cb(i) - x_b))) &
                         & - 2._wp*log(cosh(a_x*(x_b - x_a)/2._wp)))
                end do
            end do
            x%cb = x%cb*length

            x%cc(0:m) = (x%cb(0:m) + x%cb(-1:m - 1))/2._wp

            dx = minval(x%cb(0:m) - x%cb(-1:m - 1))
            print *, 'Stretched grid: min/max x grid: ', minval(x%cc(:)), maxval(x%cc(:))
            if (num_procs > 1) call s_mpi_reduce_min(dx)
        end if

        ! Grid Generation in the y-direction
        if (n == 0) return

        ! Axisymmetric cylindrical grid (r, z): half-cell offset at r=0 axis
        if (grid_geometry == 2 .and. f_approx_equal(y_domain%beg, 0.0_wp)) then
            dy = (y_domain%end - y_domain%beg)/real(2*n + 1, wp)

            y%cc(0) = y_domain%beg + 5.e-1_wp*dy
            y%cb(-1) = y_domain%beg

            do i = 1, n
                y%cc(i) = y_domain%beg + 2._wp*dy*real(i, wp)
                y%cb(i - 1) = y_domain%beg + dy*real(2*i - 1, wp)
            end do
        else
            dy = (y_domain%end - y_domain%beg)/real(n + 1, wp)

            do i = 0, n
                y%cc(i) = y_domain%beg + 5.e-1_wp*dy*real(2*i + 1, wp)
                y%cb(i - 1) = y_domain%beg + dy*real(i, wp)
            end do
        end if

        y%cb(n) = y_domain%end

        ! Hyperbolic tangent grid stretching in y-direction
        if (stretch_y) then
            length = abs(y%cb(n) - y%cb(-1))
            y%cb = y%cb/length
            y_a = y_a/length
            y_b = y_b/length

            do j = 1, loops_y
                do i = -1, n
                    y%cb(i) = y%cb(i)/a_y*(a_y + log(cosh(a_y*(y%cb(i) - y_a))) + log(cosh(a_y*(y%cb(i) - y_b))) &
                         & - 2._wp*log(cosh(a_y*(y_b - y_a)/2._wp)))
                end do
            end do

            y%cb = y%cb*length
            y%cc(0:n) = (y%cb(0:n) + y%cb(-1:n - 1))/2._wp

            dy = minval(y%cb(0:n) - y%cb(-1:n - 1))

            if (num_procs > 1) call s_mpi_reduce_min(dy)
        end if

        ! Grid Generation in the z-direction
        if (p == 0) return

        dz = (z_domain%end - z_domain%beg)/real(p + 1, wp)

        do i = 0, p
            z%cc(i) = z_domain%beg + 5.e-1_wp*dz*real(2*i + 1, wp)
            z%cb(i - 1) = z_domain%beg + dz*real(i, wp)
        end do

        z%cb(p) = z_domain%end

        ! Hyperbolic tangent grid stretching in z-direction
        if (stretch_z) then
            length = abs(z%cb(p) - z%cb(-1))
            z%cb = z%cb/length
            z_a = z_a/length
            z_b = z_b/length

            do j = 1, loops_z
                do i = -1, p
                    z%cb(i) = z%cb(i)/a_z*(a_z + log(cosh(a_z*(z%cb(i) - z_a))) + log(cosh(a_z*(z%cb(i) - z_b))) &
                         & - 2._wp*log(cosh(a_z*(z_b - z_a)/2._wp)))
                end do
            end do

            z%cb = z%cb*length
            z%cc(0:p) = (z%cb(0:p) + z%cb(-1:p - 1))/2._wp

            dz = minval(z%cb(0:p) - z%cb(-1:p - 1))

            if (num_procs > 1) call s_mpi_reduce_min(dz)
        end if

    end subroutine s_generate_serial_grid

    !> Generate a uniform or stretched rectilinear grid in parallel from user parameters.
    impure subroutine s_generate_parallel_grid

#ifdef MFC_MPI
        real(wp) :: length  !< domain lengths
        ! Locations of cell boundaries
        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb  !< Locations of cell boundaries
        character(LEN=path_len + name_len)  :: file_loc                      !< Generic string used to store the address of a file
        integer                             :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer                             :: i, j                          !< Generic loop integers

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Uniform grid: dx = (x_end - x_beg) / (m_glb + 1)
        dx = (x_domain%end - x_domain%beg)/real(m_glb + 1, wp)
        do i = 0, m_glb
            x_cb_glb(i - 1) = x_domain%beg + dx*real(i, wp)
        end do
        x_cb_glb(m_glb) = x_domain%end
        ! Hyperbolic tangent grid stretching in x-direction (parallel version)
        if (stretch_x) then
            length = abs(x_cb_glb(m_glb) - x_cb_glb(-1))

            x_cb_glb = x_cb_glb/length

            x_a = x_a/length
            x_b = x_b/length

            do j = 1, loops_x
                do i = -1, m_glb
                    x_cb_glb(i) = x_cb_glb(i)/a_x*(a_x + log(cosh(a_x*(x_cb_glb(i) - x_a))) + log(cosh(a_x*(x_cb_glb(i) - x_b))) &
                             & - 2._wp*log(cosh(a_x*(x_b - x_a)/2._wp)))
                end do
            end do

            x_cb_glb = x_cb_glb*length
        end if

        ! Grid generation in the y-direction
        if (n_glb > 0) then
            ! Axisymmetric cylindrical grid (r, z): half-cell offset at r=0 axis
            if (grid_geometry == 2 .and. f_approx_equal(y_domain%beg, 0.0_wp)) then
                dy = (y_domain%end - y_domain%beg)/real(2*n_glb + 1, wp)
                y_cb_glb(-1) = y_domain%beg
                do i = 1, n_glb
                    y_cb_glb(i - 1) = y_domain%beg + dy*real(2*i - 1, wp)
                end do
            else
                dy = (y_domain%end - y_domain%beg)/real(n_glb + 1, wp)
                do i = 0, n_glb
                    y_cb_glb(i - 1) = y_domain%beg + dy*real(i, wp)
                end do
            end if
            y_cb_glb(n_glb) = y_domain%end
            if (stretch_y) then
                length = abs(y_cb_glb(n_glb) - y_cb_glb(-1))

                y_cb_glb = y_cb_glb/length

                y_a = y_a/length
                y_b = y_b/length

                do j = 1, loops_y
                    do i = -1, n_glb
                        y_cb_glb(i) = y_cb_glb(i)/a_y*(a_y + log(cosh(a_y*(y_cb_glb(i) - y_a))) + log(cosh(a_y*(y_cb_glb(i) - y_b) &
                                 & )) - 2._wp*log(cosh(a_y*(y_b - y_a)/2._wp)))
                    end do
                end do

                y_cb_glb = y_cb_glb*length
            end if

            ! Grid generation in the z-direction
            if (p_glb > 0) then
                dz = (z_domain%end - z_domain%beg)/real(p_glb + 1, wp)
                do i = 0, p_glb
                    z_cb_glb(i - 1) = z_domain%beg + dz*real(i, wp)
                end do
                z_cb_glb(p_glb) = z_domain%end
                if (stretch_z) then
                    length = abs(z_cb_glb(p_glb) - z_cb_glb(-1))

                    z_cb_glb = z_cb_glb/length
                    z_a = z_a/length
                    z_b = z_b/length

                    do j = 1, loops_z
                        do i = -1, p_glb
                            z_cb_glb(i) = z_cb_glb(i)/a_z*(a_z + log(cosh(a_z*(z_cb_glb(i) - z_a))) + log(cosh(a_z*(z_cb_glb(i) &
                                     & - z_b))) - 2._wp*log(cosh(a_z*(z_b - z_a)/2._wp)))
                        end do
                    end do

                    z_cb_glb = z_cb_glb*length
                end if
            end if
        end if

        ! Write cell boundary locations to grid data files
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'x_cb.dat'
        data_size = m_glb + 2
        call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
        call MPI_FILE_WRITE(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
        call MPI_FILE_CLOSE(ifile, ierr)

        if (n > 0) then
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'y_cb.dat'
            data_size = n_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)

            if (p > 0) then
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'z_cb.dat'
                data_size = p_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
                call MPI_FILE_WRITE(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)
#endif

    end subroutine s_generate_parallel_grid

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_grid_module

        if (parallel_io .neqv. .true.) then
            s_generate_grid => s_generate_serial_grid
        else
            s_generate_grid => s_generate_parallel_grid
        end if

    end subroutine s_initialize_grid_module

    !> Deallocation procedures for the module
    impure subroutine s_finalize_grid_module

        s_generate_grid => null()

    end subroutine s_finalize_grid_module

end module m_grid
