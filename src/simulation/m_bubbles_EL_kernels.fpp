!>
!! @file m_bubbles_EL_kernels.f90
!! @brief Contains module m_bubbles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian bubbles
!!        in the Eulerian framework.
module m_bubbles_EL_kernels

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    implicit none

contains

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      bubbles into the Eulerian framework using different approaches.
            !! @param nBubs Number of lagrangian bubbles in the current domain
            !! @param lbk_rad Radius of the bubbles
            !! @param lbk_vel Interface velocity of the bubbles
            !! @param lbk_s Computational coordinates of the bubbles
            !! @param lbk_pos Spatial coordinates of the bubbles
            !! @param updatedvar Eulerian variable to be updated
    subroutine s_smoothfunction(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        smoothfunc:select case(lag_params%smooth_type)
        case (1)
        call s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar)
        case (2)
        call s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, updatedvar)
        end select smoothfunc

    end subroutine s_smoothfunction

    !> The purpose of this procedure contains the algorithm to use the delta kernel function to map the effect of the bubbles.
            !!      The effect of the bubbles only affects the cell where the bubble is located.
    subroutine s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, updatedvar)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        integer, dimension(3) :: cell
        real(wp) :: strength_vel, strength_vol

        real(wp) :: addFun1, addFun2, addFun3
        real(wp) :: volpart, Vol
        real(wp), dimension(3) :: s_coord
        integer :: l

        $:GPU_PARALLEL_LOOP(private='[l,s_coord,cell]')
        do l = 1, nBubs

            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            call s_get_cell(s_coord, cell)

            strength_vol = volpart
            strength_vel = 4._wp*pi*lbk_rad(l, 2)**2._wp*lbk_vel(l, 2)

            if (num_dims == 2) then
                Vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            end if

            !Update void fraction field
            addFun1 = strength_vol/Vol
            $:GPU_ATOMIC(atomic='update')
            updatedvar(1)%sf(cell(1), cell(2), cell(3)) = updatedvar(1)%sf(cell(1), cell(2), cell(3)) + real(addFun1, kind=stp)

            !Update time derivative of void fraction
            addFun2 = strength_vel/Vol
            $:GPU_ATOMIC(atomic='update')
            updatedvar(2)%sf(cell(1), cell(2), cell(3)) = updatedvar(2)%sf(cell(1), cell(2), cell(3)) + real(addFun2, kind=stp)

            !Product of two smeared functions
            !Update void fraction * time derivative of void fraction
            if (lag_params%cluster_type >= 4) then
                addFun3 = (strength_vol*strength_vel)/Vol
                $:GPU_ATOMIC(atomic='update')
                updatedvar(5)%sf(cell(1), cell(2), cell(3)) = updatedvar(5)%sf(cell(1), cell(2), cell(3)) + real(addFun3, kind=stp)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_deltafunc

    !> The purpose of this procedure contains the algorithm to use the gaussian kernel function to map the effect of the bubbles.
            !!      The effect of the bubbles affects the 3X3x3 cells that surround the bubble.
    subroutine s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        real(wp) :: stddsv
        real(wp) :: strength_vel, strength_vol

        real(wp), dimension(3) :: nodecoord
        real(wp) :: addFun1, addFun2, addFun3
        real(wp) :: func, func2, volpart
        integer, dimension(3) :: cellaux
        real(wp), dimension(3) :: s_coord
        integer :: l, i, j, k
        logical :: celloutside
        integer :: smearGrid, smearGridz

        smearGrid = mapCells - (-mapCells) + 1 ! Include the cell that contains the bubble (3+1+3)
        smearGridz = smearGrid
        if (p == 0) smearGridz = 1

        $:GPU_PARALLEL_LOOP(private='[nodecoord,l,s_coord,cell,center]', copyin='[smearGrid,smearGridz]')
        do l = 1, nBubs
            nodecoord(1:3) = 0
            center(1:3) = 0._wp
            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            center(1:2) = lbk_pos(l, 1:2, 2)
            if (p > 0) center(3) = lbk_pos(l, 3, 2)
            call s_get_cell(s_coord, cell)
            call s_compute_stddsv(cell, volpart, stddsv)

            strength_vol = volpart
            strength_vel = 4._wp*pi*lbk_rad(l, 2)**2._wp*lbk_vel(l, 2)

            $:GPU_LOOP(collapse=3,private='[cellaux,nodecoord]')
            do i = 1, smearGrid
                do j = 1, smearGrid
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !Check if the cells intended to smear the bubbles in are in the computational domain
                        !and redefine the cells for symmetric boundary
                        call s_check_celloutside(cellaux, celloutside)

                        if (.not. celloutside) then

                            nodecoord(1) = x_cc(cellaux(1))
                            nodecoord(2) = y_cc(cellaux(2))
                            if (p > 0) nodecoord(3) = z_cc(cellaux(3))
                            call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)
                            if (lag_params%cluster_type >= 4) call s_applygaussian(center, cellaux, nodecoord, stddsv, 1._wp, func2)

                            ! Relocate cells for bubbles intersecting symmetric boundaries
                            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                                call s_shift_cell_symmetric_bc(cellaux, cell)
                            end if
                        else
                            func = 0._wp
                            func2 = 0._wp
                            cellaux(1) = cell(1)
                            cellaux(2) = cell(2)
                            cellaux(3) = cell(3)
                            if (p == 0) cellaux(3) = 0
                        end if

                        !Update void fraction field
                        addFun1 = func*strength_vol
                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                            updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                            + real(addFun1, kind=stp)

                        !Update time derivative of void fraction
                        addFun2 = func*strength_vel
                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                            updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                            + real(addFun2, kind=stp)

                        !Product of two smeared functions
                        !Update void fraction * time derivative of void fraction
                        if (lag_params%cluster_type >= 4) then
                            addFun3 = func2*strength_vol*strength_vel
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun3, kind=stp)
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_gaussian

    !> The purpose of this subroutine is to apply the gaussian kernel function for each bubble (Maeda and Colonius, 2018)).
    subroutine s_applygaussian(center, cellaux, nodecoord, stddsv, strength_idx, func)
        $:GPU_ROUTINE(function_name='s_applygaussian',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: center
        integer, dimension(3), intent(in) :: cellaux
        real(wp), dimension(3), intent(in) :: nodecoord
        real(wp), intent(in) :: stddsv
        real(wp), intent(in) :: strength_idx
        real(wp), intent(out) :: func

        real(wp) :: distance
        real(wp) :: theta, dtheta, L2, dzp, Lz2
        real(wp) :: Nr, Nr_count

        distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + (center(3) - nodecoord(3))**2._wp)

        if (num_dims == 3) then
            !< 3D gaussian function
            func = exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
        else
            if (cyl_coord) then
                !< 2D cylindrical function:
                ! We smear particles in the azimuthal direction for given r
                theta = 0._wp
                Nr = ceiling(2._wp*pi*nodecoord(2)/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                dtheta = 2._wp*pi/Nr
                L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                ! Factor 2._wp is for symmetry (upper half of the 2D field (+r) is considered)
                func = dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                Nr_count = 0._wp
                do while (Nr_count < Nr - 1._wp)
                    Nr_count = Nr_count + 1._wp
                    theta = Nr_count*dtheta
                    ! trigonometric relation
                    L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                    ! nodecoord(2)*dtheta is the azimuthal width of the cell
                    func = func + &
                           dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**(3._wp*(strength_idx + 1._wp))
                end do
            else

                !< 2D cartesian function:
                ! We smear particles considering a virtual depth (lag_params%charwidth)
                theta = 0._wp
                Nr = ceiling(lag_params%charwidth/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                Nr_count = 1._wp - mapCells*1._wp
                dzp = y_cb(cellaux(2) + 1) - y_cb(cellaux(2))
                Lz2 = (center(3) - (dzp*(0.5_wp + Nr_count) - lag_params%charwidth/2._wp))**2._wp
                distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                func = dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                do while (Nr_count < Nr - 1._wp + ((mapCells - 1)*1._wp))
                    Nr_count = Nr_count + 1._wp
                    Lz2 = (center(3) - (dzp*(0.5_wp + Nr_count) - lag_params%charwidth/2._wp))**2._wp
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                    func = func + &
                           dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**(3._wp*(strength_idx + 1._wp))
                end do
            end if
        end if

    end subroutine s_applygaussian

    !> The purpose of this subroutine is to check if the current cell is outside the computational domain or not (including ghost cells).
            !! @param cellaux Tested cell to smear the bubble effect in.
            !! @param celloutside If true, then cellaux is outside the computational domain.
    subroutine s_check_celloutside(cellaux, celloutside)
        $:GPU_ROUTINE(function_name='s_check_celloutside',parallelism='[seq]', &
            & cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        logical, intent(out) :: celloutside

        celloutside = .false.

        if (num_dims == 2) then
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                celloutside = .true.
            end if
            if (cyl_coord .and. y_cc(cellaux(2)) < 0._wp) then
                celloutside = .true.
            end if
            if ((cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                celloutside = .true.
            end if
        else
            if ((cellaux(3) < -buff_size) .or. (cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(3) > p + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(1) > m + buff_size)) then
                celloutside = .true.
            end if
        end if

    end subroutine s_check_celloutside

    !> This subroutine relocates the current cell, if it intersects a symmetric boundary.
            !! @param cell Cell of the current bubble
            !! @param cellaux Cell to map the bubble effect in.
    subroutine s_shift_cell_symmetric_bc(cellaux, cell)
        $:GPU_ROUTINE(function_name='s_shift_cell_symmetric_bc', &
            & parallelism='[seq]', cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        integer, dimension(3), intent(in) :: cell

        ! x-dir
        if (bc_x%beg == BC_REFLECTIVE .and. (cell(1) <= mapCells - 1)) then
            cellaux(1) = abs(cellaux(1)) - 1
        end if
        if (bc_x%end == BC_REFLECTIVE .and. (cell(1) >= m + 1 - mapCells)) then
            cellaux(1) = cellaux(1) - (2*(cellaux(1) - m) - 1)
        end if

        !y-dir
        if (bc_y%beg == BC_REFLECTIVE .and. (cell(2) <= mapCells - 1)) then
            cellaux(2) = abs(cellaux(2)) - 1
        end if
        if (bc_y%end == BC_REFLECTIVE .and. (cell(2) >= n + 1 - mapCells)) then
            cellaux(2) = cellaux(2) - (2*(cellaux(2) - n) - 1)
        end if

        if (p > 0) then
            !z-dir
            if (bc_z%beg == BC_REFLECTIVE .and. (cell(3) <= mapCells - 1)) then
                cellaux(3) = abs(cellaux(3)) - 1
            end if
            if (bc_z%end == BC_REFLECTIVE .and. (cell(3) >= p + 1 - mapCells)) then
                cellaux(3) = cellaux(3) - (2*(cellaux(3) - p) - 1)
            end if
        end if

    end subroutine s_shift_cell_symmetric_bc

    !> Calculates the standard deviation of the bubble being smeared in the Eulerian framework.
            !! @param cell Cell where the bubble is located
            !! @param volpart Volume of the bubble
            !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, volpart, stddsv)
        $:GPU_ROUTINE(function_name='s_compute_stddsv',parallelism='[seq]', &
            & cray_inline=True)

        integer, dimension(3), intent(in) :: cell
        real(wp), intent(in) :: volpart
        real(wp), intent(out) :: stddsv

        real(wp) :: chardist, charvol
        real(wp) :: rad

        !< Compute characteristic distance
        chardist = sqrt(dx(cell(1))*dy(cell(2)))
        if (p > 0) chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1._wp/3._wp)

        !< Compute characteristic volume
        if (p > 0) then
            charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
        else
            if (cyl_coord) then
                charvol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                charvol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
            end if
        end if

        !< Compute Standard deviaton
        if (((volpart/charvol) > 0.5_wp*lag_params%valmaxvoid) .or. (lag_params%smooth_type == 1)) then
            rad = (3._wp*volpart/(4._wp*pi))**(1._wp/3._wp)
            stddsv = 1._wp*lag_params%epsilonb*max(chardist, rad)
        else
            stddsv = 0._wp
        end if

    end subroutine s_compute_stddsv

    !> The purpose of this procedure is to calculate the characteristic cell volume
            !! @param cell Computational coordinates (x, y, z)
            !! @param Charvol Characteristic volume
    subroutine s_get_char_vol(cellx, celly, cellz, Charvol)
        $:GPU_ROUTINE(function_name='s_get_char_vol',parallelism='[seq]', &
            & cray_inline=True)

        integer, intent(in) :: cellx, celly, cellz
        real(wp), intent(out) :: Charvol

        if (p > 0) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*lag_params%charwidth
            end if
        end if

    end subroutine s_get_char_vol

    !> This subroutine transforms the computational coordinates of the bubble from
            !!      real type into integer.
            !! @param s Computational coordinates of the bubble, real type
            !! @param get_cell Computational coordinates of the bubble, integer type
    subroutine s_get_cell(s_cell, get_cell)
        $:GPU_ROUTINE(function_name='s_get_cell',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: s_cell
        integer, dimension(3), intent(out) :: get_cell
        integer :: i

        get_cell(:) = int(s_cell(:))
        do i = 1, num_dims
            if (s_cell(i) < 0._wp) get_cell(i) = get_cell(i) - 1
        end do

    end subroutine s_get_cell

end module m_bubbles_EL_kernels
