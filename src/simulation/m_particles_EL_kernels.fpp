!>
!! @file m_particles_EL_kernels.f90
!! @brief Contains module m_particles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian particles
!!        in the Eulerian framework.
module m_particles_EL_kernels

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use ieee_arithmetic        !< For checking NaN

    implicit none

contains

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      particles into the Eulerian framework using different approaches.
            !! @param nParticles Number of lagrangian particles in the current domain
            !! @param lbk_rad Radius of the particles
            !! @param lbk_s Computational coordinates of the particles
            !! @param lbk_pos Spatial coordinates of the particles
            !! @param updatedvar Eulerian variable to be updated
            !! @param lbk_f_p Forces on the particles
    subroutine s_smoothfunction(nParticles, lbk_rad, lbk_s, lbk_pos, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p

        call s_trilinear_projection(nParticles, lbk_rad, lbk_s, lbk_pos, updatedvar, lbk_f_p)

    end subroutine s_smoothfunction

    !> The purpose of this procedure is to use trilinear projection to map the effect of the particles to the fluid.
    subroutine s_trilinear_projection(nParticles, lbk_rad, lbk_s, lbk_pos, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        real(wp) :: addFun1, addFun2_x, addFun2_y, addFun2_z, addFun_E
        real(wp) :: volpart, Vol
        real(wp) :: fp_x, fp_y, fp_z
        real(wp) :: vol_frac
        real(wp), dimension(3) :: s_coord
        real(wp) :: wx, wy, wz, w
        integer :: l, i, j, k, nx, ny, nz, iL, jL, kL, kL_max, cx, cy, cz

        $:GPU_PARALLEL_LOOP(private='[l,s_coord,cell,center,fp_x,fp_y,fp_z,i,j,k]')
        do l = 1, nParticles

            fp_x = -lbk_f_p(l, 1)
            fp_y = -lbk_f_p(l, 2)
            fp_z = -lbk_f_p(l, 3)

            center(1:3) = 0._wp
            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            center(1:2) = lbk_pos(l, 1:2, 2)

            if (p > 0) center(3) = lbk_pos(l, 3, 2)
            call s_get_cell(s_coord, cell)

            if (num_dims == 2) then
                Vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            end if
            vol_frac = volpart/Vol

            i = cell(1)
            j = cell(2)
            k = cell(3)

            if (center(1) >= x_cc(i)) then
                iL = i
            else
                iL = i - 1
            end if
            if (center(2) >= y_cc(j)) then
                jL = j
            else
                jL = j - 1
            end if
            kL = 0
            kL_max = -1
            if (p > 0) then
                if (center(3) >= z_cc(k)) then
                    kL = k
                    kL_max = kL
                else
                    kL = k - 1
                    kL_max = kL
                end if
            end if
            wx = (center(1) - x_cc(iL))/(x_cc(iL + 1) - x_cc(iL))
            wy = (center(2) - y_cc(jL))/(y_cc(jL + 1) - y_cc(jL))
            wz = 0._wp
            if (p > 0) wz = (center(3) - z_cc(kL))/(z_cc(kL + 1) - z_cc(kL))

            ! $:GPU_ATOMIC(atomic='update')
            ! updatedvar(1)%sf(i, j, k) = updatedvar(1)%sf(i, j, k) &
            !                             + real(vol_frac, kind=stp)

            $:GPU_LOOP(collapse=3,private='[w]')
            do nx = iL, iL + 1
                do ny = jL, jL + 1
                    do nz = kL, kL_max + 1
                        cx = nx - iL
                        cy = ny - jL
                        cz = nz - kL
                        ! Compute weight for this corner
                        w = ((1.0_wp - wx)*(1 - cx) + wx*cx)* &
                            ((1.0_wp - wy)*(1 - cy) + wy*cy)* &
                            ((1.0_wp - wz)*(1 - cz) + wz*cz)

                        addFun1 = (w*vol_frac)
                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(1)%sf(nx, ny, nz) = updatedvar(1)%sf(nx, ny, nz) &
                                                       + real(addFun1, kind=stp)

                        if (lag_params%solver_approach == 2) then
                            ! Add particle force to the grid cell
                            addFun2_x = (w*fp_x)/Vol
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(2)%sf(nx, ny, nz) = &
                                updatedvar(2)%sf(nx, ny, nz) &
                                + real(addFun2_x, kind=stp)

                            addFun2_y = (w*fp_y)/Vol
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(3)%sf(nx, ny, nz) = &
                                updatedvar(3)%sf(nx, ny, nz) &
                                + real(addFun2_y, kind=stp)

                            if (num_dims == 3) then
                                addFun2_z = (w*fp_z)/Vol
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(4)%sf(nx, ny, nz) = &
                                    updatedvar(4)%sf(nx, ny, nz) &
                                    + real(addFun2_z, kind=stp)

                                addFun_E = 0._wp*w
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(5)%sf(nx, ny, nz) = &
                                    updatedvar(5)%sf(nx, ny, nz) &
                                    + real(addFun_E, kind=stp)
                            else
                                addFun_E = 0._wp*w
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(4)%sf(nx, ny, nz) = &
                                    updatedvar(4)%sf(nx, ny, nz) &
                                    + real(addFun_E, kind=stp)
                            end if
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_trilinear_projection

    !> The purpose of this procedure contains the algorithm to use the delta kernel function to map the effect of the particles.
            !!      The effect of the particles only affects the cell where the particle is located.
    subroutine s_deltafunc(nParticles, lbk_rad, lbk_s, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        integer, dimension(3) :: cell
        real(wp) :: strength_vol

        real(wp) :: addFun1, addFun2_x, addFun2_y, addFun2_z, addFun_E
        real(wp) :: volpart, Vol
        real(wp), dimension(3) :: s_coord
        integer :: l
        real(wp) :: fp_x, fp_y, fp_z

        $:GPU_PARALLEL_LOOP(private='[l,s_coord,cell]')
        do l = 1, nParticles

            fp_x = lbk_f_p(l, 1)
            fp_y = lbk_f_p(l, 2)
            fp_z = lbk_f_p(l, 3)

            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            call s_get_cell(s_coord, cell)

            strength_vol = volpart

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

            if (lag_params%solver_approach == 2) then

                !Update x-momentum source term
                addFun2_x = -fp_x/Vol
                $:GPU_ATOMIC(atomic='update')
                updatedvar(2)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar(2)%sf(cell(1), cell(2), cell(3)) &
                    + real(addFun2_x, kind=stp)

                !Update y-momentum source term
                addFun2_y = -fp_y/Vol
                $:GPU_ATOMIC(atomic='update')
                updatedvar(3)%sf(cell(1), cell(2), cell(3)) = &
                    updatedvar(3)%sf(cell(1), cell(2), cell(3)) &
                    + real(addFun2_y, kind=stp)

                if (num_dims == 3) then
                    !Update z-momentum source term
                    addFun2_z = -fp_z/Vol
                    $:GPU_ATOMIC(atomic='update')
                    updatedvar(4)%sf(cell(1), cell(2), cell(3)) = &
                        updatedvar(4)%sf(cell(1), cell(2), cell(3)) &
                        + real(addFun2_z, kind=stp)
                    !Update energy source term
                    addFun_E = 0._wp
                    $:GPU_ATOMIC(atomic='update')
                    updatedvar(5)%sf(cell(1), cell(2), cell(3)) = &
                        updatedvar(5)%sf(cell(1), cell(2), cell(3)) &
                        + real(addFun_E, kind=stp)
                else
                    !Update energy source term
                    addFun_E = 0._wp
                    $:GPU_ATOMIC(atomic='update')
                    updatedvar(4)%sf(cell(1), cell(2), cell(3)) = &
                        updatedvar(4)%sf(cell(1), cell(2), cell(3)) &
                        + real(addFun_E, kind=stp)
                end if
            end if

        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_deltafunc

    !> The purpose of this procedure contains the algorithm to use the gaussian kernel function to map the effect of the particles.
            !!      The effect of the particles affects the 3X3x3 cells that surround the particle.
    subroutine s_gaussian(nParticles, lbk_rad, lbk_s, lbk_pos, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        real(wp) :: stddsv
        real(wp) :: strength_vol

        real(wp), dimension(3) :: nodecoord
        real(wp) :: addFun1, addFun2_x, addFun2_y, addFun2_z, addFun_E, func_sum
        real(wp) :: func, volpart, Vol
        integer, dimension(3) :: cellaux
        real(wp), dimension(3) :: s_coord
        integer :: l, i, j, k
        logical :: celloutside
        integer :: smearGrid, smearGridz
        real(wp) :: fp_x, fp_y, fp_z

        smearGrid = mapCells - (-mapCells) + 1 ! Include the cell that contains the particle (3+1+3)
        smearGridz = smearGrid
        if (p == 0) smearGridz = 1

        $:GPU_PARALLEL_LOOP(private='[cellaux,nodecoord,l,s_coord,cell,center]', copyin='[smearGrid,smearGridz]')
        do l = 1, nParticles
            nodecoord(1:3) = 0
            center(1:3) = 0._wp
            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            center(1:2) = lbk_pos(l, 1:2, 2)

            if (p > 0) center(3) = lbk_pos(l, 3, 2)
            call s_get_cell(s_coord, cell)
            call s_compute_stddsv(cell, volpart, stddsv)
            strength_vol = volpart

            if (num_dims == 2) then
                Vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            end if

            fp_x = lbk_f_p(l, 1)
            fp_y = lbk_f_p(l, 2)
            fp_z = lbk_f_p(l, 3)

            func_sum = 0._wp
            do i = 1, smearGrid
                do j = 1, smearGrid
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !Check if the cells intended to smear the particles in are in the computational domain
                        !and redefine the cells for symmetric boundary
                        call s_check_celloutside(cellaux, celloutside)

                        if (.not. celloutside) then
                            nodecoord(1) = x_cc(cellaux(1))
                            nodecoord(2) = y_cc(cellaux(2))
                            if (p > 0) nodecoord(3) = z_cc(cellaux(3))
                            call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)
                            func_sum = func_sum + func
                        end if
                    end do
                end do
            end do

            $:GPU_LOOP(collapse=3,private='[cellaux,nodecoord]')
            do i = 1, smearGrid
                do j = 1, smearGrid
                    do k = 1, smearGridz
                        cellaux(1) = cell(1) + i - (mapCells + 1)
                        cellaux(2) = cell(2) + j - (mapCells + 1)
                        cellaux(3) = cell(3) + k - (mapCells + 1)
                        if (p == 0) cellaux(3) = 0

                        !Check if the cells intended to smear the particles in are in the computational domain
                        !and redefine the cells for symmetric boundary
                        call s_check_celloutside(cellaux, celloutside)

                        if (.not. celloutside) then
                            nodecoord(1) = x_cc(cellaux(1))
                            nodecoord(2) = y_cc(cellaux(2))
                            if (p > 0) nodecoord(3) = z_cc(cellaux(3))
                            call s_applygaussian(center, cellaux, nodecoord, stddsv, 0._wp, func)
                            func = func/func_sum

                            ! Relocate cells for particles intersecting symmetric boundaries
                            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                                call s_shift_cell_symmetric_bc(cellaux, cell)
                            end if

                            !Update void fraction field
                            addFun1 = func*strength_vol !strength_vol/Vol
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun1, kind=stp)

                            if (lag_params%solver_approach == 2) then

                                !Update x-momentum source term
                                addFun2_x = -func*fp_x
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + real(addFun2_x, kind=stp)

                                !Update y-momentum source term
                                addFun2_y = -func*fp_y
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(3)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar(3)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + real(addFun2_y, kind=stp)

                                if (num_dims == 3) then
                                    !Update z-momentum source term
                                    addFun2_z = -func*fp_z
                                    $:GPU_ATOMIC(atomic='update')
                                    updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                        updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        + real(addFun2_z, kind=stp)
                                    !Update energy source term
                                    addFun_E = 0._wp
                                    $:GPU_ATOMIC(atomic='update')
                                    updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                        updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        + real(addFun_E, kind=stp)
                                else
                                    !Update energy source term
                                    addFun_E = 0._wp
                                    $:GPU_ATOMIC(atomic='update')
                                    updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                        updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                        + real(addFun_E, kind=stp)
                                end if
                            end if
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
                dzp = y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)
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
            !! @param cellaux Tested cell to smear the particle effect in.
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

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size)) then
                celloutside = .true.
            end if
        else
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size) .or. (cellaux(3) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(3) > p + buff_size)) then
                celloutside = .true.
            end if
        end if

    end subroutine s_check_celloutside

    !> This subroutine relocates the current cell, if it intersects a symmetric boundary.
            !! @param cell Cell of the current particle
            !! @param cellaux Cell to map the particle effect in.
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

    !> Calculates the standard deviation of the particle being smeared in the Eulerian framework.
            !! @param cell Cell where the particle is located
            !! @param volpart Volume of the particle
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

    !> This subroutine transforms the computational coordinates of the particle from
            !!      real type into integer.
            !! @param s Computational coordinates of the particle, real type
            !! @param get_cell Computational coordinates of the particle, integer type
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

    function f_interp_prim_trilinear(pos, cell, field_vf, field_index) result(val)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        type(scalar_field), dimension(:), intent(in) :: field_vf
        integer, intent(in) :: field_index

        real(wp) :: val
        real(wp) :: wx, wy, wz, x_lo, y_lo, z_lo
        real(wp) :: fx0, fx1, fy0, fy1
        integer :: i, j, k, iL, jL, kL

        i = cell(1)
        j = cell(2)
        k = cell(3)

        if (pos(1) >= x_cc(i)) then
            iL = i
        else
            iL = i - 1
        end if

        if (pos(2) >= y_cc(j)) then
            jL = j
        else
            jL = j - 1
        end if

        wx = (pos(1) - x_cc(iL))/(x_cc(iL + 1) - x_cc(iL))
        wy = (pos(2) - y_cc(jL))/(y_cc(jL + 1) - y_cc(jL))

        ! 2D case bilinear
        if (num_dims == 2) then

            fx0 = (1 - wx)*field_vf(field_index)%sf(iL, jL, 0) + &
                  wx*field_vf(field_index)%sf(iL + 1, jL, 0)

            fx1 = (1 - wx)*field_vf(field_index)%sf(iL, jL + 1, 0) + &
                  wx*field_vf(field_index)%sf(iL + 1, jL + 1, 0)

            val = (1 - wy)*fx0 + wy*fx1
            return
        end if

        ! 3D case trilinear
        if (pos(3) >= z_cc(k)) then
            kL = k
        else
            kL = k - 1
        end if
        wz = (pos(3) - z_cc(kL))/(z_cc(kL + 1) - z_cc(kL))

        fx0 = (1 - wx)*field_vf(field_index)%sf(iL, jL, kL) + wx*field_vf(field_index)%sf(iL + 1, jL, kL)
        fx1 = (1 - wx)*field_vf(field_index)%sf(iL, jL + 1, kL) + wx*field_vf(field_index)%sf(iL + 1, jL + 1, kL)
        fy0 = (1 - wy)*fx0 + wy*fx1

        fx0 = (1 - wx)*field_vf(field_index)%sf(iL, jL, kL + 1) + wx*field_vf(field_index)%sf(iL + 1, jL, kL + 1)
        fx1 = (1 - wx)*field_vf(field_index)%sf(iL, jL + 1, kL + 1) + wx*field_vf(field_index)%sf(iL + 1, jL + 1, kL + 1)
        fy1 = (1 - wy)*fx0 + wy*fx1

        val = (1 - wz)*fy0 + wz*fy1

    end function

    !! This function calculates the force on a particle
            !!      based on the pressure gradient, velocity, and drag model.
            !! @param pos Position of the particle
            !! @param rad Radius of the particle
            !! @param vel_p Velocity of the particle
            !! @param mass_p Particle mass
            !! @param Re Viscosity!
            !! @param rho Density of the fluid
            !! @param vol_frac Particle Volume Fraction
            !! @param cell Computational coordinates of the particle
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return a Acceleration of the particle in direction i
    subroutine f_get_particle_force(pos, rad, vel_p, mass_p, Re, gamm, vol_frac, cell, &
                                    q_prim_vf, fieldvars, force, rmass_add)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rad, mass_p, Re, gamm, vol_frac
        real(wp), dimension(3), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        real(wp), dimension(3), intent(in) :: vel_p
        ! real(stp), intent(in) :: pressure_gradient(:,:,:,:)
        ! real(stp), intent(in) :: density_gradient(:,:,:,:)
        type(scalar_field), dimension(:), intent(in) :: fieldvars
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(3), intent(out) :: force
        real(wp), intent(out) :: rmass_add

        real(wp) :: a, vol, rho_fluid, pressure_fluid
        real(wp), dimension(3) :: v_rel, dp
        real(wp), dimension(fd_order) :: xi, eta, L
        real(wp) :: particle_diam, gas_mu, vmag, cson
        real(wp) :: slip_velocity_x, slip_velocity_y, slip_velocity_z, beta
        real(wp), dimension(3) :: fluid_vel
        integer :: dir

        !Added pass params
        real(wp) :: mach, Cam, udot_grad_rho, flux_f, flux_b, div_rhou, SDrho, vgradrho, drhodt
        real(wp), dimension(3) :: rhoDuDt, grad_rho, fam
        integer, dimension(3) :: p1

        force = 0._wp
        dp = 0._wp
        grad_rho = 0._wp
        drhodt = 0._wp
        fam = 0._wp
        fluid_vel = 0._wp
        v_rel = 0._wp
        rhoDuDt = 0._wp

        !!Interpolation - either tri-linear or 0th order
        if (fd_order > 1) then
            rho_fluid = f_interp_prim_trilinear(pos, cell, q_prim_vf, 1)
            pressure_fluid = f_interp_prim_trilinear(pos, cell, q_prim_vf, E_idx)
            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = f_interp_prim_trilinear(pos, cell, fieldvars, dir)
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = f_interp_prim_trilinear(pos, cell, fieldvars, 3 + dir)
                    drhodt = drhodt + f_interp_prim_trilinear(pos, cell, fieldvars, 6 + dir)
                end if
                fluid_vel(dir) = f_interp_prim_trilinear(pos, cell, q_prim_vf, momxb + dir - 1)
            end do
        else
            rho_fluid = q_prim_vf(1)%sf(cell(1), cell(2), cell(3))
            pressure_fluid = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))
            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = fieldvars(dir)%sf(cell(1), cell(2), cell(3))
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = fieldvars(3 + dir)%sf(cell(1), cell(2), cell(3))
                    drhodt = drhodt + fieldvars(6 + dir)%sf(cell(1), cell(2), cell(3))
                end if
                fluid_vel(dir) = q_prim_vf(momxb + dir - 1)%sf(cell(1), cell(2), cell(3))
            end do
        end if

        drhodt = -drhodt

        v_rel = vel_p - fluid_vel

        if (lag_params%qs_drag_model > 0 .or. lag_params%added_mass_model > 0) then
            ! Quasi-steady Drag Force Parameters
            slip_velocity_x = fluid_vel(1) - vel_p(1)
            slip_velocity_y = fluid_vel(2) - vel_p(2)
            if (num_dims == 3) then
                slip_velocity_z = fluid_vel(3) - vel_p(3)
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y + &
                            slip_velocity_z*slip_velocity_z)
            elseif (num_dims == 2) then
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y)
            end if
            particle_diam = rad*2._wp
            cson = sqrt((gamm*pressure_fluid)/rho_fluid) !gamma*P/rho
            gas_mu = Re
        end if

        if (lag_params%added_mass_model > 0) then
            rhoDuDt = -dp
            udot_grad_rho = dot_product(fluid_vel, grad_rho)
            vgradrho = dot_product(vel_p, grad_rho)
            SDrho = (drhodt + udot_grad_rho)/(1._wp - vol_frac)
            mach = vmag/cson
        end if

        ! Step 1: Force component quasi-steady
        if (lag_params%qs_drag_model == 1) then
            beta = QS_Parmar(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 2) then
            beta = QS_Osnes(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 3) then
            beta = QS_ModifiedParmar(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 4) then
            beta = QS_Gidaspow(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else
            !No Quasi-Steady drag
        end if

        ! Step 1.1: Stokes drag
        if (lag_params%stokes_drag == 1) then ! Free slip Stokes drag
            force = force - 4._wp*pi*gas_mu*rad*v_rel
        elseif (lag_params%stokes_drag == 2) then ! No slip Stokes drag
            force = force - 6._wp*pi*gas_mu*rad*v_rel
        else
            !No stokes drag
        end if

        ! Step 2: Pressure Gradient Force
        if (lag_params%pressure_force) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            force = force - vol*dp
        end if

        ! Step 3: Gravitational Force
        if (lag_params%gravity_force) then
            force = force + (mass_p)*accel_bf
        end if

        ! Step 4: Added Mass Force
        if (lag_params%added_mass_model == 1) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            if (mach > 0.6_wp) then
                Cam = 1._wp + 1.8_wp*(0.6_wp**2) + 7.6_wp* &
                      (0.6_wp**4)
            else
                Cam = 1._wp + 1.8_wp*mach**2 + 7.6_wp*mach**4
            end if

            Cam = 0.5_wp*Cam*(1._wp + 0.68_wp*vol_frac**2)
            rmass_add = rho_fluid*vol*Cam

            fam = Cam*vol*(vel_p*SDrho + rhoDuDt + &
                           fluid_vel*(vgradrho))

            do dir = 1, num_dims
                if (.not. ieee_is_finite(fam(dir))) then
                    fam(dir) = 0._wp
                    rmass_add = 0._wp
                end if
            end do

            force = force + fam
        else
            rmass_add = 0._wp
        end if

        do dir = 1, num_dims
            if (.not. ieee_is_finite(force(dir))) then
                force(dir) = 0._wp
            end if
        end do

    end subroutine f_get_particle_force

    ! Quasi-steady force (Re_p and Ma_p corrections):
    !   Improved Drag Correlation for Spheres and Application
    !   to Shock-Tube Experiments
    !   - Parmar et al. (2010)
    !   - AIAA Journal
    !
    ! Quasi-steady force (phi corrections):
    !   The Added Mass, Basset, and Viscous Drag Coefficients
    !   in Nondilute Bubbly Liquids Undergoing Small-Amplitude
    !   Oscillatory Motion
    !   - Sangani et al. (1991)
    !   - Phys. Fluids A
    function QS_Parmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp) :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp) :: beta
        real(wp) :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1.e-14_wp) then
            rcd1 = 1.0_wp
        else
            rmacr = 0.6_wp ! Critical rmachp no
            rcd_mcr = (1.+0.15*re**(0.684)) + &
                      (re/24.0)*(0.513/(1.+483./re**(0.669)))
            if (mp <= rmacr) then
                rcd_std = (1.+0.15*re**(0.687)) + &
                          (re/24.0)*(0.42/(1.+42500./re**(1.16)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                C1 = 6.48_wp
                C2 = 9.28_wp
                C3 = 12.21_wp
                f1M = -1.884_wp + 8.422_wp*mp - 13.70_wp*mp**2 + 8.162_wp*mp**3
                f2M = -2.228_wp + 10.35_wp*mp - 16.96_wp*mp**2 + 9.840_wp*mp**3
                f3M = 4.362_wp - 16.91_wp*mp + 19.84_wp*mp**2 - 6.296_wp*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                rcd_M2 = (1.0 + 0.107*re**0.867) + &
                         (re/24.0)*0.646/(1.0 + 861.0/re**.634)
                C1 = 6.48_wp
                C2 = 8.93_wp
                C3 = 12.21_wp
                f1M = -2.963 + 4.392*mp - 1.169*mp**2 - 0.027*mp**3 &
                      - 0.233*exp((1.0 - mp)/0.011)
                f2M = -6.617 + 12.11*mp - 6.501*mp**2 + 1.182*mp**3 &
                      - 0.174*exp((1.0 - mp)/0.010)
                f3M = -5.866 + 11.57*mp - 6.665*mp**2 + 1.312*mp**3 &
                      - 0.350*exp((1.0 - mp)/0.012)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0 + 0.107*re**0.867) + &
                       (re/24.0)*0.646/(1.0 + 861.0/re**.634)
            end if ! mp
        end if    ! re

        ! Sangani's volume fraction correction for dilute random arrays
        ! Capping volume fraction at 0.5
        phi_corr = (1.0 + 5.94*min(phi, 0.5))

        cd = (24.0/re)*rcd1*phi_corr

        beta = rcd1*3.0*pi*mu_fluid*dp

        beta = beta*phi_corr

    end function QS_Parmar

    ! Quasi-steady force (Re_p and Ma_p corrections):
    !   Improved Drag Correlation for Spheres and Application
    !   to Shock-Tube Experiments
    !   - Parmar et al. (2010)
    !   - AIAA Journal
    !
    ! Quasi-steady force (phi corrections):
    !   Sangani et al. (1991) volume fraction correction overshoots
    !   the drag coefficient.
    !
    !   We adopt instead Osnes et al. (2023) volume fraction correction
    !   based on Tenneti et al. with one extra term.
    !
    !   At Mach=0, the drag coefficient from this subroutine matches very
    !   well with the one calculated using the Osnes subroutine, for various
    !   Reynolds numbers and volume fractions.
    function QS_ModifiedParmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp) :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp) :: b1, b2, b3
        real(wp) :: beta
        real(wp) :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1e-14) then
            rcd1 = 1.0
        else
            rmacr = 0.6 ! Critical rmachp no.
            rcd_mcr = (1.+0.15*re**(0.684)) + &
                      (re/24.0)*(0.513/(1.+483./re**(0.669)))
            if (mp <= rmacr) then
                rcd_std = (1.+0.15*re**(0.687)) + &
                          (re/24.0)*(0.42/(1.+42500./re**(1.16)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                C1 = 6.48
                C2 = 9.28
                C3 = 12.21
                f1M = -1.884 + 8.422*mp - 13.70*mp**2 + 8.162*mp**3
                f2M = -2.228 + 10.35*mp - 16.96*mp**2 + 9.840*mp**3
                f3M = 4.362 - 16.91*mp + 19.84*mp**2 - 6.296*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                rcd_M2 = (1.0 + 0.107*re**0.867) + &
                         (re/24.0)*0.646/(1.0 + 861.0/re**.634)
                C1 = 6.48
                C2 = 8.93
                C3 = 12.21
                f1M = -2.963 + 4.392*mp - 1.169*mp**2 - 0.027*mp**3 &
                      - 0.233*exp((1.0 - mp)/0.011)
                f2M = -6.617 + 12.11*mp - 6.501*mp**2 + 1.182*mp**3 &
                      - 0.174*exp((1.0 - mp)/0.010)
                f3M = -5.866 + 11.57*mp - 6.665*mp**2 + 1.312*mp**3 &
                      - 0.350*exp((1.0 - mp)/0.012)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0 + 0.107*re**0.867) + &
                       (re/24.0)*0.646/(1.0 + 861.0/re**.634)
            end if ! mp
        end if    ! re

        ! Osnes's volume fraction correction
        b1 = 5.81*phi/((1.0 - phi)**2) + &
             0.48*(phi**(1._wp/3._wp))/((1.0 - phi)**3)

        b2 = ((1.0 - phi)**2)*(phi**3)* &
             re*(0.95 + 0.61*(phi**3)/((1.0 - phi)*2))

        b3 = min(sqrt(20.0_wp*mp), 1.0_wp)* &
             (5.65*phi - 22.0*(phi**2) + 23.4*(phi**3))* &
             (1 + tanh((mp - (0.65 - 0.24*phi))/0.35))

        cd = (24.0/re)*rcd1

        cd = cd/(1.0 - phi) + b3 + (24.0/re)*(1.0 - phi)*(b1 + b2)

        beta = 3.0*pi*mu_fluid*dp*(re/24.0)*cd

    end function QS_ModifiedParmar

    ! QS Force calculated as a function of Re, Ma and phi
    !
    ! Use Osnes etal (2023) correlations
    ! A.N. Osnes, M. Vartdal, M. Khalloufi,
    !    J. Capecelatro, and S. Balachandar.
    ! Comprehensive quasi-steady force correlations for compressible flow
    !    through random particle suspensions.
    ! International Journal of Multiphase Flow, Vol. 165, 104485, (2023).
    ! doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
    !
    ! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura.
    ! Supersonic and hypersonic drag coefficients for a sphere.
    ! AIAA Journal, Vol. 59(8), pp. 3261-3274, (2021).
    ! doi: https://doi.org/10.2514/1.J060153.
    !
    ! NOTE: Re<45 Rarefied fomu_fluidla of Loth et al has been redefined by Balachandar
    ! to avoid singularity as Ma -> 0.
    function QS_Osnes(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rmachp, mp, phi, rep, re
        real(wp) :: Knp, fKn, CD1, s, JM, CD2, cd_loth, CM, GM, HM, b1, b2, b3, cd, sgby2, JMt
        real(wp) :: beta

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        ! Loth's correlation
        if (re <= 45.0_wp) then
            ! Rarefied-dominated regime
            Knp = sqrt(0.5_wp*pi*gamma)*mp/re
            if (Knp > 0.01) then
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/Knp)))
            else
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/0.01)))

            end if
            CD1 = (24.0_wp/re)*(1.0_wp + 0.15_wp*re**(0.687_wp))*fKn
            s = mp*sqrt(0.5_wp*gamma)
            sgby2 = sqrt(0.5_wp*gamma)
            if (mp <= 1) then
                !JMt = 2.26_wp*(mp**4) - 0.1_wp*(mp**3) + 0.14_wp*mp
                JMt = 2.26_wp*(mp**4) + 0.14_wp*mp
            else
                JMt = 1.6_wp*(mp**4) + 0.25_wp*(mp**3) + 0.11_wp*(mp**2) + 0.44_wp*mp
            end if
            !
            ! Refomu_fluidlated version of Loth et al. to avoid singularity at mp = 0
            !
            CD2 = (1.0_wp + 2.0_wp*(s**2))*exp(-s**2)*mp &
                  /((sgby2**3)*sqrt(pi)) &
                  + (4.0_wp*(s**4) + 4.0_wp*(s**2) - 1.0_wp) &
                  *erf(s)/(2.0_wp*(sgby2**4)) &
                  + (2.0_wp*(mp**3)/(3.0_wp*sgby2))*sqrt(pi)

            CD2 = CD2/(1.0_wp + (((CD2/JMt) - 1.0_wp)*sqrt(re/45.0_wp)))
            cd_loth = CD1/(1.0_wp + (mp**4)) &
                      + CD2/(1.0_wp + (mp**4))
        else
            ! Compression-dominated regime
            ! TLJ: coefficients tweaked to get continuous values
            !      on the two branches at the critical points
            if (mp < 1.5_wp) then
                CM = 1.65_wp + 0.65_wp*tanh(4_wp*mp - 3.4_wp)
            else
                !CM = 2.18_wp - 0.13_wp*tanh(0.9_wp*mp - 2.7_wp)
                CM = 2.18_wp - 0.12913149918318745_wp*tanh(0.9_wp*mp - 2.7_wp)
            end if
            if (mp < 0.8) then
                GM = 166.0_wp*(mp**3) + 3.29_wp*(mp**2) - 10.9_wp*mp + 20._wp
            else
                !GM = 5.0_wp + 40._wp*(mp**(-3))
                GM = 5.0_wp + 47.809331200000017_wp*(mp**(-3))
            end if
            if (mp < 1) then
                HM = 0.0239_wp*(mp**3) + 0.212_wp*(mp**2) &
                     - 0.074_wp*mp + 1._wp
            else
                !HM =   0.93_wp + 1.0_wp / (3.5_wp + (mp**5))
                HM = 0.93967777777777772_wp + 1.0_wp/(3.5_wp + (mp**5))
            end if

            cd_loth = (24.0_wp/re)*(1 + 0.15_wp*(re**(0.687)))*HM + &
                      0.42_wp*CM/(1 + 42500/re**(1.16*CM) + GM/sqrt(re))

        end if

        b1 = 5.81_wp*phi/((1.0 - phi)**2) + &
             0.48_wp*(phi**(1._wp/3._wp))/((1.0 - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)* &
             re*(0.95 + 0.61*(phi**3)/((1.0 - phi)*2))

        b3 = min(sqrt(20.0_wp*mp), 1.0_wp)* &
             (5.65*phi - 22.0*(phi**2) + 23.4*(phi**3))* &
             (1 + tanh((mp - (0.65 - 0.24*phi))/0.35))

        cd = cd_loth/(1.0_wp - phi) + b3 + (24.0_wp/re)*(1.0_wp - phi)*(b1 + b2)

        beta = 3.0_wp*pi*mu_fluid*dp*(re/24.0_wp)*cd

    end function QS_Osnes

! Subroutine for Quasi-Steady Drag Model of Gidaspow
!
! D. Gidaspow, Multiphase Flow and Fluidization (Academic Press, 1994)
!
! Note: Model is provided per cell volume. We convert that to per
! particle using the particle volume fraction and volume
    function QS_Gidaspow(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: cd, phifRep, phif
        real(wp) :: phi, rep, re
        real(wp) :: beta

        rep = vmag*dp*rho/mu_fluid
        phi = max(volume_fraction, 0.0001_wp)
        phif = max(1._wp - volume_fraction, 0.0001_wp)
        re = max(rep, 0.1_wp)

        phifRep = phif*re

        if (phifRep < 1000.0) then
            cd = 24.0/phifRep*(1.0 + 0.15*(phifRep)**0.687)
        else
            cd = 0.44
        end if

        if (phif < 0.8) then
            beta = 150.0*((phi**2)*mu_fluid)/(phif*dp**2) &
                   + 1.75*(rho*phi*vmag/dp)
        else
            beta = 0.75*cd*phi*rho*vmag/(dp*phif**1.65)
        end if

        beta = beta*(pi*dp**3)/(6.0*phi)

    end function QS_Gidaspow

end module m_particles_EL_kernels
