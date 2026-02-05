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

        smoothfunc:select case(lag_params%smooth_type)
        case (1)
        call s_gaussian(nParticles, lbk_rad, lbk_s, lbk_pos, updatedvar, lbk_f_p)
        case (2)
        call s_deltafunc(nParticles, lbk_rad, lbk_s, updatedvar, lbk_f_p)
        end select smoothfunc

    end subroutine s_smoothfunction

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

            fp_x = lbk_f_p(l,1)
            fp_y = lbk_f_p(l,2)
            fp_z = lbk_f_p(l,3)

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
              endif
            endif        

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
        real(wp) :: func, volpart
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
            fp_x = lbk_f_p(l,1)
            fp_y = lbk_f_p(l,2)
            fp_z = lbk_f_p(l,3)

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
                    enddo
                enddo
            enddo

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
                            addFun1 = func*strength_vol
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
                              endif
                            endif
                        endif
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

    !! This function interpolates the velocity of Eulerian field at the position
            !! of the particle.
            !! @param pos Position of the particle in directiion i
            !! @param cell Computational coordinates of the particle
            !! @param i Direction of the velocity (1: x, 2: y, 3: z)
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return v Interpolated velocity at the position of the particle
    function f_interpolate_velocity(pos, cell, i, q_prim_vf) result(v)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        integer, intent(in) :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp) :: v
        real(wp), dimension(fd_order + 1) :: xi, eta, L

        if (fd_order == 1) then
            v = q_prim_vf(momxb+i-1)%sf(cell(1), cell(2), cell(3))
        
        elseif (fd_order == 2) then
            if (i == 1) then
                xi(1) = x_cc(cell(1) - 1)
                eta(1) = q_prim_vf(momxb)%sf(cell(1) - 1, cell(2), cell(3))
                xi(2) = x_cc(cell(1))
                eta(2) = q_prim_vf(momxb)%sf(cell(1), cell(2), cell(3))
                xi(3) = x_cc(cell(1) + 1)
                eta(3) = q_prim_vf(momxb)%sf(cell(1) + 1, cell(2), cell(3))
            elseif (i == 2) then
                xi(1) = y_cc(cell(2) - 1)
                eta(1) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 1, cell(3))
                xi(2) = y_cc(cell(2))
                eta(2) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2), cell(3))
                xi(3) = y_cc(cell(2) + 1)
                eta(3) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 1, cell(3))
            elseif (i == 3) then
                xi(1) = z_cc(cell(3) - 1)
                eta(1) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 1)
                xi(2) = z_cc(cell(3))
                eta(2) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3))
                xi(3) = z_cc(cell(3) + 1)
                eta(3) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 1)
            end if

            L(1) = ((pos - xi(2))*(pos - xi(3)))/((xi(1) - xi(2))*(xi(1) - xi(3)))
            L(2) = ((pos - xi(1))*(pos - xi(3)))/((xi(2) - xi(1))*(xi(2) - xi(3)))
            L(3) = ((pos - xi(1))*(pos - xi(2)))/((xi(3) - xi(1))*(xi(3) - xi(2)))

            v = L(1)*eta(1) + L(2)*eta(2) + L(3)*eta(3)
        elseif (fd_order == 4) then
            if (i == 1) then
                xi(1) = x_cc(cell(1) - 2)
                eta(1) = q_prim_vf(momxb)%sf(cell(1) - 2, cell(2), cell(3))
                xi(2) = x_cc(cell(1) - 1)
                eta(2) = q_prim_vf(momxb)%sf(cell(1) - 1, cell(2), cell(3))
                xi(3) = x_cc(cell(1))
                eta(3) = q_prim_vf(momxb)%sf(cell(1), cell(2), cell(3))
                xi(4) = x_cc(cell(1) + 1)
                eta(4) = q_prim_vf(momxb)%sf(cell(1) + 1, cell(2), cell(3))
                xi(5) = x_cc(cell(1) + 2)
                eta(5) = q_prim_vf(momxb)%sf(cell(1) + 2, cell(2), cell(3))
            elseif (i == 2) then
                xi(1) = y_cc(cell(2) - 2)
                eta(1) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 2, cell(3))
                xi(2) = y_cc(cell(2) - 1)
                eta(2) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 1, cell(3))
                xi(3) = y_cc(cell(2))
                eta(3) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2), cell(3))
                xi(4) = y_cc(cell(2) + 1)
                eta(4) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 1, cell(3))
                xi(5) = y_cc(cell(2) + 2)
                eta(5) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 2, cell(3))
            elseif (i == 3) then
                xi(1) = z_cc(cell(3) - 2)
                eta(1) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 2)
                xi(2) = z_cc(cell(3) - 1)
                eta(2) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 1)
                xi(3) = z_cc(cell(3))
                eta(3) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3))
                xi(4) = z_cc(cell(3) + 1)
                eta(4) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 1)
                xi(5) = z_cc(cell(3) + 2)
                eta(5) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 2)
            end if

            L(1) = (pos-xi(2))*(pos-xi(3))*(pos-xi(4))*(pos-xi(5)) / &
                  ((xi(1)-xi(2))*(xi(1)-xi(3))*(xi(1)-xi(4))*(xi(1)-xi(5)))
            L(2) = (pos-xi(1))*(pos-xi(3))*(pos-xi(4))*(pos-xi(5)) / &
                  ((xi(2)-xi(1))*(xi(2)-xi(3))*(xi(2)-xi(4))*(xi(2)-xi(5)))
            L(3) = (pos-xi(1))*(pos-xi(2))*(pos-xi(4))*(pos-xi(5)) / &
                  ((xi(3)-xi(1))*(xi(3)-xi(2))*(xi(3)-xi(4))*(xi(3)-xi(5)))
            L(4) = (pos-xi(1))*(pos-xi(2))*(pos-xi(3))*(pos-xi(5)) / &
                  ((xi(4)-xi(1))*(xi(4)-xi(2))*(xi(4)-xi(3))*(xi(4)-xi(5)))
            L(5) = (pos-xi(1))*(pos-xi(2))*(pos-xi(3))*(pos-xi(4)) / &
                  ((xi(5)-xi(1))*(xi(5)-xi(2))*(xi(5)-xi(3))*(xi(5)-xi(4)))


            v = L(1)*eta(1) + L(2)*eta(2) + L(3)*eta(3) + L(4)*eta(4) + L(5)*eta(5)
        end if

    end function f_interpolate_velocity

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
            !! @param i Direction of the velocity (1: x, 2: y, 3: z)
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return a Acceleration of the particle in direction i
    function f_get_particle_force(pos, rad, vel_p, mass_p, Re, gamm, rho, vol_frac, cell, i, q_prim_vf) result(force)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pos, rad, mass_p, Re, gamm, rho, vol_frac
        integer, dimension(3), intent(in) :: cell
        real(wp), intent(in) :: vel_p(:)
        integer, intent(in) :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp) :: a, dp, vol, force
        real(wp) :: v_rel
        real(wp), dimension(fd_order) :: xi, eta, L
        real(wp) :: particle_diam, gas_mu, vmag, cson
        real(wp) :: slip_velocity_x, slip_velocity_y, slip_velocity_z, beta
        real(wp), dimension(3) :: fluid_vel
        integer :: dir

        if (fd_order <= 2) then
            if (i == 1) then
                dp = (q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3)) - &
                      q_prim_vf(E_idx)%sf(cell(1) - 1, cell(2), cell(3)))/ &
                     (x_cc(cell(1) + 1) - x_cc(cell(1) - 1))
            elseif (i == 2) then
                dp = (q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3)) - &
                      q_prim_vf(E_idx)%sf(cell(1), cell(2) - 1, cell(3)))/ &
                     (y_cc(cell(2) + 1) - y_cc(cell(2) - 1))
            elseif (i == 3) then
                dp = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 1) - &
                      q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) - 1))/ &
                     (z_cc(cell(3) + 1) - z_cc(cell(3) - 1))
            end if
        elseif (fd_order == 4) then
            if (i == 1) then
                xi(1) = x_cc(cell(1) - 1)
                eta(1) = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1) - 2, cell(2), cell(3)))/ &
                         (x_cc(cell(1)) - x_cc(cell(1) - 2))
                xi(2) = x_cc(cell(1))
                eta(2) = (q_prim_vf(E_idx)%sf(cell(1) + 1, cell(2), cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1) - 1, cell(2), cell(3)))/ &
                         (x_cc(cell(1) + 1) - x_cc(cell(1) - 1))
                xi(3) = x_cc(cell(1) + 1)
                eta(3) = (q_prim_vf(E_idx)%sf(cell(1) + 2, cell(2), cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)))/ &
                         (x_cc(cell(1) + 2) - x_cc(cell(1)))
            elseif (i == 2) then
                xi(1) = y_cc(cell(2) - 1)
                eta(1) = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2) - 2, cell(3)))/ &
                         (y_cc(cell(2)) - y_cc(cell(2) - 2))
                xi(2) = y_cc(cell(2))
                eta(2) = (q_prim_vf(E_idx)%sf(cell(1), cell(2) + 1, cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2) - 1, cell(3)))/ &
                         (y_cc(cell(2) + 1) - y_cc(cell(2) - 1))
                xi(3) = y_cc(cell(2) + 1)
                eta(3) = (q_prim_vf(E_idx)%sf(cell(1), cell(2) + 2, cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)))/ &
                         (y_cc(cell(2) + 2) - y_cc(cell(2)))
            elseif (i == 3) then
                xi(1) = z_cc(cell(3) - 1)
                eta(1) = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) - 2))/ &
                         (z_cc(cell(3)) - z_cc(cell(3) - 2))
                xi(2) = z_cc(cell(3))
                eta(2) = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 1) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) - 1))/ &
                         (z_cc(cell(3) + 1) - z_cc(cell(3) - 1))
                xi(3) = z_cc(cell(3) + 1)
                eta(3) = (q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3) + 2) - &
                          q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)))/ &
                         (z_cc(cell(3) + 2) - z_cc(cell(3)))
            end if

            L(1) = ((pos - xi(2))*(pos - xi(3)))/((xi(1) - xi(2))*(xi(1) - xi(3)))
            L(2) = ((pos - xi(1))*(pos - xi(3)))/((xi(2) - xi(1))*(xi(2) - xi(3)))
            L(3) = ((pos - xi(1))*(pos - xi(2)))/((xi(3) - xi(1))*(xi(3) - xi(2)))

            dp = L(1)*eta(1) + L(2)*eta(2) + L(3)*eta(3)
        end if

        vol = (4._wp/3._wp)*pi*(rad**3._wp)

        fluid_vel(1) = f_interpolate_velocity(pos, cell, 1, q_prim_vf)
        fluid_vel(2) = f_interpolate_velocity(pos, cell, 2, q_prim_vf)
        if (num_dims == 3) then
          fluid_vel(3) = f_interpolate_velocity(pos, cell, 3, q_prim_vf)
        endif


        do dir = 1,3
          if (abs(q_prim_vf(momxb+dir-1)%sf(cell(1), cell(2), cell(3))) < 1.e-8) then
            fluid_vel(dir) = q_prim_vf(momxb+dir-1)%sf(cell(1), cell(2), cell(3)) !0th order interpolation if fluid velocity is small
          endif
          if (fd_order>1) then
            if (.not. ieee_is_finite(fluid_vel(dir))) then
              fluid_vel(dir) =  q_prim_vf(momxb+dir-1)%sf(cell(1), cell(2), cell(3))
            endif
          endif
        enddo

        v_rel = vel_p(i) - fluid_vel(i)


        !!!!!!! Quasi-steady Drag Force Paramaeters
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
        cson = sqrt((gamm*q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3)))/rho) !gamma*P/rho
        gas_mu = Re

        if (ABS(v_rel) .le. 1.e-7_wp) then
          v_rel = 0._wp
        endif

        force = 0._wp

        ! Step 1: Force component quasi-steady
        if (lag_params%qs_drag_model == 1) then
            beta = QS_Parmar(rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 2) then
            beta = QS_Osnes(rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 3) then
            beta = QS_ModifiedParmar(rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 4) then
            beta = QS_Gidaspow(rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        end if

        if (.not. ieee_is_finite(force)) then
          force = 0._wp
        endif


        ! Step 1.1: Stokes drag?
        if (lag_params%stokes_drag == 1) then ! Free slip Stokes drag
            force = force - 4._wp*pi*gas_mu*rad*v_rel
        elseif (lag_params%stokes_drag == 2) then ! No slip Stokes drag
            force = force - 6._wp*pi*gas_mu*rad*v_rel
        else
            !No stokes drag
        end if

        ! Step 2: Pressure Gradient Force
        if (lag_params%pressure_force) then
            force = force - vol*dp
        end if

        ! Step 3: Gravitational Force
        if (lag_params%gravity_force) then
            force = force + (mass_p)*accel_bf(i)
        end if

    end function f_get_particle_force


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

        b3 = dmin1(sqrt(20.0_wp*mp), 1.0_wp)* &
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
    ! NOTE: Re<45 Rarified fomu_fluidla of Loth et al has been redefined by Balachandar
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

        b3 = dmin1(sqrt(20.0_wp*mp), 1.0_wp)* &
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
