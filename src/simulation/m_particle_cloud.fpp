!>
!! @file m_particle_cloud.fpp
!! @brief Generates particle beds: converts particle_cloud specifications into
!!        individual sphere/circle particle_cloud_ibs entries before reduction.

!> @brief Generates particle beds by converting particle_cloud patch specifications into individual immersed boundary patches before
!! domain reduction. Each rank runs the same deterministic placement so no MPI broadcast of particle positions is needed.
module m_particle_cloud

    use m_global_parameters
    use m_constants
    use m_mpi_common

    implicit none

    private

    public :: s_generate_particle_clouds

contains

    !> Generate all particle beds and fill particle_cloud_ibs. Called on all ranks before s_reduce_ib_patch_array.
    impure subroutine s_generate_particle_clouds(particle_cloud_ibs)

        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: particle_cloud_ibs
        integer                                                           :: cloud_idx, ib_idx, n_total_particles
        real(wp)                                                          :: t_start, t_end

        if (num_particle_clouds == 0) then
            allocate (particle_cloud_ibs(0))
            return
        end if

        call cpu_time(t_start)

        ! Pre-count total particles across all beds so particle_cloud_ibs can be allocated exactly once.
        n_total_particles = 0
        do cloud_idx = 1, num_particle_clouds
            n_total_particles = n_total_particles + particle_cloud(cloud_idx)%num_particles
        end do
        allocate (particle_cloud_ibs(n_total_particles))

        ib_idx = 0  ! index into particle_cloud_ibs

        do cloud_idx = 1, num_particle_clouds
            select case (particle_cloud(cloud_idx)%packing_method)
            case (1)  ! random box packing method
                call s_particle_cloud_random_box(cloud_idx, ib_idx, particle_cloud_ibs)
            case (2)  ! lattice packing method
                call s_particle_cloud_lattice(cloud_idx, ib_idx, particle_cloud_ibs)
            end select
        end do

        call cpu_time(t_end)
        if (proc_rank == 0) print '(a,i0,a,f0.3,a)', 'Particle beds placed ', ib_idx, ' particles in ', t_end - t_start, ' seconds.'

    end subroutine s_generate_particle_clouds

    !> Generates a random distributions of particles in a box with a minimum spacing
    subroutine s_particle_cloud_random_box(cloud_idx, ib_idx, particle_cloud_ibs)

        integer, intent(in)                                    :: cloud_idx
        integer, intent(inout)                                 :: ib_idx
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs
        integer                                                :: n_placed, geom, seed
        integer(8)                                             :: n_attempts, max_attempts
        real(wp)                                               :: xmin, xmax, ymin, ymax, zmin, zmax, min_dist
        real(wp)                                               :: rx, ry, rz, dist
        logical                                                :: overlaps
        real(wp), allocatable                                  :: placed(:,:)
        integer                                                :: hash_size, slot
        integer                                                :: bx, by, bz, nbx, nby, nbz
        integer                                                :: dx_b, dy_b, dz_b, dz_lo, dz_hi, j
        integer, allocatable                                   :: hash_head(:), chain_next(:)

        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
        xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
        ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
        zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing

        if (p == 0) then
            geom = 2  ! circle for 2D
            dz_lo = 0
            dz_hi = 0
        else
            geom = 8  ! sphere for 3D
            dz_lo = -1
            dz_hi = 1
        end if

        max_attempts = int(particle_cloud(cloud_idx)%num_particles, 8)*1000_8
        n_placed = 0
        n_attempts = 0
        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = 1 + cloud_idx*1013904223

        allocate (placed(3, particle_cloud(cloud_idx)%num_particles))

        ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
        ! previous occupant of its bucket.
        hash_size = max(16, 4*particle_cloud(cloud_idx)%num_particles)
        allocate (hash_head(hash_size))
        allocate (chain_next(particle_cloud(cloud_idx)%num_particles))
        hash_head = -1
        chain_next = -1

        do while (n_placed < particle_cloud(cloud_idx)%num_particles .and. n_attempts < max_attempts)
            n_attempts = n_attempts + 1

            rx = xmin + f_xorshift(seed)*(xmax - xmin)
            ry = ymin + f_xorshift(seed)*(ymax - ymin)
            if (p == 0) then
                rz = particle_cloud(cloud_idx)%z_centroid
            else
                rz = zmin + f_xorshift(seed)*(zmax - zmin)
            end if

            bx = int(floor(rx/min_dist))
            by = int(floor(ry/min_dist))
            bz = 0
            if (p /= 0) bz = int(floor(rz/min_dist))

            ! Check 3x3(x3) neighboring bins - O(1) average via hash lookup
            overlaps = .false.
            outer: do dx_b = -1, 1
                do dy_b = -1, 1
                    do dz_b = dz_lo, dz_hi
                        nbx = bx + dx_b
                        nby = by + dy_b
                        nbz = bz + dz_b
                        slot = f_bin_hash(nbx, nby, nbz, hash_size)
                        j = hash_head(slot)
                        do while (j > 0)
                            if (p == 0) then
                                dist = sqrt((rx - placed(1, j))**2 + (ry - placed(2, j))**2)
                            else
                                dist = sqrt((rx - placed(1, j))**2 + (ry - placed(2, j))**2 + (rz - placed(3, j))**2)
                            end if
                            if (dist < min_dist) then
                                overlaps = .true.
                                exit outer
                            end if
                            j = chain_next(j)
                        end do
                    end do
                end do
            end do outer

            if (.not. overlaps) then
                n_placed = n_placed + 1
                placed(1, n_placed) = rx
                placed(2, n_placed) = ry
                placed(3, n_placed) = rz

                ! Insert into hash grid as head of bucket chain
                slot = f_bin_hash(bx, by, bz, hash_size)
                chain_next(n_placed) = hash_head(slot)
                hash_head(slot) = n_placed

                call s_add_cloud_particle(cloud_idx, ib_idx, geom, rx, ry, rz, particle_cloud_ibs)
            end if
        end do

        if (n_placed < particle_cloud(cloud_idx)%num_particles) then
            call s_mpi_abort("Error :: Failed to place all particles in particle bed")
        end if

        deallocate (placed, hash_head, chain_next)

    end subroutine s_particle_cloud_random_box

    !> Places particles on the optimally dense lattice for the cloud region: a triangular lattice in 2D, a face-centered cubic
    !! lattice in 3D. The lattice spacing is set by the particle density (num_particles over the region area/volume); if that
    !! spacing falls below the required centre-to-centre distance (2*radius + min_spacing), the region is too dense and the run is
    !! aborted.
    subroutine s_particle_cloud_lattice(cloud_idx, ib_idx, particle_cloud_ibs)

        integer, intent(in)                                    :: cloud_idx
        integer, intent(inout)                                 :: ib_idx
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs
        integer                                                :: n_placed, n_target, geom
        integer                                                :: row, col, ncx, ncy, ix, jy, kz, b
        real(wp)                                               :: xmin, xmax, ymin, ymax, zmin, zmax, min_dist
        real(wp)                                               :: spacing, row_dy, cell, x0, px, py
        real(wp), dimension(4)                                 :: bx_off, by_off, bz_off

        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
        xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
        ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
        zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        n_target = particle_cloud(cloud_idx)%num_particles
        n_placed = 0

        if (p == 0) then
            geom = 2  ! circle for 2D
            ! Triangular lattice: area per particle = (sqrt(3)/2)*spacing**2.
            spacing = sqrt(2._wp*(xmax - xmin)*(ymax - ymin)/(sqrt(3._wp)*real(n_target, wp)))
        else
            geom = 8  ! sphere for 3D
            ! Face-centered cubic lattice: volume per particle = spacing**3/sqrt(2).
            spacing = (sqrt(2._wp)*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)/real(n_target, wp))**(1._wp/3._wp)
        end if

        if (spacing < min_dist) then
            call s_mpi_abort("Error :: Particle cloud is too dense for lattice packing; " &
                             & // "reduce num_particles or min_spacing, or enlarge the cloud region")
        end if

        if (p == 0) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_dy = spacing*sqrt(3._wp)/2._wp
            row = 0
            do while (n_placed < n_target)
                py = ymin + real(row, wp)*row_dy
                x0 = xmin
                if (mod(row, 2) == 1) x0 = xmin + 0.5_wp*spacing
                col = 0
                px = x0
                do while (px <= xmax .and. n_placed < n_target)
                    call s_add_cloud_particle(cloud_idx, ib_idx, geom, px, py, particle_cloud(cloud_idx)%z_centroid, &
                                              & particle_cloud_ibs)
                    n_placed = n_placed + 1
                    col = col + 1
                    px = x0 + real(col, wp)*spacing
                end do
                row = row + 1
            end do
        else
            ! Face-centered cubic lattice via the conventional cubic cell (side = spacing*sqrt(2)) and its four basis points.
            cell = spacing*sqrt(2._wp)
            bx_off = [0._wp, 0.5_wp, 0.5_wp, 0._wp]*cell
            by_off = [0._wp, 0.5_wp, 0._wp, 0.5_wp]*cell
            bz_off = [0._wp, 0._wp, 0.5_wp, 0.5_wp]*cell
            ncx = max(1, ceiling((xmax - xmin)/cell))
            ncy = max(1, ceiling((ymax - ymin)/cell))
            kz = 0
            do while (n_placed < n_target)
                do jy = 0, ncy - 1
                    do ix = 0, ncx - 1
                        do b = 1, 4
                            if (n_placed >= n_target) exit
                            call s_add_cloud_particle(cloud_idx, ib_idx, geom, xmin + real(ix, wp)*cell + bx_off(b), &
                                                      & ymin + real(jy, wp)*cell + by_off(b), zmin + real(kz, &
                                                      & wp)*cell + bz_off(b), particle_cloud_ibs)
                            n_placed = n_placed + 1
                        end do
                    end do
                end do
                kz = kz + 1
            end do
        end if

    end subroutine s_particle_cloud_lattice

    !> Writes a single placed particle into particle_cloud_ibs at the next free slot, advancing ib_idx. Shared by all packing
    !! methods so the per-particle ib_patch_parameters setup stays in one place.
    subroutine s_add_cloud_particle(cloud_idx, ib_idx, geom, px, py, pz, particle_cloud_ibs)

        integer, intent(in)                                    :: cloud_idx, geom
        integer, intent(inout)                                 :: ib_idx
        real(wp), intent(in)                                   :: px, py, pz
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs

        ib_idx = ib_idx + 1

        ! gbl_patch_id is relative within particle_cloud_ibs here; s_reduce_ib_patch_array adjusts to global indexing.
        particle_cloud_ibs(ib_idx)%gbl_patch_id = ib_idx
        particle_cloud_ibs(ib_idx)%geometry = geom
        particle_cloud_ibs(ib_idx)%x_centroid = px
        particle_cloud_ibs(ib_idx)%y_centroid = py
        particle_cloud_ibs(ib_idx)%z_centroid = pz
        particle_cloud_ibs(ib_idx)%step_x_centroid = px
        particle_cloud_ibs(ib_idx)%step_y_centroid = py
        particle_cloud_ibs(ib_idx)%step_z_centroid = pz
        particle_cloud_ibs(ib_idx)%angles(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_angles(:) = 0._wp
        particle_cloud_ibs(ib_idx)%vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%angular_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_angular_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%force(:) = 0._wp
        particle_cloud_ibs(ib_idx)%torque(:) = 0._wp
        particle_cloud_ibs(ib_idx)%centroid_offset(:) = 0._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix = 0._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(1, 1) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(2, 2) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(3, 3) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix_inverse = particle_cloud_ibs(ib_idx)%rotation_matrix
        particle_cloud_ibs(ib_idx)%radius = particle_cloud(cloud_idx)%radius
        particle_cloud_ibs(ib_idx)%mass = particle_cloud(cloud_idx)%mass
        particle_cloud_ibs(ib_idx)%moment = dflt_real
        particle_cloud_ibs(ib_idx)%moving_ibm = particle_cloud(cloud_idx)%moving_ibm
        particle_cloud_ibs(ib_idx)%slip = .false.

    end subroutine s_add_cloud_particle

    !> Xorshift PRNG. Advances seed in-place and returns a value in [0, 1).
    function f_xorshift(seed) result(rval)

        integer, intent(inout) :: seed
        real(wp)               :: rval

        seed = ieor(seed, ishft(seed, 13))
        seed = ieor(seed, ishft(seed, -17))
        seed = ieor(seed, ishft(seed, 5))

        rval = abs(real(seed, wp))/real(huge(seed), wp)

    end function f_xorshift

    !> Hash bin coordinates to a 1-indexed slot in [1, hash_size]. Uses large prime multipliers to spread bins across buckets. Hash
    !! collisions are benign: the distance check catches false neighbours.
    function f_bin_hash(bx, by, bz, hash_size) result(slot)

        integer, intent(in) :: bx, by, bz, hash_size
        integer             :: slot
        integer(8)          :: key

        key = ieor(ieor(int(bx, 8)*73856093_8, int(by, 8)*19349663_8), int(bz, 8)*83492791_8)
        slot = int(mod(abs(key), int(hash_size, 8))) + 1

    end function f_bin_hash

end module m_particle_cloud
