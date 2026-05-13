!>
!! @file m_particle_bed.fpp
!! @brief Generates particle beds: converts particle_bed specifications into
!!        individual sphere/circle patch_ib entries before MPI broadcast.

!> @brief Generates particle beds by converting particle_bed patch specifications into individual immersed boundary patches before
!! MPI broadcast.
module m_particle_bed

    use m_global_parameters
    use m_constants

    implicit none

    private

    public :: s_generate_particle_beds

contains

    !> Generate all particle beds and append the resulting particles to patch_ib. Called on rank 0 only, before
    !! s_mpi_bcast_user_inputs. Uses a spatial hash grid (cell size = min_dist) so each candidate requires only 3^dim distance
    !! checks on average instead of O(n).
    impure subroutine s_generate_particle_beds()

        integer               :: b, ib_idx, geom
        integer               :: n_placed, n_total_placed
        integer(8)            :: n_attempts, max_attempts
        real(wp)              :: xmin, xmax, ymin, ymax, zmin, zmax, min_dist
        real(wp)              :: rx, ry, rz, dist
        real(wp)              :: t_start, t_end
        integer               :: seed
        logical               :: overlaps
        real(wp), allocatable :: placed(:,:)

        ! Spatial hash grid
        integer              :: hash_size, slot
        integer              :: bx, by, bz, nbx, nby, nbz
        integer              :: dx_b, dy_b, dz_b, dz_lo, dz_hi, j
        integer, allocatable :: hash_head(:), chain_next(:)

        if (num_particle_beds == 0) return

        call cpu_time(t_start)
        n_total_placed = 0

        do b = 1, num_particle_beds
            xmin = particle_bed(b)%x_centroid - 0.5_wp*particle_bed(b)%length_x
            xmax = particle_bed(b)%x_centroid + 0.5_wp*particle_bed(b)%length_x
            ymin = particle_bed(b)%y_centroid - 0.5_wp*particle_bed(b)%length_y
            ymax = particle_bed(b)%y_centroid + 0.5_wp*particle_bed(b)%length_y
            zmin = particle_bed(b)%z_centroid - 0.5_wp*particle_bed(b)%length_z
            zmax = particle_bed(b)%z_centroid + 0.5_wp*particle_bed(b)%length_z

            min_dist = 2._wp*particle_bed(b)%radius + particle_bed(b)%min_spacing

            if (p == 0) then
                geom = 2  ! circle for 2D
                dz_lo = 0
                dz_hi = 0
            else
                geom = 8  ! sphere for 3D
                dz_lo = -1
                dz_hi = 1
            end if

            max_attempts = int(particle_bed(b)%num_particles, 8)*1000_8
            n_placed = 0
            n_attempts = 0
            seed = particle_bed(b)%seed
            if (seed == 0) seed = 1 + b*1013904223

            allocate (placed(3, particle_bed(b)%num_particles))

            ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
            ! previous occupant of its bucket.
            hash_size = max(16, 4*particle_bed(b)%num_particles)
            allocate (hash_head(hash_size))
            allocate (chain_next(particle_bed(b)%num_particles))
            hash_head = -1
            chain_next = -1

            do while (n_placed < particle_bed(b)%num_particles .and. n_attempts < max_attempts)
                n_attempts = n_attempts + 1

                rx = xmin + f_xorshift(seed)*(xmax - xmin)
                ry = ymin + f_xorshift(seed)*(ymax - ymin)
                if (p == 0) then
                    rz = particle_bed(b)%z_centroid
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

                    num_ibs = num_ibs + 1
                    ib_idx = num_ibs

                    patch_ib(ib_idx)%gbl_patch_id = ib_idx
                    patch_ib(ib_idx)%geometry = geom
                    patch_ib(ib_idx)%x_centroid = rx
                    patch_ib(ib_idx)%y_centroid = ry
                    patch_ib(ib_idx)%z_centroid = rz
                    patch_ib(ib_idx)%angles(1) = 0._wp
                    patch_ib(ib_idx)%angles(2) = 0._wp
                    patch_ib(ib_idx)%angles(3) = 0._wp
                    patch_ib(ib_idx)%vel(1) = 0._wp
                    patch_ib(ib_idx)%vel(2) = 0._wp
                    patch_ib(ib_idx)%vel(3) = 0._wp
                    patch_ib(ib_idx)%angular_vel(1) = 0._wp
                    patch_ib(ib_idx)%angular_vel(2) = 0._wp
                    patch_ib(ib_idx)%angular_vel(3) = 0._wp
                    patch_ib(ib_idx)%radius = particle_bed(b)%radius
                    patch_ib(ib_idx)%mass = particle_bed(b)%mass
                    patch_ib(ib_idx)%moving_ibm = particle_bed(b)%moving_ibm
                end if
            end do

            if (n_placed < particle_bed(b)%num_particles) then
                print *, "Error :: Failed to place all IBs ib particle bed"
                stop
            end if

            n_total_placed = n_total_placed + n_placed
            deallocate (placed, hash_head, chain_next)
        end do

        call cpu_time(t_end)
        if (proc_rank == 0) print '(a,i0,a,f0.3,a)', 'Particle beds placed ', n_total_placed, ' particles in ', t_end - t_start, &
            & ' seconds.'

    end subroutine s_generate_particle_beds

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

end module m_particle_bed
