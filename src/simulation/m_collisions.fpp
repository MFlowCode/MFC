!>
!! @file
!! @brief Contains module m_collisions

#:include 'macros.fpp'

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the
!! flow state
module m_collisions

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Definitions of the global parameters
    use m_helper
    use m_helper_basic       !< Functions to compare floating point numbers
    use m_constants
    use m_compute_levelset
    use m_ib_patches
    use m_model
    use m_mpi_proxy

    implicit none

    private; public :: s_apply_IB_collisions, s_apply_EL_collisions, s_initialize_collisions_module, &
        & s_finalize_collisions_module, f_neighborhood_ranks_own_location, collisions_active
    ! overlap distances for computing collisions
    integer, allocatable, dimension(:,:)  :: collision_lookup
    real(wp), allocatable, dimension(:,:) :: wall_overlap_distances
    real(wp)                              :: spring_stiffness, damping_parameter
    $:GPU_DECLARE(create='[spring_stiffness, damping_parameter]')
    $:GPU_DECLARE(create='[collision_lookup, wall_overlap_distances]')

    !> true when any IB-IB or IB-wall contact was detected on this rank since the last adaptive-dt computation
    logical :: collisions_active

contains

    subroutine s_initialize_collisions_module()

        real(wp) :: e
        integer  :: max_collisions

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness]')

        ! one lookup sized for the larger of the IB and EL pair counts; EL allows 12 contacts per particle (3D close packing)
        max_collisions = 0
        if (ib) max_collisions = num_local_ibs_max*27*8
        if (particles_lagrange) max_collisions = max(max_collisions, 12*particle_params%nparticles_glb)
        if (max_collisions > 0) then
            @:ALLOCATE(collision_lookup(max_collisions, 4))
            collision_lookup = 0
            $:GPU_UPDATE(device='[collision_lookup]')
        end if

        if (ib) then
            @:ALLOCATE(wall_overlap_distances(num_local_ibs_max*27, 6))
            wall_overlap_distances = 0
            $:GPU_UPDATE(device='[wall_overlap_distances]')
        end if

        collisions_active = .false.
        $:GPU_UPDATE(device='[coefficient_of_friction]')
        $:GPU_UPDATE(device='[ib_bc_x%beg, ib_bc_x%end, ib_bc_y%beg, ib_bc_y%end, ib_bc_z%beg, ib_bc_z%end]')

    end subroutine s_initialize_collisions_module

    subroutine s_apply_IB_collisions(ghost_points, num_gps, ib_markers, forces, torques)

        type(ghost_point), dimension(:), intent(in)    :: ghost_points
        integer, intent(in)                            :: num_gps
        type(integer_field), intent(in)                :: ib_markers
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: num_considered_collisions, i, l, pid1, pid2, decoded_id
        integer, dimension(3)                          :: periodic_image_1, periodic_image_2
        logical                                        :: any_wall_collision
        real(wp), dimension(3)                         :: pos, vel, angular_vel, force, torque
        real(wp), dimension(3)                         :: pos_2, vel_2, angular_vel_2, torque_2
        real(wp), dimension(6)                         :: wall_overlaps
        real(wp)                                       :: radius, mass, max_overlap

        ! return if no collisions

        if (collision_model == 0) return

        ! iterating over particles, determine if they are colliding with the wall
        max_overlap = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, pos, wall_overlaps]', reduction='[[max_overlap]]', reductionOp='[max]')
        do i = 1, num_ibs
            pos = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
            if (num_dims == 3) pos(3) = patch_ib(i)%z_centroid
            call s_detect_wall_overlaps(pos, patch_ib(i)%radius, wall_overlaps)
            wall_overlap_distances(i,:) = wall_overlaps
            max_overlap = max(max_overlap, maxval(wall_overlaps))
        end do
        $:END_GPU_PARALLEL_LOOP()
        any_wall_collision = max_overlap > 0._wp

        ! apply wall collisions
        $:GPU_PARALLEL_LOOP(private='[i, l, pos, vel, radius, angular_vel, mass, wall_overlaps, force, torque]', &
                            & copy='[forces, torques]')
        do i = 1, num_ibs
            pos = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
            if (num_dims == 3) pos(3) = patch_ib(i)%z_centroid
            vel = patch_ib(i)%vel
            angular_vel = patch_ib(i)%angular_vel
            radius = patch_ib(i)%radius
            mass = patch_ib(i)%mass
            wall_overlaps = wall_overlap_distances(i,:)

            force = 0._wp; torque = 0._wp
            select case (collision_model)
            case(1)  ! soft sphere model
                call s_soft_sphere_wall_force(pos, vel, radius, mass, wall_overlaps, force, angular_vel, torque)
            end select

            ! one thread per IB, so no atomics
            do l = 1, num_dims
                forces(i, l) = forces(i, l) + force(l)
                torques(i, l) = torques(i, l) + torque(l)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! performs a GCIBM-specific collision detection reduction step and populates into the collision lookup
        call s_detect_ib_collisions(ghost_points, ib_markers, num_gps, num_considered_collisions)

        ! accumulate across RK stages; consumed (and reset) by s_compute_dt once per time step
        collisions_active = collisions_active .or. any_wall_collision .or. (num_considered_collisions > 0)

        ! apply IB-IB collisions
        $:GPU_PARALLEL_LOOP(private='[i, l, pid1, pid2, decoded_id, periodic_image_1, periodic_image_2, pos, vel, angular_vel, &
                            & pos_2, vel_2, angular_vel_2, force, torque, torque_2]', copy='[forces, torques]')
        do i = 1, num_considered_collisions
            pid1 = collision_lookup(i, 1)
            pid2 = collision_lookup(i, 2)
            if (pid1 <= 0 .or. pid2 <= 0) cycle

            ! the encoded ids carry which periodic image of each IB is in contact
            call s_decode_patch_periodicity(collision_lookup(i, 3), decoded_id, periodic_image_1(1), periodic_image_1(2), &
                                            & periodic_image_1(3))
            call s_decode_patch_periodicity(collision_lookup(i, 4), decoded_id, periodic_image_2(1), periodic_image_2(2), &
                                            & periodic_image_2(3))
            pos = [patch_ib(pid1)%x_centroid, patch_ib(pid1)%y_centroid, 0._wp]
            pos_2 = [patch_ib(pid2)%x_centroid, patch_ib(pid2)%y_centroid, 0._wp]
            if (num_dims == 3) then
                pos(3) = patch_ib(pid1)%z_centroid
                pos_2(3) = patch_ib(pid2)%z_centroid
            end if
            do l = 1, num_dims
                pos(l) = pos(l) + real(periodic_image_1(l), wp)*(glb_bounds(l)%end - glb_bounds(l)%beg)
                pos_2(l) = pos_2(l) + real(periodic_image_2(l), wp)*(glb_bounds(l)%end - glb_bounds(l)%beg)
            end do

            ! only overlapping pairs are passed to the collision model
            if (patch_ib(pid1)%radius + patch_ib(pid2)%radius - norm2(pos_2 - pos) <= 0._wp) cycle

            vel = patch_ib(pid1)%vel
            vel_2 = patch_ib(pid2)%vel
            angular_vel = patch_ib(pid1)%angular_vel
            angular_vel_2 = patch_ib(pid2)%angular_vel

            force = 0._wp; torque = 0._wp; torque_2 = 0._wp
            select case (collision_model)
            case (1)  ! soft sphere model
                call s_perform_soft_sphere_collision(pos, vel, patch_ib(pid1)%radius, patch_ib(pid1)%mass, pos_2, vel_2, &
                                                     & patch_ib(pid2)%radius, patch_ib(pid2)%mass, force, angular_vel, &
                                                     & angular_vel_2, torque, torque_2)
            end select

            ! equal and opposite force on the second IB
            do l = 1, num_dims
                $:GPU_ATOMIC(atomic='update')
                forces(pid1, l) = forces(pid1, l) + force(l)
                $:GPU_ATOMIC(atomic='update')
                torques(pid1, l) = torques(pid1, l) + torque(l)
                $:GPU_ATOMIC(atomic='update')
                forces(pid2, l) = forces(pid2, l) - force(l)
                $:GPU_ATOMIC(atomic='update')
                torques(pid2, l) = torques(pid2, l) + torque_2(l)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_IB_collisions

    !> Add the wall and particle-particle contact accelerations of the EL particles to accelerations, dividing by the same inertia
    !! as the hydrodynamic force (mass plus added mass). Only pairs on this rank are found; there is no halo yet.
    subroutine s_apply_EL_collisions(num_particles, positions, velocities, radii, masses, added_masses, accelerations)

        integer, intent(in)                     :: num_particles
        real(wp), dimension(:,:,:), intent(in)  :: positions, velocities  !< (particle, component, RK level); level 2 is current
        real(wp), dimension(:), intent(in)      :: radii, masses, added_masses
        real(wp), dimension(:,:), intent(inout) :: accelerations
        integer                                 :: num_considered_collisions, i, l, pid1, pid2
        logical                                 :: any_wall_collision
        real(wp), dimension(3)                  :: pos, vel, pos_2, vel_2, force
        real(wp), dimension(6)                  :: wall_overlaps
        real(wp)                                :: max_overlap

        if (collision_model == 0) return

        ! apply wall collisions; one thread per particle, so no atomics
        max_overlap = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, l, pos, vel, wall_overlaps, force]', reduction='[[max_overlap]]', reductionOp='[max]')
        do i = 1, num_particles
            pos = positions(i,:,2)
            if (num_dims == 2) pos(3) = 0._wp
            vel = velocities(i,:,2)
            call s_detect_wall_overlaps(pos, radii(i), wall_overlaps)
            max_overlap = max(max_overlap, maxval(wall_overlaps))

            force = 0._wp
            select case (collision_model)
            case (1)  ! soft sphere model
                call s_soft_sphere_wall_force(pos, vel, radii(i), masses(i), wall_overlaps, force)
            end select

            do l = 1, num_dims
                accelerations(i, l) = accelerations(i, l) + force(l)/(masses(i) + added_masses(i))
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        any_wall_collision = max_overlap > 0._wp

        ! naive all-to-all search; every pair in the lookup overlaps
        call s_detect_el_collisions(num_particles, positions, radii, num_considered_collisions)

        ! accumulate across RK stages; consumed (and reset) by s_compute_dt once per time step
        collisions_active = collisions_active .or. any_wall_collision .or. (num_considered_collisions > 0)

        ! apply particle-particle collisions
        $:GPU_PARALLEL_LOOP(private='[i, l, pid1, pid2, pos, vel, pos_2, vel_2, force]')
        do i = 1, num_considered_collisions
            pid1 = collision_lookup(i, 1)
            pid2 = collision_lookup(i, 2)
            pos = positions(pid1,:,2)
            pos_2 = positions(pid2,:,2)
            if (num_dims == 2) then
                pos(3) = 0._wp
                pos_2(3) = 0._wp
            end if
            vel = velocities(pid1,:,2)
            vel_2 = velocities(pid2,:,2)

            force = 0._wp
            select case (collision_model)
            case (1)  ! soft sphere model
                call s_perform_soft_sphere_collision(pos, vel, radii(pid1), masses(pid1), pos_2, vel_2, radii(pid2), &
                                                     & masses(pid2), force)
            end select

            ! equal and opposite force on the second particle
            do l = 1, num_dims
                $:GPU_ATOMIC(atomic='update')
                accelerations(pid1, l) = accelerations(pid1, l) + force(l)/(masses(pid1) + added_masses(pid1))
                $:GPU_ATOMIC(atomic='update')
                accelerations(pid2, l) = accelerations(pid2, l) - force(l)/(masses(pid2) + added_masses(pid2))
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_EL_collisions

    !> Linear spring-dashpot with Coulomb friction (Cundall & Strack, Geotechnique 29(1), 1979) between two overlapping spheres or
    !! circles. Body 2 receives -force_1. Rotation and torques are included only when they are passed.
    subroutine s_perform_soft_sphere_collision(pos_1, vel_1, radius_1, mass_1, pos_2, vel_2, radius_2, mass_2, force_1, &
        & angular_vel_1, angular_vel_2, torque_1, torque_2)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos_1, vel_1, pos_2, vel_2
        real(wp), intent(in) :: radius_1, mass_1, radius_2, mass_2
        real(wp), dimension(3), intent(out) :: force_1
        real(wp), dimension(3), intent(in), optional :: angular_vel_1, angular_vel_2
        real(wp), dimension(3), intent(out), optional :: torque_1, torque_2
        real(wp) :: overlap_distance
        real(wp), dimension(3) :: normal_vector, contact_point
        real(wp), dimension(3) :: normal_velocity, tangential_vector, normal_force, tangential_force, radial_vector, &
             & rotation_velocity, vel1, vel2
        real(wp) :: k, eta, effective_mass  ! the spring stiffness and damping coefficient and mass of a specific interaction

        force_1 = 0._wp
        if (present(torque_1)) torque_1 = 0._wp
        if (present(torque_2)) torque_2 = 0._wp

        ! the caller guarantees overlap_distance > 0
        normal_vector = pos_2 - pos_1
        overlap_distance = radius_1 + radius_2 - norm2(normal_vector)
        normal_vector = normal_vector/norm2(normal_vector)

        ! ensure the local rank owns that collision before proceeding; owning by centroid drops or doubles pairs split across ranks
        contact_point = pos_1 + normal_vector*(radius_1 - 0.5_wp*overlap_distance)
        if (.not. f_local_rank_owns_location(contact_point, glb_bounds)) return

        ! compute constants of the collision
        effective_mass = 1.0_wp/((1.0_wp/mass_1) + (1._wp/mass_2))
        k = spring_stiffness*effective_mass
        eta = damping_parameter*effective_mass

        ! velocity of each surface at the contact point
        vel1 = vel_1
        vel2 = vel_2
        if (present(torque_1)) then
            radial_vector = normal_vector*(radius_1 - 0.5_wp*overlap_distance)
            call s_cross_product(angular_vel_1, radial_vector, rotation_velocity)
            vel1 = vel_1 + rotation_velocity
            radial_vector = normal_vector*(-1.0_wp)*(radius_2 - 0.5_wp*overlap_distance)
            call s_cross_product(angular_vel_2, radial_vector, rotation_velocity)
            vel2 = vel_2 + rotation_velocity
        end if

        normal_velocity = dot_product(vel1 - vel2, normal_vector)*normal_vector
        tangential_vector = (vel1 - vel2) - normal_velocity
        if (.not. f_approx_equal(norm2(tangential_vector), 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)

        ! compute force and torque
        normal_force = -k*overlap_distance*normal_vector - eta*normal_velocity
        tangential_force = -coefficient_of_friction*norm2(normal_force)*tangential_vector
        force_1 = normal_force + tangential_force

        ! body 2's lever arm and tangential force both flip, so its torque keeps the sign
        if (present(torque_1)) then
            call s_cross_product(normal_vector*radius_1, tangential_force, torque_1)
            torque_2 = torque_1*radius_2/radius_1
        end if

    end subroutine s_perform_soft_sphere_collision

    !> Soft-sphere contact between one sphere or circle and every domain wall it overlaps
    subroutine s_soft_sphere_wall_force(pos, vel, radius, mass, wall_overlaps, force, angular_vel, torque)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos, vel
        real(wp), intent(in) :: radius, mass
        real(wp), dimension(6), intent(in) :: wall_overlaps
        real(wp), dimension(3), intent(out) :: force
        real(wp), dimension(3), intent(in), optional :: angular_vel
        real(wp), dimension(3), intent(out), optional :: torque
        integer :: i
        real(wp), dimension(3) :: normal_force, tangential_force, normal_vector, normal_velocity, tangential_vector, wall_torque, &
             & radial_vector, rotation_velocity, velocity
        real(wp) :: k, eta  ! the spring stiffness and damping coefficient for a specific

        force = 0._wp
        if (present(torque)) torque = 0._wp

        ! ensure the local rank owns that collision before proceeding
        if (.not. f_local_rank_owns_location(pos, glb_bounds)) return

        k = spring_stiffness*mass
        eta = damping_parameter*mass

        do i = 1, num_dims*2
            ! only compute force contributions if there was an overlap
            if (f_approx_equal(wall_overlaps(i), 0._wp)) cycle

            select case (i)
            case (1)  ! x domain left
                normal_vector = [-1._wp, 0._wp, 0._wp]
            case (2)  ! x domain right
                normal_vector = [1._wp, 0._wp, 0._wp]
            case (3)  ! y domain bottom
                normal_vector = [0._wp, -1._wp, 0._wp]
            case (4)  ! y domain top
                normal_vector = [0._wp, 1._wp, 0._wp]
            case (5)  ! z domain back
                normal_vector = [0._wp, 0._wp, -1._wp]
            case (6)  ! z domain front
                normal_vector = [0._wp, 0._wp, 1._wp]
            end select

            velocity = vel
            if (present(torque)) then
                ! get the vector that points from the centroid to the point of collision
                radial_vector = normal_vector*(radius - wall_overlaps(i))
                ! convert the angular velocity to linear velocity
                call s_cross_product(angular_vel, radial_vector, rotation_velocity)
                velocity = velocity + rotation_velocity
            end if

            ! standard soft-sphere collision  with the wall
            normal_velocity = dot_product(velocity, normal_vector)*normal_vector
            tangential_vector = velocity - normal_velocity
            if (.not. f_approx_equal(norm2(tangential_vector), &
                & 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)
            normal_force = -k*wall_overlaps(i)*normal_vector - eta*normal_velocity
            tangential_force = -coefficient_of_friction*norm2(normal_force)*tangential_vector
            force = force + (normal_force + tangential_force)

            ! compute the torque if that is requested
            if (present(torque)) then
                call s_cross_product(normal_vector*radius, tangential_force, wall_torque)
                torque = torque + wall_torque
            end if
        end do

    end subroutine s_soft_sphere_wall_force

    !> uses ghost-point/image-point information to determine if it is possible if two IBs are colliding, effectively an optimized
    !! nearest neighbor search
    subroutine s_detect_ib_collisions(gps, ib_markers, num_gps, num_considered_collisions)

        type(ghost_point), dimension(num_gps), intent(in) :: gps
        type(integer_field), intent(in)                   :: ib_markers
        integer, intent(in)                               :: num_gps
        integer, intent(out)                              :: num_considered_collisions
        integer                                           :: i, j, k, z_bound, ii, jj, kk
        integer, dimension(2)                             :: decoded_pairs
        integer                                           :: gp_idx, gp_patch_id, neighbor_patch_id, local_idx
        integer                                           :: pair_idx, out_idx
        logical                                           :: already_found

        ! Temporary array to hold all detected pairs (with potential duplicates)
        integer, dimension(num_gps, 2) :: raw_pairs
        integer                        :: num_raw, local_num_raw

        num_raw = 0
        z_bound = 0; if (num_dims == 3) z_bound = 2

        $:GPU_PARALLEL_LOOP(private='[gp_idx, gp_patch_id, neighbor_patch_id, local_num_raw, i, j, k, ii, jj, kk]', &
                            & copy='[raw_pairs, num_raw]', copyin='[z_bound]')
        do gp_idx = 1, num_gps
            i = gps(gp_idx)%loc(1)
            j = gps(gp_idx)%loc(2)
            k = 0; if (num_dims == 3) k = gps(gp_idx)%loc(3)
            gp_patch_id = ib_markers%sf(i, j, k)

            ! search in a cube around the BG for Ib markers belonging to another patch
            neighbor_search: do ii = i - 2, i + 2
                do jj = j - 2, j + 2
                    do kk = k - z_bound, k + z_bound
                        neighbor_patch_id = ib_markers%sf(ii, jj, kk)

                        ! Any neighbor of a different patch is a candidate pair. Both patches record it: the rank that owns the
                        ! contact point may hold interior ghost points of only one of them, so one-sided detection can leave
                        ! that rank blind to a contact that sits within a cell of its boundary. The host pass below sorts the
                        ! pair and drops duplicates.
                        if (neighbor_patch_id /= 0 .and. neighbor_patch_id /= gp_patch_id) then
                            $:GPU_ATOMIC(atomic='capture')
                            num_raw = num_raw + 1
                            local_num_raw = num_raw
                            $:END_GPU_ATOMIC_CAPTURE()

                            raw_pairs(local_num_raw, 1) = gp_patch_id
                            raw_pairs(local_num_raw, 2) = neighbor_patch_id
                            exit neighbor_search
                        end if
                    end do
                end do
            end do neighbor_search
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Coalesce collisions unique pairs
        num_considered_collisions = 0
        collision_lookup = 0
        ! for each pair found in the raw collection
        do pair_idx = 1, num_raw
            already_found = .false.

            ! get the decoded pairs for checking if they exist, using ii,jj,kk as dummy indices
            call s_decode_patch_periodicity(raw_pairs(pair_idx, 1), decoded_pairs(1), ii, jj, kk)
            call s_decode_patch_periodicity(raw_pairs(pair_idx, 2), decoded_pairs(2), ii, jj, kk)
            call s_get_neighborhood_idx(decoded_pairs(1), local_idx)
            decoded_pairs(1) = local_idx
            call s_get_neighborhood_idx(decoded_pairs(2), local_idx)
            decoded_pairs(2) = local_idx

            ! skip self-collisions (an IB cannot collide with its own periodic image)
            if (decoded_pairs(1) == decoded_pairs(2)) cycle

            ! need to swap to guarantee the smaller decoded marker value is in index 1 and prevent double-counting
            if (decoded_pairs(2) < decoded_pairs(1)) then
                decoded_pairs(1) = decoded_pairs(1) + decoded_pairs(2)
                decoded_pairs(2) = decoded_pairs(1) - decoded_pairs(2)
                decoded_pairs(1) = decoded_pairs(1) - decoded_pairs(2)
                raw_pairs(pair_idx, 1) = raw_pairs(pair_idx, 1) + raw_pairs(pair_idx, 2)
                raw_pairs(pair_idx, 2) = raw_pairs(pair_idx, 1) - raw_pairs(pair_idx, 2)
                raw_pairs(pair_idx, 1) = raw_pairs(pair_idx, 1) - raw_pairs(pair_idx, 2)
            end if

            ! check if it is already in the list
            do out_idx = 1, num_considered_collisions
                if (collision_lookup(out_idx, 1) == decoded_pairs(1) .and. collision_lookup(out_idx, 2) == decoded_pairs(2)) then
                    already_found = .true.
                    exit
                end if
            end do

            ! and if it is not, append it to the list of pairs
            if (.not. already_found) then
                num_considered_collisions = num_considered_collisions + 1
                @:PROHIBIT(num_considered_collisions > size(collision_lookup, 1) , &
                           & "More collisions detected than memory to hold them. Consider increasing the size of the collision_lookup array")

                collision_lookup(num_considered_collisions, 1) = decoded_pairs(1)
                collision_lookup(num_considered_collisions, 2) = decoded_pairs(2)
                collision_lookup(num_considered_collisions, 3) = raw_pairs(pair_idx, 1)
                collision_lookup(num_considered_collisions, 4) = raw_pairs(pair_idx, 2)
            end if
        end do
        $:GPU_UPDATE(device='[collision_lookup]')

    end subroutine s_detect_ib_collisions

    !> Naive O(N^2) search over this rank's EL particles; appends every overlapping pair (i < j) to collision_lookup
    subroutine s_detect_el_collisions(num_particles, positions, radii, num_considered_collisions)

        integer, intent(in)                    :: num_particles
        real(wp), dimension(:,:,:), intent(in) :: positions  !< (particle, component, RK level); level 2 is current
        real(wp), dimension(:), intent(in)     :: radii
        integer, intent(out)                   :: num_considered_collisions
        real(wp)                               :: distance_squared
        integer                                :: i, j, pair_idx, num_found, max_pairs

        num_found = 0
        max_pairs = size(collision_lookup, 1)

        ! one thread per (i, j) so the GPU is filled; a per-i thread walking all j used under 1% of an A100
        $:GPU_PARALLEL_LOOP(private='[i, j, distance_squared, pair_idx]', copy='[num_found]', copyin='[max_pairs]', collapse=2)
        do i = 1, num_particles
            do j = 1, num_particles
                if (j <= i) cycle

                ! squared distances avoid norm2; d < r_i + r_j is the overlap test
                distance_squared = (positions(j, 1, 2) - positions(i, 1, 2))**2 + (positions(j, 2, 2) - positions(i, 2, 2))**2
                if (num_dims == 3) distance_squared = distance_squared + (positions(j, 3, 2) - positions(i, 3, 2))**2
                if (distance_squared >= (radii(i) + radii(j))**2) cycle

                $:GPU_ATOMIC(atomic='capture')
                num_found = num_found + 1
                pair_idx = num_found
                $:END_GPU_ATOMIC_CAPTURE()

                if (pair_idx > max_pairs) cycle
                collision_lookup(pair_idx, 1) = i
                collision_lookup(pair_idx, 2) = j
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        @:PROHIBIT(num_found > max_pairs, &
                   & "More collisions detected than memory to hold them. Consider increasing the size of the collision_lookup array")
        num_considered_collisions = num_found

    end subroutine s_detect_el_collisions

    !> Distance a sphere or circle extends past each slip or no-slip domain wall (1/2: x beg/end, 3/4: y, 5/6: z); 0 if none
    subroutine s_detect_wall_overlaps(centroid, radius, wall_overlaps)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)  :: centroid
        real(wp), intent(in)                :: radius
        real(wp), dimension(6), intent(out) :: wall_overlaps

        wall_overlaps = 0._wp
        #:for X, DIR, IDX in [('x', 1, 1), ('y', 2, 3), ('z', 3, 5)]
            if (num_dims >= ${DIR}$) then
                if (ib_bc_${X}$%beg == BC_SLIP_WALL .or. ib_bc_${X}$%beg == BC_NO_SLIP_WALL) then
                    wall_overlaps(${IDX}$) = max(0._wp, glb_bounds(${DIR}$)%beg - (centroid(${DIR}$) - radius))
                end if
                if (ib_bc_${X}$%end == BC_SLIP_WALL .or. ib_bc_${X}$%end == BC_NO_SLIP_WALL) then
                    wall_overlaps(${IDX}$ + 1) = max(0._wp, centroid(${DIR}$) + radius - glb_bounds(${DIR}$)%end)
                end if
            end if
        #:endfor

    end subroutine s_detect_wall_overlaps

    !> @brief function checks if this local MPI processor owns this specific collision
    function f_neighborhood_ranks_own_location(location) result(owns_collision)

        real(wp), dimension(3), intent(in) :: location
        logical                            :: owns_collision, periodic_owner
        real(wp)                           :: temp_neighbor_domain
        integer                            :: i

        owns_collision = .true.

#ifdef MFC_MPI
        if (num_procs > 2) then
            ! catch the edge case where th collision lies just outside the computational domain
            owns_collision = .true.
            #:for X, ID in [('x', 1), ('y', 2,), ('z', 3,)]
                if (num_dims >= ${ID}$) then
                    if (ib_bc_${X}$%beg == BC_PERIODIC .and. neighbor_domain_${X}$%beg >= neighbor_domain_${X}$%end) then
                        ! project right side to the left
                        temp_neighbor_domain = neighbor_domain_${X}$%end + (glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg)
                        periodic_owner = neighbor_domain_${X}$%beg <= location(${ID}$) .and. location(${ID}$) < temp_neighbor_domain
                        ! project the left side to the right
                        temp_neighbor_domain = neighbor_domain_${X}$%beg - (glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg)
                        periodic_owner = periodic_owner .or. (temp_neighbor_domain <= location(${ID}$) .and. location(${ID}$) &
                                                              & < neighbor_domain_${X}$%end)

                        owns_collision = owns_collision .and. periodic_owner
                    else
                        owns_collision = owns_collision .and. neighbor_domain_${X}$%beg <= location(${ID}$) .and. location(${ID}$) &
                            & < neighbor_domain_${X}$%end
                    end if
                end if
            #:endfor
        end if
#endif

    end function f_neighborhood_ranks_own_location

    subroutine s_finalize_collisions_module()

        if (particles_lagrange .or. ib) then
            @:DEALLOCATE(collision_lookup)
        end if
        if (ib) then
            @:DEALLOCATE(wall_overlap_distances)
        end if

    end subroutine s_finalize_collisions_module

end module m_collisions
