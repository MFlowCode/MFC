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

    !> One bit per unordered pair (a < b) of this rank's IBs (local + neighborhood indices, a triangular index over
    !! pair_bits_capacity): the ghost-point candidate pass claims each pair with an atomic OR, so only the first thread to find a
    !! pair records it and duplicates are dropped on the GPU. Bits are cleared after the forces are applied, so the set is empty
    !! between calls. Sized from num_ibs and regrown with it; memory is quadratic in the body count (1e4 bodies take 6 MB).
    integer, allocatable, dimension(:) :: pair_bits
    integer                            :: pair_bits_capacity  !< Number of IBs the bitset is indexed for
    $:GPU_DECLARE(create='[pair_bits]')

contains

    !> @brief the point that decides which rank owns a body pair: the contact point built from the lower global id's own centroid
    !! and the minimum-image offset to the other, so every rank holding the pair computes the same bits whichever periodic images it
    !! found the pair through
    subroutine s_pair_owner_point(pa, pb, point)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(ib_patch_parameters), intent(in) :: pa, pb
        real(wp), intent(out)                 :: point(3)
        real(wp)                              :: base(3), rel(3), rb, span

        if (pa%gbl_patch_id < pb%gbl_patch_id) then
            base = [pa%x_centroid, pa%y_centroid, 0._wp]; rel = [pb%x_centroid, pb%y_centroid, 0._wp] - base; rb = pa%radius
            if (num_dims == 3) then
                base(3) = pa%z_centroid; rel(3) = pb%z_centroid - base(3)
            end if
        else
            base = [pb%x_centroid, pb%y_centroid, 0._wp]; rel = [pa%x_centroid, pa%y_centroid, 0._wp] - base; rb = pb%radius
            if (num_dims == 3) then
                base(3) = pb%z_centroid; rel(3) = pa%z_centroid - base(3)
            end if
        end if
        #:for X, D in [('x', 1), ('y', 2), ('z', 3)]
            if (num_dims >= ${D}$ .and. ib_bc_${X}$%beg == BC_PERIODIC) then
                span = glb_bounds(${D}$)%end - glb_bounds(${D}$)%beg
                rel(${D}$) = rel(${D}$) - span*anint(rel(${D}$)/span)
            end if
        #:endfor
        point = base + rel/norm2(rel)*(rb - 0.5_wp*(pa%radius + pb%radius - norm2(rel)))

    end subroutine s_pair_owner_point

    !> @brief inverse effective tangential mass of a sphere or circle at its contact point: 1/m + r^2/I
    function f_inv_tangential_mass(patch) result(inv_mass)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(ib_patch_parameters), intent(in) :: patch
        real(wp)                              :: inv_mass

        inv_mass = 1._wp/patch%mass
        if (patch%moment > 0._wp) inv_mass = inv_mass + patch%radius**2/patch%moment

    end function f_inv_tangential_mass

    !> @brief scales the slip velocity into the friction direction, -F_t = mu*|F_n|*v_t/max(|v_t|, eps). With ib_friction_slip_eps >
    !! 0, eps is floored at the stick limit mu*|F_n|*dt/m_t, which keeps the explicit update stable. With 0 it is v_t/|v_t|, left
    !! unscaled at exactly zero slip.
    subroutine s_friction_direction(tangential_vector, normal_force_mag, inv_mass_t)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(inout) :: tangential_vector(3)
        real(wp), intent(in)    :: normal_force_mag, inv_mass_t
        real(wp)                :: slip

        slip = norm2(tangential_vector)
        if (ib_friction_slip_eps > 0._wp) then
            tangential_vector = tangential_vector/max(slip, ib_friction_slip_eps, &
                & coefficient_of_friction*normal_force_mag*dt*inv_mass_t)
        else if (.not. f_approx_equal(slip, 0._wp)) then
            tangential_vector = tangential_vector/slip
        end if

    end subroutine s_friction_direction

    subroutine s_initialize_collisions_module()

        real(wp) :: e

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness]')

        ! particles are loaded after this, so the EL part of the lookup grows on the first s_apply_EL_collisions
        call s_resize_collision_lookup(0)

        if (ib) then
            @:ALLOCATE(wall_overlap_distances(num_local_ibs_max*27, 6))
            wall_overlap_distances = 0
            $:GPU_UPDATE(device='[wall_overlap_distances]')
            call s_resize_pair_bits(num_ibs)
        end if

        collisions_active = .false.
        $:GPU_UPDATE(device='[coefficient_of_friction, ib_friction_slip_eps]')
        $:GPU_UPDATE(device='[ib_bc_x%beg, ib_bc_x%end, ib_bc_y%beg, ib_bc_y%end, ib_bc_z%beg, ib_bc_z%end]')

    end subroutine s_initialize_collisions_module

    !> Make collision_lookup hold the IB pair bound (when ib) and twice num_pairs EL pairs; it only grows. Contents are rebuilt by
    !! every detection, so they are not kept.
    subroutine s_resize_collision_lookup(num_pairs)

        integer, intent(in) :: num_pairs
        integer             :: capacity

        capacity = max(2*num_pairs, 1)
        if (ib) capacity = max(capacity, num_local_ibs_max*27*8)
        if (allocated(collision_lookup)) then
            if (size(collision_lookup, 1) >= capacity) return
            @:DEALLOCATE(collision_lookup)
        end if
        @:ALLOCATE(collision_lookup(capacity, 4))
        collision_lookup = 0
        $:GPU_UPDATE(device='[collision_lookup]')

    end subroutine s_resize_collision_lookup

    !> Size the pair bitset for num_bodies IBs, doubled so a growing body count does not reallocate every step (empty on exit)
    subroutine s_resize_pair_bits(num_bodies)

        integer, intent(in) :: num_bodies

        if (allocated(pair_bits)) then
            @:DEALLOCATE(pair_bits)
        end if

        pair_bits_capacity = min(max(2*num_bodies, 2), num_local_ibs_max)
        @:ALLOCATE(pair_bits(f_pair_bit_index(pair_bits_capacity - 1, pair_bits_capacity, pair_bits_capacity)/32 + 1))
        pair_bits = 0
        $:GPU_UPDATE(device='[pair_bits]')

    end subroutine s_resize_pair_bits

    subroutine s_apply_IB_collisions(ghost_points, num_gps, ib_markers, forces, torques)

        type(ghost_point), dimension(:), intent(in)    :: ghost_points
        integer, intent(in)                            :: num_gps
        type(integer_field), intent(in)                :: ib_markers
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: num_considered_collisions, i, l, pid1, pid2, decoded_id
        integer, dimension(3)                          :: periodic_image_1, periodic_image_2
        logical                                        :: any_wall_collision
        real(wp), dimension(3)                         :: pos, vel, angular_vel, force, torque
        real(wp), dimension(3)                         :: pos_2, vel_2, angular_vel_2, torque_2, contact_point
        real(wp), dimension(6)                         :: wall_overlaps
        real(wp)                                       :: radius, mass, max_overlap, distance, overlap_distance

        ! return if no collisions

        if (collision_model == 0) return

        ! IBs migrating onto this rank can outgrow the bitset
        if (num_ibs > pair_bits_capacity) call s_resize_pair_bits(num_ibs)

        ! iterating over particles, determine if they are colliding with the wall
        max_overlap = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, pos, wall_overlaps]', reduction='[[max_overlap]]', reductionOp='[max]')
        do i = 1, num_ibs
            pos = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
            if (num_dims == 3) pos(3) = patch_ib(i)%z_centroid
            call s_detect_wall_overlaps(pos, patch_ib(i)%radius, .true., wall_overlaps)
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
                call s_soft_sphere_wall_force(pos, vel, radius, mass, f_inv_tangential_mass(patch_ib(i)), wall_overlaps, force, &
                                              & angular_vel, torque)
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
                            & pos_2, vel_2, angular_vel_2, force, torque, torque_2, distance, overlap_distance, contact_point]', &
                            & copy='[forces, torques]')
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
            distance = norm2(pos_2 - pos)
            overlap_distance = patch_ib(pid1)%radius + patch_ib(pid2)%radius - distance
            if (overlap_distance <= 0._wp .or. distance < sgm_eps) cycle

            ! IB forces are summed across ranks afterwards, so only the rank holding the contact point applies the pair; the point
            ! is computed the same way on every rank (pid1 and the images are rank-local)
            call s_pair_owner_point(patch_ib(pid1), patch_ib(pid2), contact_point)
            if (.not. f_local_rank_owns_location(contact_point, glb_bounds)) cycle

            vel = patch_ib(pid1)%vel
            vel_2 = patch_ib(pid2)%vel
            angular_vel = patch_ib(pid1)%angular_vel
            angular_vel_2 = patch_ib(pid2)%angular_vel

            force = 0._wp; torque = 0._wp; torque_2 = 0._wp
            select case (collision_model)
            case (1)  ! soft sphere model
                call s_perform_soft_sphere_collision(pos, vel, patch_ib(pid1)%radius, patch_ib(pid1)%mass, pos_2, vel_2, &
                                                     & patch_ib(pid2)%radius, patch_ib(pid2)%mass, &
                                                     & f_inv_tangential_mass(patch_ib(pid1)) &
                                                     & + f_inv_tangential_mass(patch_ib(pid2)), force, angular_vel, &
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

        call s_clear_pair_bits(num_considered_collisions)

    end subroutine s_apply_IB_collisions

    !> Add the wall and particle-particle contact accelerations of the num_particles owned EL particles to accelerations, dividing
    !! by the same inertia as the hydrodynamic force (mass plus added mass). The num_halo copies stored after them are partners
    !! only: the rank owning each copy computes the same pair and applies the opposite force itself.
    subroutine s_apply_EL_collisions(num_particles, num_halo, positions, velocities, radii, masses, added_masses, accelerations)

        integer, intent(in)                     :: num_particles, num_halo
        real(wp), dimension(:,:,:), intent(in)  :: positions, velocities  !< (particle, component, RK level); level 2 is current
        real(wp), dimension(:), intent(in)      :: radii, masses, added_masses
        real(wp), dimension(:,:), intent(inout) :: accelerations
        integer                                 :: num_considered_collisions, i, l, pid1, pid2
        logical                                 :: any_wall_collision
        real(wp), dimension(3)                  :: pos, vel, pos_2, vel_2, force
        real(wp), dimension(6)                  :: wall_overlaps
        real(wp)                                :: max_overlap

        if (collision_model == 0) return

        ! room for 12 contacts per owned particle (3D close packing), sized from the local count
        if (12*num_particles > size(collision_lookup, 1)) call s_resize_collision_lookup(12*num_particles)

        ! apply wall collisions; one thread per particle, so no atomics
        max_overlap = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, l, pos, vel, wall_overlaps, force]', reduction='[[max_overlap]]', reductionOp='[max]')
        do i = 1, num_particles
            pos = positions(i,:,2)
            if (num_dims == 2) pos(3) = 0._wp
            vel = velocities(i,:,2)
            ! reflective faces get mirror-image halo particles instead (s_build_halo_particles)
            call s_detect_wall_overlaps(pos, radii(i), .false., wall_overlaps)
            max_overlap = max(max_overlap, maxval(wall_overlaps))

            force = 0._wp
            select case (collision_model)
            case (1)  ! soft sphere model
                call s_soft_sphere_wall_force(pos, vel, radii(i), masses(i), 1._wp/masses(i), wall_overlaps, force)
            end select

            do l = 1, num_dims
                accelerations(i, l) = accelerations(i, l) + force(l)/(masses(i) + added_masses(i))
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        any_wall_collision = max_overlap > 0._wp

        ! naive all-to-all search; every pair in the lookup overlaps
        call s_detect_el_collisions(num_particles, num_halo, positions, radii, num_considered_collisions)
        if (num_considered_collisions > size(collision_lookup, 1)) then
            ! more contacts than reserved (e.g. deep overlaps or unequal radii): grow to the count found and detect again
            call s_resize_collision_lookup(num_considered_collisions)
            call s_detect_el_collisions(num_particles, num_halo, positions, radii, num_considered_collisions)
        end if

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
                                                     & masses(pid2), 1._wp/masses(pid1) + 1._wp/masses(pid2), force)
            end select

            ! equal and opposite force on the second particle; a halo copy (pid2 > num_particles) has no added mass here and its
            ! owner applies the force itself
            do l = 1, num_dims
                $:GPU_ATOMIC(atomic='update')
                accelerations(pid1, l) = accelerations(pid1, l) + force(l)/(masses(pid1) + added_masses(pid1))
            end do
            if (pid2 <= num_particles) then
                do l = 1, num_dims
                    $:GPU_ATOMIC(atomic='update')
                    accelerations(pid2, l) = accelerations(pid2, l) - force(l)/(masses(pid2) + added_masses(pid2))
                end do
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_EL_collisions

    !> Linear spring-dashpot with Coulomb friction (Cundall & Strack, Geotechnique 29(1), 1979) between two overlapping spheres or
    !! circles. Body 2 receives -force_1. Rotation and torques are included only when they are passed.
    subroutine s_perform_soft_sphere_collision(pos_1, vel_1, radius_1, mass_1, pos_2, vel_2, radius_2, mass_2, inv_mass_t, &
        & force_1, angular_vel_1, angular_vel_2, torque_1, torque_2)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos_1, vel_1, pos_2, vel_2
        real(wp), intent(in) :: radius_1, mass_1, radius_2, mass_2
        real(wp), intent(in) :: inv_mass_t  !< inverse effective tangential mass of the contact, for the friction stick limit
        real(wp), dimension(3), intent(out) :: force_1
        real(wp), dimension(3), intent(in), optional :: angular_vel_1, angular_vel_2
        real(wp), dimension(3), intent(out), optional :: torque_1, torque_2
        real(wp) :: overlap_distance, center_distance
        real(wp), dimension(3) :: normal_vector
        real(wp), dimension(3) :: normal_velocity, tangential_vector, normal_force, tangential_force, radial_vector, &
             & rotation_velocity, vel1, vel2
        real(wp) :: k, eta, effective_mass  ! the spring stiffness and damping coefficient and mass of a specific interaction

        force_1 = 0._wp
        if (present(torque_1)) torque_1 = 0._wp
        if (present(torque_2)) torque_2 = 0._wp

        ! the caller guarantees overlap_distance > 0; coincident centers have no contact direction
        normal_vector = pos_2 - pos_1
        center_distance = norm2(normal_vector)
        if (center_distance < sgm_eps) return
        overlap_distance = radius_1 + radius_2 - center_distance
        normal_vector = normal_vector/center_distance

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

        ! compute force and torque
        normal_force = -k*overlap_distance*normal_vector - eta*normal_velocity
        call s_friction_direction(tangential_vector, norm2(normal_force), inv_mass_t)
        tangential_force = -coefficient_of_friction*norm2(normal_force)*tangential_vector
        force_1 = normal_force + tangential_force

        ! body 2's lever arm and tangential force both flip, so its torque keeps the sign
        if (present(torque_1)) then
            call s_cross_product(normal_vector*radius_1, tangential_force, torque_1)
            torque_2 = torque_1*radius_2/radius_1
        end if

    end subroutine s_perform_soft_sphere_collision

    !> Soft-sphere contact between one sphere or circle and every domain wall it overlaps. A reflective face acts as the
    !! mirror-image body, whose contact gives the same normal force and no tangential slip, so it is a frictionless wall.
    subroutine s_soft_sphere_wall_force(pos, vel, radius, mass, inv_mass_t, wall_overlaps, force, angular_vel, torque)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos, vel
        real(wp), intent(in) :: radius, mass
        real(wp), intent(in) :: inv_mass_t  !< inverse effective tangential mass, for the friction stick limit
        real(wp), dimension(6), intent(in) :: wall_overlaps
        real(wp), dimension(3), intent(out) :: force
        real(wp), dimension(3), intent(in), optional :: angular_vel
        real(wp), dimension(3), intent(out), optional :: torque
        integer :: i, face_bc
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
                face_bc = ib_bc_x%beg
            case (2)  ! x domain right
                normal_vector = [1._wp, 0._wp, 0._wp]
                face_bc = ib_bc_x%end
            case (3)  ! y domain bottom
                normal_vector = [0._wp, -1._wp, 0._wp]
                face_bc = ib_bc_y%beg
            case (4)  ! y domain top
                normal_vector = [0._wp, 1._wp, 0._wp]
                face_bc = ib_bc_y%end
            case (5)  ! z domain back
                normal_vector = [0._wp, 0._wp, -1._wp]
                face_bc = ib_bc_z%beg
            case (6)  ! z domain front
                normal_vector = [0._wp, 0._wp, 1._wp]
                face_bc = ib_bc_z%end
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
            normal_force = -k*wall_overlaps(i)*normal_vector - eta*normal_velocity
            tangential_force = 0._wp
            if (face_bc /= BC_REFLECTIVE) then
                call s_friction_direction(tangential_vector, norm2(normal_force), inv_mass_t)
                tangential_force = -coefficient_of_friction*norm2(normal_force)*tangential_vector
            end if
            force = force + (normal_force + tangential_force)

            ! compute the torque if that is requested
            if (present(torque)) then
                call s_cross_product(normal_vector*radius, tangential_force, wall_torque)
                torque = torque + wall_torque
            end if
        end do

    end subroutine s_soft_sphere_wall_force

    !> Triangular index (0-based) of the unordered pair a < b, both in 1 to num_bodies
    pure function f_pair_bit_index(a, b, num_bodies) result(idx)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: a, b, num_bodies
        integer(kind=8)     :: idx

        idx = int(a - 1, 8)*int(2*num_bodies - a, 8)/2 + int(b - a - 1, 8)

    end function f_pair_bit_index

    !> Uses ghost-point information to find IB pairs that may be colliding, effectively an optimized nearest-neighbor search: every
    !! ghost point with another patch's marker within 2 cells gives a candidate pair. Both patches record it: the rank that owns the
    !! contact point may hold interior ghost points of only one of them, so one-sided detection can leave that rank blind to a
    !! contact within a cell of its boundary. The many duplicate candidates are dropped on the GPU: each candidate is mapped to its
    !! sorted neighborhood pair (a, b) and claims that pair's bit in pair_bits with an atomic OR, and only the thread that flips the
    !! bit appends the pair. Bits are cleared by s_clear_pair_bits after the forces. A pair's periodic image is taken from the
    !! claiming candidate; two different images of one pair would need a domain about two particles across.
    subroutine s_detect_ib_collisions(gps, ib_markers, num_gps, num_considered_collisions)

        type(ghost_point), dimension(num_gps), intent(in) :: gps
        type(integer_field), intent(in)                   :: ib_markers
        integer, intent(in)                               :: num_gps
        integer, intent(out)                              :: num_considered_collisions
        integer                                           :: i, j, k, ii, jj, kk, z_bound, gp_idx, gp_patch_id, neighbor_patch_id
        integer                                           :: a, b, gbl_a, gbl_b, enc_a, enc_b, xp, yp, zp, bit, mask, old_word, slot
        integer                                           :: num_found, max_pairs, num_bodies
        integer(kind=8)                                   :: idx, word

        num_found = 0
        max_pairs = size(collision_lookup, 1)
        num_bodies = pair_bits_capacity
        z_bound = 0; if (num_dims == 3) z_bound = 2

        $:GPU_PARALLEL_LOOP(private='[gp_idx, i, j, k, ii, jj, kk, gp_patch_id, neighbor_patch_id, a, b, gbl_a, gbl_b, enc_a, &
                            & enc_b, xp, yp, zp, idx, word, bit, mask, old_word, slot]', copy='[num_found]', copyin='[z_bound, &
                            & max_pairs, num_bodies]')
        do gp_idx = 1, num_gps
            i = gps(gp_idx)%loc(1)
            j = gps(gp_idx)%loc(2)
            k = 0; if (num_dims == 3) k = gps(gp_idx)%loc(3)
            gp_patch_id = ib_markers%sf(i, j, k)

            ! first marker of another patch in the cube around the ghost point
            neighbor_patch_id = 0
            neighbor_search: do ii = i - 2, i + 2
                do jj = j - 2, j + 2
                    do kk = k - z_bound, k + z_bound
                        if (ib_markers%sf(ii, jj, kk) /= 0 .and. ib_markers%sf(ii, jj, kk) /= gp_patch_id) then
                            neighbor_patch_id = ib_markers%sf(ii, jj, kk)
                            exit neighbor_search
                        end if
                    end do
                end do
            end do neighbor_search
            if (neighbor_patch_id == 0) cycle

            ! sorted neighborhood indices, the encoded ids (which carry the periodic image) following their patch
            call s_decode_patch_periodicity(gp_patch_id, gbl_a, xp, yp, zp)
            call s_get_neighborhood_idx(gbl_a, a)
            call s_decode_patch_periodicity(neighbor_patch_id, gbl_b, xp, yp, zp)
            call s_get_neighborhood_idx(gbl_b, b)
            if (a <= 0 .or. b <= 0 .or. a == b) cycle
            enc_a = gp_patch_id
            enc_b = neighbor_patch_id
            if (b < a) then
                slot = a; a = b; b = slot
                slot = enc_a; enc_a = enc_b; enc_b = slot
            end if

            ! claim the pair: the thread whose OR flips the bit is the first to find it
            idx = f_pair_bit_index(a, b, num_bodies)
            word = idx/32 + 1
            bit = int(mod(idx, 32_8))
            mask = ishft(1, bit)
            $:GPU_ATOMIC(atomic='capture')
            old_word = pair_bits(word)
            pair_bits(word) = ior(pair_bits(word), mask)
            $:END_GPU_ATOMIC_CAPTURE()
            if (iand(old_word, mask) /= 0) cycle

            $:GPU_ATOMIC(atomic='capture')
            num_found = num_found + 1
            slot = num_found
            $:END_GPU_ATOMIC_CAPTURE()
            if (slot > max_pairs) cycle
            collision_lookup(slot, 1) = a
            collision_lookup(slot, 2) = b
            collision_lookup(slot, 3) = enc_a
            collision_lookup(slot, 4) = enc_b
        end do
        $:END_GPU_PARALLEL_LOOP()

        @:PROHIBIT(num_found > max_pairs, &
                   & "More collisions detected than memory to hold them. Consider increasing the size of the collision_lookup array")
        num_considered_collisions = num_found

    end subroutine s_detect_ib_collisions

    !> Clear the bits the last s_detect_ib_collisions set (every set bit belongs to a recorded pair, so whole words can be zeroed)
    subroutine s_clear_pair_bits(num_considered_collisions)

        integer, intent(in) :: num_considered_collisions
        integer             :: n, num_bodies

        num_bodies = pair_bits_capacity
        $:GPU_PARALLEL_LOOP(private='[n]', copyin='[num_bodies]')
        do n = 1, num_considered_collisions
            pair_bits(f_pair_bit_index(collision_lookup(n, 1), collision_lookup(n, 2), num_bodies)/32 + 1) = 0
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_clear_pair_bits

    !> Naive O(N^2) search pairing each owned EL particle i with every later owned particle or halo copy j; appends every
    !! overlapping pair to collision_lookup. Returns the full count, which exceeds the lookup size when pairs were dropped.
    subroutine s_detect_el_collisions(num_particles, num_halo, positions, radii, num_considered_collisions)

        integer, intent(in)                    :: num_particles, num_halo
        real(wp), dimension(:,:,:), intent(in) :: positions  !< (particle, component, RK level); level 2 is current
        real(wp), dimension(:), intent(in)     :: radii
        integer, intent(out)                   :: num_considered_collisions
        real(wp)                               :: distance_squared
        integer                                :: i, j, pair_idx, num_found, max_pairs, num_stored

        num_found = 0
        max_pairs = size(collision_lookup, 1)
        num_stored = num_particles + num_halo

        ! one thread per (i, j) so the GPU is filled; a per-i thread walking all j used under 1% of an A100
        $:GPU_PARALLEL_LOOP(private='[i, j, distance_squared, pair_idx]', copy='[num_found]', copyin='[max_pairs, num_stored]', &
                            & collapse=2)
        do i = 1, num_particles
            do j = 1, num_stored
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

        num_considered_collisions = num_found

    end subroutine s_detect_el_collisions

    !> Distance a sphere or circle extends past each slip or no-slip domain wall, and each reflective face when reflective_walls
    !! (1/2: x beg/end, 3/4: y, 5/6: z); 0 if none. ib_bc_* are the global faces, so interior rank faces are never walls.
    subroutine s_detect_wall_overlaps(centroid, radius, reflective_walls, wall_overlaps)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)  :: centroid
        real(wp), intent(in)                :: radius
        logical, intent(in)                 :: reflective_walls
        real(wp), dimension(6), intent(out) :: wall_overlaps

        wall_overlaps = 0._wp
        #:for X, DIR, IDX in [('x', 1, 1), ('y', 2, 3), ('z', 3, 5)]
            if (num_dims >= ${DIR}$) then
                ! with cyl_coord, a reflective y beg is the axis, not a mirror plane
                if (ib_bc_${X}$%beg == BC_SLIP_WALL .or. ib_bc_${X}$%beg == BC_NO_SLIP_WALL .or. (reflective_walls &
                    & .and. ib_bc_${X}$%beg == BC_REFLECTIVE .and. .not. (cyl_coord .and. ${DIR}$ == 2))) then
                    wall_overlaps(${IDX}$) = max(0._wp, glb_bounds(${DIR}$)%beg - (centroid(${DIR}$) - radius))
                end if
                if (ib_bc_${X}$%end == BC_SLIP_WALL .or. ib_bc_${X}$%end == BC_NO_SLIP_WALL .or. (reflective_walls &
                    & .and. ib_bc_${X}$%end == BC_REFLECTIVE)) then
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

        if (allocated(collision_lookup)) then
            @:DEALLOCATE(collision_lookup)
        end if
        if (ib) then
            @:DEALLOCATE(pair_bits)
            @:DEALLOCATE(wall_overlap_distances)
        end if

    end subroutine s_finalize_collisions_module

end module m_collisions
