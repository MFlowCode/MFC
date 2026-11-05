!>
!! @file
!! @brief Contains module m_collisions

#:include 'macros.fpp'

!> @brief Contact forces between Lagrangian bodies (IB patches, EL particles) and with walls. Detection is per kind; the force
!! models act on the shared lagrangian_body type, chosen by s_compute_pair_contact and s_compute_wall_contact.
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

    private; public :: s_apply_collision_forces, s_initialize_collisions_module, s_finalize_collisions_module, &
        & f_neighborhood_ranks_own_location, collisions_active, s_compute_pair_contact

    integer, allocatable, dimension(:,:)  :: collision_lookup  !< IB pairs: local ids (1:2) and periodic-encoded ids (3:4)
    real(wp), allocatable, dimension(:,:) :: wall_overlap_distances  !< IB overlap past each of the 6 domain walls
    real(wp)                              :: spring_stiffness, damping_parameter  !< Per unit effective mass
    integer                               :: active_collision_model  !< Device copy of collision_model
    real(wp), dimension(3)                :: kind_friction  !< Coulomb friction coefficient per particle_type_*
    $:GPU_DECLARE(create='[spring_stiffness, damping_parameter, active_collision_model, kind_friction]')
    $:GPU_DECLARE(create='[collision_lookup, wall_overlap_distances]')

    !> true when any contact was detected on this rank since the last adaptive-dt computation
    logical :: collisions_active
    logical :: collisions_initialized = .false.

contains

    !> Set the model constants once; the IB pair and wall buffers are allocated only with IBs
    subroutine s_initialize_collisions_module()

        real(wp) :: e

        if (collisions_initialized) return
        collisions_initialized = .true.

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        active_collision_model = collision_model
        kind_friction = 0._wp
        kind_friction(particle_type_ib) = ib_coefficient_of_friction
        kind_friction(particle_type_el_particle) = particle_coefficient_of_friction
        collisions_active = .false.
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness, active_collision_model, kind_friction]')

        if (ib) then
            @:ALLOCATE(collision_lookup(num_local_ibs_max * 27 * 8, 4))
            @:ALLOCATE(wall_overlap_distances(num_local_ibs_max*27, 6))
            wall_overlap_distances = 0
            $:GPU_UPDATE(device='[wall_overlap_distances]')
        end if

    end subroutine s_initialize_collisions_module

    !> IB driver: detect IB-IB and IB-wall contacts and add their forces and torques to the rank-local IB partial sums
    subroutine s_apply_collision_forces(ghost_points, num_gps, ib_markers, forces, torques)

        type(ghost_point), dimension(:), intent(in)    :: ghost_points
        integer, intent(in)                            :: num_gps
        type(integer_field), intent(in)                :: ib_markers
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: num_considered_collisions
        logical                                        :: any_wall_collision

        if (collision_model == 0) return

        call s_detect_wall_collisions(any_wall_collision)
        call s_detect_ib_collisions(ghost_points, ib_markers, num_gps, num_considered_collisions)

        ! accumulate across RK stages; consumed (and reset) by s_compute_dt once per time step
        collisions_active = collisions_active .or. any_wall_collision .or. (num_considered_collisions > 0)

        call s_apply_ib_wall_contacts(forces, torques)
        call s_apply_ib_pair_contacts(num_considered_collisions, forces, torques)

    end subroutine s_apply_collision_forces

    !> Model tree for a pair: picks the contact model from the two body kinds and collision_model. IB-EL pairs are not yet modeled.
    !! The force is on body_a; torques are nonzero only for IB bodies.
    subroutine s_compute_pair_contact(body_a, shift_a, body_b, shift_b, force_a, torque_a, torque_b, contact_point, in_contact)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(lagrangian_body), intent(in)   :: body_a, body_b
        real(wp), dimension(3), intent(in)  :: shift_a, shift_b  !< Periodic image offsets of each centroid
        real(wp), dimension(3), intent(out) :: force_a, torque_a, torque_b, contact_point
        logical, intent(out)                :: in_contact

        force_a = 0._wp; torque_a = 0._wp; torque_b = 0._wp; contact_point = 0._wp
        in_contact = .false.
        if (body_a%body_kind /= body_b%body_kind) return

        select case (active_collision_model)
        case (1)
            call s_soft_sphere_pair(body_a, shift_a, body_b, shift_b, kind_friction(body_a%body_kind), force_a, torque_a, &
                                    & torque_b, contact_point, in_contact)
        end select

    end subroutine s_compute_pair_contact

    !> Model tree for a body pressed overlap past wall wall_index (1/2: x beg/end, 3/4: y, 5/6: z)
    subroutine s_compute_wall_contact(body, wall_index, overlap, force, torque)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(lagrangian_body), intent(in)   :: body
        integer, intent(in)                 :: wall_index
        real(wp), intent(in)                :: overlap
        real(wp), dimension(3), intent(out) :: force, torque

        force = 0._wp; torque = 0._wp
        select case (active_collision_model)
        case (1)
            call s_soft_sphere_wall(body, wall_index, overlap, kind_friction(body%body_kind), force, torque)
        end select

    end subroutine s_compute_wall_contact

    !> Linear spring-dashpot with Coulomb friction (Cundall & Strack, Geotechnique 29(1), 1979) between two spheres or circles
    subroutine s_soft_sphere_pair(body_a, shift_a, body_b, shift_b, friction, force_a, torque_a, torque_b, contact_point, &
                                  & in_contact)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(lagrangian_body), intent(in) :: body_a, body_b
        real(wp), dimension(3), intent(in) :: shift_a, shift_b
        real(wp), intent(in) :: friction
        real(wp), dimension(3), intent(out) :: force_a, torque_a, torque_b, contact_point
        logical, intent(out) :: in_contact
        real(wp), dimension(3) :: centroid_a, centroid_b, normal_vector, radial_vector, rotation_velocity, vel_a, vel_b, &
             & normal_velocity, tangential_vector, normal_force, tangential_force
        real(wp) :: overlap_distance, effective_mass, spring, damping

        force_a = 0._wp; torque_a = 0._wp; torque_b = 0._wp; contact_point = 0._wp

        centroid_a = [body_a%x_centroid + shift_a(1), body_a%y_centroid + shift_a(2), 0._wp]
        centroid_b = [body_b%x_centroid + shift_b(1), body_b%y_centroid + shift_b(2), 0._wp]
        if (num_dims == 3) then
            centroid_a(3) = body_a%z_centroid + shift_a(3)
            centroid_b(3) = body_b%z_centroid + shift_b(3)
        end if

        normal_vector = centroid_b - centroid_a
        overlap_distance = body_a%radius + body_b%radius - norm2(normal_vector)
        in_contact = overlap_distance > 0._wp
        if (.not. in_contact) return

        normal_vector = normal_vector/norm2(normal_vector)
        contact_point = centroid_a + normal_vector*(body_a%radius - 0.5_wp*overlap_distance)

        effective_mass = 1.0_wp/((1.0_wp/body_a%mass) + (1._wp/(body_b%mass)))
        spring = spring_stiffness*effective_mass
        damping = damping_parameter*effective_mass

        ! contact-point velocity of each body
        radial_vector = normal_vector*(body_a%radius - 0.5_wp*overlap_distance)
        call s_cross_product(body_a%angular_vel, radial_vector, rotation_velocity)
        vel_a = body_a%vel + rotation_velocity
        radial_vector = normal_vector*(-1.0_wp)*(body_b%radius - 0.5_wp*overlap_distance)
        call s_cross_product(body_b%angular_vel, radial_vector, rotation_velocity)
        vel_b = body_b%vel + rotation_velocity

        normal_velocity = dot_product(vel_a - vel_b, normal_vector)*normal_vector
        tangential_vector = (vel_a - vel_b) - normal_velocity
        if (.not. f_approx_equal(norm2(tangential_vector), 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)

        normal_force = -spring*overlap_distance*normal_vector - damping*normal_velocity
        tangential_force = -friction*norm2(normal_force)*tangential_vector
        force_a = normal_force + tangential_force

        if (body_a%body_kind == particle_type_ib) call s_cross_product(normal_vector*body_a%radius, tangential_force, torque_a)
        if (body_b%body_kind == particle_type_ib) torque_b = torque_a*body_b%radius/body_a%radius

    end subroutine s_soft_sphere_pair

    !> Spring-dashpot with Coulomb friction between a body and a flat domain wall
    subroutine s_soft_sphere_wall(body, wall_index, overlap, friction, force, torque)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(lagrangian_body), intent(in) :: body
        integer, intent(in) :: wall_index
        real(wp), intent(in) :: overlap, friction
        real(wp), dimension(3), intent(out) :: force, torque
        real(wp), dimension(3) :: normal_vector, radial_vector, rotation_velocity, velocity, normal_velocity, tangential_vector, &
             & normal_force, tangential_force
        real(wp) :: spring, damping

        ! outward normal of the wall: odd indices face the low end of an axis, even the high end
        normal_vector = 0._wp
        normal_vector((wall_index + 1)/2) = merge(-1._wp, 1._wp, mod(wall_index, 2) == 1)

        spring = spring_stiffness*body%mass
        damping = damping_parameter*body%mass

        radial_vector = normal_vector*(body%radius - overlap)
        call s_cross_product(body%angular_vel, radial_vector, rotation_velocity)
        velocity = body%vel + rotation_velocity

        normal_velocity = dot_product(velocity, normal_vector)*normal_vector
        tangential_vector = velocity - normal_velocity
        if (.not. f_approx_equal(norm2(tangential_vector), 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)
        normal_force = -spring*overlap*normal_vector - damping*normal_velocity
        tangential_force = -friction*norm2(normal_force)*tangential_vector
        force = normal_force + tangential_force

        torque = 0._wp
        if (body%body_kind == particle_type_ib) call s_cross_product(normal_vector*body%radius, tangential_force, torque)

    end subroutine s_soft_sphere_wall

    !> Add the forces and torques of the detected IB-IB pairs whose contact point this rank owns
    subroutine s_apply_ib_pair_contacts(num_considered_collisions, forces, torques)

        integer, intent(in)                            :: num_considered_collisions
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: i, l, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2
        real(wp), dimension(3)                         :: shift_1, shift_2, force_a, torque_a, torque_b, contact_point
        type(lagrangian_body)                          :: body_a, body_b
        logical                                        :: in_contact

        if (num_considered_collisions == 0) return

        $:GPU_PARALLEL_LOOP(private='[i, l, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2, shift_1, &
                            & shift_2, force_a, torque_a, torque_b, contact_point, body_a, body_b, in_contact]', copy='[forces, torques]')
        do i = 1, num_considered_collisions
            encoded_pid1 = collision_lookup(i, 3)
            encoded_pid2 = collision_lookup(i, 4)
            call s_decode_patch_periodicity(encoded_pid1, pid1, xp1, yp1, zp1)
            call s_decode_patch_periodicity(encoded_pid2, pid2, xp2, yp2, zp2)
            pid1 = collision_lookup(i, 1)
            pid2 = collision_lookup(i, 2)

            if (pid1 <= 0 .or. pid2 <= 0) cycle

            shift_1 = 0._wp; shift_2 = 0._wp
            shift_1(1) = real(xp1, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
            shift_1(2) = real(yp1, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
            shift_2(1) = real(xp2, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
            shift_2(2) = real(yp2, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
            if (num_dims == 3) then
                shift_1(3) = real(zp1, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)
                shift_2(3) = real(zp2, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)
            end if

            body_a = patch_ib(pid1)%lagrangian_body
            body_b = patch_ib(pid2)%lagrangian_body
            call s_compute_pair_contact(body_a, shift_1, body_b, shift_2, force_a, torque_a, torque_b, contact_point, in_contact)

            ! pid1 is a rank-local index, so owning the pair by its centroid drops or doubles pairs split across ranks
            if (in_contact) then
                if (f_local_rank_owns_location(contact_point, glb_bounds)) then
                    do l = 1, num_dims
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid1, l) = forces(pid1, l) + force_a(l)
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid1, l) = torques(pid1, l) + torque_a(l)
                        ! equal and opposite force on the second IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid2, l) = forces(pid2, l) - force_a(l)
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid2, l) = torques(pid2, l) + torque_b(l)
                    end do
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_ib_pair_contacts

    !> Add the forces and torques of IBs pressed into the domain walls, on the rank that owns each IB's centroid
    subroutine s_apply_ib_wall_contacts(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: patch_id, i, l
        real(wp), dimension(3)                         :: collision_location, force, torque
        type(lagrangian_body)                          :: body

        $:GPU_PARALLEL_LOOP(private='[patch_id, i, l, collision_location, force, torque, body]', copy='[forces, torques]', &
                            & collapse=2)
        do patch_id = 1, num_ibs
            do i = 1, num_dims*2
                if (f_approx_equal(wall_overlap_distances(patch_id, i), 0._wp)) cycle

                collision_location = [patch_ib(patch_id)%x_centroid, patch_ib(patch_id)%y_centroid, 0._wp]
                if (num_dims == 3) collision_location(3) = patch_ib(patch_id)%z_centroid
                if (f_local_rank_owns_location(collision_location, glb_bounds)) then
                    body = patch_ib(patch_id)%lagrangian_body
                    call s_compute_wall_contact(body, i, wall_overlap_distances(patch_id, i), force, torque)
                    do l = 1, num_dims
                        $:GPU_ATOMIC(atomic='update')
                        forces(patch_id, l) = forces(patch_id, l) + force(l)
                        $:GPU_ATOMIC(atomic='update')
                        torques(patch_id, l) = torques(patch_id, l) + torque(l)
                    end do
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_ib_wall_contacts

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

    !> @brief uses boundary conditions and particle locations to check for wall conditions
    subroutine s_detect_wall_collisions(any_wall_collision)

        logical, intent(out) :: any_wall_collision
        integer              :: gp_idx, i, j, k, patch_id
        real(wp)             :: edge_location, overlap_distance, max_overlap

        max_overlap = 0._wp
        $:GPU_PARALLEL_LOOP(private='[patch_id, edge_location, overlap_distance]', reduction='[[max_overlap]]', reductionOp='[max]')
        do patch_id = 1, num_ibs
            #:for X, DIR, IDX in [('x', 1, 1), ('y', 2, 3), ('z', 3, 5)]
                ! check if the boundaries are either of the two conditions we should compute collisions with
                if (ib_bc_${X}$%beg == BC_SLIP_WALL .or. ib_bc_${X}$%beg == BC_NO_SLIP_WALL) then
                    ! get the location of the true IB surface towards the domain boundary
                    edge_location = patch_ib(patch_id)%${X}$_centroid - patch_ib(patch_id)%radius
                    ! check if that edge actually extends out of the comutational domain
                    if (edge_location < glb_bounds(${DIR}$)%beg) then
                        ! the distance that the IB extends out of the domain
                        overlap_distance = glb_bounds(${DIR}$)%beg - edge_location
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$) = overlap_distance
                    max_overlap = max(max_overlap, overlap_distance)
                end if

                if (ib_bc_${X}$%end == BC_SLIP_WALL .or. ib_bc_${X}$%end == BC_NO_SLIP_WALL) then
                    edge_location = patch_ib(patch_id)%${X}$_centroid + patch_ib(patch_id)%radius
                    if (edge_location > glb_bounds(${DIR}$)%end) then
                        overlap_distance = edge_location - glb_bounds(${DIR}$)%end
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$ + 1) = overlap_distance
                    max_overlap = max(max_overlap, overlap_distance)
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

        any_wall_collision = max_overlap > 0._wp

    end subroutine s_detect_wall_collisions

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

        if (.not. collisions_initialized) return
        collisions_initialized = .false.
        if (ib) then
            @:DEALLOCATE(collision_lookup)
            @:DEALLOCATE(wall_overlap_distances)
        end if

    end subroutine s_finalize_collisions_module

end module m_collisions
