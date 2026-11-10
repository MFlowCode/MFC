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

    private; public :: s_apply_collision_forces, s_initialize_collisions_module, s_finalize_collisions_module, &
        & f_neighborhood_ranks_own_location, collisions_active
    ! overlap distances for computing collisions
    integer, allocatable, dimension(:,:)  :: collision_lookup
    real(wp), allocatable, dimension(:,:) :: wall_overlap_distances
    real(wp)                              :: spring_stiffness, damping_parameter
    $:GPU_DECLARE(create='[spring_stiffness, damping_parameter]')
    $:GPU_DECLARE(create='[collision_lookup, wall_overlap_distances]')

    !> true when any IB-IB or IB-wall contact was detected on this rank since the last adaptive-dt computation
    logical :: collisions_active

    !> One bit per unordered pair (a < b) of this rank's IBs (local + neighborhood indices, a triangular index over
    !! num_local_ibs_max): the ghost-point candidate pass claims each pair with an atomic OR, so only the first thread to find a
    !! pair records it and duplicates are dropped on the GPU. Bits are cleared after the forces are applied, so the set is empty
    !! between calls. num_local_ibs_max = 30000 gives 4.5e8 bits (56 MB).
    integer, allocatable, dimension(:) :: pair_bits
    $:GPU_DECLARE(create='[pair_bits]')

contains

    subroutine s_initialize_collisions_module()

        real(wp) :: e

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness]')

        @:ALLOCATE(collision_lookup(num_local_ibs_max * 27 * 8, 4))
        @:ALLOCATE(wall_overlap_distances(num_local_ibs_max*27, 6))

        @:ALLOCATE(pair_bits(f_pair_bit_index(num_local_ibs_max - 1, num_local_ibs_max)/32 + 1))
        pair_bits = 0
        $:GPU_UPDATE(device='[pair_bits]')

        wall_overlap_distances = 0
        collisions_active = .false.
        $:GPU_UPDATE(device='[wall_overlap_distances]')
        $:GPU_UPDATE(device='[ib_coefficient_of_friction]')

    end subroutine s_initialize_collisions_module

    subroutine s_apply_collision_forces(ghost_points, num_gps, ib_markers, forces, torques)

        type(ghost_point), dimension(:), intent(in)    :: ghost_points
        integer, intent(in)                            :: num_gps
        type(integer_field), intent(in)                :: ib_markers
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: num_considered_collisions
        logical                                        :: any_wall_collision

        ! return if no collisions

        if (collision_model == 0) return

        ! get is distance used in the force calculation with each IB and each wall
        call s_detect_wall_collisions(any_wall_collision)
        call s_detect_ib_collisions(ghost_points, ib_markers, num_gps, num_considered_collisions)

        ! accumulate across RK stages; consumed (and reset) by s_compute_dt once per time step
        collisions_active = collisions_active .or. any_wall_collision .or. (num_considered_collisions > 0)

        select case (collision_model)
        case (1)  ! soft sphere model
            call s_apply_wall_collision_forces_soft_sphere(forces, torques)
            call s_apply_ib_collision_forces_soft_sphere(num_considered_collisions, forces, torques)
        end select

        call s_clear_pair_bits(num_considered_collisions)

    end subroutine s_apply_collision_forces

    !> @brief applies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_ib_collision_forces_soft_sphere(num_considered_collisions, forces, torques)

        integer, intent(in) :: num_considered_collisions
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer :: i, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2, l  ! iterators and patch IDs
        real(wp) :: overlap_distance
        real(wp), dimension(3) :: normal_vector, centroid_1, centroid_2, contact_point
        real(wp), dimension(3) :: normal_velocity, tangential_vector, normal_force, tangential_force, torque, radial_vector, &
             & rotation_velocity, vel1, vel2
        real(wp) :: k, eta, effective_mass  ! the spring stiffness and damping coefficient and mass of a specific interaction

        if (num_considered_collisions == 0) return

        ! Iterate over all collisions detected
        $:GPU_PARALLEL_LOOP(private='[i, l, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2, centroid_1, &
                            & centroid_2, normal_vector, overlap_distance, effective_mass, k, eta, normal_velocity, &
                            & tangential_vector, normal_force, tangential_force, torque, radial_vector, rotation_velocity, vel1, &
                            & vel2, contact_point]', copy='[forces, torques]')
        do i = 1, num_considered_collisions
            encoded_pid1 = collision_lookup(i, 3)
            encoded_pid2 = collision_lookup(i, 4)
            call s_decode_patch_periodicity(encoded_pid1, pid1, xp1, yp1, zp1)
            call s_decode_patch_periodicity(encoded_pid2, pid2, xp2, yp2, zp2)
            pid1 = collision_lookup(i, 1)
            pid2 = collision_lookup(i, 2)

            if (pid1 <= 0 .or. pid2 <= 0) cycle

            centroid_1(1) = patch_ib(pid1)%x_centroid + real(xp1, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
            centroid_1(2) = patch_ib(pid1)%y_centroid + real(yp1, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
            centroid_1(3) = 0._wp
            centroid_2(1) = patch_ib(pid2)%x_centroid + real(xp2, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
            centroid_2(2) = patch_ib(pid2)%y_centroid + real(yp2, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
            centroid_2(3) = 0._wp
            if (num_dims == 3) then
                centroid_1(3) = patch_ib(pid1)%z_centroid + real(zp1, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)
                centroid_2(3) = patch_ib(pid2)%z_centroid + real(zp2, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)
            end if

            normal_vector = centroid_2 - centroid_1
            overlap_distance = patch_ib(pid1)%radius + patch_ib(pid2)%radius - norm2(normal_vector)
            if (overlap_distance > 0._wp) then  ! if the two patches are close enough to collide
                normal_vector = normal_vector/norm2(normal_vector)
                ! pid1 is a rank-local index, so owning the pair by its centroid drops or doubles pairs split across ranks
                contact_point = centroid_1 + normal_vector*(patch_ib(pid1)%radius - 0.5_wp*overlap_distance)
                if (f_local_rank_owns_location(contact_point, glb_bounds)) then
                    ! compute constants of the collision
                    effective_mass = 1.0_wp/((1.0_wp/patch_ib(pid1)%mass) + (1._wp/(patch_ib(pid2)%mass)))
                    k = spring_stiffness*effective_mass
                    eta = damping_parameter*effective_mass

                    ! Get the vectors and velcoities
                    radial_vector = normal_vector*(patch_ib(pid1)%radius - 0.5_wp*overlap_distance)
                    call s_cross_product(patch_ib(pid1)%angular_vel, radial_vector, rotation_velocity)
                    vel1 = patch_ib(pid1)%vel + rotation_velocity
                    radial_vector = normal_vector*(-1.0_wp)*(patch_ib(pid2)%radius - 0.5_wp*overlap_distance)
                    call s_cross_product(patch_ib(pid2)%angular_vel, radial_vector, rotation_velocity)
                    vel2 = patch_ib(pid2)%vel + rotation_velocity

                    normal_velocity = dot_product(vel1 - vel2, normal_vector)*normal_vector
                    tangential_vector = (vel1 - vel2) - normal_velocity
                    if (.not. f_approx_equal(norm2(tangential_vector), &
                        & 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)

                    ! compute force and torque
                    normal_force = -k*overlap_distance*normal_vector - eta*normal_velocity
                    tangential_force = -ib_coefficient_of_friction*norm2(normal_force)*tangential_vector
                    call s_cross_product(normal_vector*patch_ib(pid1)%radius, tangential_force, torque)

                    do l = 1, num_dims
                        ! update the first IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid1, l) = forces(pid1, l) + (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid1, l) = torques(pid1, l) + torque(l)

                        ! apply equal and opposite force/torque to second IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid2, l) = forces(pid2, l) - (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid2, l) = torques(pid2, l) + torque(l)*patch_ib(pid2)%radius/patch_ib(pid1)%radius
                    end do
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_ib_collision_forces_soft_sphere

    !> @brief applies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_wall_collision_forces_soft_sphere(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer :: patch_id, i, l
        real(wp), dimension(3) :: normal_force, tangential_force, normal_vector, normal_velocity, tangential_vector, &
             & collision_location, torque, radial_vector, rotation_velocity, velocity
        real(wp) :: k, eta  ! the spring stiffness and damping coefficient for a specific IB

        $:GPU_PARALLEL_LOOP(private='[patch_id, i, l, collision_location, normal_vector, k, eta, normal_velocity, &
                            & tangential_vector, normal_force, tangential_force, torque, radial_vector, rotation_velocity, &
                            & velocity]', copy='[forces, torques]', collapse=2)
        do patch_id = 1, num_ibs
            do i = 1, num_dims*2
                ! only compute force contributions if there was an overlap
                if (f_approx_equal(wall_overlap_distances(patch_id, i), 0._wp)) cycle

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

                ! ensure the local rank owns that collision before proceeding
                collision_location = [patch_ib(patch_id)%x_centroid, patch_ib(patch_id)%y_centroid, 0._wp]
                if (num_dims == 3) collision_location(3) = patch_ib(patch_id)%z_centroid
                if (f_local_rank_owns_location(collision_location, glb_bounds)) then
                    k = spring_stiffness*patch_ib(patch_id)%mass
                    eta = damping_parameter*patch_ib(patch_id)%mass

                    ! get the vector that points from the centroid to the point of collision
                    radial_vector = normal_vector*(patch_ib(patch_id)%radius - wall_overlap_distances(patch_id, i))
                    ! convert the angular velocity to linear velocity
                    call s_cross_product(patch_ib(patch_id)%angular_vel, radial_vector, rotation_velocity)
                    velocity = patch_ib(patch_id)%vel + rotation_velocity

                    ! standard soft-sphere collision  with the wall
                    normal_velocity = dot_product(velocity, normal_vector)*normal_vector
                    tangential_vector = velocity - normal_velocity
                    if (.not. f_approx_equal(norm2(tangential_vector), &
                        & 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)
                    normal_force = -k*wall_overlap_distances(patch_id, i)*normal_vector - eta*normal_velocity
                    tangential_force = -ib_coefficient_of_friction*norm2(normal_force)*tangential_vector
                    call s_cross_product(normal_vector*patch_ib(patch_id)%radius, tangential_force, torque)

                    do l = 1, num_dims
                        $:GPU_ATOMIC(atomic='update')
                        forces(patch_id, l) = forces(patch_id, l) + (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(patch_id, l) = torques(patch_id, l) + torque(l)
                    end do
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_wall_collision_forces_soft_sphere

    !> Triangular index (0-based) of the unordered pair a < b, both in 1 ... num_local_ibs_max
    pure function f_pair_bit_index(a, b) result(idx)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: a, b
        integer(kind=8)     :: idx

        idx = int(a - 1, 8)*int(2*num_local_ibs_max - a, 8)/2 + int(b - a - 1, 8)

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
        type(integer_field), intent(in) :: ib_markers
        integer, intent(in) :: num_gps
        integer, intent(out) :: num_considered_collisions
        integer :: i, j, k, ii, jj, kk, z_bound, gp_idx, gp_patch_id, neighbor_patch_id
        integer :: a, b, gbl_a, gbl_b, enc_a, enc_b, xp, yp, zp, word, bit, mask, old_word, slot
        integer :: num_found, max_pairs
        integer(kind=8) :: idx

        num_found = 0
        max_pairs = size(collision_lookup, 1)
        z_bound = 0; if (num_dims == 3) z_bound = 2

        $:GPU_PARALLEL_LOOP(private='[gp_idx, i, j, k, ii, jj, kk, gp_patch_id, neighbor_patch_id, a, b, gbl_a, gbl_b, enc_a, &
                            & enc_b, xp, yp, zp, idx, word, bit, mask, old_word, slot]', copy='[num_found]', copyin='[z_bound, max_pairs]')
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
            idx = f_pair_bit_index(a, b)
            word = int(idx/32) + 1
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
        integer             :: n

        $:GPU_PARALLEL_LOOP(private='[n]')
        do n = 1, num_considered_collisions
            pair_bits(int(f_pair_bit_index(collision_lookup(n, 1), collision_lookup(n, 2))/32) + 1) = 0
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_clear_pair_bits

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

        @:DEALLOCATE(pair_bits)
        @:DEALLOCATE(collision_lookup)
        @:DEALLOCATE(wall_overlap_distances)

    end subroutine s_finalize_collisions_module

end module m_collisions
