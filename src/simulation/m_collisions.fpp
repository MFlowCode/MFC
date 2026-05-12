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

    implicit none

    private; public :: s_apply_collision_forces, s_initialize_collisions_module, s_finalize_collisions_module
    ! overlap distances for computing collisions
    integer, allocatable, dimension(:,:)  :: collision_lookup
    real(wp), allocatable, dimension(:,:) :: wall_overlap_distances
    real(wp)                              :: spring_stiffness, damping_parameter
    $:GPU_DECLARE(create='[spring_stiffness, damping_parameter]')
    $:GPU_DECLARE(create='[collision_lookup, wall_overlap_distances]')

contains

    subroutine s_initialize_collisions_module()

        real(wp) :: e

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness]')

        @:ALLOCATE(collision_lookup(num_ibs * 8, 4))
        @:ALLOCATE(wall_overlap_distances(num_ibs, 6))

        wall_overlap_distances = 0
        $:GPU_UPDATE(device='[wall_overlap_distances]')
        $:GPU_UPDATE(device='[ib_coefficient_of_friction]')

    end subroutine s_initialize_collisions_module

    subroutine s_apply_collision_forces(ghost_points, num_gps, ib_markers, forces, torques)

        type(ghost_point), dimension(:), intent(in)    :: ghost_points
        integer, intent(in)                            :: num_gps
        type(integer_field), intent(in)                :: ib_markers
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer                                        :: num_considered_collisions

        ! return if no collisions

        if (collision_model == 0) return

        ! get is distance used in the force calculation with each IB and each wall
        call s_detect_wall_collisions()
        call s_detect_ib_collisions(ghost_points, ib_markers, num_gps, num_considered_collisions)

        select case (collision_model)
        case (1)  ! soft sphere model
            call s_apply_wall_collision_forces_soft_sphere(forces, torques)
            call s_apply_ib_collision_forces_soft_sphere(num_considered_collisions, forces, torques)
        end select

    end subroutine s_apply_collision_forces

    !> @brief applies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_ib_collision_forces_soft_sphere(num_considered_collisions, forces, torques)

        integer, intent(in) :: num_considered_collisions
        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer :: i, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2, l  ! iterators and patch IDs
        real(wp) :: overlap_distance
        real(wp), dimension(3) :: normal_vector, centroid_1, centroid_2
        real(wp), dimension(3) :: normal_velocity, tangental_vector, normal_force, tangental_force, torque, radial_vector, &
             & rotation_velocity, vel1, vel2
        real(wp) :: k, eta, effective_mass  ! the spring stiffness and damping coefficient and mass of a specific interaction

        if (num_considered_collisions == 0) return

        ! print *, "Checking Collisions: ", num_considered_collisions, " on rank ", proc_rank

        ! Iterate over all collisions detected
        $:GPU_PARALLEL_LOOP(private='[i, l, encoded_pid1, encoded_pid2, xp1, xp2, yp1, yp2, zp1, zp2, pid1, pid2, centroid_1, &
                            & centroid_2, normal_vector, overlap_distance, effective_mass, k, eta, normal_velocity, &
                            & tangental_vector, normal_force, tangental_force, torque, radial_vector, rotation_velocity, vel1, &
                            & vel2]', copy='[forces, torques]')
        do i = 1, num_considered_collisions
            encoded_pid1 = collision_lookup(i, 3)
            encoded_pid2 = collision_lookup(i, 4)
            call s_decode_patch_periodicity(encoded_pid1, pid1, xp1, yp1, zp1)
            call s_decode_patch_periodicity(encoded_pid2, pid2, xp2, yp2, zp2)

            centroid_1(1) = patch_ib(pid1)%x_centroid + real(xp1, wp)*(x_domain%end - x_domain%beg)
            centroid_1(2) = patch_ib(pid1)%y_centroid + real(yp1, wp)*(y_domain%end - y_domain%beg)
            centroid_1(3) = 0._wp
            centroid_2(1) = patch_ib(pid2)%x_centroid + real(xp2, wp)*(x_domain%end - x_domain%beg)
            centroid_2(2) = patch_ib(pid2)%y_centroid + real(yp2, wp)*(y_domain%end - y_domain%beg)
            centroid_2(3) = 0._wp
            if (num_dims == 3) then
                centroid_1(3) = patch_ib(pid1)%z_centroid + real(zp1, wp)*(z_domain%end - z_domain%beg)
                centroid_2(3) = patch_ib(pid2)%z_centroid + real(zp2, wp)*(z_domain%end - z_domain%beg)
            end if

            normal_vector = centroid_2 - centroid_1
            overlap_distance = patch_ib(pid1)%radius + patch_ib(pid2)%radius - norm2(normal_vector)
            if (overlap_distance > 0._wp) then  ! if the two patches are close enough to collide
                normal_vector = normal_vector/norm2(normal_vector)
                if (f_local_rank_owns_collision(centroid_1)) then
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
                    tangental_vector = (vel1 - vel2) - normal_velocity
                    if (.not. f_approx_equal(norm2(tangental_vector), &
                        & 0._wp)) tangental_vector = tangental_vector/norm2(tangental_vector)

                    ! compute force and torque
                    normal_force = -k*overlap_distance*normal_vector - eta*normal_velocity
                    tangental_force = -ib_coefficient_of_friction*norm2(normal_force)*tangental_vector
                    call s_cross_product(normal_vector*patch_ib(pid1)%radius, tangental_force, torque)

                    do l = 1, num_dims
                        ! update the first IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid1, l) = forces(pid1, l) + (normal_force(l) + tangental_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid1, l) = torques(pid1, l) + torque(l)

                        ! apply equal and opposite force/torque to second IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid2, l) = forces(pid2, l) - (normal_force(l) + tangental_force(l))
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
        real(wp), dimension(3) :: normal_force, tangental_force, normal_vector, normal_velocity, tangental_vector, &
             & collision_location, torque, radial_vector, rotation_velocity, velocity
        real(wp) :: k, eta  ! the spring stiffness and damping coefficient for a specific IB

        $:GPU_PARALLEL_LOOP(private='[patch_id, i, l, collision_location, normal_vector, k, eta, normal_velocity, &
                            & tangental_vector, normal_force, tangental_force, torque, radial_vector, rotation_velocity, &
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
                if (f_local_rank_owns_collision(collision_location)) then
                    k = spring_stiffness*patch_ib(patch_id)%mass
                    eta = damping_parameter*patch_ib(patch_id)%mass

                    ! get the vector that points from the centroid to the point of collision
                    radial_vector = normal_vector*(patch_ib(patch_id)%radius - wall_overlap_distances(patch_id, i))
                    ! convert the angular velocity to linear velocity
                    call s_cross_product(patch_ib(patch_id)%angular_vel, radial_vector, rotation_velocity)
                    velocity = patch_ib(patch_id)%vel + rotation_velocity

                    ! standard soft-sphere collision  with the wall
                    normal_velocity = dot_product(velocity, normal_vector)*normal_vector
                    tangental_vector = velocity - normal_velocity
                    if (.not. f_approx_equal(norm2(tangental_vector), &
                        & 0._wp)) tangental_vector = tangental_vector/norm2(tangental_vector)
                    normal_force = -k*wall_overlap_distances(patch_id, i)*normal_vector - eta*normal_velocity
                    tangental_force = -ib_coefficient_of_friction*norm2(normal_force)*tangental_vector
                    call s_cross_product(normal_vector*patch_ib(patch_id)%radius, tangental_force, torque)

                    do l = 1, num_dims
                        $:GPU_ATOMIC(atomic='update')
                        forces(patch_id, l) = forces(patch_id, l) + (normal_force(l) + tangental_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(patch_id, l) = torques(patch_id, l) + torque(l)
                    end do
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_wall_collision_forces_soft_sphere

    !> uses ghost-point/image-point information to determine if it is possible if two IBs are colliding, effectively an optimized
    !! nearest neighbor search
    subroutine s_detect_ib_collisions(gps, ib_markers, num_gps, num_considered_collisions)

        type(ghost_point), dimension(num_gps), intent(in) :: gps
        type(integer_field), intent(in)                   :: ib_markers
        integer, intent(in)                               :: num_gps
        integer, intent(out)                              :: num_considered_collisions
        integer                                           :: i, j, k, z_bound, ii, jj, kk
        integer, dimension(2)                             :: decoded_pairs
        integer                                           :: gp_idx, gp_patch_id, neighbor_patch_id
        integer                                           :: pair_idx, out_idx
        logical                                           :: already_found

        ! Temporary array to hold all detected pairs (with potential duplicates)
        integer, dimension(num_gps, 2) :: raw_pairs
        integer                        :: num_raw, local_num_raw

        num_raw = 0
        z_bound = 0; if (num_dims == 3) z_bound = 1

        $:GPU_PARALLEL_LOOP(private='[gp_idx, gp_patch_id, neighbor_patch_id, local_num_raw, i, j, k, ii, jj, kk]', &
                            & copy='[raw_pairs, num_raw]', copyin='[z_bound]')
        do gp_idx = 1, num_gps
            i = gps(gp_idx)%loc(1)
            j = gps(gp_idx)%loc(2)
            k = 0; if (num_dims == 3) k = gps(gp_idx)%loc(3)
            gp_patch_id = ib_markers%sf(i, j, k)

            ! search in a cube around the BG for Ib markers belonging to another patch
            neighbor_search: do ii = i - 1, i + 1
                do jj = j - 1, j + 1
                    do kk = k - z_bound, k + z_bound
                        neighbor_patch_id = ib_markers%sf(ii, jj, kk)

                        ! If any neighbors are of a different/higher marker value, we consider it for possible collision
                        if (gp_patch_id < neighbor_patch_id) then
                            $:GPU_ATOMIC(atomic='capture')
                            num_raw = num_raw + 1
                            local_num_raw = num_raw
                            $:END_GPU_ATOMIC_CAPTURE()

                            ! Store with smaller ID first for consistent ordering
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
                collision_lookup(num_considered_collisions, 1) = decoded_pairs(1)
                collision_lookup(num_considered_collisions, 2) = decoded_pairs(2)
                collision_lookup(num_considered_collisions, 3) = raw_pairs(pair_idx, 1)
                collision_lookup(num_considered_collisions, 4) = raw_pairs(pair_idx, 2)
            end if
        end do
        $:GPU_UPDATE(device='[collision_lookup]')

    end subroutine s_detect_ib_collisions

    subroutine s_detect_ib_collisions_n2(num_considered_collisions)

        integer, intent(out)   :: num_considered_collisions
        integer                :: pid1, pid2, encoded_pid1, encoded_pid2, current_collisions
        real(wp), dimension(3) :: centroid_1, centroid_2, distance_vec

        num_considered_collisions = 0

        $:GPU_PARALLEL_LOOP(private='[pid1, pid2, centroid_1, centroid_2, distance_vec, current_collisions]', &
                            & copy='[num_considered_collisions]')
        do pid1 = 1, num_ibs - 1
            do pid2 = pid1 + 1, num_ibs
                centroid_1 = [patch_ib(pid1)%x_centroid, patch_ib(pid1)%y_centroid, 0._wp]
                centroid_2 = [patch_ib(pid2)%x_centroid, patch_ib(pid2)%y_centroid, 0._wp]
                if (num_dims == 3) then
                    centroid_1(3) = patch_ib(pid1)%z_centroid
                    centroid_2(3) = patch_ib(pid2)%z_centroid
                end if
                distance_vec = centroid_2 - centroid_1

                if (norm2(distance_vec) < patch_ib(pid1)%radius + patch_ib(pid2)%radius) then
                    $:GPU_ATOMIC(atomic='capture')
                    num_considered_collisions = num_considered_collisions + 1
                    current_collisions = num_considered_collisions
                    $:END_GPU_ATOMIC_CAPTURE()

                    collision_lookup(current_collisions, 1) = pid1
                    collision_lookup(current_collisions, 2) = pid2
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_detect_ib_collisions_n2

    !> @brief uses boundary conditions and particle locations to check for wall conditions
    subroutine s_detect_wall_collisions()

        integer  :: gp_idx, i, j, k, patch_id
        real(wp) :: edge_location, overlap_distance

        ! iterate over all ghost points to detect the one that is most-overlapping in each direction

        $:GPU_PARALLEL_LOOP(private='[patch_id, edge_location, overlap_distance]')
        do patch_id = 1, num_ibs
            #:for X, IDX in [('x', 1), ('y', 3), ('z', 5)]
                ! check if the boundaries are either of the two conditions we should compute collisions with
                if (ib_bc_${X}$%beg == BC_SLIP_WALL .or. ib_bc_${X}$%beg == BC_NO_SLIP_WALL) then
                    ! get the location of the true IB surface towards the domain boundary
                    edge_location = patch_ib(patch_id)%${X}$_centroid - patch_ib(patch_id)%radius
                    ! check if that edge actually extends out of the comutational domain
                    if (edge_location < ${X}$_domain%beg) then
                        overlap_distance = ${X}$_domain%beg - edge_location  ! the distance that the IB extends out of the domain
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$) = overlap_distance
                end if

                if (ib_bc_${X}$%end == BC_SLIP_WALL .or. ib_bc_${X}$%end == BC_NO_SLIP_WALL) then
                    edge_location = patch_ib(patch_id)%${X}$_centroid + patch_ib(patch_id)%radius
                    if (edge_location > ${X}$_domain%end) then
                        overlap_distance = edge_location - ${X}$_domain%end
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$ + 1) = overlap_distance
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_detect_wall_collisions

    !> @brief function checks if this local MPI processor owns this specific collision
    function f_local_rank_owns_collision(collision_location) result(owns_collision)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: collision_location
        logical                            :: owns_collision
        real(wp), dimension(3)             :: projected_location

        #:if defined('MFC_MPI')
            if (num_procs == 1) then
                owns_collision = .true.
            else
                projected_location(:) = collision_location(:)

                ! catch the edge case where th collision lies just outside the computational domain
                #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
                    if (num_dims >= ${ID}$) then
                        if (ib_bc_${X}$%beg /= BC_PERIODIC) then
                            ! if it is outside the domain in one direction, project it somewhere inside so at least one rank owns it
                            if (collision_location(${ID}$) < ${X}$_domain%beg) then
                                projected_location(${ID}$) = ${X}$_domain%beg
                            else if (${X}$_domain%end < collision_location(${ID}$)) then
                                projected_location(${ID}$) = ${X}$_domain%end - 1.0e-10_wp
                            end if
                        end if
                    end if
                #:endfor

                ! the object that contains the collision location owns the collisions
                owns_collision = x%cb(-1) <= projected_location(1) .and. projected_location(1) < x%cb(m)
                owns_collision = owns_collision .and. y%cb(-1) <= projected_location(2) .and. projected_location(2) < y%cb(n)
                if (num_dims == 3) owns_collision = owns_collision .and. z%cb(-1) <= projected_location(3) &
                    & .and. projected_location(3) < z%cb(p)
            end if
        #:else
            owns_collision = .true.
        #:endif

    end function f_local_rank_owns_collision

    subroutine s_finalize_collisions_module()

        @:DEALLOCATE(collision_lookup)
        @:DEALLOCATE(wall_overlap_distances)

    end subroutine s_finalize_collisions_module

end module m_collisions
