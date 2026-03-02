!>
!! @file
!! @brief Contains module m_collisions

#:include 'macros.fpp'

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the flow state
module m_collisions

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_constants

    use m_compute_levelset

    use m_ib_patches

    use m_model

    implicit none

    private; public;

contains

    subroutine s_detect_ib_wall_collisions(gps)

        type(ghost_point), dimension(num_gps), intent(in) :: gps
        real(wp), dimension(num_ibs, (num_dims*2)), intent(out) :: overlap_distances

        integer :: gp_idx, i, j, k, patch_id
        real(wp) :: edge_location, overlap_distance

        collision_distances = 0._wp

        do gp_idx = 1, num_gps
            i = gps(gp_idx)%loc(1)
            j = gps(gp_idx)%loc(2)
            if (p /= 0) k = gps(gp_idx)%loc(3)
            patch_id = gps(gp_idx)%ib_patch_id

            if (bc_x%beg == BC_SLIP_WALL .or. bc_x%beg == BC_NO_SLIP_WALL) then
                edge_location = x_cc(i) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(1))
                if (edge_location < x_domain%beg) then
                    overlap_distance = x_domain%beg - edge_location
                    overlap_distances(patch_id, 1) = max(overlap_distances(patch_id, 1), overlap_distance)
                end if
            end if

            if (bc_x%end == BC_SLIP_WALL .or. bc_x%end == BC_NO_SLIP_WALL) then
                edge_location = x_cc(i) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(1))
                if (edge_location > x_domain%end) then
                    overlap_distance = edge_location - x_domain%end
                    overlap_distances(patch_id, 2) = max(overlap_distances(patch_id, 2), overlap_distance)
                end if
            end if

            if (bc_y%beg == BC_SLIP_WALL .or. bc_y%beg == BC_NO_SLIP_WALL) then
                edge_location = y_cc(j) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(2))
                if (edge_location < y_domain%beg) then
                    overlap_distance = y_domain%beg - edge_location
                    overlap_distances(patch_id, 3) = max(overlap_distances(patch_id, 3), overlap_distance)
                end if
            end if

            if (bc_y%end == BC_SLIP_WALL .or. bc_y%end == BC_NO_SLIP_WALL) then
                edge_location = y_cc(j) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(2))
                if (edge_location > y_domain%end) then
                    overlap_distance = edge_location - y_domain%end
                    overlap_distances(patch_id, 4) = max(overlap_distances(patch_id, 4), overlap_distance)
                end if
            end if

            if (p > 0) then
                if (bc_z%beg == BC_SLIP_WALL .or. bc_z%beg == BC_NO_SLIP_WALL) then
                    edge_location = z_cc(k) - abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(3))
                    if (edge_location < z_domain%beg) then
                        overlap_distance = z_domain%beg - edge_location
                        overlap_distances(patch_id, 6) = max(overlap_distances(patch_id, 5), overlap_distance)
                    end if
                end if

                if (bc_z%end == BC_SLIP_WALL .or. bc_z%end == BC_NO_SLIP_WALL) then
                    edge_location = z_cc(k) + abs(gps(gp_idx)%levelset * gps(gp_idx)%levelset_norm(3))
                    if (edge_location > z_domain%end) then
                        overlap_distance = edge_location - z_domain%end
                        overlap_distances(patch_id, 6) = max(overlap_distances(patch_id, 6), overlap_distance)
                    end if
                end if
            end if
        end do

    end subroutine s_detect_ib_wall_collisions

end module m_collisions
