!>
!! @file
!! @brief Contains module m_ibm

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Per-block fine-grid immersed-boundary state for AMR. A separate module from m_ibm: declaring these derived-type
!! allocatables inside m_ibm corrupts CCE's OpenMP declare-target descriptor table for m_ibm's ghost_points (lib-4425 on its
!! @:ALLOCATE, even for non-AMR runs). Do not fold them back into m_ibm.
module m_ibm_fine

    use m_derived_types

    implicit none

    !> Per-block fine IB state: each AMR slot keeps a host-side copy of its markers, computed at fine resolution and copied into
    !! m_ibm's declare-target ib_markers for the fine advance (that array is never reallocated or pointer-swapped, since that
    !! corrupts the Cray present table; the swap syncs via GPU_UPDATE). Ghost-point lists park on-device in gp_park; num_gps is the
    !! slot's ghost-point count across a swap, off its gp_park offset and cap the entries reserved there (0 = none).
    type ib_fine_state
        type(integer_field) :: markers
        integer             :: num_gps
        integer(kind=8)     :: off, cap
    end type ib_fine_state
    type(ib_fine_state), allocatable :: ib_fine(:)
    integer                          :: num_gps_save

    !> Device-resident park for the coarse and fine ghost-point lists across a fine swap, copied to/from ghost_points by on-device
    !! kernels. One packed pool: entries 1..park_coarse hold the coarse list; each slot holding ghost points (a block containing a
    !! body) gets a segment after it, so the pool scales with the ghost points parked, not blocks x the longest list. Mapped via
    !! move_alloc + GPU_ENTER_DATA (the amr_cg idiom; a direct @:ALLOCATE aborts with lib-4425 on CCE OpenMP).
    type(ghost_point), allocatable :: gp_park(:)

    !> The coarse ghost points and markers park (host copies) in an extra ib_fine slot, ib_coarse_slot (= nslots+1): a new
    !! module-level derived-type allocatable here would corrupt a sibling's descriptor on CCE OpenMP (lib-4425).
    integer         :: ib_coarse_slot = 0
    integer(kind=8) :: park_coarse = 0_8  !< gp_park entries reserved for the coarse list
    integer(kind=8) :: park_top = 0_8     !< gp_park entries in use (coarse segment + fine segments)

    !> Fine-block ghost-point capacity (= buffered fine-block cell count), set by s_ibm_alloc_fine; 0 when AMR-IB is inactive.
    integer(kind=8) :: fine_gps_cap = 0_8

    !> Bounds of ib_markers, enclosing both coarse and fine blocks so it never needs reallocating. Set by s_ibm_alloc_fine.
    integer :: mkr_lo(3) = 0, mkr_hi(3) = 0
end module m_ibm_fine

!> @brief Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the
!! flow state
module m_ibm

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_compile_specific, only: s_create_directory
    use m_delay_file_access, only: s_delay_file_access
    use m_variables_conversion
    use m_eos
    use m_helper
    use m_helper_basic
    use m_constants
    use m_compute_levelset
    use m_ib_patches
    use m_model
    use m_phase_timing, only: s_mem_note, MEM_IB, MEM_IBHOST
    use m_patch_geometries
    use m_collisions
    use m_ib_wall_column
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf, flux_hatR_rsx_vf
    use m_thermochem, only: num_species, gas_constant, molecular_weights, get_mixture_molecular_weight, get_mixture_energy_mass, &
        & get_mixture_thermal_conductivity_mixavg, get_species_mass_diffusivities_mixavg, get_mixture_specific_heat_cp_mass

    use m_surface_thermochem, only: get_surface_net_production_rates, get_surface_reaction_heat_flux
    use m_thermochem_state, only: get_species_enthalpies_mass

    ! Fine-IB AMR state lives in m_ibm_fine (see there)
    use m_ibm_fine

    implicit none

    private :: s_compute_image_points, s_compute_interpolation_coeffs, s_interpolate_image_point, s_find_ghost_points, &
        & s_find_num_ghost_points, s_compute_ghost_point_pressure, s_compute_ghost_point_velocity, s_interpolation_weights, &
        & s_wall_fill, f_tsc_apply, f_quad_fill
    ; public :: s_initialize_ibm_module, s_ibm_setup, s_ibm_correct_state, s_finalize_ibm_module, s_report_ibm_surface, &
        & s_write_ib_surface, s_ib_flux_wall_prepare, s_update_ib_temperatures, s_ibm_alloc_fine, s_ibm_setup_fine, &
        & s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_load_fine_markers, s_ibm_reset_fine_park, s_ibm_note_fine_memory, &
        & num_gps, f_ibm_fine_gps, s_ibm_tally_gps_by_body, s_ibm_check_bodies_on_one_rank, s_conserve_ib_mass, &
        & s_add_swept_cells, s_fill_uncovered_cells, s_ibm_set_conduction_temperature, s_ibm_restore_conduction_temperature, &
        & s_write_ib_conduction_temperature, s_read_ib_conduction_temperature, s_ibm_advance_wall_columns, s_wall_column_tcond, &
        & f_wall_column_dt

    !> Ghost points at which the reacting-surface Newton solve did not reach its tolerance, so the point fell back to a chemically
    !! inert wall (keeping a prescribed Twall). Counted because that fallback is otherwise indistinguishable from a surface
    !! mechanism that simply does nothing, and a run can look converged while the surface chemistry never engaged.
    integer :: n_surface_not_converged = 0

    !> Ghost points whose levelset distance was not positive, leaving the surface gradient (X_IP - X_s)/d undefined. Separate from
    !! the counter above because this is a grid/geometry degeneracy, not a kinetics failure.
    integer :: n_surface_ill_posed = 0

    !> Quadratic wall-normal fill of one ghost point (ib_wall_order = 2): TSC stencils of two probes along the surface normal and
    !! the Lagrange weights that carry the profile through them back to the ghost (see s_wall_fill).
    type :: wall_fill
        integer                   :: tier      !< 1 quadratic, 2 linear mirror (the outer probe lacks fluid), 3 the wall value
        real(wp), dimension(3)    :: lw        !< Lagrange weights on (phi_w, phi_1, phi_2)
        real(wp), dimension(2)    :: llin      !< Linear mirror on (phi_w, phi_1), which bounds the quadratic
        integer, dimension(3)     :: ic1, ic2  !< TSC stencil centres of the probes
        real(wp), dimension(3, 3) :: wd1, wd2  !< TSC weights of the probes
    end type wall_fill

    !> Ghost-point updates per fill tier and AMR level (0 = the coarse grid), summed over the run and reported at its end, because a
    !! tight packing can sit mostly in tier 2 and a second-order claim is then not true of that surface.
    integer, allocatable        :: n_fill_tier(:,:)
    type(integer_field), public :: ib_markers
    $:GPU_DECLARE(create='[ib_markers]')

    !> Body markers in the fine-advance frame, read by s_compute_rhs's body-cell RHS zeroing while amr_in_fine_advance: the advanced
    !! block's fine markers (for a batch, every member's at its slab offset), loaded by s_ibm_load_fine_markers before each fine RHS
    !! pass. ib_markers holds the coarse markers during a fine RHS. Sized to the allocation extents plus one zeroed layer, which the
    !! face loops read (bodies sit inside their blocks).
    type(integer_field), public :: ib_markers_fine
    $:GPU_DECLARE(create='[ib_markers_fine]')

    !> Layout of the last s_ibm_load_fine_markers (member slots, stacking dimension and stride), for the fine conduction temperature
    integer              :: fine_lay_nb = 0, fine_lay_sd = 1, fine_lay_w = 0
    integer, allocatable :: fine_lay_slots(:)

    !> 1 at ghost points whose image point reaches the fluid, 2 at those with a buried image point, 0 elsewhere. The latter average
    !! over the former, and never over each other.
    type(integer_field) :: corrected_gps
    $:GPU_DECLARE(create='[corrected_gps]')

    type(ghost_point), dimension(:), allocatable :: ghost_points

    !> New primitive state per ghost point, staged so it is published only after every image-point stencil has been read
    !! (MFlowCode/MFC#1901). Writing q_prim inside the same loop that reads it at neighboring image-point stencils makes a ghost's
    !! result depend on whether its neighbors were updated first. Inverse-distance weights mask solid cells, so the race is latent
    !! there; a stencil that reads ghost cells (ib_ip_interp = 2) needs it.
    real(wp), dimension(:,:), allocatable :: gp_prim_new
    $:GPU_DECLARE(create='[ghost_points]')
    integer(kind=8), parameter :: tcond_magic = 7164990_8  !< first header word of restart_data/ib_tcond files
    integer(kind=8), parameter :: tcond_version = 1_8      !< format version of the ib_tcond header
    integer, parameter         :: tcond_chunk = 65536      !< records per block read or written of an ib_tcond file

    !> Packed (i, j, k, T_cond) per ghost for the T_cond restart file, so only these four reals per ghost cross to the host
    real(wp), dimension(:,:), allocatable :: gp_tc_rec
    !> per-(i,j)-column and per-i-row ghost-point counts/offsets for the ordered gather
    integer, allocatable, dimension(:,:) :: gp_col
    integer, allocatable, dimension(:)   :: gp_slab
    $:GPU_DECLARE(create='[gp_col, gp_slab]')

    !> Surface record per ghost point from the latest s_ibm_correct_state, kept when ib_surface_wrt: (1) area weight, (2) wall
    !! temperature, (3) gasified mass flux. Only ghost points within two cells of the surface carry a weight, so summing weight*flux
    !! over them integrates over the surface; see s_record_gp_surface.
    integer, parameter                    :: ib_surf_nvars = 3
    real(wp), allocatable, dimension(:,:) :: gp_surf
    $:GPU_DECLARE(create='[gp_surf]')

    !> ib_flux_wall: per ghost point with a fluid face neighbor of a body, from its surface solve (zero for an inert body), the
    !! species mass fluxes into the gas (1:num_species), the energy flux (ibw_e), the blowing velocity (ibw_u:ibw_u + 2) and the
    !! wall area per face area, 1/sum|n_i| (ibw_a, 0 off a wall); ib_wall_theta: each fluid cell's positivity factor on its wall
    !! faces' rates
    integer, parameter, public :: ibw_e = num_species + 1, ibw_u = num_species + 2
    integer, parameter, public :: ibw_a = num_species + 5
    ! Not allocatable: the AMD OpenMP kernels of m_riemann_solver_hllc assume every allocatable they name is present, and these
    ! are only filled under ib_flux_wall (their sf stay unassociated otherwise, as alf_sum's do in m_rhs)
    type(scalar_field), dimension(ibw_a), public :: ib_wall_src
    type(scalar_field), public                   :: ib_wall_theta
    $:GPU_DECLARE(create='[ib_wall_src, ib_wall_theta]')

    !> gp_surf is filled when it is written (ib_surface_wrt) or integrated (ib_exchange).
    logical         :: record_surface
    logical, public :: ib_exchange = .false.   !< ib_heat is accumulated: some IB has thermal_bc = 3, ib_heat_wrt or ib_mass_relax
    logical         :: ib_mass_only = .false.  !< only the mass is accumulated, from fixed bodies without AMR (s_accumulate_ib_mass)
    logical, public :: ib_wall_col = .false.   !< Some IB has wall_column (m_ib_wall_column)
    logical         :: wc_ready = .false.      !< Columns initialized (from the gas, at the first step of a run or restart)
    $:GPU_DECLARE(create='[record_surface]')

    !> Per patch_ib slot (the bodies this rank holds, not global ids), for the thermal_bc = 3 update, ib_heat_wrt and ib_mass_relax:
    !! (1) surface area and (4) gasification rate, from gp_surf; (2) energy into the solid and (3) mass out of it, both summed over
    !! the faces between the body's cells and fluid cells with the RK stage weights, so over a step they are exactly the energy and
    !! mass the fluid exchanged with the body; (5) the gas mass of the body's ib_mass_relax shell; (6) a copy of the body's bank
    !! (patch_ib%mass_bank) and (7) the share of it released this step, for the multi-rank exchange. Rows 1:5 are this rank's
    !! partial sums until s_ib_exchange_heat adds them over the ranks that hold the body, through its owner.
    integer, parameter                    :: ib_nheat = 7
    real(wp), allocatable, dimension(:,:) :: ib_heat
    $:GPU_DECLARE(create='[ib_heat]')

    !> Per patch_ib slot and shard: (1) energy and (2) mass from the face fluxes and (3) shell mass of ib_mass_relax, before they
    !! are summed into ib_heat(2:3) and ib_heat(5). Every cell of a body's surface adds to the same few sums, and on AMD a device
    !! atomic on a double is a compare-and-swap loop that serializes by address, so the adds are spread over ib_nshard sums.
    real(wp), allocatable, dimension(:,:,:) :: ib_part
    integer, parameter                      :: ib_nshard = 32  !< a power of two (f_ib_shard)
    $:GPU_DECLARE(create='[ib_part]')

    !> ib_mass_relax: the bank of a body (the gas mass it gained beyond its gasification and not yet returned) is
    !! patch_ib%mass_bank, held by the owner. ib_shell: patch_ib slot of each fluid cell within ib_shell_width cells of a body, else
    !! 0 (valid until the step's handoff). ib_markers_prev: markers at the last placement of a moving body (s_add_swept_cells).
    type(integer_field) :: ib_shell, ib_markers_prev
    integer, parameter  :: ib_shell_width = 2
    $:GPU_DECLARE(create='[ib_shell, ib_markers_prev]')
    integer :: ib_slot_cap = 0  !< slots of ib_heat and ib_part: a bound on the bodies one rank holds
    integer :: num_gps          !< Number of ghost points
    $:GPU_DECLARE(create='[num_gps]')
    logical :: moving_immersed_boundary_flag
    !> Cells per direction in the image-point interpolation stencil: 2 (inverse distance, bilinear) or 3 (TSC)
    integer :: ib_ip_width
    !> a buried-image-point ghost exists on some rank (set with the interpolation coefficients): the fallback average then needs the
    !! corrected ghosts of the neighbouring ranks, and the halos of corrected_gps and the state are filled before it
    logical :: ib_any_buried = .false.
    !> corrected_gps in the real fields the halo exchange carries
    type(scalar_field), dimension(1) :: corrected_gps_r
    $:GPU_DECLARE(create='[corrected_gps_r]')
    logical :: ib_keep_ghost  !< The interpolation or the wall probes read solid cells, so ghost cells keep the last correction
    logical :: ib_wall_quad   !< ib_wall_order = 2: quadratic wall-normal ghost velocity at no-slip walls
    $:GPU_DECLARE(create='[ib_ip_width, ib_keep_ghost, ib_wall_quad]')
    logical :: centroid_offsets_active = .false.  !< some patch is a moving airfoil or STL, which carries a centroid_offset

    ! IB force reduction (s_communicate_ib_forces): distinct neighborhood ranks in ascending order, and buffers allocated once
    integer, parameter    :: ib_rec_len = 7                           !< gbl_patch_id, force(1:3), torque(1:3)
    integer               :: n_ib_nbrs = 0
    integer, allocatable  :: ib_nbrs(:), ib_send_counts(:), ib_reqs(:), ib_stats(:,:)
    integer, allocatable  :: ib_reply_slots(:,:), ib_reply_counts(:)  !< s_ib_exchange_heat: bodies each neighbor sent partials of
    real(wp), allocatable :: ib_send_bufs(:,:), ib_recv_bufs(:,:)
    private :: ib_rec_len, n_ib_nbrs, ib_nbrs, ib_send_counts, ib_reqs, ib_stats, ib_send_bufs, ib_recv_bufs, ib_reply_slots, &
        & ib_reply_counts

    ! s_handoff_ib_ownership buffers, allocated on its first call
    character(len=1), allocatable :: ho_send_buf(:), ho_recv_bufs(:,:)
    integer, allocatable          :: ho_reqs(:), ho_recv_nbrs(:), ho_gain_gid(:), ho_gain_from(:), ho_gain_slot(:)
    private :: ho_send_buf, ho_recv_bufs, ho_reqs, ho_recv_nbrs, ho_gain_gid, ho_gain_from, ho_gain_slot

    !> Force (1:3) and torque (4:6) on each local IB from the fluid-body face fluxes of the last RHS evaluation
    real(wp), allocatable, dimension(:,:) :: ib_face_ft
    $:GPU_DECLARE(create='[ib_face_ft]')

contains

    !> Allocates memory for the variables in the IBM module
    impure subroutine s_initialize_ibm_module()

        allocate (n_fill_tier(3,0:amr_max_level)); n_fill_tier = 0
        if (amr .and. ib) then
            ! size the declare-target ib_markers for the deepest fine level (never reallocated on Cray); at amr_max_level = 1 the
            ! bounds reduce to the coarse extents
            call s_ibm_marker_bounds()
            @:ALLOCATE(ib_markers%sf(mkr_lo(1):mkr_hi(1), mkr_lo(2):mkr_hi(2), mkr_lo(3):mkr_hi(3)))
            ! corrected_gps spans the same bounds
            @:ALLOCATE(corrected_gps%sf(mkr_lo(1):mkr_hi(1), mkr_lo(2):mkr_hi(2), mkr_lo(3):mkr_hi(3)))
            @:ALLOCATE(ib_markers_fine%sf(-1:m_alloc + 1, -min(n, 1):n_alloc + min(n, 1), -min(p, 1):p_alloc + min(p, 1)))
            @:ACC_SETUP_SFs(ib_markers_fine)
            ib_markers_fine%sf = 0
            $:GPU_UPDATE(device='[ib_markers_fine%sf]')
        else if (p > 0) then
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
            @:ALLOCATE(corrected_gps%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
            @:ALLOCATE(corrected_gps_r(1)%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, -buff_size:p+buff_size))
        else
            @:ALLOCATE(ib_markers%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
            @:ALLOCATE(corrected_gps%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
            @:ALLOCATE(corrected_gps_r(1)%sf(-buff_size:m+buff_size, -buff_size:n+buff_size, 0:0))
        end if

        @:ACC_SETUP_SFs(ib_markers)
        @:ACC_SETUP_SFs(corrected_gps)
        if (associated(corrected_gps_r(1)%sf)) then
            @:ACC_SETUP_SFs(corrected_gps_r(1))
        end if

        $:GPU_ENTER_DATA(copyin='[num_gps]')
        ! Defaults (-1): TSC image points and the quadratic no-slip ghost, except with ib_second_order_vel (own weights and mirror)
        if (ib_ip_interp < 0) ib_ip_interp = merge(1, 2, ib_second_order_vel)
        if (ib_wall_order < 0) ib_wall_order = merge(1, 2, ib_second_order_vel)
        ib_ip_width = merge(3, 2, ib_ip_interp == 2)
        ib_wall_quad = ib_wall_order == 2
        ib_keep_ghost = ib_ip_interp == 2 .or. ib_wall_quad
        $:GPU_UPDATE(device='[ib_second_order_vel, ib_ip_min_dist, ib_ip_interp, ib_ip_width, ib_keep_ghost, ib_wall_quad, ib_slip_ccst]')

    end subroutine s_initialize_ibm_module

    !> Initializes the values of various IBM variables, such as ghost points and image points.
    impure subroutine s_ibm_setup()

        integer                        :: i, j, k, gid, gbl_id
        integer(kind=8)                :: n_need
        real(wp)                       :: t_stage  !< initial time for prescribed kinematics and analytic velocities
        real(wp)                       :: max_num_gps_rank
        integer(kind=8)                :: max_num_gps, fine_gps_est, total_gps, n_pool
        integer(kind=8)                :: n_lumped
        type(ghost_point), allocatable :: tmp_park(:)

        call nvtxStartRange("SETUP-IBM-MODULE")

        ! GPU routines require updated cell centers and boundaries
        $:GPU_UPDATE(device='[num_ibs, num_gbl_ibs, x_cc, y_cc, x_cb, y_cb, dx, dy, ib_bc_x%beg, ib_bc_x%end, ib_bc_y%beg, ib_bc_y%end]')
        if (p /= 0) then
            $:GPU_UPDATE(device='[z_cc, z_cb, dz, ib_bc_z%beg, ib_bc_z%end]')
        end if
        $:GPU_UPDATE(device='[patch_ib(1:num_ibs), glb_bounds]')

        ! do all set up for moving immersed boundaries; prescribed kinematics and analytic velocities are evaluated at the
        ! initial time so the first stage already sees the correct body state, velocity and angular velocity
        t_stage = mytime
        $:GPU_PARALLEL_LOOP(private='[i, gbl_id]', copyin='[t_stage]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                if (patch_ib(i)%kin_model > 0) call s_prescribed_kinematics(i, t_stage)
                if (patch_ib(i)%moving_ibm == 1 .and. patch_ib(i)%kin_model <= 0) then
                    @:mib_analytical()
                end if
            end if
            call s_update_ib_rotation_matrix(i)
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')

        ! allocate some arrays for MPI communication, if required by this simulation
#ifdef MFC_MPI
        if (num_procs > 1) call s_setup_ib_force_comm()
#endif

        @:ALLOCATE(ib_face_ft(6, num_local_ibs_max))
        ib_face_ft = 0._wp
        $:GPU_UPDATE(device='[ib_face_ft]')

        call s_update_ib_lookup()

        ! recompute the new ib_patch locations
        ib_markers%sf = 0._wp
        corrected_gps%sf = 0
        $:GPU_UPDATE(device='[ib_markers%sf, corrected_gps%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        call s_check_every_patch_marked()
        ! One reduction decides whether any patch needs an offset, so cases without one (particle clouds have thousands of
        ! patches) skip the per-patch collectives. Then loop over global ids: ranks hold different patches.
        call s_mpi_allreduce_integer_sum(int(count([(f_needs_centroid_offset(patch_ib(i)), i=1, num_ibs)]), 8), n_need)
        centroid_offsets_active = n_need > 0_8
        if (centroid_offsets_active) then
            do gid = 1, num_gbl_ibs
                call s_get_neighborhood_idx(gid, i)
                call s_compute_centroid_offset(gid, i)
            end do
            call s_restore_centroid_offsets(t_stage)
        end if
        do i = 1, num_ibs
            $:GPU_UPDATE(device='[patch_ib(i)]')
        end do

        ! after the markers and centroid offsets: bodies without an analytic moment take it from their cells' second moment
        call s_sum_ib_inertia_tensors()
        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) call s_compute_moment_of_inertia(patch_ib(i), patch_ib(i)%angular_vel, &
                & patch_ib(i)%moment)
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')

        ! find the number of ghost points and set them to be the maximum total across ranks
        call s_find_num_ghost_points(num_gps)
        if (moving_immersed_boundary_flag) then
            ! Every rank allocates the same size, which must hold any one rank's peak as the body moves: twice the global
            ! total, but at least eight times the largest per-rank count (a thin body crossing few of many ranks)
            call s_mpi_allreduce_integer_sum(int(num_gps, 8), max_num_gps)
            call s_mpi_allreduce_max(real(num_gps, wp), max_num_gps_rank)
            max_num_gps = max(max_num_gps*2_8, 8_8*int(max_num_gps_rank, 8))
            max_num_gps = min(max_num_gps, int(m + 1, 8)*int(n + 1, 8)*int(p + 1, 8))
        else
            max_num_gps = int(num_gps, 8)
        end if
        ! AMR-IB: fine blocks swap their (larger) ghost-point list into this declare-target ghost_points, allocated once (a
        ! realloc corrupts the Cray present table), so size it for the fine list too: fine_gps_cap (buffered cell count) is a volume
        ! bound; refine it by a surface estimate, the whole-body coarse count (global: a seam-spanning block's marker region
        ! reaches neighbor-owned surface) times amr_ref_ratio**(level*(num_dims-1)), x4 margin, floored. The overflow PROHIBITs at
        ! the swap/rebuild sites are the backstop.
        if (allocated(ib_fine)) then
            call s_mpi_allreduce_integer_sum(int(num_gps, 8), total_gps)
            fine_gps_est = min(fine_gps_cap, 4_8*total_gps*int(amr_ref_ratio, 8)**(amr_max_level*(num_dims - 1)) + 4096_8)
            max_num_gps = max(max_num_gps, fine_gps_est)
        end if

        ! set the size of the ghost point arrays to be the amount of points total, plus a factor of 2 buffer
        $:GPU_UPDATE(device='[num_gps]')
        ! no explicit copyin (the device pipeline writes the contents); a dynamic map on top of the declared entry corrupts
        ! CCE-OMP's descriptor (lib-4425)
        @:ALLOCATE(ghost_points(1:max_num_gps))
        @:ALLOCATE(gp_prim_new(1:sys_size, 1:max_num_gps))
        $:GPU_ENTER_DATA(copyin='[gp_prim_new]')
        if (chemistry .and. chem_params%diffusion) then
            @:ALLOCATE(gp_tc_rec(4, 1:max_num_gps))
        end if
        ! Global, not per rank: s_update_ib_temperatures is collective, and a rank whose neighborhood holds no thermal_bc = 3
        ! body must still take part in its allreduce.
        call s_mpi_allreduce_integer_sum(int(count(patch_ib(1:num_ibs)%thermal_bc == 3), 8), n_lumped)
        ib_exchange = n_lumped > 0 .or. ib_heat_wrt .or. ib_mass_relax > 0
        ib_mass_only = ib_mass_relax > 0 .and. n_lumped == 0 .and. .not. (ib_heat_wrt .or. ib_state_wrt .or. ib_force_wrt &
            & .or. moving_immersed_boundary_flag .or. amr .or. hypo_nc_mode == hypo_nc_mode_dual_pass)
        record_surface = ib_surface_wrt .or. ib_exchange
        $:GPU_UPDATE(device='[record_surface]')
        if (record_surface) then
            @:ALLOCATE(gp_surf(ib_surf_nvars, 1:max_num_gps))
            gp_surf = 0._wp
            $:GPU_UPDATE(device='[gp_surf]')
        end if
        if (ib_exchange) then
            ! One slot per body this rank holds (patch_ib slots), with room for the count to change as bodies move between ranks
            call s_mpi_allreduce_max(real(num_ibs, wp), max_num_gps_rank)
            ib_slot_cap = int(max(1_8, min(int(num_gbl_ibs, 8), 2_8*int(max_num_gps_rank, 8) + 64_8)))
            @:ALLOCATE(ib_heat(ib_nheat, ib_slot_cap), ib_part(3, ib_nshard, ib_slot_cap))
            ib_heat = 0._wp
            ib_part = 0._wp
            $:GPU_UPDATE(device='[ib_heat, ib_part]')
        end if
        if (ib_mass_relax > 0) then
            @:ALLOCATE(ib_shell%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(ib_shell)
        end if
        if ((ib_exchange .or. ib_uncover_fill) .and. moving_immersed_boundary_flag) then
            @:ALLOCATE(ib_markers_prev%sf(lbound(ib_markers%sf,1):ubound(ib_markers%sf,1), lbound(ib_markers%sf, &
                       & 2):ubound(ib_markers%sf,2), lbound(ib_markers%sf,3):ubound(ib_markers%sf,3)))
            @:ACC_SETUP_SFs(ib_markers_prev)
            ib_markers_prev%sf = ib_markers%sf
            $:GPU_UPDATE(device='[ib_markers_prev%sf]')
        end if
        @:ALLOCATE(gp_col(0:max(m, mkr_hi(1)) + 1, 0:max(n, mkr_hi(2)) + 1))
        @:ALLOCATE(gp_slab(0:max(m, mkr_hi(1)) + 1))
        ! AMR-IB: device-resident park for the coarse/fine ghost-point swap (move_alloc + GPU_ENTER_DATA, the amr_cg idiom)
        if (allocated(ib_fine)) then
            ! the coarse list plus every fine list a rank can hold. Fixed bodies: lists keep their length between resets, so
            ! segments are exact; size by the surface estimate of all bodies' fine ghost points (block buffers overlap, so x8
            ! not x4), capped at every slot full. Moving bodies: lists change length every stage, so each slot holding any
            ! gets a full-length segment (a body's ghost points fall in one block per level; x2 margin)
            park_coarse = max_num_gps
            if (moving_immersed_boundary_flag) then
                n_pool = park_coarse*(1_8 + min(int(ib_coarse_slot - 1, 8), max(2_8*int(num_gbl_ibs, 8)*int(amr_max_level, 8), &
                                      & 8_8)))
            else
                n_pool = park_coarse + min(int(ib_coarse_slot - 1, 8)*fine_gps_cap, 8_8*total_gps*int(amr_ref_ratio, &
                                           & 8)**(amr_max_level*(num_dims - 1)) + 4096_8)
            end if
            park_top = park_coarse
            allocate (tmp_park(1:n_pool))
            call move_alloc(tmp_park, gp_park)
            $:GPU_ENTER_DATA(create='[gp_park]')
        end if
        ! Ghost-cell IBM, Tseng & Ferziger JCP (2003), Mittal & Iaccarino ARFM (2005)
        call s_find_ghost_points()
        call s_apply_levelset(ghost_points, num_gps)

        call s_compute_image_points()
        call s_compute_interpolation_coeffs()
        if (ib_flux_wall) call s_ib_flux_wall_setup()
        call s_wall_column_setup()

        call nvtxEndRange

    end subroutine s_ibm_setup

    !> Velocity of the body's material point at r from its centroid: translation plus rotation about the centroid. Zero for a fixed
    !! body.
    function f_wall_velocity(gp_patch_id, r) result(vel_w)

        $:GPU_ROUTINE(function_name='f_wall_velocity', parallelism='[seq]', cray_inline=True)
        integer, intent(in)                :: gp_patch_id
        real(wp), dimension(3), intent(in) :: r
        real(wp), dimension(3)             :: vel_w

        vel_w = 0._wp
        if (patch_ib(gp_patch_id)%moving_ibm /= 0) then
            call s_cross_product(patch_ib(gp_patch_id)%angular_vel, r, vel_w)
            vel_w = vel_w + patch_ib(gp_patch_id)%vel
        end if

    end function f_wall_velocity

    !> The wall point nearest a ghost point, relative to the body's centroid: the ghost point is |levelset| inside the wall along
    !! the outward normal.
    function f_wall_point(gp, radial_vector) result(r_w)

        $:GPU_ROUTINE(function_name='f_wall_point', parallelism='[seq]', cray_inline=True)
        type(ghost_point), intent(in)      :: gp
        real(wp), dimension(3), intent(in) :: radial_vector  !< Centroid to ghost point
        real(wp), dimension(3)             :: r_w

        r_w = radial_vector + abs(real(gp%levelset, kind=wp))*gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)

    end function f_wall_point

    !> sum_i kappa_i u_i^2 over the principal directions of the wall, for a velocity u_t tangent to it: the centripetal acceleration
    !! of fluid sliding along a curved wall with outward normal n. Known for a circle and a sphere (1/R in every direction) and the
    !! curved side of a cylinder (none along its axis, the first column of its rotation matrix); zero for any other shape, which
    !! keeps the zero-curvature condition there.
    function f_curvature_u2(gp_patch_id, n, u_t) result(k_u2)

        $:GPU_ROUTINE(function_name='f_curvature_u2', parallelism='[seq]', cray_inline=True)
        integer, intent(in)                :: gp_patch_id
        real(wp), dimension(3), intent(in) :: n, u_t
        real(wp)                           :: k_u2
        real(wp), dimension(3)             :: axis

        k_u2 = 0._wp
        select case (patch_ib(gp_patch_id)%geometry)
        case (ib_circle, ib_sphere)
            k_u2 = sum(u_t**2)/patch_ib(gp_patch_id)%radius
        case (ib_cylinder)
            axis = patch_ib(gp_patch_id)%rotation_matrix(:,1)
            if (abs(dot_product(n, axis)) < 0.5_wp) k_u2 = (sum(u_t**2) - dot_product(u_t, axis)**2)/patch_ib(gp_patch_id)%radius
        end select

    end function f_curvature_u2

    !> Ghost pressure from the wall-normal momentum balance, dp/dn = -rho (a_w . n) + rho sum_i kappa_i u_i^2 with n out of the
    !! body. a_w is the acceleration of the wall point: the body's translation (force/mass for a two-way coupled body, else the
    !! stage-differenced one), alpha x r_w from its angular acceleration, and the centripetal acceleration omega x (omega x r_w) of
    !! its rotation about its centroid. A fluid held at an accelerating wall's velocity needs this pressure gradient to stay with
    !! it; mirroring the image-point pressure instead leaves the wall layer unbalanced. The second term is the centripetal
    !! acceleration of fluid sliding along a curved slip wall, from the image point's velocity relative to the wall; a no-slip wall
    !! holds its fluid at the wall velocity, so it has none. Clamped both ways: the linearization holds only while the correction is
    !! order one, and an unbounded one drives the ghost state to vacuum.
    subroutine s_compute_ghost_point_pressure(gp, gp_patch_id, radial_vector, vel_IP, rho, pres_IP, pres_GP)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in)      :: gp
        integer, intent(in)                :: gp_patch_id
        real(wp), dimension(3), intent(in) :: radial_vector  !< Centroid to ghost point
        real(wp), dimension(3), intent(in) :: vel_IP
        real(wp), intent(in)               :: rho, pres_IP
        real(wp), intent(out)              :: pres_GP
        real(wp), dimension(3)             :: accel, norm, r_w, w_x_r, w_x_w_x_r, u_rel
        real(wp)                           :: a_n

        norm = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        accel = patch_ib(gp_patch_id)%accel
        if (patch_ib(gp_patch_id)%moving_ibm > 1) accel = patch_ib(gp_patch_id)%force/patch_ib(gp_patch_id)%mass
        r_w = f_wall_point(gp, radial_vector)
        call s_cross_product(patch_ib(gp_patch_id)%angular_vel, r_w, w_x_r)
        call s_cross_product(patch_ib(gp_patch_id)%angular_vel, w_x_r, w_x_w_x_r)
        call s_cross_product(patch_ib(gp_patch_id)%angular_accel, r_w, w_x_r)
        a_n = dot_product(accel + w_x_w_x_r + w_x_r, norm)

        if (gp%slip) then
            u_rel = vel_IP - f_wall_velocity(gp_patch_id, r_w)
            u_rel = u_rel - dot_product(u_rel, norm)*norm
            a_n = a_n - f_curvature_u2(gp_patch_id, norm, u_rel)
        end if

        pres_GP = pres_IP/min(max(1._wp - (abs(gp%levelset) + gp%ip_dist)*rho/pres_IP*a_n, 5.e-1_wp), 2._wp)

    end subroutine s_compute_ghost_point_pressure

    subroutine s_compute_ghost_point_velocity(gp, gp_patch_id, radial_vector, vel_IP, pres_IP, dh, vel_GP, q_prim_vf, wf)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in) :: gp
        integer, intent(in) :: gp_patch_id
        real(wp), dimension(3), intent(in) :: radial_vector, vel_IP
        real(wp), intent(in) :: pres_IP
        real(wp), intent(in) :: dh  !< h_IP - h_GP of a slip ghost (ib_slip_ccst); 0 where it does not apply
        real(wp), dimension(3), intent(out) :: vel_GP
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf  !< Probed by the quadratic fill
        type(wall_fill), intent(in) :: wf
        real(wp), dimension(3) :: norm, vel_norm_IP, vel_w, u_t
        real(wp) :: ut2, ut2_g
        real(wp) :: buf, v_blow_eff
        integer :: q  !< Iterator variable

        ! Calculate velocity of ghost cell; the body's motion is evaluated at the wall point, not the ghost point
        vel_w = f_wall_velocity(gp_patch_id, f_wall_point(gp, radial_vector))
        if (gp%slip) then
            norm(1:3) = gp%levelset_norm
            buf = sqrt(sum(norm**2))
            norm = norm/buf
            vel_norm_IP = sum(vel_IP*norm)*norm
            vel_GP = vel_IP - vel_norm_IP
            ! add only the component of the IB's motion that is normal to the surface
            vel_GP = vel_GP + sum(vel_w*norm)*norm
            if (ib_slip_ccst .and. dh /= 0._wp) then
                ! Total enthalpy held across the wall (ghost normal velocity zero relative to it) sets the tangential speed; the
                ! scaling is bounded so a nearly stagnant image point is not given a large velocity
                u_t = vel_IP - vel_w
                u_t = u_t - sum(u_t*norm)*norm
                ut2 = sum(u_t**2)
                ut2_g = ut2 + sum((vel_IP - vel_w)*norm)**2 + 2._wp*dh
                if (ut2 > sgm_eps .and. ut2_g > 0._wp) vel_GP = vel_w + u_t*min(max(sqrt(ut2_g/ut2), 0.5_wp), 2._wp)
            end if
        else
            vel_GP = vel_w
        end if

        ! Burning/injecting surface: superimpose wall-normal (outward) blowing on the
        ! ghost velocity so the immersed surface transpires/injects gas into the flow.
        if (patch_ib(gp_patch_id)%v_blow > 0._wp) then
            v_blow_eff = patch_ib(gp_patch_id)%v_blow
            ! Pressure-coupled burn rate (Vieille's law r_dot ~ p^n)
            if (patch_ib(gp_patch_id)%burn_rate_pref > 0._wp) then
                ! max(pres_IP, 0) guards the fractional power against a transient negative
                v_blow_eff = v_blow_eff*(max(pres_IP, &
                                         & 0._wp)/patch_ib(gp_patch_id)%burn_rate_pref)**patch_ib(gp_patch_id)%burn_rate_exp
            end if
            norm(1:3) = gp%levelset_norm
            buf = sqrt(sum(norm**2))
            if (buf > 0._wp) vel_GP = vel_GP + v_blow_eff*norm/buf
        end if

        if (wf%tier < 3) then
            ! Quadratic profile through the wall value and the two probes. The ghost value is what enforces no-slip in the solve,
            ! so its order sets the order of the wall shear and the drag.
            do q = 1, num_dims
                vel_GP(q) = f_quad_fill(vel_GP(q), f_tsc_apply(q_prim_vf(eqn_idx%mom%beg + q - 1), wf%ic1, wf%wd1), &
                       & f_tsc_apply(q_prim_vf(eqn_idx%mom%beg + q - 1), wf%ic2, wf%wd2), wf%lw, wf%llin)
            end do
        else if (ib_second_order_vel) then
            ! extrapolate linearly through the boundary intercept, Mittal et al. (2008)
            vel_GP = vel_GP - abs(gp%levelset)/max(gp%ip_dist, sgm_eps)*(vel_IP - vel_GP)
        end if

    end subroutine s_compute_ghost_point_velocity

    !> Update the conservative variables at the ghost points
    subroutine s_ibm_correct_state(q_cons_vf, q_prim_vf, pb_in, mv_in, lvl)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf  !< Primitive Variables
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf  !< Primitive Variables
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in) :: lvl  !< AMR level of the frame (0 = the coarse grid), for the fill-tier tally
        integer :: i, j, k, l, q, r  !< Iterator variables
        integer :: jj, kk, ll  !< Neighbor iterators
        integer :: rad, rad_z, num_nbrs  !< Neighbor stencil radius and population
        integer :: patch_id, patch_id_temp  !< Patch ID of ghost point
        real(wp) :: rho, gamma, pi_inf, dyn_pres  !< Mixture variables
        real(wp) :: vel_sum_g, E_ghost  !< Ghost-point velocity magnitude and energy
        real(wp), dimension(2) :: Re_K
        real(wp) :: qv_K
        real(wp) :: pres_IP, pres_GP
        real(wp), dimension(3) :: vel_IP
        real(wp) :: c_IP
        type(wall_fill) :: wf  !< Quadratic wall-normal fill of the current ghost point
        integer :: nf1, nf2, nf3  !< Per-call fill-tier tallies

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_IP, alpha_IP, alpha_rho_GP
            real(wp), dimension(3)  :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(18) :: nmom_IP
            real(wp), dimension(12) :: presb_IP, massv_IP
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_IP, alpha_IP, alpha_rho_GP
            real(wp), dimension(nb)         :: r_IP, v_IP, pb_IP, mv_IP
            real(wp), dimension(nb*nmom)    :: nmom_IP
            real(wp), dimension(nb*nnode)   :: presb_IP, massv_IP
        #:endif
        real(wp), dimension(${NUM_SPECIES}$) :: Ys_IP, Ys_g, Ys_s, W_species
        real(wp) :: alpha_q, alpha_rho_q, e_q
        real(wp) :: T_IP, mw_IP, e_IP  !< Image-point temperature, mixture MW, and mass-specific internal energy (chemistry)
        real(wp) :: T_s, T_g, mw_s, mw_g, rho_s, rho_g, mdot_s, v_stefan, d
        real(wp) :: T_adv  !< Temperature of the advected ghost state: T_s, or the windowed mirror T_g on a moving wall
        real(wp), dimension(${NUM_SPECIES}$) :: mdot_k, h_w  !< ib_flux_wall: species mass fluxes into the gas; species enthalpies
        real(wp) :: e_s, nrm_l  !< ib_flux_wall: energy flux into the gas; length of the ghost's normal
        real(wp) :: h_IP, dh_ccst, p_ratio  !< Enthalpy at the image point, its excess over the slip ghost, (p_GP + B)/(p_IP + B)
        logical :: use_ccst  !< ib_slip_ccst applies to this ghost point
        logical :: surface_converged
        integer :: n_not_converged, n_ill_posed  !< Per-call reacting-surface failure tallies (see the module-level counters)
        real(wp), dimension(3) :: norm  !< Levelset normal at the ghost point, normalized below with buf
        ! Primitive variables at the image point associated with a ghost point, interpolated from surrounding fluid cells.

        real(wp), dimension(3) :: physical_loc   !< Physical loc of GP
        real(wp), dimension(3) :: vel_g          !< Velocity of GP
        real(wp), dimension(3) :: radial_vector  !< vector from centroid to ghost point
        real(wp)               :: nbub
        real(wp)               :: buf, buf_prim
        type(ghost_point)      :: gp

        ! set the Moving IBM interior conservative variables. A solid cell that is a ghost, its image point buried or not, keeps the
        ! value of the last correction when the interpolation reads solid cells (ib_keep_ghost): resetting it to the body velocity,
        ! and without chemistry to a pressure of 1 Pa, would feed those values into the image-point state. Inverse-distance weights
        ! drop solid cells, so they never see it.
        $:GPU_PARALLEL_LOOP(private='[i, j, k, patch_id, rho, patch_id_temp]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    patch_id = ib_markers%sf(j, k, l)
                    if (patch_id /= 0) then
                        call s_decode_patch_periodicity(patch_id, patch_id_temp)
                        call s_get_neighborhood_idx(patch_id_temp, patch_id)
                        if (patch_id > 0 .and. .not. (ib_keep_ghost .and. corrected_gps%sf(j, k, l) /= 0)) then
                            ! skip pressure correction with chemistry to prevent unphysical pressure
                            if (.not. chemistry) q_prim_vf(eqn_idx%E)%sf(j, k, l) = 1._wp
                            rho = 0._wp
                            do i = 1, num_fluids
                                rho = rho + q_prim_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l)
                            end do

                            ! Sets the momentum
                            do i = 1, num_dims
                                q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)*rho
                                q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = patch_ib(patch_id)%vel(i)
                            end do
                        end if  ! patch_id > 0
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        n_not_converged = 0
        n_ill_posed = 0
        nf1 = 0; nf2 = 0; nf3 = 0

        if (num_gps > 0) then
            ! Seed the staging buffer, so a component no ghost point writes survives the flush unchanged.
            $:GPU_PARALLEL_LOOP(private='[i, q, j, k, l, gp]', present='[ghost_points, gp_prim_new]')
            do i = 1, num_gps
                gp = ghost_points(i)
                j = gp%loc(1); k = gp%loc(2); l = gp%loc(3)
                do q = 1, sys_size
                    gp_prim_new(q, i) = q_prim_vf(q)%sf(j, k, l)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            $:GPU_PARALLEL_LOOP(private='[i, physical_loc, dyn_pres, alpha_rho_IP, alpha_IP, alpha_rho_GP, pres_IP, pres_GP, &
                                & vel_IP, vel_g, r_IP, v_IP, pb_IP, mv_IP, nmom_IP, presb_IP, massv_IP, rho, gamma, pi_inf, Re_K, &
                                & gp, norm, buf, radial_vector, j, k, l, q, qv_K, c_IP, nbub, patch_id, Ys_IP, W_species, T_IP, &
                                & mw_IP, e_IP, Ys_g, Ys_s, T_s, T_g, T_adv, mw_s, mw_g, rho_s, rho_g, mdot_s, v_stefan, d, &
                                & surface_converged, vel_sum_g, E_ghost, alpha_q, alpha_rho_q, e_q, mdot_k, h_w, e_s, nrm_l, wf, &
                                & h_IP, dh_ccst, p_ratio, use_ccst]', reduction='[[n_not_converged, n_ill_posed, nf1, nf2, &
                                & nf3]]', reductionOp='[+]', present='[ghost_points]')
            do i = 1, num_gps
                if (record_surface) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, ib_surf_nvars
                        gp_surf(q, i) = 0._wp
                    end do
                end if
                ghost_points(i)%T_cond = T_cond_unset
                gp = ghost_points(i)
                if (.not. gp%interp_valid) cycle
                j = gp%loc(1)
                k = gp%loc(2)
                l = gp%loc(3)
                patch_id = ghost_points(i)%ib_patch_id

                ! Calculate physical location of GP
                physical_loc = [x_cc(j), y_cc(k), 0._wp]
                if (num_dims == 3) physical_loc(3) = z_cc(l)

                ! Interpolate primitive variables at image point associated w/ GP
                if (bubbles_euler .and. .not. qbmm) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP)
                else if (qbmm .and. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP, nmom_IP)
                else if (qbmm .and. .not. polytropic) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, &
                                                   & pb_IP, mv_IP, nmom_IP, pb_in, mv_in, presb_IP, massv_IP)
                else if (chemistry) then
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, Ys_IP=Ys_IP)
                else
                    call s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP)
                end if

                ! Injecting (burning) surface: replace the mirrored ghost composition with pure
                ! injected fuel at the local pressure and the ambient (image-point) temperature.
                ! Setting a consistent injected density here (rather than reusing the heavy ambient
                ! rho) keeps the light fuel at a physical temperature and feeds the surface flame.
                if (chemistry .and. patch_ib(patch_id)%inj_species > 0) then
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    T_IP = pres_IP*mw_IP/(alpha_rho_IP(1)*gas_constant)
                    Ys_IP = 0._wp
                    Ys_IP(patch_ib(patch_id)%inj_species) = 1._wp
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    alpha_rho_IP(1) = pres_IP*mw_IP/(T_IP*gas_constant)
                end if

                ! Thermal and heterogeneous reacting-surface boundary conditions.
                v_stefan = 0._wp
                surface_converged = .false.

                if (chemistry .and. patch_ib(patch_id)%inj_species == 0) then
                    ! Intrinsic gas state at the image point:     rho_IP = (alpha*rho)_IP / alpha_IP
                    call get_mixture_molecular_weight(Ys_IP, mw_IP)
                    T_IP = pres_IP*mw_IP*alpha_IP(1)/(alpha_rho_IP(1)*gas_constant)

                    if (patch_ib(patch_id)%surface_reaction == 1) then
                        ! Heterogeneous reacting surface.
                        d = abs(real(gp%levelset, kind=wp))

                        W_species(1:num_species) = molecular_weights(:)

                        ! d is the gas-side gradient length for every surface flux, so a ghost point sitting exactly on the
                        ! immersed surface (grid-aligned rectangle patch, coincident STL facet) leaves the whole balance
                        ! undefined. Left alone the 1/d produces NaN residuals, and a NaN then defeats the pivot test in
                        ! s_solve_surface_linear_system -- NaN <= epsilon is false -- so the solve reports success and writes
                        ! NaN into the ghost state. Such a point is not solvable, so it is recorded and skipped rather than
                        ! having its distance rescaled to something it is not.
                        if (d > 0._wp) then
                            call s_solve_surface(pres_IP, T_IP, patch_ib(patch_id)%Twall, d, Ys_IP, W_species, merge(1, &
                                                 & patch_ib(patch_id)%thermal_bc, f_holds_twall(patch_ib(patch_id)%thermal_bc)), &
                                                 & Ys_s, T_s, mdot_s, surface_converged)
                            if (.not. surface_converged) n_not_converged = n_not_converged + 1
                        else
                            n_ill_posed = n_ill_posed + 1
                        end if

                        if (surface_converged) then
                            call get_mixture_molecular_weight(Ys_s, mw_s)

                            ! Intrinsic gas density at the reacting surface.
                            rho_s = pres_IP*mw_s/(gas_constant*T_s)
                            if (rho_s > 0._wp) v_stefan = mdot_s/rho_s
                        end if

                        ! ib_flux_wall: the rates this ghost's wall faces carry (zero where the solve failed: an inert, impermeable
                        ! face)
                        if (ib_flux_wall) then
                            if (ib_wall_src(ibw_a)%sf(j, k, l) > 0._wp) then
                                mdot_k = 0._wp
                                e_s = 0._wp
                                if (surface_converged .and. rho_s > 0._wp) then
                                    ! W_k omega_k and the energy they carry, sum_k mdot_k (h_k(T_s) + v_n^2/2), per unit wall area
                                    call get_surface_net_production_rates(rho_s, T_s, Ys_s, mdot_k)
                                    call get_species_enthalpies_mass(T_s, h_w)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, num_species
                                        mdot_k(q) = molecular_weights(q)*mdot_k(q)
                                        e_s = e_s + mdot_k(q)*(h_w(q) + 0.5_wp*v_stefan**2)
                                    end do
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, num_species
                                    ib_wall_src(q)%sf(j, k, l) = mdot_k(q)
                                end do
                                ib_wall_src(ibw_e)%sf(j, k, l) = e_s
                                nrm_l = max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, 3
                                    ib_wall_src(ibw_u + q - 1)%sf(j, k, l) = v_stefan*gp%levelset_norm(q)/nrm_l
                                end do
                            end if
                        end if
                    end if

                    ! Inert surface, or a reacting surface whose solve failed: zero species flux, so the surface composition
                    ! is the image-point composition. Only the chemistry is dropped on failure; a prescribed Twall is kept,
                    ! and thermal_bc = 2, which has no wall temperature without the energy balance, takes a zero gradient.
                    if (.not. surface_converged) then
                        Ys_s(1:num_species) = Ys_IP(1:num_species)
                        T_s = T_IP
                        if (f_holds_twall(patch_ib(patch_id)%thermal_bc)) T_s = patch_ib(patch_id)%Twall
                    end if

                    if (record_surface) call s_record_gp_surface(i, gp, T_s, mdot_s, surface_converged)

                    ! The mirror T_g goes to the wall conduction (T_cond) with only the face bound, no window. The advected ghost is
                    ! wall gas at T_s; a moving wall's lab-frame faces advect it into the fluid, where wall gas cut the heat flux
                    ! below that of the same body held fixed, so it carries the windowed mirror.
                    call s_blend_ghost_temperature(T_IP, T_s, -huge(1._wp), huge(1._wp), f_face_theta(q_prim_vf, ib_markers, gp, &
                                                   & T_IP, T_s), T_g)
                    T_adv = T_s
                    if (patch_ib(patch_id)%moving_ibm /= 0) T_adv = min(max(T_g, T_surface_min), T_surface_max)
                    call s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

                    call get_mixture_molecular_weight(Ys_g, mw_g)
                    alpha_rho_IP(1) = alpha_IP(1)*pres_IP*mw_g/(gas_constant*T_adv)
                end if

                ! Isothermal wall without chemistry: keep p_IP, set the ghost density from the mirrored temperature.
                if (.not. chemistry .and. patch_ib(patch_id)%thermal_bc == 1) then
                    call s_phase_temperature(alpha_rho_IP(1)/alpha_IP(1), pres_IP, 1, T_IP)
                    call s_blend_ghost_temperature(T_IP, patch_ib(patch_id)%Twall, 0._wp, huge(1._wp), 1._wp, T_g)
                    call s_phase_density(pres_IP, T_g, 1, rho_g)
                    alpha_rho_IP(1) = alpha_IP(1)*rho_g
                end if

                call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_K, alpha_IP, alpha_rho_IP, Re_K)

                if (surface_tension) gp_prim_new(eqn_idx%c, i) = c_IP

                ! get the vector that points from the centroid to the ghost
                radial_vector(1) = physical_loc(1) - (patch_ib(patch_id)%x_centroid + real(ghost_points(i)%x_periodicity, &
                              & wp)*(glb_bounds(1)%end - glb_bounds(1)%beg))
                radial_vector(2) = physical_loc(2) - (patch_ib(patch_id)%y_centroid + real(ghost_points(i)%y_periodicity, &
                              & wp)*(glb_bounds(2)%end - glb_bounds(2)%beg))
                radial_vector(3) = 0._wp
                if (num_dims == 3) radial_vector(3) = physical_loc(3) - (patch_ib(patch_id)%z_centroid &
                    & + real(ghost_points(i)%z_periodicity, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg))

                ! set the pressure and density. A fixed no-slip wall keeps the image-point pressure
                dh_ccst = 0._wp
                if (patch_ib(patch_id)%moving_ibm == 0 .and. .not. gp%slip) then
                    pres_GP = pres_IP
                    alpha_rho_GP = alpha_rho_IP
                else
                    call s_compute_ghost_point_pressure(gp, patch_id, radial_vector, vel_IP, rho, pres_IP, pres_GP)

                    ! The adiabatic wall condition T_GP = T_IP the correction is derived from also
                    ! fixes the ghost density: p + B = (n - 1)*cv*rho*T at both points under the one
                    ! temperature leaves each partial density carrying the pressure ratio. The volume
                    ! fractions are untouched, so only rho and qv move with it.
                    h_IP = ((gamma + 1._wp)*pres_IP + pi_inf + qv_K)/rho
                    use_ccst = ib_slip_ccst .and. gp%slip
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, num_fluids
                        p_ratio = (pres_GP + isentrope_B(q))/(pres_IP + isentrope_B(q))
                        alpha_rho_GP(q) = alpha_rho_IP(q)*p_ratio
                        ! Slip wall, ib_slip_ccst: the ghost density follows the isentrope instead; a non-positive p + B (a gap
                        ! closing towards vacuum) keeps the pressure scaling and the copied tangential speed
                        if (use_ccst .and. p_ratio > 0._wp) then
                            alpha_rho_GP(q) = alpha_rho_IP(q)*p_ratio**(1._wp/isentrope_n(q))
                        else
                            use_ccst = .false.
                        end if
                    end do
                    call s_compute_mixture_coefficients(alpha_rho_GP, alpha_IP, rho, gamma, pi_inf, qv_K)
                    if (use_ccst) dh_ccst = h_IP - ((gamma + 1._wp)*pres_GP + pi_inf + qv_K)/rho
                end if
                gp_prim_new(eqn_idx%E, i) = pres_GP

                ! Set q_prim_vf params at GP
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    gp_prim_new(q, i) = alpha_rho_GP(q)
                    gp_prim_new(eqn_idx%adv%beg + q - 1, i) = alpha_IP(q)
                end do

                ! Calculate velocity of ghost cell
                wf%tier = 3
                if (ib_wall_quad .and. .not. gp%slip) then
                    call s_wall_fill(gp, wf)
                    if (wf%tier == 1) nf1 = nf1 + 1
                    if (wf%tier == 2) nf2 = nf2 + 1
                    if (wf%tier == 3) nf3 = nf3 + 1
                end if
                call s_compute_ghost_point_velocity(gp, patch_id, radial_vector, vel_IP, pres_IP, dh_ccst, vel_g, q_prim_vf, wf)

                if (chemistry .and. patch_ib(patch_id)%inj_species == 0 .and. patch_ib(patch_id)%surface_reaction == 1 &
                    & .and. surface_converged) then
                    norm(1:3) = gp%levelset_norm
                    buf = sqrt(sum(norm**2))
                    if (buf > 0._wp) vel_g = vel_g + v_stefan*norm/buf
                end if

                ! Set momentum
                vel_sum_g = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = eqn_idx%mom%beg, eqn_idx%mom%end
                    q_cons_vf(q)%sf(j, k, l) = rho*vel_g(q - eqn_idx%mom%beg + 1)
                    gp_prim_new(q, i) = vel_g(q - eqn_idx%mom%beg + 1)
                    vel_sum_g = vel_sum_g + vel_g(q - eqn_idx%mom%beg + 1)**2._wp
                end do
                dyn_pres = 5.e-1_wp*rho*vel_sum_g

                ! Set continuity and adv vars
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    q_cons_vf(q)%sf(j, k, l) = alpha_rho_GP(q)
                    q_cons_vf(eqn_idx%adv%beg + q - 1)%sf(j, k, l) = alpha_IP(q)
                end do

                ! Set color function
                if (surface_tension) q_cons_vf(eqn_idx%c)%sf(j, k, l) = c_IP

                ! Set Energy
                if (chemistry) then
                    ! Use the reconstructed thermal/species ghost state for an inert
                    ! thermal surface or a converged heterogeneous reacting surface.
                    if (patch_ib(patch_id)%inj_species == 0 .and. (patch_ib(patch_id)%surface_reaction == 0 &
                        & .or. surface_converged)) then

                        call get_mixture_energy_mass(T_adv, Ys_g, e_IP)
                        ghost_points(i)%T_cond = T_g
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            q_cons_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l) = rho*Ys_g(q)
                        end do
                    else
                        ! Ordinary chemistry/injection, or fallback after a failed
                        ! heterogeneous surface solve: retain the image-point state.
                        mw_IP = 0._wp
                        call get_mixture_molecular_weight(Ys_IP, mw_IP)
                        T_IP = pres_GP*mw_IP/(rho*gas_constant)
                        call get_mixture_energy_mass(T_IP, Ys_IP, e_IP)
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            q_cons_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l) = rho*Ys_IP(q)
                        end do
                    end if

                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = rho*e_IP + dyn_pres
                else
                    call s_compute_energy(pres_GP, alpha_rho_GP, alpha_IP, vel_sum_g, E_ghost)
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = E_ghost
                end if
                ! Set bubble vars
                if (bubbles_euler .and. .not. qbmm) then
                    call s_comp_n_from_prim(alpha_IP(1), r_IP, nbub, weight)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*2)%sf(j, k, l) = nbub*r_IP(q)
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*2 + 1)%sf(j, k, l) = nbub*v_IP(q)
                        if (.not. polytropic) then
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4)%sf(j, k, l) = nbub*r_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 1)%sf(j, k, l) = nbub*v_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 2)%sf(j, k, l) = nbub*pb_IP(q)
                            q_cons_vf(eqn_idx%bub%beg + (q - 1)*4 + 3)%sf(j, k, l) = nbub*mv_IP(q)
                        end if
                    end do
                end if

                if (qbmm) then
                    nbub = nmom_IP(1)
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb*nmom
                        q_cons_vf(eqn_idx%bub%beg + q - 1)%sf(j, k, l) = nbub*nmom_IP(q)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        q_cons_vf(eqn_idx%bub%beg + (q - 1)*nmom)%sf(j, k, l) = nbub
                    end do

                    if (.not. polytropic) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = 1, nnode
                                pb_in(j, k, l, r, q) = presb_IP((q - 1)*nnode + r)
                                mv_in(j, k, l, r, q) = massv_IP((q - 1)*nnode + r)
                            end do
                        end do
                    end if
                end if

                if (model_eqns == model_eqns_6eq) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = eqn_idx%int_en%beg, eqn_idx%int_en%end
                        alpha_q = alpha_IP(q - eqn_idx%int_en%beg + 1)
                        alpha_rho_q = alpha_rho_GP(q - eqn_idx%int_en%beg + 1)
                        call s_phase_internal_energy(pres_GP, alpha_q, alpha_rho_q, q - eqn_idx%int_en%beg + 1, e_q)
                        q_cons_vf(q)%sf(j, k, l) = e_q
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! Publish the staged state now that every stencil has been read, so the buried-point average below sees it.
            $:GPU_PARALLEL_LOOP(private='[i, q, j, k, l, gp]', present='[ghost_points, gp_prim_new]')
            do i = 1, num_gps
                gp = ghost_points(i)
                j = gp%loc(1); k = gp%loc(2); l = gp%loc(3)
                do q = 1, sys_size
                    q_prim_vf(q)%sf(j, k, l) = gp_prim_new(q, i)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! The ghosts corrected above and their flags, on the neighbouring ranks, are what the average below reads; collective, so
        ! outside the num_gps test
        if (ib_any_buried) call s_ibm_exchange_gp_halos(q_cons_vf, q_prim_vf)

        if (num_gps > 0) then
            ! A ghost point whose image point is buried in a neighboring IB has no fluid to mirror, so
            ! it takes the average of the ghost points corrected above, growing the stencil until it
            ! reaches one. Those neighbors already carry the wall condition, so the average does too.
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, q, r, jj, kk, ll, gp, rad, rad_z, num_nbrs, buf, buf_prim]', &
                                & present='[ghost_points]')
            do i = 1, num_gps
                gp = ghost_points(i)
                if (gp%interp_valid) cycle
                j = gp%loc(1)
                k = gp%loc(2)
                l = gp%loc(3)

                num_nbrs = 0
                rad = 0
                rad_z = 0
                do while (num_nbrs == 0 .and. rad < gp_layers)
                    rad = rad + 1
                    if (p /= 0) rad_z = rad
                    do jj = j - rad, j + rad
                        do kk = k - rad, k + rad
                            do ll = l - rad_z, l + rad_z
                                if (corrected_gps%sf(jj, kk, ll) == 1) num_nbrs = num_nbrs + 1
                            end do
                        end do
                    end do
                end do

                ! Fully enclosed: nothing to average from, so leave the interior placeholder in place
                if (num_nbrs == 0) cycle

                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, sys_size
                    buf = 0._wp
                    buf_prim = 0._wp
                    do jj = j - rad, j + rad
                        do kk = k - rad, k + rad
                            do ll = l - rad_z, l + rad_z
                                if (corrected_gps%sf(jj, kk, ll) == 1) then
                                    buf = buf + q_cons_vf(q)%sf(jj, kk, ll)
                                    buf_prim = buf_prim + q_prim_vf(q)%sf(jj, kk, ll)
                                end if
                            end do
                        end do
                    end do
                    q_cons_vf(q)%sf(j, k, l) = buf/real(num_nbrs, wp)
                    q_prim_vf(q)%sf(j, k, l) = buf_prim/real(num_nbrs, wp)
                end do

                if (qbmm .and. .not. polytropic) then
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = 1, nnode
                            buf = 0._wp
                            buf_prim = 0._wp
                            do jj = j - rad, j + rad
                                do kk = k - rad, k + rad
                                    do ll = l - rad_z, l + rad_z
                                        if (corrected_gps%sf(jj, kk, ll) == 1) then
                                            buf = buf + pb_in(jj, kk, ll, r, q)
                                            buf_prim = buf_prim + mv_in(jj, kk, ll, r, q)
                                        end if
                                    end do
                                end do
                            end do
                            pb_in(j, k, l, r, q) = buf/real(num_nbrs, wp)
                            mv_in(j, k, l, r, q) = buf_prim/real(num_nbrs, wp)
                        end do
                    end do
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        n_surface_not_converged = n_surface_not_converged + n_not_converged
        n_surface_ill_posed = n_surface_ill_posed + n_ill_posed
        n_fill_tier(:,lvl) = n_fill_tier(:,lvl) + [nf1, nf2, nf3]

    end subroutine s_ibm_correct_state

    !> Wall columns (patch_ib%wall_column): give each first-layer ghost of such a body its gas faces' weights, area times (l(s_f) -
    !! l(s_g))/h with l the steady radial conduction coordinate (f_wc_ell) and s the signed wall distance of the face's gas cell and
    !! of the ghost. Over a circle these sum to its length to 1e-5 at 5-160 cells/R, where the face normals err by 0.3-5 %. The
    !! ghost list of fixed bodies is built once, so the columns stay with their ghosts (and ranks).
    impure subroutine s_wall_column_setup()

        integer(kind=8)        :: n_col
        integer                :: i, c, f, dir, sgn
        integer, dimension(3)  :: nb
        real(wp)               :: kap, sg, sn, area, h, sm
        real(wp), dimension(3) :: cb, top, nrm
        logical                :: off_halo
        type(ghost_point)      :: gp

        call s_mpi_allreduce_integer_sum(int(count(patch_ib(1:num_ibs)%wall_column), 8), n_col)
        ib_wall_col = n_col > 0_8
        if (.not. ib_wall_col) return
        @:PROHIBIT(moving_immersed_boundary_flag, "patch_ib%wall_column requires fixed immersed boundaries")

        ! The columns are held for the ghosts of column bodies only (wc_gp lists them), not for every ghost of the rank
        @:ALLOCATE(wc_slot(max(num_gps, 1)))
        $:GPU_PARALLEL_LOOP(private='[i]', present='[ghost_points]')
        do i = 1, max(num_gps, 1)
            wc_slot(i) = 0
            if (i <= num_gps) then
                if (patch_ib(ghost_points(i)%ib_patch_id)%wall_column .and. ghost_points(i)%interp_valid) wc_slot(i) = 1
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[wc_slot]')
        n_wc = 0
        do i = 1, num_gps
            if (wc_slot(i) == 0) cycle
            n_wc = n_wc + 1
            wc_slot(i) = n_wc
        end do
        @:ALLOCATE(wc_gp(max(n_wc, 1)))
        wc_gp = 0
        do i = 1, num_gps
            if (wc_slot(i) > 0) wc_gp(wc_slot(i)) = i
        end do
        $:GPU_UPDATE(device='[wc_slot, wc_gp]')
        @:ALLOCATE(wc_T(wc_n, max(n_wc, 1)), wc_w(6, max(n_wc, 1)), wc_g(6, max(n_wc, 1)), wc_ext(6, max(n_wc, 1)))
        @:ALLOCATE(wc_q(max(n_wc, 1)), wc_cut(max(n_wc, 1)), wc_on(max(n_wc, 1)))
        off_halo = .false.
        $:GPU_PARALLEL_LOOP(private='[i, c, f, dir, sgn, nb, kap, sg, sn, area, h, sm, cb, top, nrm, gp]', copy='[off_halo]', &
                            & present='[ghost_points]')
        do c = 1, max(n_wc, 1)
            wc_on(c) = -1
            wc_q(c) = 0._wp
            wc_cut(c) = 0._wp
            do f = 1, 6
                wc_w(f, c) = 0._wp
                wc_g(f, c) = 0._wp
                wc_ext(f, c) = 0._wp
            end do
            i = wc_gp(c)
            if (i == 0) cycle
            gp = ghost_points(i)
            kap = 1._wp/patch_ib(gp%ib_patch_id)%radius
            cb = f_wc_body_centre(gp)
            sg = f_wc_dist(gp, cb, f_cell_centre(gp%loc))
            if (1._wp + kap*sg <= 0._wp) cycle
            do dir = 1, num_dims
                do sgn = -1, 1, 2
                    nb = gp%loc
                    nb(dir) = nb(dir) + sgn
                    if (ib_markers%sf(nb(1), nb(2), nb(3)) /= 0) cycle
                    call s_wc_face_metric(gp%loc, dir, nb, area, h)
                    sn = f_wc_dist(gp, cb, f_cell_centre(nb))
                    f = 2*dir - 1 + (sgn + 1)/2
                    wc_w(f, c) = area*(f_wc_ell(kap, sn) - f_wc_ell(kap, sg))/h
                    wc_g(f, c) = area/h
                    wc_on(c) = 0
                end do
            end do
            ! The column-top probe's TSC stencil has to lie inside the halo
            if (wc_on(c) == 0) then
                call s_wc_axis(gp, sm, nrm, top)
                top = top + sm*nrm
                if (top(1) < x_cc(-buff_size + 2) .or. top(1) > x_cc(m + buff_size - 2) .or. top(2) < y_cc(-buff_size + 2) &
                    & .or. top(2) > y_cc(n + buff_size - 2)) off_halo = .true.
                if (p > 0) then
                    if (top(3) < z_cc(-buff_size + 2) .or. top(3) > z_cc(p + buff_size - 2)) off_halo = .true.
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()
        @:DEALLOCATE(wc_slot)
        @:PROHIBIT(off_halo, "patch_ib%wall_column: a column top lies beyond the halo; shorten wall_column_len or use fewer ranks")

    end subroutine s_wall_column_setup

    !> Advance every wall column over dt from the state of the step's first RHS, and set the columns' guard targets. T and p at the
    !! top s_m come through the Lagrange probe (s_lag3_weights); where its stencil has less than half its weight in the gas (a
    !! narrow gap) the ghost keeps the mirror this step. lambda and c_p are tabulated per column, at the ghost's (image-point)
    !! composition, between the extremes of the column. On the first step of a run, restart included, the columns start from the gas
    !! field sampled along their normals.
    impure subroutine s_ibm_advance_wall_columns(q_cons_vf, q_prim_vf, q_T_sf, dt_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        type(scalar_field), intent(in)                      :: q_T_sf
        real(wp), intent(in)                                :: dt_step
        real(wp), dimension(wc_n)                           :: Tc
        real(wp), dimension(wc_ntab)                        :: lam_t, cp_t
        real(wp), dimension(${NUM_SPECIES}$)                :: Ys
        real(wp), dimension(3)                              :: cb, nrm, foot, x, wdt
        real(wp), dimension(3, 3)                           :: wd
        integer, dimension(3)                               :: ic, nb
        real(wp)                                            :: Tw, Tm, pm, mw, kap, sm, ff, Tlo, dTt, prop, qv, acc, rho_g
        integer                                             :: i, c, j, k, f, dir, nsub, a1, a2, a3
        logical                                             :: init
        type(ghost_point)                                   :: gp

        if (n_wc == 0) return
        init = .not. wc_ready
        wc_ready = .true.
        $:GPU_PARALLEL_LOOP(private='[i, c, j, k, f, dir, nsub, a1, a2, a3, Tc, lam_t, cp_t, Ys, cb, nrm, foot, x, wdt, wd, ic, &
                            & nb, Tw, Tm, pm, mw, kap, sm, ff, Tlo, dTt, prop, qv, acc, rho_g, gp]', copyin='[init, dt_step]', present='[ghost_points]')
        do c = 1, n_wc
            if (wc_on(c) < 0) cycle
            i = wc_gp(c)
            gp = ghost_points(i)
            call s_wc_axis(gp, sm, nrm, foot)
            call s_lag3_weights(ib_markers, foot + sm*nrm, ic, wd, ff)
            wc_on(c) = 0
            wc_q(c) = 0._wp
            if (ff < 0.5_wp) cycle
            wc_on(c) = 1
            Tw = patch_ib(gp%ib_patch_id)%Twall
            kap = 1._wp/patch_ib(gp%ib_patch_id)%radius
            Tm = f_tsc_apply(q_T_sf, ic, wd)
            pm = f_tsc_apply(q_prim_vf(eqn_idx%E), ic, wd)
            rho_g = 0._wp
            do k = eqn_idx%cont%beg, eqn_idx%cont%end
                rho_g = rho_g + q_cons_vf(k)%sf(gp%loc(1), gp%loc(2), gp%loc(3))
            end do
            do k = 1, num_species
                Ys(k) = q_cons_vf(eqn_idx%species%beg + k - 1)%sf(gp%loc(1), gp%loc(2), gp%loc(3))/rho_g
            end do
            call get_mixture_molecular_weight(Ys, mw)

            do k = 1, wc_n
                if (init) then
                    call s_tsc_weights(ib_markers, foot + f_wc_node(k, sm)*nrm, ic, wd, ff)
                    Tc(k) = Tw
                    if (ff > 0._wp) Tc(k) = f_tsc_fluid_mean(q_T_sf, ib_markers, ic, wd)
                else
                    Tc(k) = wc_T(k, c)
                end if
            end do

            ! Property tables over the column's range
            Tlo = min(Tw, Tm)
            dTt = max(Tw, Tm)
            do k = 1, wc_n
                Tlo = min(Tlo, Tc(k))
                dTt = max(dTt, Tc(k))
            end do
            dTt = max(dTt - Tlo, 1.e-3_wp)/real(wc_ntab - 1, wp)
            do j = 1, wc_ntab
                call get_mixture_thermal_conductivity_mixavg(Tlo + real(j - 1, wp)*dTt, Ys, prop)
                lam_t(j) = prop
                call get_mixture_specific_heat_cp_mass(Tlo + real(j - 1, wp)*dTt, Ys, prop)
                cp_t(j) = prop
            end do

            call s_wc_step(Tc, Tw, Tm, pm*mw/gas_constant, kap, sm, dt_step, Tlo, dTt, lam_t, cp_t, qv)
            wc_q(c) = qv
            do k = 1, wc_n
                wc_T(k, c) = Tc(k)
            end do

            ! Guard targets: the column averaged over each gas face cell, on 4^d sub-points of its true wall distance
            cb = f_wc_body_centre(gp)
            nsub = 4
            do f = 1, 2*num_dims
                if (wc_g(f, c) <= 0._wp) cycle
                dir = (f + 1)/2
                nb = gp%loc
                nb(dir) = nb(dir) + 2*mod(f + 1, 2) - 1
                wdt = 0._wp
                wdt(1) = dx(nb(1)); wdt(2) = dy(nb(2))
                if (p > 0) wdt(3) = dz(nb(3))
                acc = 0._wp
                do a1 = 1, nsub
                    do a2 = 1, nsub
                        do a3 = 1, merge(nsub, 1, p > 0)
                            x = f_cell_centre(nb)
                            x(1) = x(1) + ((real(a1, wp) - 0.5_wp)/real(nsub, wp) - 0.5_wp)*wdt(1)
                            x(2) = x(2) + ((real(a2, wp) - 0.5_wp)/real(nsub, wp) - 0.5_wp)*wdt(2)
                            x(3) = x(3) + ((real(a3, wp) - 0.5_wp)/real(nsub, wp) - 0.5_wp)*wdt(3)
                            acc = acc + f_wc_target(Tc, Tw, sm, f_wc_dist(gp, cb, x))
                        end do
                    end do
                end do
                wc_ext(f, c) = acc/real(nsub**num_dims, wp)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_advance_wall_columns

    !> Conduction temperature of each active column ghost that makes its gas faces carry the column's wall flux, less a one-sided
    !! guard: sum_f G_f (T_f - T_cond) = sum_f (q w_f - cut_f), G_f = lambda_f area/h with the face conductivity the diffusion flux
    !! forms (each side at its own state, the ghost at T_w), and cut_f = g lambda (T_ext - T_f)^+ area/h. The guard stops the drain
    !! of a gas cell colder than the column over its own area: heat of gas the staircase has no room for, zero in steady state. The
    !! body still books q w_f (s_accumulate_ib_face_fluxes adds the cut back). Called after each ghost correction, so the gas values
    !! lag one stage as the image point's do.
    impure subroutine s_wall_column_tcond(q_cons_vf, q_prim_vf, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        type(scalar_field), intent(in)                      :: q_T_sf
        real(wp), dimension(${NUM_SPECIES}$)                :: Ys_g, Ys_f
        integer, dimension(3)                               :: nb
        real(wp)                                            :: rho, Tw, lam_g, lam_f, lam_e, Tf, Te, G, num, den, drain, cut
        integer                                             :: i, c, f, q, dir
        type(ghost_point)                                   :: gp

        if (n_wc == 0) return
        $:GPU_PARALLEL_LOOP(private='[i, c, f, q, dir, nb, Ys_g, Ys_f, rho, Tw, lam_g, lam_f, lam_e, Tf, Te, G, num, den, drain, &
                            & cut, gp]', present='[ghost_points]')
        do c = 1, n_wc
            if (wc_on(c) /= 1) cycle
            i = wc_gp(c)
            gp = ghost_points(i)
            rho = 0._wp
            do q = eqn_idx%cont%beg, eqn_idx%cont%end
                rho = rho + q_cons_vf(q)%sf(gp%loc(1), gp%loc(2), gp%loc(3))
            end do
            do q = 1, num_species
                Ys_g(q) = q_cons_vf(eqn_idx%species%beg + q - 1)%sf(gp%loc(1), gp%loc(2), gp%loc(3))/rho
            end do
            Tw = patch_ib(gp%ib_patch_id)%Twall
            call get_mixture_thermal_conductivity_mixavg(Tw, Ys_g, lam_g)
            num = 0._wp; den = 0._wp; drain = 0._wp; cut = 0._wp
            do f = 1, 2*num_dims
                if (wc_g(f, c) <= 0._wp) cycle
                dir = (f + 1)/2
                nb = gp%loc
                nb(dir) = nb(dir) + 2*mod(f + 1, 2) - 1
                Tf = q_T_sf%sf(nb(1), nb(2), nb(3))
                do q = 1, num_species
                    Ys_f(q) = q_prim_vf(eqn_idx%species%beg + q - 1)%sf(nb(1), nb(2), nb(3))
                end do
                call get_mixture_thermal_conductivity_mixavg(Tf, Ys_f, lam_f)
                G = 0.5_wp*(lam_f + lam_g)*wc_g(f, c)
                num = num + G*Tf
                den = den + G
                drain = drain + wc_q(c)*wc_w(f, c)
                Te = wc_ext(f, c)
                if (Te > Tf) then
                    call get_mixture_thermal_conductivity_mixavg(Te, Ys_f, lam_e)
                    cut = cut + wc_guard*0.5_wp*(lam_f + lam_e)*(Te - Tf)*wc_g(f, c)
                end if
            end do
            ghost_points(i)%T_cond = (num - drain + cut)/den
            wc_cut(c) = cut
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_wall_column_tcond

    !> A wall column's length s_m = max(3h, L) (L = wall_column_len, default 0.15 R), unit normal and wall foot
    subroutine s_wc_axis(gp, sm, nrm, foot)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(ghost_point), intent(in)       :: gp
        real(wp), intent(out)               :: sm
        real(wp), dimension(3), intent(out) :: nrm, foot
        real(wp)                            :: L

        L = patch_ib(gp%ib_patch_id)%wall_column_len
        if (.not. L > 0._wp) L = wc_len_frac*patch_ib(gp%ib_patch_id)%radius
        sm = max(3._wp*f_cell_scale(gp), L)
        nrm = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        foot = f_cell_centre(gp%loc) + abs(real(gp%levelset, kind=wp))*nrm

    end subroutine s_wc_axis

    !> Centre of the periodic image of a ghost point's circle or sphere
    function f_wc_body_centre(gp) result(cb)

        $:GPU_ROUTINE(function_name='f_wc_body_centre', parallelism='[seq]', cray_inline=True)
        type(ghost_point), intent(in) :: gp
        real(wp), dimension(3)        :: cb

        cb(1) = patch_ib(gp%ib_patch_id)%x_centroid + real(gp%x_periodicity, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
        cb(2) = patch_ib(gp%ib_patch_id)%y_centroid + real(gp%y_periodicity, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
        cb(3) = 0._wp
        if (p > 0) cb(3) = patch_ib(gp%ib_patch_id)%z_centroid + real(gp%z_periodicity, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)

    end function f_wc_body_centre

    !> Signed wall distance of point x from a ghost point's circle or sphere (centre cb)
    function f_wc_dist(gp, cb, x) result(s)

        $:GPU_ROUTINE(function_name='f_wc_dist', parallelism='[seq]', cray_inline=True)
        type(ghost_point), intent(in)      :: gp
        real(wp), dimension(3), intent(in) :: cb, x
        real(wp)                           :: s

        s = sqrt(sum((x - cb)**2)) - patch_ib(gp%ib_patch_id)%radius

    end function f_wc_dist

    !> Centre of cell c (z = 0 in 2D)
    function f_cell_centre(c) result(x)

        $:GPU_ROUTINE(function_name='f_cell_centre', parallelism='[seq]', cray_inline=True)
        integer, dimension(3), intent(in) :: c
        real(wp), dimension(3)            :: x

        x(1) = x_cc(c(1))
        x(2) = y_cc(c(2))
        x(3) = 0._wp
        if (p > 0) x(3) = z_cc(c(3))

    end function f_cell_centre

    !> Area of the face between cells c and nb (neighbours along dir) and the distance between their centres
    subroutine s_wc_face_metric(c, dir, nb, area, h)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, dimension(3), intent(in) :: c, nb
        integer, intent(in)               :: dir
        real(wp), intent(out)             :: area, h

        select case (dir)
        case (1)
            area = dy(c(2))
            h = abs(x_cc(nb(1)) - x_cc(c(1)))
        case (2)
            area = dx(c(1))
            h = abs(y_cc(nb(2)) - y_cc(c(2)))
        case default
            area = dx(c(1))*dy(c(2))
            h = abs(z_cc(nb(3)) - z_cc(c(3)))
        end select
        if (p > 0 .and. dir /= 3) area = area*dz(c(3))

    end subroutine s_wc_face_metric

    !> Wall-column time-step limit: a step's column drain takes at most wall_column_dt_frac of a near-wall gas cell's heat span, rho
    !! c_v V |T_top - T_w| (rho, c_v at the cell's state). At a step start the column's flux is the true singular one, about 10x the
    !! mirror's, and at the acoustic step it removed 40 % of a cell's heat per wall face in one stage: the collapsing pressure sent
    !! an impulsively started flow non-finite. From the last step's flux; before the first, from the semi-infinite start, q dt = 2
    !! lambda dT sqrt(dt/(pi alpha_w)). Idle in steady state, where the drain per step is a small part of the span.
    impure function f_wall_column_dt(q_prim_vf, q_T_sf) result(dt_lim)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(in)                      :: q_T_sf
        real(wp)                                            :: dt_lim
        real(wp), dimension(${NUM_SPECIES}$)                :: Ys
        integer, dimension(3)                               :: nb
        real(wp)                                            :: Tw, Tf, pf, mw, cp_f, cp_w, lam_f, lam_w, rcv, sw, span, vol, lim
        integer                                             :: i, c, f, q, dir
        real(wp)                                            :: frac
        logical                                             :: ready
        type(ghost_point)                                   :: gp

        dt_lim = huge(1._wp)
        if (n_wc == 0) return
        ready = wc_ready
        frac = wall_column_dt_frac
        $:GPU_PARALLEL_LOOP(private='[i, c, f, q, dir, nb, Ys, Tw, Tf, pf, mw, cp_f, cp_w, lam_f, lam_w, rcv, sw, span, vol, lim, &
                            & gp]', copyin='[ready, frac]', reduction='[[dt_lim]]', reductionOp='[min]', present='[ghost_points]')
        do c = 1, n_wc
            if (wc_on(c) < 0 .or. (ready .and. wc_on(c) /= 1)) cycle
            i = wc_gp(c)
            gp = ghost_points(i)
            ! the column's total face weight, and the gas cell of its last gas face
            sw = 0._wp
            nb = gp%loc
            do f = 1, 2*num_dims
                if (wc_g(f, c) <= 0._wp) cycle
                sw = sw + wc_w(f, c)
                dir = (f + 1)/2
                nb = gp%loc
                nb(dir) = nb(dir) + 2*mod(f + 1, 2) - 1
            end do
            Tf = q_T_sf%sf(nb(1), nb(2), nb(3))
            pf = q_prim_vf(eqn_idx%E)%sf(nb(1), nb(2), nb(3))
            do q = 1, num_species
                Ys(q) = q_prim_vf(eqn_idx%species%beg + q - 1)%sf(nb(1), nb(2), nb(3))
            end do
            call get_mixture_molecular_weight(Ys, mw)
            call get_mixture_specific_heat_cp_mass(Tf, Ys, cp_f)
            Tw = patch_ib(gp%ib_patch_id)%Twall
            vol = dx(nb(1))*dy(nb(2))
            if (p > 0) vol = vol*dz(nb(3))
            rcv = pf*mw/(gas_constant*Tf)*(cp_f - gas_constant/mw)*vol
            if (ready) then
                span = abs(wc_T(wc_n, c) - Tw)
                lim = frac*rcv*span/max(abs(wc_q(c))*sw, tiny(1._wp))
            else
                ! q dt = 2 lambda_bar dT sqrt(dt/(pi alpha_w)) <= frac rcv dT / sw, alpha_w = lambda_w R T_w/(p W c_p,w)
                call get_mixture_thermal_conductivity_mixavg(Tf, Ys, lam_f)
                call get_mixture_thermal_conductivity_mixavg(Tw, Ys, lam_w)
                call get_mixture_specific_heat_cp_mass(Tw, Ys, cp_w)
                lim = (frac*rcv*sqrt(pi*lam_w*gas_constant*Tw/(pf*mw*cp_w))/((lam_f + lam_w)*max(sw, tiny(1._wp))))**2
            end if
            dt_lim = min(dt_lim, lim)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end function f_wall_column_dt

    !> Smallest cell width at a ghost point's cell
    function f_cell_scale(gp) result(h)

        $:GPU_ROUTINE(function_name='f_cell_scale', parallelism='[seq]', cray_inline=True)
        type(ghost_point), intent(in) :: gp
        real(wp)                      :: h

        h = min(dx(gp%loc(1)), dy(gp%loc(2)))
        if (p > 0) h = min(h, dz(gp%loc(3)))

    end function f_cell_scale

    !> Put the mirrored ghost temperature T_cond from s_ibm_correct_state in q_T_sf, where the wall conduction takes its face
    !! difference. Only that difference may read it: the mirror of a cold wall in hot gas lies below 0 K, so the diffusion
    !! properties use each ghost's own state, and s_ibm_restore_conduction_temperature undoes this once the fluxes are formed. Call
    !! after the conservative-to-primitive conversion, before the halo exchange.
    subroutine s_ibm_set_conduction_temperature(q_T_sf)

        type(scalar_field), intent(inout) :: q_T_sf
        integer                           :: i

        if (amr_in_fine_advance) then
            call s_ibm_fine_conduction_temperature(q_T_sf, .true.)
            return
        end if
        if (num_gps == 0) return
        $:GPU_PARALLEL_LOOP(private='[i]', present='[ghost_points]')
        do i = 1, num_gps
            if (ghost_points(i)%T_cond > T_cond_unset) then
                ghost_points(i)%T_save = q_T_sf%sf(ghost_points(i)%loc(1), ghost_points(i)%loc(2), ghost_points(i)%loc(3))
                q_T_sf%sf(ghost_points(i)%loc(1), ghost_points(i)%loc(2), ghost_points(i)%loc(3)) = ghost_points(i)%T_cond
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_set_conduction_temperature

    !> Undo s_ibm_set_conduction_temperature on this rank's ghosts, so the reaction source and the split chemistry substep see each
    !! ghost's own temperature
    subroutine s_ibm_restore_conduction_temperature(q_T_sf)

        type(scalar_field), intent(inout) :: q_T_sf
        integer                           :: i

        if (amr_in_fine_advance) then
            call s_ibm_fine_conduction_temperature(q_T_sf, .false.)
            return
        end if
        if (num_gps == 0) return
        $:GPU_PARALLEL_LOOP(private='[i]', present='[ghost_points]')
        do i = 1, num_gps
            if (ghost_points(i)%T_cond > T_cond_unset) then
                q_T_sf%sf(ghost_points(i)%loc(1), ghost_points(i)%loc(2), ghost_points(i)%loc(3)) = ghost_points(i)%T_save
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_restore_conduction_temperature

    !> The set (set = T) or restore of the conduction temperature on an AMR fine frame: the coarse ghost list stays active there, so
    !! each member block's own parked ghosts (gp_park) are used, at its slab offset from s_ibm_load_fine_markers.
    impure subroutine s_ibm_fine_conduction_temperature(q_T_sf, set)

        type(scalar_field), intent(inout) :: q_T_sf
        logical, intent(in)               :: set
        integer                           :: ibm, a, n_f, o(3), o1, o2, o3, j, k, l
        integer(kind=8)                   :: off

        do ibm = 1, fine_lay_nb
            n_f = ib_fine(fine_lay_slots(ibm))%num_gps
            if (n_f == 0) cycle
            off = ib_fine(fine_lay_slots(ibm))%off
            o = 0; o(fine_lay_sd) = (ibm - 1)*fine_lay_w
            o1 = o(1); o2 = o(2); o3 = o(3)
            $:GPU_PARALLEL_LOOP(private='[a, j, k, l]', firstprivate='[n_f, off, o1, o2, o3, set]')
            do a = 1, n_f
                if (gp_park(off + a)%T_cond > T_cond_unset) then
                    j = gp_park(off + a)%loc(1) + o1
                    k = gp_park(off + a)%loc(2) + o2
                    l = gp_park(off + a)%loc(3) + o3
                    if (set) then
                        gp_park(off + a)%T_save = q_T_sf%sf(j, k, l)
                        q_T_sf%sf(j, k, l) = gp_park(off + a)%T_cond
                    else
                        q_T_sf%sf(j, k, l) = gp_park(off + a)%T_save
                    end if
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_ibm_fine_conduction_temperature

    !> Pack every ghost's global cell index (i, j, k) and T_cond into gp_tc_rec(1:4, 1:num_gps) on the device and bring only those
    !! four reals per ghost to the host, instead of the whole ghost_point (about 380 B).
    impure subroutine s_pack_ghost_tcond()

        integer :: i, s(3), s1, s2, s3

        if (num_gps == 0) return
        s = 0
        s(1:num_dims) = start_idx(1:num_dims)
        s1 = s(1); s2 = s(2); s3 = s(3)
        $:GPU_PARALLEL_LOOP(private='[i]', firstprivate='[s1, s2, s3]', present='[ghost_points, gp_tc_rec]')
        do i = 1, num_gps
            gp_tc_rec(1, i) = real(ghost_points(i)%loc(1) + s1, wp)
            gp_tc_rec(2, i) = real(ghost_points(i)%loc(2) + s2, wp)
            gp_tc_rec(3, i) = real(ghost_points(i)%loc(3) + s3, wp)
            gp_tc_rec(4, i) = ghost_points(i)%T_cond
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_UPDATE(host='[gp_tc_rec(:, 1:num_gps)]')

    end subroutine s_pack_ghost_tcond

    !> Send gp_tc_rec(4, :) to the T_cond of each ghost: a component copy, never the whole ghost_point.
    impure subroutine s_unpack_ghost_tcond()

        integer :: i

        if (num_gps == 0) return
        $:GPU_UPDATE(device='[gp_tc_rec(:, 1:num_gps)]')
        $:GPU_PARALLEL_LOOP(private='[i]', present='[ghost_points, gp_tc_rec]')
        do i = 1, num_gps
            ghost_points(i)%T_cond = gp_tc_rec(4, i)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_unpack_ghost_tcond

    !> Path of the T_cond file of a save: one shared file, or one file per rank with file_per_process (as ib_state).
    impure function f_tcond_file(t_step) result(file_loc)

        integer, intent(in)                  :: t_step
        character(len=path_len + 2*name_len) :: file_loc

        if (file_per_process) then
            write (file_loc, '(A,I0,A,I0,A,I7.7,A)') trim(case_dir) // '/restart_data/lustre_', t_step, '/ib_tcond_', t_step, &
                   & '_', proc_rank, '.dat'
        else
            write (file_loc, '(A,I0,A)') trim(case_dir) // '/restart_data/ib_tcond_', t_step, '.dat'
        end if

    end function f_tcond_file

    !> Write every ghost's T_cond to restart_data/ib_tcond_<t_step>.dat (lustre_<t_step>/ib_tcond_<t_step>_<rank>.dat with
    !! file_per_process). The file is a header of 64-bit integers (magic, format version, bytes per real, number of slabs n_w, then
    !! per slab its record count and the bounding box of its global cell indices) followed by one slab per writing rank, each record
    !! four reals: the global cell index (i, j, k) and T_cond of a ghost that has one. Keyed by cell, so the file does not depend on
    !! the order the ghosts were found in or on the rank count. Each rank writes its own slab at an offset from an exclusive scan of
    !! the counts: no rank holds the global records, and only 7 integers per rank reach rank 0 for the header. T_cond is formed from
    !! the primitive state at the start of the last stage, so it cannot be rebuilt from the saved conservative state. Collective.
    impure subroutine s_write_ib_conduction_temperature(t_step)

        integer, intent(in)                  :: t_step
        character(len=path_len + 2*name_len) :: file_loc
        character(len=10)                    :: t_step_string
        integer(kind=8)                      :: info(7), off
        integer(kind=8), allocatable         :: tab(:,:), hdr(:)
        integer                              :: i, n_loc, n_w
        logical                              :: lead

#ifdef MFC_MPI
        integer(kind=MPI_OFFSET_KIND)       :: base, rb
        integer                             :: ifile, ierr, k0, nc
        integer, dimension(MPI_STATUS_SIZE) :: status
#else
        integer :: file_unit
#endif

        call s_pack_ghost_tcond()
        n_loc = 0
        info = 0
        do i = 1, num_gps
            if (gp_tc_rec(4, i) <= T_cond_unset) cycle
            n_loc = n_loc + 1
            gp_tc_rec(:,n_loc) = gp_tc_rec(:,i)
            if (n_loc == 1) then
                info(2:4) = nint(gp_tc_rec(1:3,i), 8)
                info(5:7) = info(2:4)
            else
                info(2:4) = min(info(2:4), nint(gp_tc_rec(1:3,i), 8))
                info(5:7) = max(info(5:7), nint(gp_tc_rec(1:3,i), 8))
            end if
        end do
        info(1) = n_loc

        n_w = merge(1, num_procs, file_per_process)
        lead = file_per_process .or. proc_rank == 0
        allocate (tab(7, merge(n_w, 1, lead)))
        off = 0
        tab(:,1) = info
#ifdef MFC_MPI
        if (.not. file_per_process) then
            call MPI_GATHER(info, 7, MPI_INTEGER8, tab, 7, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
            call MPI_EXSCAN(int(n_loc, 8), off, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) off = 0
        end if
#endif

        if (file_per_process) then
            call s_int_to_str(t_step, t_step_string)
            if (proc_rank == 0) call s_create_directory(trim(case_dir) // '/restart_data/lustre_' // trim(t_step_string))
        else if (proc_rank == 0) then
            call s_create_directory(trim(case_dir) // '/restart_data')
        end if
        call s_mpi_barrier()
        if (file_per_process) call s_delay_file_access(proc_rank)
        file_loc = f_tcond_file(t_step)

        allocate (hdr(4 + 7*n_w))
        if (lead) then
            hdr(1:4) = [tcond_magic, tcond_version, int(storage_size(0._wp)/8, 8), int(n_w, 8)]
            hdr(5:) = reshape(tab, [7*n_w])
        end if

#ifdef MFC_MPI
        call MPI_FILE_OPEN(merge(MPI_COMM_SELF, MPI_COMM_WORLD, file_per_process), file_loc, ior(MPI_MODE_WRONLY, &
                           & MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
        call s_check_mpi_file_open(ierr, file_loc)
        call MPI_FILE_SET_SIZE(ifile, 0_MPI_OFFSET_KIND, ierr)
        rb = int(storage_size(0._wp)/8, MPI_OFFSET_KIND)*4
        base = size(hdr, kind=MPI_OFFSET_KIND)*8 + int(off, MPI_OFFSET_KIND)*rb
        if (lead) call MPI_FILE_WRITE_AT(ifile, 0_MPI_OFFSET_KIND, hdr, size(hdr), MPI_INTEGER8, status, ierr)
        do k0 = 1, n_loc, tcond_chunk
            nc = min(tcond_chunk, n_loc - k0 + 1)
            call MPI_FILE_WRITE_AT(ifile, base + int(k0 - 1, MPI_OFFSET_KIND)*rb, gp_tc_rec(:,k0:k0 + nc - 1), 4*nc, mpi_p, &
                                   & status, ierr)
        end do
        call MPI_FILE_CLOSE(ifile, ierr)
#else
        open (newunit=file_unit, file=trim(file_loc), form='unformatted', access='stream', status='replace')
        write (file_unit) hdr
        write (file_unit) gp_tc_rec(:,1:n_loc)
        close (file_unit)
#endif

    end subroutine s_write_ib_conduction_temperature

    !> Restore T_cond of this rank's ghosts from the file of s_write_ib_conduction_temperature, after s_ibm_setup has found them.
    !! Rank 0 reads the header and broadcasts it (with file_per_process each rank reads its own); a rank then reads only the slabs
    !! whose bounding box meets its cells, in blocks of tcond_chunk records, keeping the records it owns. On the rank count that
    !! wrote the file that is its own slab: O(local ghosts) reads plus O(n_w) header entries; on another rank count a few
    !! neighbouring slabs. Nothing is proportional to the global ghost count. A file from before T_cond was saved is absent: the
    !! first stage then conducts through each ghost's own temperature, as it always did.
    impure subroutine s_read_ib_conduction_temperature(t_step)

        integer, intent(in)                     :: t_step
        real(wp), allocatable, dimension(:,:,:) :: tc
        real(wp), allocatable, dimension(:,:)   :: buf
        character(len=path_len + 2*name_len)    :: file_loc
        integer(kind=8)                         :: head(4), pre, hbytes, rb, fsize, k0, cnt
        integer(kind=8), allocatable            :: tab(:,:)
        integer                                 :: i, w, r, nc, ios, file_unit, c(3), lo(3), hi(3)
        logical                                 :: file_exist, rd

#ifdef MFC_MPI
        integer :: ierr
#endif

        file_loc = f_tcond_file(t_step)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            if (proc_rank == 0) print '(A)', &
                & 'Warning: ' // trim(file_loc) // ' not found; the first stage conducts through the ghost temperature'
            return
        end if

        open (newunit=file_unit, file=trim(file_loc), form='unformatted', access='stream', status='old', action='read')
        rd = file_per_process .or. proc_rank == 0
        head = 0
        ios = 0
        if (rd) read (file_unit, iostat=ios, pos=1) head
        if (ios /= 0) head = 0
#ifdef MFC_MPI
        if (.not. file_per_process) call MPI_BCAST(head, 4, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
#endif
        rb = int(storage_size(0._wp)/8, 8)*4
        if (head(1) /= tcond_magic .or. head(2) /= tcond_version .or. head(3)*4 /= rb .or. head(4) < 1) then
            call s_mpi_abort(trim(file_loc) // ' is not an ib_tcond file of this format (an old file without a header, a ' &
                             & // 'different real precision, or a corrupt file): delete it to restart without the wall ' &
                             & // 'conduction temperature, or rerun the save that wrote it')
        end if

        allocate (tab(7, head(4)))
        if (rd) read (file_unit, iostat=ios, pos=33) tab
        if (rd .and. ios /= 0) call s_mpi_abort('Error reading the header of ' // trim(file_loc))
#ifdef MFC_MPI
        if (.not. file_per_process) call MPI_BCAST(tab, 7*int(head(4)), MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
#endif
        hbytes = 8*(4 + 7*head(4))
        inquire (UNIT=file_unit, SIZE=fsize)
        if (fsize /= hbytes + sum(tab(1,:))*rb) call s_mpi_abort(trim(file_loc) // ' is truncated or does not match its header')

        lo = 0
        lo(1:num_dims) = start_idx(1:num_dims)
        hi = lo + [m, n, p]
        allocate (tc(0:m,0:n,0:p), buf(4, tcond_chunk))
        tc = T_cond_unset
        pre = 0
        do w = 1, int(head(4))
            cnt = tab(1, w)
            if (cnt > 0 .and. all(tab(2:4,w) <= hi) .and. all(tab(5:7,w) >= lo)) then
                do k0 = 0, cnt - 1, tcond_chunk
                    nc = int(min(int(tcond_chunk, 8), cnt - k0))
                    read (file_unit, iostat=ios, pos=hbytes + (pre + k0)*rb + 1) buf(:,1:nc)
                    if (ios /= 0) call s_mpi_abort('Error reading ' // trim(file_loc))
                    do r = 1, nc
                        c = nint(buf(1:3,r)) - lo
                        if (all(c >= 0) .and. all(c <= [m, n, p])) tc(c(1), c(2), c(3)) = buf(4, r)
                    end do
                end do
            end if
            pre = pre + cnt
        end do
        close (file_unit)

        call s_pack_ghost_tcond()
        do i = 1, num_gps
            c = nint(gp_tc_rec(1:3,i)) - lo
            gp_tc_rec(4, i) = tc(c(1), c(2), c(3))
        end do
        call s_unpack_ghost_tcond()

    end subroutine s_read_ib_conduction_temperature

    !> One line at the end of a run if the reacting surface ever failed to solve. Silence means every ghost point on every step
    !! reached the surface-balance tolerance, i.e. the surface chemistry was actually applied everywhere it was asked for.
    impure subroutine s_report_ibm_surface

        integer               :: n_not_converged_glb, n_ill_posed_glb, q, lvl
        integer, dimension(3) :: ft  !< Global fill-tier tallies of one level
        character(len=16)     :: tag

        call s_mpi_reduce_int_sum(n_surface_not_converged, n_not_converged_glb)
        call s_mpi_reduce_int_sum(n_surface_ill_posed, n_ill_posed_glb)
        if (proc_rank == 0) then
            if (n_not_converged_glb > 0) then
                print '(A,I0,A)', ' Immersed-boundary surface chemistry: the Newton solve did not converge at ', &
                    & n_not_converged_glb, &
                    & ' ghost-point updates, which fell back to a chemically inert wall (Twall kept for thermal_bc = 1 or 3, ' &
                    & // 'zero temperature gradient for thermal_bc = 2).'
            end if

            if (n_ill_posed_glb > 0) then
                print '(A,I0,A)', ' Immersed-boundary surface chemistry: ', n_ill_posed_glb, &
                    & ' ghost-point updates had a zero levelset distance, leaving the surface balance undefined; those points ' &
                    & // 'fell back to a chemically inert wall (Twall kept for thermal_bc = 1 or 3, zero temperature gradient for ' // 'thermal_bc = 2). Check the immersed geometry against the grid.'
            end if
        end if

        ! One line per AMR level; coarse ghost points under a fine block are counted on level 0 too
        if (.not. ib_wall_quad) return
        do lvl = 0, amr_max_level
            do q = 1, 3
                call s_mpi_reduce_int_sum(n_fill_tier(q, lvl), ft(q))
            end do
            tag = ''
            if (lvl > 0) write (tag, '(A,I0)') ', AMR level ', lvl
            if (proc_rank == 0 .and. sum(ft) > 0) print '(A,3(F6.1,A))', &
                & ' [WALLFILL] ghost fill tiers' // trim(tag) // ': quadratic ', 1.e2_wp*ft(1)/sum(ft), '%, linear ', &
                & 1.e2_wp*ft(2)/sum(ft), '%, wall value ', 1.e2_wp*ft(3)/sum(ft), '%'
        end do

    end subroutine s_report_ibm_surface

    !> Store ghost point i's surface record in gp_surf. The weight is the surface area the point stands for. Ghost points within two
    !! cell sizes h = dV^(1/d) of the surface fill a band of volume ~2hA, so dV/(2h) each sums to the area A (a length in 2D)
    !! whatever the surface's orientation to the grid. The band lies inside the solid, where a layer at depth s has area A (1 -
    !! s/R)^(d-1) on a curved surface; for circles, spheres and cylinders the weight is scaled back by (R/(R - s))^(d-1), which
    !! removes that bias (sampled on random lattice offsets: < 0.1% mean, 0.3% spread for a sphere at R = 12h). The factor is
    !! limited to R/(R - s) <= 2, so a body of only a few cells cannot give an infinite or negative weight.
    subroutine s_record_gp_surface(i, gp, T_s, mdot_s, reacting)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)           :: i
        type(ghost_point), intent(in) :: gp
        real(wp), intent(in)          :: T_s, mdot_s
        logical, intent(in)           :: reacting
        real(wp)                      :: d, dV, h, R, curv, axial
        real(wp), dimension(3)        :: axis

        d = abs(real(gp%levelset, kind=wp))
        dV = dx(gp%loc(1))*dy(gp%loc(2))
        if (num_dims == 3) dV = dV*dz(gp%loc(3))
        h = dV**(1._wp/real(num_dims, wp))
        if (.not. (d > 0._wp .and. d <= 2._wp*h)) return

        gp_surf(1, i) = dV/(2._wp*h)
        R = patch_ib(gp%ib_patch_id)%radius
        curv = 1._wp
        if (R > 0._wp) curv = R/max(R - d, 0.5_wp*R)
        select case (patch_ib(gp%ib_patch_id)%geometry)
        case (ib_circle, ib_sphere)
            gp_surf(1, i) = gp_surf(1, i)*curv**(num_dims - 1)
        case (ib_cylinder)
            ! A cylinder's axis is its body-frame x, which rotation_matrix maps to its first column in the lab frame. A normal
            ! along the axis is a flat cap, which needs no correction.
            axis(:) = patch_ib(gp%ib_patch_id)%rotation_matrix(:,1)
            axial = abs(dot_product(axis, gp%levelset_norm))
            if (axial < 0.5_wp) gp_surf(1, i) = gp_surf(1, i)*curv
        end select

        ! ib_flux_wall: a reacting body's records come only from the ghosts that set its face fluxes, with their faces' wall area
        if (ib_flux_wall) then
            if (patch_ib(gp%ib_patch_id)%surface_reaction == 1) gp_surf(1, i) = f_flux_wall_record_area(ib_wall_src, &
                & ib_wall_theta, ib_markers, gp%loc)
        end if

        gp_surf(2, i) = T_s
        gp_surf(3, i) = 0._wp
        if (reacting) gp_surf(3, i) = mdot_s

    end subroutine s_record_gp_surface

    !> Write this rank's surface records to D/ib_surface_<rank>_<save>.dat, one line per weighted ghost point.
    impure subroutine s_write_ib_surface(save_count)

        integer, intent(in)                  :: save_count
        character(LEN=path_len + 2*name_len) :: file_loc
        real(wp)                             :: x(3)
        integer                              :: i, g, unit
        integer, allocatable                 :: order(:)

        if (.not. allocated(gp_surf)) return

        if (num_gps > 0) then
            $:GPU_UPDATE(host='[gp_surf(:, 1:num_gps), ghost_points(1:num_gps)]')
        end if

        call s_sort_gps_by_location(order)

        write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/D/ib_surface_', proc_rank, '_', save_count, '.dat'
        open (newunit=unit, file=trim(file_loc), status='replace', action='write')
        do i = 1, num_gps
            g = order(i)
            if (.not. gp_surf(1, g) > 0._wp) cycle
            x = 0._wp
            x(1) = x_cc(ghost_points(g)%loc(1))
            x(2) = y_cc(ghost_points(g)%loc(2))
            if (num_dims == 3) x(3) = z_cc(ghost_points(g)%loc(3))
            write (unit, '(3ES16.8,1X,I0,3ES15.6,3ES16.8)') x, patch_ib(ghost_points(g)%ib_patch_id)%gbl_patch_id, &
                   & ghost_points(g)%levelset_norm, gp_surf(:,g)
        end do
        close (unit)

    end subroutine s_write_ib_surface

    !> The surface temperature is held at Twall: prescribed (thermal_bc = 1) or evolved by the heat balance (thermal_bc = 3).
    pure logical function f_holds_twall(thermal_bc)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: thermal_bc

        f_holds_twall = thermal_bc == 1 .or. thermal_bc == 3

    end function f_holds_twall

    !> The ib_part shard of cell or face (j, k, l); neighbouring cells, which a body's surface visits together, get different ones
    pure integer function f_ib_shard(j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: j, k, l

        f_ib_shard = 1 + iand(j + 5*k + 11*l, ib_nshard - 1)

    end function f_ib_shard

    !> Solid volume of a circle (per unit depth), sphere or cylinder IB.
    pure real(wp) function f_ib_volume(ib_patch)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: ib_patch

        select case (ib_patch%geometry)
        case (ib_circle)
            f_ib_volume = pi*ib_patch%radius**2
        case (ib_sphere)
            f_ib_volume = 4._wp/3._wp*pi*ib_patch%radius**3
        case default
            f_ib_volume = pi*ib_patch%radius**2*max(ib_patch%length_x, ib_patch%length_y, ib_patch%length_z)
        end select

    end function f_ib_volume

    !> Advance the temperature of each thermal_bc = 3 IB over one step, treating it as one lumped body (small Biot number). With
    !! E_solid = m c_s (T - T_ref), the step's face energy in from s_accumulate_ib_face_fluxes and the gasification rate from the
    !! surface solve, m c_s dT/dt = Q_in + c_s (T - T_ref) Mdot_out + heat_power - emissivity sigma A (T^4 - T_rad^4). Mdot_out is
    !! the solid's mass loss: the face mass sum also carries the ghost cells' non-conservation (a net inflow on a reacting
    !! cylinder). T_ref = 298.15 K is the reference of the formation enthalpies, so the gas energy the solid's lost carbon carries
    !! is consistent. The step is explicit: the body's thermal time m c_s/(h A) is many orders of magnitude longer than a flow step.
    !! For D/ib_heat.dat every body also keeps its gasification rate, face mass out and fluid heat in, Q_in (+ the c_s term for
    !! thermal_bc = 3), in patch_ib%mdot, %mdot_face and %qdot.
    impure subroutine s_update_ib_temperatures(dt_step, q_cons_vf)

        real(wp), intent(in) :: dt_step
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), parameter :: sigma_sb = 5.670374419e-8_wp  !< Stefan-Boltzmann constant [W/m^2/K^4]
        real(wp), parameter :: T_ref = 298.15_wp
        real(wp) :: q_in, bank, rel
        integer :: i, g, n_relax

        ! Area and gasification rate from the surface records of the last stage

        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            ib_heat(1, i) = 0._wp
            ib_heat(4, i) = 0._wp
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (num_gps > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, g]')
            do i = 1, num_gps
                if (gp_surf(1, i) > 0._wp) then
                    g = ghost_points(i)%ib_patch_id
                    $:GPU_ATOMIC(atomic='update')
                    ib_heat(1, g) = ib_heat(1, g) + gp_surf(1, i)
                    $:GPU_ATOMIC(atomic='update')
                    ib_heat(4, g) = ib_heat(4, g) + gp_surf(1, i)*gp_surf(3, i)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (ib_mass_relax > 0) call s_ib_shell_mass(q_cons_vf)
        $:GPU_PARALLEL_LOOP(private='[g, i]')
        do g = 1, num_ibs
            ib_heat(5, g) = 0._wp
            do i = 1, ib_nshard
                ib_heat(2, g) = ib_heat(2, g) + ib_part(1, i, g)
                ib_heat(3, g) = ib_heat(3, g) + ib_part(2, i, g)
                ib_heat(5, g) = ib_heat(5, g) + ib_part(3, i, g)
                ib_part(1, i, g) = 0._wp
                ib_part(2, i, g) = 0._wp
                ib_part(3, i, g) = 0._wp
            end do
            ib_heat(6, g) = patch_ib(g)%mass_bank
        end do
        $:END_GPU_PARALLEL_LOOP()
        ! A host round trip is a device sync, so one rank, which holds every body, skips it
        if (num_procs > 1) then
            $:GPU_UPDATE(host='[ib_heat(1:ib_nheat, 1:num_ibs)]')
            call s_ib_exchange_heat(dt_step)
            $:GPU_UPDATE(device='[ib_heat(1:ib_nheat, 1:num_ibs)]')
            if (ib_mass_relax > 0) then
                $:GPU_PARALLEL_LOOP(private='[i]')
                do i = 1, num_ibs
                    patch_ib(i)%mass_bank = ib_heat(6, i)
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        ! On the device, where patch_ib is current (moving IBs update it there). Twall advances only for the bodies this rank holds,
        ! which stays consistent because thermal_bc = 3 requires moving_ibm = 0 (case_validator)
        $:GPU_PARALLEL_LOOP(private='[i, q_in]', copyin='[dt_step]')
        do i = 1, num_ibs
            q_in = ib_heat(2, i)
            if (patch_ib(i)%thermal_bc == 3) q_in = q_in + patch_ib(i)%cp_solid*(patch_ib(i)%Twall - T_ref)*ib_heat(4, i)
            patch_ib(i)%mdot = ib_heat(4, i)
            patch_ib(i)%mdot_face = ib_heat(3, i)
            patch_ib(i)%qdot = q_in
            if (patch_ib(i)%thermal_bc == 3) then
                q_in = q_in + patch_ib(i)%heat_power - patch_ib(i)%emissivity*sigma_sb*ib_heat(1, &
                                       & i)*(patch_ib(i)%Twall**4 - patch_ib(i)%T_rad**4)
                patch_ib(i)%Twall = min(max(patch_ib(i)%Twall &
                         & + dt_step*q_in/(patch_ib(i)%rho_solid*patch_ib(i)%cp_solid*f_ib_volume(patch_ib(i))), T_surface_min), &
                         & T_surface_max)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! ib_mass_relax with one rank (several: s_ib_exchange_heat): the gas should gain exactly the gasified mass; bank the rest
        ! and release 1/ib_mass_relax of the bank
        if (ib_mass_relax > 0 .and. num_procs == 1) then
            n_relax = ib_mass_relax
            $:GPU_PARALLEL_LOOP(private='[i, bank, rel]', copyin='[dt_step, n_relax]')
            do i = 1, num_ibs
                bank = patch_ib(i)%mass_bank
                call s_ib_bank_step(bank, rel, (ib_heat(3, i) - ib_heat(4, i))*dt_step, n_relax)
                patch_ib(i)%mass_bank = bank
                ib_heat(7, i) = rel
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! The next step accumulates afresh
        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            ib_heat(2, i) = 0._wp
            ib_heat(3, i) = 0._wp
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_update_ib_temperatures

    !> One step of an ib_mass_relax bank: add the mass gained beyond gasification, dm, and release 1/ib_mass_relax of the sum.
    pure subroutine s_ib_bank_step(bank, rel, dm, n_relax)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(inout) :: bank
        real(wp), intent(out)   :: rel
        real(wp), intent(in)    :: dm
        integer, intent(in)     :: n_relax  !< ib_mass_relax, passed because a device routine does not see the global

        bank = bank + dm
        rel = bank/real(n_relax, wp)
        bank = bank - rel

    end subroutine s_ib_bank_step

    !> ib_mass_relax: tag the fluid cells within ib_shell_width cells of each body (ib_shell) and add their gas mass to ib_part(3)
    !! (summed into ib_heat(5) by s_update_ib_temperatures). A cell within reach of two bodies goes to the nearer one.
    impure subroutine s_ib_shell_mass(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp) :: rho, vol
        integer :: i, j, k, l, jj, kk, ll, wz, mk, ib_tmp, ib_id, g, xp, yp, zp, sh, c, d2, d2min

        wz = 0
        if (num_dims == 3) wz = ib_shell_width
        ! Shell membership from the markers, which carry the halo, so a cell next to another rank's body is found too
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, jj, kk, ll, mk, ib_tmp, ib_id, g, xp, yp, zp, rho, vol, sh, c, d2, d2min]', &
                            & collapse=3, copyin='[wz]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    g = 0
                    if (ib_markers%sf(j, k, l) == 0) then
                        ! Nearest marker, ties to the smaller patch id (not the encoded marker, whose periodic offset is
                        ! axis-ordered): independent of the scan order, so of orientation
                        mk = 0
                        d2min = huge(d2min)
                        do ll = -wz, wz
                            do kk = -ib_shell_width, ib_shell_width
                                do jj = -ib_shell_width, ib_shell_width
                                    c = ib_markers%sf(j + jj, k + kk, l + ll)
                                    d2 = jj*jj + kk*kk + ll*ll
                                    if (c /= 0 .and. (d2 < d2min .or. (d2 == d2min .and. mod(c - 1, &
                                        & num_gbl_ibs + 1) < mod(mk - 1, num_gbl_ibs + 1)))) then
                                        mk = c
                                        d2min = d2
                                    end if
                                end do
                            end do
                        end do
                        if (mk /= 0) then
                            call s_decode_patch_periodicity(mk, ib_tmp, xp, yp, zp)
                            call s_get_neighborhood_idx(ib_tmp, ib_id)
                            g = max(ib_id, 0)  ! always assigned: a conditional store left g at 0 on NVHPC GPU
                        end if
                    end if
                    ib_shell%sf(j, k, l) = g
                    if (g > 0) then
                        rho = 0._wp
                        do i = eqn_idx%cont%beg, eqn_idx%cont%end
                            rho = rho + q_cons_vf(i)%sf(j, k, l)
                        end do
                        vol = dx(j)*dy(k)
                        if (num_dims == 3) vol = vol*dz(l)
                        sh = f_ib_shard(j, k, l)
                        $:GPU_ATOMIC(atomic='update')
                        ib_part(3, sh, g) = ib_part(3, sh, g) + rho*vol
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_shell_mass

    !> ib_mass_relax: take this step's release of each body's bank (s_update_ib_temperatures) out of the gas within ib_shell_width
    !! cells of the body, scaling the partial densities, momentum, energy and species of those cells by one factor (velocity,
    !! temperature and composition kept). The ghost-cell treatment does not conserve mass at a moving or gasifying wall (face fluxes
    !! into overwritten ghost cells against cells changing side); this returns the gas to exactly the gasified mass, lagged by about
    !! ib_mass_relax steps.
    impure subroutine s_conserve_ib_mass(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp)                                               :: f
        integer                                                :: i, j, k, l, g

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, g, f]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    g = ib_shell%sf(j, k, l)
                    if (g > 0) then
                        if (ib_heat(5, g) > 0._wp) then
                            f = 1._wp - ib_heat(7, g)/ib_heat(5, g)
                            do i = eqn_idx%cont%beg, eqn_idx%E
                                q_cons_vf(i)%sf(j, k, l) = f*q_cons_vf(i)%sf(j, k, l)
                            end do
                            if (chemistry) then
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    q_cons_vf(i)%sf(j, k, l) = f*q_cons_vf(i)%sf(j, k, l)
                                end do
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_conserve_ib_mass

    !> Moving bodies: book the mass and the body-frame energy E' = E - rho u.V_b + rho |V_b|^2/2 of each cell a body covered (+) or
    !! uncovered (-) at this stage's placement, into ib_heat(3) (mass out of the body, as a rate) and ib_heat(2) (energy in). The
    !! grid faces do not move with the wall, so the face fluxes alone miss the gas swept between face and wall. Called at every
    !! placement, after the stage update and before the ghost correction overwrites the covered cells: a cell's content enters the
    !! step's result with weight a_w (the product of the later stages' rk_coef(:,1)/rk_coef(:,4)) and the step-start state with 1 -
    !! a_w, so it is booked as a_w q + (1 - a_w) q_start. V_b = V + omega x r at the cell, r the minimum-image lever arm to the
    !! centroid. A cell passing directly from one body to another is not counted.
    impure subroutine s_add_swept_cells(q_cons_vf, q_start_vf, a_w, dt_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_start_vf
        real(wp), intent(in)                                :: a_w, dt_step
        real(wp), dimension(3)                              :: r, vb, box_len
        real(wp)                                            :: rho, e_b, vol, b_w
        integer                                             :: i, j, k, l, mk, sgn, ib_tmp, ib_id, g, xp, yp, zp

        b_w = 1._wp - a_w
        box_len = 0._wp
        do i = 1, num_dims
            box_len(i) = glb_bounds(i)%end - glb_bounds(i)%beg
        end do
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, mk, sgn, ib_tmp, ib_id, g, xp, yp, zp, r, vb, rho, e_b, vol]', collapse=3, &
                            & copyin='[dt_step, a_w, b_w, box_len]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if ((ib_markers_prev%sf(j, k, l) == 0) .neqv. (ib_markers%sf(j, k, l) == 0)) then
                        sgn = merge(1, -1, ib_markers%sf(j, k, l) /= 0)
                        mk = merge(ib_markers%sf(j, k, l), ib_markers_prev%sf(j, k, l), sgn == 1)
                        call s_decode_patch_periodicity(mk, ib_tmp, xp, yp, zp)
                        call s_get_neighborhood_idx(ib_tmp, ib_id)
                        if (ib_id > 0) then
                            r = 0._wp
                            r(1) = x_cc(j) - patch_ib(ib_id)%x_centroid
                            r(2) = y_cc(k) - patch_ib(ib_id)%y_centroid
                            if (num_dims == 3) r(3) = z_cc(l) - patch_ib(ib_id)%z_centroid
                            do i = 1, num_dims
                                r(i) = r(i) - box_len(i)*anint(r(i)/box_len(i))
                            end do
                            call s_cross_product(patch_ib(ib_id)%angular_vel, r, vb)
                            vb = vb + patch_ib(ib_id)%vel
                            rho = 0._wp
                            do i = eqn_idx%cont%beg, eqn_idx%cont%end
                                rho = rho + a_w*q_cons_vf(i)%sf(j, k, l) + b_w*q_start_vf(i)%sf(j, k, l)
                            end do
                            e_b = a_w*q_cons_vf(eqn_idx%E)%sf(j, k, l) + b_w*q_start_vf(eqn_idx%E)%sf(j, k, &
                                                & l) + 0.5_wp*rho*(vb(1)**2 + vb(2)**2 + vb(3)**2)
                            do i = 1, num_dims
                                e_b = e_b - (a_w*q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                             & l) + b_w*q_start_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l))*vb(i)
                            end do
                            vol = dx(j)*dy(k)
                            if (num_dims == 3) vol = vol*dz(l)
                            g = ib_id
                            $:GPU_ATOMIC(atomic='update')
                            ib_heat(2, g) = ib_heat(2, g) + real(sgn, wp)*e_b*vol/dt_step
                            $:GPU_ATOMIC(atomic='update')
                            ib_heat(3, g) = ib_heat(3, g) - real(sgn, wp)*rho*vol/dt_step
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_advance_markers_prev()

    end subroutine s_add_swept_cells

    !> ib_markers_prev = ib_markers over the whole halo
    impure subroutine s_advance_markers_prev()

        integer :: j, k, l, lo(3), hi(3)

        lo = lbound(ib_markers%sf)
        hi = ubound(ib_markers%sf)
        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[lo, hi]')
        do l = lo(3), hi(3)
            do k = lo(2), hi(2)
                do j = lo(1), hi(1)
                    ib_markers_prev%sf(j, k, l) = ib_markers%sf(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_advance_markers_prev

    !> ib_uncover_fill: a cell the body left at this placement (solid at the last one, fluid now) holds the cold ghost state it last
    !! had, which the stage update then carries as gas. Give it the mean density, species and pressure of the neighbours in its
    !! 3x3(x3) that were fluid at both placements, the state of the flow it moved into; with chemistry the energy is rebuilt from
    !! the mean state's temperature. The velocity stays the ghost's: the mean of the neighbours' makes the drag of a moving body
    !! disagree more with its Galilean twin. A cell with no such neighbour is left as it is. The step-start state q_start, which the
    !! later stages blend in, gets the same state; q_prim_vf (the stage-start primitives) is updated for the image-point
    !! interpolation. The halo is filled as well (see the loop bounds). Call before s_add_swept_cells so the booking sees the filled
    !! state.
    impure subroutine s_fill_uncovered_cells(q_cons_vf, q_start_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_start_vf, q_prim_vf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho_f, alpha_f
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_f, alpha_f
        #:endif
        real(wp), dimension(${NUM_SPECIES}$) :: Ys_f
        real(wp)                             :: s, rho_f, p_f, vsq, ysum, mw_f, T_f, e_f, en_f
        integer                              :: i, j, k, l, jj, kk, ll, cnt, pd, lo(3), hi(3)

        pd = merge(1, 0, p > 0)
        ! The halo is filled too, from the halo's own valid data, so a rank holds what its neighbour filled in the cells that the
        ! image-point interpolation reads: the same cells and arithmetic as on one rank, without a halo exchange.
        lo = lbound(ib_markers%sf)
        hi = ubound(ib_markers%sf)
        do i = 1, num_dims
            lo(i) = lo(i) + 1
            hi(i) = hi(i) - 1
        end do
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, jj, kk, ll, cnt, s, rho_f, p_f, vsq, ysum, mw_f, T_f, e_f, en_f, alpha_rho_f, &
                            & alpha_f, Ys_f]', collapse=3, copyin='[pd, lo, hi]')
        do l = lo(3), hi(3)
            do k = lo(2), hi(2)
                do j = lo(1), hi(1)
                    if (ib_markers_prev%sf(j, k, l) /= 0 .and. ib_markers%sf(j, k, l) == 0) then
                        cnt = 0
                        do ll = -pd, pd
                            do kk = -1, 1
                                do jj = -1, 1
                                    if (ib_markers_prev%sf(j + jj, k + kk, l + ll) == 0 .and. ib_markers%sf(j + jj, k + kk, &
                                        & l + ll) == 0) cnt = cnt + 1
                                end do
                            end do
                        end do
                        if (cnt > 0) then
                            ! mean of every primitive variable but the velocity over the fluid neighbours
                            do i = 1, sys_size
                                if (i >= eqn_idx%mom%beg .and. i <= eqn_idx%mom%end) cycle
                                s = 0._wp
                                do ll = -pd, pd
                                    do kk = -1, 1
                                        do jj = -1, 1
                                            if (ib_markers_prev%sf(j + jj, k + kk, l + ll) == 0 .and. ib_markers%sf(j + jj, &
                                                & k + kk, l + ll) == 0) s = s + q_prim_vf(i)%sf(j + jj, k + kk, l + ll)
                                        end do
                                    end do
                                end do
                                q_prim_vf(i)%sf(j, k, l) = s/real(cnt, wp)
                            end do
                            rho_f = 0._wp
                            do i = 1, num_fluids
                                alpha_rho_f(i) = q_prim_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l)
                                alpha_f(i) = q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l)
                                rho_f = rho_f + alpha_rho_f(i)
                            end do
                            vsq = 0._wp
                            do i = 1, num_dims
                                vsq = vsq + q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)**2
                                q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = rho_f*q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, &
                                          & k, l)
                            end do
                            p_f = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                            do i = 1, num_fluids
                                q_cons_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l) = alpha_rho_f(i)
                                q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = alpha_f(i)
                            end do
                            if (chemistry) then
                                ysum = 0._wp
                                do i = 1, num_species
                                    ysum = ysum + q_prim_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l)
                                end do
                                do i = 1, num_species
                                    Ys_f(i) = q_prim_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l)/ysum
                                    q_prim_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l) = Ys_f(i)
                                    q_cons_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l) = rho_f*Ys_f(i)
                                end do
                                call get_mixture_molecular_weight(Ys_f, mw_f)
                                T_f = p_f*mw_f/(rho_f*gas_constant)
                                call get_mixture_energy_mass(T_f, Ys_f, e_f)
                                en_f = rho_f*e_f + 0.5_wp*rho_f*vsq
                            else
                                call s_compute_energy(p_f, alpha_rho_f, alpha_f, vsq, en_f)
                            end if
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = en_f
                            do i = 1, sys_size
                                q_start_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (.not. ib_exchange) call s_advance_markers_prev()  ! s_add_swept_cells, which also does it, is not called

    end subroutine s_fill_uncovered_cells

    !> Ghost-point indices ordered by grid location, so per-point output is sorted whatever order the list holds
    !! (s_find_ghost_points gathers in location order, but later list edits need not keep it). Heapsort: a particle array carries
    !! O(10^5) ghost points.
    pure subroutine s_sort_gps_by_location(order)

        integer, allocatable, intent(out) :: order(:)
        integer(8), allocatable           :: key(:)
        integer                           :: i, root, g
        integer(8)                        :: k

        allocate (order(max(num_gps, 1)), key(max(num_gps, 1)))
        do i = 1, num_gps
            order(i) = i
            ! Lexicographic (i, j, k) as one key; loc is a buffer-inclusive index, so offset it to stay non-negative.
            k = int(ghost_points(i)%loc(1) + buff_size, 8)
            k = k*int(n + 1 + 2*buff_size, 8) + int(ghost_points(i)%loc(2) + buff_size, 8)
            key(i) = k*int(p + 1 + 2*buff_size, 8) + int(ghost_points(i)%loc(3) + buff_size, 8)
        end do

        do root = num_gps/2, 1, -1
            call s_sift_gp(order, key, root, num_gps)
        end do
        do i = num_gps, 2, -1
            g = order(1); order(1) = order(i); order(i) = g
            call s_sift_gp(order, key, 1, i - 1)
        end do

    end subroutine s_sort_gps_by_location

    !> Restore the max-heap rooted at root over order(1:last), comparing by key.
    pure subroutine s_sift_gp(order, key, root, last)

        integer, intent(inout) :: order(:)
        integer(8), intent(in) :: key(:)
        integer, intent(in)    :: root, last
        integer                :: r, child, g

        r = root
        do
            child = 2*r
            if (child > last) exit
            if (child < last) then
                if (key(order(child + 1)) > key(order(child))) child = child + 1
            end if
            if (key(order(child)) <= key(order(r))) exit
            g = order(r); order(r) = order(child); order(child) = g
            r = child
        end do

    end subroutine s_sift_gp

    !> Compute the image points for each ghost point
    impure subroutine s_compute_image_points()

        real(wp)                        :: dist
        real(wp)                        :: min_cell_width
        real(wp)                        :: temp_loc
        real(wp), pointer, dimension(:) :: s_cc => null()
        integer                         :: bound
        type(ghost_point)               :: gp
        integer                         :: q, dim      !< Iterator variables
        integer                         :: i, j, k, l  !< Location indexes
        integer                         :: patch_id    !< IB Patch ID
        integer                         :: dir
        integer                         :: index
        logical                         :: bounds_error

        bounds_error = .false.

        $:GPU_PARALLEL_LOOP(private='[q, gp, i, j, k, patch_id, dist, min_cell_width, dim, bound, dir, index, temp_loc, s_cc]', &
                            & copy='[bounds_error]', present='[ghost_points]')
        do q = 1, num_gps
            gp = ghost_points(q)
            i = gp%loc(1)
            j = gp%loc(2)
            k = gp%loc(3)

            ! Calculate and store the precise location of the image point
            patch_id = gp%ib_patch_id
            dist = abs(real(gp%levelset, kind=wp))
            ghost_points(q)%ip_dist = dist
            ! keep the image point stencil off the wall
            if (ib_second_order_vel) then
                min_cell_width = min(dx(i), dy(j))
                if (p > 0) min_cell_width = min(min_cell_width, dz(k))
                ghost_points(q)%ip_dist = max(dist, ib_ip_min_dist*min_cell_width)
            end if
            ! Per-component scalar stores: the array-valued form of this assignment was intermittently overwritten by another
            ! thread on AMD GPUs (a few ip_loc(2) of ~4000 per call), which sent the image-point search past the halo
            ghost_points(q)%ip_loc(1) = x_cc(i) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(1)
            ghost_points(q)%ip_loc(2) = y_cc(j) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(2)
            ghost_points(q)%ip_loc(3) = 0._wp
            if (p > 0) ghost_points(q)%ip_loc(3) = z_cc(k) + (dist + ghost_points(q)%ip_dist)*gp%levelset_norm(3)

            ! Find the closest grid point to the image point
            do dim = 1, num_dims
                ! s_cc points to the dim array we need
                if (dim == 1) then
                    s_cc => x_cc
                    bound = m + buff_size - 1
                else if (dim == 2) then
                    s_cc => y_cc
                    bound = n + buff_size - 1
                else
                    s_cc => z_cc
                    bound = p + buff_size - 1
                end if

                if (f_approx_equal(gp%levelset_norm(dim), 0._wp)) then
                    ! if the ghost point is almost equal to a cell location, we set it equal and continue
                    ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim)
                else
                    if (gp%levelset_norm(dim) > 0) then
                        dir = 1
                    else
                        dir = -1
                    end if

                    index = ghost_points(q)%loc(dim)
                    temp_loc = ghost_points(q)%ip_loc(dim)
                    do while ((temp_loc < s_cc(index) .or. temp_loc > s_cc(index + 1)) .and. (.not. bounds_error))
                        index = index + dir
                        if (index < -buff_size .or. index > bound) then
#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP)
                            print *, "A required image point is not located in this computational domain."
                            print *, "Ghost Point is located at :"
                            if (p == 0) then
                                print *, [x_cc(i), y_cc(j)]
                            else
                                print *, [x_cc(i), y_cc(j), z_cc(k)]
                            end if
                            print *, "We are searching in dimension ", dim, " for image point at ", ghost_points(q)%ip_loc(:)
                            print *, "Domain size: "
                            print *, "x: ", x_cc(-buff_size), " to: ", x_cc(m + buff_size - 1)
                            print *, "y: ", y_cc(-buff_size), " to: ", y_cc(n + buff_size - 1)
                            if (p /= 0) print *, "z: ", z_cc(-buff_size), " to: ", z_cc(p + buff_size - 1)
                            print *, "Image point is located approximately ", &
                                & (ghost_points(q)%loc(dim) - ghost_points(q) %ip_loc(dim))/(s_cc(1) - s_cc(0)), " grid cells away"
                            print *, "Levelset ", dist, " and Norm: ", gp%levelset_norm
                            print *, &
                                & "A short term fix may include increasing buff_size further in m_helper_basic (currently set to a minimum of 10)"
#endif
                            bounds_error = .true.
                        end if
                    end do

                    ghost_points(q)%ip_grid(dim) = index
                    if (ghost_points(q)%DB(dim) == -1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) + 1
                    else if (ghost_points(q)%DB(dim) == 1) then
                        ghost_points(q)%ip_grid(dim) = ghost_points(q)%loc(dim) - 1
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        @:PROHIBIT(bounds_error, "Ghost Point and Image Point on Different Processors. Exiting")

    end subroutine s_compute_image_points

    !> Count the number of ghost points for memory allocation
    subroutine s_find_num_ghost_points(num_gps_out)

        integer, intent(out) :: num_gps_out
        integer              :: i, j, k, ii, jj, kk, gp_layers_z  !< Iterator variables
        integer              :: num_gps_local                     !< local copies of the gp count to support GPU compute
        logical              :: is_gp

        num_gps_local = 0
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp]', copy='[num_gps_local]', firstprivate='[gp_layers_z]', &
                            & copyin='[ib_markers%sf]', collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search

                        if (is_gp) then
                            $:GPU_ATOMIC(atomic='update')
                            num_gps_local = num_gps_local + 1
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        num_gps_out = num_gps_local

    end subroutine s_find_num_ghost_points

    !> Find the ghost points (solid cells within gp_layers of a fluid cell) in lexicographic (i,j,k) order on every backend, so the
    !! list and anything reduced over it do not depend on the thread schedule (#2039): count them per (i,j) column, scan the counts
    !! into column offsets (per i-row in parallel, then the row totals), then fill each column in k order.
    subroutine s_find_ghost_points()

        integer :: i, j, k, ii, jj, kk, gp_layers_z, slab_n, run, c
        integer :: patch_id, encoded_patch_id, neighborhood_patch_id, xp, yp, zp
        integer :: min_id  !< smallest neighborhood id seen (see s_check_gp_patch_ids)
        logical :: is_gp

        min_id = huge(1)
        gp_layers_z = gp_layers
        if (p == 0) gp_layers_z = 0
        slab_n = m

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp, c]', copyin='[ib_markers%sf]', firstprivate='[gp_layers_z]', &
                            & present='[gp_col]', collapse=2)
        do i = 0, m
            do j = 0, n
                c = 0
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        is_gp = .true.
                                        exit marker_search
                                    end if
                                end do
                            end do
                        end do marker_search
                        if (is_gp) c = c + 1
                    end if
                end do
                gp_col(i, j) = c
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! exclusive scan of each row's column counts, the row total into gp_slab; then of the row totals (and their sum), one thread
        $:GPU_PARALLEL_LOOP(private='[i, j, run, c]', present='[gp_col, gp_slab]')
        do i = 0, m
            run = 0
            do j = 0, n
                c = gp_col(i, j)
                gp_col(i, j) = run
                run = run + c
            end do
            gp_slab(i) = run
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(private='[i, run, c]', firstprivate='[slab_n]', present='[gp_slab]')
        do k = 1, 1
            run = 0
            do i = 0, slab_n
                c = gp_slab(i)
                gp_slab(i) = run
                run = run + c
            end do
            gp_slab(slab_n + 1) = run
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[i, j, k, ii, jj, kk, is_gp, c, patch_id, encoded_patch_id, neighborhood_patch_id, xp, yp, &
                            & zp]', copyin='[ib_markers%sf, glb_bounds]', firstprivate='[gp_layers_z]', reduction='[[min_id]]', &
                            & reductionOp='[min]', present='[ghost_points, gp_slab, gp_col]', collapse=2)
        do i = 0, m
            do j = 0, n
                c = gp_slab(i) + gp_col(i, j)
                ! skip a column the count pass found empty: its end is the next column's (or row's) start
                if (j < n) then
                    if (gp_slab(i) + gp_col(i, j + 1) == c) cycle
                else
                    if (gp_slab(i + 1) == c) cycle
                end if
                do k = 0, p
                    if (ib_markers%sf(i, j, k) /= 0) then
                        is_gp = .false.
                        marker_search2: do ii = i - gp_layers, i + gp_layers
                            do jj = j - gp_layers, j + gp_layers
                                do kk = k - gp_layers_z, k + gp_layers_z
                                    if (ib_markers%sf(ii, jj, kk) == 0) then
                                        is_gp = .true.
                                        exit marker_search2
                                    end if
                                end do
                            end do
                        end do marker_search2
                        if (is_gp) then
                            c = c + 1
                            ghost_points(c)%loc = [i, j, k]
                            ghost_points(c)%T_cond = T_cond_unset
                            encoded_patch_id = ib_markers%sf(i, j, k)
                            call s_decode_patch_periodicity(encoded_patch_id, patch_id, xp, yp, zp)
                            call s_get_neighborhood_idx(patch_id, neighborhood_patch_id)
                            min_id = min(min_id, neighborhood_patch_id)
                            neighborhood_patch_id = max(neighborhood_patch_id, 1)
                            ghost_points(c)%ib_patch_id = neighborhood_patch_id
                            ghost_points(c)%x_periodicity = xp
                            ghost_points(c)%y_periodicity = yp
                            ghost_points(c)%z_periodicity = zp
                            ghost_points(c)%slip = patch_ib(neighborhood_patch_id)%slip

                            #:for X, ID, IDX in [('x', 1, 'i'), ('y', 2, 'j'), ('z', 3, 'k')]
                                ghost_points(c)%DB(${ID}$) = 0
                                if (${ID}$ <= num_dims) then  ! Onlyrun the numeher of dimensions present
                                    if (ib_bc_${X}$%beg /= BC_PERIODIC) then
                                        if ((${X}$_cc(${IDX}$) - d${X}$(${IDX}$)) < glb_bounds(${ID}$)%beg) then
                                            ! if the grid cell is in a wall on the "left"
                                            ghost_points(c)%DB(${ID}$) = -1
                                        else if ((${X}$_cc(${IDX}$) + d${X}$(${IDX}$)) > glb_bounds(${ID}$)%end) then
                                            ! if the grid cell is in a wall on the "right"
                                            ghost_points(c)%DB(${ID}$) = 1
                                        end if
                                    end if
                                end if
                            #:endfor
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_check_gp_patch_ids(min_id)

    end subroutine s_find_ghost_points

    !> Abort if a ghost point's body was not in the rank's neighborhood list (s_get_neighborhood_idx returned -1). The setters clamp
    !! the id to 1 so no device read of patch_ib(id <= 0) happens before this check.
    subroutine s_check_gp_patch_ids(min_id)

        integer, intent(in) :: min_id

        if (min_id <= 0) call s_mpi_abort('IBM: a ghost point marks a body that is not in this rank''s IB neighborhood list')

    end subroutine s_check_gp_patch_ids

    !> Compute the interpolation coefficients for image points
    subroutine s_compute_interpolation_coeffs()

        real(wp), dimension(2, 2, 2) :: dist
        real(wp), dimension(2, 2, 2) :: fluid_mask
        real(wp), dimension(2, 2, 2) :: interp_coeffs
        real(wp), dimension(3, 3, 3) :: stencil_coeffs          !< Coefficients on the stencil that starts at ip_stencil
        real(wp), dimension(3, 3)    :: wd                      !< 1D TSC weights per direction
        real(wp)                     :: fluid_w                 !< Fraction of the TSC weight carried by fluid cells
        integer, dimension(3)        :: ic                      !< Cell centring the TSC stencil
        type(ghost_point)            :: gp
        integer                      :: q, i, j, k, ii, jj, kk  !< Grid indexes and iterators
        integer                      :: n_buried                !< Ghost points of this rank without fluid to mirror
        integer(kind=8)              :: n_buried_loc, n_buried_glb

        $:GPU_PARALLEL_LOOP(private='[q, i, j, k, ii, jj, kk, dist, gp, interp_coeffs, fluid_mask, stencil_coeffs, wd, fluid_w, &
                            & ic]', copyin='[ib_markers%sf]', present='[ghost_points]')
        do q = 1, num_gps
            gp = ghost_points(q)
            stencil_coeffs = 0._wp
            ic = 0
            if (ib_ip_interp == 2) then
                call s_tsc_weights(ib_markers, gp%ip_loc, ic, wd, fluid_w)
                stencil_coeffs = f_tsc_stencil(wd)
                ! An injecting body's ghost cells hold the imposed injectant state, not a mirror of the gas, so they must not feed
                ! its own image point: keep the fluid cells only, as the inverse-distance stencil does
                if (patch_ib(gp%ib_patch_id)%inj_species > 0) then
                    fluid_w = 0._wp
                    do kk = 1, merge(3, 1, p /= 0)
                        do jj = 1, 3
                            do ii = 1, 3
                                if (ib_markers%sf(ic(1) + ii - 2, ic(2) + jj - 2, ic(3) + (kk - 2)*min(p, &
                                    & 1)) /= 0) stencil_coeffs(ii, jj, kk) = 0._wp
                                fluid_w = fluid_w + stencil_coeffs(ii, jj, kk)
                            end do
                        end do
                    end do
                    if (fluid_w > 0._wp) stencil_coeffs = stencil_coeffs/fluid_w
                end if
                ghost_points(q)%interp_valid = fluid_w > 0._wp
            else
                ! Get the interpolation points
                i = gp%ip_grid(1)
                j = gp%ip_grid(2)
                if (p /= 0) then
                    k = gp%ip_grid(3)
                else
                    k = 0
                end if

                ! get the distance to a cell in each direction
                dist = 0._wp
                do ii = 0, 1
                    do jj = 0, 1
                        if (p == 0) then
                            dist(1 + ii, 1 + jj, 1) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2)
                        else
                            do kk = 0, 1
                                dist(1 + ii, 1 + jj, &
                                     & 1 + kk) = sqrt((x_cc(i + ii) - gp%ip_loc(1))**2 + (y_cc(j + jj) - gp%ip_loc(2))**2 &
                                     & + (z_cc(k + kk) - gp%ip_loc(3))**2)
                            end do
                        end if
                    end do
                end do

                ! flag the fluid cells of the stencil
                fluid_mask = 0._wp
                do kk = 0, min(p, 1)
                    do jj = 0, 1
                        do ii = 0, 1
                            if (ib_markers%sf(i + ii, j + jj, k + kk) == 0) fluid_mask(ii + 1, jj + 1, kk + 1) = 1._wp
                        end do
                    end do
                end do

                call s_interpolation_weights(gp, i, j, k, dist, fluid_mask, interp_coeffs)

                stencil_coeffs(1:2,1:2,1:2) = interp_coeffs
                ghost_points(q)%interp_valid = any(interp_coeffs > 0._wp)
            end if

            ! An image point buried in a neighboring IB gets no weights at all: its ghost point is
            ! averaged from its own neighbors instead, in s_ibm_correct_state.
            do kk = 1, 3
                do jj = 1, 3
                    do ii = 1, 3
                        ghost_points(q)%interp_coeffs(ii, jj, kk) = stencil_coeffs(ii, jj, kk)
                    end do
                end do
                ! The stencil's lower corner: one cell below the TSC centre, or the image point's grid cell
                ghost_points(q)%ip_stencil(kk) = merge(ic(kk) - 1, gp%ip_grid(kk), ib_ip_interp == 2)
            end do
            if (p == 0) ghost_points(q)%ip_stencil(3) = 0
            corrected_gps%sf(gp%loc(1), gp%loc(2), gp%loc(3)) = merge(1, 2, ghost_points(q)%interp_valid)
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Collective: whether any rank holds a ghost point without fluid to mirror (the AMR frames and non-polytropic QBMM
        ! keep the rank-local average)
        n_buried = 0
        $:GPU_PARALLEL_LOOP(private='[q]', reduction='[[n_buried]]', reductionOp='[+]', present='[ghost_points]')
        do q = 1, num_gps
            if (.not. ghost_points(q)%interp_valid) n_buried = n_buried + 1
        end do
        $:END_GPU_PARALLEL_LOOP()
        ! not with AMR: fine frames call this per block, on some ranks only, so it must not reach a collective
        ib_any_buried = .false.
        if (.not. amr .and. .not. (qbmm .and. .not. polytropic)) then
            n_buried_loc = int(n_buried, kind=8)
            call s_mpi_allreduce_integer_sum(n_buried_loc, n_buried_glb)
            ib_any_buried = n_buried_glb > 0
        end if

    end subroutine s_compute_interpolation_coeffs

    !> Fill the halo cells of the corrected ghosts: the state (both forms) and the flags of corrected_gps, from the neighbouring
    !! rank at a processor edge or from the periodic image, so the buried-point average sees the same ghosts as a single rank
    impure subroutine s_ibm_exchange_gp_halos(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf
        integer                                                :: j, k, l
        integer                                                :: lo1, hi1, lo2, hi2, lo3, hi3

        lo1 = lbound(corrected_gps%sf, 1); hi1 = ubound(corrected_gps%sf, 1)
        lo2 = lbound(corrected_gps%sf, 2); hi2 = ubound(corrected_gps%sf, 2)
        lo3 = lbound(corrected_gps%sf, 3); hi3 = ubound(corrected_gps%sf, 3)

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[lo1, hi1, lo2, hi2, lo3, hi3]')
        do l = lo3, hi3
            do k = lo2, hi2
                do j = lo1, hi1
                    corrected_gps_r(1)%sf(j, k, l) = real(corrected_gps%sf(j, k, l), kind=stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_ibm_fill_halos(q_cons_vf, sys_size)
        call s_ibm_fill_halos(q_prim_vf, sys_size)
        call s_ibm_fill_halos(corrected_gps_r, 1)

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[lo1, hi1, lo2, hi2, lo3, hi3]')
        do l = lo3, hi3
            do k = lo2, hi2
                do j = lo1, hi1
                    corrected_gps%sf(j, k, l) = nint(corrected_gps_r(1)%sf(j, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_exchange_gp_halos

    !> Halo cells of q(1:nv) along each direction, from the neighbouring rank at a processor edge (the exchange of the RHS) and from
    !! the periodic image where the direction has one rank; other boundaries keep their halos
    impure subroutine s_ibm_fill_halos(q, nv)

        type(scalar_field), dimension(1:), intent(inout) :: q
        integer, intent(in)                              :: nv
        integer                                          :: edge, loc

        #:for X, D in [('x', 1), ('y', 2), ('z', 3)]
            if (num_dims >= ${D}$) then
                do loc = -1, 1, 2
                    edge = merge(bc_${X}$%beg, bc_${X}$%end, loc == -1)
                    if (edge >= 0) then
                        call s_mpi_sendrecv_variables_buffers(q, ${D}$, loc, nv)
                    else if (edge == BC_PERIODIC) then
                        call s_ibm_wrap_halo(q, nv, ${D}$, loc)
                    end if
                end do
            end if
        #:endfor

    end subroutine s_ibm_fill_halos

    !> Copy the opposite interior cells into the halo of direction d at side loc (-1 beg, 1 end). Earlier directions are already
    !! filled, so the corner cells come out right.
    impure subroutine s_ibm_wrap_halo(q, nv, d, loc)

        type(scalar_field), dimension(1:), intent(inout) :: q
        integer, intent(in)                              :: nv, d, loc
        integer                                          :: i, j, k, l
        integer                                          :: jb, je, kb, ke, lb, le, shj, shk, shl, ext

        jb = 0; je = m; kb = 0; ke = n; lb = 0; le = p
        if (d > 1) then
            jb = -buff_size; je = m + buff_size
        end if
        if (d > 2) then
            kb = -buff_size; ke = n + buff_size
        end if
        shj = 0; shk = 0; shl = 0
        ext = merge(1, -1, loc == -1)*(merge(m, merge(n, p, d == 2), d == 1) + 1)
        if (d == 1) then
            shj = ext
            jb = merge(-buff_size, m + 1, loc == -1); je = merge(-1, m + buff_size, loc == -1)
        else if (d == 2) then
            shk = ext
            kb = merge(-buff_size, n + 1, loc == -1); ke = merge(-1, n + buff_size, loc == -1)
        else
            shl = ext
            lb = merge(-buff_size, p + 1, loc == -1); le = merge(-1, p + buff_size, loc == -1)
        end if

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4, copyin='[jb, je, kb, ke, lb, le, shj, shk, shl]')
        do l = lb, le
            do k = kb, ke
                do j = jb, je
                    do i = 1, nv
                        q(i)%sf(j, k, l) = q(i)%sf(j + shj, k + shk, l + shl)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_wrap_halo

    !> Triangular-shaped-cloud (quadratic B-spline) weights of a point on the 3x3(x3) stencil centred on the cell nearest to it
    !! (Ezra, Peles & Kozak 2025, Int. J. Heat Mass Transf. 252:127314), and the fraction of that weight carried by fluid cells. The
    !! kernel is a partition of unity and smooth across cell boundaries, so an interpolated wall quantity converges one order faster
    !! than from inverse-distance or bilinear weights, which are only continuous. The weights use the fractional position in the
    !! bracketing cell at integer offsets, which keeps the partition of unity exactly on a stretched grid. The stencil is kept
    !! inside the halo. Solid cells in it carry the ghost values of the last correction (see ib_keep_ghost); a point buried in a
    !! body has no fluid weight.
    subroutine s_tsc_weights(markers, loc, stencil_center, axis_weights, fluid_frac)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(integer_field), intent(in) :: markers  !< ib_markers, passed so the device routine reads the caller's present copy
        real(wp), dimension(3), intent(in) :: loc
        integer, dimension(3), intent(out) :: stencil_center  !< Cell nearest to loc, centring the stencil
        real(wp), dimension(3, 3), intent(out) :: axis_weights  !< 1D weights on stencil_center - 1, +0, +1, per direction
        real(wp), intent(out) :: fluid_frac  !< Fraction of the stencil weight on fluid cells
        integer, dimension(3) :: containing_cell  !< Cell whose boundaries hold loc
        real(wp), dimension(3) :: cell_coord  !< Computational coordinates from s_locate_cell (unused)
        real(wp) :: cell_offset  !< Offset of loc from the stencil's centre cell, in center-to-center spacings
        real(wp) :: weight_sum, fluid_weight, cell_weight
        integer :: ii, jj, kk, n_solid

        stencil_center = 0
        axis_weights = 0._wp
        ! Binary search, so the lookup holds on any stretched grid
        call s_locate_cell(loc, containing_cell, cell_coord)
        #:for D, CC, N in [(1, 'x_cc', 'm'), (2, 'y_cc', 'n'), (3, 'z_cc', 'p')]
            if (${D}$ < 3 .or. num_dims == 3) then
                ! Lower of the two cell centres bracketing loc, kept inside the halo
                stencil_center(${D}$) = containing_cell(${D}$)
                if (loc(${D}$) < ${CC}$(stencil_center(${D}$))) stencil_center(${D}$) = stencil_center(${D}$) - 1
                stencil_center(${D}$) = max(min(stencil_center(${D}$), ${N}$ + buff_size - 2), -buff_size + 2)
                cell_offset = (loc(${D}$) - ${CC}$(stencil_center(${D}$)))/(${CC}$(stencil_center(${D}$) + 1) &
                               & - ${CC}$(stencil_center(${D}$)))
                ! Centre the stencil on the nearer of the two
                if (cell_offset > 0.5_wp) then
                    stencil_center(${D}$) = stencil_center(${D}$) + 1
                    cell_offset = cell_offset - 1._wp
                end if
                axis_weights(1, ${D}$) = 0.5_wp*(1.5_wp - abs(cell_offset + 1._wp))**2
                axis_weights(2, ${D}$) = 0.75_wp - cell_offset*cell_offset
                axis_weights(3, ${D}$) = 0.5_wp*(1.5_wp - abs(cell_offset - 1._wp))**2
            end if
        #:endfor
        if (p == 0) axis_weights(1, 3) = 1._wp  ! 2D: the consumers index this slot as k - k1 + 1 = 1

        weight_sum = 0._wp; fluid_weight = 0._wp; n_solid = 0
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, merge(3, 1, p /= 0)
                    cell_weight = axis_weights(ii, 1)*axis_weights(jj, 2)*axis_weights(kk, 3)
                    weight_sum = weight_sum + cell_weight
                    if (markers%sf(stencil_center(1) + ii - 2, stencil_center(2) + jj - 2, stencil_center(3) + (kk - 2)*min(p, &
                        & 1)) == 0) then
                        fluid_weight = fluid_weight + cell_weight
                    else if (cell_weight > 0._wp) then
                        n_solid = n_solid + 1  ! a solid cell in a zero-weight plane (offset exactly 1/2) takes no weight
                    end if
                end do
            end do
        end do
        ! Exactly 1 for an all-fluid stencil: a compiler that contracts only one of the two sums into an FMA leaves the ratio an ulp
        ! below 1, and s_lag3_weights would then switch to TSC on round-off (#2443)
        fluid_frac = 1._wp
        if (n_solid > 0) fluid_frac = fluid_weight/weight_sum

    end subroutine s_tsc_weights

    !> Apply TSC weights from s_tsc_weights to one field.
    function f_tsc_apply(sf, ic, wd) result(val)

        $:GPU_ROUTINE(function_name='f_tsc_apply', parallelism='[seq]', cray_inline=True)
        type(scalar_field), intent(in)        :: sf
        integer, dimension(3), intent(in)     :: ic
        real(wp), dimension(3, 3), intent(in) :: wd
        real(wp)                              :: val
        integer                               :: ii, jj, kk

        val = 0._wp
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, merge(3, 1, p /= 0)
                    val = val + wd(ii, 1)*wd(jj, 2)*wd(kk, 3)*sf%sf(ic(1) + ii - 2, ic(2) + jj - 2, ic(3) + (kk - 2)*min(p, 1))
                end do
            end do
        end do

    end function f_tsc_apply

    !> TSC weights from s_tsc_weights applied to the gas cells only, renormalized
    function f_tsc_fluid_mean(sf, markers, ic, wd) result(val)

        $:GPU_ROUTINE(function_name='f_tsc_fluid_mean', parallelism='[seq]', cray_inline=True)
        type(scalar_field), intent(in)        :: sf
        type(integer_field), intent(in)       :: markers
        integer, dimension(3), intent(in)     :: ic
        real(wp), dimension(3, 3), intent(in) :: wd
        real(wp)                              :: val, w, wsum
        integer                               :: ii, jj, kk, kz

        val = 0._wp
        wsum = 0._wp
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, merge(3, 1, p /= 0)
                    kz = ic(3) + (kk - 2)*min(p, 1)
                    if (markers%sf(ic(1) + ii - 2, ic(2) + jj - 2, kz) /= 0) cycle
                    w = wd(ii, 1)*wd(jj, 2)*wd(kk, 3)
                    val = val + w*sf%sf(ic(1) + ii - 2, ic(2) + jj - 2, kz)
                    wsum = wsum + w
                end do
            end do
        end do
        val = val/max(wsum, tiny(1._wp))

    end function f_tsc_fluid_mean

    !> Lagrange (biquadratic, triquadratic in 3D) weights of point loc on the stencil s_tsc_weights centres there, from the cells'
    !! own centres: an interpolant, where TSC smooths a sample by about (h^2/4) lap phi. Where a stencil cell is solid the TSC
    !! weights are kept.
    subroutine s_lag3_weights(markers, loc, ic, wd, fluid_frac)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(integer_field), intent(in)        :: markers
        real(wp), dimension(3), intent(in)     :: loc
        integer, dimension(3), intent(out)     :: ic
        real(wp), dimension(3, 3), intent(out) :: wd
        real(wp), intent(out)                  :: fluid_frac
        real(wp)                               :: xm, x0, xp
        integer                                :: dim

        call s_tsc_weights(markers, loc, ic, wd, fluid_frac)
        if (fluid_frac < 1._wp) return
        do dim = 1, num_dims
            if (dim == 1) then
                xm = x_cc(ic(1) - 1); x0 = x_cc(ic(1)); xp = x_cc(ic(1) + 1)
            else if (dim == 2) then
                xm = y_cc(ic(2) - 1); x0 = y_cc(ic(2)); xp = y_cc(ic(2) + 1)
            else
                xm = z_cc(ic(3) - 1); x0 = z_cc(ic(3)); xp = z_cc(ic(3) + 1)
            end if
            wd(1, dim) = (loc(dim) - x0)*(loc(dim) - xp)/((xm - x0)*(xm - xp))
            wd(2, dim) = (loc(dim) - xm)*(loc(dim) - xp)/((x0 - xm)*(x0 - xp))
            wd(3, dim) = (loc(dim) - xm)*(loc(dim) - x0)/((xp - xm)*(xp - x0))
        end do

    end subroutine s_lag3_weights

    !> Ghost value from the quadratic wall-normal profile, bounded by the linear mirror (minmod). The two extrapolations agree
    !! exactly for a linear profile, so this limits only the curvature correction, and to the size the linear mirror already
    !! produces. Unlimited, the quadratic amplifies that correction with ghost depth and can leave unphysical near-wall states.
    function f_quad_fill(phi_w, phi_1, phi_2, lw, llin) result(phi_g)

        $:GPU_ROUTINE(function_name='f_quad_fill', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in)               :: phi_w, phi_1, phi_2
        real(wp), dimension(3), intent(in) :: lw
        real(wp), dimension(2), intent(in) :: llin
        real(wp)                           :: phi_g, d_quad, d_lin

        d_quad = lw(1)*phi_w + lw(2)*phi_1 + lw(3)*phi_2 - phi_w
        d_lin = llin(1)*phi_w + llin(2)*phi_1 - phi_w
        if (d_quad*d_lin <= 0._wp) then
            phi_g = phi_w + d_lin  ! opposite signs: the curvature term has reversed the extrapolation
        else if (abs(d_quad) <= abs(d_lin)) then
            phi_g = phi_w + d_quad
        else
            phi_g = phi_w + d_lin
        end if

    end function f_quad_fill

    !> Probes and Lagrange weights of the quadratic wall-normal fill of a ghost point. The wall is at s = 0, the probes at s1 and s2
    !! = 2 s1 along the surface normal from the wall foot, the ghost at s = -d (d = |levelset|):
    !!
    !!     phi_g = L0 phi_w + L1 phi_1 + L2 phi_2
    !!     L0 = (d + s1)(d + s2)/(s1 s2),  L1 = -d (d + s2)/(s1 (s2 - s1)),  L2 = d (d + s1)/(s2 (s2 - s1))
    !!
    !! (d = ghost_depth, s1 = probe1_dist, s2 = probe2_dist below.) L0 + L1 + L2 = 1 and at d = 0 they reduce to (1, 0, 0), so the
    !! wall value is exact on the wall. s1 = max(h, d) bounds the weights by (3, -3, 1) at any ghost depth. Tier 1 uses both probes;
    !! tier 2 is the linear mirror where the outer probe has too little fluid around it (a narrow gap, a rim, a neighbouring body)
    !! or lies beyond the halo (a ghost several cells deep); tier 3 leaves the wall value (no quadratic fill).
    subroutine s_wall_fill(gp, wf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(ghost_point), intent(in) :: gp
        type(wall_fill), intent(out)  :: wf
        real(wp), parameter           :: FLUID_FRAC_MIN = 0.5_wp  !< A probe needs this much fluid weight to be used
        real(wp), dimension(3)        :: wall_normal, wall_loc  !< Unit outward normal, and the surface point nearest the ghost
        real(wp)                      :: ghost_depth  !< Distance from the ghost cell centre to the surface
        real(wp)                      :: probe1_dist, probe2_dist  !< Distances of the probes from the wall along wall_normal
        real(wp)                      :: probe1_fluid_frac, probe2_fluid_frac  !< Fraction of each probe's TSC weight on fluid cells
        real(wp)                      :: cell_vol
        real(wp), dimension(3)        :: probe_loc  !< Named, not an expression argument (see the Silent-Failure Traps)
        real(wp)                      :: reach  !< Largest probe offset from the ghost, in its cell widths, inside the halo
        real(wp), dimension(3)        :: width  !< Cell widths at the ghost
        logical                       :: in1, in2  !< Each probe's stencil lies inside the halo
        integer                       :: d

        cell_vol = dx(gp%loc(1))*dy(gp%loc(2))
        if (p > 0) cell_vol = cell_vol*dz(gp%loc(3))
        ghost_depth = abs(real(gp%levelset, kind=wp))
        probe1_dist = max(cell_vol**(1._wp/real(num_dims, wp)), ghost_depth)
        probe2_dist = 2._wp*probe1_dist

        wall_normal = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        wall_loc = 0._wp
        wall_loc(1) = x_cc(gp%loc(1)); wall_loc(2) = y_cc(gp%loc(2))
        if (p > 0) wall_loc(3) = z_cc(gp%loc(3))
        wall_loc = wall_loc + ghost_depth*wall_normal
        probe_loc = wall_loc + probe1_dist*wall_normal
        call s_tsc_weights(ib_markers, probe_loc, wf%ic1, wf%wd1, probe1_fluid_frac)
        probe_loc = wall_loc + probe2_dist*wall_normal
        call s_tsc_weights(ib_markers, probe_loc, wf%ic2, wf%wd2, probe2_fluid_frac)

        ! A deep ghost's probe can lie beyond the halo, where its stencil is clamped and the fill would depend on the domain
        ! decomposition. Decided from the ghost's own geometry, so every decomposition drops the same ghosts a tier.
        reach = real(buff_size, wp) - 2.5_wp
        width(1) = dx(gp%loc(1)); width(2) = dy(gp%loc(2)); width(3) = 1._wp
        if (p > 0) width(3) = dz(gp%loc(3))
        in1 = .true.; in2 = .true.
        do d = 1, num_dims
            in1 = in1 .and. abs((ghost_depth + probe1_dist)*wall_normal(d)) <= reach*width(d)
            in2 = in2 .and. abs((ghost_depth + probe2_dist)*wall_normal(d)) <= reach*width(d)
        end do

        wf%llin = [1._wp + ghost_depth/probe1_dist, -ghost_depth/probe1_dist]
        wf%lw = [wf%llin(1), wf%llin(2), 0._wp]
        wf%tier = 3
        if (probe1_fluid_frac < FLUID_FRAC_MIN .or. .not. in1) return
        wf%tier = 2
        if (probe2_fluid_frac < FLUID_FRAC_MIN .or. .not. in2) return
        wf%lw = [(ghost_depth + probe1_dist)*(ghost_depth + probe2_dist)/(probe1_dist*probe2_dist), &
                  & -ghost_depth*(ghost_depth + probe2_dist)/(probe1_dist*(probe2_dist - probe1_dist)), &
                  & ghost_depth*(ghost_depth + probe1_dist)/(probe2_dist*(probe2_dist - probe1_dist))]
        wf%tier = 1

    end subroutine s_wall_fill

    !> Tensor product of the 1D TSC weights: the coefficient on each cell of the 3x3x3 stencil.
    function f_tsc_stencil(wd) result(c)

        $:GPU_ROUTINE(function_name='f_tsc_stencil', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(3, 3), intent(in) :: wd
        real(wp), dimension(3, 3, 3)          :: c
        integer                               :: ii, jj, kk

        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, 3
                    c(ii, jj, kk) = wd(ii, 1)*wd(jj, 2)*wd(kk, 3)
                end do
            end do
        end do

    end function f_tsc_stencil

    !> Fluid-masked image-point weights: bilinear for second-order IBM, else inverse distance squared
    subroutine s_interpolation_weights(gp, i, j, k, dist, fluid_mask, interp_coeffs)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ghost_point), intent(in)             :: gp
        integer, intent(in)                       :: i, j, k      !< lower corner of the stencil
        real(wp), dimension(2, 2, 2), intent(in)  :: dist, fluid_mask
        real(wp), dimension(2, 2, 2), intent(out) :: interp_coeffs
        real(wp), dimension(2, 2, 2)              :: weights
        real(wp), dimension(2, 3)                 :: axis_weight  !< (stencil side, dimension) linear weights
        real(wp)                                  :: weight_sum
        integer                                   :: ii, jj, kk, num_z_layers

        num_z_layers = min(p, 1) + 1

        if (ib_second_order_vel) then
            ! bilinear weights, Mittal et al. (2008)
            axis_weight(2, 1) = (gp%ip_loc(1) - x_cc(i))/(x_cc(i + 1) - x_cc(i))
            axis_weight(2, 2) = (gp%ip_loc(2) - y_cc(j))/(y_cc(j + 1) - y_cc(j))
            axis_weight(2, 3) = 0._wp
            if (p /= 0) axis_weight(2, 3) = (gp%ip_loc(3) - z_cc(k))/(z_cc(k + 1) - z_cc(k))
            axis_weight(2,:) = min(max(axis_weight(2,:), 0._wp), 1._wp)  ! clamp to the cell edges
            axis_weight(1,:) = 1._wp - axis_weight(2,:)
            do kk = 1, num_z_layers
                do jj = 1, 2
                    do ii = 1, 2
                        weights(ii, jj, kk) = axis_weight(ii, 1)*axis_weight(jj, 2)*axis_weight(kk, 3)
                    end do
                end do
            end do
        else
            ! inverse distance squared, Chaudhuri et al. JCP (2011)
            if (minval(dist(:,:,1:num_z_layers)) <= 1.e-16_wp) then
                weights(:,:,1:num_z_layers) = merge(1._wp, 0._wp, dist(:,:,1:num_z_layers) <= 1.e-16_wp)
            else
                weights(:,:,1:num_z_layers) = 1._wp/dist(:,:,1:num_z_layers)**2
            end if
        end if

        interp_coeffs = 0._wp
        weights(:,:,1:num_z_layers) = fluid_mask(:,:,1:num_z_layers)*weights(:,:,1:num_z_layers)
        weight_sum = sum(weights(:,:,1:num_z_layers))
        if (weight_sum > 0._wp) interp_coeffs(:,:,1:num_z_layers) = weights(:,:,1:num_z_layers)/weight_sum

    end subroutine s_interpolation_weights

    !> Interpolate primitive variables to a ghost point's image point with the precomputed interp_coeffs: inverse-distance weights
    !! over the 2x2(x2) cell (ib_ip_interp = 1) or TSC weights over 3x3(x3) cells (ib_ip_interp = 2)
    subroutine s_interpolate_image_point(q_prim_vf, gp, alpha_rho_IP, alpha_IP, pres_IP, vel_IP, c_IP, r_IP, v_IP, pb_IP, mv_IP, &
                                         & nmom_IP, pb_in, mv_in, presb_IP, massv_IP, Ys_IP)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf  !< Primitive Variables
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_in, mv_in
        type(ghost_point), intent(in) :: gp
        real(wp), intent(inout) :: pres_IP
        real(wp), dimension(3), intent(inout) :: vel_IP
        real(wp), intent(inout) :: c_IP
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(inout) :: alpha_IP, alpha_rho_IP
        #:else
            real(wp), dimension(num_fluids), intent(inout) :: alpha_IP, alpha_rho_IP
        #:endif
        real(wp), optional, dimension(:), intent(inout) :: r_IP, v_IP, pb_IP, mv_IP
        real(wp), optional, dimension(:), intent(inout) :: nmom_IP
        real(wp), optional, dimension(:), intent(inout) :: presb_IP, massv_IP
        real(wp), optional, dimension(:), intent(inout) :: Ys_IP  !< Interpolated species mass fractions (chemistry)
        integer                                         :: i, j, k, l, q  !< Iterator variables
        integer                                         :: i1, i2, j1, j2, k1, k2  !< Iterator variables
        real(wp)                                        :: coeff

        i1 = gp%ip_stencil(1); i2 = i1 + ib_ip_width - 1
        j1 = gp%ip_stencil(2); j2 = j1 + ib_ip_width - 1
        k1 = gp%ip_stencil(3); k2 = k1 + ib_ip_width - 1

        if (p == 0) then
            k1 = 0
            k2 = 0
        end if

        alpha_rho_IP = 0._wp
        alpha_IP = 0._wp
        pres_IP = 0._wp
        vel_IP = 0._wp

        if (chemistry) Ys_IP = 0._wp

        if (surface_tension) c_IP = 0._wp

        if (bubbles_euler) then
            r_IP = 0._wp
            v_IP = 0._wp
            if (.not. polytropic) then
                mv_IP = 0._wp
                pb_IP = 0._wp
            end if
        end if

        if (qbmm) then
            nmom_IP = 0._wp
            if (.not. polytropic) then
                presb_IP = 0._wp
                massv_IP = 0._wp
            end if
        end if

        $:GPU_LOOP(parallelism='[seq]')
        do i = i1, i2
            $:GPU_LOOP(parallelism='[seq]')
            do j = j1, j2
                $:GPU_LOOP(parallelism='[seq]')
                do k = k1, k2
                    coeff = gp%interp_coeffs(i - i1 + 1, j - j1 + 1, k - k1 + 1)

                    pres_IP = pres_IP + coeff*q_prim_vf(eqn_idx%E)%sf(i, j, k)

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = eqn_idx%mom%beg, eqn_idx%mom%end
                        vel_IP(q + 1 - eqn_idx%mom%beg) = vel_IP(q + 1 - eqn_idx%mom%beg) + coeff*q_prim_vf(q)%sf(i, j, k)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do l = eqn_idx%cont%beg, eqn_idx%cont%end
                        alpha_rho_IP(l) = alpha_rho_IP(l) + coeff*q_prim_vf(l)%sf(i, j, k)
                        alpha_IP(l) = alpha_IP(l) + coeff*q_prim_vf(eqn_idx%adv%beg + l - 1)%sf(i, j, k)
                    end do

                    if (surface_tension) then
                        c_IP = c_IP + coeff*q_prim_vf(eqn_idx%c)%sf(i, j, k)
                    end if

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, num_species
                            Ys_IP(q) = Ys_IP(q) + coeff*q_prim_vf(eqn_idx%species%beg + q - 1)%sf(i, j, k)
                        end do
                    end if

                    if (bubbles_euler .and. .not. qbmm) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do l = 1, nb
                            if (polytropic) then
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + (l - 1)*2)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 1 + (l - 1)*2)%sf(i, j, k)
                            else
                                r_IP(l) = r_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + (l - 1)*4)%sf(i, j, k)
                                v_IP(l) = v_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 1 + (l - 1)*4)%sf(i, j, k)
                                pb_IP(l) = pb_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 2 + (l - 1)*4)%sf(i, j, k)
                                mv_IP(l) = mv_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg + 3 + (l - 1)*4)%sf(i, j, k)
                            end if
                        end do
                    end if

                    if (qbmm) then
                        do l = 1, nb*nmom
                            nmom_IP(l) = nmom_IP(l) + coeff*q_prim_vf(eqn_idx%bub%beg - 1 + l)%sf(i, j, k)
                        end do
                        if (.not. polytropic) then
                            do q = 1, nb
                                do l = 1, nnode
                                    presb_IP((q - 1)*nnode + l) = presb_IP((q - 1)*nnode + l) + coeff*real(pb_in(i, j, k, l, q), &
                                             & kind=wp)
                                    massv_IP((q - 1)*nnode + l) = massv_IP((q - 1)*nnode + l) + coeff*real(mv_in(i, j, k, l, q), &
                                             & kind=wp)
                                end do
                            end do
                        end if
                    end if
                end do
            end do
        end do

    end subroutine s_interpolate_image_point

    !> Resets the current indexes of immersed boundaries and replaces them after updating
    !> the position of each moving immersed boundary
    !> Recalculate the moving bodies' rotation matrices from their current angles (device).
    subroutine s_update_ib_rotation_matrices(num_ibs)

        integer, intent(in) :: num_ibs
        integer             :: i

        $:GPU_PARALLEL_LOOP(private='[i]')
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm /= 0) then
                call s_update_ib_rotation_matrix(i)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_update_ib_rotation_matrices

    impure subroutine s_update_mib(num_ibs)

        integer, intent(in) :: num_ibs
        integer             :: i, j, k, z_buff_size

        call nvtxStartRange("UPDATE-MIBM")

        ! Clears the existing immersed boundary indices over the same halo extent s_get_bounding_indices draws them
        z_buff_size = 0; if (p /= 0) z_buff_size = buff_size
        $:GPU_PARALLEL_LOOP(private='[i, j, k]')
        do i = -buff_size, m + buff_size; do j = -buff_size, n + buff_size; do k = -z_buff_size, p + z_buff_size
            ib_markers%sf(i, j, k) = 0._wp
            corrected_gps%sf(i, j, k) = 0
        end do; end do; end do
        $:END_GPU_PARALLEL_LOOP()

        call s_update_ib_rotation_matrices(num_ibs)

        ! recompute the new ib_patch locations and broadcast them.
        call nvtxStartRange("APPLY-IB-PATCHES")
        call s_apply_ib_patches(ib_markers)
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-GHOST-POINTS")
        ! recalculate the ghost point locations and coefficients
        call s_find_num_ghost_points(num_gps)
        ! the ghost_points capacity (a setup-time heuristic) can be outgrown when the surface's discrete cell count increases; the
        ! fill below has no bound check. size(ghost_points) is the active array's capacity (the fine slot's when swapped in).
        @:PROHIBIT(num_gps > size(ghost_points), &
                   & "moving IB: the ghost-point count outgrew the ghost-point array capacity; the body's surface-cell count increased beyond the setup-time sizing")
        ! num_gps is GPU_DECLARE'd and the on-device sort reads the device copy: mirror the init path's update
        $:GPU_UPDATE(device='[num_gps]')
        call s_find_ghost_points()
        call nvtxEndRange

        call nvtxStartRange("COMPUTE-IMAGE-POINTS")
        call s_apply_levelset(ghost_points, num_gps)
        call s_compute_image_points()
        call s_compute_interpolation_coeffs()
        call nvtxEndRange

        call nvtxEndRange

    end subroutine s_update_mib

    !> log(cosh(x)) without overflow for large |x|
    pure function f_log_cosh(x) result(y)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: x
        real(wp)             :: y
        y = abs(x) + log(0.5_wp*(1._wp + exp(-2._wp*abs(x))))

    end function f_log_cosh

    !> Prescribed kinematics for patch i at time t: hinged flapping (kin_model = 1) or the Eldredge pitch ramp (kin_model = 2).
    !! Flapping rolls phi about the lab x axis through the hinge and pitches theta about the body spanwise (y) axis through it; the
    !! ramp holds phi at zero and drives theta alone. Both compose as R = Rx(phi) Ry(theta) (MFC's angle convention) and set the
    !! angles, centroid, velocity and lab-frame angular velocity; nothing is integrated, so restarts are exact.
    subroutine s_prescribed_kinematics(i, t)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)    :: i
        real(wp), intent(in)   :: t
        real(wp)               :: tau, amp, damp, omega, phi, theta, phid, thetad, arg_phi, arg_th, t_p, a_s
        real(wp), dimension(3) :: r, c, w, v

        tau = t - patch_ib(i)%kin_t0
        omega = 2._wp*pi*patch_ib(i)%kin_freq
        amp = 0._wp
        damp = 0._wp
        if (tau >= 0._wp) then
            amp = 1._wp
            if (patch_ib(i)%kin_ramp > 0._wp .and. tau < patch_ib(i)%kin_ramp) then
                amp = 0.5_wp*(1._wp - cos(pi*tau/patch_ib(i)%kin_ramp))
                damp = 0.5_wp*pi/patch_ib(i)%kin_ramp*sin(pi*tau/patch_ib(i)%kin_ramp)
            end if
        end if
        arg_phi = omega*tau
        arg_th = omega*tau + patch_ib(i)%kin_phase

        if (patch_ib(i)%kin_model == 2) then
            ! Eldredge smoothed linear pitch ramp and hold (AIAA canonical pitch-up): with a = kin_smooth, nominal rate
            ! Omega = kin_pitch_rate and pitch time t_p = theta0/Omega, G = log[cosh(a tau)/cosh(a (tau - t_p))] rises from
            ! -a t_p to +a t_p, so theta = mean + theta0/2 [1 + G/(a t_p)] and theta' = Omega/2 [tanh(a tau) - tanh(a (tau - t_p))]
            phi = 0._wp
            phid = 0._wp
            t_p = patch_ib(i)%kin_theta0/patch_ib(i)%kin_pitch_rate
            a_s = patch_ib(i)%kin_smooth
            theta = patch_ib(i)%kin_theta_mean + 0.5_wp*patch_ib(i)%kin_theta0*(1._wp + (f_log_cosh(a_s*tau) &
                             & - f_log_cosh(a_s*(tau - t_p)))/(a_s*t_p))
            thetad = 0.5_wp*patch_ib(i)%kin_pitch_rate*(tanh(a_s*tau) - tanh(a_s*(tau - t_p)))
        else
            phi = amp*patch_ib(i)%kin_phi0*sin(arg_phi)
            phid = damp*patch_ib(i)%kin_phi0*sin(arg_phi) + amp*patch_ib(i)%kin_phi0*omega*cos(arg_phi)
            theta = patch_ib(i)%kin_theta_mean + amp*patch_ib(i)%kin_theta0*sin(arg_th)
            thetad = damp*patch_ib(i)%kin_theta0*sin(arg_th) + amp*patch_ib(i)%kin_theta0*omega*cos(arg_th)
        end if

        ! centroid relative to the hinge: Rx(phi) Ry(theta) offset. kin_offset points to the case-file centroid; a
        ! centroid moved to the centre of mass (geometries 4, 5, 11, 12) sits centroid_offset short of it
        r = patch_ib(i)%kin_offset - patch_ib(i)%centroid_offset
        c(1) = cos(theta)*r(1) + sin(theta)*r(3)
        c(2) = r(2)
        c(3) = -sin(theta)*r(1) + cos(theta)*r(3)
        v(1) = c(1)
        v(2) = cos(phi)*c(2) - sin(phi)*c(3)
        v(3) = sin(phi)*c(2) + cos(phi)*c(3)
        c = v

        ! lab-frame angular velocity: phi' e_x + theta' Rx(phi) e_y
        w(1) = phid
        w(2) = thetad*cos(phi)
        w(3) = thetad*sin(phi)
        call s_cross_product(w, c, v)

        patch_ib(i)%angles(1) = phi
        patch_ib(i)%angles(2) = theta
        patch_ib(i)%angles(3) = 0._wp
        patch_ib(i)%x_centroid = patch_ib(i)%kin_hinge(1) + c(1)
        patch_ib(i)%y_centroid = patch_ib(i)%kin_hinge(2) + c(2)
        patch_ib(i)%z_centroid = patch_ib(i)%kin_hinge(3) + c(3)
        patch_ib(i)%vel = v
        patch_ib(i)%angular_vel = w

    end subroutine s_prescribed_kinematics

    !> Face flux of component c. The dual-pass HLLD RHS takes the hat_R set at the left face of a cell (hatR) and hat_L at its
    !! right.
    function f_ib_face_flux(hatR, j, k, l, c) result(f)

        $:GPU_ROUTINE(function_name='f_ib_face_flux', parallelism='[seq]', cray_inline=True)
        logical, intent(in) :: hatR
        integer, intent(in) :: j, k, l, c
        real(wp)            :: f

        if (hatR) then
            f = flux_hatR_rsx_vf(j, k, l, c)
        else
            f = flux_rsx_vf(j, k, l, c)
        end if

    end function f_ib_face_flux

    !> Mass flux through face (j, k, l)
    function f_ib_face_mass(hatR, j, k, l) result(f_m)

        $:GPU_ROUTINE(function_name='f_ib_face_mass', parallelism='[seq]', cray_inline=True)
        logical, intent(in) :: hatR
        integer, intent(in) :: j, k, l
        real(wp)            :: f_m
        integer             :: i

        f_m = 0._wp
        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            f_m = f_m + f_ib_face_flux(hatR, j, k, l, i)
        end do

    end function f_ib_face_mass

    !> Momentum density of component i in cell (j, k, l). The momentum halos of q_cons_vf are not exchanged between ranks, so a halo
    !! cell (rank seam, periodic edge) takes rho*u from q_prim_vf, whose halos are.
    function f_ib_cell_momentum(q_cons_vf, q_prim_vf, i, j, k, l) result(mom)

        $:GPU_ROUTINE(function_name='f_ib_cell_momentum', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        integer, intent(in)                                 :: i, j, k, l
        real(wp)                                            :: mom, rho
        integer                                             :: c

        mom = q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)
        if (igr .or. (j >= 0 .and. j <= m .and. k >= 0 .and. k <= n .and. l >= 0 .and. l <= p)) return
        rho = 0._wp
        do c = eqn_idx%cont%beg, eqn_idx%cont%end
            rho = rho + q_prim_vf(c)%sf(j, k, l)
        end do
        mom = rho*q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)

    end function f_ib_cell_momentum

    !> Area of the face (j, k, l) normal to direction idir
    pure real(wp) function f_ib_face_area(idir, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: idir, j, k, l

        if (idir == 1) f_ib_face_area = dy(k)
        if (idir == 2) f_ib_face_area = dx(j)
        if (idir == 3) f_ib_face_area = dx(j)*dy(k)
        if (num_dims == 3 .and. idir /= 3) f_ib_face_area = f_ib_face_area*dz(l)

    end function f_ib_face_area

    !> Weight of RK stage `stage` in the step (SSP-RK: 1; 1/2, 1/2; 1/6, 1/6, 2/3)
    pure real(wp) function f_ib_stage_weight(stage)
        integer, intent(in)               :: stage
        real(wp), dimension(3), parameter :: w_rk3 = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

        f_ib_stage_weight = 1._wp
        if (time_stepper == time_stepper_rk2) f_ib_stage_weight = 0.5_wp
        if (time_stepper == time_stepper_rk3) f_ib_stage_weight = w_rk3(stage)

    end function f_ib_stage_weight

    !> The mass part of s_accumulate_ib_face_fluxes alone: ib_mass_relax on fixed bodies with no forces, heat output, lumped body,
    !! AMR or dual-pass HLLD. A kernel's launch cost grows with the arrays it names (about 0.5 ms on an MI300A for the general one),
    !! so this one names only what the mass sum needs. It reads flux_rsx_vf itself: f_ib_face_mass, which NVHPC did not inline here,
    !! returned the fluxes of the first RHS evaluation at every step on an A100 (#2414).
    subroutine s_accumulate_ib_mass(idir, stage)

        integer, intent(in) :: idir, stage
        real(wp)            :: w, area, f_m
        integer             :: i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g, sh, j0, k0, l0

        w = f_ib_stage_weight(stage)
        j0 = 0; k0 = 0; l0 = 0
        if (idir == 1) j0 = -1
        if (idir == 2) k0 = -1
        if (idir == 3) l0 = -1

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g, sh, &
                            & area, f_m]', firstprivate='[w]')
        do l = l0, p
            do k = k0, n
                do j = j0, m
                    jn = j; kn = k; ln = l
                    if (idir == 1) jn = j + 1
                    if (idir == 2) kn = k + 1
                    if (idir == 3) ln = l + 1
                    side = 0
                    if (jn <= m .and. kn <= n .and. ln <= p) then
                        if (ib_markers%sf(j, k, l) == 0 .and. ib_markers%sf(jn, kn, ln) /= 0) then
                            if (j >= 0 .and. k >= 0 .and. l >= 0) side = 1
                        else if (ib_markers%sf(j, k, l) /= 0 .and. ib_markers%sf(jn, kn, ln) == 0) then
                            side = -1
                        end if
                    else if (ib_markers%sf(j, k, l) == 0 .and. ib_markers%sf(jn, kn, ln) /= 0) then
                        side = 1
                    end if
                    if (side /= 0) then
                        ci = jn; cj = kn; ck = ln
                        if (side == -1) then
                            ci = j; cj = k; ck = l
                        end if
                        call s_decode_patch_periodicity(ib_markers%sf(ci, cj, ck), ib_tmp, xp, yp, zp)
                        call s_get_neighborhood_idx(ib_tmp, ib_id)
                        if (ib_id > 0) then
                            area = f_ib_face_area(idir, j, k, l)
                            f_m = 0._wp
                            do i = eqn_idx%cont%beg, eqn_idx%cont%end
                                f_m = f_m + flux_rsx_vf(j, k, l, i)
                            end do
                            g = ib_id
                            sh = f_ib_shard(j, k, l)
                            $:GPU_ATOMIC(atomic='update')
                            ib_part(2, sh, g) = ib_part(2, sh, g) - real(side, wp)*w*area*f_m
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_accumulate_ib_mass

    !> Add direction idir's fluid-body face fluxes to each IB. These are the fluxes the RHS differences, Riemann (flux_rsx_vf) plus
    !! viscous, capillary, conductive and diffusive (flux_src_rsx_vf), so each body receives exactly what the fluid loses. A face
    !! counts on the rank owning its fluid cell. With IB forces: force (1:3) and torque about the centroid (4:6) of each local IB
    !! into ib_face_ft, from the last RHS evaluation (idir = 1 starts afresh); a moving body's momentum flux is taken through a face
    !! moving with its wall (face momentum times the wall's normal velocity is removed). With ib_exchange: energy in and mass out of
    !! each lumped body (every body with ib_heat_wrt) into ib_heat(2:3), weighted so a step's stages sum to its exchange (SSP-RK: 1;
    !! 1/2, 1/2; 1/6, 1/6, 2/3).
    subroutine s_accumulate_ib_face_fluxes(idir, stage, q_cons_vf, q_prim_vf, markers)

        integer, intent(in) :: idir, stage
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        type(integer_field), intent(in) :: markers  !< ib_markers, or ib_markers_fine on an AMR fine frame
        real(wp), dimension(3) :: f, r, tq, vb
        real(wp) :: area, w, f_E, f_m, f_mom
        integer :: i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g
        integer :: j0, k0, l0, fj, fk, fl, ibm, jl, kl, ll, sd, bw, cov, nc, sh, c
        logical :: has_src, has_src_E, do_ft, all_ibs, dual_pass, hatR, fine, heat, covd

        do_ft = moving_immersed_boundary_flag .or. ib_state_wrt .or. ib_force_wrt
        if (num_ibs == 0 .or. .not. (do_ft .or. ib_exchange)) return
        if (ib_mass_only) then
            call s_accumulate_ib_mass(idir, stage)
            return
        end if
        all_ibs = ib_heat_wrt .or. ib_mass_relax > 0
        ! AMR: each face counts on the finest level covering its fluid cell. The coarse pass of a stage starts the sums and skips
        ! covered cells (cov = 1, the cell map); a fine frame adds its faces, skipping those under its own children (cov = 2, their
        ! boxes). A fine frame takes interior faces only (its markers have no ghosts; bodies sit inside their blocks) and no
        ! heat exchange (thermal walls are rejected with amr). Subcycled levels are not taken: their passes span several stages.
        fine = amr_in_fine_advance
        heat = ib_exchange .and. .not. fine
        cov = 0
        if (amr .and. .not. amr_subcycle) cov = merge(2, 1, fine)
        nc = merge(amr_cov_n, 0, cov == 2); sd = amr_bat_sd; bw = max(amr_bat_w, 1)
        if (fine .and. .not. (do_ft .and. cov == 2)) return

        if (do_ft .and. idir == 1 .and. .not. fine) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', collapse=2)
            do i = 1, num_ibs
                do l = 1, 6
                    ib_face_ft(l, i) = 0._wp
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! The conditions under which s_compute_additional_physics_rhs differences flux_src for the momentum and the energy
        has_src = viscous .or. surface_tension .or. heat_conduction
        has_src_E = has_src .or. chem_params%diffusion
        dual_pass = hypo_nc_mode == hypo_nc_mode_dual_pass
        w = f_ib_stage_weight(stage)
        j0 = 0; k0 = 0; l0 = 0
        if (.not. fine) then
            if (idir == 1) j0 = -1
            if (idir == 2) k0 = -1
            if (idir == 3) l0 = -1
        end if

        ! Face (j, k, l) lies between that cell and the next one along idir
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jn, kn, ln, ci, cj, ck, side, ib_id, ib_tmp, xp, yp, zp, g, area, &
                            & f, r, tq, vb, f_E, f_m, f_mom, hatR, fj, fk, fl, ibm, jl, kl, ll, covd, sh]', firstprivate='[fine, &
                            & heat, cov, nc, sd, bw]')
        do l = l0, p
            do k = k0, n
                do j = j0, m
                    jn = j; kn = k; ln = l
                    if (idir == 1) jn = j + 1
                    if (idir == 2) kn = k + 1
                    if (idir == 3) ln = l + 1
                    ! side = +1: body on the far side, the face flux enters it; -1: body on the near side, the flux leaves it
                    side = 0
                    if (jn <= m .and. kn <= n .and. ln <= p) then
                        if (markers%sf(j, k, l) == 0 .and. markers%sf(jn, kn, ln) /= 0) then
                            if (j >= 0 .and. k >= 0 .and. l >= 0) side = 1
                        else if (markers%sf(j, k, l) /= 0 .and. markers%sf(jn, kn, ln) == 0) then
                            side = -1
                        end if
                    else if (.not. fine) then
                        if (markers%sf(j, k, l) == 0 .and. markers%sf(jn, kn, ln) /= 0) side = 1
                    end if
                    if (side /= 0) then
                        hatR = dual_pass .and. side == -1
                        ci = jn; cj = kn; ck = ln
                        fj = j; fk = k; fl = l
                        if (side == -1) then
                            ci = j; cj = k; ck = l
                            fj = jn; fk = kn; fl = ln
                        end if
                        ! the fluid cell (fj, fk, fl) under a finer level: that level takes the face
                        covd = .false.
                        if (cov == 1) then
                            covd = amr_cfl_rf(fj, fk, fl) > 1._wp
                        else
                            do i = 1, nc
                                if (fj >= amr_cov_box(1, i) .and. fj <= amr_cov_box(4, i) .and. fk >= amr_cov_box(2, &
                                    & i) .and. fk <= amr_cov_box(5, i) .and. fl >= amr_cov_box(3, i) .and. fl <= amr_cov_box(6, &
                                    & i)) covd = .true.
                            end do
                        end if
                        if (covd) side = 0
                    end if
                    if (side /= 0) then
                        call s_decode_patch_periodicity(markers%sf(ci, cj, ck), ib_tmp, xp, yp, zp)
                        call s_get_neighborhood_idx(ib_tmp, ib_id)
                        if (ib_id > 0) then
                            area = f_ib_face_area(idir, j, k, l)
                            f_m = f_ib_face_mass(hatR, j, k, l)
                            if (do_ft) then
                                ! Lever arm: face centre minus the centroid of the body's periodic image
                                r = 0._wp
                                if (fine) then
                                    ! a batch member reads the leader's coordinates: its own are those at its local index, shifted
                                    ibm = 1 + merge(j, merge(k, l, sd == 2), sd == 1)/bw
                                    jl = j - amr_bat_off(1, ibm); kl = k - amr_bat_off(2, ibm); ll = l - amr_bat_off(3, ibm)
                                    r(1) = merge(x_cb(jl), x_cc(jl), idir == 1) + amr_bat_xsh(1, &
                                      & ibm) - (patch_ib(ib_id)%x_centroid + real(xp, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg))
                                    r(2) = merge(y_cb(kl), y_cc(kl), idir == 2) + amr_bat_xsh(2, &
                                      & ibm) - (patch_ib(ib_id)%y_centroid + real(yp, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg))
                                    if (num_dims == 3) r(3) = merge(z_cb(ll), z_cc(ll), idir == 3) + amr_bat_xsh(3, &
                                        & ibm) - (patch_ib(ib_id)%z_centroid + real(zp, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg))
                                else
                                    r(1) = merge(x_cb(j), x_cc(j), idir == 1) - (patch_ib(ib_id)%x_centroid + real(xp, &
                                      & wp)*(glb_bounds(1)%end - glb_bounds(1)%beg))
                                    r(2) = merge(y_cb(k), y_cc(k), idir == 2) - (patch_ib(ib_id)%y_centroid + real(yp, &
                                      & wp)*(glb_bounds(2)%end - glb_bounds(2)%beg))
                                    if (num_dims == 3) r(3) = merge(z_cb(l), z_cc(l), &
                                        & idir == 3) - (patch_ib(ib_id)%z_centroid + real(zp, &
                                        & wp)*(glb_bounds(3)%end - glb_bounds(3)%beg))
                                end if
                                ! Wall velocity at the face, V_b = V + omega x r
                                call s_cross_product(patch_ib(ib_id)%angular_vel, r, vb)
                                vb = vb + patch_ib(ib_id)%vel
                                f = 0._wp
                                do i = 1, num_dims
                                    f(i) = f_ib_face_flux(hatR, j, k, l, eqn_idx%mom%beg + i - 1)
                                    if (has_src) f(i) = f(i) + flux_src_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                    ! A moving wall's force is the momentum flux through a face moving with it, (rho u_i)(u_n - V_n)
                                    ! + p n_i: the lab flux carries rho u_i u_n, so remove the face momentum (mean of the fluid
                                    ! and body cell) times the wall's normal velocity V_n
                                    if (patch_ib(ib_id)%moving_ibm > 0) then
                                        f(i) = f(i) - 0.5_wp*vb(idir)*(f_ib_cell_momentum(q_cons_vf, q_prim_vf, i, fj, fk, &
                                          & fl) + f_ib_cell_momentum(q_cons_vf, q_prim_vf, i, ci, cj, ck))
                                    end if
                                    f(i) = real(side, wp)*area*f(i)
                                end do
                                call s_cross_product(r, f, tq)
                                do i = 1, 3
                                    $:GPU_ATOMIC(atomic='update')
                                    ib_face_ft(i, ib_id) = ib_face_ft(i, ib_id) + f(i)
                                    $:GPU_ATOMIC(atomic='update')
                                    ib_face_ft(i + 3, ib_id) = ib_face_ft(i + 3, ib_id) + tq(i)
                                end do
                            end if
                            if (heat .and. (all_ibs .or. patch_ib(ib_id)%thermal_bc == 3)) then
                                f_E = f_ib_face_flux(hatR, j, k, l, eqn_idx%E)
                                if (has_src_E) f_E = f_E + flux_src_rsx_vf(j, k, l, eqn_idx%E)
                                if (patch_ib(ib_id)%moving_ibm > 0) then
                                    ! Energy flux in the body frame, f_E - V_b.f_mom + |V_b|^2 f_m/2: the lab flux also carries the
                                    ! wall's pressure/viscous work and the energy of the gas it sweeps (that is booked with the
                                    ! swept cells), neither of which is heat into the body
                                    do i = 1, num_dims
                                        f_mom = f_ib_face_flux(hatR, j, k, l, eqn_idx%mom%beg + i - 1)
                                        if (has_src) f_mom = f_mom + flux_src_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                        f_E = f_E - vb(i)*f_mom
                                    end do
                                    f_E = f_E + 0.5_wp*(vb(1)**2 + vb(2)**2 + vb(3)**2)*f_m
                                end if
                                g = ib_id
                                sh = f_ib_shard(j, k, l)
                                $:GPU_ATOMIC(atomic='update')
                                ib_part(1, sh, g) = ib_part(1, sh, g) + real(side, wp)*w*area*f_E
                                $:GPU_ATOMIC(atomic='update')
                                ib_part(2, sh, g) = ib_part(2, sh, g) - real(side, wp)*w*area*f_m
                            end if
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Wall columns: the body books the column's flux, i.e. the gas faces' flux plus the guard's cut of the gas drain (once per
        ! stage)
        if (ib_wall_col .and. all_ibs .and. idir == 1 .and. n_wc > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, c, g, sh]', copyin='[w]', present='[ghost_points]')
            do c = 1, n_wc
                ! unconditional: a private scalar set only under an if is undefined on NVHPC GPUs
                i = max(wc_gp(c), 1)
                g = max(ghost_points(i)%ib_patch_id, 0)  ! local body slot, as ib_part is indexed since #2330
                sh = f_ib_shard(ghost_points(i)%loc(1), ghost_points(i)%loc(2), ghost_points(i)%loc(3))
                if (wc_on(c) == 1) then
                    $:GPU_ATOMIC(atomic='update')
                    ib_part(1, sh, g) = ib_part(1, sh, g) + w*wc_cut(c)
                end if
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_accumulate_ib_face_fluxes

    !> Set each IB's force and torque: the fluid's share from the face fluxes of the last RHS evaluation
    !! (s_accumulate_ib_face_fluxes), plus collisions and body forces
    subroutine s_compute_ib_forces()

        integer                         :: i, l
        real(wp), dimension(num_ibs, 3) :: forces, torques

        call nvtxStartRange("COMPUTE-IB-FORCES")

        forces = 0._wp
        torques = 0._wp

        ! no kernel over zero-size mapped arrays: a rank can hold no patch after a handoff
        if (num_ibs > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', copy='[forces, torques]')
            do i = 1, num_ibs
                do l = 1, 3
                    forces(i, l) = ib_face_ft(l, i)
                    torques(i, l) = ib_face_ft(l + 3, i)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call nvtxStartRange("APPLY-COLLISION-FORCES")
        call s_apply_IB_collisions(ghost_points, num_gps, ib_markers, forces, torques)
        call nvtxEndRange

        ! reduce the forces across local neighborhood ranks
        call nvtxStartRange("COMMUNICATE-IB-FORCES")
        call s_communicate_ib_forces(forces, torques)
        call nvtxEndRange

        ! Body forces after reducing, to avoid double counting. Only a two-way body feels them; a fixed or prescribed body keeps
        ! the fluid and collision force, and its mass may be unset (-1e6)
        do i = 1, num_ibs
            if (patch_ib(i)%moving_ibm >= 2) then
                if (bf_x) forces(i, 1) = forces(i, 1) + accel_bf(1)*patch_ib(i)%mass
                if (bf_y) forces(i, 2) = forces(i, 2) + accel_bf(2)*patch_ib(i)%mass
                if (bf_z) forces(i, 3) = forces(i, 3) + accel_bf(3)*patch_ib(i)%mass
            end if
        end do

        ! apply the summed forces
        if (num_ibs > 0) then
            $:GPU_PARALLEL_LOOP(private='[i, l]', copyin='[forces, torques]')
            do i = 1, num_ibs
                do l = 1, 3
                    patch_ib(i)%force(l) = forces(i, l)
                    patch_ib(i)%torque(l) = torques(i, l)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call nvtxEndRange

    end subroutine s_compute_ib_forces

    !> Computes the center of mass for IB patch types where we are unable to determine their center of mass analytically.
    !> These patches include things like NACA airfoils and STL models
    !> Move a moving airfoil/STL patch's centroid to its centre of mass, keeping the difference in centroid_offset. Collective:
    !! every rank calls it for every global id, contributing zeros for patches it does not hold.
    subroutine s_compute_centroid_offset(gid, ib_marker)

        integer, intent(in)      :: gid        !< global patch id
        integer, intent(in)      :: ib_marker  !< local index on this rank; <= 0 if not held
        integer                  :: i, j, k, num_cells_local, decoded_gbl_id
        integer(kind=8)          :: num_cells, needs_loc, needs_glb
        real(wp), dimension(1:3) :: center_of_mass, center_of_mass_local

        needs_loc = 0_8
        if (ib_marker > 0) then
            if (f_needs_centroid_offset(patch_ib(ib_marker))) needs_loc = 1_8
        end if
        call s_mpi_allreduce_integer_sum(needs_loc, needs_glb)
        if (needs_glb == 0_8) then
            if (ib_marker > 0) patch_ib(ib_marker)%centroid_offset(:) = 0._wp
            return
        end if

        center_of_mass_local = [0._wp, 0._wp, 0._wp]
        num_cells_local = 0
        if (ib_marker > 0) then
            ! get the summed mass distribution and number of cells to divide by
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) /= 0) then
                            call s_decode_patch_periodicity(ib_markers%sf(i, j, k), decoded_gbl_id)
                            if (decoded_gbl_id == gid) then
                                num_cells_local = num_cells_local + 1
                                center_of_mass_local = center_of_mass_local + [x_cc(i), y_cc(j), 0._wp]
                                if (num_dims == 3) center_of_mass_local(3) = center_of_mass_local(3) + z_cc(k)
                            end if
                        end if
                    end do
                end do
            end do
        end if

        ! reduce the mass contribution over all MPI ranks and compute COM
        call s_mpi_allreduce_integer_sum(int(num_cells_local, 8), num_cells)
        call s_mpi_allreduce_sum(center_of_mass_local(1), center_of_mass(1))
        call s_mpi_allreduce_sum(center_of_mass_local(2), center_of_mass(2))
        call s_mpi_allreduce_sum(center_of_mass_local(3), center_of_mass(3))
        if (ib_marker <= 0) return
        if (num_cells == 0) then
            patch_ib(ib_marker)%centroid_offset = [0._wp, 0._wp, 0._wp]
            return
        end if
        center_of_mass = center_of_mass/real(num_cells, wp)

        ! assign the centroid offset as a vector pointing from the true COM to the "centroid" in the input file and replace the
        ! current centroid
        patch_ib(ib_marker)%centroid_offset = [patch_ib(ib_marker)%x_centroid, patch_ib(ib_marker)%y_centroid, &
                 & patch_ib(ib_marker)%z_centroid] - center_of_mass
        patch_ib(ib_marker)%x_centroid = center_of_mass(1)
        patch_ib(ib_marker)%y_centroid = center_of_mass(2)
        patch_ib(ib_marker)%z_centroid = center_of_mass(3)

        ! rotate the centroid offset back into the local coords of the IB
        patch_ib(ib_marker)%centroid_offset = matmul(patch_ib(ib_marker)%rotation_matrix_inverse, &
                 & patch_ib(ib_marker)%centroid_offset)

    end subroutine s_compute_centroid_offset

    !> A moving airfoil or STL (geometries 4, 5, 11, 12) moves its centroid to the centre of mass and keeps a centroid_offset
    pure logical function f_needs_centroid_offset(patch)

        type(ib_patch_parameters), intent(in) :: patch

        f_needs_centroid_offset = patch%moving_ibm /= 0 .and. any(patch%geometry == [ib_airfoil_2d, ib_model_2d, ib_airfoil_3d, &
            & ib_model_3d])

    end function f_needs_centroid_offset

    !> On restart, replace the re-measured centroid offsets with the ones the run was using (restart_data/ib_offset_<step>.dat,
    !! written with each checkpoint), and re-place kinematics-driven bodies about them. No file: the measured offsets stand.
    impure subroutine s_restore_centroid_offsets(t_init)

        real(wp), intent(in)                 :: t_init
        character(len=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: gid, i, ios, file_unit, step
        real(wp), dimension(3)               :: off

        step = t_step_start
        if (cfl_dt) step = n_start
        if (step == 0) return
        write (file_loc, '(A,I0,A)') trim(case_dir) // '/restart_data/ib_offset_', step, '.dat'
        inquire (file=trim(file_loc), exist=file_exist)
        if (.not. file_exist) return
        open (newunit=file_unit, file=trim(file_loc), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (file_unit, *, iostat=ios) gid, off
            if (ios /= 0) exit
            call s_get_neighborhood_idx(gid, i)
            if (i > 0) then
                patch_ib(i)%centroid_offset = off
                if (patch_ib(i)%moving_ibm /= 0 .and. patch_ib(i)%kin_model > 0) call s_prescribed_kinematics(i, t_init)
            end if
        end do
        close (file_unit)

    end subroutine s_restore_centroid_offsets

    !> Computes the moment of inertia for an immersed boundary
    subroutine s_compute_moment_of_inertia(patch, axis, moment)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(ib_patch_parameters), intent(in) :: patch
        real(wp), dimension(3), intent(in)    :: axis
        real(wp), intent(out)                 :: moment
        real(wp), dimension(3)                :: normal_axis, n_b

        ! if the IB is in 2D or a 3D sphere, we can compute this exactly
        if (patch%geometry == ib_circle) then  ! circle
            moment = 0.5_wp*patch%mass*(patch%radius)**2
        else if (patch%geometry == ib_rectangle) then  ! rectangle
            moment = patch%mass*(patch%length_x**2 + patch%length_y**2)/12._wp
        else if (patch%geometry == ib_ellipse) then  ! ellipse
            moment = 0.0625_wp*patch%mass*(patch%length_x**2 + patch %length_y**2)
        else if (patch%geometry == ib_sphere) then  ! sphere
            moment = 0.4_wp*patch%mass*(patch%radius)**2
        else  ! no analytic moment: from the body-frame second moment summed at set-up (s_sum_ib_inertia_tensors)
            if (p == 0) then
                normal_axis = [0._wp, 0._wp, 1._wp]
            else if (sqrt(sum(axis**2)) < sgm_eps) then
                ! if the object is not actually rotating at this time, return a dummy value and exit
                moment = 1._wp
                return
            else
                normal_axis = axis/sqrt(sum(axis**2))
            end if
            ! sum over cells of the squared distance to the axis: trace(S) - n.S.n, with n in the body frame
            n_b(1) = dot_product(patch%rotation_matrix(:,1), normal_axis)
            n_b(2) = dot_product(patch%rotation_matrix(:,2), normal_axis)
            n_b(3) = dot_product(patch%rotation_matrix(:,3), normal_axis)
            moment = patch%mass*(patch%inertia_b(1)*(1._wp - n_b(1)**2) + patch%inertia_b(2)*(1._wp - n_b(2)**2) &
                                 & + patch%inertia_b(3)*(1._wp - n_b(3)**2) - 2._wp*(n_b(1)*n_b(2)*patch%inertia_b(4) + n_b(1) &
                                 & *n_b(3)*patch%inertia_b(5) + n_b(2)*n_b(3)*patch%inertia_b(6)))
            ! a body too small to hold more than one cell keeps the dummy value
            if (moment <= 0._wp) moment = 1._wp
        end if

    end subroutine s_compute_moment_of_inertia

    !> Second moment about the centroid, per unit mass, of each moving body without an analytic moment, in its body frame. Every
    !! rank sums its own cells (periodic images measured from their own centroid); the owner adds the partials and returns the
    !! totals (s_communicate_ib_forces), so every rank holding the body gets the same tensor. Once at set-up: the rotating body then
    !! needs no per-step cell sum, only n.R S R^T.n.
    impure subroutine s_sum_ib_inertia_tensors()

        real(wp), dimension(num_ibs, 3) :: s_diag, s_off, cnt, dummy
        real(wp), dimension(3)          :: r, centroid
        real(wp), dimension(3, 3)       :: s_w, s_b
        integer                         :: i, j, k, l, gbl_id, xp, yp, zp

        s_diag = 0._wp
        s_off = 0._wp
        cnt = 0._wp
        dummy = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib_markers%sf(j, k, l) == 0) cycle
                    call s_decode_patch_periodicity(ib_markers%sf(j, k, l), gbl_id, xp, yp, zp)
                    call s_get_neighborhood_idx(gbl_id, i)
                    if (i <= 0) cycle
                    if (.not. f_needs_cell_moment(patch_ib(i))) cycle
                    centroid(1) = patch_ib(i)%x_centroid + real(xp, wp)*(glb_bounds(1)%end - glb_bounds(1)%beg)
                    centroid(2) = patch_ib(i)%y_centroid + real(yp, wp)*(glb_bounds(2)%end - glb_bounds(2)%beg)
                    r = [x_cc(j) - centroid(1), y_cc(k) - centroid(2), 0._wp]
                    if (num_dims == 3) then
                        centroid(3) = patch_ib(i)%z_centroid + real(zp, wp)*(glb_bounds(3)%end - glb_bounds(3)%beg)
                        r(3) = z_cc(l) - centroid(3)
                    end if
                    s_diag(i,:) = s_diag(i,:) + r**2
                    s_off(i,:) = s_off(i,:) + [r(1)*r(2), r(1)*r(3), r(2)*r(3)]
                    cnt(i, 1) = cnt(i, 1) + 1._wp
                end do
            end do
        end do

        call s_communicate_ib_forces(s_diag, s_off)
        call s_communicate_ib_forces(cnt, dummy)

        do i = 1, num_ibs
            patch_ib(i)%inertia_b = 0._wp
            if (.not. f_needs_cell_moment(patch_ib(i)) .or. cnt(i, 1) < 1._wp) cycle
            s_w = reshape([s_diag(i, 1), s_off(i, 1), s_off(i, 2), s_off(i, 1), s_diag(i, 2), s_off(i, 3), s_off(i, 2), s_off(i, &
                          & 3), s_diag(i, 3)], [3, 3])/cnt(i, 1)
            ! rotation_matrix maps body-frame vectors to the lab frame, so S_body = R^T S_lab R
            s_b = matmul(transpose(patch_ib(i)%rotation_matrix), matmul(s_w, patch_ib(i)%rotation_matrix))
            patch_ib(i)%inertia_b = [s_b(1, 1), s_b(2, 2), s_b(3, 3), s_b(1, 2), s_b(1, 3), s_b(2, 3)]
        end do
        do i = 1, num_ibs
            $:GPU_UPDATE(device='[patch_ib(i)]')
        end do

    end subroutine s_sum_ib_inertia_tensors

    !> A moving body whose moment of inertia is not analytic (everything but circle, rectangle, ellipse and sphere)
    pure logical function f_needs_cell_moment(patch)

        type(ib_patch_parameters), intent(in) :: patch

        f_needs_cell_moment = patch%moving_ibm /= 0 .and. .not. any(patch%geometry == [ib_circle, ib_rectangle, ib_ellipse, &
            & ib_sphere])

    end function f_needs_cell_moment

    !> Wrap immersed boundary positions across periodic domain boundaries
    subroutine s_wrap_periodic_ibs()

        integer :: patch_id

        $:GPU_PARALLEL_LOOP(private='[patch_id]')
        do patch_id = 1, num_ibs
            ! check domain wraps in x, y,
            #:for X, DIR in [('x', 1), ('y', 2), ('z', 3)]
                if (num_dims >= ${DIR}$) then
                    ! check for periodicity
                    if (ib_bc_${X}$%beg == BC_PERIODIC) then
                        ! check if the boundary has left the domain, and then correct
                        if (patch_ib(patch_id)%${X}$_centroid < glb_bounds(${DIR}$)%beg) then
                            ! if the boundary exited "left", wrap it back around to the "right"
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid + (glb_bounds(${DIR}$)%end &
                                     & - glb_bounds(${DIR}$)%beg)
                        else if (patch_ib(patch_id)%${X}$_centroid >= glb_bounds(${DIR}$)%end) then
                            ! "right": the last rank owns [beg, end), so a centroid on the end belongs to no rank until wrapped
                            patch_ib(patch_id)%${X}$_centroid = patch_ib(patch_id)%${X}$_centroid - (glb_bounds(${DIR}$)%end &
                                     & - glb_bounds(${DIR}$)%beg)
                        end if
                    end if
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_wrap_periodic_ibs

    !> Builds the ascending list of distinct neighborhood ranks used by s_communicate_ib_forces and allocates its buffers once. The
    !! ascending order is the fixed order in which an owner adds its neighbors' contributions.
    impure subroutine s_setup_ib_force_comm()

#ifdef MFC_MPI
        integer                            :: i, j, r, max_recs
        integer, allocatable, dimension(:) :: flat

        if (num_procs == 1 .or. allocated(ib_nbrs)) return
        @:PROHIBIT(num_gbl_ibs >= min(2_8**digits(0._wp), int(huge(0), 8)), &
                   & "Too many IBs to encode their ids exactly in real(wp) force messages")

        ! a rank can fill several table slots (periodicity, few ranks) or be its own neighbor; keep each distinct rank once
        flat = reshape(ib_neighbor_ranks, [size(ib_neighbor_ranks)])
        allocate (ib_nbrs(size(flat)))
        n_ib_nbrs = 0
        do i = 1, size(flat)
            r = flat(i)
            if (r < 0 .or. r == proc_rank) cycle
            if (any(ib_nbrs(1:n_ib_nbrs) == r)) cycle
            ! insert in ascending order
            j = n_ib_nbrs
            do while (j > 0)
                if (ib_nbrs(j) < r) exit
                ib_nbrs(j + 1) = ib_nbrs(j)
                j = j - 1
            end do
            ib_nbrs(j + 1) = r
            n_ib_nbrs = n_ib_nbrs + 1
        end do

        ! an owner never receives more records than it owns, and a broadcast never carries more than its sender owns
        max_recs = min(size(patch_ib), num_local_ibs_max)
        allocate (ib_send_counts(n_ib_nbrs), ib_reqs(2*n_ib_nbrs), ib_stats(MPI_STATUS_SIZE, 2*n_ib_nbrs))
        allocate (ib_reply_slots(max_recs, max(n_ib_nbrs, 1)), ib_reply_counts(max(n_ib_nbrs, 1)))
        allocate (ib_send_bufs(ib_rec_len*max_recs, max(n_ib_nbrs, 1)), ib_recv_bufs(ib_rec_len*max_recs, max(n_ib_nbrs, 1)))
#endif

    end subroutine s_setup_ib_force_comm

    !> Position of rank in ib_nbrs
    function f_ib_nbr_slot(rank) result(slot)

        integer, intent(in) :: rank
        integer             :: slot

        do slot = 1, n_ib_nbrs
            if (ib_nbrs(slot) == rank) return
        end do
        call s_mpi_abort('IB owner is not a neighborhood rank')

    end function f_ib_nbr_slot

    !> Adds ib_heat(1:5) (the host copy: this rank's partial sums per slot) over the ranks that hold each body, through its owner
    !! (patch_ib%owner_rank), and, with ib_mass_relax, steps the owners' banks. Phase 1: each rank sends its nonzero partials of
    !! bodies it does not own to the owner, which adds them in ascending source-rank order. Phase 2: each owner returns the totals
    !! and the bank's release (row 7) to the ranks that sent, the only ranks with cells of the body. Both are point-to-point with
    !! neighborhood ranks and carry only bodies a rank shares, so the volume follows the bodies on rank boundaries, not the body
    !! count.
    impure subroutine s_ib_exchange_heat(dt_step)

        real(wp), intent(in) :: dt_step

#ifdef MFC_MPI
        integer, parameter :: nf = 6  !< gbl_patch_id, ib_heat(1:5)
        integer            :: i, j, s, r, slot, pos, nvals, ierr

        if (n_ib_nbrs > 0) then
            do s = 1, n_ib_nbrs
                call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 612, MPI_COMM_WORLD, ib_reqs(s), ierr)
            end do
            ib_send_counts = 0
            do i = 1, num_ibs
                if (patch_ib(i)%owner_rank == proc_rank) cycle
                if (all(ib_heat(1:5,i) == 0._wp)) cycle
                slot = f_ib_nbr_slot(patch_ib(i)%owner_rank)
                pos = nf*ib_send_counts(slot)
                ib_send_bufs(pos + 1, slot) = real(patch_ib(i)%gbl_patch_id, wp)
                ib_send_bufs(pos + 2:pos + nf,slot) = ib_heat(1:5,i)
                ib_send_counts(slot) = ib_send_counts(slot) + 1
            end do
            do s = 1, n_ib_nbrs
                call MPI_ISEND(ib_send_bufs(:,s), nf*ib_send_counts(s), mpi_p, ib_nbrs(s), 612, MPI_COMM_WORLD, &
                               & ib_reqs(n_ib_nbrs + s), ierr)
            end do
            call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

            ! own partial first, then each neighbor in ascending rank order, whatever order the messages arrived in
            ib_reply_counts = 0
            do s = 1, n_ib_nbrs
                call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
                do r = 0, nvals/nf - 1
                    pos = nf*r
                    call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                    @:ASSERT(j > 0, 'IB partial for an IB this rank does not hold')
                    @:ASSERT(patch_ib(j)%owner_rank == proc_rank, 'IB partial for an IB this rank does not own')
                    ib_heat(1:5,j) = ib_heat(1:5,j) + ib_recv_bufs(pos + 2:pos + nf,s)
                    ib_reply_counts(s) = ib_reply_counts(s) + 1
                    ib_reply_slots(ib_reply_counts(s), s) = j
                end do
            end do
        end if

        ib_heat(7,1:num_ibs) = 0._wp
        if (ib_mass_relax > 0) then
            do i = 1, num_local_ibs
                j = local_ib_patch_ids(i)
                call s_ib_bank_step(ib_heat(6, j), ib_heat(7, j), (ib_heat(3, j) - ib_heat(4, j))*dt_step, ib_mass_relax)
            end do
        end if
        if (n_ib_nbrs == 0) return

        do s = 1, n_ib_nbrs
            ib_reqs(s) = MPI_REQUEST_NULL
            ib_reqs(n_ib_nbrs + s) = MPI_REQUEST_NULL
            if (ib_send_counts(s) > 0) then
                call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 613, MPI_COMM_WORLD, ib_reqs(s), ierr)
            end if
            if (ib_reply_counts(s) == 0) cycle
            do r = 0, ib_reply_counts(s) - 1
                j = ib_reply_slots(r + 1, s)
                pos = ib_rec_len*r
                ib_send_bufs(pos + 1, s) = real(patch_ib(j)%gbl_patch_id, wp)
                ib_send_bufs(pos + 2:pos + 6,s) = ib_heat(1:5,j)
                ib_send_bufs(pos + 7, s) = ib_heat(7, j)
            end do
            call MPI_ISEND(ib_send_bufs(:,s), ib_rec_len*ib_reply_counts(s), mpi_p, ib_nbrs(s), 613, MPI_COMM_WORLD, &
                           & ib_reqs(n_ib_nbrs + s), ierr)
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

        do s = 1, n_ib_nbrs
            if (ib_send_counts(s) == 0) cycle
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            do r = 0, nvals/ib_rec_len - 1
                pos = ib_rec_len*r
                call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                @:ASSERT(j > 0, 'IB total for an IB this rank does not hold')
                ib_heat(1:5,j) = ib_recv_bufs(pos + 2:pos + 6,s)
                ib_heat(7, j) = ib_recv_bufs(pos + 7, s)
            end do
        end do
#endif

    end subroutine s_ib_exchange_heat

    !> At a handoff, move the ib_mass_relax banks of the bodies that changed owner from the old owner to the new one. The new owner
    !! knows which neighbors to expect from (gain_from); the old owner queued its records in s_handoff_ib_ownership.
    impure subroutine s_ib_migrate_banks(n_gain, gain_gid, gain_from, num_entries)

        integer, intent(in) :: n_gain, gain_gid(:), gain_from(:), num_entries

#ifdef MFC_MPI
        integer :: s, r, j, nvals, ierr, n_exp(max(n_ib_nbrs, 1))

        if (n_ib_nbrs == 0) return
        n_exp = 0
        do r = 1, n_gain
            s = f_ib_nbr_slot(gain_from(r))
            n_exp(s) = n_exp(s) + 1
        end do
        do s = 1, n_ib_nbrs
            ib_reqs(s) = MPI_REQUEST_NULL
            ib_reqs(n_ib_nbrs + s) = MPI_REQUEST_NULL
            if (n_exp(s) > 0) then
                call MPI_IRECV(ib_recv_bufs(:,s), 2*n_exp(s), mpi_p, ib_nbrs(s), 620, MPI_COMM_WORLD, ib_reqs(s), ierr)
            end if
            if (ib_send_counts(s) > 0) then
                call MPI_ISEND(ib_send_bufs(:,s), 2*ib_send_counts(s), mpi_p, ib_nbrs(s), 620, MPI_COMM_WORLD, &
                               & ib_reqs(n_ib_nbrs + s), ierr)
            end if
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)
        do s = 1, n_ib_nbrs
            if (n_exp(s) == 0) cycle
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            @:ASSERT(nvals == 2*n_exp(s), 'IB bank handoff: record count mismatch')
            do r = 0, n_exp(s) - 1
                call s_get_neighborhood_idx(nint(ib_recv_bufs(2*r + 1, s)), j, num_entries)
                @:ASSERT(j > 0, 'IB bank handoff for an IB this rank does not hold')
                patch_ib(j)%mass_bank = ib_recv_bufs(2*r + 2, s)
            end do
        end do
#endif

    end subroutine s_ib_migrate_banks

    !> Reduces forces and torques so every rank holding an IB ends with bit-identical values. Each IB has one owner, the rank whose
    !! subdomain holds its centroid (patch_ib%owner_rank, stamped at handoff). Phase 1: each rank sends its nonzero partials to the
    !! owner, which adds them to its own partial in ascending source-rank order once all have arrived. Phase 2: each owner sends its
    !! totals to every neighborhood rank, which overwrite their partials. One nonblocking message per neighbor per phase.
    subroutine s_communicate_ib_forces(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques

#ifdef MFC_MPI
        integer :: i, j, s, r, slot, pos, nvals, ierr

        if (num_procs == 1 .or. n_ib_nbrs == 0) return

        ! Phase 1: send each partial to the owner of its IB
        do s = 1, n_ib_nbrs
            call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 610, MPI_COMM_WORLD, ib_reqs(s), ierr)
        end do

        ib_send_counts = 0
        do i = 1, num_ibs
            if (patch_ib(i)%owner_rank == proc_rank) cycle
            if (all(forces(i,:) == 0._wp) .and. all(torques(i,:) == 0._wp)) cycle
            slot = f_ib_nbr_slot(patch_ib(i)%owner_rank)
            pos = ib_rec_len*ib_send_counts(slot)
            ib_send_bufs(pos + 1, slot) = real(patch_ib(i)%gbl_patch_id, wp)
            ib_send_bufs(pos + 2:pos + 4,slot) = forces(i,1:3)
            ib_send_bufs(pos + 5:pos + 7,slot) = torques(i,1:3)
            ib_send_counts(slot) = ib_send_counts(slot) + 1
        end do

        ! every neighbor gets a message, possibly empty, so every posted receive completes
        do s = 1, n_ib_nbrs
            call MPI_ISEND(ib_send_bufs(:,s), ib_rec_len*ib_send_counts(s), mpi_p, ib_nbrs(s), 610, MPI_COMM_WORLD, &
                           & ib_reqs(n_ib_nbrs + s), ierr)
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

        ! owner sum: own partial first, then each neighbor in ascending rank order, whatever order the messages arrived in
        do s = 1, n_ib_nbrs
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            do r = 0, nvals/ib_rec_len - 1
                pos = ib_rec_len*r
                call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                @:ASSERT(j > 0, 'IB force contribution for an IB this rank does not hold')
                @:ASSERT(patch_ib(j)%owner_rank == proc_rank, 'IB force contribution for an IB this rank does not own')
                forces(j,1:3) = forces(j,1:3) + ib_recv_bufs(pos + 2:pos + 4,s)
                torques(j,1:3) = torques(j,1:3) + ib_recv_bufs(pos + 5:pos + 7,s)
            end do
        end do

        ! Phase 2: each owner sends its totals to every neighbor
        do s = 1, n_ib_nbrs
            call MPI_IRECV(ib_recv_bufs(:,s), size(ib_recv_bufs, 1), mpi_p, ib_nbrs(s), 611, MPI_COMM_WORLD, ib_reqs(s), ierr)
        end do

        do i = 1, num_local_ibs
            j = local_ib_patch_ids(i)
            pos = ib_rec_len*(i - 1)
            ib_send_bufs(pos + 1, 1) = real(patch_ib(j)%gbl_patch_id, wp)
            ib_send_bufs(pos + 2:pos + 4,1) = forces(j,1:3)
            ib_send_bufs(pos + 5:pos + 7,1) = torques(j,1:3)
        end do

        do s = 1, n_ib_nbrs
            call MPI_ISEND(ib_send_bufs(:,1), ib_rec_len*num_local_ibs, mpi_p, ib_nbrs(s), 611, MPI_COMM_WORLD, &
                           & ib_reqs(n_ib_nbrs + s), ierr)
        end do
        call MPI_WAITALL(2*n_ib_nbrs, ib_reqs, ib_stats, ierr)

        ! holders take the owner's bits
        do s = 1, n_ib_nbrs
            call MPI_GET_COUNT(ib_stats(:,s), mpi_p, nvals, ierr)
            do r = 0, nvals/ib_rec_len - 1
                pos = ib_rec_len*r
                call s_get_neighborhood_idx(nint(ib_recv_bufs(pos + 1, s)), j)
                if (j <= 0) cycle
                @:ASSERT(patch_ib(j)%owner_rank == ib_nbrs(s), 'IB force total from a rank that does not own the IB')
                forces(j,1:3) = ib_recv_bufs(pos + 2:pos + 4,s)
                torques(j,1:3) = ib_recv_bufs(pos + 5:pos + 7,s)
            end do
        end do
#endif

    end subroutine s_communicate_ib_forces

    !> @brief Swaps ownership of IBs and passes ownership of IBs to neighbor processors
    subroutine s_handoff_ib_ownership()

        integer                     :: i, j, k, output_idx, local_output_idx
        integer                     :: num_ibs_old, num_ibs_pre
        integer                     :: recv_count
        integer                     :: pack_pos, unpack_pos, buf_size, patch_bytes
        integer                     :: send_neighbor, recv_neighbor, ierr
        integer                     :: dx, dy, dz, tag, nbr_idx, nreqs
        real(wp), dimension(3)      :: centroid
        type(ib_patch_parameters)   :: tmp_patch
        integer, dimension(num_ibs) :: old_to_new  ! old patch_ib slot -> slot after compaction
        ! (2R+1)^num_dims - 1 neighbors in the radius-R neighborhood; each gets its own recv buffer
        integer :: max_nbrs
        integer :: n_gain, slot

#ifdef MFC_MPI
        call nvtxStartRange("HANDOFF-IB-OWNERSHIP")

        if (num_procs > 1) then
            patch_bytes = storage_size(tmp_patch)/8
            ! a rank can hand off at most every global patch, so size by num_gbl_ibs, not num_local_ibs_max
            buf_size = storage_size(0)/8 + patch_bytes*max(1, min(num_local_ibs_max, num_gbl_ibs))
            max_nbrs = (2*ib_neighborhood_radius + 1)**num_dims - 1
            if (.not. allocated(ho_send_buf)) then
                allocate (ho_send_buf(buf_size), ho_recv_bufs(buf_size, max_nbrs), ho_reqs(2*max_nbrs), ho_recv_nbrs(max_nbrs))
                allocate (ho_gain_gid(num_local_ibs_max), ho_gain_from(num_local_ibs_max), ho_gain_slot(num_local_ibs_max))
            end if

            ! Sync GPU-updated fields (angles, angular_vel, centroids) to host before
            ! compaction and MPI packing, which read from host memory.
            $:GPU_UPDATE(host='[patch_ib(1:num_ibs)]')

            ! delete any particles that no longer need to be tracked and coalesce the array
            num_ibs_old = num_ibs
            output_idx = 0
            local_output_idx = 0
            n_gain = 0
            do i = 1, num_ibs
                old_to_new(i) = -1
                centroid = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, 0._wp]
                if (num_dims == 3) centroid(3) = patch_ib(i)%z_centroid

                ! delete if not in neighborhood
                if (f_neighborhood_ranks_own_location(centroid)) then
                    output_idx = output_idx + 1
                    old_to_new(i) = output_idx
                    if (i /= output_idx) then
                        patch_ib(output_idx) = patch_ib(i)
                    end if

                    ! check if in local domain
                    if (f_local_rank_owns_location(centroid, glb_bounds)) then
                        local_output_idx = local_output_idx + 1
                        @:PROHIBIT(local_output_idx > num_local_ibs_max, &
                                   & "Too many IBs on a single processor rank. Modify case file or increase num_local_ibs_max to resolve.")
                        local_ib_patch_ids(local_output_idx) = output_idx
                        if (patch_ib(output_idx)%owner_rank /= proc_rank) then
                            ! gained from a neighbor, which sends the ib_mass_relax bank (s_ib_migrate_banks)
                            n_gain = n_gain + 1
                            ho_gain_gid(n_gain) = patch_ib(output_idx)%gbl_patch_id
                            ho_gain_from(n_gain) = patch_ib(output_idx)%owner_rank
                            ho_gain_slot(n_gain) = output_idx
                        end if
                        patch_ib(output_idx)%owner_rank = proc_rank
                    end if
                end if
            end do
            num_ibs = output_idx
            num_local_ibs = local_output_idx
            ! num_ibs shrinks here, so refresh it with patch_ib: s_get_neighborhood_idx bounds its search by it on the device
            $:GPU_UPDATE(device='[patch_ib(1:num_ibs), num_ibs]')
            call s_compact_ib_lookup(old_to_new, num_ibs_old)

            ! Broadcast newly-owned patches (the ones gained above, in local order) to all neighborhood neighbors
            pack_pos = 0
            call MPI_PACK(n_gain, 1, MPI_INTEGER, ho_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, ierr)
            do i = 1, n_gain
                call MPI_PACK(patch_ib(ho_gain_slot(i)), patch_bytes, MPI_BYTE, ho_send_buf, buf_size, pack_pos, MPI_COMM_WORLD, &
                              & ierr)
            end do

            ! Post all receives first, then sends
            nreqs = 0
            nbr_idx = 0
            do dz = merge(-ib_neighborhood_radius, 0, num_dims == 3), merge(ib_neighborhood_radius, 0, num_dims == 3)
                do dy = -ib_neighborhood_radius, ib_neighborhood_radius
                    do dx = -ib_neighborhood_radius, ib_neighborhood_radius
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        nbr_idx = nbr_idx + 1
                        ! one tag per offset in the (2R+1)^3 neighbourhood: a radix-3 encoding collides once R > 1
                        tag = 200 + ((dx + ib_neighborhood_radius)*(2*ib_neighborhood_radius + 1) + (dy + ib_neighborhood_radius)) &
                                     & *(2*ib_neighborhood_radius + 1) + (dz + ib_neighborhood_radius)
                        recv_neighbor = ib_neighbor_ranks(-dx, -dy, -dz)
                        ho_recv_nbrs(nbr_idx) = MPI_PROC_NULL
                        if (recv_neighbor < 0) cycle
                        ho_recv_nbrs(nbr_idx) = recv_neighbor
                        nreqs = nreqs + 1
                        call MPI_IRECV(ho_recv_bufs(:,nbr_idx), buf_size, MPI_PACKED, recv_neighbor, tag, MPI_COMM_WORLD, &
                                       & ho_reqs(nreqs), ierr)
                    end do
                end do
            end do

            do dz = merge(-ib_neighborhood_radius, 0, num_dims == 3), merge(ib_neighborhood_radius, 0, num_dims == 3)
                do dy = -ib_neighborhood_radius, ib_neighborhood_radius
                    do dx = -ib_neighborhood_radius, ib_neighborhood_radius
                        if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                        ! one tag per offset in the (2R+1)^3 neighbourhood: a radix-3 encoding collides once R > 1
                        tag = 200 + ((dx + ib_neighborhood_radius)*(2*ib_neighborhood_radius + 1) + (dy + ib_neighborhood_radius)) &
                                     & *(2*ib_neighborhood_radius + 1) + (dz + ib_neighborhood_radius)
                        send_neighbor = ib_neighbor_ranks(dx, dy, dz)
                        if (send_neighbor < 0) cycle
                        nreqs = nreqs + 1
                        call MPI_ISEND(ho_send_buf, pack_pos, MPI_PACKED, send_neighbor, tag, MPI_COMM_WORLD, ho_reqs(nreqs), ierr)
                    end do
                end do
            end do

            call MPI_WAITALL(nreqs, ho_reqs, MPI_STATUSES_IGNORE, ierr)

            ! Unpack all received buffers
            num_ibs_pre = num_ibs
            ib_send_counts = 0
            do nbr_idx = 1, max_nbrs
                if (ho_recv_nbrs(nbr_idx) == MPI_PROC_NULL) cycle
                unpack_pos = 0
                call MPI_UNPACK(ho_recv_bufs(:,nbr_idx), buf_size, unpack_pos, recv_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                do i = 1, recv_count
                    call MPI_UNPACK(ho_recv_bufs(:,nbr_idx), buf_size, unpack_pos, tmp_patch, patch_bytes, MPI_BYTE, &
                                    & MPI_COMM_WORLD, ierr)
                    call s_get_neighborhood_idx(tmp_patch%gbl_patch_id, j, num_ibs_pre)
                    ! A rank filling several table slots (periodicity, few ranks) sends each arrival once per slot, and the map
                    ! holds only the num_ibs_pre older entries: look among this call's arrivals too
                    if (j < 0) then
                        do k = num_ibs_pre + 1, num_ibs
                            if (patch_ib(k)%gbl_patch_id == tmp_patch%gbl_patch_id) j = k
                        end do
                    end if
                    if (j < 0) then
                        num_ibs = num_ibs + 1
                        @:ASSERT(num_ibs <= size(patch_ib), 'patch_ib overflow in neighborhood handoff')
                        @:PROHIBIT(num_ibs > num_local_ibs_max, &
                                   & "Owned plus neighboring IBs on a rank exceed num_local_ibs_max, which sizes ib_box and ib_face_ft. Use more ranks or increase num_local_ibs_max.")
                        patch_ib(num_ibs) = tmp_patch
                    else
                        ! an IB we already hold changed owner: record the new owner, keep our own copy of its state; the old
                        ! owner queues its bank for the new one
                        if (ib_mass_relax > 0 .and. patch_ib(j)%owner_rank == proc_rank) then
                            slot = f_ib_nbr_slot(tmp_patch%owner_rank)
                            ib_send_counts(slot) = ib_send_counts(slot) + 1
                            ib_send_bufs(2*ib_send_counts(slot) - 1, slot) = real(tmp_patch%gbl_patch_id, wp)
                            ib_send_bufs(2*ib_send_counts(slot), slot) = patch_ib(j)%mass_bank
                        end if
                        patch_ib(j)%owner_rank = tmp_patch%owner_rank
                    end if
                end do
            end do

            @:PROHIBIT(num_ibs > num_gbl_ibs, "IB handoff duplicated a patch: num_ibs exceeds num_gbl_ibs.")
            if (ib_mass_relax > 0) call s_ib_migrate_banks(n_gain, ho_gain_gid, ho_gain_from, num_ibs_pre)
            $:GPU_UPDATE(device='[patch_ib(1:num_ibs), num_ibs]')
            call s_merge_ib_lookup(num_ibs_pre)
            @:PROHIBIT(ib_slot_cap > 0 .and. num_ibs > ib_slot_cap, &
                       & "A rank holds more IBs than its ib_heat slots; clustering of IBs exceeds the 2x margin set at start-up")
        end if

        call nvtxEndRange()
#endif

    end subroutine s_handoff_ib_ownership

    !> Abort if any namelist immersed boundary marks no cell anywhere. The owning rank is chosen by the patch centroid, so an STL
    !! placed with model_translate while its centroid is left elsewhere is handed to ranks that do not hold it, and it vanishes
    !! silently. Cloud particles are skipped: a particle's centroid is its body, so it cannot vanish this way.
    impure subroutine s_check_every_patch_marked()

        integer(kind=8), allocatable :: cnt_loc(:), cnt_glb(:)
        integer                      :: gid, i, j, k, num_checked

        num_checked = num_gbl_ibs - sum(particle_cloud(1:num_particle_clouds)%num_particles)
        if (num_checked <= 0) return
        allocate (cnt_loc(num_checked), cnt_glb(num_checked))
        cnt_loc = 0_8
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (ib_markers%sf(i, j, k) /= 0) then
                        call s_decode_patch_periodicity(ib_markers%sf(i, j, k), gid)
                        if (gid <= num_checked) cnt_loc(gid) = cnt_loc(gid) + 1_8
                    end if
                end do
            end do
        end do
        call s_mpi_allreduce_integer_sum_vec(cnt_loc, cnt_glb)
        @:PROHIBIT(any(cnt_glb == 0_8), &
                   & "An immersed boundary marked no cell anywhere: its centroid decides which rank owns it, so a body placed " &
                   & // "elsewhere with model_translate is handed to a rank that does not hold it. Set patch_ib%x/y/z_centroid " &
                   & // "to where the body actually is.")
        deallocate (cnt_loc, cnt_glb)

    end subroutine s_check_every_patch_marked

    !> Compute the deepest-level marker-field bounds into mkr_lo/mkr_hi: the coarse block (m/n/p with ghosts) or the deepest fine
    !! block a rank can own, whichever is larger (a level-l block has amr_ref_ratio**l * base_ext - 1 interior cells per active dim,
    !! base_ext = the user footprint in coarse cells). Called from s_initialize_ibm_module and s_ibm_alloc_fine.
    impure subroutine s_ibm_marker_bounds()

        mkr_lo(1) = -buff_size
        mkr_hi(1) = max(m, amr_ref_ratio**amr_max_level*(amr_block_end(1) - amr_block_beg(1) + 1) - 1) + buff_size
        mkr_lo(2) = -buff_size
        if (n_glb > 0) then
            mkr_hi(2) = max(n, amr_ref_ratio**amr_max_level*(amr_block_end(2) - amr_block_beg(2) + 1) - 1) + buff_size
        else
            mkr_hi(2) = n + buff_size
        end if
        if (p > 0) then
            mkr_lo(3) = -buff_size
            mkr_hi(3) = max(p, amr_ref_ratio**amr_max_level*(amr_block_end(3) - amr_block_beg(3) + 1) - 1) + buff_size
        else
            mkr_lo(3) = 0; mkr_hi(3) = 0
        end if

    end subroutine s_ibm_marker_bounds

    !> Allocate the per-slot fine-IB marker fields, sized to the max buffered fine extents; ghost-point lists start empty. No-op
    !! unless amr .and. ib.
    impure subroutine s_ibm_alloc_fine(nslots, f1_lo, f1_hi, f2_lo, f2_hi, f3_lo, f3_hi)

        integer, intent(in) :: nslots, f1_lo, f1_hi, f2_lo, f2_hi, f3_lo, f3_hi
        integer             :: islot

        call s_ibm_marker_bounds()
        @:PROHIBIT(f1_hi > mkr_hi(1) .or. f2_hi > mkr_hi(2) .or. (p > 0 .and. f3_hi > mkr_hi(3)), &
                   & "AMR fine IB: fine block extent exceeds the deepest-level ib_markers bounds; the copy-based fine-marker swap needs ib_markers sized to enclose the fine block")

        ! ghost-point capacity: the deepest fine block's buffered cell count, from the widened mkr bounds
        fine_gps_cap = int(mkr_hi(1) - mkr_lo(1) + 1, 8)*int(mkr_hi(2) - mkr_lo(2) + 1, 8)*int(max(mkr_hi(3) - mkr_lo(3) + 1, 1), 8)

        ! Extra slot (nslots+1) parks the coarse markers during a fine swap (see ib_coarse_slot). markers%sf is host-only park
        ! storage, allocated on a slot's first swap (s_ibm_fine_markers_ensure); fine ghost-point lists park on-device in gp_park.
        ib_coarse_slot = nslots + 1
        allocate (ib_fine(1:ib_coarse_slot))
        do islot = 1, ib_coarse_slot
            ib_fine(islot)%num_gps = 0
            ib_fine(islot)%off = 0_8
            ib_fine(islot)%cap = 0_8
        end do

    end subroutine s_ibm_alloc_fine

    !> Allocate slot islot's host marker copy (zeroed) on first use: a copy per possible slot would be amr_max_blocks + 1 full
    !! deepest-level fields per rank, nearly all unused (#2219).
    impure subroutine s_ibm_fine_markers_ensure(islot)

        integer, intent(in) :: islot

        if (associated(ib_fine(islot)%markers%sf)) return
        allocate (ib_fine(islot)%markers%sf(mkr_lo(1):mkr_hi(1),mkr_lo(2):mkr_hi(2),mkr_lo(3):mkr_hi(3)))
        ib_fine(islot)%markers%sf = 0

    end subroutine s_ibm_fine_markers_ensure

    !> Note the fine IB state's memory for the end-of-run AMR table.
    impure subroutine s_ibm_note_fine_memory()

        integer :: islot

        if (.not. allocated(ib_fine)) return
        call s_mem_note(MEM_IB, storage_size(gp_park, kind=8)/8_8*size(gp_park, kind=8))
        do islot = 1, ib_coarse_slot
            if (.not. associated(ib_fine(islot)%markers%sf)) cycle
            call s_mem_note(MEM_IBHOST, storage_size(ib_fine(islot)%markers%sf, kind=8)/8_8*size(ib_fine(islot)%markers%sf, kind=8))
        end do

    end subroutine s_ibm_note_fine_memory

    !> Forget every slot's parked fine list (before the fine lists are rebuilt for the current block set).
    impure subroutine s_ibm_reset_fine_park()

        if (.not. allocated(ib_fine)) return
        ib_fine(:)%num_gps = 0
        ib_fine(:)%off = 0_8
        ib_fine(:)%cap = 0_8
        park_top = park_coarse

    end subroutine s_ibm_reset_fine_park

    !> Ghost points stored for fine slot islot: zero for a block that holds no immersed body.
    impure integer function f_ibm_fine_gps(islot) result(n)

        integer, intent(in) :: islot

        n = ib_fine(islot)%num_gps

    end function f_ibm_fine_gps

    !> Add to cnt(1:num_gbl_ibs) the ghost points each body has in the active list (the fine slot's when swapped in); counted on the
    !! device, since ghost_points never crosses to the host.
    impure subroutine s_ibm_tally_gps_by_body(cnt)

        integer(kind=8), intent(inout) :: cnt(:)
        integer                        :: dcnt(size(cnt))
        integer                        :: a, gid, ng

        if (num_gps == 0) return
        dcnt = 0
        ng = num_gps
        $:GPU_PARALLEL_LOOP(private='[a, gid]', copy='[dcnt]', firstprivate='[ng]')
        do a = 1, ng
            gid = patch_ib(ghost_points(a)%ib_patch_id)%gbl_patch_id
            $:GPU_ATOMIC(atomic='update')
            dcnt(gid) = dcnt(gid) + 1
        end do
        $:END_GPU_PARALLEL_LOOP()
        cnt = cnt + int(dcnt, 8)

    end subroutine s_ibm_tally_gps_by_body

    !> Abort if any one body has ghost points on more than one rank: the fine image-point stencil is exact only when a body's whole
    !! fine state lives on one rank. Counts are per global body, so bodies on different ranks are fine.
    impure subroutine s_ibm_check_bodies_on_one_rank(cnt)

        integer(kind=8), intent(in)  :: cnt(:)
        integer(kind=8), allocatable :: flag(:), tot(:)
        integer                      :: ierr, b

        allocate (flag(size(cnt)), tot(size(cnt)))
        flag = merge(1_8, 0_8, cnt > 0_8)
        tot = flag
#ifdef MFC_MPI
        call MPI_ALLREDUCE(flag, tot, size(flag), MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        do b = 1, size(tot)
            if (tot(b) > 1_8) then
                if (proc_rank == 0) write (0, '(A,I0,A,I0,A)') ' [amr] immersed body ', b, ' has ghost points on ', tot(b), ' ranks'
                call s_mpi_abort('amr with ib: an immersed body straddles a rank boundary, where the fine-IB ' &
                                 & // 'image-point stencil is not yet decomposition-exact; keep each body within a single ' &
                                 & // 'rank subdomain (use fewer ranks or reposition it).')
            end if
        end do

    end subroutine s_ibm_check_bodies_on_one_rank

    !> Load ib_markers_fine for the next fine RHS pass: member ibm's fine markers (interior 0:mext(:, ibm)) land at offset (ibm -
    !! 1)*w along dimension sd; the rest of the installed frame and its first ghost layer read 0, so padding and gaps zero nothing.
    !! Per-block callers pass nb = 1, w = 0. Grid globals must already be swapped.
    impure subroutine s_ibm_load_fine_markers(nb, slots, mext, sd, w)

        integer, intent(in) :: nb, sd, w
        integer, intent(in) :: slots(nb), mext(3, nb)
        integer             :: ibm, o(3), e2, e3

        e2 = min(n, 1); e3 = min(p, 1)
        ib_markers_fine%sf(-1:m + 1,-e2:n + e2,-e3:p + e3) = 0
        do ibm = 1, nb
            if (.not. associated(ib_fine(slots(ibm))%markers%sf)) cycle  ! never swapped in: no body markers
            o = 0; o(sd) = (ibm - 1)*w
            ib_markers_fine%sf(o(1):o(1) + mext(1, ibm),o(2):o(2) + mext(2, ibm),o(3):o(3) + mext(3, &
                               & ibm)) = ib_fine(slots(ibm))%markers%sf(0:mext(1, ibm),0:mext(2, ibm),0:mext(3, ibm))
        end do
        $:GPU_UPDATE(device='[ib_markers_fine%sf(:, :, -e3:p + e3)]')
        fine_lay_nb = nb; fine_lay_sd = sd; fine_lay_w = w
        fine_lay_slots = slots

    end subroutine s_ibm_load_fine_markers

    !> Swap the module IB globals (ib_markers/ghost_points/num_gps) to fine slot islot's state; the coarse state parks in the save
    !! slot. Pair with s_ibm_restore_from_fine; grid globals must already be swapped to the fine block.
    impure subroutine s_ibm_swap_to_fine(islot, gps_on_device)

        integer, intent(in) :: islot
        logical, intent(in) :: gps_on_device  !< fine ghost points already present on device (per-stage correct path)
        integer             :: n_c, n_f, a
        integer(kind=8)     :: off

        ! ib_markers is device-resident and never pointer-swapped (corrupts the Cray present table): park the coarse markers as a
        ! host copy; on the correct/moving path copy this slot's fine markers in and push. On the setup path s_ibm_setup_fine
        ! rebuilds them in place.

        call s_ibm_fine_markers_ensure(ib_coarse_slot)
        call s_ibm_fine_markers_ensure(islot)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        ib_fine(ib_coarse_slot)%markers%sf = ib_markers%sf
        if (gps_on_device) then
            ib_markers%sf = ib_fine(islot)%markers%sf
            $:GPU_UPDATE(device='[ib_markers%sf]')
        end if

        ! ghost_points and gp_park are device-resident and never reallocated after setup. Park the coarse list into gp_park's coarse
        ! segment, then (correct/moving path) pull this slot's fine list in, via on-device kernels; a whole-array GPU_UPDATE of
        ! ghost_points aborts amdflang (per-element custom mapper, ROCm HSA OUT_OF_RESOURCES). On the setup path s_ibm_setup_fine
        ! fills ghost_points in place.
        n_c = num_gps
        @:PROHIBIT(int(n_c, 8) > park_coarse, "AMR fine IB: coarse ghost-point count exceeds the gp_park capacity")
        $:GPU_PARALLEL_LOOP(private='[a]', firstprivate='[n_c]')
        do a = 1, n_c
            gp_park(a) = ghost_points(a)
        end do
        $:END_GPU_PARALLEL_LOOP()
        num_gps_save = num_gps
        num_gps = ib_fine(islot)%num_gps
        if (gps_on_device) then
            n_f = num_gps
            off = ib_fine(islot)%off
            @:PROHIBIT(int(n_f, 8) > ib_fine(islot)%cap, "AMR fine IB: a slot with ghost points has no gp_park segment")
            $:GPU_PARALLEL_LOOP(private='[a]', firstprivate='[n_f, off]')
            do a = 1, n_f
                ghost_points(a) = gp_park(off + a)
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_ibm_rebuild_corrected_gps(n_f)
        end if
        $:GPU_UPDATE(device='[num_gps]')

    end subroutine s_ibm_swap_to_fine

    !> Restore the coarse IB globals saved by s_ibm_swap_to_fine, parking the (possibly updated) fine state back in slot islot.
    impure subroutine s_ibm_restore_from_fine(islot)

        integer, intent(in) :: islot
        integer             :: n_c, n_f, a
        integer(kind=8)     :: off

        ! Mirror s_ibm_swap_to_fine: save this slot's fine markers to its host store, copy the parked coarse markers back to the
        ! device ib_markers, and restore the ghost-point lists between gp_park and ghost_points on-device.

        call s_ibm_fine_markers_ensure(islot)
        $:GPU_UPDATE(host='[ib_markers%sf]')
        ib_fine(islot)%markers%sf = ib_markers%sf
        ib_markers%sf = ib_fine(ib_coarse_slot)%markers%sf
        $:GPU_UPDATE(device='[ib_markers%sf]')

        n_f = num_gps
        if (int(n_f, 8) > ib_fine(islot)%cap) call s_ibm_reserve_park(islot, n_f)
        off = ib_fine(islot)%off
        $:GPU_PARALLEL_LOOP(private='[a]', firstprivate='[n_f, off]')
        do a = 1, n_f
            gp_park(off + a) = ghost_points(a)
        end do
        $:END_GPU_PARALLEL_LOOP()
        ib_fine(islot)%num_gps = num_gps
        num_gps = num_gps_save
        n_c = num_gps
        $:GPU_PARALLEL_LOOP(private='[a]', firstprivate='[n_c]')
        do a = 1, n_c
            ghost_points(a) = gp_park(a)
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_ibm_rebuild_corrected_gps(n_c)
        $:GPU_UPDATE(device='[num_gps]')

    end subroutine s_ibm_restore_from_fine

    !> Give slot islot a gp_park segment for n_f ghost points at the pool top: exact for fixed bodies (the list keeps its length
    !! until s_ibm_reset_fine_park), full-length for moving bodies (so a slot never needs a second one, and a vacated one is
    !! reused).
    impure subroutine s_ibm_reserve_park(islot, n_f)

        integer, intent(in) :: islot, n_f
        integer             :: j

        if (moving_immersed_boundary_flag) then
            ! take over the segment of a slot the body has left
            do j = 1, ib_coarse_slot - 1
                if (j == islot .or. ib_fine(j)%cap == 0_8 .or. ib_fine(j)%num_gps > 0) cycle
                ib_fine(islot)%off = ib_fine(j)%off
                ib_fine(islot)%cap = ib_fine(j)%cap
                ib_fine(j)%cap = 0_8
                return
            end do
        end if
        ib_fine(islot)%off = park_top
        ib_fine(islot)%cap = merge(park_coarse, int(n_f, 8), moving_immersed_boundary_flag)
        park_top = park_top + ib_fine(islot)%cap
        @:PROHIBIT(park_top > size(gp_park, kind=8), &
                   & "AMR fine IB: the fine blocks hold more ghost points than gp_park was sized for (see s_ibm_setup)")

    end subroutine s_ibm_reserve_park

    !> Rebuild corrected_gps from the ghost-point list now in ghost_points: s_compute_interpolation_coeffs only sets entries, so
    !! clear the whole field and re-mark the ghost points on-device.
    impure subroutine s_ibm_rebuild_corrected_gps(n_g)

        integer, intent(in) :: n_g
        integer             :: i, j, k, a, lo1, hi1, lo2, hi2, lo3, hi3

        lo1 = lbound(corrected_gps%sf, 1); hi1 = ubound(corrected_gps%sf, 1)
        lo2 = lbound(corrected_gps%sf, 2); hi2 = ubound(corrected_gps%sf, 2)
        lo3 = lbound(corrected_gps%sf, 3); hi3 = ubound(corrected_gps%sf, 3)
        $:GPU_PARALLEL_LOOP(private='[i, j, k]', firstprivate='[lo1, hi1, lo2, hi2, lo3, hi3]', collapse=3)
        do k = lo3, hi3; do j = lo2, hi2; do i = lo1, hi1
            corrected_gps%sf(i, j, k) = 0
        end do; end do; end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[a]', firstprivate='[n_g]')
        do a = 1, n_g
            corrected_gps%sf(ghost_points(a)%loc(1), ghost_points(a)%loc(2), ghost_points(a)%loc(3)) = merge(1, 2, &
                             & ghost_points(a)%interp_valid)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ibm_rebuild_corrected_gps

    !> Compute the fine-grid IB state (markers, ghost points, levelset, image points, interpolation coeffs) for the current block,
    !! mirroring the static-body portion of s_ibm_setup. Grid globals must be swapped to the fine block and IB globals to this slot
    !! (s_ibm_swap_to_fine) first.
    impure subroutine s_ibm_setup_fine()

        integer :: n_g  !< num_gps copied to a local before the rebuild call below; see the comment there

        ib_markers%sf = 0
        $:GPU_UPDATE(device='[ib_markers%sf]')
        call s_apply_ib_patches(ib_markers)
        $:GPU_UPDATE(host='[ib_markers%sf]')

        call s_find_num_ghost_points(num_gps)
        $:GPU_UPDATE(device='[num_gps]')
        ! ghost_points is allocated once in s_ibm_setup (sized to the fine cell-count cap) and filled in place; never reallocated
        @:PROHIBIT(int(num_gps, 8) > size(ghost_points, kind=8), &
                   & "AMR fine IB: ghost-point count exceeds the ghost_points capacity sized at s_ibm_setup")

        call s_find_ghost_points()
        call s_apply_levelset(ghost_points, num_gps)
        call s_compute_image_points()
        call s_compute_interpolation_coeffs()
        ! s_compute_interpolation_coeffs only marks entries; drop the coarse pattern still in corrected_gps.
        ! num_gps goes through a local: NVHPC's two-pass IPO inlines s_ibm_rebuild_corrected_gps here, which would put firstprivate
        ! on a declare-target symbol (NVFORTRAN-S-0155 "Load of NULL symbol").
        n_g = num_gps
        call s_ibm_rebuild_corrected_gps(n_g)

    end subroutine s_ibm_setup_fine

    !> Finalize the IBM module
    impure subroutine s_finalize_ibm_module()

        integer :: i

        deallocate (n_fill_tier)
        @:DEALLOCATE(ib_markers%sf)
        @:DEALLOCATE(corrected_gps%sf)
        if (associated(corrected_gps_r(1)%sf)) then
            @:DEALLOCATE(corrected_gps_r(1)%sf)
        end if
        if (associated(ib_wall_theta%sf)) then
            do i = 1, ibw_a
                @:DEALLOCATE(ib_wall_src(i)%sf)
            end do
            @:DEALLOCATE(ib_wall_theta%sf)
        end if
        @:DEALLOCATE(ib_face_ft)
        if (associated(ib_markers_fine%sf)) then
            @:DEALLOCATE(ib_markers_fine%sf)
        end if
        do i = 1, num_ib_airfoils_max
            if (allocated(ib_airfoil_grids(i)%upper)) then
                @:DEALLOCATE(ib_airfoil_grids(i)%upper)
                @:DEALLOCATE(ib_airfoil_grids(i)%lower)
            end if
        end do
        if (allocated(models)) then
            @:DEALLOCATE(models)
        end if
        if (allocated(ghost_points)) then
            @:DEALLOCATE(ghost_points)
            @:DEALLOCATE(gp_prim_new)
        end if
        if (allocated(gp_tc_rec)) then
            @:DEALLOCATE(gp_tc_rec)
        end if
        if (allocated(gp_surf)) then
            @:DEALLOCATE(gp_surf)
        end if
        ! gp_park is device-mapped (GPU_ENTER_DATA); @:DEALLOCATE unmaps and frees it
        if (allocated(gp_park)) then
            @:DEALLOCATE(gp_park)
        end if
        if (allocated(ib_fine)) then
            do i = 1, size(ib_fine)
                ! markers%sf is host-only parking storage: plain deallocate
                if (associated(ib_fine(i)%markers%sf)) then
                    deallocate (ib_fine(i)%markers%sf)
                end if
            end do
            deallocate (ib_fine)
        end if
        if (allocated(ib_heat)) then
            @:DEALLOCATE(ib_heat, ib_part)
        end if
        if (associated(ib_shell%sf)) then
            @:DEALLOCATE(ib_shell%sf)
        end if
        if (associated(ib_markers_prev%sf)) then
            @:DEALLOCATE(ib_markers_prev%sf)
        end if
        if (allocated(wc_T)) then
            @:DEALLOCATE(wc_T, wc_w, wc_g, wc_ext, wc_q, wc_cut, wc_on, wc_gp)
        end if
        n_wc = 0
        wc_ready = .false.
        call s_finalize_ib_patches()
#ifdef MFC_MPI
        if (allocated(ib_nbrs)) deallocate (ib_nbrs, ib_send_counts, ib_reqs, ib_stats, ib_send_bufs, ib_recv_bufs, &
            & ib_reply_slots, ib_reply_counts)
        if (allocated(ho_send_buf)) deallocate (ho_send_buf, ho_recv_bufs, ho_reqs, ho_recv_nbrs, ho_gain_gid, ho_gain_from, &
            & ho_gain_slot)
#endif

    end subroutine s_finalize_ibm_module

    !> ib_flux_wall: allocate the wall fields, give each body's ghost points with a fluid face neighbor their wall area per face
    !! area, and report the wall area the faces carry. A planar wall of area A has sum|n_i| A of staircase faces, so each face
    !! carries 1/sum|n_i| of its own area.
    impure subroutine s_ib_flux_wall_setup()

        integer  :: i, j, k, l, q, d, sgn, pid
        logical  :: wet
        real(wp) :: area, area_glb, vol

        do q = 1, ibw_a
            @:ALLOCATE(ib_wall_src(q)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(ib_wall_src(q))
        end do
        @:ALLOCATE(ib_wall_theta%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(ib_wall_theta)

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    ib_wall_theta%sf(j, k, l) = 1._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, ibw_a
                        ib_wall_src(q)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, d, sgn, wet, pid]', present='[ghost_points]')
        do i = 1, num_gps
            j = ghost_points(i)%loc(1); k = ghost_points(i)%loc(2); l = ghost_points(i)%loc(3)
            wet = .false.
            do d = 1, num_dims
                do sgn = -1, 1, 2
                    if (ib_markers%sf(j + sgn*merge(1, 0, d == 1), k + sgn*merge(1, 0, d == 2), l + sgn*merge(1, 0, &
                        & d == 3)) == 0) wet = .true.
                end do
            end do
            ! A reacting body, or a no-slip, non-blowing inert one (its zero source carries no mass); all are fixed (case validator)
            pid = ghost_points(i)%ib_patch_id
            if (wet .and. patch_ib(pid)%inj_species == 0 .and. sum(abs(ghost_points(i)%levelset_norm)) > 0._wp &
                & .and. (patch_ib(pid)%surface_reaction == 1 .or. (.not. ghost_points(i)%slip .and. patch_ib(pid)%v_blow <= 0._wp) &
                & )) then
                ib_wall_src(ibw_a)%sf(j, k, l) = sqrt(sum(ghost_points(i)%levelset_norm**2))/sum(abs(ghost_points(i)%levelset_norm))
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_ib_wall_src_buffers()

        area = 0._wp
        $:GPU_PARALLEL_LOOP(private='[j, k, l, vol]', collapse=3, reduction='[[area]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib_markers%sf(j, k, l) == 0) then
                        vol = dx(j)*dy(k)
                        if (num_dims == 3) vol = vol*dz(l)
                        area = area + vol*f_wall_face_sum(ib_wall_src, j, k, l, 0)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mpi_allreduce_sum(area, area_glb)
        if (proc_rank == 0) print '(A, ES14.6)', ' ib_flux_wall: wall area carried by the faces = ', area_glb

    end subroutine s_ib_flux_wall_setup

    !> ib_flux_wall, before each Riemann solve (primitive buffers filled): fill ib_wall_src's buffer cells and set the positivity
    !! factor of every fluid cell in the domain and its first buffer layer (a wall face on a rank seam is solved on both sides), the
    !! largest theta <= 1 for which its wall faces take at most ib_wall_safety of any species in one forward-Euler step of size dt
    !! (each SSP-RK stage is a convex combination of such steps). Conservative: the body books the scaled flux.
    impure subroutine s_ib_flux_wall_prepare(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), parameter                                 :: ib_wall_safety = 0.9_wp
        real(wp), parameter                                 :: ib_wall_share = 1.e-6_wp  !< trace consumption not limited
        real(wp)                                            :: dt_loc, use, th, rhoY, tot
        integer, dimension(3)                               :: lo, hi
        integer                                             :: j, k, l, q

        call s_ib_wall_src_buffers()
        dt_loc = dt
        lo = 0; hi = [m, n, p]
        do q = 1, num_dims
            lo(q) = -1; hi(q) = hi(q) + 1
        end do
        $:GPU_PARALLEL_LOOP(private='[j, k, l, q, use, th, rhoY, tot]', collapse=3, copyin='[dt_loc, lo, hi]')
        do l = lo(3), hi(3)
            do k = lo(2), hi(2)
                do j = lo(1), hi(1)
                    th = 1._wp
                    if (ib_markers%sf(j, k, l) == 0) then
                        if (f_wall_face_sum(ib_wall_src, j, k, l, 0) > 0._wp) then
                            tot = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, num_species
                                tot = tot + f_wall_face_sum(ib_wall_src, j, k, l, q)
                            end do
                            ! A species with a negligible share of the consumption (a radical the cell holds almost none of)
                            ! would otherwise set theta ~ 0 and starve the whole wall
                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, num_species
                                use = f_wall_face_sum(ib_wall_src, j, k, l, q)
                                if (.not. use > ib_wall_share*tot) cycle
                                use = dt_loc*use
                                rhoY = q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l)*q_prim_vf(eqn_idx%species%beg + q - 1)%sf(j, k, l)
                                ! Only a consumed species divides (use > 0 here): 0/0 and min/max with a NaN are processor dependent
                                if (rhoY <= 0._wp) then
                                    th = 0._wp
                                else if (use > ib_wall_safety*rhoY) then
                                    th = min(th, ib_wall_safety*rhoY/use)
                                end if
                            end do
                        end if
                    end if
                    if (.not. (th >= 0._wp)) th = 0._wp
                    ib_wall_theta%sf(j, k, l) = th
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_ib_flux_wall_prepare

    !> ib_flux_wall: sum over the wall faces of fluid cell (j, k, l) of the wall area per cell volume, times species q's consumption
    !! rate max(-mdot_q, 0) (q > 0) or 1 (q = 0)
    function f_wall_face_sum(src, j, k, l, q) result(s)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(ibw_a), intent(in) :: src  !< ib_wall_src, passed so the device routine reads the present copy
        integer, intent(in)                              :: j, k, l, q
        real(wp)                                         :: s

        s = (f_wall_face_term(src, j - 1, k, l, q) + f_wall_face_term(src, j + 1, k, l, q))/dx(j) + (f_wall_face_term(src, j, &
             & k - 1, l, q) + f_wall_face_term(src, j, k + 1, l, q))/dy(k)
        if (num_dims == 3) s = s + (f_wall_face_term(src, j, k, l - 1, q) + f_wall_face_term(src, j, k, l + 1, q))/dz(l)

    end function f_wall_face_sum

    !> ib_flux_wall: wall area per face area times species q's consumption rate (q > 0) or 1 (q = 0) of the body cell (j, k, l)
    function f_wall_face_term(src, j, k, l, q) result(s)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(ibw_a), intent(in) :: src
        integer, intent(in)                              :: j, k, l, q
        real(wp)                                         :: s

        s = src(ibw_a)%sf(j, k, l)
        if (q > 0) s = s*max(-src(q)%sf(j, k, l), 0._wp)

    end function f_wall_face_term

    !> ib_flux_wall: the surface area a reacting body's ghost point (j, k, l) stands for in its surface record, the wall area its
    !! fluid faces carry (a theta times their area, as the Riemann solver takes them); 0 off the faces
    function f_flux_wall_record_area(src, theta, marks, loc) result(w)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(ibw_a), intent(in) :: src
        type(scalar_field), intent(in)                   :: theta
        type(integer_field), intent(in)                  :: marks  !< ib_markers, passed like src
        integer, dimension(3), intent(in)                :: loc
        real(wp)                                         :: w
        real(wp), dimension(3)                           :: h
        integer, dimension(3)                            :: c
        integer                                          :: d, sgn

        w = 0._wp
        h = [dx(loc(1)), dy(loc(2)), 1._wp]
        if (num_dims == 3) h(3) = dz(loc(3))
        do d = 1, num_dims
            do sgn = -1, 1, 2
                c = loc
                c(d) = c(d) + sgn
                if (marks%sf(c(1), c(2), c(3)) == 0) w = w + src(ibw_a)%sf(loc(1), loc(2), loc(3))*theta%sf(c(1), c(2), &
                    & c(3))*h(1)*h(2)*h(3)/h(d)
            end do
        end do

    end function f_flux_wall_record_area

    !> ib_flux_wall: fill ib_wall_src's buffer cells from the neighbouring ranks or the periodic image; they stay zero (no wall) at
    !! other domain boundaries
    impure subroutine s_ib_wall_src_buffers()

        integer, dimension(3) :: lo, hi, c, dst, src_c
        integer               :: d, loc, edge, nd, j, k, l, q

        do d = 1, num_dims
            do loc = -1, 1, 2
                select case (d)
                case (1)
                    edge = merge(bc_x%beg, bc_x%end, loc == -1); nd = m
                case (2)
                    edge = merge(bc_y%beg, bc_y%end, loc == -1); nd = n
                case default
                    edge = merge(bc_z%beg, bc_z%end, loc == -1); nd = p
                end select
                if (edge >= 0) then
                    call s_mpi_sendrecv_variables_buffers(ib_wall_src, d, loc, ibw_a)
                else if (edge == BC_PERIODIC) then
                    lo = [idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg]
                    hi = [idwbuff(1)%end, idwbuff(2)%end, idwbuff(3)%end]
                    lo(d) = 1; hi(d) = buff_size
                    $:GPU_PARALLEL_LOOP(private='[j, k, l, q, c, dst, src_c]', collapse=3, copyin='[lo, hi, d, loc, nd]')
                    do l = lo(3), hi(3)
                        do k = lo(2), hi(2)
                            do j = lo(1), hi(1)
                                c = [j, k, l]
                                dst = c; src_c = c
                                if (loc == -1) then
                                    dst(d) = -c(d); src_c(d) = nd + 1 - c(d)
                                else
                                    dst(d) = nd + c(d); src_c(d) = c(d) - 1
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do q = 1, ibw_a
                                    ib_wall_src(q)%sf(dst(1), dst(2), dst(3)) = ib_wall_src(q)%sf(src_c(1), src_c(2), src_c(3))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end do
        end do

    end subroutine s_ib_wall_src_buffers

    !> Species flux residual at a heterogeneous reacting surface.
    subroutine s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_s, d
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        real(wp), intent(out) :: R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: mw_IP, mw_s, rho_s, sum_BG
        real(wp)              :: Xs_IP(num_species), Xs_s(num_species)
        real(wp)              :: D_s(num_species), B_s(num_species), G_s(num_species)
        integer               :: k

        call get_mixture_molecular_weight(Ys_IP, mw_IP)
        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        do k = 1, num_species
            Xs_IP(k) = Ys_IP(k)*mw_IP/W_species(k)
            Xs_s(k) = Ys_s(k)*mw_s/W_species(k)
        end do

        call get_species_mass_diffusivities_mixavg(pres, T_s, Ys_s, D_s)
        call get_surface_net_production_rates(rho_s, T_s, Ys_s, omega_s)

        mdot_s = 0._wp
        do k = 1, num_species
            mdot_s = mdot_s + W_species(k)*omega_s(k)
        end do

        sum_BG = 0._wp
        do k = 1, num_species
            B_s(k) = rho_s*D_s(k)*W_species(k)/mw_s
            G_s(k) = (Xs_IP(k) - Xs_s(k))/d
            sum_BG = sum_BG + B_s(k)*G_s(k)
        end do

        do k = 1, num_species
            R_species(k) = -B_s(k)*G_s(k) + Ys_s(k)*(sum_BG + mdot_s) - W_species(k)*omega_s(k)
        end do

    end subroutine s_surface_species_residual

    !> Surface energy residual: gas-side conduction balances heterogeneous reaction heat. Radiation and solid-side conduction are
    !! omitted.
    subroutine s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d
        real(wp), intent(in)  :: Ys_s(num_species)
        real(wp), intent(out) :: R_energy
        real(wp)              :: mw_s, rho_s, k_s, q_rxn

        call get_mixture_molecular_weight(Ys_s, mw_s)
        rho_s = pres*mw_s/(gas_constant*T_s)

        call get_mixture_thermal_conductivity_mixavg(T_s, Ys_s, k_s)
        call get_surface_reaction_heat_flux(rho_s, T_s, Ys_s, q_rxn)

        R_energy = k_s*(T_s - T_IP)/d - q_rxn

    end subroutine s_surface_energy_residual

    !> Assemble the Newton residual for Ns species, with temperature appended only when it is solved. k_bath selects which species
    !! balance the sum constraint displaces; see s_solve_surface for why it must be the most abundant one.
    subroutine s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                  & energy_scale, R, R_species, omega_s, mdot_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_s, d, flux_scale, energy_scale
        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species), W_species(num_species)
        integer, intent(in)   :: k_bath
        logical, intent(in)   :: solve_temperature
        real(wp), intent(out) :: R(num_species + 1), R_species(num_species), omega_s(num_species), mdot_s
        real(wp)              :: R_energy
        integer               :: k

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)

        R = 0._wp
        do k = 1, num_species
            R(k) = R_species(k)/flux_scale
        end do
        R(k_bath) = sum(Ys_s) - 1._wp

        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            R(num_species + 1) = R_energy/energy_scale
        end if

    end subroutine s_surface_residual

    !> Largest temperature mirror weight, at most 1, that keeps the ghost consistent with every fluid cell it shares a face with.
    !! Conduction across that face is the two-point difference (T_n - T_g)/h, so the line through the ghost (a behind the wall) and
    !! the cell (b in front) must not reach the wall beyond T_s: T_g - T_s <= (a/b)(T_s - T_n) on the mirror's side, and T_g = T_s
    !! once the cell is past T_s. The image point sits about a cell out, so without this an under-resolved layer leaves the mirror
    !! set by the second gas cell and drives the first past the wall temperature. A resolved, convex profile puts the mirror inside
    !! the bound, so it leaves the mirror untouched there. A cell closer than a quarter of its spacing to the wall plane is left
    !! out: its temperature is within noise of T_s, and the bound, whose slope there is a/b, then switched theta between 0 and 1
    !! from step to step and kept a steady conduction annulus from converging.
    function f_face_theta(q_prim_vf, markers, gp, T_IP, T_s) result(theta)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(integer_field), intent(in) :: markers  !< ib_markers, passed so the device routine reads the present copy
        type(ghost_point), intent(in) :: gp
        real(wp), intent(in) :: T_IP, T_s
        real(wp) :: theta, a, b, h, T_n, mw_n, r
        real(wp), dimension(3) :: nrm
        real(wp), dimension(${NUM_SPECIES}$) :: Ys_n
        integer, dimension(3) :: c
        integer :: dir, sgn, q

        theta = 1._wp
        if (f_approx_equal(T_s, T_IP)) return
        a = abs(real(gp%levelset, kind=wp))
        nrm = gp%levelset_norm/max(sqrt(sum(gp%levelset_norm**2)), sgm_eps)
        do dir = 1, num_dims
            do sgn = -1, 1, 2
                c = gp%loc
                c(dir) = c(dir) + sgn
                if (markers%sf(c(1), c(2), c(3)) /= 0) cycle
                select case (dir)
                case (1)
                    h = x_cc(c(1)) - x_cc(gp%loc(1))
                case (2)
                    h = y_cc(c(2)) - y_cc(gp%loc(2))
                case default
                    h = z_cc(c(3)) - z_cc(gp%loc(3))
                end select
                ! Wall distance of the neighbor from the ghost's tangent plane (it can come out <= 0 at an edge)
                b = h*nrm(dir) - a
                if (b < 0.25_wp*abs(h)) cycle
                do q = 1, num_species
                    Ys_n(q) = q_prim_vf(eqn_idx%species%beg + q - 1)%sf(c(1), c(2), c(3))
                end do
                call get_mixture_molecular_weight(Ys_n, mw_n)
                T_n = q_prim_vf(eqn_idx%E)%sf(c(1), c(2), c(3))*mw_n/(q_prim_vf(eqn_idx%cont%beg)%sf(c(1), c(2), c(3))*gas_constant)
                r = (T_s - T_n)/(T_s - T_IP)
                if (r <= 0._wp) then
                    theta = 0._wp
                else
                    theta = min(theta, (a/b)*r)
                end if
            end do
        end do

    end function f_face_theta

    !> Ghost state for an immersed surface, from the image-point state and the surface state.
    !!
    !! The natural closure is the linear mirror phi_g = 2*phi_s - phi_IP, which is what makes the midpoint of the ghost/image
    !! pair reproduce the surface value. Its -1 coefficient on phi_IP is unconditional, though, so it drives strictly positive
    !! quantities negative whenever the surface value sits far below the image-point value: a cold wall in hot gas gives a
    !! negative ghost temperature and hence a negative ghost density, and a species the surface consumes faster than half the
    !! free-stream value gives a negative mass fraction. Rather than clamp the result -- which would break sum(Y) = 1 and hide
    !! the excursion -- the ghost state is the largest convex blend of the mirror back toward the surface value that stays
    !! physical:
    !!
    !!     phi_g(theta) = phi_s + theta*(phi_s - phi_IP),   theta in [0, 1]
    !!
    !! theta = 1 is the full second-order mirror; theta = 0 is the first-order Dirichlet ghost phi_g = phi_s, which is
    !! Gibou et al. (JCP 176:205, 2002) Eq. 17 and is what they likewise fall back to where the linear form is ill-behaved.
    !! Each species and the temperature get their own theta (s_blend_ghost_species keeps sum(Y) = 1 through the bath species).
    !! Sharing one would let a trace radical the surface consumes to ~1e-9 (O and OH at a burning carbon wall) drag every other
    !! mirror down with it, turning a 2100 K ghost into 1205 K over a species whose own excursion was 1e-9. The cost of limiting
    !! is small: the mirror already delivers only first-order wall flux (Ezra et al., Int. J. Heat Mass Transfer, 2025),
    !! and a boundary closure one order below the interior scheme retains the interior convergence rate (Gustafsson,
    !! Math. Comp. 29:396, 1975).
    !!
    !! This bounds the ghost state; it does not make the boundary strictly conservative. Baskaya et al. (Computers & Fluids
    !! 270:106134, 2024) trace the same cold-wall, large-gradient regime in ablation to ghost-cell mass conservation error
    !! that surfaces as spurious blowing, and resolve it only by moving to a flux-based cut-cell boundary -- a different
    !! discretization from this one, not a tuning of it.
    subroutine s_blend_ghost_temperature(T_IP, T_s, T_lo, T_hi, theta_max, T_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: T_IP, T_s, T_lo, T_hi
        real(wp), intent(in)  :: theta_max  !< Largest mirror weight the caller allows (f_face_theta), at most 1
        real(wp), intent(out) :: T_g
        ! Stop short of the bound so roundoff cannot cross it.
        real(wp), parameter :: blend_safety = 0.9_wp
        real(wp)            :: theta_T

        ! Largest theta <= theta_max keeping T_g in [T_lo, T_hi]. The bounds are compared before dividing, so +-huge() cannot
        ! overflow.
        theta_T = theta_max
        if (T_s < T_IP .and. blend_safety*(T_s - T_lo) < theta_T*(T_IP - T_s)) then
            theta_T = blend_safety*(T_s - T_lo)/(T_IP - T_s)
        else if (T_s > T_IP .and. blend_safety*(T_hi - T_s) < theta_T*(T_s - T_IP)) then
            theta_T = blend_safety*(T_hi - T_s)/(T_s - T_IP)
        end if
        theta_T = max(theta_T, 0._wp)

        T_g = T_s + theta_T*(T_s - T_IP)

    end subroutine s_blend_ghost_temperature

    !> Species half of the ghost-state blend: each species takes its own positivity-limited theta and the bath species (the most
    !! abundant at the surface) absorbs the difference, so sum(Y_g) = sum(Y_s). One theta shared by every species collapsed to
    !! ~0.006 at a burning carbon wall, where the surface consumes the trace radicals (O, OH), and made every species ghost first
    !! order. The shared theta is the fallback if the bath species would go negative.
    subroutine s_blend_ghost_species(Ys_IP, Ys_s, Ys_g)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: Ys_IP(num_species), Ys_s(num_species)
        real(wp), intent(out) :: Ys_g(num_species)
        real(wp), parameter   :: blend_safety = 0.9_wp
        real(wp)              :: theta_Y
        integer               :: k, k_bath

        call s_pick_bath_species(Ys_s, k_bath)
        Ys_g(k_bath) = Ys_s(k_bath)
        do k = 1, num_species
            if (k == k_bath) cycle
            theta_Y = 1._wp
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(1._wp, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
            Ys_g(k_bath) = Ys_g(k_bath) + Ys_s(k) - Ys_g(k)
        end do
        if (Ys_g(k_bath) >= 0._wp) return

        theta_Y = 1._wp
        do k = 1, num_species
            if (Ys_s(k) < Ys_IP(k)) theta_Y = min(theta_Y, blend_safety*Ys_s(k)/(Ys_IP(k) - Ys_s(k)))
        end do
        do k = 1, num_species
            Ys_g(k) = Ys_s(k) + theta_Y*(Ys_s(k) - Ys_IP(k))
        end do

    end subroutine s_blend_ghost_species

    !> Index of the most abundant species, i.e. the balance the sum constraint displaces.
    !!
    !! Closing sum(Y) = 1 by dropping one species' flux balance is standard (Surface CHEMKIN), but the dropped species absorbs
    !! the roundoff of every other balance, so it must be the bath gas. Dropping a fixed index instead -- the last species in
    !! the mechanism -- lands that error on whatever the mechanism happens to list last, which for the reduced GRI mechanism
    !! shipped with the reacting-surface example is H2O2, a trace radical whose own surface balance then goes unenforced.
    !! Cantera's solveSP re-scans for the largest species each iteration (evalSurfLarge) for this reason.
    subroutine s_pick_bath_species(Ys_s, k_bath)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: Ys_s(num_species)
        integer, intent(out) :: k_bath
        integer              :: k

        k_bath = 1
        do k = 2, num_species
            if (Ys_s(k) > Ys_s(k_bath)) k_bath = k
        end do

    end subroutine s_pick_bath_species

    !> Newton solve for a reacting surface. thermal_bc=0: zero-normal-gradient T; 1: prescribed T; 2: energy balance.
    subroutine s_solve_surface(pres, T_IP, T_wall, d, Ys_IP, W_species, thermal_bc, Ys_s, T_s, mdot_s, converged)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: pres, T_IP, T_wall, d
        real(wp), intent(in)  :: Ys_IP(num_species), W_species(num_species)
        integer, intent(in)   :: thermal_bc
        real(wp), intent(out) :: Ys_s(num_species), T_s, mdot_s
        logical, intent(out)  :: converged
        integer, parameter    :: max_iter = 30, max_backtrack = 20
        real(wp), parameter   :: fd_eps_Y = 1.e-7_wp, fd_eps_T = 1.e-6_wp
        ! The Jacobian is a forward difference, so the Newton update carries O(sqrt(epsilon)) noise and the iteration turns
        ! linear near the root. The residual itself is evaluated exactly, so a tighter tolerance is not unreachable, but it is
        ! reached by stalling rather than converging -- which the caller then reports as a failed surface. Cantera's solveSP
        ! runs the same 1e-7 forward-difference step against a 1e-4 relative tolerance; 1e-6 here stays well inside that.
        real(wp), parameter :: tol = 1.e-6_wp
        real(wp), parameter :: Y_tol = 100._wp*epsilon(1._wp)
        real(wp)            :: A(num_species + 1, num_species + 1), rhs(num_species + 1), delta(num_species + 1)
        real(wp)            :: R(num_species + 1), R_pert(num_species + 1), R_trial(num_species + 1)
        real(wp)            :: R_species(num_species), R_species_pert(num_species), R_species_trial(num_species)
        real(wp)            :: omega_s(num_species), omega_pert(num_species), omega_trial(num_species)
        real(wp)            :: Ys_pert(num_species), Ys_trial(num_species)
        real(wp)            :: mdot_pert, mdot_trial, T_pert, T_trial
        real(wp)            :: flux_scale, energy_scale, R_energy, dx, lambda, norm_R, norm_trial
        logical             :: solve_temperature, linear_success, accepted
        integer             :: nsolve, iter, j, iback, k_bath, k_bath_prev

        converged = .false.
        Ys_s = Ys_IP

        select case (thermal_bc)
        case (0)
            T_s = T_IP
            solve_temperature = .false.
        case (1)
            T_s = T_wall
            solve_temperature = .false.
        case (2)
            T_s = min(max(T_IP, T_surface_min), T_surface_max)
            solve_temperature = .true.
        case default
            T_s = T_IP
            omega_s = 0._wp
            mdot_s = 0._wp
            return
        end select

        nsolve = num_species + merge(1, 0, solve_temperature)

        call s_pick_bath_species(Ys_s, k_bath)

        call s_surface_species_residual(pres, T_s, d, Ys_IP, Ys_s, W_species, R_species, omega_s, mdot_s)
        flux_scale = max(maxval(abs(R_species)), 1.e-12_wp)
        energy_scale = 1._wp
        if (solve_temperature) then
            call s_surface_energy_residual(pres, T_IP, T_s, d, Ys_s, R_energy)
            energy_scale = max(abs(R_energy), 1._wp)
        end if

        call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, energy_scale, &
                                & R, R_species, omega_s, mdot_s)
        norm_R = maxval(abs(R(1:nsolve)))
        if (norm_R < tol) then
            converged = .true.
            return
        end if

        do iter = 1, max_iter
            do j = 1, num_species
                Ys_pert = Ys_s
                T_pert = T_s
                dx = fd_eps_Y*max(abs(Ys_s(j)), 1._wp)
                if (Ys_s(j) + dx > 1._wp) dx = -dx
                Ys_pert(j) = Ys_pert(j) + dx

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,j) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end do

            if (solve_temperature) then
                Ys_pert = Ys_s
                dx = fd_eps_T*max(abs(T_s), 1._wp)
                T_pert = T_s + dx
                if (T_pert > T_surface_max) then
                    dx = -dx
                    T_pert = T_s + dx
                end if

                call s_surface_residual(pres, T_IP, T_pert, d, Ys_IP, Ys_pert, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R_pert, R_species_pert, omega_pert, mdot_pert)
                A(1:nsolve,nsolve) = (R_pert(1:nsolve) - R(1:nsolve))/dx
            end if

            rhs(1:nsolve) = -R(1:nsolve)
            call s_solve_surface_linear_system(A, rhs, delta, nsolve, linear_success)
            if (.not. linear_success) return

            lambda = 1._wp
            accepted = .false.
            do iback = 1, max_backtrack
                Ys_trial = Ys_s + lambda*delta(1:num_species)
                T_trial = T_s
                if (solve_temperature) T_trial = T_s + lambda*delta(nsolve)

                if (minval(Ys_trial) < -Y_tol .or. maxval(Ys_trial) > 1._wp + Y_tol .or. T_trial < T_surface_min &
                    & .or. T_trial > T_surface_max) then
                    lambda = 0.5_wp*lambda
                    cycle
                end if

                where (Ys_trial < 0._wp) Ys_trial = 0._wp
                where (Ys_trial > 1._wp) Ys_trial = 1._wp

                call s_surface_residual(pres, T_IP, T_trial, d, Ys_IP, Ys_trial, W_species, k_bath, solve_temperature, &
                                        & flux_scale, energy_scale, R_trial, R_species_trial, omega_trial, mdot_trial)
                norm_trial = maxval(abs(R_trial(1:nsolve)))
                if (norm_trial < norm_R) then
                    accepted = .true.
                    exit
                end if
                lambda = 0.5_wp*lambda
            end do

            if (.not. accepted) return

            Ys_s = Ys_trial
            T_s = T_trial
            R = R_trial
            R_species = R_species_trial
            omega_s = omega_trial
            mdot_s = mdot_trial
            norm_R = norm_trial

            ! The composition has moved, so the species the sum constraint displaces may no longer be the most abundant one.
            ! Re-closing changes what R means, so the residual is rebuilt on the rare steps where the choice actually changes.
            k_bath_prev = k_bath
            call s_pick_bath_species(Ys_s, k_bath)
            if (k_bath /= k_bath_prev) then
                call s_surface_residual(pres, T_IP, T_s, d, Ys_IP, Ys_s, W_species, k_bath, solve_temperature, flux_scale, &
                                        & energy_scale, R, R_species, omega_s, mdot_s)
                norm_R = maxval(abs(R(1:nsolve)))
            end if

            if (norm_R < tol) then
                converged = .true.
                return
            end if
        end do

    end subroutine s_solve_surface

    !> Small dense linear solve with partial pivoting for the local surface Newton system.
    subroutine s_solve_surface_linear_system(A, b, x, nsolve, success)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)     :: nsolve
        real(wp), intent(inout) :: A(num_species + 1, num_species + 1)
        real(wp), intent(inout) :: b(num_species + 1)
        real(wp), intent(out)   :: x(num_species + 1)
        logical, intent(out)    :: success
        real(wp)                :: factor, pivot_value, tmp, row_tmp(num_species + 1)
        integer                 :: i, j, k, pivot

        success = .true.
        x = 0._wp

        do k = 1, nsolve - 1
            pivot = k
            pivot_value = abs(A(k, k))
            do i = k + 1, nsolve
                if (abs(A(i, k)) > pivot_value) then
                    pivot = i
                    pivot_value = abs(A(i, k))
                end if
            end do
            if (pivot_value <= epsilon(1._wp)) then
                success = .false.
                return
            end if

            if (pivot /= k) then
                row_tmp(1:nsolve) = A(k,1:nsolve)
                A(k,1:nsolve) = A(pivot,1:nsolve)
                A(pivot,1:nsolve) = row_tmp(1:nsolve)
                tmp = b(k)
                b(k) = b(pivot)
                b(pivot) = tmp
            end if

            do i = k + 1, nsolve
                factor = A(i, k)/A(k, k)
                A(i, k) = 0._wp
                do j = k + 1, nsolve
                    A(i, j) = A(i, j) - factor*A(k, j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(A(nsolve, nsolve)) <= epsilon(1._wp)) then
            success = .false.
            return
        end if

        x(nsolve) = b(nsolve)/A(nsolve, nsolve)
        do i = nsolve - 1, 1, -1
            tmp = b(i)
            do j = i + 1, nsolve
                tmp = tmp - A(i, j)*x(j)
            end do
            if (abs(A(i, i)) <= epsilon(1._wp)) then
                success = .false.
                return
            end if
            x(i) = tmp/A(i, i)
        end do

    end subroutine s_solve_surface_linear_system

end module m_ibm
