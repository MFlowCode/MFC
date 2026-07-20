!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic
    use m_constants, only: recon_type_weno, recon_type_muscl, muscl_order_first_order, time_stepper_rk3, BC_RIEMANN_EXTRAP

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        call s_check_inputs_time_stepping

        @:PROHIBIT(chemistry .and. chem_params%reaction_substeps < 0, &
                   & "chem_params%reaction_substeps must be >= 0 (0 = reaction source in the flow RHS; > 0 = operator-split sub-stepping)")

        @:PROHIBIT(chemistry .and. igr .and. chem_params%reaction_substeps > 0, &
                   & "operator-split reaction sub-stepping (reaction_substeps > 0) is not supported with igr: the reactor reads the post-flow (rho, e, T) state, which the IGR update path does not guarantee")

        @:PROHIBIT(chemistry .and. chem_params%adap_substeps .and. chem_params%reaction_substeps < 1, &
                   & "chem_params%adap_substeps requires reaction_substeps >= 1 (the operator-split floor)")

        @:PROHIBIT(chemistry .and. chem_params%adap_substeps &
                   & .and. chem_params%reaction_substeps_max < chem_params%reaction_substeps, &
                   & "chem_params%reaction_substeps_max must be >= reaction_substeps when adap_substeps = T")

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")
        @:PROHIBIT(many_ib_patch_parallelism .and. .not. ib, "many_ib_patch_parallelism requires ib to be enabled")

        if (active_box) then
            @:PROHIBIT(recon_type /= recon_type_weno, "active_box requires WENO reconstruction")
            @:PROHIBIT(ib, "active_box is incompatible with immersed boundaries")
            @:PROHIBIT(acoustic_source, "active_box is incompatible with acoustic sources")
            @:PROHIBIT(bodyForces, "active_box is incompatible with body forces")
            @:PROHIBIT(bubbles_lagrange, "active_box is incompatible with Lagrangian bubbles")
            @:PROHIBIT(relax, "active_box is incompatible with phase change")
            @:PROHIBIT(igr, "active_box is incompatible with the IGR solver")
            @:PROHIBIT(time_stepper /= time_stepper_rk3, "active_box requires time_stepper = 3 (SSP-RK3)")
            @:PROHIBIT(viscous, "active_box is incompatible with viscous (no finite domain of dependence for the frozen exterior)")
            @:PROHIBIT(surface_tension, &
                       & "active_box is incompatible with surface_tension (nonlocal curvature coupling violates the static-uniform-exterior assumption)")
            @:PROHIBIT(cyl_coord, &
                       & "active_box is incompatible with cyl_coord (geometric source terms are nonzero for uniform flow; exterior is not static)")
            @:PROHIBIT(hypoelasticity, &
                       & "active_box is incompatible with hypoelasticity (stress source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(hyperelasticity, &
                       & "active_box is incompatible with hyperelasticity (stress source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(mhd, &
                       & "active_box is incompatible with mhd (magnetic field source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(chemistry, &
                       & "active_box is incompatible with chemistry (reactive source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(bubbles_euler, &
                       & "active_box is incompatible with bubbles_euler (cell-local bubble sources in a non-equilibrium ambient violate the static-uniform-exterior assumption)")
        end if

        @:PROHIBIT(sfc_partition_wrt .and. partition_tile_size < 1, "partition_tile_size must be >= 1")
        @:PROHIBIT(load_balance .and. .not. parallel_io, "load_balance requires parallel_io = T")
        @:PROHIBIT(load_balance .and. num_procs == 1, "load_balance requires more than one MPI rank")
        @:PROHIBIT(load_balance .and. file_per_process, &
                   & "load_balance is incompatible with file_per_process (per-rank restart files are sized for the equal decomposition; rereading them at rebalanced extents corrupts the fields)")

        if (amr) then
            @:PROHIBIT((.not. igr) .and. recon_type /= recon_type_weno, "amr requires WENO reconstruction (or the IGR solver)")
            @:PROHIBIT(time_stepper /= time_stepper_rk3, "amr requires time_stepper = 3 (SSP-RK3)")
            ! 6-equation support: the internal-energy equations prolong/restrict on the generic
            ! conservative path and the per-stage pressure relaxation (cell-local) also runs on
            ! each fine block, mirroring the coarse stage order.
            @:PROHIBIT(model_eqns /= 2 .and. model_eqns /= 3, "amr requires model_eqns = 2 (5-equation) or 3 (6-equation)")
            ! Non-polytropic QBMM carries a pb/mv quadrature side-state whose coarse<->fine coupling (prolong +
            ! restriction) is distributed across ranks for SINGLE-LEVEL AMR via the pb/mv P2P gather/scatter
            ! (s_amr_gather_coarse_patch_pbmv / s_amr_scatter_pbmv, mirroring the q_cons distribution). The MULTI-level
            ! (parent-side) pb/mv coupling is not yet distributed, so amr_max_level > 1 stays fail-closed at np>=2.
            @:PROHIBIT(qbmm .and. (.not. polytropic) .and. amr_max_level > 1 .and. num_procs > 1, &
                       & "amr with non-polytropic QBMM on more than one MPI rank is only supported at amr_max_level = 1: the multi-level (parent-side) pb/mv quadrature side-state coarse/fine coupling is not yet distributed. Run multi-level non-polytropic QBMM on a single rank.")
            ! Riemann-extrapolation BCs modify the WENO coefficient rows near the domain boundary;
            ! the fine advance reuses (or block-locally recomputes) those arrays, and neither form
            ! can carry the coarse boundary special-casing onto an interior block correctly
            @:PROHIBIT(bc_x%beg == BC_RIEMANN_EXTRAP .or. bc_x%end == BC_RIEMANN_EXTRAP .or. (n_glb > 0 &
                       & .and. (bc_y%beg == BC_RIEMANN_EXTRAP .or. bc_y%end == BC_RIEMANN_EXTRAP)) .or. (p_glb > 0 &
                       & .and. (bc_z%beg == BC_RIEMANN_EXTRAP .or. bc_z%end == BC_RIEMANN_EXTRAP)), &
                       & "amr does not support Riemann-extrapolation boundary conditions (bc = -4): they alter the WENO coefficient rows near the boundary, which the fine-block reconstruction cannot inherit correctly")
            @:PROHIBIT(num_fluids > 1 .and. (.not. mpp_lim) .and. (.not. bubbles_lagrange), &
                       & "amr with num_fluids > 1 requires mpp_lim (its volume-fraction clamp+renormalize maintains coarse/fine alpha consistency); Lagrangian bubbles are exempt (their alphas sum to the local liquid fraction and prolong without the sum-to-one closure)")
            @:PROHIBIT(surface_tension, &
                       & "amr does not support surface_tension: the capillary force depends on the interface normal (grad-c direction), which the conservative-linearly-prolonged fine ghost color cannot reproduce consistently with the coarse solver across a 2:1 coarse/fine boundary - a growing spurious seam current results")
            ! hypoelasticity is supported: stress components prolong via the generic conservative-linear
            ! path and the swap/restore recomputes the spacing-dependent FD coefficients per grid
            @:PROHIBIT(hyperelasticity, "amr does not support hyperelasticity")
            ! MHD stays gated ON MEASURED EVIDENCE (not just caution): B/psi ride the generic
            ! conservative machinery structurally, but the per-component prolongation/reflux is
            ! not divergence-preserving - on a magnetized 2D Brio-Wu the c/f seam acts as a
            ! continuous O(1) monopole source that GLM cleaning spreads but cannot remove
            ! (max|divB| 0.53 in the block interior and 0.36 far-field vs the no-AMR run's
            ! 1.4e-3 cleaning background; HLLD, which has no GLM coupling at all, NaNs
            ! outright). Supporting MHD needs divergence-preserving (constrained-transport
            ! class) prolongation and reflux for B.
            ! 1D MHD/RMHD is exempt: div(B) = d(Bx)/dx and 1D evolves only By/Bz (Bx is the
            ! uniform Bx0 parameter), so div(B) is IDENTICALLY zero - the failure mode above
            ! is structurally absent and By/Bz reflux/restrict as ordinary conserved scalars.
            @:PROHIBIT(mhd .and. n > 0, &
                       & "amr with mhd is 1D-only (the coarse/fine seam is not divergence-preserving for B; in 1D div(B) = 0 by construction)")
            ! IGR is supported with restriction-only coarse/fine coupling (stage 1): the fine
            ! block runs its own fixed-iteration sigma solve seeded and Dirichlet-bounded by the
            ! converged coarse sigma; the Berger-Colella reflux is not yet captured from the
            ! fused IGR flux kernels, so seam conservation is truncation-order, not exact
            @:PROHIBIT(igr .and. amr_subcycle, "amr_subcycle with the IGR solver is not yet supported (lockstep only)")
            ! Lagrangian bubbles are supported with the cloud EXCLUDED from fine blocks (two-way
            ! coupling lives on the coarse grid): regrid suppresses tags and clips boxes around
            ! the cloud's padded bbox, and a per-stage guard aborts if the cloud reaches a block
            ! 2D axisymmetric is supported: the geometric sources read the live grid arrays the fine
            ! swap replaces, and the axis-singularity viscous treatment is skipped on fine blocks
            ! (blocks cannot touch the axis - the domain-edge clamp keeps them buff_size inside).
            ! 3D cylindrical stays gated: its per-stage azimuthal Fourier filter is a global
            ! operation incompatible with the block-local fine advance.
            @:PROHIBIT(cyl_coord .and. p > 0, &
                       & "amr with cyl_coord supports 2D axisymmetric only: the 3D cylindrical azimuthal Fourier filter is a global operation incompatible with the block-local fine advance")
            ! 2D axisymmetric conservation (radius-weighted restriction + area-weighted reflux) is implemented for the L0/L1 coarse
            ! frame only. Multi-level folds/refluxes in the PARENT-FINE frame (host-only per-block coords) and non-polytropic QBMM
            ! carries a pb/mv side-state whose fold-back is not radius-weighted - both stay fail-closed under cyl_coord.
            @:PROHIBIT(cyl_coord .and. amr_max_level > 1, &
                       & "amr with cyl_coord supports amr_max_level = 1 only: multi-level axisymmetric restriction/reflux in the parent-fine frame is not yet radius-weighted (conservation would drift)")
            @:PROHIBIT(cyl_coord .and. qbmm .and. (.not. polytropic), &
                       & "amr with cyl_coord and non-polytropic QBMM is not supported: the pb/mv quadrature side-state fold-back is not radius-weighted (conservation would drift)")
            ! non-polytropic QBMM: each block carries its own pb/mv quadrature side-state (prolonged
            ! piecewise-constant to preserve CHyQMOM realizability, advanced with the block's own rhs
            ! scratch, restricted back with the moments). The subcycle time-lerps the pb/mv ghost
            ! shell like the conservative ghosts, and the regrid bounces the side-state through the
            ! pb/mv_stor arrays exactly like q_cons - dynamic regrid and subcycling are supported.
            ! static-body IB AMR (SP20) + prescribed-motion moving bodies (SP21): fixed or analytically-moving
            ! (moving_ibm==1) bodies resolved on a static fine block. Multi-body (num_ibs>1) is supported - every body
            ! shares the one static block and reuses the multi-body-capable core IB setup. Force/torque-driven motion
            ! (moving_ibm==2) and STL geometry remain gated (unvalidated).
            @:PROHIBIT(ib .and. any(patch_ib(1:num_ibs)%moving_ibm == 2), &
                       & "amr with ib supports static or prescribed-motion (moving_ibm=1) bodies only; force-driven moving IB (moving_ibm=2) under amr is not yet validated")
            @:PROHIBIT(ib .and. any(patch_ib(1:num_ibs)%geometry == 12), &
                       & "amr with ib does not support STL-model geometry (not yet validated)")
            ! dynamic regrid with bodies (static or prescribed-motion): candidate boxes expand to
            ! fully contain every body at its LIVE position (partial coverage is an untested
            ! regime), overlapping expansions merge, and the fine IB state is rebuilt from the
            ! geometry after every regrid. Between regrids a moving body's containment is
            ! guarded per substage (abort if it reaches the block boundary).
            ! active_box is supported (np=1 by active_box's own gate): blocks must sit strictly
            ! inside the monotonically-growing active window (init check + regrid clamp; the
            ! windowed coarse update would drop reflux corrections at faces outside it), and
            ! the fine advance disables the coarse-indexed windowing on the swapped block grid
            ! no acoustic_source gate here: acoustic sources act on the coarse grid only (their spatial support is precomputed as
            ! coarse cell indices). A startup check aborts if the support overlaps the user-placed
            ! initial block; the dynamic regrid keeps its own boxes clear of the support (tags are
            ! suppressed over it and candidate boxes are clipped), so the source region stays coarse.
            @:PROHIBIT(any(amr_block_beg(1:num_dims) < 0), "amr_block_beg must be >= 0")
            @:PROHIBIT(amr_block_end(1) > m_glb .or. (num_dims >= 2 .and. amr_block_end(2) > n_glb) .or. (num_dims >= 3 &
                       & .and. amr_block_end(3) > p_glb), "amr_block_end must be <= global cell max per axis")
            @:PROHIBIT(any(amr_block_end(1:num_dims) <= amr_block_beg(1:num_dims)), &
                       & "amr_block_end must exceed amr_block_beg on each active axis")
            @:PROHIBIT(ref_ratio*(amr_block_end(1) - amr_block_beg(1) + 1) - 1 > m_glb, &
                       & "amr fine x-extent exceeds the base grid (module scratch is sized to the base)")
            @:PROHIBIT(num_dims >= 2 .and. ref_ratio*(amr_block_end(2) - amr_block_beg(2) + 1) - 1 > n_glb, &
                       & "amr fine y-extent exceeds the base grid")
            @:PROHIBIT(num_dims >= 3 .and. ref_ratio*(amr_block_end(3) - amr_block_beg(3) + 1) - 1 > p_glb, &
                       & "amr fine z-extent exceeds the base grid")
            @:PROHIBIT(amr_regrid_int < 0, "amr_regrid_int must be >= 0")
            @:PROHIBIT(amr_regrid_int > 0 .and. amr_tag_eps <= 0._wp, "amr_tag_eps must be > 0 when regridding")
            @:PROHIBIT(amr_regrid_int > 0 .and. amr_buf < 1, "amr_buf must be >= 1 when regridding")
            @:PROHIBIT(amr_max_blocks < 1, "amr_max_blocks must be >= 1")
            @:PROHIBIT(amr_max_level < 1, "amr_max_level must be >= 1")
            @:PROHIBIT(amr_max_level > 1 .and. amr_max_blocks < 2, &
                       & "multi-level AMR (amr_max_level > 1) needs amr_max_blocks >= 2 (at least one level-1 block plus one nested level-2 block); a run-time abort catches the tiled case where even more blocks are required")
            @:PROHIBIT(amr_max_level > 1 .and. ib .and. num_procs > 1, &
                       & "multi-level AMR (amr_max_level > 1) with immersed boundaries is only supported at num_procs = 1 (the fine-IB image-point stencil is not decomposition-exact across a rank seam)")
            @:PROHIBIT(amr_max_level > 1 .and. ib .and. any(patch_ib(1:num_ibs)%moving_ibm /= 0), &
                       & "multi-level AMR (amr_max_level > 1) with a MOVING immersed body is not yet supported; use a static body")
            @:PROHIBIT(amr_regrid_int == 0 .and. amr_max_level > 2, &
                       & "static multi-level AMR (amr_regrid_int = 0) nests exactly one level-2 block in block 1, so it supports at most amr_max_level = 2; use amr_regrid_int > 0 for deeper or multi-block nesting")
            @:PROHIBIT(amr_cluster_eff <= 0._wp .or. amr_cluster_eff > 1._wp, &
                       & "amr_cluster_eff must satisfy 0 < amr_cluster_eff <= 1")
            @:PROHIBIT(ref_ratio /= 2 .and. ref_ratio /= 4, "ref_ratio must be 2 or 4")
            @:PROHIBIT(ref_ratio /= 2 .and. (amr_max_level > 1 .or. amr_subcycle), &
                       & "ref_ratio /= 2 is only supported at amr_max_level = 1 without subcycling (v1)")
        end if
        @:PROHIBIT(.not. amr .and. amr_regrid_int > 0, "amr_regrid_int requires amr")
        @:PROHIBIT(amr_subcycle .and. .not. amr, "amr_subcycle requires amr")
        @:PROHIBIT(amr_subcycle .and. cfl_dt, "amr_subcycle requires a fixed dt (cfl_dt not supported)")
        ! Subcycled fine advance at np>1 needs the block-to-block fine-fine seam halo (s_amr_fine_fine_halo) run PER SUBSTEP:
        ! max_grid_size TILING can split a feature into ADJACENT same-level sub-blocks, and the halo overwrites their shared-face
        ! ghosts with the neighbour's fine interior so both sides compute a MATCHING seam flux (else mass leaks at the seam).
        ! s_amr_advance_fine_subcycle_all advances all LEVEL-1 blocks stage-by-stage in lockstep with the halo interposed, so
        ! single-level subcycle np>1 is conservation-safe. The level-2 children still advance per-block (s_amr_advance_children),
        ! so L2-L2 seams are not yet reconciled - keep multi-level (amr_max_level > 1) subcycle gated at np>1 until the recursive
        ! per-substep L2 halo lands. np=1 never tiles into adjacent blocks (halo skipped there, byte-identical to before).
        @:PROHIBIT(amr_subcycle .and. amr_regrid_int > 0 .and. num_procs > 1 .and. amr_max_level > 1, &
                   & "multi-level (amr_max_level > 1) amr_subcycle with dynamic regrid is not yet conservation-safe at num_procs > 1: the level-2 seam halo is per-block, not lockstep (single-level subcycling IS supported at np > 1). Use amr_subcycle = F (lock-step) for multi-level dynamic multi-rank runs")

        if (num_particle_clouds > 0) then
            call s_check_inputs_particle_clouds
        end if

        if (synthetic_turbulence) then
            call s_check_inputs_synthetic_turbulence
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(muscl_order == muscl_order_first_order .and. int_comp > 0, &
                   & "int_comp requires muscl_order >= 2 (muscl_order=1 leaves the reconstruction workspace uninitialised)")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    !> Checks that each active particle cloud has a valid packing_method specified
    impure subroutine s_check_inputs_particle_clouds

        integer          :: i
        character(len=5) :: idxStr

        do i = 1, num_particle_clouds
            call s_int_to_str(i, idxStr)
            @:PROHIBIT(particle_cloud(i)%packing_method == dflt_int, &
                       & "particle_cloud("//trim(idxStr) &
                       & //")%packing_method must be specified (1 = rejection sampling, 2 = lattice)")
            @:PROHIBIT(particle_cloud(i)%packing_method /= 1 .and. particle_cloud(i)%packing_method /= 2, &
                       & "particle_cloud("//trim(idxStr) //")%packing_method must be 1 (rejection sampling) or 2 (lattice)")
        end do

    end subroutine s_check_inputs_particle_clouds

    !> Checks that each active synthetic-turbulence forcing zone has a fully specified position and a positive size in every active
    !! dimension
    impure subroutine s_check_inputs_synthetic_turbulence

        integer          :: i, d
        character(len=5) :: idxStr

        @:PROHIBIT(num_turbulent_sources <= 0, "num_turbulent_sources must be > 0 when synthetic_turbulence is enabled")

        do i = 1, num_turbulent_sources
            call s_int_to_str(i, idxStr)
            do d = 1, num_dims
                @:PROHIBIT(f_is_default(turb_pos(i, d)), &
                           & "turb_pos("//trim(idxStr) &
                           & //",:) must be specified for all num_dims when synthetic_turbulence is enabled")
                @:PROHIBIT(f_is_default(synth_L(i, d)) .or. synth_L(i, d) <= 0._wp, &
                           & "synth_L("//trim(idxStr)//",:) must be positive for all num_dims when synthetic_turbulence is enabled")
            end do
        end do

    end subroutine s_check_inputs_synthetic_turbulence

end module m_checker
