!>
!! @file
!! @brief Contains module m_projection

#:include 'case.fpp'
#:include 'macros.fpp'

#! Multigrid row at flat index idx of a ghosted level (y and z strides sy, sz): dg = diagonal, nb = conductance-weighted
! neighbor sum, with the values read sh past idx (a block with the level's layout stored elsewhere in mg_e). Each cell stores the
! conductance of its low faces; boundary faces carry zero unless they are rank or periodic
#
#! seams, whose neighbor values sit in the ghost layer
#:def MG_ROW(sh='0')
    dg = mg_d(idx) + mg_kx(idx) + mg_kx(idx + 1)
    nb = mg_kx(idx)*mg_e(idx + ${sh}$ - 1) + mg_kx(idx + 1)*mg_e(idx + ${sh}$ + 1)
    if (num_dims > 1) then
        dg = dg + mg_ky(idx) + mg_ky(idx + sy)
        nb = nb + mg_ky(idx)*mg_e(idx + ${sh}$ - sy) + mg_ky(idx + sy)*mg_e(idx + ${sh}$ + sy)
    end if
    if (num_dims > 2) then
        dg = dg + mg_kz(idx) + mg_kz(idx + sz)
        nb = nb + mg_kz(idx)*mg_e(idx + ${sh}$ - sz) + mg_kz(idx + sz)*mg_e(idx + ${sh}$ + sz)
    end if
#:enddef

#! First (lo) and last (hi) fine child of coarse index ic along a direction with nc coarse and nf fine cells: an odd last child
#! folds into the last coarse cell, and ic = nc (the ghost past the end) maps to the fine high face
#:def MG_CHILDREN(lo, hi, ic, nc, nf)
    ${lo}$ = 2*${ic}$
    ${hi}$ = merge(${nf}$ - 1, 2*${ic}$ + 1, ${ic}$ == ${nc}$ - 1)
    if (${ic}$ >= ${nc}$) ${hi}$ = 2*${ic}$
#:enddef

#! Flat index of cell (i, j, k) on the level whose metadata is in off, ex, ey, gx, gy, gz
#:def MG_IX(i, j, k)
    off + ((${k}$ + gz)*ey + ${j}$ + gy)*ex + ${i}$ + gx + 1
#:enddef

!> All-Mach pressure projection (Fuster & Popinet, JCP 374, 2018). Advection uses a persistent face velocity; the pressure then
!! solves p - rho*c^2*tau^2*div(rho_f^-1 grad p) = p_adv - rho*c^2*tau*div(u*_f) with div, grad and the Laplacian all taken on
!! faces, so they compose exactly and the projected face velocity satisfies the discrete pressure equation. The solve is PCG
!! preconditioned by a geometric multigrid V-cycle coupled across ranks: each level exchanges one ghost layer with its neighbors
!! (one aggregated message per neighbor), coarsens to one cell per rank, and every rank continues the hierarchy on that gathered
!! rank-level problem: solved directly when small, else by one host V-cycle.
!> @brief All-Mach pressure projection
module m_projection

#ifdef MFC_MPI
    use mpi
#endif
    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_boundary_common
    use m_eos
    use m_body_forces, only: s_compute_acceleration
    use m_riemann_state, only: Re_avg_rsx_vf, vel_src_rsx_vf, s_compute_interface_reynolds
    use m_nvtx

    implicit none

    private; public :: s_initialize_projection_module, s_projection_rhs, s_projection_face_props, s_projection_heat, &
        & s_projection_apply, s_finalize_projection_module

    integer, parameter :: mg_maxlev = 24
    real(wp), parameter :: res_floor = 1.e2_wp*epsilon(1._wp)  !< residual round-off floor, relative to the right-hand side
    real(wp), allocatable, dimension(:,:,:,:) :: uf            !< face velocity; index j is the face between cells j and j+1
    real(wp), allocatable, dimension(:,:,:,:) :: uf0           !< the face velocity that transported the current stage
    real(wp), allocatable, dimension(:,:,:) :: divu, rhs_p     !< div of uf, and the pressure transport rate
    real(wp), allocatable, dimension(:,:,:) :: pflx            !< upwind pressure flux on the faces of one direction
    real(wp), allocatable, dimension(:,:,:) :: p_stage, p_step0
    real(wp), allocatable, dimension(:,:,:) :: rhoc            !< star density with one ghost layer
    real(wp), allocatable, dimension(:,:,:) :: dcoef, bvec     !< SPD system: D_c and right-hand side
    real(wp), allocatable, dimension(:,:,:) :: xs, rs, zs, qs  !< PCG vectors
    real(stp), allocatable, dimension(:,:,:), target :: pk     !< search direction, then solution, with ghosts for the halo
    type(scalar_field), dimension(1) :: pk_sf
    $:GPU_DECLARE(create='[pk_sf]')
    $:GPU_DECLARE(create='[uf, uf0, divu, rhs_p, pflx, p_stage, p_step0, rhoc, dcoef, bvec, xs, rs, zs, qs, pk]')
    !> Well-balanced surface tension, with one ghost layer: curvature (1) and |grad c| (2), both zero outside the interface band
    real(wp), allocatable, dimension(:,:,:,:) :: kap
    real(wp), allocatable, dimension(:,:,:,:) :: gnd  !< grad(alpha_1) (1:num_dims) and its magnitude (0), two ghost layers
    $:GPU_DECLARE(create='[kap, gnd]')

    !> Multigrid hierarchy, flattened, with one ghost layer in each active direction; level lv starts after mg_off(lv), x fastest.
    !! Every level coarsens by two (an odd size folds its last cell into the last coarse cell) down to the bottom level, the first
    !! with at most mg_bottom_max cells in all (or one cell per rank). PCG's search direction follows the levels in mg_e, at
    !! mg_poff, laid out as level 1
    integer, parameter                  :: mg_bottom_max = 128
    integer                             :: mg_nlev, mg_gx, mg_gy, mg_gz, mg_poff
    integer                             :: mg_blev             !< Level this solve stops at (s_mg_omega): mg_nlev, or the truncation
    real(wp), dimension(mg_maxlev)      :: mg_om               !< Coarse-correction scale into each level (s_mg_omega)
    real(wp), parameter                 :: mg_om_r0 = 0.02_wp  !< Helmholtz-to-Laplacian ratio that halves the scale's excess
    integer, dimension(mg_maxlev)       :: mg_nx, mg_ny, mg_nz, mg_off
    real(wp), allocatable, dimension(:) :: mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r
    $:GPU_DECLARE(create='[mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r]')
    !> K-cycle vectors of the coarse levels (indexed as mg_e; allocated with mg_kfull > 0): the right-hand side (kr), the first
    !! step's preconditioned direction (kc) and its image under the level operator (kv)
    real(wp), allocatable, dimension(:) :: mg_kr, mg_kc, mg_kv
    $:GPU_DECLARE(create='[mg_kr, mg_kc, mg_kv]')
    real(wp), parameter :: mg_kt = 0.25_wp  !< K-cycle skips its second step once the residual falls by this

    !> Ghost-layer exchange: side q = 1..6 is (x, y, z) x (low, high); mg_nbr(q) is the neighbor rank, this rank for a periodic
    !! seam, or -1. Each exchange sends one message per distinct neighbor, holding all its sides
    integer, dimension(6)               :: mg_nbr
    real(wp), allocatable, dimension(:) :: mg_sbuf, mg_rbuf
    $:GPU_DECLARE(create='[mg_sbuf, mg_rbuf]')
    !> Side q's action (0 none, 1 pack for another rank, 3 periodic self copy) and its offset in mg_sbuf (1) and mg_rbuf (2) by
    !! level; message i goes to rank mg_nlist(i), from mg_sbeg + 1 for mg_slen values, and arrives at mg_rbeg + 1 for mg_rlen
    integer                             :: mg_nmsg
    integer, dimension(6)               :: mg_smode, mg_nlist
    integer, dimension(6, mg_maxlev, 2) :: mg_boff
    integer, dimension(6, mg_maxlev)    :: mg_sbeg, mg_slen, mg_rbeg, mg_rlen
    $:GPU_DECLARE(create='[mg_smode, mg_boff]')

    !> Bottom level, gathered to every rank: crs_n cells in all, rank r's crs_cnt(r + 1) of them (a crs_sz(:, r + 1) block) numbered
    !! from crs_disp(r + 1) + 1, x fastest. crs_n is at most mg_bottom_max, or num_procs past that, when it is the rank grid itself.
    !! Every rank continues the hierarchy on it: level lv holds crs_ln(lv) cells from crs_loff(lv) + 1, at rank-grid coordinates
    !! crs_co, as a stencil (diagonal crs_ad, conductance crs_ak(q, i) to cell crs_aj(q, i), 0 for none); a 2x2x2 box aggregates to
    !! cell crs_agg(i) of the next level, and the last level (at most mg_bottom_max cells) is factored densely in crs_l. With one
    !! level the bottom solve is exact; past mg_bottom_max cells it is one V-cycle (an exact solve costs every rank more the more
    !! ranks there are), unless crs_exact keeps the bottom whole and CG solves it (crs_cg)
    integer :: crs_n, crs_nl, crs_tot
    logical :: crs_exact  !< bottom solved exactly (proj_mg_bottom, or by rank count)
    logical :: crs_cg  !< exactly, by CG, as it is too large to factor densely
    real(wp), parameter :: crs_tol = 1.e-12_wp  !< CG bottom tolerance, relative to its right-hand side
    integer :: mg_kfull  !< levels between K-cycle steps when a solve reaches the bottom
    integer :: mg_kspace  !< levels between K-cycle steps in this solve: mg_kfull, or V-cycles once truncated (s_mg_omega)
    integer, dimension(mg_maxlev) :: crs_ln, crs_loff
    integer, dimension(3, mg_maxlev) :: crs_ld
    integer, allocatable, dimension(:) :: crs_cnt, crs_disp, crs_agg
    integer, allocatable, dimension(:,:) :: crs_sz, crs_aj, crs_co
    real(wp), allocatable, dimension(:,:) :: crs_l, crs_ak
    real(wp), allocatable, dimension(:) :: crs_ad, crs_x, crs_b
    !> Boundary code of each domain face, low (1) and high (2) side, copied from bc_type: walls and boundary patches per face
    integer, allocatable, dimension(:,:,:) :: fbx, fby, fbz
    $:GPU_DECLARE(create='[fbx, fby, fbz]')
    !> An extrapolation (-3) face is a pressure outlet held at bc_[x,y,z]%pres_out where that is set: its conductance couples the
    !! boundary cell to a ghost at pout, so flow leaves where the solve would otherwise see a wall
    logical, dimension(3)  :: has_pout
    real(wp), dimension(3) :: pout
    logical, dimension(3)  :: seam_lo, seam_hi    !< that side couples to another rank or, periodically, to this one
    logical                :: faces_ready         !< uf has been seeded from the cell velocities
    logical                :: wb_st               !< well-balanced surface tension
    integer                :: gk0, gk1, gl0, gl1  !< y and z extents including one ghost layer where those directions exist

contains

    impure subroutine s_initialize_projection_module()

        integer  :: lv, tot, kl, ku
        real(wp) :: rlev

#ifdef MFC_MIXED_PRECISION
        call s_mpi_abort('proj_method needs stp = wp; mixed precision is not supported')
#endif

        ! Automatic choices by rank count; the thresholds are machine-dependent (docs: "Projection method iterative solve tuning")
        mg_kfull = proj_mg_kcycle
        if (mg_kfull < 0) mg_kfull = merge(2, 0, num_procs >= proj_mg_k_ranks)
        mg_kspace = mg_kfull
        ! A K-cycle reaches the bottom several times per cycle, which multiplies the cost of an exact solve there
        crs_exact = proj_mg_bottom == 1 .or. (proj_mg_bottom == -1 .and. mg_kfull == 0 .and. num_procs <= proj_mg_cg_ranks)

        @:ALLOCATE(uf(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims), uf0(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims))
        @:ALLOCATE(divu(0:m, 0:n, 0:p), rhs_p(0:m, 0:n, 0:p), p_stage(0:m, 0:n, 0:p), p_step0(0:m, 0:n, 0:p))
        @:ALLOCATE(pflx(-1:m + 1, -1:n + 1, -1:p + 1), rhoc(-1:m + 1, -1:n + 1, -1:p + 1))
        @:ALLOCATE(dcoef(0:m, 0:n, 0:p), bvec(0:m, 0:n, 0:p))
        @:ALLOCATE(xs(0:m, 0:n, 0:p), rs(0:m, 0:n, 0:p), zs(0:m, 0:n, 0:p), qs(0:m, 0:n, 0:p))
        @:ALLOCATE(pk(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        pk_sf(1)%sf => pk
        $:GPU_ENTER_DATA(copyin='[pk_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[pk_sf(1)%sf]')

        ! Halve every direction (rounding down) until each rank holds one cell; all ranks take the global level count, so a rank
        ! that reaches one cell early keeps it, and neighbors, which share their tangential sizes, stay aligned. The hierarchy then
        ! stops at the first level small enough for the dense bottom solve
        mg_gx = 1; mg_gy = merge(1, 0, n > 0); mg_gz = merge(1, 0, p > 0)
        mg_nx(1) = m + 1; mg_ny(1) = n + 1; mg_nz(1) = p + 1
        lv = 1
        do while (max(mg_nx(lv), mg_ny(lv), mg_nz(lv)) > 1)
            mg_nx(lv + 1) = max(mg_nx(lv)/2, 1); mg_ny(lv + 1) = max(mg_ny(lv)/2, 1); mg_nz(lv + 1) = max(mg_nz(lv)/2, 1)
            lv = lv + 1
        end do
        call s_mpi_allreduce_max(real(lv, wp), rlev)
        mg_nlev = nint(rlev)
        do lv = lv + 1, mg_nlev
            mg_nx(lv) = 1; mg_ny(lv) = 1; mg_nz(lv) = 1
        end do
        do lv = 1, mg_nlev
            call s_mpi_allreduce_max(real(mg_nx(lv)*mg_ny(lv)*mg_nz(lv), wp), rlev)
            if (rlev*num_procs <= max(mg_bottom_max, num_procs)) exit
        end do
        mg_nlev = min(lv, mg_nlev)
        tot = 0
        do lv = 1, mg_nlev
            mg_off(lv) = tot
            tot = tot + (mg_nx(lv) + 2*mg_gx)*(mg_ny(lv) + 2*mg_gy)*(mg_nz(lv) + 2*mg_gz)
        end do
        mg_poff = tot
        @:ALLOCATE(mg_d(tot), mg_kx(tot), mg_ky(tot), mg_kz(tot), mg_f(tot), mg_r(tot))
        @:ALLOCATE(mg_e(tot + (mg_nx(1) + 2*mg_gx)*(mg_ny(1) + 2*mg_gy)*(mg_nz(1) + 2*mg_gz)))
        kl = mg_off(min(2, mg_nlev)) + 1; ku = merge(mg_poff, kl - 1, mg_kfull > 0)
        @:ALLOCATE(mg_kr(kl:ku), mg_kc(kl:ku), mg_kv(kl:ku))
        tot = 2*((n + 1)*(p + 1) + (m + 1)*(p + 1) + (m + 1)*(n + 1))
        @:ALLOCATE(mg_sbuf(tot), mg_rbuf(tot))
        @:PIN_HOST(mg_sbuf, mg_rbuf)

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            seam_lo(${D}$) = bc_${XYZ}$%beg >= 0 .or. bc_${XYZ}$%beg == BC_PERIODIC
            seam_hi(${D}$) = bc_${XYZ}$%end >= 0 .or. bc_${XYZ}$%end == BC_PERIODIC
            mg_nbr(2*${D}$ - 1) = f_seam_rank(bc_${XYZ}$%beg, num_dims >= ${D}$)
            mg_nbr(2*${D}$) = f_seam_rank(bc_${XYZ}$%end, num_dims >= ${D}$)
        #:endfor
        call s_mg_exchange_layout()
        call s_mg_bottom_layout()

        #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            has_pout(${D}$) = .not. f_is_default(bc_${XYZ}$%pres_out) .and. num_dims >= ${D}$
            pout(${D}$) = merge(bc_${XYZ}$%pres_out, 0._wp, has_pout(${D}$))
        #:endfor
        @:ALLOCATE(fbx(0:n, 0:p, 1:2), fby(0:merge(m, 0, n > 0), 0:p, 1:2), fbz(0:merge(m, 0, p > 0), 0:merge(n, 0, p > 0), 1:2))
        faces_ready = .false.
        gk0 = merge(-1, 0, n > 0); gk1 = merge(n + 1, n, n > 0)
        gl0 = merge(-1, 0, p > 0); gl1 = merge(p + 1, p, p > 0)

        wb_st = surface_tension .and. surface_tension_model == surface_tension_model_well_balanced
        if (wb_st) then
            @:ALLOCATE(kap(-1:m + 1, gk0:gk1, gl0:gl1, 1:2))
            @:ALLOCATE(gnd(-2:m + 2, 2*gk0:gk1 - gk0, 2*gl0:gl1 - gl0, 0:3))
        else
            @:ALLOCATE(kap(0:0, 0:0, 0:0, 1:2), gnd(0:0, 0:0, 0:0, 0:3))
        end if

    end subroutine s_initialize_projection_module

    !> A wall face: its normal velocity vanishes
    pure logical function f_wall_face(code)

        $:GPU_ROUTINE(function_name='f_wall_face', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: code

        f_wall_face = code == BC_REFLECTIVE .or. code == BC_SLIP_WALL .or. code == BC_NO_SLIP_WALL

    end function f_wall_face

    !> Copy the boundary code of each domain face from bc_type, on the host where pre-process left it
    impure subroutine s_projection_face_codes(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: e

        do e = 1, 2
            fbx(:,:,e) = int(bc_type(1, e)%sf(0,:,:))
            if (num_dims > 1) fby(:,:,e) = int(bc_type(2, e)%sf(0:m,0,:))
            if (num_dims > 2) fbz(:,:,e) = int(bc_type(3, e)%sf(0:m,0:n,0))
        end do
        $:GPU_UPDATE(device='[fbx, fby, fbz]')

    end subroutine s_projection_face_codes

    !> Face velocity from the average of the cell velocities either side (primitive velocities, with ghosts)
    subroutine s_projection_init_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: j, k, l

        #:for D, IP1, LB, KB, JB in [(1, 'j + 1, k, l', 0, 0, -1), (2, 'j, k + 1, l', 0, -1, 0), (3, 'j, k, l + 1', -1, 0, 0)]
            if (num_dims >= ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            uf(j, k, l, ${D}$) = 0.5_wp*(real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), &
                               & wp) + real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${IP1}$), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

    end subroutine s_projection_init_faces

    !> A Dirichlet (-17) face, a whole boundary or a patch such as a nozzle set into a wall, carries the normal velocity of its
    !! ghost cell, which s_dirichlet sets from the boundary buffer and its inflow ramp: prescribed, it is neither predicted nor
    !! pressure-corrected (the pressure ghost is a zero-gradient copy) and the solve sees it as a source of volume
    subroutine s_dirichlet_faces(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: j, k, l

        #:set UB = {'j': 'm', 'k': 'n', 'l': 'p'}
        ! The low face and its ghost cell share index -1; the high face is m (n, p), its ghost m + 1
        #:for D, NV, LO, HI, GH, FB, FI in [(1, 'j', '-1, k, l', 'm, k, l', 'm + 1, k, l', 'fbx', 'k, l'), &
            (2, 'k', 'j, -1, l', 'j, n, l', 'j, n + 1, l', 'fby', 'j, l'), (3, 'l', 'j, k, -1', 'j, k, p', 'j, k, p + 1', 'fbz', &
             & 'j, k')]
            #:set TV = [v for v in ['l', 'k', 'j'] if v != NV]
            if (num_dims >= ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do ${TV[0]}$ = 0, ${UB[TV[0]]}$
                    do ${TV[1]}$ = 0, ${UB[TV[1]]}$
                        if (${FB}$(${FI}$, 1) == BC_DIRICHLET) uf(${LO}$, &
                            & ${D}$) = real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${LO}$), wp)
                        if (${FB}$(${FI}$, 2) == BC_DIRICHLET) uf(${HI}$, &
                            & ${D}$) = real(q_prim_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${GH}$), wp)
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_dirichlet_faces

    !> Normal velocity vanishes on solid walls
    subroutine s_zero_wall_faces()

        integer :: j, k, l

        #:set UB = {'j': 'm', 'k': 'n', 'l': 'p'}
        #:for D, NV, LO, HI, FB, FI in [(1, 'j', '-1, k, l', 'm, k, l', 'fbx', 'k, l'), (2, 'k', 'j, -1, l', 'j, n, l', 'fby', 'j, l'), &
            (3, 'l', 'j, k, -1', 'j, k, p', 'fbz', 'j, k')]
            #:set TV = [v for v in ['l', 'k', 'j'] if v != NV]
            if (num_dims >= ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do ${TV[0]}$ = 0, ${UB[TV[0]]}$
                    do ${TV[1]}$ = 0, ${UB[TV[1]]}$
                        if (f_wall_face(${FB}$(${FI}$, 1))) uf(${LO}$, ${D}$) = 0._wp
                        if (f_wall_face(${FB}$(${FI}$, 2))) uf(${HI}$, ${D}$) = 0._wp
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_zero_wall_faces

    !> Advective right-hand side of one direction sweep. Every quantity is carried by the projected face velocity and upwinded on
    !! its sign; the momentum flux is the summed partial-density flux times the upwind velocity, so mass and momentum move with one
    !! operator. Energy is left at zero here and rebuilt from the equation of state after the pressure solve.
    subroutine s_projection_rhs(id, qfl_rs, qfr_rs, q_prim_vf, flux_vf, rhs_vf, bc_type)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(in)                                 :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_vf, rhs_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                          :: bc_type
        real(wp)                                                                            :: vf, a_up, ar_up, fm
        logical                                                                             :: up_l
        integer                                                                             :: i, j, k, l

        if (id == 1) then
            if (.not. faces_ready) then
                call s_projection_face_codes(bc_type)
                call s_projection_init_faces(q_prim_vf)
                faces_ready = .true.
            end if
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu(j, k, l) = f_div_uf(j, k, l)
                        p_stage(j, k, l) = real(q_prim_vf(eqn_idx%E)%sf(j, k, l), wp)
                        ! Transport sources q*div(u), which make alpha and p advect rather than compress
                        rhs_p(j, k, l) = p_stage(j, k, l)*divu(j, k, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, sys_size
                            rhs_vf(i)%sf(j, k, l) = 0._stp
                        end do
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            rhs_vf(i)%sf(j, k, l) = real(real(q_prim_vf(i)%sf(j, k, l), wp)*divu(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_dirichlet_faces(q_prim_vf)
        end if

        #:for D, SV, COORDS, JB, KB, LB, DXV in [(1, 'j', '{SI}, k, l', -1, 0, 0, 'dx'), &
            (2, 'k', 'j, {SI}, l', 0, -1, 0, 'dy'), (3, 'l', 'j, k, {SI}', 0, 0, -1, 'dz')]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${D}$) then
                ! Face fluxes. Left state of face j is the right edge of cell j, right state the left edge of cell j+1
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, vf, a_up, ar_up, fm, up_l]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            vf = uf(j, k, l, ${D}$)
                            up_l = vf >= 0._wp
                            fm = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                ! Reconstructed partial densities upwinded directly. Taking the phase density alpha_rho/alpha from
                                ! the upwind cell instead diverges where a phase is vanishing: both are round-off there, and their
                                ! ratio (seen at 6e10) times the face alpha flux injects mass
                                if (up_l) then
                                    a_up = qfl_rs(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    ar_up = qfl_rs(${SF('')}$, i)
                                else
                                    a_up = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                    ar_up = qfr_rs(${SF(' + 1')}$, i)
                                end if
                                flux_vf(eqn_idx%adv%beg + i - 1)%sf(${SF('')}$) = real(a_up*vf, stp)
                                flux_vf(i)%sf(${SF('')}$) = real(ar_up*vf, stp)
                                fm = fm + ar_up*vf
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                if (up_l) then
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(${SF('')}$) = real(fm*qfl_rs(${SF('')}$, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                else
                                    flux_vf(eqn_idx%mom%beg + i - 1)%sf(${SF('')}$) = real(fm*qfr_rs(${SF(' + 1')}$, &
                                            & eqn_idx%mom%beg + i - 1), stp)
                                end if
                            end do
                            if (up_l) then
                                pflx(j, k, l) = vf*qfl_rs(${SF('')}$, eqn_idx%E)
                            else
                                pflx(j, k, l) = vf*qfr_rs(${SF(' + 1')}$, eqn_idx%E)
                            end if
                            ! The color function is not reconstructed; upwind its cell value
                            if (surface_tension) then
                                if (up_l) then
                                    flux_vf(eqn_idx%c)%sf(${SF('')}$) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(${SF('')}$), wp), stp)
                                else
                                    flux_vf(eqn_idx%c)%sf(${SF('')}$) = real(vf*real(q_prim_vf(eqn_idx%c)%sf(${SF(' + 1')}$), &
                                            & wp), stp)
                                end if
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, sys_size
                                if (i /= eqn_idx%E) then
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + real((real(flux_vf(i)%sf(${SF(' - 1')}$), &
                                           & wp) - real(flux_vf(i)%sf(j, k, l), wp))/${DXV}$(${SV}$), stp)
                                end if
                            end do
                            rhs_p(j, k, l) = rhs_p(j, k, l) + (pflx(${SF(' - 1')}$) - pflx(j, k, l))/${DXV}$(${SV}$)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_rhs

    !> Face data that the viscous and capillary source fluxes otherwise take from a Riemann solve: interface Reynolds numbers, the
    !! mean face velocity (read only by their energy terms, which the projection rebuilds from the EOS), and for surface tension the
    !! face velocity whose divergence makes the color function advect
    subroutine s_projection_face_props(id, qfl_rs, qfr_rs, flux_src_vf)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: qfl_rs, qfr_rs
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_src_vf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: al, ar
        #:else
            real(wp), dimension(num_fluids) :: al, ar
        #:endif
        real(wp), dimension(2) :: re_l, re_r
        integer                :: i, j, k, l, rs1, rs2

        rs1 = Re_size(1); rs2 = Re_size(2)

        #:for D, SV, COORDS, JB, KB, LB in [(1, 'j', '{SI}, k, l', -1, 0, 0), (2, 'k', 'j, {SI}, l', 0, -1, 0), &
            (3, 'l', 'j, k, {SI}', 0, 0, -1)]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, al, ar, re_l, re_r]', firstprivate='[rs1, rs2]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    al(i) = qfl_rs(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    ar(i) = qfr_rs(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end do
                                call s_compute_interface_reynolds(al, re_l, rs1, rs2)
                                call s_compute_interface_reynolds(ar, re_r, rs1, rs2)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rsx_vf(j, k, l, i) = 2._wp/(1._wp/re_l(i) + 1._wp/re_r(i))
                                end do
                            end if
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_src_rsx_vf(j, k, l, i) = 0.5_wp*(qfl_rs(${SF('')}$, &
                                               & eqn_idx%mom%beg + i - 1) + qfr_rs(${SF(' + 1')}$, eqn_idx%mom%beg + i - 1))
                            end do
                            if (surface_tension) flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = real(uf(j, k, l, ${D}$), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_face_props

    !> Conduction heating of one direction sweep from the energy source flux fe (-k dT/dn on faces). At fixed alpha and phase
    !! densities d(rho e) = Gamma dp, so the pressure transport rate gains -div(q)/Gamma, and the energy rebuilt from p carries it
    subroutine s_projection_heat(id, q_prim_vf, fe)

        integer, intent(in)                                 :: id
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(in)                      :: fe
        real(wp)                                            :: rho, gam, pinf, qv

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: ar, al
        #:else
            real(wp), dimension(num_fluids) :: ar, al
        #:endif
        integer :: i, j, k, l

        #:for D, SV, COORDS, DXV in [(1, 'j', '{SI}, k, l', 'dx'), (2, 'k', 'j, {SI}, l', 'dy'), (3, 'l', 'j, k, {SI}', 'dz')]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${D}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, ar, al]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                ar(i) = real(q_prim_vf(i)%sf(j, k, l), wp)
                                al(i) = real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                            end do
                            call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                            rhs_p(j, k, l) = rhs_p(j, k, l) + (real(fe%sf(${SF(' - 1')}$), wp) - real(fe%sf(j, k, l), &
                                  & wp))/(${DXV}$(${SV}$)*gam)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_heat

    !> Divergence of the face velocity in cell (j, k, l): the one operator the transport sources and the pressure equation share
    function f_div_uf(j, k, l) result(dv)

        $:GPU_ROUTINE(function_name='f_div_uf', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: j, k, l
        real(wp)            :: dv

        dv = (uf(j, k, l, 1) - uf(j - 1, k, l, 1))/dx(j)
        if (num_dims > 1) dv = dv + (uf(j, k, l, 2) - uf(j, k - 1, l, 2))/dy(k)
        if (num_dims > 2) dv = dv + (uf(j, k, l, 3) - uf(j, k, l - 1, 3))/dz(l)

    end function f_div_uf

    !> Face conductance A_f/(rho_f*d_f) with the arithmetic face density: the inertia of a face volume straddling an interface,
    !! which a heavy phase resting on a light one needs to stay at rest
    pure function f_cond(ra, rb, area, dist) result(kf)

        $:GPU_ROUTINE(function_name='f_cond', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ra, rb, area, dist
        real(wp)             :: kf

        kf = 2._wp*area/(max(ra + rb, sgm_eps)*dist)

    end function f_cond

    !> Pressure correction of a face velocity, tau*grad(p)/rho_f with the operator's own conductance (area 1), so div(uf) matches
    !! the solved pressure exactly
    pure function f_face_gf(ra, rb, dist, pa, pb, tau) result(gf)

        $:GPU_ROUTINE(function_name='f_face_gf', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ra, rb, dist, pa, pb, tau
        real(wp)             :: gf

        gf = tau*f_cond(ra, rb, 1._wp, dist)*(pb - pa)

    end function f_face_gf

    !> Well-balanced (Brackbill CSF) capillary acceleration of a face: sigma*kappa_f*(c_b - c_a)/(d_f*rho_f), with rho_f the same
    !! arithmetic face density as the pressure operator, so a constant curvature is balanced exactly by a pressure jump. kappa_f is
    !! the |grad c|-weighted mean of the adjacent cells' curvature (weights wa, wb)
    pure function f_capillary_accel(ka, kb, wa, wb, ca, cb, ra, rb, dist) result(acc)

        $:GPU_ROUTINE(function_name='f_capillary_accel', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ka, kb, wa, wb, ca, cb, ra, rb, dist
        real(wp)             :: acc

        acc = 0._wp
        if (wa + wb > 0._wp) acc = 2._wp*sigma*(wa*ka + wb*kb)/(wa + wb)*(cb - ca)/(dist*max(ra + rb, sgm_eps))

    end function f_capillary_accel

    !> Interface normal component, grad_d(c)/|grad(c)|; zero outside the interface band
    pure function f_normal(gd, g) result(nd)

        $:GPU_ROUTINE(function_name='f_normal', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: gd, g
        real(wp)             :: nd

        nd = 0._wp
        if (g > capillary_cutoff) nd = gd/g

    end function f_normal

    !> Curvature kappa = -div(n) in the interface band, with n = grad(alpha_1)/|grad(alpha_1)|. The CSF force sigma*kappa*grad(c) is
    !! unchanged under c -> 1 - c, so the conservatively transported, compression-sharpened volume fraction serves as the indicator;
    !! the color function, upwinded for the stress-tensor model, smears and loses its interior value over time
    subroutine s_compute_curvature(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp)                                            :: kv
        integer                                             :: j, k, l, k0, k1, l0, l1, dk, dl, ia

        k0 = gk0; k1 = gk1; l0 = gl0; l1 = gl1
        dk = -k0; dl = -l0
        ia = eqn_idx%adv%beg
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = l0 - dl, l1 + dl
            do k = k0 - dk, k1 + dk
                do j = -2, m + 2
                    gnd(j, k, l, 1) = real(q_cons_vf(ia)%sf(j + 1, k, l) - q_cons_vf(ia)%sf(j - 1, k, l), &
                        & wp)/(x_cc(j + 1) - x_cc(j - 1))
                    gnd(j, k, l, 2) = 0._wp; gnd(j, k, l, 3) = 0._wp
                    if (num_dims > 1) gnd(j, k, l, 2) = real(q_cons_vf(ia)%sf(j, k + 1, l) - q_cons_vf(ia)%sf(j, k - 1, l), &
                        & wp)/(y_cc(k + 1) - y_cc(k - 1))
                    if (num_dims > 2) gnd(j, k, l, 3) = real(q_cons_vf(ia)%sf(j, k, l + 1) - q_cons_vf(ia)%sf(j, k, l - 1), &
                        & wp)/(z_cc(l + 1) - z_cc(l - 1))
                    gnd(j, k, l, 0) = sqrt(gnd(j, k, l, 1)**2 + gnd(j, k, l, 2)**2 + gnd(j, k, l, 3)**2)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, kv]')
        do l = l0, l1
            do k = k0, k1
                do j = -1, m + 1
                    kap(j, k, l, 1) = 0._wp
                    kap(j, k, l, 2) = 0._wp
                    if (gnd(j, k, l, 0) > capillary_cutoff) then
                        #:set NRM = lambda d, idx: f"f_normal(gnd({idx}, {d}), gnd({idx}, 0))"
                        kv = -(${NRM(1, 'j + 1, k, l')}$ - ${NRM(1, 'j - 1, k, l')}$)/(x_cc(j + 1) - x_cc(j - 1))
                        if (num_dims > 1) kv = kv - (${NRM(2, 'j, k + 1, l')}$ - ${NRM(2, 'j, k - 1, l')}$)/(y_cc(k + 1) - y_cc(k &
                            & - 1))
                        if (num_dims > 2) kv = kv - (${NRM(3, 'j, k, l + 1')}$ - ${NRM(3, 'j, k, l - 1')}$)/(z_cc(l + 1) - z_cc(l &
                            & - 1))
                        kap(j, k, l, 1) = kv
                        kap(j, k, l, 2) = gnd(j, k, l, 0)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_curvature

    !> Pressure solve and correction on the blended (star) state of one RK stage
    impure subroutine s_projection_apply(q_cons_vf, bc_type, pb_in, mv_in, q_T_sf, rkc1, rkc2, rkc3, rkc4, stage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), intent(inout) :: q_T_sf
        real(wp), intent(in) :: rkc1, rkc2, rkc3, rkc4
        integer, intent(in) :: stage
        real(wp) :: tau, rho, gam, pinf, qv, rc2, dv
        real(wp) :: vol, ke, ga, gf
        real(wp), dimension(3) :: acc
        logical :: wbl, fix

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: ar, al
        #:else
            real(wp), dimension(num_fluids) :: ar, al
        #:endif
        integer :: i, j, k, l, k0, k1, l0, l1

        tau = rkc3*dt/rkc4
        k0 = gk0; k1 = gk1; l0 = gl0; l1 = gl1

        call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)

        acc = 0._wp
        if (bodyForces) then
            call s_compute_acceleration(mytime)
            #:for D, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (bf_${XYZ}$) acc(${D}$) = accel_bf(${D}$)
            #:endfor
        end if

        ! Keep the face velocity this stage was transported with, for s_mass_lag below
        if (proj_lag_corr > 0) then
            $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
            do i = 1, num_dims
                do l = l0, l1
                    do k = k0, k1
                        do j = -1, m + 1
                            uf0(j, k, l, i) = uf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Star density with its ghosts, and the face predictor from the star cell velocities
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
        do l = l0, l1
            do k = k0, k1
                do j = -1, m + 1
                    rhoc(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rhoc(j, k, l) = rhoc(j, k, l) + real(q_cons_vf(i)%sf(j, k, l), wp)
                    end do
                    rhoc(j, k, l) = max(rhoc(j, k, l), sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (wb_st) call s_compute_curvature(q_cons_vf)
        wbl = wb_st

        ! Body forces and well-balanced surface tension enter on faces, where they meet the pressure gradient that balances them.
        ! A Dirichlet face keeps its prescribed velocity
        #:for D, DXV, SV, UB, IP1, LB, KB, JB, FB, FI in [(1, 'dx', 'j', 'm', 'j + 1, k, l', 0, 0, -1, 'fbx', 'k, l'), &
            (2, 'dy', 'k', 'n', 'j, k + 1, l', 0, -1, 0, 'fby', 'j, l'), (3, 'dz', 'l', 'p', 'j, k, l + 1', -1, 0, 0, 'fbz', &
             & 'j, k')]
            if (num_dims >= ${D}$) then
                ga = tau*acc(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, fix]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            fix = .false.
                            if (${SV}$ == -1) fix = ${FB}$(${FI}$, 1) == BC_DIRICHLET
                            if (${SV}$ == ${UB}$) fix = ${FB}$(${FI}$, 2) == BC_DIRICHLET
                            if (.not. fix) then
                                uf(j, k, l, ${D}$) = 0.5_wp*(real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), wp)/rhoc(j, &
                                   & k, l) + real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(${IP1}$), wp)/rhoc(${IP1}$)) + ga
                                if (wbl) uf(j, k, l, ${D}$) = uf(j, k, l, ${D}$) + tau*f_capillary_accel(kap(j, k, l, 1), &
                                    & kap(${IP1}$, 1), kap(j, k, l, 2), kap(${IP1}$, 2), real(q_cons_vf(eqn_idx%adv%beg)%sf(j, k, &
                                    & l), wp), real(q_cons_vf(eqn_idx%adv%beg)%sf(${IP1}$), wp), rhoc(j, k, l), rhoc(${IP1}$), &
                                    & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)))
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        ! SPD system D_c p + sum_f K_f (p - p_nb) = b, the Helmholtz row scaled by V_c/(rho c^2 tau^2)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, rc2, dv, vol, ar, al]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (stage == 1) p_step0(j, k, l) = p_stage(j, k, l)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ar(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        al(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                    ! Allaire's model advects alpha, so gamma_mix and pi_inf_mix are advected and Dp/Dt = -K div(u) with K the
                    ! mixture bulk modulus (not Wood's, which belongs to the Kapila model)
                    rc2 = max(f_bulk_modulus(p_stage(j, k, l), gam, pinf), sgm_eps)
                    dv = f_div_uf(j, k, l)
                    vol = dx(j)
                    if (num_dims > 1) vol = vol*dy(k)
                    if (num_dims > 2) vol = vol*dz(l)
                    dcoef(j, k, l) = vol/(rc2*tau*tau)
                    bvec(j, k, l) = dcoef(j, k, l)*((rkc1*p_stage(j, k, l) + rkc2*p_step0(j, k, l) + rkc3*dt*rhs_p(j, k, &
                         & l))/rkc4 - rc2*tau*dv)
                    xs(j, k, l) = p_stage(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (stage == 1) then
            proj_pcg_iters = 0
            proj_pcg_res = 0._wp
        end if
        call nvtxStartRange("PROJ-PCG")
        call s_pcg_solve(bc_type)
        call nvtxEndRange

        if (proj_lag_corr == 2) call s_mass_lag(q_cons_vf, tau, .true.)

        ! Face correction with the operator's own conductance (area 1), so div(uf) matches the solved pressure exactly. Cells take
        ! the mean of their faces' net acceleration (body force less pressure gradient, zero on walls): a hydrostatic balance on
        ! the faces then leaves the cells at rest too, and for uniform density this is the centered pressure gradient
        #:for D, DXV, SV, UB, IP1, IM1, LB, KB, JB, FB, FI in [(1, 'dx', 'j', 'm', 'j + 1, k, l', 'j - 1, k, l', 0, 0, -1, 'fbx', &
            & 'k, l'), (2, 'dy', 'k', 'n', 'j, k + 1, l', 'j, k - 1, l', 0, -1, 0, 'fby', 'j, l'), (3, 'dz', 'l', 'p', 'j, k, l + 1', &
            & 'j, k, l - 1', -1, 0, 0, 'fbz', 'j, k')]
            if (num_dims >= ${D}$) then
                ga = tau*acc(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, gf]')
                do l = ${LB}$, p
                    do k = ${KB}$, n
                        do j = ${JB}$, m
                            gf = f_face_gf(rhoc(j, k, l), rhoc(${IP1}$), 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)), &
                                           & real(pk(j, k, l), wp), real(pk(${IP1}$), wp), tau)
                            uf(j, k, l, ${D}$) = uf(j, k, l, ${D}$) - gf
                            pflx(j, k, l) = ga - gf
                            if (wbl) pflx(j, k, l) = pflx(j, k, l) + tau*f_capillary_accel(kap(j, k, l, 1), kap(${IP1}$, 1), &
                                & kap(j, k, l, 2), kap(${IP1}$, 2), real(q_cons_vf(eqn_idx%adv%beg)%sf(j, k, l), wp), &
                                & real(q_cons_vf(eqn_idx%adv%beg)%sf(${IP1}$), wp), rhoc(j, k, l), rhoc(${IP1}$), &
                                & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)))
                            if (${SV}$ == -1) then
                                if (f_wall_face(${FB}$(${FI}$, 1)) .or. ${FB}$(${FI}$, 1) == BC_DIRICHLET) pflx(j, k, l) = 0._wp
                            else if (${SV}$ == ${UB}$) then
                                if (f_wall_face(${FB}$(${FI}$, 2)) .or. ${FB}$(${FI}$, 2) == BC_DIRICHLET) pflx(j, k, l) = 0._wp
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, &
                                      & l) = real(real(q_cons_vf(eqn_idx%mom%beg + ${D}$ - 1)%sf(j, k, l), wp) + rhoc(j, k, &
                                      & l)*0.5_wp*(pflx(${IM1}$) + pflx(j, k, l)), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor
        call s_zero_wall_faces()

        if (proj_lag_corr > 0) call s_mass_lag(q_cons_vf, tau, .false.)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho, gam, pinf, qv, ke, ar, al]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ar(i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                        al(i) = real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                    end do
                    call s_compute_mixture_coefficients(ar, al, rho, gam, pinf, qv)
                    ke = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        ke = ke + 0.5_wp*real(q_cons_vf(i)%sf(j, k, l), wp)*(real(q_cons_vf(i)%sf(j, k, l), wp)/rho)
                    end do
                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = real(gam*real(pk(j, k, l), wp) + pinf + qv + ke, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_apply

    !> The stage transported the partial densities, volume fractions and momentum with uf0, while the pressure solve compressed p
    !! with the corrected face velocity: left alone, density and pressure disagree by tau*rho*div(u_f - uf0) every stage, an entropy
    !! error that accumulates wherever that divergence persists (small scales, walls, low acoustic CFL), and the moved mass would
    !! not carry its momentum. Move them by the difference too, conservatively (alpha in its advective form), upwinded on the
    !! corrected face velocity. Momentum (mom) is moved before the correction, which forms that velocity here as the correction
    !! will, so every upwind value, ghosts included, is the pre-correction state on both sides of a rank seam; the densities and
    !! alpha after it, where uf is final and the correction has read the transported alpha. Each component is summed into rhs_p,
    !! free once the pressure system is built, then applied.
    subroutine s_mass_lag(q_cons_vf, tau, mom)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: tau
        logical, intent(in)                                    :: mom
        real(wp)                                               :: c, dd, du, vn, tl
        integer                                                :: iq, q, j, k, l, nalpha, i1, i2
        logical                                                :: pre

        tl = tau; pre = mom  ! a dummy may alias a host variable, which a device kernel must not reference
        ! A lone fluid's alpha is 1, which the difference leaves unchanged
        nalpha = merge(0, num_fluids, num_fluids == 1)
        i1 = merge(num_fluids + nalpha + 1, 1, mom)
        i2 = merge(num_fluids + nalpha + num_dims, num_fluids + nalpha, mom)
        do iq = i1, i2
            if (iq <= num_fluids) then
                q = iq
            else if (iq <= num_fluids + nalpha) then
                q = eqn_idx%adv%beg + iq - num_fluids - 1
            else
                q = eqn_idx%mom%beg + iq - num_fluids - nalpha - 1
            end if
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, c, dd, du, vn]', firstprivate='[q, tl, pre]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c = 0._wp; dd = 0._wp
                        #:for D, DXV, SV, UB, IP1, IM1, FB, FI in [(1, 'dx', 'j', 'm', 'j + 1, k, l', 'j - 1, k, l', 'fbx', 'k, l'), &
                            (2, 'dy', 'k', 'n', 'j, k + 1, l', 'j, k - 1, l', 'fby', 'j, l'), &
                            (3, 'dz', 'l', 'p', 'j, k, l + 1', 'j, k, l - 1', 'fbz', 'j, k')]
                            if (num_dims >= ${D}$) then
                                ! High face: between this cell and the next
                                vn = uf(j, k, l, ${D}$)
                                ! The correction
                                if (pre) vn = vn - f_face_gf(rhoc(j, k, l), rhoc(${IP1}$), &
                                    & 0.5_wp*(${DXV}$(${SV}$) + ${DXV}$(${SV}$ + 1)), real(pk(j, k, l), wp), real(pk(${IP1}$), &
                                    & wp), tl)
                                if (${SV}$ == ${UB}$) then
                                    if (f_wall_face(${FB}$(${FI}$, 2))) vn = 0._wp
                                end if
                                du = vn - uf0(j, k, l, ${D}$)
                                c = c + du*real(merge(q_cons_vf(q)%sf(j, k, l), q_cons_vf(q)%sf(${IP1}$), vn >= 0._wp), &
                                                & wp)/${DXV}$(${SV}$)
                                dd = dd + du/${DXV}$(${SV}$)
                                ! Low face: between the previous cell and this one
                                vn = uf(${IM1}$, ${D}$)
                                if (pre) vn = vn - f_face_gf(rhoc(${IM1}$), rhoc(j, k, l), &
                                    & 0.5_wp*(${DXV}$(${SV}$ - 1) + ${DXV}$(${SV}$)), real(pk(${IM1}$), wp), real(pk(j, k, l), &
                                    & wp), tl)
                                if (${SV}$ == 0) then
                                    if (f_wall_face(${FB}$(${FI}$, 1))) vn = 0._wp
                                end if
                                du = vn - uf0(${IM1}$, ${D}$)
                                c = c - du*real(merge(q_cons_vf(q)%sf(${IM1}$), q_cons_vf(q)%sf(j, k, l), vn >= 0._wp), &
                                                & wp)/${DXV}$(${SV}$)
                                dd = dd - du/${DXV}$(${SV}$)
                            end if
                        #:endfor
                        if (iq > num_fluids .and. iq <= num_fluids + nalpha) c = c - real(q_cons_vf(q)%sf(j, k, l), wp)*dd
                        rhs_p(j, k, l) = c
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', firstprivate='[q, tl]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(q)%sf(j, k, l) = real(real(q_cons_vf(q)%sf(j, k, l), wp) - tl*rhs_p(j, k, l), stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_mass_lag

    !> PCG on the SPD pressure system, preconditioned by one multigrid V-cycle. The solution is left in pk with filled ghosts.
    impure subroutine s_pcg_solve(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: bnorm, rnorm, rtol, rz, rz_new, zq, dq, alpha, beta, r0
        integer                                                    :: it, j, k, l, idx, off, ex, ey, gx, gy, gz, sh

        call nvtxStartRange("PROJ-MG-BUILD")
        call s_mg_build()
        call nvtxEndRange
        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy; sh = mg_poff

        ! The search direction's wall ghosts stay zero: those faces carry no conductance, and the ghosts must not hold garbage
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_e(${MG_IX('j', 'k', 'l')}$ + sh) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (any(has_pout)) call s_pout_rhs()
        call s_set_dir(xs)
        call s_apply_operator()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rs(j, k, l) = bvec(j, k, l) - qs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Relative to the initial residual, i.e. to the pressure change being solved for; relative to b the tolerance would
        ! depend on the ambient pressure, which b carries in full
        bnorm = sqrt(f_dot(bvec, bvec))
        rnorm = sqrt(f_dot(rs, rs))
        rtol = max(proj_tol*rnorm, res_floor*bnorm)
        r0 = rnorm

        if (rnorm > rtol) then
            call nvtxStartRange("PROJ-VCYCLE")
            call s_mg_precond()
            call nvtxEndRange
            call s_set_dir(zs)
            rz = f_dot(rs, zs)

            do it = 1, proj_max_iters
                proj_pcg_iters = proj_pcg_iters + 1
                call s_apply_operator()
                dq = f_dot_dir(qs)
                alpha = rz/dq
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            idx = ${MG_IX('j', 'k', 'l')}$ + sh
                            xs(j, k, l) = xs(j, k, l) + alpha*mg_e(idx)
                            rs(j, k, l) = rs(j, k, l) - alpha*qs(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                rnorm = sqrt(f_dot(rs, rs))
                if (rnorm <= rtol) exit

                call nvtxStartRange("PROJ-VCYCLE")
                call s_mg_precond()
                call nvtxEndRange
                ! Flexible (Polak-Ribiere) beta, z.(r - r_old)/rz with r_old - r = alpha*q, as a K-cycle preconditioner varies
                call s_dot_rz_zq(rz_new, zq)
                beta = -alpha*zq/rz
                rz = rz_new
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            idx = ${MG_IX('j', 'k', 'l')}$ + sh
                            mg_e(idx) = zs(j, k, l) + beta*mg_e(idx)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            ! Inside the branch: a solve already below its floor takes no iterations and reports none
            proj_pcg_res = max(proj_pcg_res, rnorm/r0)
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    pk(j, k, l) = real(xs(j, k, l), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_populate_F_igr_buffers(bc_type, pk_sf)
        if (any(has_pout)) call s_pout_ghosts()

    end subroutine s_pcg_solve

    !> A pressure outlet's fixed ghost enters the right-hand side as conductance times pout, its search-direction ghost being zero
    subroutine s_pout_rhs()

        integer  :: j, k, l, off, ex, ey, gx, gy, gz
        real(wp) :: po

        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy
        #:set UB = {'j': 'm', 'k': 'n', 'l': 'p'}
        #:for D, NV, KC, LO, HI, HF, FB, FI in [(1, 'j', 'mg_kx', '0, k, l', 'm, k, l', 'm + 1, k, l', 'fbx', 'k, l'), &
            (2, 'k', 'mg_ky', 'j, 0, l', 'j, n, l', 'j, n + 1, l', 'fby', 'j, l'), &
            (3, 'l', 'mg_kz', 'j, k, 0', 'j, k, p', 'j, k, p + 1', 'fbz', 'j, k')]
            #:set TV = [v for v in ['l', 'k', 'j'] if v != NV]
            if (has_pout(${D}$)) then
                po = pout(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[po]')
                do ${TV[0]}$ = 0, ${UB[TV[0]]}$
                    do ${TV[1]}$ = 0, ${UB[TV[1]]}$
                        if (${FB}$(${FI}$, &
                            & 1) == BC_GHOST_EXTRAP) bvec(${LO}$) = bvec(${LO}$) + ${KC}$(${MG_IX(*LO.split(', '))}$)*po
                        if (${FB}$(${FI}$, &
                            & 2) == BC_GHOST_EXTRAP) bvec(${HI}$) = bvec(${HI}$) + ${KC}$(${MG_IX(*HF.split(', '))}$)*po
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_pout_rhs

    !> A pressure outlet's ghost holds pout, so the face correction pushes flow out against it
    subroutine s_pout_ghosts()

        integer  :: j, k, l
        real(wp) :: po

        #:set UB = {'j': 'm', 'k': 'n', 'l': 'p'}
        #:for D, NV, LO, HI, FB, FI in [(1, 'j', '-1, k, l', 'm + 1, k, l', 'fbx', 'k, l'), (2, 'k', 'j, -1, l', 'j, n + 1, l', 'fby', &
            & 'j, l'), (3, 'l', 'j, k, -1', 'j, k, p + 1', 'fbz', 'j, k')]
            #:set TV = [v for v in ['l', 'k', 'j'] if v != NV]
            if (has_pout(${D}$)) then
                po = pout(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[po]')
                do ${TV[0]}$ = 0, ${UB[TV[0]]}$
                    do ${TV[1]}$ = 0, ${UB[TV[1]]}$
                        if (${FB}$(${FI}$, 1) == BC_GHOST_EXTRAP) pk(${LO}$) = real(po, stp)
                        if (${FB}$(${FI}$, 2) == BC_GHOST_EXTRAP) pk(${HI}$) = real(po, stp)
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_pout_ghosts

    !> Search direction (interior) = v
    subroutine s_set_dir(v)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: v
        integer                                      :: j, k, l, off, ex, ey, gx, gy, gz, sh

        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy; sh = mg_poff
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_e(${MG_IX('j', 'k', 'l')}$ + sh) = v(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_set_dir

    !> qs = A d for the search direction d, held in level 1's layout: A is level 1's operator, whose boundary faces carry no
    !! conductance unless they are rank or periodic seams, so one exchanged ghost layer is all it needs
    impure subroutine s_apply_operator()

        integer  :: j, k, l, idx, off, ex, ey, gx, gy, gz, sy, sz, sh
        real(wp) :: dg, nb

        call nvtxStartRange("PROJ-DIR-EXCHANGE")
        call s_mg_exchange(1, mg_poff)
        call nvtxEndRange
        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy; sy = ex; sz = ex*ey
        sh = mg_poff
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, idx, dg, nb]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    idx = ${MG_IX('j', 'k', 'l')}$
                    @:MG_ROW(sh)
                    qs(j, k, l) = dg*mg_e(idx + sh) - nb
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_operator

    !> Global r.z and q.z over the interior, in one reduction
    impure subroutine s_dot_rz_zq(rz, zq)

        real(wp), intent(out)     :: rz, zq
        real(wp), dimension(2, 1) :: sl, sg
        real(wp)                  :: a, b
        integer                   :: j, k, l

        a = 0._wp; b = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[a, b]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    a = a + rs(j, k, l)*zs(j, k, l)
                    b = b + qs(j, k, l)*zs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        sl(:,1) = [a, b]
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_vectors_sum(sl, sg, 2, 1)
        call nvtxEndRange
        rz = sg(1, 1); zq = sg(2, 1)

    end subroutine s_dot_rz_zq

    !> Global inner product over the interior
    impure function f_dot(a, b) result(res)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: a, b
        real(wp)                                     :: res, loc
        integer                                      :: j, k, l

        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + a(j, k, l)*b(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_sum(loc, res)
        call nvtxEndRange

    end function f_dot

    !> Global inner product of the search direction with a vector
    impure function f_dot_dir(a) result(res)

        real(wp), dimension(0:m,0:n,0:p), intent(in) :: a
        real(wp)                                     :: res, loc
        integer                                      :: j, k, l, off, ex, ey, gx, gy, gz, sh

        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy; sh = mg_poff
        loc = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[loc]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    loc = loc + mg_e(${MG_IX('j', 'k', 'l')}$ + sh)*a(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call nvtxStartRange("PROJ-ALLREDUCE")
        call s_mpi_allreduce_sum(loc, res)
        call nvtxEndRange

    end function f_dot_dir

    !> Level 1 from the fine system, with this rank's boundary faces dropped, then Galerkin coarsening for piecewise-constant
    !! aggregation: a coarse diagonal sums its children, a coarse face sums the fine faces lying on it. Exact for any coefficient
    !! jump, which rediscretizing an averaged density would not be.
    !> Neighbor across a domain side from its bc value: a rank, this rank for a periodic seam it owns alone, or -1
    pure function f_seam_rank(bc, active) result(rank)

        integer, intent(in) :: bc
        logical, intent(in) :: active
        integer             :: rank

        rank = -1
        if (.not. active) return
        if (bc >= 0) rank = bc
        if (bc == BC_PERIODIC) rank = proc_rank

    end function f_seam_rank

    !> Galerkin hierarchy. The fine level holds D_c and each cell's low-face conductances, seam faces included (the high seam face
    !! sits in the ghost cell past the last); a coarse cell sums its children's diagonals, and a coarse face the fine faces on it.
    !! The bottom level, one cell per rank, is gathered and factored.
    impure subroutine s_mg_build()

        integer  :: lv, nx, ny, nz, cnx, cny, cnz, coff, ii, jj, kk, a, b, c, idx, cidx, sy, sz, off, ex, ey, gx, gy, gz, cex, cey
        integer  :: fi, fj, fk, a0, a1, b0, b1, c0, c1
        real(wp) :: area, sd, skx, sky, skz
        logical  :: sl1, sh1, sl2, sh2, sl3, sh3, po1, po2, po3

        sl1 = seam_lo(1); sh1 = seam_hi(1); sl2 = seam_lo(2); sh2 = seam_hi(2); sl3 = seam_lo(3); sh3 = seam_hi(3)
        po1 = has_pout(1); po2 = has_pout(2); po3 = has_pout(3)
        gx = mg_gx; gy = mg_gy; gz = mg_gz
        off = mg_off(1); ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = -gz, p + gz
            do jj = -gy, n + gy
                do ii = -gx, m + gx
                    idx = ${MG_IX('ii', 'jj', 'kk')}$
                    mg_d(idx) = 0._wp; mg_kx(idx) = 0._wp; mg_ky(idx) = 0._wp; mg_kz(idx) = 0._wp
                    mg_e(idx) = 0._wp; mg_f(idx) = 0._wp; mg_r(idx) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, area]', firstprivate='[po1, po2, po3]')
        do kk = 0, p + gz
            do jj = 0, n + gy
                do ii = 0, m + gx
                    idx = ${MG_IX('ii', 'jj', 'kk')}$
                    if (ii <= m .and. jj <= n .and. kk <= p) mg_d(idx) = dcoef(ii, jj, kk)
                    if (jj <= n .and. kk <= p .and. ((ii > 0 .and. ii <= m) .or. (ii == 0 .and. (sl1 .or. (po1 .and. fbx(jj, kk, &
                        & 1) == BC_GHOST_EXTRAP))) .or. (ii == m + 1 .and. (sh1 .or. (po1 .and. fbx(jj, kk, &
                        & 2) == BC_GHOST_EXTRAP))))) then
                        area = 1._wp
                        if (num_dims > 1) area = dy(jj)
                        if (num_dims > 2) area = area*dz(kk)
                        mg_kx(idx) = f_cond(rhoc(ii, jj, kk), rhoc(ii - 1, jj, kk), area, 0.5_wp*(dx(ii - 1) + dx(ii)))
                    end if
                    if (num_dims > 1) then
                        if (ii <= m .and. kk <= p .and. ((jj > 0 .and. jj <= n) .or. (jj == 0 .and. (sl2 .or. (po2 .and. fby(ii, &
                            & kk, 1) == BC_GHOST_EXTRAP))) .or. (jj == n + 1 .and. (sh2 .or. (po2 .and. fby(ii, kk, &
                            & 2) == BC_GHOST_EXTRAP))))) then
                            area = dx(ii)
                            if (num_dims > 2) area = area*dz(kk)
                            mg_ky(idx) = f_cond(rhoc(ii, jj, kk), rhoc(ii, jj - 1, kk), area, 0.5_wp*(dy(jj - 1) + dy(jj)))
                        end if
                    end if
                    if (num_dims > 2) then
                        if (ii <= m .and. jj <= n .and. ((kk > 0 .and. kk <= p) .or. (kk == 0 .and. (sl3 .or. (po3 .and. fbz(ii, &
                            & jj, 1) == BC_GHOST_EXTRAP))) .or. (kk == p + 1 .and. (sh3 .or. (po3 .and. fbz(ii, jj, &
                            & 2) == BC_GHOST_EXTRAP))))) then
                            mg_kz(idx) = f_cond(rhoc(ii, jj, kk), rhoc(ii, jj, kk - 1), dx(ii)*dy(jj), 0.5_wp*(dz(kk - 1) + dz(kk)))
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do lv = 1, mg_nlev - 1
            nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv); ex = nx + 2*gx; ey = ny + 2*gy
            cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1); coff = mg_off(lv + 1)
            cex = cnx + 2*gx; cey = cny + 2*gy
            ! Coarse cells and faces, ghosts included (a ghost holds only its low face, the level's high boundary face)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, idx, cidx, sd, skx, sky, skz, fi, fj, fk, a0, a1, b0, &
                                & b1, c0, c1]')
            do kk = -gz, cnz - 1 + gz
                do jj = -gy, cny - 1 + gy
                    do ii = -gx, cnx - 1 + gx
                        sd = 0._wp; skx = 0._wp; sky = 0._wp; skz = 0._wp
                        if (ii >= 0 .and. jj >= 0 .and. kk >= 0) then
                            fi = merge(nx, 2*ii, ii == cnx); fj = merge(ny, 2*jj, jj == cny); fk = merge(nz, 2*kk, kk == cnz)
                            @:MG_CHILDREN(a0, a1, ii, cnx, nx)
                            @:MG_CHILDREN(b0, b1, jj, cny, ny)
                            @:MG_CHILDREN(c0, c1, kk, cnz, nz)
                            $:GPU_LOOP(parallelism='[seq]')
                            do c = c0, c1
                                $:GPU_LOOP(parallelism='[seq]')
                                do b = b0, b1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do a = a0, a1
                                        if (ii < cnx .and. jj < cny .and. kk < cnz) sd = sd + mg_d(${MG_IX('a', 'b', 'c')}$)
                                        if (a == a0 .and. jj < cny .and. kk < cnz) skx = skx + mg_kx(${MG_IX('fi', 'b', 'c')}$)
                                        if (b == b0 .and. ii < cnx .and. kk < cnz) sky = sky + mg_ky(${MG_IX('a', 'fj', 'c')}$)
                                        if (c == c0 .and. ii < cnx .and. jj < cny) skz = skz + mg_kz(${MG_IX('a', 'b', 'fk')}$)
                                    end do
                                end do
                            end do
                        end if
                        cidx = coff + ((kk + gz)*cey + jj + gy)*cex + ii + gx + 1
                        mg_d(cidx) = sd; mg_kx(cidx) = skx; mg_ky(cidx) = sky; mg_kz(cidx) = skz
                        mg_e(cidx) = 0._wp; mg_f(cidx) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_mg_omega()
        if (mg_blev == mg_nlev) call s_mg_bottom_build()

    end subroutine s_mg_build

    !> Coarse-correction scale of each level. Piecewise-constant aggregation under-corrects the Laplacian but represents the
    !! Helmholtz diagonal exactly, so the correction from level lv + 1 is scaled by 1 + (proj_mg_omega - 1)/(1 + r/mg_om_r0), r
    !! being that level's Helmholtz-to-Laplacian diagonal ratio: the full scale for Poisson-like (low Mach) levels, none once
    !! compressibility dominates. The form and mg_om_r0 are fitted to TGV iteration counts from Mach 0.1 to 0.001
    impure subroutine s_mg_omega()

        real(wp), dimension(2, mg_maxlev) :: sums, sums_glb
        real(wp)                          :: sd, sk, dg, nb, r1
        integer                           :: lv, ii, jj, kk, idx, off, nx, ny, nz, ex, ey, gx, gy, gz, sy, sz

        gx = mg_gx; gy = mg_gy; gz = mg_gz
        do lv = 1, mg_nlev
            nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv); ex = nx + 2*gx; ey = ny + 2*gy; sy = ex; sz = ex*ey
            sd = 0._wp; sk = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[sd, sk]]', reductionOp='[+]')
            do kk = 0, nz - 1
                do jj = 0, ny - 1
                    do ii = 0, nx - 1
                        idx = ${MG_IX('ii', 'jj', 'kk')}$
                        @:MG_ROW()
                        sd = sd + mg_d(idx)
                        sk = sk + dg - mg_d(idx)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            sums(:,lv) = [sd, sk]
        end do
        call s_mpi_allreduce_vectors_sum(sums(:,1:mg_nlev), sums_glb(:,1:mg_nlev), 2, mg_nlev)
        do lv = 2, mg_nlev
            mg_om(lv - 1) = 1._wp + (proj_mg_omega - 1._wp)/(1._wp + sums_glb(1, lv)/(mg_om_r0*max(sums_glb(2, lv), tiny(1._wp))))
        end do

        ! Stop at the first level whose spacing spans proj_mg_trunc screening lengths. The pressure equation screens over
        ! l = h/sqrt(6 r1), r1 the fine level's Helmholtz-to-Laplacian diagonal ratio (l is about c*dt, the acoustic CFL in cells);
        ! error longer than l has eigenvalues near the Helmholtz diagonal, which CG resolves cheaply, so the coarser, rank-spanning
        ! levels and the gathered bottom add little. The coarse ratios cannot show this: summing aggregates doubles a level's
        ! ratio, where the spacing squared quadruples
        mg_blev = mg_nlev
        if (proj_mg_trunc > 0._wp) then
            r1 = sums_glb(1, 1)/max(sums_glb(2, 1), tiny(1._wp))
            do lv = 2, mg_nlev
                if (6._wp*r1*4._wp**(lv - 1) >= proj_mg_trunc**2) then
                    mg_blev = lv
                    exit
                end if
            end do
        end if

        ! K-cycles pay off by keeping iterations flat as the hierarchy deepens with the grid; a truncated hierarchy stops at a
        ! fixed depth set by the screening length, so with proj_mg_kcycle = -1 it keeps V-cycles, cheaper and fewer iterations
        mg_kspace = mg_kfull
        if (proj_mg_kcycle < 0 .and. mg_blev < mg_nlev) mg_kspace = 0

    end subroutine s_mg_omega

    !> Every rank's bottom-level block size, for the global numbering of the bottom cells, and the bottom hierarchy's levels
    impure subroutine s_mg_bottom_layout()

        integer, dimension(3) :: pc
        integer               :: r, lv, i, ierr

        allocate (crs_cnt(num_procs), crs_disp(num_procs), crs_sz(3, num_procs))
        crs_sz(:,proc_rank + 1) = [mg_nx(mg_nlev), mg_ny(mg_nlev), mg_nz(mg_nlev)]
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, crs_sz, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        crs_n = 0
        do r = 1, num_procs
            crs_cnt(r) = product(crs_sz(:,r)); crs_disp(r) = crs_n; crs_n = crs_n + crs_cnt(r)
        end do

        crs_nl = 1; crs_ln(1) = crs_n; crs_loff(1) = 0; crs_tot = crs_n
        crs_ld(:,1) = [num_procs_x, num_procs_y, num_procs_z]
        do while (crs_ln(crs_nl) > mg_bottom_max .and. .not. crs_exact)
            crs_ld(:,crs_nl + 1) = (crs_ld(:,crs_nl) + 1)/2
            crs_loff(crs_nl + 1) = crs_tot; crs_ln(crs_nl + 1) = product(crs_ld(:,crs_nl + 1))
            crs_nl = crs_nl + 1; crs_tot = crs_tot + crs_ln(crs_nl)
        end do
        allocate (crs_ad(crs_tot), crs_ak(6, crs_tot), crs_aj(6, crs_tot), crs_agg(crs_tot), crs_co(3, crs_tot))
        crs_cg = crs_ln(crs_nl) > mg_bottom_max
        allocate (crs_x(crs_tot), crs_b(crs_tot))
        if (.not. crs_cg) allocate (crs_l(crs_ln(crs_nl), crs_ln(crs_nl)))
        if (crs_nl == 1) return

        if (any(crs_cnt /= 1) .or. product(crs_ld(:,1)) /= crs_n) then
            call s_mpi_abort('projection: a bottom level past mg_bottom_max cells must hold one cell per rank')
        end if
        pc = 0; pc(1:num_dims) = proc_coords(1:num_dims)
        crs_co(:,proc_rank + 1) = pc
#ifdef MFC_MPI
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, crs_co, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
        do lv = 1, crs_nl - 1
            do i = crs_loff(lv) + 1, crs_loff(lv) + crs_ln(lv)
                pc = crs_co(:,i)/2
                crs_agg(i) = crs_loff(lv + 1) + 1 + pc(1) + crs_ld(1, lv + 1)*(pc(2) + crs_ld(2, lv + 1)*pc(3))
                crs_co(:,crs_agg(i)) = pc
            end do
        end do

    end subroutine s_mg_bottom_layout

    !> Global bottom index of the neighbor of this rank's bottom cell (i, j, k) across side q, or 0 where a boundary has no seam
    integer function f_crs_nbr(q, i, j, k) result(g)

        integer, intent(in)   :: q, i, j, k
        integer, dimension(3) :: ic
        integer               :: nd, r

        ic = [i, j, k]; nd = (q + 1)/2; r = proc_rank; g = 0
        if (mod(q, 2) == 1) then
            ic(nd) = ic(nd) - 1
            if (ic(nd) < 0) then
                r = mg_nbr(q); if (r < 0) return
                ic(nd) = crs_sz(nd, r + 1) - 1
            end if
        else
            ic(nd) = ic(nd) + 1
            if (ic(nd) >= crs_sz(nd, proc_rank + 1)) then
                r = mg_nbr(q); if (r < 0) return
                ic(nd) = 0
            end if
        end if
        g = crs_disp(r + 1) + (ic(3)*crs_sz(2, r + 1) + ic(2))*crs_sz(1, r + 1) + ic(1) + 1

    end function f_crs_nbr

    !> Gather every bottom cell's row (its diagonal sum and its face conductances with their global neighbors) into the global
    !! bottom stencil, coarsen it as the device levels are, and factor the last level; a face whose two sides fall in one cell (a
    !! periodic seam across a one-cell width, or a face inside an aggregated box) couples the cell to itself and cancels
    impure subroutine s_mg_bottom_build()

        real(wp), dimension(13, crs_cnt(proc_rank + 1)) :: row
        real(wp), dimension(13, crs_n)                  :: rows
        integer                                         :: i, j, k, c, q, g, jr, idx, off, ex, ey, gx, gy, gz, nt, ierr
        integer                                         :: lv, o, nn, ci, cj

        off = mg_off(mg_nlev); gx = mg_gx; gy = mg_gy; gz = mg_gz
        ex = mg_nx(mg_nlev) + 2*gx; ey = mg_ny(mg_nlev) + 2*gy; nt = ex*ey*(mg_nz(mg_nlev) + 2*gz)
        $:GPU_UPDATE(host='[mg_d(off + 1:off + nt), mg_kx(off + 1:off + nt), mg_ky(off + 1:off + nt), mg_kz(off + 1:off + nt)]')
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; idx = ${MG_IX('i', 'j', 'k')}$
                    row(:,c) = 0._wp
                    row(1, c) = mg_d(idx)
                    row(3, c) = mg_kx(idx); row(5, c) = mg_kx(idx + 1)
                    if (num_dims > 1) then
                        row(7, c) = mg_ky(idx); row(9, c) = mg_ky(idx + ex)
                    end if
                    if (num_dims > 2) then
                        row(11, c) = mg_kz(idx); row(13, c) = mg_kz(idx + ex*ey)
                    end if
                    do q = 1, 2*num_dims
                        row(2*q, c) = real(f_crs_nbr(q, i, j, k), wp)
                    end do
                end do
            end do
        end do

        rows(:,1:c) = row
#ifdef MFC_MPI
        call MPI_ALLGATHERV(row, 13*c, mpi_p, rows, 13*crs_cnt, 13*crs_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif

        do i = 1, crs_n
            crs_ad(i) = rows(1, i)
            do q = 1, 6
                g = nint(rows(2*q, i))
                ! No neighbour: a domain face, whose conductance (a pressure outlet's; zero on a wall) couples to a fixed ghost
                if (g == 0) crs_ad(i) = crs_ad(i) + rows(2*q + 1, i)
                if (g == i) g = 0
                crs_aj(q, i) = max(g, 0)
                crs_ak(q, i) = merge(rows(2*q + 1, i), 0._wp, g > 0)
                crs_ad(i) = crs_ad(i) + crs_ak(q, i)
            end do
        end do

        do lv = 1, crs_nl - 1
            o = crs_loff(lv + 1); nn = crs_ln(lv + 1)
            crs_ad(o + 1:o + nn) = 0._wp; crs_ak(:,o + 1:o + nn) = 0._wp; crs_aj(:,o + 1:o + nn) = 0
            do i = crs_loff(lv) + 1, crs_loff(lv) + crs_ln(lv)
                ci = crs_agg(i); crs_ad(ci) = crs_ad(ci) + crs_ad(i)
                do q = 1, 6
                    if (crs_aj(q, i) == 0) cycle
                    cj = crs_agg(crs_aj(q, i))
                    if (cj == ci) then
                        crs_ad(ci) = crs_ad(ci) - crs_ak(q, i)
                    else
                        crs_aj(q, ci) = cj; crs_ak(q, ci) = crs_ak(q, ci) + crs_ak(q, i)
                    end if
                end do
            end do
        end do

        if (crs_cg) return
        o = crs_loff(crs_nl); nn = crs_ln(crs_nl)
        crs_l = 0._wp
        do i = 1, nn
            crs_l(i, i) = crs_ad(o + i)
            do q = 1, 6
                j = crs_aj(q, o + i)
                if (j > 0) crs_l(i, j - o) = crs_l(i, j - o) - crs_ak(q, o + i)
            end do
        end do

        ! In-place Cholesky, lower triangle
        do jr = 1, nn
            crs_l(jr, jr) = sqrt(crs_l(jr, jr) - sum(crs_l(jr,1:jr - 1)**2))
            do i = jr + 1, nn
                crs_l(i, jr) = (crs_l(i, jr) - sum(crs_l(i,1:jr - 1)*crs_l(jr,1:jr - 1)))/crs_l(jr, jr)
            end do
        end do

    end subroutine s_mg_bottom_build

    !> Bottom solve: gather the ranks' right-hand sides, run the bottom hierarchy's V-cycle on them, keep this rank's block
    impure subroutine s_mg_bottom_solve()

        real(wp), dimension(crs_cnt(proc_rank + 1)) :: v
        integer                                     :: i, j, k, c, idx, off, ex, ey, gx, gy, gz, nt, ierr

        off = mg_off(mg_nlev); gx = mg_gx; gy = mg_gy; gz = mg_gz
        ex = mg_nx(mg_nlev) + 2*gx; ey = mg_ny(mg_nlev) + 2*gy; nt = ex*ey*(mg_nz(mg_nlev) + 2*gz)
        $:GPU_UPDATE(host='[mg_f(off + 1:off + nt)]')
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; v(c) = mg_f(${MG_IX('i', 'j', 'k')}$)
                end do
            end do
        end do
        crs_b(1:c) = v
#ifdef MFC_MPI
        call MPI_ALLGATHERV(v, c, mpi_p, crs_b, crs_cnt, crs_disp, mpi_p, MPI_COMM_WORLD, ierr)
#endif
        if (crs_cg) then
            call s_crs_cg()
        else
            call s_crs_cycle(1)
        end if

        ! The ghosts stay zero: nothing reads the bottom level's
        mg_e(off + 1:off + nt) = 0._wp
        c = 0
        do k = 0, mg_nz(mg_nlev) - 1
            do j = 0, mg_ny(mg_nlev) - 1
                do i = 0, mg_nx(mg_nlev) - 1
                    c = c + 1; mg_e(${MG_IX('i', 'j', 'k')}$) = crs_x(crs_disp(proc_rank + 1) + c)
                end do
            end do
        end do
        $:GPU_UPDATE(device='[mg_e(off + 1:off + nt)]')

    end subroutine s_mg_bottom_solve

    !> One V-cycle on bottom level lv, crs_b into crs_x: proj_mg_sweeps forward Gauss-Seidel sweeps, the coarser levels, the same
    !! sweeps backward; the last level is solved densely. It stays linear (no K-cycle), as a bottom solve that varies with its
    !! right-hand side slows the outer flexible CG markedly. Every rank runs it on the same data, so all agree
    recursive subroutine s_crs_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i, s, o, nn

        o = crs_loff(lv); nn = crs_ln(lv)
        if (lv == crs_nl) then
            do i = 1, nn
                crs_x(o + i) = (crs_b(o + i) - sum(crs_l(i,1:i - 1)*crs_x(o + 1:o + i - 1)))/crs_l(i, i)
            end do
            do i = nn, 1, -1
                crs_x(o + i) = (crs_x(o + i) - sum(crs_l(i + 1:nn,i)*crs_x(o + i + 1:o + nn)))/crs_l(i, i)
            end do
            return
        end if
        crs_x(o + 1:o + nn) = 0._wp
        crs_b(crs_loff(lv + 1) + 1:crs_loff(lv + 1) + crs_ln(lv + 1)) = 0._wp
        do s = 1, proj_mg_sweeps
            call s_crs_gs(lv, 1)
        end do
        do i = o + 1, o + nn
            crs_b(crs_agg(i)) = crs_b(crs_agg(i)) + crs_b(i) - f_crs_ax(i)
        end do
        call s_crs_cycle(lv + 1)
        do i = o + 1, o + nn
            crs_x(i) = crs_x(i) + crs_x(crs_agg(i))
        end do
        do s = 1, proj_mg_sweeps
            call s_crs_gs(lv, -1)
        end do

    end subroutine s_crs_cycle

    !> Exact bottom solve past mg_bottom_max cells, too many to factor densely: Jacobi-preconditioned CG on the level-1 stencil,
    !! crs_b into crs_x, from zero to crs_tol relative to crs_b. Every rank runs it on the same data, so all agree
    subroutine s_crs_cg()

        real(wp), dimension(crs_n) :: r, z, d, q
        real(wp)                   :: rz, rz_new, alpha, bnorm
        integer                    :: it, i, f

        crs_x(1:crs_n) = 0._wp; r = crs_b(1:crs_n)
        bnorm = sqrt(sum(r*r))
        if (bnorm == 0._wp) return
        z = r/crs_ad(1:crs_n); d = z; rz = sum(r*z)
        do it = 1, crs_n
            do i = 1, crs_n
                q(i) = crs_ad(i)*d(i)
                do f = 1, 6
                    if (crs_aj(f, i) > 0) q(i) = q(i) - crs_ak(f, i)*d(crs_aj(f, i))
                end do
            end do
            alpha = rz/sum(d*q)
            crs_x(1:crs_n) = crs_x(1:crs_n) + alpha*d; r = r - alpha*q
            if (sqrt(sum(r*r)) <= crs_tol*bnorm) exit
            z = r/crs_ad(1:crs_n); rz_new = sum(r*z)
            d = z + (rz_new/rz)*d; rz = rz_new
        end do

    end subroutine s_crs_cg

    !> (A x)_i of bottom-hierarchy cell i
    pure real(wp) function f_crs_ax(i) result(ax)
        integer, intent(in) :: i
        integer             :: q

        ax = crs_ad(i)*crs_x(i)
        do q = 1, 6
            if (crs_aj(q, i) > 0) ax = ax - crs_ak(q, i)*crs_x(crs_aj(q, i))
        end do

    end function f_crs_ax

    !> One Gauss-Seidel sweep of crs_x on bottom level lv, in cell order (dir = 1) or reversed (dir = -1)
    subroutine s_crs_gs(lv, dir)

        integer, intent(in) :: lv, dir
        integer             :: i, q, i0, i1
        real(wp)            :: s

        i0 = crs_loff(lv) + 1; i1 = crs_loff(lv) + crs_ln(lv)
        if (dir < 0) then
            i0 = i1; i1 = crs_loff(lv) + 1
        end if
        do i = i0, i1, dir
            s = crs_b(i)
            do q = 1, 6
                if (crs_aj(q, i) > 0) s = s + crs_ak(q, i)*crs_x(crs_aj(q, i))
            end do
            crs_x(i) = s/crs_ad(i)
        end do

    end subroutine s_crs_gs

    !> Fill the ghost layer of the block of mg_e at offset off, laid out as level lv, from the neighbors: one message per distinct
    !! neighbor rank holding all the sides it shares with this one, and periodic seams a rank owns alone copied in place
    impure subroutine s_mg_exchange(lv, off)

        integer, intent(in)    :: lv, off
        integer, dimension(12) :: req
        integer                :: i, nt, ierr

        call s_mg_sides(lv, off, 1)
        if (mg_nmsg == 0) return
#ifdef MFC_MPI
        if (rdma_mpi) then
            $:GPU_WAIT()
            #:call GPU_HOST_DATA(use_device_addr='[mg_sbuf, mg_rbuf]')
                do i = 1, mg_nmsg
                    call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), &
                                   & ierr)
                    call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                                   & req(mg_nmsg + i), ierr)
                end do
                call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            #:endcall GPU_HOST_DATA
        else
            nt = mg_sbeg(mg_nmsg, lv) + mg_slen(mg_nmsg, lv)
            $:GPU_UPDATE(host='[mg_sbuf(1:nt)]')
            do i = 1, mg_nmsg
                call MPI_IRECV(mg_rbuf(mg_rbeg(i, lv) + 1), mg_rlen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, req(i), ierr)
                call MPI_ISEND(mg_sbuf(mg_sbeg(i, lv) + 1), mg_slen(i, lv), mpi_p, mg_nlist(i), 7101, MPI_COMM_WORLD, &
                               & req(mg_nmsg + i), ierr)
            end do
            call MPI_WAITALL(2*mg_nmsg, req, MPI_STATUSES_IGNORE, ierr)
            nt = mg_rbeg(mg_nmsg, lv) + mg_rlen(mg_nmsg, lv)
            $:GPU_UPDATE(device='[mg_rbuf(1:nt)]')
        end if
#endif
        call s_mg_sides(lv, off, 2)

    end subroutine s_mg_exchange

    !> Message layout of every level, fixed by the neighbor ranks: this rank sends its sides in its own side order, and receives
    !! each neighbor's sides in that neighbor's side order, so both ends agree without further communication
    subroutine s_mg_exchange_layout()

        integer :: lv, q, qq, i, tot

        mg_nmsg = 0
        do q = 1, 6
            mg_smode(q) = 0
            if (mg_nbr(q) < 0) cycle
            mg_smode(q) = merge(3, 1, mg_nbr(q) == proc_rank)
            if (mg_smode(q) == 3 .or. any(mg_nlist(1:mg_nmsg) == mg_nbr(q))) cycle
            mg_nmsg = mg_nmsg + 1; mg_nlist(mg_nmsg) = mg_nbr(q)
        end do
        mg_boff = 0
        do lv = 1, mg_nlev
            tot = 0
            do i = 1, mg_nmsg
                mg_sbeg(i, lv) = tot
                do q = 1, 6
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 1) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_slen(i, lv) = tot - mg_sbeg(i, lv)
            end do
            tot = 0
            do i = 1, mg_nmsg
                mg_rbeg(i, lv) = tot
                do qq = 1, 6
                    q = f_opposite(qq)
                    if (mg_nbr(q) /= mg_nlist(i)) cycle
                    mg_boff(q, lv, 2) = tot; tot = tot + f_side_size(lv, q)
                end do
                mg_rlen(i, lv) = tot - mg_rbeg(i, lv)
            end do
        end do
        $:GPU_UPDATE(device='[mg_smode, mg_boff]')

    end subroutine s_mg_exchange_layout

    !> Side q of the opposite end of the same direction
    pure integer function f_opposite(q)

        integer, intent(in) :: q

        f_opposite = q + merge(1, -1, mod(q, 2) == 1)

    end function f_opposite

    !> Cells in the layer of side q on level lv
    pure integer function f_side_size(lv, q)

        integer, intent(in) :: lv, q

        select case ((q + 1)/2)
        case (1); f_side_size = mg_ny(lv)*mg_nz(lv)
        case (2); f_side_size = mg_nx(lv)*mg_nz(lv)
        case default; f_side_size = mg_nx(lv)*mg_ny(lv)
        end select

    end function f_side_size

    !> All six sides of the block at offset off laid out as level lv, in one kernel. Phase 1 copies each periodic seam a rank owns
    !! alone into its ghost layer and packs each side facing another rank into mg_sbuf; phase 2 unpacks mg_rbuf into those ghosts
    impure subroutine s_mg_sides(lv, off, phase)

        integer, intent(in) :: lv, off, phase
        integer             :: nx, ny, nz, ex, ey, base, nmax, q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx, lvl, ph

        ! Dummies may alias host variables, which a device kernel must not reference

        lvl = lv; ph = phase
        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); ex = nx + 2*mg_gx; ey = ny + 2*mg_gy
        base = off + (mg_gz*ey + mg_gy)*ex + mg_gx + 1
        nmax = max(nx, ny, nz)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[q, a, b, md, n1, n2, s1, s2, st, nn, lay, gho, idx, bidx]')
        do q = 1, 6
            do b = 0, nmax - 1
                do a = 0, nmax - 1
                    ! Phase 1 acts on self seams (3) and packs (1), phase 2 unpacks (becomes 2)
                    md = mg_smode(q)
                    if (ph == 2) md = merge(2, 0, md == 1)
                    if ((q + 1)/2 == 1) then
                        nn = nx; st = 1; n1 = ny; s1 = ex; n2 = nz; s2 = ex*ey
                    else if ((q + 1)/2 == 2) then
                        nn = ny; st = ex; n1 = nx; s1 = 1; n2 = nz; s2 = ex*ey
                    else
                        nn = nz; st = ex*ey; n1 = nx; s1 = 1; n2 = ny; s2 = ex
                    end if
                    if (md > 0 .and. a < n1 .and. b < n2) then
                        lay = merge(0, nn - 1, mod(q, 2) == 1)
                        gho = merge(-1, nn, mod(q, 2) == 1)
                        idx = base + a*s1 + b*s2
                        bidx = mg_boff(q, lvl, merge(1, 2, md == 1)) + b*n1 + a + 1
                        if (md == 1) then
                            mg_sbuf(bidx) = mg_e(idx + lay*st)
                        else if (md == 2) then
                            mg_e(idx + gho*st) = mg_rbuf(bidx)
                        else
                            mg_e(idx + gho*st) = mg_e(idx + (nn - 1 - lay)*st)
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_sides

    !> Multigrid preconditioner on rs, returned in zs, coupled across ranks
    impure subroutine s_mg_precond()

        integer :: j, k, l, off, ex, ey, gx, gy, gz

        off = mg_off(1); gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = mg_nx(1) + 2*gx; ey = mg_ny(1) + 2*gy
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -gx, m + gx
                    mg_e(${MG_IX('j', 'k', 'l')}$) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_f(${MG_IX('j', 'k', 'l')}$) = rs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mg_cycle(1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    zs(j, k, l) = mg_e(${MG_IX('j', 'k', 'l')}$)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_precond

    !> One cycle on level lv, from mg_f(lv) into mg_e(lv), which must start at zero, ghosts included. Each smoothing phase freezes
    !! the ghost layer for all its sweeps; pre-sweeps run red then black and post-sweeps black then red, so a V-cycle is a symmetric
    !! operator. The coarse problem gets one cycle, or every mg_kspace levels two flexible-CG steps (s_mg_kstep)
    recursive impure subroutine s_mg_cycle(lv)

        integer, intent(in) :: lv
        integer             :: i

        if (lv == mg_nlev) then
            call nvtxStartRange("PROJ-MG-BOTTOM")
            call s_mg_bottom_solve()
            call nvtxEndRange
            return
        end if
        ! Two exchanges per level: pre-smoothing starts from zero, so its frozen ghosts are already current; the residual and the
        ! post-smoothing each take one exchange
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 0); call s_mg_smooth(lv, 1)
        end do
        call nvtxStartRange("PROJ-MG-EXCHANGE")
        call s_mg_exchange(lv, mg_off(lv))
        call nvtxEndRange
        ! At a truncation level (mg_blev < mg_nlev) the cycle is its two smoothing phases alone, still symmetric; CG takes the
        ! screened error longer than the level's spacing
        if (lv < mg_blev) then
            call s_mg_restrict(lv)
            if (mg_kspace > 0 .and. lv + 1 < mg_blev .and. mod(lv, max(mg_kspace, 1)) == 0) then
                call s_mg_kstep(lv + 1)
            else
                call s_mg_cycle(lv + 1)
            end if
            call s_mg_prolong(lv)
            call nvtxStartRange("PROJ-MG-EXCHANGE")
            call s_mg_exchange(lv, mg_off(lv))
            call nvtxEndRange
        end if
        do i = 1, proj_mg_sweeps
            call s_mg_smooth(lv, 1); call s_mg_smooth(lv, 0)
        end do

    end subroutine s_mg_cycle

    !> K-cycle on coarse level c (Notay & Vassilevski, NLAA 15, 2008): two flexible-CG steps on its system, mg_f(c) into mg_e(c),
    !! each preconditioned by one cycle; the second is skipped once the first has cut the residual by mg_kt. The dot products of a
    !! step are fused into one reduction
    recursive impure subroutine s_mg_kstep(c)

        integer, intent(in)       :: c
        real(wp), dimension(5, 1) :: sl, sg
        real(wp)                  :: s1, s2, s3, s4, s5, a1, rho2, cc, cd, dg, nb
        integer                   :: ii, jj, kk, idx, nx, ny, nz, off, ex, ey, gx, gy, gz, sy, sz

        nx = mg_nx(c); ny = mg_ny(c); nz = mg_nz(c); off = mg_off(c)
        gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = nx + 2*gx; ey = ny + 2*gy; sy = ex; sz = ex*ey

        ! Step 1: c1 = cycle(r), v = A c1
        call s_mg_cycle(c)
        call s_mg_exchange(c, off)
        s1 = 0._wp; s2 = 0._wp; s3 = 0._wp; s4 = 0._wp; s5 = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[s1, s2, s3, s4, s5]]', reductionOp='[+]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = ${MG_IX('ii', 'jj', 'kk')}$
                    @:MG_ROW()
                    mg_kr(idx) = mg_f(idx); mg_kc(idx) = mg_e(idx); mg_kv(idx) = dg*mg_e(idx) - nb
                    s1 = s1 + mg_e(idx)*mg_kv(idx); s2 = s2 + mg_e(idx)*mg_f(idx); s3 = s3 + mg_f(idx)*mg_f(idx)
                    s4 = s4 + mg_kv(idx)*mg_f(idx); s5 = s5 + mg_kv(idx)*mg_kv(idx)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        sl(:,1) = [s1, s2, s3, s4, s5]
        call s_mpi_allreduce_vectors_sum(sl, sg, 5, 1)
        s1 = sg(1, 1); s2 = sg(2, 1); s3 = sg(3, 1); s4 = sg(4, 1); s5 = sg(5, 1)
        a1 = s2/s1
        cc = a1; cd = 0._wp

        ! Step 2, unless |r - a1 v| <= mg_kt |r| already: d = cycle(r - a1 v), from zero
        if (s3 - 2._wp*a1*s4 + a1*a1*s5 > mg_kt*mg_kt*s3) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
            do kk = -gz, nz - 1 + gz
                do jj = -gy, ny - 1 + gy
                    do ii = -gx, nx - 1 + gx
                        idx = ${MG_IX('ii', 'jj', 'kk')}$
                        mg_e(idx) = 0._wp
                        if (ii >= 0 .and. ii < nx .and. jj >= 0 .and. jj < ny .and. kk >= 0 .and. kk < nz) then
                            mg_f(idx) = mg_kr(idx) - a1*mg_kv(idx)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mg_cycle(c)
            call s_mg_exchange(c, off)
            s2 = 0._wp; s4 = 0._wp; s5 = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]', reduction='[[s2, s4, s5]]', reductionOp='[+]')
            do kk = 0, nz - 1
                do jj = 0, ny - 1
                    do ii = 0, nx - 1
                        idx = ${MG_IX('ii', 'jj', 'kk')}$
                        @:MG_ROW()
                        s2 = s2 + mg_e(idx)*mg_kv(idx); s4 = s4 + mg_e(idx)*(dg*mg_e(idx) - nb); s5 = s5 + mg_e(idx)*mg_f(idx)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            sl(1:3,1) = [s2, s4, s5]
            call s_mpi_allreduce_vectors_sum(sl(1:3,:), sg(1:3,:), 3, 1)
            ! gamma = d.v, beta = d.A d, alpha2 = d.(r - a1 v); rho2 = beta - gamma^2/rho1
            rho2 = sg(2, 1) - sg(1, 1)**2/s1
            if (rho2 > 0._wp) then
                cd = sg(3, 1)/rho2; cc = a1 - sg(1, 1)*cd/s1
            end if
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = ${MG_IX('ii', 'jj', 'kk')}$
                    mg_e(idx) = cc*mg_kc(idx) + cd*mg_e(idx)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_kstep

    !> One red-black Gauss-Seidel half-sweep of the given color on level lv, one thread per cell of that color
    impure subroutine s_mg_smooth(lv, color)

        integer, intent(in) :: lv, color
        integer             :: nx, ny, nz, off, ex, ey, gx, gy, gz, sy, sz, ih, ii, jj, kk, idx, cl
        real(wp)            :: dg, nb

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = nx + 2*gx; ey = ny + 2*gy; sy = ex; sz = ex*ey; cl = color
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ih, ii, jj, kk, idx, dg, nb]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ih = 0, (nx - 1)/2
                    ii = 2*ih + mod(jj + kk + cl, 2)
                    if (ii < nx) then
                        idx = ${MG_IX('ii', 'jj', 'kk')}$
                        @:MG_ROW()
                        mg_e(idx) = (mg_f(idx) + nb)/dg
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_smooth

    !> Residual of level lv (ghosts current) summed onto the coarse right-hand side; restriction is the transpose of prolongation
    impure subroutine s_mg_restrict(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, off, ex, ey, gx, gy, gz, sy, sz, cnx, cny, cnz, coff, cidx, ii, jj, kk, idx, a, b, c
        integer             :: a0, a1, b0, b1, c0, c1
        real(wp)            :: dg, nb

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = nx + 2*gx; ey = ny + 2*gy; sy = ex; sz = ex*ey
        cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1); coff = mg_off(lv + 1)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, idx, dg, nb]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    idx = ${MG_IX('ii', 'jj', 'kk')}$
                    @:MG_ROW()
                    mg_r(idx) = mg_f(idx) - (dg*mg_e(idx) - nb)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! The coarse correction starts from zero, ghosts included, so its first sweep needs no exchange
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk]')
        do kk = -gz, cnz - 1 + gz
            do jj = -gy, cny - 1 + gy
                do ii = -gx, cnx - 1 + gx
                    mg_e(coff + ((kk + gz)*(cny + 2*gy) + jj + gy)*(cnx + 2*gx) + ii + gx + 1) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Gather per coarse cell rather than scatter with atomics, so the sum order is fixed and results are reproducible
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk, a, b, c, cidx, dg, a0, a1, b0, b1, c0, c1]')
        do kk = 0, cnz - 1
            do jj = 0, cny - 1
                do ii = 0, cnx - 1
                    dg = 0._wp
                    @:MG_CHILDREN(a0, a1, ii, cnx, nx)
                    @:MG_CHILDREN(b0, b1, jj, cny, ny)
                    @:MG_CHILDREN(c0, c1, kk, cnz, nz)
                    $:GPU_LOOP(parallelism='[seq]')
                    do c = c0, c1
                        $:GPU_LOOP(parallelism='[seq]')
                        do b = b0, b1
                            $:GPU_LOOP(parallelism='[seq]')
                            do a = a0, a1
                                dg = dg + mg_r(${MG_IX('a', 'b', 'c')}$)
                            end do
                        end do
                    end do
                    cidx = coff + ((kk + gz)*(cny + 2*gy) + jj + gy)*(cnx + 2*gx) + ii + gx + 1
                    mg_f(cidx) = dg
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_restrict

    !> Add each coarse correction, scaled by mg_om(lv) (see s_mg_omega), to its children; a scale in (0, 2) keeps the V-cycle
    !! symmetric positive definite
    impure subroutine s_mg_prolong(lv)

        integer, intent(in) :: lv
        integer             :: nx, ny, nz, off, ex, ey, gx, gy, gz, cnx, cny, cnz, coff, ii, jj, kk
        real(wp)            :: om

        nx = mg_nx(lv); ny = mg_ny(lv); nz = mg_nz(lv); off = mg_off(lv)
        gx = mg_gx; gy = mg_gy; gz = mg_gz; ex = nx + 2*gx; ey = ny + 2*gy
        cnx = mg_nx(lv + 1); cny = mg_ny(lv + 1); cnz = mg_nz(lv + 1); coff = mg_off(lv + 1); om = mg_om(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ii, jj, kk]')
        do kk = 0, nz - 1
            do jj = 0, ny - 1
                do ii = 0, nx - 1
                    mg_e(${MG_IX('ii', 'jj', 'kk')}$) = mg_e(${MG_IX('ii', 'jj', 'kk')}$) + om*mg_e(coff + ((min(kk/2, &
                         & cnz - 1) + gz)*(cny + 2*gy) + min(jj/2, cny - 1) + gy)*(cnx + 2*gx) + min(ii/2, cnx - 1) + gx + 1)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_mg_prolong

    impure subroutine s_finalize_projection_module()

        $:GPU_EXIT_DATA(detach='[pk_sf(1)%sf]')
        @:DEALLOCATE(uf, uf0, divu, rhs_p, p_stage, p_step0, pflx, rhoc, dcoef, bvec, xs, rs, zs, qs, pk, kap, gnd, fbx, fby, fbz)
        @:DEALLOCATE(mg_d, mg_kx, mg_ky, mg_kz, mg_e, mg_f, mg_r, mg_kr, mg_kc, mg_kv)
        if (allocated(crs_l)) deallocate (crs_l)
        deallocate (crs_x, crs_b, crs_ad, crs_ak, crs_aj, crs_agg, crs_co, crs_cnt, crs_disp, crs_sz)
        @:UNPIN_HOST(mg_sbuf, mg_rbuf)
        @:DEALLOCATE(mg_sbuf, mg_rbuf)

    end subroutine s_finalize_projection_module

end module m_projection
