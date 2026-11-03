!>
!!@file
!!@brief Contains module m_amr_distribution

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Rank decomposition, box arithmetic and block ownership (SFC cut, owner lookups, owned-block lists).
module m_amr_distribution

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_box, only: f_morton
    use m_global_parameters
    use m_constants, only: K_ib
    use m_mpi_proxy, only: s_mpi_abort
    use m_ibm, only: ib_markers
    use m_amr_state

    implicit none

    private
    public :: f_amr_overlap_count, f_amr_rank_overlaps, s_amr_adopt_block_owners, s_amr_assign_block_owners, s_amr_body_bbox, &
        & s_amr_box_isect, s_amr_compute_isect, s_amr_expand_box_over_bodies, s_amr_rank_coarse_range, s_amr_rank_decomp, &
        & s_amr_rank_interior, s_amr_ranks_overlapping, s_amr_refresh_lists, s_amr_refresh_my_blocks, s_amr_sfc_cut, &
        & s_amr_sort_by_key, s_amr_tile_box, s_amr_validate_decomp, s_amr_validate_owner

contains

    !> Rank r's coarse-grid decomposition (start_idx + local extent m/n/p), computed O(1) from its cartesian coords. The domain cart
    !! is MPI_CART_CREATE(reorder=.false.), so ranks keep MPI_COMM_WORLD order and coords are row-major in r. The per-dim split is
    !! the integer split-with-remainder of s_mpi_decompose_computational_domain.
    pure subroutine s_amr_rank_decomp(r, sidx, ext)

        integer, intent(in)  :: r
        integer, intent(out) :: sidx(3), ext(3)
        integer              :: coords(3), gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        ! row-major: r = cx*(py*pz) + cy*pz + cz  (py=pz=1 when that dim is not decomposed)
        coords(1) = r/(pd(2)*pd(3))
        coords(2) = mod(r, pd(2)*pd(3))/pd(3)
        coords(3) = mod(r, pd(3))
        sidx = 0; ext = 0
        do d = 1, num_dims
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            ext(d) = base - 1 + merge(1, 0, coords(d) < rem)
            sidx(d) = coords(d)*base + min(coords(d), rem)
        end do

    end subroutine s_amr_rank_decomp

    !> Guard: the computed decomposition must reproduce this rank's actual one (each rank checks itself, so all r are covered).
    impure subroutine s_amr_validate_decomp()

        integer :: sidx(3), ext(3)

        call s_amr_rank_decomp(proc_rank, sidx, ext)
        if (any(sidx /= amr_sidx .or. ext /= amr_ext)) then
            call s_mpi_abort('s_amr_rank_decomp does not reproduce this rank''s decomposition - computed split disagrees with ' &
                             & // 's_mpi_decompose_computational_domain')
        end if

    end subroutine s_amr_validate_decomp

    !> Closed-form inverse of the per-dim split in s_amr_rank_decomp: the cart coord owning global coarse cell g, with
    !! base=(gd+1)/pd, rem=mod(gd+1,pd). Exact for g in [0,gd]; callers clamp for out-of-range g (ghost reach).
    pure integer function f_amr_cell_coord(g, base, rem) result(c)

        integer, intent(in) :: g, base, rem

        if (g < rem*(base + 1)) then
            c = g/(base + 1)
        else
            c = rem + (g - rem*(base + 1))/base
        end if

    end function f_amr_cell_coord

    !> Per-dim contiguous rank-coord range [clo:chi] whose owned coarse slab meets box [blo:bhi], clamped to [0,pd-1]. The boundary
    !! coords' coarse_range is extended by buff_size at the domain edge (s_amr_rank_coarse_range), and clamping reproduces both the
    !! interior (scatter) and coarse_range (gather) sets. Collapsed dims contribute coord 0.
    pure subroutine s_amr_coord_range(blo, bhi, clo, chi)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: clo(3), chi(3)
        integer              :: gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        clo = 0; chi = 0
        do d = 1, num_dims
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            clo(d) = min(max(f_amr_cell_coord(blo(d), base, rem), 0), pd(d) - 1)
            chi(d) = min(max(f_amr_cell_coord(bhi(d), base, rem), 0), pd(d) - 1)
        end do

    end subroutine s_amr_coord_range

    !> Ascending rank list overlapping coarse box [blo:bhi]: the coord brick is enumerated cx->cy->cz so r = cx*(Py*Pz)+cy*Pz+cz
    !! ascends, reproducing the 0..num_procs-1 scan order. Owner not excluded. Caller sizes ranks(:) >= f_amr_overlap_count.
    pure subroutine s_amr_ranks_overlapping(blo, bhi, ranks, nr)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: ranks(:)
        integer, intent(out) :: nr
        integer              :: clo(3), chi(3), cx, cy, cz, pyz

        call s_amr_coord_range(blo, bhi, clo, chi)
        pyz = num_procs_y*num_procs_z
        nr = 0
        do cx = clo(1), chi(1)
            do cy = clo(2), chi(2)
                do cz = clo(3), chi(3)
                    nr = nr + 1
                    ranks(nr) = cx*pyz + cy*num_procs_z + cz
                end do
            end do
        end do

    end subroutine s_amr_ranks_overlapping

    !> Does rank r's subdomain overlap coarse box [blo:bhi]? Answers membership without enumerating the overlap set, which is
    !! O(num_procs) for a box spanning the machine.
    pure logical function f_amr_rank_overlaps(blo, bhi, r) result(hit)

        integer, intent(in) :: blo(3), bhi(3), r
        integer             :: clo(3), chi(3), c(3)

        call s_amr_coord_range(blo, bhi, clo, chi)
        c(1) = r/(num_procs_y*num_procs_z)  ! same coord -> rank map as the enumeration
        c(2) = mod(r/num_procs_z, num_procs_y)
        c(3) = mod(r, num_procs_z)
        hit = all(c >= clo) .and. all(c <= chi)

    end function f_amr_rank_overlaps

    !> Overlap count only (allocation sizing), = product of the per-dim coord-range widths.
    pure integer function f_amr_overlap_count(blo, bhi) result(nr)

        integer, intent(in) :: blo(3), bhi(3)
        integer             :: clo(3), chi(3)

        call s_amr_coord_range(blo, bhi, clo, chi)
        nr = (chi(1) - clo(1) + 1)*(chi(2) - clo(2) + 1)*(chi(3) - clo(3) + 1)

    end function f_amr_overlap_count

    !> Rank r's contiguous owned coarse-cell range per dim: the interior plus buff_size ghosts where the subdomain touches the
    !! domain edge. One span, so box intersections identify contributors without a per-cell scan.
    pure subroutine s_amr_rank_coarse_range(r, crlo, crhi)

        integer, intent(in)  :: r
        integer, intent(out) :: crlo(3), crhi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        crlo = merge(merge(-buff_size, sidx, sidx == 0), 0, amr_dim)
        crhi = merge(sidx + ext + merge(buff_size, 0, sidx + ext == [m_glb, n_glb, p_glb]), 0, amr_dim)

    end subroutine s_amr_rank_coarse_range

    !> Per-dim intersection of two global boxes [alo:ahi] and [blo:bhi] -> [olo:ohi] (empty when olo > ohi in some dim).
    pure subroutine s_amr_box_isect(alo, ahi, blo, bhi, olo, ohi)

        integer, intent(in)  :: alo(3), ahi(3), blo(3), bhi(3)
        integer, intent(out) :: olo(3), ohi(3)

        olo = max(alo, blo); ohi = min(ahi, bhi)

    end subroutine s_amr_box_isect

    !> Per-block cost over each block's level-0 footprint (1 per cell + K_ib per IB-marked cell), summed over each rank's owned
    !! coarse cells and replicated by one MPI_ALLREDUCE(SUM).
    impure subroutine s_amr_block_cost(cost)

        real(wp), intent(out) :: cost(:)
        integer               :: k, j, kk, l, lo(3), hi(3)
        real(wp)              :: c

#ifdef MFC_MPI
        integer :: ierr
#endif

        if (ib) then
            $:GPU_UPDATE(host='[ib_markers%sf]')
        end if
        do k = 1, amr_num_blocks
            ! block footprint /\ this rank's coarse subdomain, in local interior indices (empty -> no-trip loops)
            lo = merge(max(amr_region_lo_all(:,k) - amr_sidx, 0), 0, amr_dim)
            hi = merge(min(amr_region_hi_all(:,k) - amr_sidx, amr_ext), 0, amr_dim)
            c = 0._wp
            do l = lo(3), hi(3)
                do kk = lo(2), hi(2)
                    do j = lo(1), hi(1)
                        c = c + 1._wp
                        if (ib) then
                            if (ib_markers%sf(j, kk, l) /= 0) c = c + K_ib
                        end if
                    end do
                end do
            end do
            cost(k) = c
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, cost, amr_num_blocks, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_block_cost

    !> Cost-weighted SFC partition of n items into num_procs contiguous Morton-key ranges: walk the key-sorted items (stable merge
    !! sort) and advance the owner rank when the running weight crosses the next even share. Writes owner(1:n) and
    !! cut(0:num_procs-1) (running Morton upper bound per rank; ranks with no item inherit the predecessor's bound, so cut is
    !! non-decreasing).
    subroutine s_amr_sfc_cut(keys, wt, n, cut, owner)

        integer, intent(in)          :: n
        integer(kind=8), intent(in)  :: keys(n)
        real(wp), intent(in)         :: wt(n)
        integer(kind=8), intent(out) :: cut(0:num_procs - 1)
        integer, intent(out)         :: owner(n)
        integer                      :: ord(n), k, r
        real(wp)                     :: total, cum, tgt, tol

        cut = -1_8
        if (n < 1) return

        call s_amr_sort_by_key(keys, n, ord)

        total = 0._wp
        do k = 1, n
            total = total + wt(k)
        end do

        ! chains-on-chains in SFC order; replicated weights in a fixed order, so every rank computes the identical assignment
        r = 0; cum = 0._wp
        do k = 1, n
            tgt = real(r + 1, wp)*total/real(num_procs, wp)
            ! tolerance = accumulated rounding bound: at an exact share boundary cum and tgt differ by rounding, and a fractional
            ! cost term could otherwise split identical weights unevenly (e.g. 5/3 instead of 4/4)
            tol = spacing(tgt)*real(n, wp)
            if (cum >= tgt - tol .and. r < num_procs - 1) r = r + 1
            owner(ord(k)) = r
            cut(r) = keys(ord(k))  ! items visited in ascending Morton key => running upper bound for rank r
            cum = cum + wt(ord(k))
        end do
        ! ranks with no item inherit the predecessor's bound so the search never lands on them
        do r = 1, num_procs - 1
            if (cut(r) < cut(r - 1)) cut(r) = cut(r - 1)
        end do

    end subroutine s_amr_sfc_cut

    !> ord(1:n) = item indices in ascending key order. Stable iterative bottom-up merge sort, a pure function of the input so every
    !! rank derives identical orders from replicated keys (differing assignments would deadlock or abort).
    pure subroutine s_amr_sort_by_key(keys, n, ord)

        integer, intent(in)         :: n
        integer(kind=8), intent(in) :: keys(n)
        integer, intent(out)        :: ord(n)
        integer                     :: mrg(n), k, width, lo_m, mid_m, hi_m, i_m, j_m, t_m

        do k = 1, n
            ord(k) = k
        end do
        width = 1
        do while (width < n)
            lo_m = 1
            do while (lo_m <= n - width)
                mid_m = lo_m + width - 1
                hi_m = min(lo_m + 2*width - 1, n)
                i_m = lo_m; j_m = mid_m + 1; t_m = lo_m
                do while (i_m <= mid_m .and. j_m <= hi_m)
                    ! <= keeps the left run first on ties: stability
                    if (keys(ord(i_m)) <= keys(ord(j_m))) then
                        mrg(t_m) = ord(i_m); i_m = i_m + 1
                    else
                        mrg(t_m) = ord(j_m); j_m = j_m + 1
                    end if
                    t_m = t_m + 1
                end do
                do while (i_m <= mid_m); mrg(t_m) = ord(i_m); i_m = i_m + 1; t_m = t_m + 1; end do
                do while (j_m <= hi_m); mrg(t_m) = ord(j_m); j_m = j_m + 1; t_m = t_m + 1; end do
                ord(lo_m:hi_m) = mrg(lo_m:hi_m)
                lo_m = lo_m + 2*width
            end do
            width = 2*width
        end do

    end subroutine s_amr_sort_by_key

    !> Refresh the owned-block list if any owner write has happened since the last refresh. Lazy, since s_amr_assign_block_owners is
    !! not the only writer of amr_block_owner.
    impure subroutine s_amr_refresh_my_blocks()

        integer :: b

        if (.not. amr_myblk_dirty) return
        if (allocated(amr_my_blk)) then
            if (size(amr_my_blk) < amr_max_blocks) deallocate (amr_my_blk)
        end if
        if (.not. allocated(amr_my_blk)) allocate (amr_my_blk(amr_max_blocks))
        amr_n_my = 0
        do b = 1, amr_num_blocks
            if (amr_block_owner(b) /= proc_rank) cycle
            amr_n_my = amr_n_my + 1
            amr_my_blk(amr_n_my) = b
        end do
        amr_myblk_dirty = .false.

    end subroutine s_amr_refresh_my_blocks

    !> Rebuild the epoch-keyed block lists when the mesh epoch has moved (one O(global blocks) walk per regrid, not per RK stage).
    !! Also caches the parent index (f_amr_parent_block is an O(global blocks) scan) and the children adjacency used by
    !! s_amr_sibling_face_weights.
    impure subroutine s_amr_refresh_lists()

        integer :: b, rlo(3), rhi(3), milo(3), mihi(3), bl(3), bh(3)
        integer :: crlo(3), crhi(3), plo(3), phi(3), pl(3), ph(3), mar, pblk, acc, nxt

        if (amr_l1r_epoch == amr_mesh_epoch .and. amr_l1r_nblk == amr_num_blocks) return
        #:for A in ['amr_l1r_blk', 'amr_l1p_blk', 'amr_fch_blk', 'amr_own_blk', 'amr_parent_blk', 'amr_child_idx']
            if (allocated(${A}$)) then
                if (size(${A}$) < amr_max_blocks) deallocate (${A}$)
            end if
            if (.not. allocated(${A}$)) allocate (${A}$ (amr_max_blocks))
        #:endfor
        if (allocated(amr_child_ptr)) then
            if (size(amr_child_ptr) < amr_max_blocks + 1) deallocate (amr_child_ptr)
        end if
        if (.not. allocated(amr_child_ptr)) allocate (amr_child_ptr(0:amr_max_blocks))
        call s_amr_rank_interior(proc_rank, milo, mihi)
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        mar = amr_cpat_mar
        amr_n_l1p = 0
        amr_n_l1r = 0
        amr_n_fch = 0
        amr_n_own = 0
        amr_child_ptr(0:amr_num_blocks) = 0
        do b = 1, amr_num_blocks
            if (amr_owns_all(b)) then
                amr_n_own = amr_n_own + 1
                amr_own_blk(amr_n_own) = b
            end if
            amr_parent_blk(b) = 0
            if (amr_block_level(b) >= 2) then
                pblk = f_amr_parent_block(b)
                amr_parent_blk(b) = pblk
                ! pblk == 0 is impossible on a nested mesh, but counting it would land in amr_child_ptr(0) and corrupt the lists
                if (pblk > 0) then
                    amr_child_ptr(pblk) = amr_child_ptr(pblk) + 1
                    if (amr_block_owner(pblk) == proc_rank .and. amr_block_owner(b) /= proc_rank) then
                        amr_n_fch = amr_n_fch + 1
                        amr_fch_blk(amr_n_fch) = b
                    end if
                end if
                cycle
            end if
            if (amr_block_level(b) /= 1) cycle
            if (amr_block_owner(b) == proc_rank) cycle
            call s_amr_region_box(b, rlo, rhi)
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (all(bl <= bh)) then
                amr_n_l1r = amr_n_l1r + 1
                amr_l1r_blk(amr_n_l1r) = b
            end if
            call s_amr_patch_box(b, plo, phi)
            call s_amr_box_isect(plo, phi, crlo, crhi, pl, ph)
            if (all(pl <= ph)) then
                amr_n_l1p = amr_n_l1p + 1
                amr_l1p_blk(amr_n_l1p) = b
            end if
        end do
        ! counts -> exclusive prefix, then the ascending fill bumps each amr_child_ptr(p) back to p's inclusive end
        acc = 0
        do b = 1, amr_num_blocks
            nxt = amr_child_ptr(b); amr_child_ptr(b) = acc; acc = acc + nxt
        end do
        do b = 1, amr_num_blocks
            if (amr_block_level(b) < 2) cycle
            pblk = amr_parent_blk(b)
            if (pblk == 0) cycle  ! uncounted above; a bump here would clobber a real entry
            amr_child_ptr(pblk) = amr_child_ptr(pblk) + 1
            amr_child_idx(amr_child_ptr(pblk)) = b
        end do
        amr_l1r_epoch = amr_mesh_epoch
        amr_l1r_nblk = amr_num_blocks

    end subroutine s_amr_refresh_lists

    !> Coarse-side overlap per rank of block k's region grown by g coarse cells: level 1 against the level-0 rank bricks, level >= 2
    !! against the owners of the level-(l-1) blocks. Adds into vol(0:) and lists touched ranks in rk(1:nr); the caller zeroes
    !! vol(rk(1:nr)) after use.
    impure subroutine s_amr_place_overlap(k, g, vol, rk, nr)

        integer, intent(in)     :: k, g
        real(wp), intent(inout) :: vol(0:)
        integer, intent(out)    :: rk(:), nr
        integer                 :: lo(3), hi(3), ilo(3), ihi(3), ol(3), oh(3), j, q

        lo = amr_region_lo_all(:,k) - merge(g, 0, amr_dim); hi = amr_region_hi_all(:,k) + merge(g, 0, amr_dim)
        nr = 0
        if (amr_block_level(k) == 1) then
            call s_amr_ranks_overlapping(lo, hi, rk, nr)
            do j = 1, nr
                call s_amr_rank_interior(rk(j), ilo, ihi)
                call s_amr_box_isect(lo, hi, ilo, ihi, ol, oh)
                vol(rk(j)) = product(real(max(oh - ol + 1, 0), wp))
            end do
            return
        end if
        do j = 1, amr_num_blocks
            if (amr_block_level(j) /= amr_block_level(k) - 1) cycle
            call s_amr_box_isect(lo, hi, amr_region_lo_all(:,j), amr_region_hi_all(:,j), ol, oh)
            if (any(ol > oh)) cycle
            q = amr_block_owner(j)
            if (vol(q) <= 0._wp) then
                nr = nr + 1; rk(nr) = q
            end if
            vol(q) = vol(q) + product(real(oh - ol + 1, wp))
        end do

    end subroutine s_amr_place_overlap

    !> Owner-aware placement (amr_place_owner) of one level's na blocks (slots aidx, weights wt, Morton keys key): each block goes
    !! to the rank owning most of the coarse cells under it. A rank above (1 + place_tol) x the mean sheds its last blocks in Morton
    !! order, first to the least-loaded rank next to the block, then to the least-loaded rank anywhere, never pushing the receiver
    !! above the cap. Replicated inputs in a fixed order, so every rank computes the same owners.
    impure subroutine s_amr_place_owner(key, wt, aidx, na, own)

        integer, intent(in)         :: na, aidx(na)
        integer(kind=8), intent(in) :: key(na)
        real(wp), intent(in)        :: wt(na)
        integer, intent(inout)      :: own(na)
        real(wp), parameter         :: place_tol = 0.05_wp
        real(wp), allocatable       :: vol(:), load(:)
        integer, allocatable        :: rk(:), ord(:)
        integer                     :: a, b, j, r, m, nr, pass
        real(wp)                    :: cap

        allocate (vol(0:num_procs - 1), load(0:num_procs - 1), rk(num_procs), ord(na))
        vol = 0._wp; load = 0._wp
        do a = 1, na
            call s_amr_place_overlap(aidx(a), 0, vol, rk, nr)
            do j = 1, nr  ! nr = 0 (no coarse owner found) keeps the SFC owner
                if (j == 1) own(a) = rk(1)
                if (vol(rk(j)) > vol(own(a)) .or. (vol(rk(j)) == vol(own(a)) .and. rk(j) < own(a))) own(a) = rk(j)
            end do
            vol(rk(1:nr)) = 0._wp
            load(own(a)) = load(own(a)) + wt(a)
        end do

        cap = (1._wp + place_tol)*sum(load)/real(num_procs, wp)
        call s_amr_sort_by_key(key, na, ord)
        do pass = 1, 2  ! 1: to a coarse neighbour, 2: anywhere
            do r = 0, num_procs - 1
                do b = na, 1, -1
                    if (load(r) <= cap) exit
                    a = ord(b)
                    if (own(a) /= r) cycle
                    if (pass == 1) then
                        call s_amr_place_overlap(aidx(a), 1, vol, rk, nr)
                        vol(rk(1:nr)) = 0._wp
                        m = -1
                        do j = 1, nr
                            if (rk(j) == r) cycle
                            if (m < 0) then
                                m = rk(j)
                            else if (load(rk(j)) < load(m) .or. (load(rk(j)) == load(m) .and. rk(j) < m)) then
                                m = rk(j)
                            end if
                        end do
                    else
                        m = minloc(load, 1) - 1
                    end if
                    if (m < 0) cycle
                    if (load(m) + wt(a) > cap) cycle
                    own(a) = m; load(r) = load(r) - wt(a); load(m) = load(m) + wt(a)
                end do
            end do
        end do
        deallocate (vol, load, rk, ord)

    end subroutine s_amr_place_owner

    impure subroutine s_amr_assign_block_owners()

        integer :: k, a, lev, maxlev, na
        ! heap, not stack: O(global boxes)
        integer, allocatable         :: aidx(:), aown(:)
        integer(kind=8), allocatable :: key(:), akey(:)
        real(wp), allocatable        :: wt(:), cost(:), awt(:)

        if (amr_num_blocks < 1) return

        allocate (aidx(amr_num_blocks), aown(amr_num_blocks), key(amr_num_blocks), akey(amr_num_blocks), wt(amr_num_blocks), &
                  & cost(amr_num_blocks), awt(amr_num_blocks))

        call s_amr_block_cost(cost)

        ! per-block fine-work weight = footprint cost x amr_ref_ratio**(level*active dims)
        do k = 1, amr_num_blocks
            wt(k) = cost(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (n_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (p_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            key(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
        end do

        ! Per-level distribution: one cut per level on each block's own weight, so a deep tower does not pin its subtree to one
        ! rank. Not one mixed cut: same-level boxes have distinct Morton keys, which f_amr_owner's cut search needs. amr_fine_cut(:,
        ! lev) holds each level's cut.
        maxlev = maxval(amr_block_level(1:amr_num_blocks))
        do lev = 1, maxlev
            na = 0
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                na = na + 1
                akey(na) = key(k); awt(na) = wt(k); aidx(na) = k
            end do
            if (na < 1) cycle
            call s_amr_sfc_cut(akey, awt, na, amr_fine_cut(:,lev), aown)
            if (amr_place_owner) call s_amr_place_owner(akey, awt, aidx, na, aown)
            do a = 1, na
                amr_block_owner(aidx(a)) = aown(a); amr_myblk_dirty = .true.
            end do
        end do

        call s_amr_validate_owner()
        call s_amr_report_balance(wt, maxlev)

        deallocate (aidx, aown, key, akey, wt, cost, awt)

    end subroutine s_amr_assign_block_owners

    !> Install a recorded owner table instead of recomputing it (restart at the writer's rank count): the cost weights read live
    !! state (body markers) that differs between the regrid before the write and the restart. Each level's cut is the running max
    !! Morton key of each rank's blocks; s_amr_validate_owner rejects a table not contiguous in that order (unless amr_place_owner).
    impure subroutine s_amr_adopt_block_owners(own)

        integer, intent(in) :: own(:)
        integer             :: k, r, lev, maxlev
        integer(kind=8)     :: key

        maxlev = maxval(amr_block_level(1:amr_num_blocks))
        do lev = 1, maxlev
            amr_fine_cut(:,lev) = -1_8
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                key = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
                amr_fine_cut(own(k), lev) = max(amr_fine_cut(own(k), lev), key)
                amr_block_owner(k) = own(k); amr_myblk_dirty = .true.
            end do
            do r = 1, num_procs - 1
                if (amr_fine_cut(r, lev) < amr_fine_cut(r - 1, lev)) amr_fine_cut(r, lev) = amr_fine_cut(r - 1, lev)
            end do
        end do
        call s_amr_validate_owner()

    end subroutine s_amr_adopt_block_owners

    !> Per-level and total load-balance report (rank 0, rank_time_wrt): max/mean assigned block weight over ranks (1 = perfect).
    !! no_blocks_ranks is the granularity floor of a level with fewer boxes than ranks, not idleness.
    impure subroutine s_amr_report_balance(wt, maxlev)

        real(wp), intent(in)  :: wt(:)
        integer, intent(in)   :: maxlev
        real(wp), allocatable :: rw(:), tw(:), rc(:)
        real(wp)              :: mx, mean, cmx, cmean
        integer               :: k, lev, nb, empty

        if (.not. rank_time_wrt .or. proc_rank /= 0) return

        ! heap: num_procs long
        allocate (rw(0:num_procs - 1), tw(0:num_procs - 1), rc(0:num_procs - 1))
        tw = 0._wp
        do lev = 1, maxlev
            rw = 0._wp; rc = 0._wp
            nb = 0
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                rw(amr_block_owner(k)) = rw(amr_block_owner(k)) + wt(k)
                rc(amr_block_owner(k)) = rc(amr_block_owner(k)) + 1._wp
                nb = nb + 1
            end do
            if (nb == 0) cycle
            tw = tw + rw
            mx = maxval(rw); mean = sum(rw)/real(num_procs, wp)
            empty = count(rw <= 0._wp)
            ! box-count imbalance beside weight: a block's advance cost has a large fixed component
            cmx = maxval(rc); cmean = sum(rc)/real(num_procs, wp)
            if (mean > 0._wp) print '(A,I0,A,I0,A,I0,A,F8.3,A,F8.3,A,I0,A,I0)', ' [amr-balance] level ', lev, ': boxes ', nb, &
                & '/ranks ', num_procs, ' max/mean ', mx/mean, ' boxes_max/mean ', cmx/cmean, ' no_blocks_ranks ', empty, ' of ', &
                & num_procs
        end do
        ! per-rank weight (fine cells)
        write (*, '(A)', advance='no') ' [amr-balance] per-rank fine_work :'
        do k = 0, num_procs - 1; write (*, '(I12)', advance='no') nint(tw(k), kind=8); end do
        write (*, '(A)') ''
        mean = sum(tw)/real(num_procs, wp)
        ! fine_work = sum of assigned weights = fine cells advanced per step
        if (mean > 0._wp) print '(A,F8.3,A,I0,A,I0)', ' [amr-balance] TOTAL   : max/mean ', maxval(tw)/mean, &
            & ' ranks_with_no_fine_block ', count(tw <= 0._wp), ' fine_work ', nint(sum(tw), kind=8)
        deallocate (rw, tw, rc)

    end subroutine s_amr_report_balance

    !> Owner rank of block k by binary search of its Morton key in its level's cut amr_fine_cut(:, level); reproduces
    !! s_amr_assign_block_owners' assignment exactly.
    pure integer function f_amr_owner(k) result(r)

        integer, intent(in) :: k
        integer(kind=8)     :: mk, cut(0:num_procs - 1)
        integer             :: lo, hi, mid

        cut = amr_fine_cut(:,amr_block_level(k))
        mk = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
        lo = 0; hi = num_procs - 1
        do while (lo < hi)
            mid = (lo + hi)/2
            if (mk <= cut(mid)) then
                hi = mid
            else
                lo = mid + 1
            end if
        end do
        r = lo

    end function f_amr_owner

    !> Guard: the SFC cut-point accessor must reproduce the stored owner table exactly.
    impure subroutine s_amr_validate_owner()

        integer :: k

        do k = 1, amr_num_blocks
            if (amr_place_owner) cycle  ! placed owners are not SFC-contiguous
            if (f_amr_owner(k) /= amr_block_owner(k)) &
                & call s_mpi_abort('SFC cut-point owner disagrees with amr_block_owner - cut capture or search is wrong')
        end do

    end subroutine s_amr_validate_owner

    !> Set this rank's intersection of box lo:hi with its subdomain (amr_isect_lo/hi) and whether it holds fine cells
    !! (amr_rank_owns_block: nonempty in all active dims). Call with the coarse grid state in m/n/p, never inside the fine advance.
    impure subroutine s_amr_compute_isect(lo, hi)

        integer, intent(in) :: lo(3), hi(3)

        amr_isect_lo = max(lo, amr_sidx); amr_isect_hi = min(hi, amr_sidx + amr_ext)
        amr_rank_owns_block = all(amr_isect_lo <= amr_isect_hi .or. .not. amr_dim)

    end subroutine s_amr_compute_isect

    !> Rank r's coarse interior box (global, no ghosts) from s_amr_rank_decomp; restriction targets are identified by interior
    !! overlap alone.
    pure subroutine s_amr_rank_interior(r, ilo, ihi)

        integer, intent(in)  :: r
        integer, intent(out) :: ilo(3), ihi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        ilo = merge(sidx, 0, amr_dim); ihi = merge(sidx + ext, 0, amr_dim)

    end subroutine s_amr_rank_interior

    !> Margin-padded global coarse-index bounding box of immersed body i (supported analytic geometries only; aborts on others).
    !! Reads the body's live centroid, so a moving body's box tracks it.
    impure subroutine s_amr_body_bbox(i, mrg, blo, bhi)

        integer, intent(in)  :: i, mrg
        integer, intent(out) :: blo(3), bhi(3)
        real(wp)             :: c(3), half(3)

        c = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, patch_ib(i)%z_centroid]
        select case (patch_ib(i)%geometry)
        case (2, 8, 10)  ! circle, sphere, cylinder: radius-bounded (cylinder length adds below)
            half = patch_ib(i)%radius
            if (patch_ib(i)%geometry == 10) then
                half(1) = max(half(1), 0.5_wp*patch_ib(i)%length_x)
                half(2) = max(half(2), 0.5_wp*patch_ib(i)%length_y)
                half(3) = max(half(3), 0.5_wp*patch_ib(i)%length_z)
            end if
        case (3, 9)  ! rectangle, box
            half = 0.5_wp*[patch_ib(i)%length_x, patch_ib(i)%length_y, patch_ib(i)%length_z]
        case default
            call s_mpi_abort('amr dynamic regrid with ib: unsupported body geometry for the ' &
                             & // 'containment bounding box (supported: circle/rectangle/sphere/box/cylinder)')
        end select
        call s_amr_phys_to_cells(c - half, c + half, mrg, blo, bhi)

    end subroutine s_amr_body_bbox

    !> Expand a candidate regrid box (global indices) to fully contain every immersed body it overlaps, with a margin; aborts if a
    !! body exceeds the per-rank block cap.
    impure subroutine s_amr_expand_box_over_bodies(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: i, blo(3), bhi(3), mrg

        ! margin: the IB image-point stencil reaches a few cells past the surface (buff_size would exceed the block cap). For
        ! amr_max_level > 1 add amr_cpat_mar per nesting level so the finest C/F boundary stays a stencil off the surface.

        mrg = max(amr_buf, 4) + max(0, amr_max_level - 1)*amr_cpat_mar

        do i = 1, num_ibs
            call s_amr_body_bbox(i, mrg, blo, bhi)
            ! blocks stay buff_size inside the domain: fail with a named message instead of a clipped body
            if (any(amr_dim .and. (blo < buff_size .or. bhi > [m_glb, n_glb, p_glb] - buff_size))) then
                call s_mpi_abort('amr dynamic regrid with ib: the immersed body plus its containment ' &
                                 & // 'margin does not fit inside the refinable domain interior (blocks stay buff_size off the edges)')
            end if
            if (.not. f_amr_boxes_overlap(lo, hi, blo, bhi)) cycle
            lo = min(lo, blo); hi = max(hi, bhi)
            if (any(amr_dim .and. hi - lo + 1 > amr_maxc_fit)) then
                call s_mpi_abort('amr dynamic regrid with ib: containing the immersed body plus margin ' &
                                 & // 'exceeds the per-rank block size cap; use fewer ranks or a larger amr_maxc_fit')
            end if
        end do

    end subroutine s_amr_expand_box_over_bodies

    !> Split box [lo:hi] into an even grid of adjacent sub-boxes each <= amr_maxc_fit per dim (the max a rank can whole-own),
    !! appended to out(nt+1:); sets capped=1 and stops if the cap is hit. Collapsed dims stay [0:0].
    pure subroutine s_amr_tile_box(lo, hi, out, nt, cap, capped, tsz)

        integer, intent(in)        :: lo(3), hi(3), cap
        type(t_box), intent(inout) :: out(:)
        integer, intent(inout)     :: nt, capped
        !> per-dim tile size (default amr_maxc_fit)
        integer, intent(in), optional :: tsz(3)
        integer                       :: ntl(3), s(3), t1, t2, t3, qlo(3), qhi(3), tc(3)

        tc = amr_maxc_fit; if (present(tsz)) tc = tsz
        tc = max(tc, 1)  ! a 0 tile size would divide by zero; 1 cell is the floor
        ntl = (hi - lo + tc)/tc; s = (hi - lo + ntl)/ntl  ! a collapsed dim ([0:0]) tiles once with size 1
        do t3 = 0, ntl(3) - 1
            qlo(3) = lo(3) + t3*s(3); qhi(3) = min(lo(3) + (t3 + 1)*s(3) - 1, hi(3))
            do t2 = 0, ntl(2) - 1
                qlo(2) = lo(2) + t2*s(2); qhi(2) = min(lo(2) + (t2 + 1)*s(2) - 1, hi(2))
                do t1 = 0, ntl(1) - 1
                    if (nt >= cap) then; capped = 1; return; end if
                    qlo(1) = lo(1) + t1*s(1); qhi(1) = min(lo(1) + (t1 + 1)*s(1) - 1, hi(1))
                    nt = nt + 1; out(nt)%lo = qlo; out(nt)%hi = qhi
                end do
            end do
        end do

    end subroutine s_amr_tile_box

end module m_amr_distribution
