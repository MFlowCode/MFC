# Automatic `ib_neighborhood_radius` when ranks are not cubes

`ib_neighborhood_radius` is a count of **rank hops**. A rank keeps an immersed-boundary patch only while the
patch's centroid lies inside its own subdomain grown outward by that many hops
(`s_get_neighbor_bounds`, `f_neighborhood_ranks_own_location`), and drops it otherwise. So the radius has to be
large enough that the body is reachable within that many hops **in every direction** — and the hop that needs
the most is the one stepping across the *thinnest* rank.

When the radius is not set in the case file, MFC chooses it from the body's half-extent divided by a rank
width. The width it used was assembled the wrong way round:

```fortran
local_rank_width = -1._wp
do each direction:
    local_rank_width = max(local_rank_width, <this rank's extent in that direction>)
call s_mpi_allreduce_min(local_rank_width, min_rank_width)
```

Each rank reports its **widest** extent, and the minimum is taken over ranks. On a decomposition where every
rank is long in one direction and thin in another, the reported width is the long one, the radius comes out too
small, and ranks that should have kept the patch drop it.

## The case

A thin plate in a long, narrow channel: 20 chords by 8 by 6, with 400 x 50 x 50 cells chosen so MFC's topology
search settles on **16 x 2 x 2** at 64 ranks. That is an ordinary shape for a wake, a jet or a channel — the
flow direction resolved far more finely than the cross-stream ones — and it makes the ranks strongly
anisotropic:

| direction | ranks | extent per rank |
| --- | --- | --- |
| x | 16 | **1.250** |
| y | 2 | 4.000 |
| z | 2 | 3.000 |

The plate is 1.0 x 2.4 x 0.1, so `s_get_ib_bound` (geometry 9, the cuboid's half-diagonal) returns **1.3010**.
Crossing that at 1.250 per hop needs `ceil(1.1 * 1.3010 / 1.250) = 2` hops. The old width of 4.000 gives
`ceil(1.1 * 1.3010 / 4.000) = 1`.

## Running it

```
./mfc.sh run examples/3D_ibm_neighborhood_radius/case.py -n 64
```

and read the line MFC prints at start-up:

```
Automatic choice of ib_neighborhood_radius selected:  N
```

| | printed radius |
| --- | --- |
| before | **1** |
| after | **2** |

Both runs complete; the case is 1 M cells and takes a few minutes on two CPU nodes. `SUMMARY=1 python3
case.py` prints the half-extent, the rank extents and the arithmetic above without running anything.

## It is not only this case

The same two production grids that motivated the fix, measured from their own `lustre_*_cb.dat`:

| case | topology | old width | old radius | new width | new radius |
| --- | --- | --- | --- | --- | --- |
| gust encounter, 128 ranks | 16 x 2 x 4 | 2.051 | 1 | 1.052 | **2** |
| flapping wing, 128 ranks | 8 x 4 x 4 | 1.745 | 1 | 1.027 | **2** |

Both pick 1 where 2 is required.

## Scope

The new width is never larger than the old one, so the chosen radius never decreases: the change can only make
the neighbourhood more conservative, at the cost of more hops in the force reduction. Cases that set
`ib_neighborhood_radius` explicitly are untouched, and so is any decomposition whose ranks are close to cubic,
where the widest and narrowest extents coincide — which is why a uniform grid with a balanced topology shows no
difference.
