#!/usr/bin/env python3
"""Generate a per-block cost-weights file for MFC_SFC_WEIGHTS (weighted SFC partition).

Per-cell cost model, following PR MFlowCode/MFC#1628's m_load_weight and adding the cost
of the fused reconstruction kernels:

    w = base * (1 + K*I_stencil + K_bub*void)

  base        1 everywhere, or --box G: 1 inside the grown active box around non-ambient
              cells and 0 outside (the #1628 active-box membership indicator; pairs with
              the runtime MFC_ACTIVE_BOX window, which makes skipped cells cost ~0).
  I_stencil   1 where any five-point alpha stencil along x, y or z is nonconstant. The
              fused reconstruction kernels shortcut constant stencils, so only these
              cells carry the reconstruction cost.
  void        EL-bubble void fraction (1 - alpha_liquid), the per-cell bubble-ODE count
              proxy of #1628's K_bub term; --kbub only, for bubbles_lagrange cases.

Usage:
  gen_sfc_weights.py <restart_data_dir> <out.weights> <K> [block] --grid G1,G2,G3
                     --alpha-var N [--box G] [--kbub KBUB] [--sample DIR ...]

Reads restart_data/lustre_0.dat (parallel-IO output, var-major, x fastest, float64) and
writes the file format consumed by s_sfc_read_weights in src/common/m_mpi_common.fpp:

  header: g1 g2 g3 block snap1 snap2 snap3
  body:   nblk(1)*nblk(2)*nblk(3) block weight sums, x fastest

K is the stencil-cost coefficient: the excess fused-work cost of a stencil-nonconstant
cell over a constant one. Measure it as the work ratio W0/W1 of two ranks with cell
counts C0, C1 and stencil-nonconstant counts E0, E1:
K = (W0*C1 - W1*C0) / (W1*E0 - W0*E1), e.g. from a two-rank split of the grid.

The partition is fixed for the whole run, so the quantity it should balance is the
time-integrated cost. Pass restart dumps sampled across the run as repeated --sample
DIR (the positional dir is one more sample). The per-cell weight fields are averaged
over the samples with equal weight. Without --sample the weights come from the single
positional directory.
"""

import argparse
import sys

import numpy as np


def nonconstant_stencil(a, radius=2):
    """Cells whose (2*radius+1)-point stencil along any axis varies (edge-clamped)."""
    out = np.zeros_like(a, dtype=bool)
    for axis in range(3):
        if a.shape[axis] <= 1:
            continue  # degenerate axis (1-D/2-D grids): no stencil to vary along it
        ap = np.pad(a, radius, mode="edge")
        av = np.moveaxis(ap, axis, 0)
        nc = np.zeros(av.shape, dtype=bool)
        for off in range(-radius, radius + 1):
            if off != 0:
                nc |= av != np.roll(av, -off, axis=0)
        res = np.moveaxis(nc, 0, axis)
        sl = tuple(slice(radius, s - radius) for s in res.shape)
        out |= res[sl]
    return out


def active_box(a, growth, tol=1e-12, floor_frac=1e-4):
    """Grown bounding box of cells deviating (relative tol) from the far-corner ambient state,
    matching the runtime window of m_active_box.fpp (PR MFlowCode/MFC#1628). Components whose
    ambient value is ~0 (quiescent momentum) borrow the largest ambient component's scale, as
    in the runtime predicate -- otherwise their threshold collapses to ~1e-22 absolute and
    every computed cell's roundoff drift registers as activity."""
    amb = a[:, -1, -1, -1]
    thresh = (tol*np.maximum(np.abs(amb), floor_frac*np.abs(amb).max()))[:, None, None, None]
    dev = np.any(np.abs(a - amb[:, None, None, None]) > thresh, axis=0)
    idx = np.argwhere(dev)
    if idx.size == 0:
        return None, 0.0
    lo = np.maximum(idx.min(axis=0) - growth, 0)
    hi = np.minimum(idx.max(axis=0) + growth, np.array(a.shape[1:]) - 1)
    box = np.zeros(a.shape[1:], dtype=bool)
    box[lo[0] : hi[0] + 1, lo[1] : hi[1] + 1, lo[2] : hi[2] + 1] = True
    return box, box.mean()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("restart_data_dir")
    p.add_argument("out")
    p.add_argument("K", type=float)
    p.add_argument("block", nargs="?", type=int, default=3)
    p.add_argument("--box", type=int, default=None, metavar="G",
                   help="active-box base: weight 1 inside the box grown by G cells, 0 outside")
    p.add_argument("--kbub", type=float, default=0.0, metavar="KBUB",
                   help="PR #1628 K_bub coefficient on the EL-bubble void fraction")
    p.add_argument("--sample", action="append", default=[], metavar="DIR",
                   help="additional restart_data dir to time-average over (repeatable); "
                        "the positional dir is one more sample")
    p.add_argument("--grid", required=True, metavar="G1,G2,G3",
                   help="cell counts of the restart file's grid")
    p.add_argument("--alpha-var", type=int, required=True, metavar="N",
                   help="1-based variable index of alpha(1) in the restart file")
    args = p.parse_args()
    rdirs = [args.restart_data_dir] + args.sample
    out, block, K, growth, kbub = (
        args.out, args.block, args.K, args.box, args.kbub
    )

    g = tuple(int(v) for v in args.grid.split(","))
    w = np.zeros((g[2], g[1], g[0]), dtype=float)  # [z, y, x], the restart reader's layout
    for i, rdir in enumerate(rdirs):
        arr = np.fromfile(f"{rdir}/lustre_0.dat", dtype="<f8")
        nvar = arr.size // (g[0] * g[1] * g[2])
        assert nvar * g[0] * g[1] * g[2] == arr.size, "unexpected lustre_0.dat layout"
        arr = arr.reshape(nvar, g[2], g[1], g[0])  # [var, z, y, x]

        alpha1 = arr[args.alpha_var - 1]  # alpha(1)
        ex = nonconstant_stencil(alpha1)

        print(f"[{rdir}] stencil-nonconstant cells: {ex.sum()} "
              f"({100 * ex.sum() / np.prod(g):.2f}%), K = {K:.6f}")

        base = np.ones((g[2], g[1], g[0]), dtype=float)
        if growth is not None:
            box, frac = active_box(arr, growth)
            if box is None:
                sys.exit("no non-ambient cells: an all-zero weights file is meaningless")
            base = box.astype(float)
            print(f"[{rdir}] active box (growth {growth}): {100 * frac:.2f}% of cells")

        # PR #1628 m_load_weight: EL-bubble void proxy. The gas volume fraction of a resolved
        # two-fluid region is NOT bubble-ODE cost, so kbub models cost only for
        # bubbles_lagrange cases.
        w += base * (1.0 + K * ex + kbub * (1.0 - alpha1))
    w /= len(rdirs)

    # Per-axis block counts with zero padding to full chunks, so degenerate axes
    # (1-D/2-D grids) block cleanly.
    nbs = [(s_ + block - 1) // block for s_ in (g[2], g[1], g[0])]      # (z, y, x) blocks
    chs = [(s_ + n_ - 1) // n_ for s_, n_ in zip((g[2], g[1], g[0]), nbs)]
    pad = np.zeros((nbs[0]*chs[0], nbs[1]*chs[1], nbs[2]*chs[2]))
    pad[:g[2], :g[1], :g[0]] = w
    blk = pad.reshape(nbs[0], chs[0], nbs[1], chs[1], nbs[2], chs[2]).sum(axis=(1, 3, 5))
    wblk = blk.transpose(2, 1, 0).reshape(-1)  # x fastest

    with open(out, "w") as f:
        f.write(f"{g[0]} {g[1]} {g[2]} {block} 1 1 1\n")
        f.write("\n".join(f"{v:.17g}" for v in wblk))
        f.write("\n")
    print(f"wrote {out}: {tuple(reversed(nbs))} blocks of {block}^3 cells, total weight {wblk.sum():.0f}"
          + (f", averaged over {len(rdirs)} samples" if len(rdirs) > 1 else ""))


if __name__ == "__main__":
    main()
