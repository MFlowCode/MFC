#!/usr/bin/env python3
"""
Binary STL of the flat ring r_in <= r <= r_out (counter-clockwise triangles, z = 0): the solid outside an annular gas gap.

    ring_stl.py OUT.stl R_IN R_OUT [NSEG]

`ring.stl` in this directory is `ring_stl.py ring.stl 2 5 512` (lengths in units of the inner radius R). Binary, 50 bytes per
triangle, so the file stays small; 512 segments put the chord error of the inner edge below 4e-5 R.
"""

import math
import struct
import sys


def write_ring(path, r_in, r_out, nseg=512):
    pts = [[(r * math.cos(t), r * math.sin(t), 0.0) for t in (2.0 * math.pi * k / nseg for k in range(nseg))] for r in (r_in, r_out)]
    tris = []
    for k in range(nseg):
        a, b, c, d = pts[0][k], pts[0][(k + 1) % nseg], pts[1][k], pts[1][(k + 1) % nseg]
        tris += [(a, c, d), (a, d, b)]  # both counter-clockwise
    with open(path, "wb") as f:
        f.write(b"binary STL ring".ljust(80, b" "))  # must not start with 'solid' (that selects the ASCII reader)
        f.write(struct.pack("<I", len(tris)))
        for tri in tris:
            f.write(struct.pack("<3f", 0.0, 0.0, 1.0))
            for v in tri:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


if __name__ == "__main__":
    write_ring(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 512)
