#!/usr/bin/env python3
"""
Write the binary STL of a flat 2D ring r_in <= r <= r_out (counter-clockwise triangles in z = 0): the solid outside a circular gas gap.

    ring_stl.py OUT.stl R_IN R_OUT [NSEG]

ring.stl in this directory is `ring_stl.py ring.stl 2 5 1024`.
"""

import math
import struct
import sys


def write_ring(path, r_in, r_out, nseg=1024):
    th = [2.0 * math.pi * k / nseg for k in range(nseg)]
    inner, outer = ([(r * math.cos(t), r * math.sin(t), 0.0) for t in th] for r in (r_in, r_out))
    tris = []
    for k in range(nseg):
        a, b, c, d = inner[k], inner[(k + 1) % nseg], outer[k], outer[(k + 1) % nseg]
        tris += [(a, c, d), (a, d, b)]
    with open(path, "wb") as f:
        f.write(b"ring".ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        for tri in tris:
            f.write(struct.pack("<12fH", 0.0, 0.0, 1.0, *(x for v in tri for x in v), 0))


if __name__ == "__main__":
    write_ring(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 1024)
