# Rising-bubble benchmark reference data

From the FeatFlow benchmark page (http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html),
accompanying Hysing et al., "Quantitative benchmark computations of two-dimensional bubble
dynamics", International Journal for Numerical Methods in Fluids 60:1259-1288 (2009). Names follow the source, c#g#l#:
test case, group, refinement level. Each group's finest level is kept:

| file         | group             |
|--------------|-------------------|
| c1g1l7, c2g1l8 | TU Dortmund (TP2D) |
| c1g2l3, c2g2l3 | EPFL Lausanne (FreeLIFE) |
| c1g3l4, c2g3l4 | Uni Magdeburg (MooNMD) |

`c#g#l#.txt` columns: t, bubble area, circularity, centroid height y_c, rise velocity v_c.
Thinned to one row per 0.005 in t (the source stores every time step).

`c#g#l#s.txt`: bubble contour points (x, y) at t = 3, thinned to at most 2000 points.
