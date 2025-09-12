# Scaling and Performance test

The scaling case can exercise both weak- and strong-scaling. It
adjusts itself depending on the number of requested ranks.

This directory also contains a collection of scripts used to test strong and weak
scaling on OLCF Frontier.

## Weak Scaling

Pass `--scaling weak`. The `--memory` option controls (approximately) how much
memory each rank should use, in Gigabytes. The number of cells in each dimension
is then adjusted according to the number of requested ranks and an approximation
for the relation between cell count and memory usage. The problem size increases
linearly with the number of ranks.

## Strong Scaling

Pass `--scaling strong`. The `--memory` option controls (approximately) how much
memory should be used in total during simulation, across all ranks, in Gigabytes.
The problem size remains constant as the number of ranks increases.

## Example

For example, to run a weak-scaling test that uses ~4GB of GPU memory per rank
on 8 2-rank nodes with case optimization, one could:

```shell
./mfc.sh run examples/scaling/case.py -t pre_process simulation                    \
             -e batch -p mypartition -N 8 -n 2 -w "01:00:00" -# "MFC Weak Scaling" \
             --case-optimization -j 32 -- --scaling weak --memory 4
```
