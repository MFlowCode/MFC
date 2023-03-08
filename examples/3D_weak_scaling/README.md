# 3D Weak Scaling

The [**3D_weak_scaling**](case.py) case depends on two parameters:

- **The number of MPI ranks** (_procs_): As _procs_ increases, the problem
size per rank remains constant. _procs_ is determined using information provided
to the case file by `mfc.sh run`.

- **GPU memory usage per rank** (_gbpp_): As _gbpp_ increases, the problem
size per rank increases and the number of timesteps decreases so that wall times
consistent. _gbpp_ is a user-defined optional argument to the [case.py](case.py) 
file. It can be specified right after the case filepath when invoking `mfc.sh run`.

Weak scaling benchmarks can be produced by keeping _gbpp_ constant and varying _procs_.

For example, to run a weak scaling test that uses ~4GB of GPU memory per rank
on 8 2-rank nodes with case optimization, one could:

```console
./mfc.sh run examples/3D_weak_scaling/case.py 4 -t pre_process simulation          \
             -e batch -p mypartition -N 8 -n 2 -w "01:00:00" -# "MFC Weak Scaling" \
             --case-optimization -j 32
```

