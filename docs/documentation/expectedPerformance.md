# Performance Results

MFC has been extensively benchmarked on both CPUs and GPUs. A summary of these results follow.

## Expected time-steps/hour

The following table outlines expected performance in terms of number of time-steps per hour 
(rounded to the nearest hundred) for various problem sizes and hardware for a inviscid, 6-equation,
3D simulation. CPU results utilize an entire die.

| Hardware             | # Ranks | 1M Cells       | 4M Cells       | 8M Cells     | Compiler    | Computer      |
| ---:                 | :----:  |    :----:      |  :---:         | :---:        | :----:      | :---          |
| Nvidia V100          | 1       | 88.5k          | 18.7k          | N/A          | NVHPC 22.11 | PACE Phoenix  |
| Nvidia A100          | 1       | 114.4k         | 34.6k          | 16.5k        | NVHPC 23.5  | Wingtip       |
| AMD MI250x           | 1       | 77.5k          | 22.3k          | 11.2k        | CCE 16.0.1  | OLCF Frontier |
| Intel Xeon Gold 6226 | 12      | 2.5k           | 0.7k           | 0.4k         | GNU 10.3.0  | Pace Phoenix  |
| Apple Silicon M2     | 6       | 2.8k           | 0.6k           | 0.2k         | GNU 13.2.0  | N/A           |

If `'model_eqns' : 3` is replaced by `'model_eqns' : 2`, an inviscid 5-equation model is used.
The following table outlines expected performance in terms of number of time-steps per hour 
(rounded to the nearest hundred) for various problem sizes and hardware for a inviscid, 5-equation,
3D simulation. CPU results utilize an entire die.

| Hardware             | # Ranks | 1M Cells       | 4M Cells       | 8M Cells     | Compiler    | Computer      |
| ---:                 | :----:  |    :----:      |  :---:         | :---:        | :----:      | :---          |
| Nvidia V100          | 1       | 113.4k         | 26.2k          | N/A          | NVHPC 22.11 | PACE Phoenix  |
| Nvidia A100          | 1       | 153.5k         | 48.0k          | 22.5k        | NVHPC 23.5  | Wingtip       |
| AMD MI250x           | 1       | 104.2k         | 31.0k          | 14.8k        | CCE 16.0.1  | OLCF Frontier |
| Intel Xeon Gold 6226 | 12      | 5.4k           | 1.6k           | 0.8k         | GNU 10.3.0  | Pace Phoenix  |
| Apple Silicon M2     | 6       | 3.7k           | 11.0k          | 0.3k         | GNU 13.2.0  | N/A           |

## Weak scaling

Strong scaling results are obtained by increasing the problem size with the number of processes
so that work per process remains constant.

### AMD MI250X GPU
MFC weask scales to 65,536 AMD MI250X GPUs on OLCF Frontier with 96% efficiency.

<img src="../res/weakScaling/frontier.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

### Nvidia V100 GPU
MFC weak scales to 13,824 V100 Nvidia V100 GPUs on OLCF Summit with 97% efficiency.

<img src="../res/weakScaling/summit.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

### IMB Power9 CPU
MFC Weak scales to 13,824 Power9 CPU cores on OLCF Summit with 1% of ideal scaling.

<img src="../res/weakScaling/cpuScaling.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

## Strong scaling

Strong scaling results are obtained by keeping the problem size constant and increasing
the number of process so that work per process decreases.

### Nvidia V100 GPU

For these tests, the base case utilizes 8 GPUs with one MPI process per GPU. The performance
is analyzed at two different problem sizes of 16 and 64M grid points, with the base case using
2 and 8M grid points per process.

#### 16M Grid Points
<img src="../res/strongScaling/strongScaling16.svg" style="width: 50%; border-radius: 10pt"/>

#### 64M Grid Points
<img src="../res/strongScaling/strongScaling64.svg" style="width: 50%; border-radius: 10pt"/>

### IBM Power 9 CPU

CPU strong scaling tests are done with problem sizes of 16, 32, and 64M grid points, with the
base case using 2, 4, and 8M cells per process.

<img src="../res/strongScaling/cpuStrongScaling.svg" style="width: 50%; border-radius: 10pt"/>