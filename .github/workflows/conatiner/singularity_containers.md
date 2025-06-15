# Singularity/Apptainer Containerization for MFC

This guide provides Singularity definition files for four scenarios:

- **CPU**: Standard CPU build  
- **GPU**: GPU-enabled build (NVIDIA HPC SDK)  
- **CPU Benchmark**: CPU build optimized for a specific benchmark  
- **GPU Benchmark**: GPU build optimized for a specific benchmark  

## Usage

```sh
# Build example (CPU)
sudo singularity build mfc_cpu.sif Singularity.cpu

# Run example
singularity exec mfc_cpu.sif ./mfc.sh --help
```

## Files

- `Singularity.cpu` — Standard CPU build
- `Singularity.gpu` — GPU-enabled build
- `Singularity.cpu_bench` — CPU benchmark-optimized
- `Singularity.gpu_bench` — GPU benchmark-optimized

---
