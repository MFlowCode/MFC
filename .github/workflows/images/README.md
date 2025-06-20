# MFC Apptainer/Singularity Container Images

This directory contains Apptainer/Singularity recipe files for building containerized versions of MFC (Multi-phase Flow Code).

## Container Images

### 1. CPU Image (`Singularity.cpu`)
- Standard CPU-only build of MFC
- Ubuntu 24.04 base
- Includes all necessary dependencies for CPU simulations
- Optimized cache configuration

### 2. CPU Benchmark Image (`Singularity.cpu_bench`)
- CPU build with additional benchmarking tools
- Includes performance monitoring utilities (htop, iotop, sysstat, perf)
- Same base as CPU image with benchmarking focus

### 3. GPU Image (`Singularity.gpu`)
- GPU-enabled build with CUDA support
- Ubuntu 24.04 base with CUDA 12.6 toolkit
- Includes NVIDIA drivers and GPU libraries
- Larger cache size for GPU computations

### 4. GPU Benchmark Image (`Singularity.gpu_bench`)
- GPU build using NVIDIA HPC SDK
- Based on NVIDIA's official HPC container
- Includes NVIDIA compilers (nvc, nvc++, nvfortran)
- Pre-built with benchmarking suite
- Includes GPU profiling tools (nsight-systems, nsight-compute)

## Features

All images include:
- **Fakeroot support**: Allows running containers without root privileges
- **Cache configuration**: Optimized cache directories for better performance
- **MPI support**: OpenMPI for parallel computations
- **Pre-built MFC**: MFC is built during image creation for immediate use
- **Help documentation**: Use `--help` flag with any image for usage examples

## Building Images

### Using GitHub Actions (Automated)
Images are automatically built when changes are pushed to the repository.

### Building Locally
1. Install Apptainer: https://apptainer.org/docs/admin/main/installation.html

2. Enable fakeroot for your user:
   ```bash
   sudo apptainer config fakeroot --enable $(whoami)
   ```

3. Use the provided build script:
   ```bash
   cd .github/workflows/images
   ./build-local.sh
   ```

   Or build individual images:
   ```bash
   apptainer build --fakeroot mfc_cpu.sif Singularity.cpu
   apptainer build --fakeroot mfc_gpu.sif Singularity.gpu
   ```

## Using the Container Images

### CPU Image
```bash
# Run MFC with CPU image
apptainer run --fakeroot mfc_cpu.sif run examples/2D_shockbubble/case.py -n 4

# Run tests
apptainer run --fakeroot mfc_cpu.sif test -j 8

# Interactive shell
apptainer shell --fakeroot mfc_cpu.sif
```

### GPU Image
```bash
# Run MFC with GPU acceleration (note the --nv flag)
apptainer run --nv --fakeroot mfc_gpu.sif run examples/2D_shockbubble/case.py -n 4 --gpu

# GPU profiling
apptainer run --nv --fakeroot mfc_gpu.sif run case.py --nsys
apptainer run --nv --fakeroot mfc_gpu.sif run case.py --ncu
```

### Benchmark Images
```bash
# CPU benchmarking
apptainer run --fakeroot mfc_cpu_bench.sif bench -o bench.yaml

# GPU benchmarking
apptainer run --nv --fakeroot mfc_gpu_bench.sif bench -o bench.yaml
```

## Cache Configuration

All images are configured with optimized cache settings:
- **Apptainer cache**: `/tmp/apptainer-cache`
- **Singularity cache**: `/tmp/singularity-cache`
- **CUDA cache** (GPU images): `/tmp/cuda-cache`
- **NVIDIA compiler cache** (GPU bench): `/tmp/nvcompiler-cache`

Cache sizes:
- CPU images: Standard system cache
- GPU image: 1GB CUDA cache
- GPU benchmark: 2GB CUDA cache

## Mounting External Directories

To work with files outside the container:
```bash
# Mount current directory
apptainer run --fakeroot --bind $(pwd):/work mfc_cpu.sif run /work/case.py

# Mount multiple directories
apptainer run --fakeroot --bind /data:/data,/results:/results mfc_cpu.sif run case.py
```

## Troubleshooting

### Fakeroot Issues
If you encounter permission errors:
```bash
# Check if fakeroot is enabled
apptainer config fakeroot --show

# Enable for your user
sudo apptainer config fakeroot --enable $(whoami)
```

### GPU Not Detected
- Ensure NVIDIA drivers are installed on the host
- Use the `--nv` flag when running GPU containers
- Check GPU availability: `nvidia-smi`

### Cache Permission Errors
- Clear cache directories: `rm -rf /tmp/*-cache`
- Use `--no-cache` flag during build if needed

## Performance Tips

1. **Use appropriate image**: CPU for CPU-only systems, GPU for NVIDIA GPUs
2. **Bind mount for I/O**: Mount data directories to avoid copying large files
3. **Adjust cache size**: Modify cache environment variables for your workload
4. **Use benchmark images**: For performance testing and optimization