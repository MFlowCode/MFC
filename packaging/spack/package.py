# Copyright 2013-2025 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Mfc(CMakePackage):
    """MFC (Multicomponent Flow Code) is an exascale multiphase/multiphysics
    compressible flow solver. It scales ideally to 43K+ GPUs on leadership-class
    supercomputers and supports high-order WENO/TENO schemes, immersed boundary
    methods, phase change, surface tension, and MHD."""

    homepage = "https://mflowcode.github.io/"
    url = "https://github.com/MFlowCode/MFC/archive/refs/tags/v5.1.0.tar.gz"
    git = "https://github.com/MFlowCode/MFC.git"

    maintainers("sbryngelson")

    license("MIT")

    version("master", branch="master")
    version("5.1.0", sha256="4684bee6a529287f243f8929fb7edb0dfebbb04df7c1806459761c9a6c9261cf")

    # Build options
    variant("mpi", default=True, description="Build with MPI support")
    variant("openacc", default=False, description="Build with OpenACC GPU support")
    variant("openmp", default=False, description="Build with OpenMP GPU support")
    variant(
        "precision",
        default="double",
        values=("single", "double"),
        description="Floating point precision",
    )
    variant("post_process", default=True, description="Build post-processing tool")
    variant("chemistry", default=False, description="Enable thermochemistry via Pyrometheus/Cantera")

    # GPU architecture variants for performance optimization
    variant(
        "cuda_arch",
        default="none",
        values=("none", "60", "70", "75", "80", "90"),
        multi=False,
        description="CUDA architecture for NVIDIA GPUs (60=P100, 70=V100, 80=A100, 90=H100)"
    )
    variant(
        "amdgpu_target",
        default="none",
        values=("none", "gfx908", "gfx90a", "gfx940", "gfx941", "gfx942"),
        multi=False,
        description="AMD GPU architecture (gfx908=MI100, gfx90a=MI250X, gfx940/941/942=MI300)"
    )

    # Required dependencies
    depends_on("cmake@3.20:", type="build")
    depends_on("py-fypp", type="build")
    depends_on("python@3:", type="build")

    # Chemistry dependencies
    depends_on("cantera+python", type="build", when="+chemistry")
    # Note: py-pyrometheus may not be in Spack yet; will be added to PYTHONPATH via resource
    resource(
        name="pyrometheus",
        url="https://files.pythonhosted.org/packages/21/77/1e48bef25dfef5d9e35c1ab3a3a2ea1c82adb59aceb82b18d13b3d6c8a2b/pyrometheus-1.0.5.tar.gz",
        sha256="a572ab6db954f4a850d1292bb1ef6d6055916784a894d149d657996fa98d0367",
        when="+chemistry",
        placement="pydeps/pyrometheus",
        expand=True
    )

    # Runtime dependencies
    depends_on("fftw@3:", when="~mpi")
    depends_on("fftw@3:+mpi", when="+mpi")
    depends_on("lapack")

    # Optional dependencies
    depends_on("mpi", when="+mpi")
    depends_on("silo+hdf5", when="+post_process~mpi")
    depends_on("silo+hdf5+mpi", when="+post_process+mpi")

    # GPU dependencies
    depends_on("cuda", when="+openacc %nvhpc")
    depends_on("cuda", when="+openmp %nvhpc")
    depends_on("hip", when="+openacc %cce")
    depends_on("hip", when="+openmp %cce")

    # Compiler requirements
    conflicts("%gcc@:4.999", msg="MFC requires GCC 5.0 or newer")
    conflicts("%nvhpc@:21.6.999", msg="MFC requires NVHPC 21.7 or newer")
    conflicts("%apple-clang", msg="MFC does not support Apple Clang")
    conflicts("+openacc", when="%gcc", msg="OpenACC requires NVHPC or Cray compilers")
    conflicts(
        "+openacc", when="+openmp", msg="OpenACC and OpenMP GPU offload are mutually exclusive"
    )
    conflicts(
        "+openmp", when="+openacc", msg="OpenACC and OpenMP GPU offload are mutually exclusive"
    )

    def cmake_args(self):
        args = [
            self.define_from_variant("MFC_MPI", "mpi"),
            self.define_from_variant("MFC_OpenACC", "openacc"),
            self.define_from_variant("MFC_OpenMP", "openmp"),
            self.define("MFC_PRE_PROCESS", True),
            self.define("MFC_SIMULATION", True),
            self.define_from_variant("MFC_POST_PROCESS", "post_process"),
            self.define_from_variant("MFC_CHEMISTRY", "chemistry"),
        ]

        if self.spec.variants["precision"].value == "single":
            args.append(self.define("MFC_SINGLE_PRECISION", True))

        # GPU architecture targeting for optimal performance
        if "+openacc" in self.spec or "+openmp" in self.spec:
            # NVIDIA GPU architecture
            if self.spec.variants["cuda_arch"].value != "none":
                cuda_arch = self.spec.variants["cuda_arch"].value
                if self.spec.satisfies("%nvhpc"):
                    args.append(self.define("CMAKE_CUDA_ARCHITECTURES", cuda_arch))
                    # Add NVHPC-specific GPU flags
                    args.append(self.define("CMAKE_Fortran_FLAGS", 
                                          f"-gpu=cc{cuda_arch} -acc -Minfo=accel"))
            
            # AMD GPU architecture
            if self.spec.variants["amdgpu_target"].value != "none":
                amdgpu = self.spec.variants["amdgpu_target"].value
                if self.spec.satisfies("%cce") or self.spec.satisfies("%rocmcc"):
                    args.append(self.define("CMAKE_HIP_ARCHITECTURES", amdgpu))
                    # Add Cray/ROCm-specific GPU flags
                    args.append(self.define("CMAKE_Fortran_FLAGS",
                                          f"-fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march={amdgpu}"))

        return args

    def setup_build_environment(self, env):
        # Fypp is required for preprocessing
        env.prepend_path("PATH", self.spec["py-fypp"].prefix.bin)
        
        # Make vendored Pyrometheus importable when chemistry is enabled
        if "+chemistry" in self.spec:
            env.prepend_path("PYTHONPATH", os.path.join(self.stage.source_path, "pydeps", "pyrometheus"))
        
        # Cray system-specific configuration
        # On Cray systems, use compiler wrappers that automatically include system libraries
        if self.spec.satisfies("%cce") or self.spec.satisfies("^cray-mpich"):
            env.set("CC", "cc")
            env.set("CXX", "CC")
            env.set("FC", "ftn")
            # Cray wrappers handle MPI automatically
            if "+mpi" in self.spec:
                env.set("MPI_CC", "cc")
                env.set("MPI_CXX", "CC")
                env.set("MPI_FC", "ftn")
        
        # Enable GPU-aware MPI if available
        if "+mpi" in self.spec and ("+openacc" in self.spec or "+openmp" in self.spec):
            # NVIDIA GPU-aware MPI
            if self.spec.satisfies("^cuda"):
                env.set("MPICH_GPU_SUPPORT_ENABLED", "1")
                env.set("OMPI_MCA_opal_cuda_support", "true")
            # AMD GPU-aware MPI
            if self.spec.satisfies("^hip") or self.spec.satisfies("^rocm"):
                env.set("MPICH_GPU_SUPPORT_ENABLED", "1")
                env.set("HSA_ENABLE_SDMA", "0")  # Better performance for small messages
        
        # Set parallel build jobs based on available resources
        # HPC nodes often have many cores
        import multiprocessing
        env.set("CMAKE_BUILD_PARALLEL_LEVEL", str(min(multiprocessing.cpu_count(), 64)))
