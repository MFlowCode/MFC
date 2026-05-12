ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ARG TARGET
ARG CC_COMPILER
ARG CXX_COMPILER
ARG FC_COMPILER
ARG COMPILER_PATH
ARG COMPILER_LD_LIBRARY_PATH
ARG AFAR_VERSION
ARG OLCF_AFAR_ROOT=""

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update -y && \
    apt-get install -y software-properties-common ca-certificates gnupg wget && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update -y && \
    if [ "$TARGET" = "cpu" ]; then \
        apt-get install -y \
            build-essential git make cmake gcc g++ gfortran bc \
            python3.12 python3.12-venv python3-pip \
            openmpi-bin libopenmpi-dev libfftw3-dev \
            mpich libmpich-dev; \
    elif [ "$TARGET" = "gpu" ]; then \
        apt-get install -y \
            build-essential git make cmake bc \
            python3.12 python3.12-venv python3-pip \
            libfftw3-dev \
            openmpi-bin libopenmpi-dev; \
    elif [ "$TARGET" = "amd" ]; then \
        apt-get install -y \
            build-essential git make gcc g++ gfortran bc \
            python3.12 python3.12-venv python3-pip \
            libfftw3-dev libnuma1 libdrm2 libdrm-amdgpu1; \
    fi && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 2 && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# AMD: download AFAR, install cmake 3.28 (Ubuntu 22.04 ships 3.22 which doesn't
# recognize amdflang as LLVMFlang), then build MPICH with amdflang so the
# generated mpi.mod is compiler-compatible.
RUN if [ "$TARGET" = "amd" ] && [ -n "$AFAR_VERSION" ]; then \
        OLCF_AFAR_ROOT="/opt/${AFAR_VERSION}" && \
        wget -q "https://repo.radeon.com/rocm/misc/flang/${AFAR_VERSION}-ubuntu.tar.bz2" -O /tmp/afar.tar.bz2 && \
        tar -xjf /tmp/afar.tar.bz2 -C /opt/ && \
        rm /tmp/afar.tar.bz2 && \
        CMAKE_VER=3.28.6 && \
        wget -q "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-x86_64.sh" \
            -O /tmp/cmake-install.sh && \
        chmod +x /tmp/cmake-install.sh && \
        /tmp/cmake-install.sh --prefix=/usr/local --skip-license --exclude-subdir && \
        rm /tmp/cmake-install.sh && \
        printf '#!/bin/bash\nargs=()\nwhile [ "$#" -gt 0 ]; do\n    if [ "$1" = "-soname" ]; then\n        args+=("-Wl,-soname,$2"); shift 2\n    else\n        args+=("$1"); shift\n    fi\ndone\nexec '"${OLCF_AFAR_ROOT}"'/bin/amdflang "${args[@]}"\n' \
            > /usr/local/bin/amdflang-ld-wrap && \
        chmod +x /usr/local/bin/amdflang-ld-wrap && \
        MPICH_VER=3.4.3 && \
        wget -q "https://www.mpich.org/static/downloads/${MPICH_VER}/mpich-${MPICH_VER}.tar.gz" \
            -O /tmp/mpich.tar.gz && \
        mkdir -p /tmp/mpich-src && \
        tar -xzf /tmp/mpich.tar.gz -C /tmp/mpich-src --strip-components=1 && \
        cd /tmp/mpich-src && \
        FC=/usr/local/bin/amdflang-ld-wrap CC=gcc CXX=g++ \
            ./configure --prefix=/opt/mpich --enable-shared --disable-static \
            --with-device=ch3 2>&1 && \
        make -j$(nproc) 2>&1 && \
        make install 2>&1 && \
        cd / && \
        rm -rf /tmp/mpich-src /tmp/mpich.tar.gz; \
    fi

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV HYDRA_LAUNCHER=fork
ENV PATH="/opt/MFC:$PATH"

COPY ../ /opt/MFC

ENV CC=${CC_COMPILER}
ENV CXX=${CXX_COMPILER}
ENV FC=${FC_COMPILER}
ENV OLCF_AFAR_ROOT=${OLCF_AFAR_ROOT}
ENV PATH="${COMPILER_PATH}:/opt/mpich/bin:$PATH"
ENV LD_LIBRARY_PATH="${COMPILER_LD_LIBRARY_PATH}:/opt/mpich/lib:${LD_LIBRARY_PATH:-}"

# Pre-install numpy into the venv before mfc.sh runs, as it's required at
# build time by several dependencies (pandas, cantera, matplotlib, etc.) that
# may need to compile from source on architectures without pre-built wheels
RUN python3.12 -m venv /opt/MFC/build/venv && \
    /opt/MFC/build/venv/bin/pip install --upgrade pip && \
    /opt/MFC/build/venv/bin/pip install numpy

RUN echo "TARGET=$TARGET CC=$CC_COMPILER FC=$FC_COMPILER" && \
    cd /opt/MFC && \
    if [ "$TARGET" = "gpu" ]; then \
      ./mfc.sh build --gpu acc -j $(nproc); \
    elif [ "$TARGET" = "amd" ]; then \
      ./mfc.sh build --gpu mp -j $(nproc); \
    else \
      ./mfc.sh build -j $(nproc); \
    fi

RUN cd /opt/MFC && \
    if [ "$TARGET" = "gpu" ]; then \
      ./mfc.sh test -a --dry-run --gpu acc -j $(nproc); \
    elif [ "$TARGET" = "amd" ]; then \
      ./mfc.sh test -a --dry-run --gpu mp -j $(nproc); \
    else \
      ./mfc.sh test -a --dry-run -j $(nproc); \
    fi

WORKDIR /opt/MFC
ENTRYPOINT ["tail", "-f", "/dev/null"]