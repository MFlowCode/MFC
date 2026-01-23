ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ARG TARGET
ARG CC_COMPILER
ARG CXX_COMPILER
ARG FC_COMPILER
ARG COMPILER_PATH
ARG COMPILER_LD_LIBRARY_PATH

RUN apt-get update -y && \
    if [ "$TARGET" != "gpu" ]; then \
        apt-get install -y \
            build-essential git make cmake gcc g++ gfortran bc\
            python3 python3-venv python3-pip \
            openmpi-bin libopenmpi-dev libfftw3-dev \
            mpich libmpich-dev; \
    else \
        apt-get install -y \
            build-essential git make cmake bc\
            python3 python3-venv python3-pip \
            libfftw3-dev \
            openmpi-bin libopenmpi-dev; \
    fi && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV PATH="/opt/MFC:$PATH"

COPY ../ /opt/MFC

ENV CC=${CC_COMPILER}
ENV CXX=${CXX_COMPILER}
ENV FC=${FC_COMPILER}
ENV PATH="${COMPILER_PATH}:$PATH"
ENV LD_LIBRARY_PATH="${COMPILER_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"

RUN echo "TARGET=$TARGET CC=$CC_COMPILER FC=$FC_COMPILER" && \
    cd /opt/MFC && \
    if [ "$TARGET" = "gpu" ]; then \
      ./mfc.sh build --gpu -j $(nproc); \
    else \
      ./mfc.sh build -j $(nproc); \
    fi

RUN cd /opt/MFC && \
    if [ "$TARGET" = "gpu" ]; then \
      ./mfc.sh test -a --dry-run --gpu -j $(nproc); \
    else \
      ./mfc.sh test -a --dry-run -j $(nproc); \
    fi

WORKDIR /opt/MFC
ENTRYPOINT ["tail", "-f", "/dev/null"]