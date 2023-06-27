FROM       ubuntu:jammy
MAINTAINER Henry Le Berre <hberre3@gatech.edu>

ENV PS1='[\u@mfc] \e[1;32m\w\e[0m \$ '
ENV CMAKE_GENERATOR=Ninja

RUN apt update && apt upgrade -y

RUN apt install -y git tar wget make cmake ninja-build gcc g++  \
                   python3 python3-dev python3-pip python3-venv \
                   'openmpi-*' libopenmpi-dev sudo vim nano     \
                   hdf5-tools

RUN groupadd --gid 1000 me && \
     useradd --uid 1000 --gid 1000 -m me

RUN echo me ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/me && \
    chmod 0440 /etc/sudoers.d/me

USER me

RUN pip3 install fypp rich pyyaml

CMD ["/bin/bash"]
