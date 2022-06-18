FROM archlinux:latest

# Install packages/dependencies
RUN pacman-key --init
RUN pacman -Syyu --noconfirm
RUN pacman -S --noconfirm          \
           git ninja cmake wget    \
           python3 openmpi vim     \
           coreutils gcc-fortran   \
           openssh base-devel tree

# Generic (non-root) user because MPI dislikes running as root.
# This user is allowed to run "sudo" commands without a password.
RUN     useradd --create-home --shell /bin/bash me
RUN     passwd -d me
RUN     echo "me ALL=(ALL:ALL) ALL" >> /etc/sudoers
USER    me
WORKDIR /home/me

# Env
ENV PS1='[\u@mfc] \e[1;32m\w\e[0m \$ '
ENV CMAKE_GENERATOR=Ninja

# Bashrc
RUN echo "export PS1='[\u@mfc] \e[1;32m\w\e[0m \$ '" >> /home/me/.bashrc

# Run bash as the default command
CMD ["/bin/bash"]
