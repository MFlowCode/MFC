FROM archlinux:latest

# Install packages/dependencies
RUN pacman-key --init
RUN pacman -Syyu --noconfirm
RUN pacman -S --noconfirm         \
           base-devel coreutils   \
           git ninja gcc-fortran  \
           cmake openmpi python3  \
           python-pip openssh     \
           python-virtualenv vim  \
           wget tree

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
RUN echo $' \
export PS1="[\u@mfc] \e[1;32m\w\e[0m \$ " \n\
echo " " \n\
echo "        Welcome to the MFC Docker Environment!" \n\
echo " " \n\
if [ -d "/home/me/MFC" ]; then \n\
    echo "  Please wait a few seconds while MFC configures directory permissions." \n\
    echo "    > sudo chown -R me:me /home/me/MFC" \n\
    sudo chown -R me:me /home/me/MFC \n\
    echo "    > Complete." \n\
    echo "  Follow the README\'s instructions to get started with MFC." \n\
else \n\
    echo "  Error: The MFC mount isn\'t available." \n\
fi \n\
echo "" \n\
cd /home/me/MFC \n\
' >> /home/me/.bashrc

# Run bash as the default command
CMD ["/bin/bash"]
