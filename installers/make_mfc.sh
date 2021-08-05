#!/bin/bash

cd ~

# Clean first
rm -rf MFC-Caltech packages/silo-4.10.2 packages/hdf5

user=$(whoami)

# Make sure we have appropriate modules and prepend path for Python
echo 'module purge' >> ~/.bash_profile
echo 'module load openmpi-1.8/gcc' >> ~/.bash_profile
echo 'export PATH=/usr/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile

mkdir ~/packages
cd ~/packages

# Download HDF5 binaries
wget https://support.hdfgroup.org/ftp/HDF5/current18/bin/hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz 
tar -zxf hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz

# Download and make Silo
wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
tar -zxf silo-4.10.2.tar.gz
cd silo-4.10.2 
./configure --prefix=/home/$user/packages/silo-4.10.2 --enable-pythonmodule --enable-optimization --disable-hzip --disable-fpzip FC=mpif90 F77=mpif77 -with-hdf5=/home/$user/packages/hdf5/include,/home/$user/packages/hdf5/lib --disable-silex
make
make install

# Paths
cd ~
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/packages/silo-4.10.2/lib:$HOME/packages/hdf5/lib' >> ~/.bash_profile
source ~/.bash_profile

# Make MFC
git clone https://github.com/ComputationalFlowPhysics/MFC-Caltech.git
cd MFC-Caltech/installers
./install_fftw.sh
./install_lapack.sh
cd ..

echo "silo_lib_dir = $HOME/packages/silo-4.10.2/lib" >> Makefile.user
echo "silo_include_dir = $HOME/packages/silo-4.10.2/include" >> Makefile.user
make
make test
