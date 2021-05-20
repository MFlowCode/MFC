#!/bin/bash

cd ../lib
lib_dir=$(pwd)

tar -xvf hdf5-1.8.20.tar
tar -zxvf silo-4.10.2.tar.gz

mv hdf5-1.8.20 hdf5
mv silo-4.10.2 silo

cp config* hdf5/bin/
cp config* silo/config/

### INSTALLING HDF5
cd hdf5
./configure --prefix=$lib_dir/hdf5
make
make install

### INSTALLING SILO
cd ../silo
./configure --prefix=$lib_dir/silo --enable-pythonmodule --enable-optimization --disable-hzip --disable-fpzip --enable-portable-binary FC=mpif90 F77=mpif77 --with-hdf5=$lib_dir/hdf5/include,$lib_dir/hdf5/lib --disable-silex
make 
make install
