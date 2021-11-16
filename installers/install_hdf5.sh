#!/bin/bash

cd ../lib

library_name=hdf5
library_version=1.8.21

wget -O $library_name".tar.gz" "https://support.hdfgroup.org/ftp/HDF5/current18/src/"$library_name"-"$library_version".tar.gz"

mkdir -p $library_name

cd $library_name
    tar --strip-components 1 -xf "../"$library_name".tar.gz"
cd ..

rm $library_name".tar.gz"

cd $library_name
    ./configure --prefix=$(pwd) --enable-parallel --enable-fortran --enable-cxx
    make "$@"
    make install
cd ..