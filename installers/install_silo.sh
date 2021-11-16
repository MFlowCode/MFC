#!/bin/bash

cd ../
homedir=$(pwd)
cd lib

library_name=silo
library_version=4.11

wget -O "$library_name.tgz" "https://wci.llnl.gov/sites/wci/files/2021-09/"$library_name"-"$library_version".tgz"

mkdir -p $library_name

cd $library_name
    tar --strip-components 1 -xf "../"$library_name".tgz"
cd ..

rm $library_name".tgz"

cd $library_name
    echo $(pwd)
    export PYTHON_CPPFLAGS=$(python3-config --cflags)
    ./configure --prefix=$(pwd) --enable-pythonmodule --enable-optimization --disable-hzip --disable-fpzip FC=mpif90 F77=mpif77 -with-hdf5=$(pwd)/../hdf5/hdf5/include,$(pwd)/../hdf5/hdf5/lib --disable-silex
    make "$@"
    make install
cd ..
