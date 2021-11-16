#!/bin/bash

cd ../lib

library_name=fftw
library_version=3.3.10

wget -O $library_name".tar.gz" "http://www.fftw.org/"$library_name"-"$library_version".tar.gz"

mkdir -p $library_name

cd $library_name
    tar --strip-components 1 -xf "../"$library_name".tar.gz"
cd ..

rm $library_name".tar.gz"

cd $library_name
    ./configure --prefix=$(pwd) --enable-threads --enable-mpi
    make "$@"
    make install
cd ..
