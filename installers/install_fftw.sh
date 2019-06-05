#!/bin/bash

cd ../lib

tar -zxf fftw-3.3.8.tar.gz
cd fftw-3.3.8
    lib_dir=$(pwd)
    ./configure --prefix=$lib_dir
    make 
    make install
cd ..
