#!/bin/bash

cd ..
homedir=$(pwd)
cd lib

library_name=lapack
library_version=3.4.2

wget -O $library_name".tar.gz" "http://www.netlib.org/lapack/"$library_name"-"$library_version".tgz"

mkdir -p $library_name

cd $library_name
    tar --strip-components 1 -xf "../"$library_name".tar.gz"
cd ..

rm $library_name".tar.gz"

cd $library_name
    cat ../make_lapack.inc | sed -e s+XXMAKEDIR+$homedir+g > make.inc
    make "$@" blaslib
    make "$@" lapacklib
    cp librefblas.a libblas.a
cd ..
