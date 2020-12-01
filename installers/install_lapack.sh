#!/bin/bash

cd ../
homedir=$(pwd)
cd lib

library=lapack-3.4.2
tar -zxf $library'.tgz'
cd $library
    # cp ../make_lapack.inc make.inc
    cat ../make_lapack.inc | sed -e s+XXMAKEDIR+$homedir+g > make.inc
    make blaslib 
    make lapacklib 
    make install
cd ..
