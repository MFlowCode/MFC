#!/bin/bash

declare -a dependencies

dependencies[0]="FFTW3|http://www.fftw.org/fftw-3.3.10.tar.gz"
dependencies[1]="LAPACK|http://www.netlib.org/lapack/lapack-3.4.2.tgz"
dependencies[2]="HDF5|https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.21.tar.gz"
dependencies[3]="SILO|https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11.tgz"

echo Please ensure you have an appropriate build environment
echo Note: This file is meant to be called from within dependencies/
echo 
echo "|------------------------------------------------------"
echo "|--------------- Fetching Dependencies ----------------"
echo "|------------------------------------------------------"
echo "|"

# For each dependency
for dependency in "${dependencies[@]}"
do
    # Extract name and download link
    IFS="|" read -r -a arr <<< "${dependency}"

    name="${arr[0]}"
    link="${arr[1]}"

    # If we haven't downloaded it before (the directory named $name doesn't exist)
    if [ ! -d $(pwd)"/"$name ]
    then
        archive_filename=$name".tar.gz"
        
        echo "+ Fetching "$name" @ "$link

        # Download Archive
        wget -O $archive_filename -q $link

        mkdir -p $name

        cd $name
            echo "| Uncompressing "$archive_filename

            tar --strip-components 1 -xf "../"$archive_filename
        cd ..

        echo "| Removing "$archive_filename
        
        rm $archive_filename
    else
        echo "+ Skipping "$name" because it has already been downloaded"
    fi

    echo "|"
done

echo -------------------------------------------------------
echo ---------------- Building Dependencies ----------------
echo -------------------------------------------------------
echo 

# FFTW3
cd FFTW3
    ./configure --prefix=$(pwd) --enable-threads --enable-mpi

    make "$@"
    make install
cd ..

# LAPACK
cd LAPACK
    cat ../make_lapack.inc | sed -e s+XXMAKEDIR+$(pwd/..)+g > make.inc

    make "$@" blaslib
    make "$@" lapacklib
    cp librefblas.a libblas.a
cd ..

# HDF5
cd HDF5
    ./configure --prefix=$(pwd) --enable-parallel --enable-fortran --enable-cxx

    make "$@"
    make install
cd ..

# SILO
cd SILO
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --cflags)"

    ./configure --prefix=$(pwd) --enable-pythonmodule --enable-optimization    \
                --disable-hzip --disable-fpzip FC=mpif90 F77=mpif77            \
                -with-hdf5=$(pwd)/../HDF5/hdf5/include,$(pwd)/../HDF5/hdf5/lib \
                --disable-silex

    make "$@"
    make install
cd ..

echo "|"
echo "|------------------------------------------------------"
echo "|--------------------- Completed ----------------------"
echo "|------------------------------------------------------"
echo 
echo If a build failed, please ensure you have installed
echo the necessary dependencies to build MFC\'s dependencies
echo 
