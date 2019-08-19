#!/bin/bash

set -e
set -u
set -x

sourcedir="$(pwd)"
builddir=build-contraction

cd ..

###############################################################################
#                              Install Packages                               #
###############################################################################

ubuntu_packages=(
    cmake
    libhdf5-dev 
    hdf5-tools
    libeigen3-dev
    libboost-filesystem-dev libboost-system-dev libboost-program-options-dev libboost-graph-dev
    libgtest-dev
    openmpi-bin openmpi-common libopenmpi-dev
    libopenblas-base libopenblas-dev
)
sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) main universe restricted multiverse"
#sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y "${ubuntu_packages[@]}"




###############################################################################
#                               Install C-LIME                                #
###############################################################################

git clone https://github.com/usqcd-software/c-lime.git
pushd c-lime
./autogen.sh
./configure
make -j $(nproc)
sudo make install
limedir=$(pwd)
popd

###############################################################################
#                               Install tmLQCD                                #
###############################################################################
git clone https://github.com/etmc/tmlqcd.git tmlqcd
pushd tmlqcd
autoconf
CC=mpicc \
CFLAGS="-O3 -std=c99" \
./configure --disable-omp --enable-mpi --with-mpidimension=4 \
  --disable-sse2 --disable-sse3 \
  --enable-halfspinor --enable-gaugecopy \
  --enable-alignment=32 \
  --with-limedir=${limedir} \
  --with-lapack="-lblas -llapack"
make -j $(nproc)
tmlqcd_dir=$(pwd)
popd

#################################################################################
###                              Build Google Test                              #
#################################################################################
##
### https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/
##
##mkdir "gtest"
##pushd "gtest"
##cmake /usr/src/gtest
##make -j $(nproc)
##sudo cp *.a /usr/lib/
##popd
##
#################################################################################
###                          Build sLapH-contractions                           #
#################################################################################
##
##rm -rf "$builddir"
##mkdir -p "$builddir"
##pushd "$builddir"
##
##CXX=$(which g++)
##
### Compile gtest
### Modified from https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/
##pushd /usr/src/gtest
##sudo cmake CMakeLists.txt -DCMAKE_CXX_COMPILER="$CXX"
##sudo make -j $(nproc)
##sudo cp *.a /usr/lib
##popd
##
##cmake "$sourcedir" -DCMAKE_MODULE_PATH=../cmake-module -DCMAKE_CXX_COMPILER="$CXX" \
##  -DLIME_INCLUDE_DIRS="$limedir/include" -DLIME_LIBRARIES="-L$limedir/libs -llime"
##make -j $(nproc) || make VERBOSE=1
##
##ctest --output-on-failure