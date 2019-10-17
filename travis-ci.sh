#!/bin/bash

set -e
set -u
set -x

sourcedir="$(pwd)"

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
    libopenblas-base libopenblas-dev
    libmpich-dev mpich
    numdiff
)
    #openmpi-bin openmpi-common libopenmpi-dev
sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) main universe restricted multiverse"
sudo apt-get update
sudo apt-get install -y "${ubuntu_packages[@]}"

# Installation path for the various libraries to be installed below
mkdir -p local
install_prefix="$(realpath local)"

###############################################################################
#                               Install C-LIME                                #
###############################################################################
git clone https://github.com/usqcd-software/c-lime.git
pushd c-lime
./autogen.sh
CFLAGS="-fPIC -O0 -g" \
  ./configure --prefix="$install_prefix"
make -j $(nproc)
make install
popd


###############################################################################
#                               Install LEMON                                 #
###############################################################################
git clone https://github.com/etmc/lemon lemon
pushd lemon
lemon_srcdir="$(pwd)"
autoreconf -vfi
popd

mkdir -p lemon_builddir
pushd lemon_builddir
lemon_builddir="$(pwd)"
CC=$(which mpicc) \
  CFLAGS="-fPIC -O0 -g" \
  "$lemon_srcdir/configure" --prefix="$install_prefix"
make -j $(nproc)
make install
popd

###############################################################################
#                               Install tmLQCD                                #
###############################################################################
git clone https://github.com/etmc/tmlqcd.git tmlqcd -b quda_work
pushd tmlqcd
tmlqcd_srcdir="$(pwd)"
autoconf
popd

rm -f tmlqcd_builddir
mkdir -p tmlqcd_builddir
pushd tmlqcd_builddir
CC=$(which mpicc) \
CFLAGS="-O0 -g -std=c99 -fPIC" \
"$tmlqcd_srcdir"/configure --disable-omp --enable-mpi --with-mpidimension=4 \
  --disable-sse2 --disable-sse3 \
  --enable-halfspinor --enable-gaugecopy \
  --enable-alignment=32 \
  --with-limedir="$install_prefix" \
  --with-lemondir="$install_prefix" \
  --with-lapack="-lblas -llapack"
make -j $(nproc)
tmlqcd_builddir="$(pwd)"
popd

###############################################################################
#                         Fetch and install HighFive lib                      #
###############################################################################
# note we use v2.0 as 'known-good'
highfive_builddir=highfive_builddir
git clone https://github.com/BlueBrain/HighFive.git -b v2.0
mkdir -p "$highfive_builddir"
pushd "$highfive_builddir"
cmake \
  -DCMAKE_CXX_FLAGS="-O0 -g" \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DHIGHFIVE_UNIT_TESTS=FALSE \
  -DHIGHFIVE_EXAMPLES=FALSE \
  ../HighFive 
make -j $(nproc)
make install
popd

###############################################################################
#                         Fetch and install yaml-cpp                          #
###############################################################################
yamlcpp_builddir=yamlcpp_builddir
git clone https://github.com/jbeder/yaml-cpp
pushd yaml-cpp
git checkout e0e01d53c27ffee6c86153fa41e7f5e57d3e5c90
popd
mkdir -p "$yamlcpp_builddir"
pushd "$yamlcpp_builddir"
cmake \
  -DCMAKE_CXX_FLAGS="-fPIC -O0 -g" \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DYAML_CPP_BUILD_TESTS=OFF \
  ../yaml-cpp
make -j $(nproc)
make install
popd

###############################################################################
#                         Build correlators executable                        #
###############################################################################

builddir="$(pwd)/cvc_depgraph_builddir"
rm -rf "$builddir"
mkdir -p "$builddir"
pushd "$builddir"
CXX=$(which mpicxx)

cmake \
  -DCMAKE_CXX_FLAGS="-O0 -g -fPIC" \
  -DCMAKE_PREFIX_PATH="$install_prefix" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DLIME_HOME="$install_prefix" \
  -DLEMON_HOME="$install_prefix" \
  -DTMLQCD_SRC="$tmlqcd_srcdir" \
  -DTMLQCD_BUILD="$tmlqcd_builddir" \
  -DPARALLEL_LEVEL=TXYZ \
  "$sourcedir"

make -j $(nproc) correlators || VERBOSE=1 make correlators

###############################################################################
#                         run integration test(s)                             #
###############################################################################
ctest --output-on-failure

popd
