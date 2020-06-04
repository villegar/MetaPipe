#!/usr/bin/env sh

# check if OpenMPI is cached from previous build
if [ -f "openmpi/bin/mpirun"]; then
  echo "Using cached OpenMPI"
else
  export VERSION=4.0.3
  echo "Downloading OpenMPI $VERSION source"
  wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-$VERSION.tar.gz
  tar xfz openmpi-$VERSION.tar.gz
  rm openmpi-$VERSION.tar.gz
  echo "Configuring and building openmpi"
  cd openmpi-$VERSION
  ./configure --prefix=$(pwd)
  make -j 4 all
  make install
  export PATH=$PATH:$(pwd)/bin
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/lib
  cd ..
fi