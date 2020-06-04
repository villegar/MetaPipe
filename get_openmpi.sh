#!/usr/bin/env sh

export VERSION=4.0.3
# check if OpenMPI is cached from previous build
if [ -f "openmpi-$VERSION/bin/mpirun"]; then
  echo "Using cached for OpenMPI $VERSION"
  export MPI_DIR=`pwd`
  export PATH=$PATH:`pwd`/openmpi-$VERSION/bin
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/openmpi-$VERSION/lib
  export LIBRARY_PATH=$LIBRARY_PATH:`pwd`/openmpi-$VERSION/lib
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:`pwd`/openmpi-$VERSION/include
  export CPPFLAGS="-I`pwd`/openmpi-$VERSION/include"
  export LDFLAGS="-L`pwd`/openmpi-$VERSION/lib"
else
  echo "Downloading OpenMPI $VERSION source"
  rm -rf openmpi-$VERSION
  wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-$VERSION.tar.gz
  tar xfz openmpi-$VERSION.tar.gz
  rm openmpi-$VERSION.tar.gz
  echo "Configuring and building OpenMPI $VERSION"
  cd openmpi-$VERSION
  ./configure --prefix=`pwd` &> log.config
  make -j 4 all &> log.make
  make install
  export MPI_DIR=`pwd`
  export PATH=$PATH:`pwd`/bin
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/lib
  export LIBRARY_PATH=$LIBRARY_PATH:`pwd`/lib
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:`pwd`/include
  export CPPFLAGS="-I`pwd`/include"
  export LDFLAGS="-L`pwd`/lib"
  echo $CPPFLAGS
  ls $CPPFLAGS
  which mpicc
  echo "Done installing OpenMPI $VERSION"
  cd ..
fi