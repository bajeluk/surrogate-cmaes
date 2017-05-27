#!/bin/bash
#
# bobyqa_install.sh -- Compile and install the BOBYQA algorithm
#
# BOBYQA is a quadratic-approximation optimizer by J. Powell's (2009).
# This is a installation script which downloads the Dlib C++ machine-learning
# library (which contains BOBYQA) and compiles it (after a minor patch). After
# that, Matlab MEX file dlib_bobyqa.cpp is compiled via GCC in order to
# be used in Matlab.
#
# Requirements (Debian/Ubuntu packages in brackets)
# - GCC (gcc)
# - Fortran95/GCC (gfortran)
# - CMake (cmake)
# - BLAS library (libblas3)
# - dlib_bobyqa.cpp Matlab/C++ MEX file
#
# (CC) Lukas Bajer, 2017

DLIB_VER=19.4
DLIB_INSTALL="http://dlib.net/files/dlib-$DLIB_VER.tar.bz2"
DLIB_DEST_DIR="exp/vendor/dlib"

last_addr=`pwd`
# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "$CWD/.."

if [ ! -d "$CWD/../../../$DLIB_DEST_DIR" ]; then
  # Download the Dlib's installation package
  if [ ! -f "dlib-$DLIB_VER.tar.bz2" ]; then
    echo "Trying to download Dlib library..."
    wget "$DLIB_INSTALL"
    if [ $? -ne 0 ]; then
      echo "Error: Dlib cannot be downloaded nor it exists. Exitting."
      exit 1
    fi
  fi

  # Extract only the directory 'dlib/' form Dlib's package
  tar xjf dlib-${DLIB_VER}.tar.bz2 --wildcards "dlib-${DLIB_VER}/dlib/"
  mv dlib-${DLIB_VER}/dlib .
  rm -rf dlib-${DLIB_VER}
else
  echo "Directory ${DLIB_DEST_DIR} already exists."
fi

cp bobyqa/dlib_bobyqa.cpp dlib/matlab/
patch -N -r- -p1 < bobyqa/optimization_bobyqa.h.patch
if ! grep -q dlib_bobyqa dlib/matlab/CMakeLists.txt; then
  sed -i 's/^[^#].*example_mex.*/# \0/' dlib/matlab/CMakeLists.txt
  echo "add_mex_function(dlib_bobyqa dlib)" >> dlib/matlab/CMakeLists.txt
fi

# Make Dlib
mkdir -p dlib/matlab/build
cd dlib/matlab/build
cmake .. && cmake --build . --config release

# Deploy BOBYQA to the right place
if [ $? -eq 0 ]; then
  cd "$CWD/.."
  cp dlib/matlab/build/dlib_bobyqa.mex* bobyqa/
else
  echo ""
  echo "Dlib or BOBYQA is not sucessfully compiled. Exitting with error."
  exit 1
fi
