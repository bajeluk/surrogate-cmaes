#!/bin/bash
#
# bobyqa_install.sh -- Install the BOBYQA
#
# Requirements (Debian/Ubuntu packages in brackets)
# - GCC (gcc)
# - Fortran95/GCC (gfortran)
# - CMake (cmake)
# - BLAS library (libblas3)

DLIB_VER=19.4
DLIB_INSTALL="http://dlib.net/files/dlib-$DLIB_VER.tar.bz2"
DLIB_DEST_DIR="exp/vendor/dlib"


last_addr=`pwd`
# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "$CWD/.."

if [ ! -d "$CWD/../../../$DLIB_DEST_DIR" ]; then
  # Download the installation
  if [ ! -f "dlib-$DLIB_VER.tar.bz2" ]; then
    wget "$DLIB_INSTALL"
  fi

  # Extract Dlib
  tar xjf dlib-${DLIB_VER}.tar.bz2 --wildcards "dlib-${DLIB_VER}/dlib/"
  mv dlib-${DLIB_VER}/dlib .
  rm -rf dlib-${DLIB_VER}
else
  echo "Directory ${DLIB_DEST_DIR} already exists."
  exit 0
fi

cp bobyqa/dlib_bobyqa.cpp dlib/matlab/
patch -p1 < bobyqa/optimization_bobyqa.h.patch
if ! grep -q dlib_bobyqa dlib/matlab/CMakeLists.txt; then
  echo "add_mex_function(dlib_bobyqa dlib)" >> dlib/matlab/CMakeLists.txt
fi

# Make Dlib
mkdir dlib/matlab/build
cd dlib/matlab/build
cmake ..
cmake --build . --config release

# Install BOBYQA
# TODO

# Finalize
cd "$last_dir"
