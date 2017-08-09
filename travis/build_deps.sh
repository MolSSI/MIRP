#!/bin/bash

set -eu

GMP_VER=6.1.2
MPFR_VER=3.1.5
FLINT_VER=2.5.2
ARB_VER=2.11.1

GMP_DIR="gmp-${GMP_VER}"
MPFR_DIR="mpfr-${MPFR_VER}"
FLINT_DIR="flint-${FLINT_VER}"
ARB_DIR="arb-${ARB_VER}"

GMP_FILE="gmp-${GMP_VER}.tar.bz2"
MPFR_FILE="mpfr-${MPFR_VER}.tar.bz2"
FLINT_FILE="flint-${FLINT_VER}.tar.gz"
ARB_FILE="arb-${ARB_VER}.tar.gz"

GMP_URL="https://ftp.gnu.org/gnu/gmp/${GMP_FILE}"
MPFR_URL="https://ftp.gnu.org/gnu/mpfr/${MPFR_FILE}"
FLINT_URL="http://flintlib.org/${FLINT_FILE}"
ARB_URL="https://github.com/fredrik-johansson/arb/archive/${ARB_VER}.tar.gz"

PREFIX="$(pwd)/mirp_deps"

mkdir -p deps_build
cd deps_build

# Download the file
wget -O "${GMP_FILE}" "${GMP_URL}"
wget -O "${MPFR_FILE}" "${MPFR_URL}"
wget -O "${FLINT_FILE}" "${FLINT_URL}"
wget -O "${ARB_FILE}" "${ARB_URL}"

rm -Rf "${GMP_DIR}"
rm -Rf "${MPFR_DIR}"
rm -Rf "${FLINT_DIR}"
rm -Rf "${ARB_DIR}"
tar -xf "${GMP_FILE}"
tar -xf "${MPFR_FILE}"
tar -xf "${FLINT_FILE}"
tar -xf "${ARB_FILE}"

mkdir -p "${GMP_DIR}/build"
mkdir -p "${MPFR_DIR}/build"

###################
# Common flags
###################
CFLAGS="-fomit-frame-pointer -O2 -march=native -m64"
CXXFLAGS="$CFLAGS"

###################
# Build GMP
###################
cd "${GMP_DIR}/build"
../configure CFLAGS="${CFLAGS}" CXXFLAGS="${CXXFLAGS}" \
             --disable-static \
             --prefix=${PREFIX}
make
#make check
make install
cd ../../

###################
# Build MPFR
###################
cd "${MPFR_DIR}/build"
../configure --with-gmp=${PREFIX} \
             CFLAGS="-fomit-frame-pointer -O2 -march=nehalem -m64" \
             CXXFLAGS="-fomit-frame-pointer -O2 -march=nehalem -m64" \
             --disable-static \
             --prefix=${PREFIX}
make
#make check
make install
cd ../../

###################
# Build FLINT
###################
cd "${FLINT_DIR}"
./configure --with-gmp=${PREFIX} --with-mpfr=${PREFIX} \
            CFLAGS="${CFLAGS}" \
            --disable-static \
            --prefix=${PREFIX}
make
#make check
make install
cd ../

###################
# Build ARB
###################
cd "${ARB_DIR}"
./configure --with-gmp=${PREFIX} --with-mpfr=${PREFIX} --with-flint=${PREFIX} \
            CFLAGS="${CFLAGS}" \
             --disable-static \
            --prefix=${PREFIX}
make
#make check
make install
cd ../

# Cleanup
rm -Rf "${GMP_DIR}"   "${GMP_FILE}"
rm -Rf "${MPFR_DIR}"  "${MPFR_FILE}"
rm -Rf "${FLINT_DIR}" "${FLINT_FILE}"
rm -Rf "${ARB_DIR}"   "${ARB_FILE}"

cd ../
