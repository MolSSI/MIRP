#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"
MIRP_DIR="$(dirname ${MYDIR})"

COMP_MOD=$1
CMAKE_MOD=$2
DEPS_ARCH=$3
OMP=$4
BTYPE=$5

module load ${COMP_MOD}
module load ${CMAKE_MOD}

# Determine CC and CXX from the module name/version
COMP_NAME="${COMP_MOD%/*}"   # Erase the slash and everything after
COMP_VER="${COMP_MOD#*/}"    # Erase everything up to the slash

if [[ "${COMP_NAME}" == "gcc" ]]
then
  CC="$(which gcc)"
  CXX="$(which g++)"
elif [[ "${COMP_NAME}" == "clang" ]]
then
  CC="$(which clang)"
  CXX="$(which clang++)"
else
  echo "Unknown compiler: ${COMP_NAME}"
  exit 1
fi
 
CMAKE="$(which cmake)"

FULL_DEPS_DIR="${DEPS_DIR}/mirp_deps_${DEPS_ARCH}"
BUILD_DIR="$(mktemp -d -p /tmp)"
CURDIR="$(pwd)"

echo "    Build dir: ${BUILD_DIR}"
echo " Dependencies: ${FULL_DEPS_DIR}"
echo "    Deps Arch: ${DEPS_ARCH}"
echo "           CC: ${CC}"
echo "          CXX: ${CXX}"
echo "        CMAKE: ${CMAKE}"
echo "   BUILD TYPE: ${BTYPE}"
echo "       OPENMP: ${OMP}"

cd "${BUILD_DIR}"
${CMAKE} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} \
         -DCMAKE_BUILD_TYPE=${BTYPE} \
         -DMIRP_OPENMP=${OMP} \
         -DCMAKE_PREFIX_PATH=${FULL_DEPS_DIR} \
         ${MIRP_DIR}

make
ctest

cd "${CURDIR}"
rm -Rf "${BUILD_DIR}"
