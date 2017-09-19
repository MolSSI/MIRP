#!/bin/bash

set -eu

COMP_MOD=$1
CMAKE_MOD=$2
DEPS_ARCH=$3
OMP=$4
BTYPE=$5

OUTFILE="${COMP_MOD}_${CMAKE_MOD}_deps-${DEPS_ARCH}_omp-${OMP}-${BTYPE}.log"
OUTFILE=${OUTFILE//\//_}  # Replace all slashes with underscore
OUTFILE="$(pwd)/${OUTFILE}"
rm -f ${OUTFILE}

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
  echo "Unknown compiler: ${COMP_NAME}" >> ${OUTFILE}
  exit 1
fi
 
CMAKE="$(which cmake)"

FULL_DEPS_DIR="${DEPS_DIR}/mirp_deps_${DEPS_ARCH}"
BUILD_DIR="$(mktemp -d -p /tmp)"

echo "    Build dir: ${BUILD_DIR}"     >> ${OUTFILE}
echo " Dependencies: ${FULL_DEPS_DIR}" >> ${OUTFILE}
echo "    Deps Arch: ${DEPS_ARCH}"     >> ${OUTFILE}
echo "           CC: ${CC}"            >> ${OUTFILE}
echo "          CXX: ${CC}"            >> ${OUTFILE}
echo "        CMAKE: ${CMAKE}"         >> ${OUTFILE}
echo "   BUILD TYPE: ${BTYPE}"         >> ${OUTFILE}
echo "       OPENMP: ${OMP}"           >> ${OUTFILE}

cd "${BUILD_DIR}"
${CMAKE} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} \
         -DCMAKE_BUILD_TYPE=${BTYPE} \
         -DMIRP_OPENMP=${OMP} \
         -DCMAKE_PREFIX_PATH=${FULL_DEPS_DIR} \
         ${MIRP_DIR} &>> ${OUTFILE}

make &>> ${OUTFILE}
ctest &>> ${OUTFILE}
