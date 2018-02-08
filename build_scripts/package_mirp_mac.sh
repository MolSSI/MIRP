#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"
MIRP_DIR="$(dirname ${MYDIR})"

MIRP_VER=$(cat "${MIRP_DIR}/VERSION")

ARCH=$1

CC=`which clang`
CXX=`which clang++`

CURDIR="$(pwd)"

FULL_DEPS_DIR="${DEPS_DIR}/mirp_deps_v${MIRP_VER}_${ARCH}"
PREFIX="${CURDIR}/mirp_v${MIRP_VER}_${ARCH}"

# Copy the dependencies to the install prefix
rm -Rf "${PREFIX}"
cp -R "${FULL_DEPS_DIR}" "${PREFIX}"

BUILD_DIR="$(mktemp -d /tmp/mirp_pkg.${ARCH}.XXXXXX)"
cd "${BUILD_DIR}"
cmake -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} \
      -DCMAKE_C_FLAGS="-march=${ARCH}" \
      -DCMAKE_CXX_FLAGS="-march=${ARCH}" \
      -DCMAKE_BUILD_TYPE=Release \
      -DMIRP_OPENMP=False \
      -DMIRP_STATIC=False \
      -DCMAKE_PREFIX_PATH=${PREFIX} \
      -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
      ${MIRP_DIR}

make
ctest
make install

cd "${CURDIR}"
rm -Rf "${BUILD_DIR}"

# Fix the rpaths (if we have install_name_tool)
if [[ $(command -v install_name_tool 2>&1) ]]
then
    for I in ${PREFIX}/bin/*
    do
        RP1=`otool -l $I | grep -A3 LC_RPATH | tail -n 1 | awk '{print $2}'`
        install_name_tool -add_rpath '@executable_path/../lib' "$I"
        RP2=`otool -l $I | grep -A3 LC_RPATH | tail -n 1 | awk '{print $2}'`
        echo "${I}: RPATH changed from \"${RP1}\" to \"${RP2}\""

        # Change all the dependency entries to use rpath
        otool -L ${I} | grep ${PREFIX} | awk '{print $1}' | while read F
        do
            FNAME=`basename $F`
            install_name_tool -change "$F" "@rpath/$FNAME" "$I"
        done
    done
else
    echo
    echo "!!! install_name_tool not installed. Skipping fixing RPATHS !!!"
    echo
fi

mv ${PREFIX}/README ${PREFIX}/README.dependencies

# Copy documentation
cd "${PREFIX}"
cp -R ${MIRP_DIR}/doc/html html_doc
tar -cvjf html_doc.tar.bz2 html_doc
rm -R html_doc
cd "${CURDIR}"

# Create the readme file
COMPILER_VER=$(${CC} --version | head -n 1)
BUILD_DATE=$(date +%Y-%m-%d)
cp "${MIRP_DIR}/LICENSE"                      "${PREFIX}"
cp "${MYDIR}/mirp_README.in"                  "${PREFIX}/README"
sed -i "" "s/MIRP_VER/${MIRP_VER}/g"          "${PREFIX}/README"
sed -i "" "s/COMPILER_VER/${COMPILER_VER}/g"  "${PREFIX}/README"
sed -i "" "s/ARCH/${ARCH}/g"                  "${PREFIX}/README"
sed -i "" "s/BUILD_DATE/${BUILD_DATE}/g"      "${PREFIX}/README"
