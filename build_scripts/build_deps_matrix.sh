#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"
OS=$1
MATFILE=$2

for ARCH in `cat ${MATFILE}`
do
    echo "Building ${ARCH}"
    bash "${MYDIR}/build_deps_${OS}.sh" ${ARCH} &> ${OS}_${ARCH}.build.log
done
