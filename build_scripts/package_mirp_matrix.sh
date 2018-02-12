#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"
OS=$1
MATFILE=$2

for ARCH in `cat ${MATFILE}`
do
    echo "Building ${ARCH}"
    bash "${MYDIR}/package_mirp_${OS}.sh" ${ARCH} &> mirp_${ARCH}_package.log
done
