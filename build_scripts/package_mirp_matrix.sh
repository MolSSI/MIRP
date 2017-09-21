#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"

for M in `cat $1`
do
    echo "Building ${M}"
    bash "${MYDIR}/package_mirp.sh" ${M} &> mirp_${M}_package.log
done
