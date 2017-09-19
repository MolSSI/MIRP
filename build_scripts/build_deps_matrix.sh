#!/bin/bash

set -eu

MYDIR="$(cd "$(dirname "$0")" && pwd)"

for M in `cat $1`
do
    echo "Building ${M}"
    bash "${MYDIR}/build_deps.sh" ${M} &> $M.build.log
done
