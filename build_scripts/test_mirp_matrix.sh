#!/bin/bash

set -eu

export MYDIR="$(cd "$(dirname "$0")" && pwd)"
export OS=$1
MATFILE=$2

run_build() {
    OUTFILE="${OS}_${1}_${2}_deps-${3}_omp-${4}-${5}.log"
    OUTFILE=${OUTFILE//\//_}  # Replace all slashes with underscore

    bash "${MYDIR}/test_mirp_${OS}.sh" $@ &> ${OUTFILE}

    if [[ $? == 0 ]]
    then
      echo "passed : $@"
    else
      echo "FAILED : $@"
    fi
}

export -f run_build

cat ${MATFILE} | parallel -j1 --colsep ' ' run_build
