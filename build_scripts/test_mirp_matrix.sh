#!/bin/bash

set -eu

export MYDIR="$(cd "$(dirname "$0")" && pwd)"

run_build() {
    OUTFILE="${1}_${2}_deps-${3}_omp-${4}-${5}.log"
    OUTFILE=${OUTFILE//\//_}  # Replace all slashes with underscore

    bash "${MYDIR}/test_mirp.sh" $@ &> ${OUTFILE}

    if [[ $? == 0 ]]
    then
      echo "passed : $@"
    else
      echo "FAILED : $@"
    fi
}

export -f run_build

cat $1 | parallel --colsep ' ' run_build

