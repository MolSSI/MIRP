#!/bin/bash

set -eu

export MYDIR="$(cd "$(dirname "$0")" && pwd)"
export MIRP_DIR="$(dirname ${MYDIR})"

run_build() {
    MIRP_DIR=${MIRP_DIR} bash "${MYDIR}/build_mirp.sh" $@
    if [[ $? == 0 ]]
    then
      echo "passed : $@"
    else
      echo "FAILED : $@"
    fi
}

export -f run_build


if [[ "${DEPS_DIR}" == "" ]]
then
  echo "DEPS_DIR variable is not set"
  exit 1
fi

cat $1 | parallel --colsep ' ' run_build

