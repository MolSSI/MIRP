#!/bin/bash

rm -f matrix.txt

ALL_DEPS_ARCH="haswell"
ALL_COMP_MOD="gcc/7.2.0 gcc/7.1.0 gcc/6.4.0 gcc/5.4.0 gcc/4.9.4 clang/4.0.1 clang/5.0.0"
ALL_CMAKE_MOD="cmake/3.9.2 cmake/3.8.2 cmake/3.7.2 cmake/3.6.3 cmake/3.5.2 cmake/3.4.3 cmake/3.3.2 cmake/3.2.3"
ALL_OMP="True False"
ALL_BTYPES="Debug Release"

for COMP_MOD in $ALL_COMP_MOD; do
  for CMAKE_MOD in $ALL_CMAKE_MOD; do
    for DEPS_ARCH in $ALL_DEPS_ARCH; do
      for OMP in $ALL_OMP; do
        for BTYPE in $ALL_BTYPES; do
          echo "${COMP_MOD} ${CMAKE_MOD} ${DEPS_ARCH} ${OMP} ${BTYPE}" >> matrix.txt
        done
      done
    done
  done
done
