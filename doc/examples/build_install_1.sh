#!/bin/bash

# Abort if an error is encountered or on undefined variables
set -eu


# Download MIRP sources
git clone https://github.com/MolSSI/MIRP.git

# Build dependencies in a separate directory
mkdir mirp_deps
cd mirp_deps
bash ../MIRP/build_scripts/build_deps.sh native
cd ../

# mirp_deps directory should now contain
# mirp_deps_native subdirectory

# Now build MIRP
mkdir mirp_build
cd mirp_build
cmake -DCMAKE_PREFIX_PATH="$(pwd)/../mirp_deps/mirp_deps_native" \
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
      -DCMAKE_INSTALL_PREFIX="$(pwd)/../mirp_install" \
      -DCMAKE_BUILD_TYPE=Release \
      ../MIRP

make install
