#!/bin/bash
# Build solve_pf against the already-built GridKit in ~/gridkit/build/.
# Run once from any directory:
#   bash ~/gridkit/uq-usecase/pf-solver/build.sh
#
# On success: ~/gridkit/uq-usecase/pf-solver/build/solve_pf

set -e

echo "=== Building solve_pf ==="

GRIDKIT_DIR=~/gridkit
DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6

PF_SRC_DIR=$GRIDKIT_DIR/uq-usecase/pf-solver
PF_BUILD_DIR=$PF_SRC_DIR/build
GRIDKIT_BUILD_DIR=$GRIDKIT_DIR/build

# Load the same modules used to build GridKit
module load gcc-stdalone/13.1.0
module load cmake/3.29.2

# Derive GCC 13 toolchain root (clang needs it for libstdc++ headers)
GCC13_ROOT=$(dirname $(dirname $(which g++)))
GCC13_LIBDIR=$(${GCC13_ROOT}/bin/g++ --print-file-name=libstdc++.so | xargs dirname)

echo "LLVM prefix: $LLVM_PREFIX"
echo "GCC13 root:  $GCC13_ROOT"
echo "GCC13 libs:  $GCC13_LIBDIR"
echo "GridKit build: $GRIDKIT_BUILD_DIR"

# Clean and reconfigure
rm -rf "$PF_BUILD_DIR"
mkdir -p "$PF_BUILD_DIR"
cd "$PF_BUILD_DIR"

cmake "$PF_SRC_DIR" \
  -DCMAKE_C_COMPILER="$LLVM_PREFIX/bin/clang" \
  -DCMAKE_CXX_COMPILER="$LLVM_PREFIX/bin/clang++" \
  -DCMAKE_CXX_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_C_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,${GCC13_LIBDIR}" \
  -DCMAKE_PREFIX_PATH="$GRIDKIT_BUILD_DIR" \
  -DGridKit_DIR="$GRIDKIT_BUILD_DIR/GridKit" \
  -DCMAKE_BUILD_TYPE=Release

make -j4

echo ""
echo "=== Done ==="
echo "Binary: $PF_BUILD_DIR/solve_pf"
ls -lh "$PF_BUILD_DIR/solve_pf"
