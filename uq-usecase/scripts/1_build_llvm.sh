#!/bin/bash
# LLVM 16.0.6 — built from source with LLVM_ENABLE_RTTI=ON
# RTTI is required by Enzyme; the Kestrel llvm/16 module was built without it.

echo "=== Building LLVM 16.0.6 ==="

DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
SRC_DIR=$DEPS_BASE/src
mkdir -p "$SRC_DIR" "$LLVM_PREFIX"
cd $SRC_DIR

module load gcc-stdalone/13.1.0
module load cmake/3.29.2

GCC13_ROOT=$(dirname $(dirname $(which g++)))

# Download LLVM 16.0.6
wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-16.0.6.tar.gz
tar -xzf llvmorg-16.0.6.tar.gz
cd llvm-project-llvmorg-16.0.6

rm -rf build
mkdir build && cd build
cmake ../llvm \
  -DCMAKE_INSTALL_PREFIX=$LLVM_PREFIX \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DLLVM_ENABLE_PROJECTS="clang;openmp" \
  -DLLVM_ENABLE_RTTI=ON \
  -DLLVM_ENABLE_EH=ON \
  -DLLVM_TARGETS_TO_BUILD="X86" \
  -DLLVM_BUILD_EXAMPLES=OFF \
  -DLLVM_BUILD_TESTS=OFF \
  -DLLVM_INCLUDE_TESTS=OFF \
  -DLLVM_INCLUDE_BENCHMARKS=OFF

make -j$(nproc)
make install

echo ""
echo "LLVM 16.0.6 installed to: $LLVM_PREFIX"
echo "Verify: $LLVM_PREFIX/bin/clang --version"
