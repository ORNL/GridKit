#!/bin/bash
# Enzyme v0.0.206 — AD plugin for LLVM 16
# Must be built against the same LLVM 16 built in 1_build_llvm.sh (RTTI=ON).

echo "=== Building Enzyme v0.0.206 ==="

DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
ENZYME_PREFIX=$DEPS_BASE/enzyme/0.0.206
SRC_DIR=$DEPS_BASE/src
mkdir -p "$SRC_DIR" "$ENZYME_PREFIX"
cd $SRC_DIR

module load gcc-stdalone/13.1.0
module load cmake/3.29.2

GCC13_ROOT=$(dirname $(dirname $(which g++)))

export PATH=$LLVM_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_PREFIX/lib:$LD_LIBRARY_PATH

# Download Enzyme v0.0.206
wget https://github.com/EnzymeAD/Enzyme/archive/refs/tags/v0.0.206.tar.gz
tar -xzf v0.0.206.tar.gz
cd Enzyme-0.0.206/enzyme

rm -rf build
mkdir build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$ENZYME_PREFIX \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=$LLVM_PREFIX/bin/clang \
  -DCMAKE_CXX_COMPILER=$LLVM_PREFIX/bin/clang++ \
  -DCMAKE_CXX_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_C_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DLLVM_DIR=$LLVM_PREFIX/lib/cmake/llvm \
  -DENZYME_CLANG=ON

make -j$(nproc)
make install

echo ""
echo "Enzyme installed to: $ENZYME_PREFIX/lib/LLVMEnzyme-16.so"
ls -la $ENZYME_PREFIX/lib/LLVMEnzyme-16.so
