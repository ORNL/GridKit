#!/bin/bash
# SUNDIALS develop branch — required by GridKit for adjoint memory fix.
# Must use develop branch, NOT a release tarball.
# Built with KLU (from SuiteSparse) and OpenMP support.

echo "=== Building SUNDIALS (develop branch) ==="

DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
SUITESPARSE_PREFIX=$DEPS_BASE/suitesparse/7.8.2
SUNDIALS_PREFIX=$DEPS_BASE/sundials/develop
SRC_DIR=$DEPS_BASE/src
mkdir -p "$SRC_DIR" "$SUNDIALS_PREFIX"
cd $SRC_DIR

module load gcc-stdalone/13.1.0
module load cmake/3.29.2

GCC13_ROOT=$(dirname $(dirname $(which g++)))

export PATH=$LLVM_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_PREFIX/lib:$SUITESPARSE_PREFIX/lib64:$LD_LIBRARY_PATH
export CC=$LLVM_PREFIX/bin/clang
export CXX=$LLVM_PREFIX/bin/clang++

# Clone SUNDIALS develop branch (or update if already cloned)
if [ ! -d sundials-src ]; then
  git clone --branch develop --depth 1 https://github.com/LLNL/sundials.git sundials-src
else
  git -C sundials-src pull
fi
cd sundials-src

rm -rf build
mkdir build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$SUNDIALS_PREFIX \
  -DCMAKE_C_COMPILER=$LLVM_PREFIX/bin/clang \
  -DCMAKE_CXX_COMPILER=$LLVM_PREFIX/bin/clang++ \
  -DCMAKE_C_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_CXX_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_KLU=ON \
  -DKLU_INCLUDE_DIR=$SUITESPARSE_PREFIX/include/suitesparse \
  -DKLU_LIBRARY_DIR=$SUITESPARSE_PREFIX/lib64 \
  -DENABLE_OPENMP=ON \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_STATIC_LIBS=OFF \
  -DEXAMPLES_ENABLE_C=OFF \
  -DEXAMPLES_ENABLE_CXX=OFF

make -j$(nproc)
make install

echo ""
echo "SUNDIALS installed to: $SUNDIALS_PREFIX"
echo "CMake config location:"
find $SUNDIALS_PREFIX -name "SUNDIALSConfig.cmake" -o -name "sundials-config.cmake" 2>/dev/null | head -3
echo "Build completed: $(date)"
