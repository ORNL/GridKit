#!/bin/bash
# SuiteSparse 7.8.2 — provides KLU sparse direct solver used by SUNDIALS.
# KLU is built by default; only SPEX is disabled (requires GMP/MPFR, not on Kestrel).

echo "=== Building SuiteSparse 7.8.2 ==="

DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
SUITESPARSE_PREFIX=$DEPS_BASE/suitesparse/7.8.2
SRC_DIR=$DEPS_BASE/src
mkdir -p "$SRC_DIR" "$SUITESPARSE_PREFIX"
cd $SRC_DIR

module load gcc-stdalone/13.1.0
module load cmake/3.29.2
module load netlib-lapack/3.11.0-gcc

export PATH=$LLVM_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_PREFIX/lib:$LD_LIBRARY_PATH
export CC=$LLVM_PREFIX/bin/clang
export CXX=$LLVM_PREFIX/bin/clang++

echo "Using Clang 16 for SuiteSparse build:"
$CC --version

wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.8.2.tar.gz
tar -xzf v7.8.2.tar.gz
cd SuiteSparse-7.8.2

rm -rf build
mkdir build && cd build
# SPEX is disabled: requires GMP and MPFR (arbitrary-precision libs) not on Kestrel;
# not needed — KLU is what SUNDIALS uses.
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$SUITESPARSE_PREFIX \
  -DCMAKE_C_COMPILER=$LLVM_PREFIX/bin/clang \
  -DCMAKE_CXX_COMPILER=$LLVM_PREFIX/bin/clang++ \
  -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;klu;umfpack;ldl;rbio" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSUITESPARSE_USE_OPENMP=ON

make -j$(nproc)
make install

echo ""
echo "SuiteSparse components installed:"
ls -la $SUITESPARSE_PREFIX/lib64/libamd*
ls -la $SUITESPARSE_PREFIX/lib64/libklu*
ls -la $SUITESPARSE_PREFIX/lib64/libsuitesparseconfig*
echo "Build completed: $(date)"
