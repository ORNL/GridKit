#!/bin/bash
# GridKit Build with Enzyme Support

echo "=== Building GridKit with Enzyme support ==="

# Set up environment
GRIDKIT_DIR=~/gridkit
DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
ENZYME_PREFIX=$DEPS_BASE/enzyme/0.0.206
SUITESPARSE_PREFIX=$DEPS_BASE/suitesparse/7.8.2
SUNDIALS_PREFIX=$DEPS_BASE/sundials/develop
IPOPT_PREFIX=$DEPS_BASE/ipopt/3.14.14
cd $GRIDKIT_DIR

# Load modules
module load gcc-stdalone/13.1.0
module load cmake/3.29.2
module load netlib-lapack/3.11.0-gcc
# NETLIB_LAPACK_ROOT_DIR is set by the module; expose it to the linker at build time
export LIBRARY_PATH=${NETLIB_LAPACK_ROOT_DIR}/lib64:${LIBRARY_PATH}

# Derive GCC 13 toolchain root so clang uses its libstdc++ headers
GCC13_ROOT=$(dirname $(dirname $(which g++)))

# Set up complete environment with all dependencies
export PATH=$LLVM_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_PREFIX/lib:$SUITESPARSE_PREFIX/lib64:$IPOPT_PREFIX/lib:$LD_LIBRARY_PATH
export CC=$LLVM_PREFIX/bin/clang
export CXX=$LLVM_PREFIX/bin/clang++

# Export OpenMP flags for GridKit
export OpenMP_C_FLAGS="-fopenmp"
export OpenMP_CXX_FLAGS="-fopenmp"
export OpenMP_C_LIB_NAMES="omp"
export OpenMP_CXX_LIB_NAMES="omp"
export OpenMP_omp_LIBRARY=$LLVM_PREFIX/lib/libomp.so

echo "GridKit build environment:"
echo "CC: $CC"
echo "CXX: $CXX"
echo "OpenMP lib: $OpenMP_omp_LIBRARY"

# GridKit repo is already at $GRIDKIT_DIR (no clone needed)
cd $GRIDKIT_DIR

# Ensure we are on our personal branch
git checkout isatkaus/uq-usecase
echo "GridKit branch and recent commits:"
git log --oneline -5
git submodule update --init third-party/

# Derive GCC 13 lib dir for rpath (so binaries find libstdc++ at runtime without module load)
GCC13_LIBDIR=$(${GCC13_ROOT}/bin/g++ --print-file-name=libstdc++.so | xargs dirname)

# Configure GridKit build (clean any stale cache)
rm -rf build
mkdir -p build && cd build
cmake .. \
  -DCMAKE_CXX_COMPILER=$LLVM_PREFIX/bin/clang++ \
  -DCMAKE_C_COMPILER=$LLVM_PREFIX/bin/clang \
  -DCMAKE_CXX_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_C_FLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  -DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,${GCC13_LIBDIR} -Wl,-rpath,${NETLIB_LAPACK_ROOT_DIR}/lib64" \
  -DCMAKE_SHARED_LINKER_FLAGS="-Wl,-rpath,${GCC13_LIBDIR} -Wl,-rpath,${NETLIB_LAPACK_ROOT_DIR}/lib64" \
  -DCMAKE_PREFIX_PATH="$IPOPT_PREFIX;$ENZYME_PREFIX;$SUITESPARSE_PREFIX;$LLVM_PREFIX" \
  -DSUNDIALS_DIR=$SUNDIALS_PREFIX/lib64/cmake/sundials \
  -DGridKit_ENABLE_SUNDIALS=ON \
  -DIPOPT_DIR=$IPOPT_PREFIX \
  -DGridKit_ENABLE_IPOPT=ON \
  -DENZYME_DIR=$ENZYME_PREFIX/lib/cmake/Enzyme \
  -DGridKit_ENABLE_ENZYME=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_VERBOSE_MAKEFILE=ON \
  -DOpenMP_C_FLAGS="$OpenMP_C_FLAGS" \
  -DOpenMP_CXX_FLAGS="$OpenMP_CXX_FLAGS" \
  -DOpenMP_omp_LIBRARY="$OpenMP_omp_LIBRARY"

# Build GridKit
echo "Building GridKit with $(nproc) cores..."
make -j$(nproc) -k

# Test results
echo ""
echo "=== BUILD RESULTS ==="
echo "Working Enzyme examples:"
find . -path "*/Enzyme/*" -type f -executable | head -5

echo "Phasor dynamics libraries:"
find . -path "*/PhasorDynamics/*" -name "*.so" | head -5

echo "Available make targets for phasor dynamics:"
make help | grep -i classical

echo ""
echo "GridKit build completed: $(date)"
echo ""
echo "To run all tests:"
echo "  ctest --test-dir $GRIDKIT_DIR/build --output-on-failure"
echo ""
echo "To test Enzyme only: ./examples/Enzyme/Standalone/EnzymeStandaloneScalarCheck"
