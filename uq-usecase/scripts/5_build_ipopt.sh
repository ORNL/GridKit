#!/bin/bash
# IPOPT 3.14.14 with MUMPS (open-source linear solver via ThirdParty-Mumps).
# HSL solvers (MA27/MA57) are NOT used — libhsl.so does not exist on Kestrel.
# netlib-lapack module provides BLAS/LAPACK for MUMPS.

echo "=== Building IPOPT 3.14.14 + MUMPS ==="

DEPS_BASE=/nopt/nrel/apps/cpu_stack/software/gridkit
LLVM_PREFIX=$DEPS_BASE/llvm/16.0.6
IPOPT_PREFIX=$DEPS_BASE/ipopt/3.14.14
SRC_DIR=$DEPS_BASE/src
mkdir -p "$SRC_DIR" "$IPOPT_PREFIX"
cd $SRC_DIR

module load gcc-stdalone/13.1.0
module load cmake/3.29.2
module load netlib-lapack/3.11.0-gcc

GCC13_ROOT=$(dirname $(dirname $(which g++)))

export PATH=$LLVM_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_PREFIX/lib:$IPOPT_PREFIX/lib:$LD_LIBRARY_PATH
export CC=$LLVM_PREFIX/bin/clang
export CXX=$LLVM_PREFIX/bin/clang++
export FC=gfortran

# --- ThirdParty-Mumps (open-source MUMPS for IPOPT) ---
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps
./get.Mumps
./configure --prefix=$IPOPT_PREFIX \
  CC=$CC CXX=$CXX FC=$FC \
  CFLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  CXXFLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  LDFLAGS="-Wl,-rpath,${NETLIB_LAPACK_ROOT_DIR}/lib64"
make -j$(nproc)
make install
cd $SRC_DIR

# --- IPOPT 3.14.14 ---
wget https://github.com/coin-or/Ipopt/archive/refs/tags/releases/3.14.14.tar.gz
tar -xzf 3.14.14.tar.gz
cd Ipopt-releases-3.14.14
rm -rf build
mkdir build && cd build
../configure --prefix=$IPOPT_PREFIX \
  CC=$CC CXX=$CXX FC=$FC \
  CFLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  CXXFLAGS="--gcc-toolchain=${GCC13_ROOT}" \
  LDFLAGS="-Wl,-rpath,${NETLIB_LAPACK_ROOT_DIR}/lib64" \
  --with-mumps-cflags="-I$IPOPT_PREFIX/include/coin-or/mumps" \
  --with-mumps-lflags="-L$IPOPT_PREFIX/lib -lcoinmumps"

make -j$(nproc)
make install

echo ""
echo "IPOPT installed to: $IPOPT_PREFIX"
echo "Verify:"
ls -la $IPOPT_PREFIX/lib/libipopt*
ls -la $IPOPT_PREFIX/lib/libcoinmumps*
echo "Build completed: $(date)"
