# syntax=docker/dockerfile:1.7

FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

ARG LLVM_VERSION=16
ARG BUILD_JOBS=4
# GridKit currently requires SUNDIALS develop. Replace this with the full,
# tested commit SHA after the first successful end-to-end build.
ARG SUNDIALS_REF=develop
ARG IPOPT_REF=releases/3.14.19
# GridKit requires Enzyme >= 0.0.206 and warns that some versions are
# incompatible. Start with this candidate and pin its full SHA after testing.
ARG ENZYME_REF=v0.0.206

ENV PATH=/usr/lib/llvm-${LLVM_VERSION}/bin:${PATH} \
    CC=clang-${LLVM_VERSION} \
    CXX=clang++-${LLVM_VERSION} \
    FC=gfortran \
    CMAKE_PREFIX_PATH=/opt/gridkit-deps/sundials:/opt/gridkit-deps/ipopt:/opt/gridkit-deps/enzyme \
    PKG_CONFIG_PATH=/opt/gridkit-deps/ipopt/lib/pkgconfig

USER root

# Development tools, the LLVM 16 toolchain, SuiteSparse/KLU for SUNDIALS,
# and the open-source MUMPS linear solver used by Ipopt.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      git \
      build-essential \
      cmake \
      ninja-build \
      pkg-config \
      gfortran \
      clang-${LLVM_VERSION} \
      clangd-${LLVM_VERSION} \
      lldb-${LLVM_VERSION} \
      llvm-${LLVM_VERSION} \
      llvm-${LLVM_VERSION}-dev \
      llvm-${LLVM_VERSION}-tools \
      libclang-${LLVM_VERSION}-dev \
      libomp-${LLVM_VERSION}-dev \
      libzstd-dev \
      libsuitesparse-dev \
      libmumps-seq-dev \
      libmetis-dev \
      liblapack-dev \
      libopenblas-dev \
      gdb \
      python3 \
      python3-pip \
      python3-venv; \
    rm -rf /var/lib/apt/lists/*

# Build GridKit's required SUNDIALS revision against Ubuntu's SuiteSparse/KLU.
RUN set -eux; \
    git init /tmp/sundials; \
    git -C /tmp/sundials remote add origin https://github.com/LLNL/sundials.git; \
    git -C /tmp/sundials fetch --depth 1 origin "${SUNDIALS_REF}"; \
    git -C /tmp/sundials checkout --detach FETCH_HEAD; \
    cmake -S /tmp/sundials -B /tmp/sundials-build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=clang-${LLVM_VERSION} \
      -DCMAKE_CXX_COMPILER=clang++-${LLVM_VERSION} \
      -DCMAKE_INSTALL_PREFIX=/opt/gridkit-deps/sundials \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_STATIC_LIBS=OFF \
      -DSUNDIALS_ENABLE_KLU=ON \
      -DKLU_INCLUDE_DIR=/usr/include/suitesparse \
      -DSUNDIALS_ENABLE_KLU_CHECKS=ON \
      -DSUNDIALS_ENABLE_MPI=OFF \
      -DSUNDIALS_ENABLE_FORTRAN=OFF \
      -DSUNDIALS_ENABLE_C_EXAMPLES=OFF \
      -DSUNDIALS_ENABLE_CXX_EXAMPLES=OFF \
      -DSUNDIALS_ENABLE_FORTRAN_EXAMPLES=OFF \
      -DSUNDIALS_ENABLE_EXAMPLES_INSTALL=OFF; \
    cmake --build /tmp/sundials-build --target install --parallel "${BUILD_JOBS}"; \
    rm -rf /tmp/sundials /tmp/sundials-build

# Ubuntu 24.04's Ipopt package is too old for GridKit, so build a pinned
# Ipopt 3.14 release and link it to Ubuntu's sequential MUMPS package.
RUN set -eux; \
    git clone --depth 1 --branch "${IPOPT_REF}" \
      https://github.com/coin-or/Ipopt.git /tmp/ipopt; \
    mkdir -p /tmp/ipopt-build; \
    cd /tmp/ipopt-build; \
    CC=clang-${LLVM_VERSION} \
    CXX=clang++-${LLVM_VERSION} \
    FC=gfortran \
    F77=gfortran \
      /tmp/ipopt/configure \
        --prefix=/opt/gridkit-deps/ipopt \
        --with-mumps-cflags="-I/usr/include/mumps_seq" \
        --with-mumps-lflags="-ldmumps_seq" \
        --without-hsl \
        --without-asl \
        --disable-java \
        --enable-shared \
        --disable-static; \
    make -j"${BUILD_JOBS}"; \
    make install; \
    rm -rf /tmp/ipopt /tmp/ipopt-build

# Build Enzyme against the exact LLVM/Clang 16 toolchain used by GridKit.
RUN set -eux; \
    git clone --depth 1 --branch "${ENZYME_REF}" \
      https://github.com/EnzymeAD/Enzyme.git /tmp/enzyme; \
    cmake -S /tmp/enzyme/enzyme -B /tmp/enzyme-build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=clang-${LLVM_VERSION} \
      -DCMAKE_CXX_COMPILER=clang++-${LLVM_VERSION} \
      -DCMAKE_INSTALL_PREFIX=/opt/gridkit-deps/enzyme \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DLLVM_DIR=/usr/lib/llvm-${LLVM_VERSION}/lib/cmake/llvm \
      -DClang_DIR=/usr/lib/llvm-${LLVM_VERSION}/lib/cmake/clang \
      -DLLVM_EXTERNAL_LIT=/usr/lib/llvm-${LLVM_VERSION}/build/utils/lit/lit.py \
      -DENZYME_ENABLE_PLUGINS=ON \
      -DENZYME_CLANG=ON \
      -DENZYME_MLIR=OFF \
      -DENZYME_FORTRAN=OFF \
      -DENZYME_ENABLE_BENCHMARKS=OFF; \
    cmake --build /tmp/enzyme-build --target check-enzyme --parallel "${BUILD_JOBS}"; \
    cmake --build /tmp/enzyme-build --target install --parallel "${BUILD_JOBS}"; \
    printf '%s\n' 'define void @enzyme_smoke() { ret void }' \
      > /tmp/enzyme-smoke.ll; \
    /usr/lib/llvm-${LLVM_VERSION}/bin/opt \
      -load-pass-plugin=/opt/gridkit-deps/enzyme/lib/LLVMEnzyme-${LLVM_VERSION}.so \
      -passes=enzyme \
      -disable-output \
      /tmp/enzyme-smoke.ll; \
    clang-${LLVM_VERSION} \
      -fplugin=/opt/gridkit-deps/enzyme/lib/ClangEnzyme-${LLVM_VERSION}.so \
      -x c /dev/null -c -o /tmp/enzyme-smoke.o; \
    rm -rf /tmp/enzyme /tmp/enzyme-build /tmp/enzyme-smoke.ll /tmp/enzyme-smoke.o

# Make the source-built shared libraries available without LD_LIBRARY_PATH.
RUN set -eux; \
    printf '%s\n' \
      /opt/gridkit-deps/sundials/lib \
      /opt/gridkit-deps/ipopt/lib \
      /opt/gridkit-deps/enzyme/lib \
      > /etc/ld.so.conf.d/gridkit-deps.conf; \
    ldconfig

USER vscode
