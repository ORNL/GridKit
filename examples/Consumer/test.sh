#!/bin/bash

# This script gets executed in `make test_install`. 

# GridKit install prefix used for CMake dependency finding
export GridKit_DIR=${1}
echo ${GridKit_DIR}

# Locate source of the consumer test app
export INSTALL_BUILD_CONSUME=${GridKit_DIR}/share/gridkit/examples/GridKitConsumer
echo ${INSTALL_BUILD_CONSUME}

# Create build directory
mkdir -p ${INSTALL_BUILD_CONSUME}/build

rm -rf ${INSTALL_BUILD_CONSUME}/build/* &&

# Configure consumer project
cmake -B ${INSTALL_BUILD_CONSUME}/build \
    -S ${INSTALL_BUILD_CONSUME} \
    -DGridKit_DIR=${GridKit_DIR} &&

# Build and install
cmake --build ${INSTALL_BUILD_CONSUME}/build &&

cmake --install ${INSTALL_BUILD_CONSUME}/build &&

# Check
cd ${INSTALL_BUILD_CONSUME}/build &&

ctest --output-on-failure

exit $?
