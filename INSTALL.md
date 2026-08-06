# GridKit Installation Guide

GridKit has been built and tested on Linux and macOS. Windows users are
encouraged to use [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install),
which provides a Linux environment on Windows and is the easiest path to
a working GridKit build on that platform. Native Windows builds (Cygwin or
MSVC toolchain) are possible but not regularly tested.

## Contents

- [Requirements](#requirements)
- [Getting the source](#getting-the-source)
- [Building GridKit](#building-gridkit)
- [CMake options](#cmake-options)
- [Installing dependencies](#installing-dependencies)
  - [SUNDIALS](#sundials)
  - [Ipopt and HSL](#ipopt-and-hsl)
  - [Enzyme](#enzyme)
- [Using GridKit](#using-gridkit)
- [Building with Spack](#building-with-spack)

---

## Requirements

| Requirement | Version | Notes |
|---|---|---|
| CMake | >= 3.13 | |
| C++ compiler | C++20 | Clang or GCC |
| SUNDIALS | `develop` branch | Optional; disabled by default |
| SuiteSparse (KLU) | >= 7.x | Optional; needed for sparse solvers in SUNDIALS |
| Ipopt | >= 3.14 | Optional; disabled by default |
| HSL | >= 2015 | Optional; required by Ipopt for efficient linear solvers |
| gfortran | F90 | Optional; required when building with Ipopt/HSL |
| Enzyme | >= 0.0.206 | Optional; disabled by default; experimental |
| LLVM | >= 15 (tested with 16) | Optional; required when building with Enzyme |
| Doxygen | any | Optional; needed to build API documentation |

---

## Getting the source

```sh
git clone git@github.com:ORNL/GridKit.git
cd GridKit
git submodule update --init third-party/
```

The second step initializes bundled third-party headers (`magic-enum`,
`nlohmann-json`). Re-run it after pulling updates that change submodule
revisions.

---

## Building GridKit

The build directory must be outside the GridKit source tree. A sibling
directory is a common convention:

```sh
mkdir GridKit-build && cd GridKit-build
cmake ../GridKit \
  -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build . -j$(nproc)
cmake --install .
```

If dependencies are not in standard locations, point CMake to them:

```sh
cmake ../GridKit \
  -DSUNDIALS_DIR=/path/to/sundials/install \
  -DIPOPT_DIR=/path/to/ipopt/install \
  -DCMAKE_INSTALL_PREFIX=/path/to/install \
  -DGridKit_ENABLE_SUNDIALS=ON \
  -DGridKit_ENABLE_IPOPT=ON
```

### Running tests

On Linux:
```sh
cd build
ctest --output-on-failure
```

On Windows:
```powershell
cd build
ctest -C Debug --output-on-failure
```

### Testing the installation

After installing, you can verify that the installed CMake configuration is
correct and that GridKit can be consumed as a dependency by an external
project:

```sh
cd build
make test_install
```

`test_install` configures, builds, and runs a small consumer project located
in the GridKit install tree (`share/gridkit/examples/GridKitConsumer`). It
uses `find_package(GridKit)` against the installation prefix, so it exercises
the same code path a downstream user would follow. The build and test results
are placed under that consumer directory.

This target must be run **after** `cmake --install` because it depends on the
installed headers, libraries, and CMake config files.

On MSVC, instead use `cmake --build . --target test_install`.

---

## CMake options

| Option | Default | Description |
|---|---|---|
| `GridKit_ENABLE_SUNDIALS` | `OFF` | Build with SUNDIALS DAE integrators |
| `GridKit_ENABLE_IPOPT` | `OFF` | Build with Ipopt optimization solver |
| `GridKit_ENABLE_ENZYME` | `OFF` | Build with Enzyme automatic differentiation |
| `GridKit_ENABLE_ASAN` | `OFF` | Enable address sanitizer |
| `GridKit_ENABLE_UBSAN` | `OFF` | Enable undefined behavior sanitizer |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries |
| `CMAKE_INSTALL_PREFIX` | system default | Installation root |

Options may also be spelled with the `GridKit_ENABLE_*` prefix.

Dependency root directories:

| Variable | Points to |
|---|---|
| `SUNDIALS_DIR` | SUNDIALS install prefix |
| `IPOPT_DIR` | Ipopt install prefix |
| `ENZYME_DIR` | Enzyme install prefix |
| `SUITESPARSE_DIR` | SuiteSparse install prefix |

---

## Installing dependencies

### SUNDIALS

GridKit requires the SUNDIALS `develop` branch because it includes an adjoint
memory fix needed by the dynamic optimization examples. For sparse linear
algebra support (recommended), SUNDIALS must be built with KLU enabled, which
requires SuiteSparse.

**Install SuiteSparse first** (via the system package manager, or from source):

```sh
# Ubuntu/Debian
sudo apt install libsuitesparse-dev

# macOS with Homebrew
brew install suite-sparse
```

**Build SUNDIALS with KLU:**

```sh
git clone https://github.com/LLNL/sundials.git sundials-src
git -C sundials-src checkout develop
cmake -B sundials-build -S sundials-src \
  -DCMAKE_INSTALL_PREFIX=/path/to/sundials/install \
  -DENABLE_KLU=ON \
  -DKLU_INCLUDE_DIR=/path/to/suitesparse/include \
  -DKLU_LIBRARY_DIR=/path/to/suitesparse/lib \
  -DENABLE_MPI=OFF
cmake --build sundials-build -j$(nproc)
cmake --install sundials-build
```

Without KLU, omit the `ENABLE_KLU` and `KLU_*` flags.

---

### Ipopt and HSL

Ipopt requires a Fortran compiler for the HSL linear solver library. The
HSL source must be obtained separately from the
[HSL website](https://licences.stfc.ac.uk/product/coin-hsl) (free for
academic use).

Because HSL is Fortran library, it requires a feature complete Fortran
compiler. We were not able to build HSL with Flang, so to enable both, Ipopt
and Enzyme support, we used Clang as C/C++ and `gfortran` as Fortran compiler.
The compilers are set per-language via the `CC`, `CXX`, and `FC` environment
variables.

**Step 1 — Build HSL via ThirdParty-HSL:**

```sh
git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
cd ThirdParty-HSL
# Place the HSL source tarball in the same directory where you plan to
# install Ipopt and follow the instructions in the ThirdParty-HSL README
# (get_HSL.sh or manual unpack).
./configure CC=clang FC=gfortran --prefix=/path/to/hsl/install
make -j$(nproc)
make install
```

**Step 2 — Build Ipopt:**

```sh
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
./configure \
  CC=clang CXX=clang++ FC=gfortran \
  --prefix=/path/to/ipopt/install \
  --with-hsl-lflags="-L/path/to/hsl/install/lib -lcoinhsl" \
  --with-hsl-cflags="-I/path/to/hsl/install/include"
make -j$(nproc)
make install
```

Note: Ipopt's built-in tests will fail with mixed Clang/gfortran compilers.
This is because Ipopt build system does not add `-lgfortran` flag in that
configuration. Ignore these failed tests and cross your fingers.

**Step 3 — Build GridKit with Ipopt:**

```sh
cmake -B build -S . \
  -DGridKit_ENABLE_IPOPT=ON \
  -DIPOPT_DIR=/path/to/ipopt/install \
  -DSUNDIALS_DIR=/path/to/sundials/install \
  -DCMAKE_INSTALL_PREFIX=/path/to/gridkit/install
cmake --build build -j$(nproc)
cmake --install build
```

> **Note on the Fortran linker:** GridKit's CMakeLists sets
> `LINKER_LANGUAGE Fortran` for executables that use Ipopt+HSL. This
> causes `gfortran` to drive the final link step, which automatically adds
> the Fortran runtime libraries (`libgfortran`, `libquadmath`). No manual
> `-lgfortran` flags should be needed.

---

### Enzyme

Enzyme support is experimental. There is a known bug in Enzyme
that causes errors when compiler optimization is turned off with
`-O0` flag. A workaround is to use `Release` or `RelWithDebInfo`
build types in CMake.

Enzyme requires a matching LLVM installation (currently tested with LLVM 16)
and must be compiled with the same Clang used to build GridKit.

For instructions how to acquire Enzyme and (optionally) LLVM check Enzyme
[installation instructions](https://enzyme.mit.edu/Installation/).

After you installed Enzyme, you can build GridKit as

```sh
# Build GridKit with Enzyme
CC=clang CXX=clang++ FC=gcc cmake -B build -S . \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DGridKit_ENABLE_ENZYME=ON \
  -DENZYME_DIR=/path/to/enzyme/install \
  -DSUNDIALS_DIR=/path/to/sundials/install
cmake --build build -j 4
```

---

## Using GridKit

GridKit installs CMake config files that allow downstream projects to find and
link against GridKit with `find_package`. The config file automatically
re-finds GridKit's own dependencies (e.g. SUNDIALS) so consuming projects do
not need to locate them manually.

### Minimal CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.13)
project(MyApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)     # GridKit public headers require C++20
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GridKit REQUIRED)

add_executable(myapp myapp.cpp)
target_link_libraries(myapp GridKit::solvers_dyn GridKit::phasor_dynamics_components)
```

### Locating the GridKit installation

If GridKit was not installed to a system-wide prefix, tell CMake where to
find it via `GridKit_DIR` or `CMAKE_PREFIX_PATH`:

```sh
# Option A — point directly to the install prefix
cmake -B build -S . -DGridKit_DIR=/path/to/gridkit/install

# Option B — add the prefix to the search path
cmake -B build -S . -DCMAKE_PREFIX_PATH=/path/to/gridkit/install
```

### Available CMake targets

All targets are under the `GridKit::` namespace.

| Target | Description |
|---|---|
| `GridKit::solvers_dyn` | DAE integrators (requires SUNDIALS) |
| `GridKit::solvers_opt` | Optimization solvers (requires Ipopt) |
| `GridKit::phasor_dynamics_components` | Phasor dynamics model components |
| `GridKit::phasor_dynamics_core` | Phasor dynamics core infrastructure |
| `GridKit::phasor_dynamics_signal` | Signal node models |
| `GridKit::DenseMatrix` | Dense matrix utilities |
| `GridKit::sparse_matrix` | Sparse matrix utilities |
| `GridKit::Utilities` | General-purpose utility types |
| `GridKit::definitions` | Compile-time definitions and configuration |

Individual model component targets (generators, buses, loads, exciters,
governors, stabilizers) follow the pattern `GridKit::<component_name>` and
are also available; refer to the GridKit headers under `GridKit/Model/` for
the full list.

### C++20 requirement

GridKit's public headers use C++20 features. Consuming projects must enable
C++20 or later:

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
```

---

## Building with Spack

GridKit ships a Spack repository under `buildsystem/spack_repo` and includes
Spack as a submodule. This is the approach used for CI builds.

```sh
# Initialize submodules (if not done already)
git submodule update --init buildsystem/spack

# Set up Spack
source ./buildsystem/spack/share/spack/setup-env.sh

# Create and activate an environment
spack env create gridkit
spack env activate -p gridkit

# Register the GridKit Spack repo
spack repo add buildsystem/spack_repo/gridkit

# Tell Spack to use the local source tree
spack develop --path=$(pwd) gridkit@develop

# Detect available compilers
spack compiler find

# Add GridKit with desired variants, then build
spack add gridkit+sundials+ipopt+klu ^sundials@develop
spack concretize -f
spack install
spack env deactivate
```

Available Spack variants: `sundials`, `ipopt`, `klu`, `enzyme`, `asan`, `ubsan`.
Note that `+klu` requires `+sundials`.
